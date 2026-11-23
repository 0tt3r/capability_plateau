"""Parallel PBC via the open-chain units of Litinski (2019), Sec. 5.

Each unit prepares one width-w T layer; U units execute U-1 layers per wave.
Unit preparation and correction readout overlap across waves. Square padding
supports the published Bell preparation schedule; endpoint units are charged
as full interior units. Factories are unchanged black boxes, followed by an
explicit batch conversion to post-corrected states. Delivery remains an exposed
interface assumption, as in the sequential model. This is a resource model,
not a compiled factory/decoder layout. See docs/parallel_pbc_audit.md.
"""
from math import ceil, log1p

import numpy as np

from . import surface
from .config import PBC_UNITS, WEEK_SECONDS
from .factories import menu
from .layouts import ceil_sqrt, unit_data_tiles
from .quops import ENSEMBLES, layer_counts
from .surface_pbc import batch_schedule


def pipeline_cycles(layers, unit_cycles, units, reaction_cycles, final_cycles):
    """Cold pipeline, bounded correction storage and an explicit final flush."""
    if layers < 1 or int(layers) != layers or np.any(np.asarray(units) < 2):
        raise ValueError('Need positive layer count and at least two units')
    slots = np.asarray(units)-1
    waves = (layers-1)//slots
    last = (layers-1)%slots+1
    return unit_cycles + waves*np.maximum(unit_cycles, slots*reaction_cycles) + last*reaction_cycles + final_cycles


def _validate(w, ensemble, ch, a):
    if not isinstance(w, int) or w < 2 or w % 2 or ensemble not in ENSEMBLES:
        raise ValueError('Need even integer QUOPS width')
    if a.schedule != 'pbc_units' or a.layout != 'fast_columns' or a.supply != 'buffered':
        raise ValueError('Units require fast_columns and buffered supply')
    if (not a.units or any(int(u) != u or u < 2 for u in a.units)
            or not np.isfinite(a.reaction_cycles) or a.reaction_cycles <= 0
            or a.pair_preparation_steps < 1 or a.correction_storage_factor < 1
            or a.delivery_cycles < 0 or a.feedforward_cycles < 0
            or ch.t_per_injection_steps < 1 or ch.executions < 1 or ch.cycle_seconds <= 0):
        raise ValueError('Invalid unit schedule settings')


def _costs(w, ensemble, hw, ch, a, d, f, nt, eps, counts, units):
    rotations = w*layer_counts(w, ensemble)[0]
    layers = layer_counts(w, ensemble)[0]*nt
    # Serial conversion of one batch at each independent factory interface.
    latency = f.batch_cycles + a.pair_preparation_steps*d
    service = ch.t_per_injection_steps*d+a.delivery_cycles
    consume = batch_schedule(w, counts, f.states_per_batch, latency, service, 0)
    bell_cycles = (ceil_sqrt(w)+3)*d
    unit_cycles = consume+bell_cycles
    slots = units-1
    waves = (layers-1)//slots
    last = (layers-1)%slots+1
    reaction = a.reaction_cycles+a.feedforward_cycles
    cycles = unit_cycles + waves*np.maximum(unit_cycles, slots*reaction) + last*reaction + w*d
    data = unit_data_tiles(w)
    correction = 2*w*a.correction_storage_factor
    # Per batch output: magic state, conversion/correction qubit and access tile.
    buffers = 3*counts*f.states_per_batch
    tiles = units*(data+correction+buffers)
    patches = tiles*cycles*(-log1p(-surface.p_sc(d,hw,ch)))
    magic = rotations*nt*(-log1p(-f.error))
    synthesis = ch.synthesis_scale*rotations*(eps**2 if ch.synthesis_charge == 'squared' else eps)
    if ch.synthesis_budget is None:
        ok = patches+magic+synthesis <= ch.threshold
    else:
        ok = (patches+magic <= ch.threshold) & (synthesis <= ch.synthesis_budget)
    runs = np.minimum(np.floor(WEEK_SECONDS/(cycles*ch.cycle_seconds)),ch.executions)
    copies = np.ceil(ch.executions/np.maximum(runs,1))
    per_copy = tiles*2*d*d+units*counts*f.qubits
    return dict(qubits=copies*per_copy,qubits_per_copy=per_copy,cycles=cycles,
                unit_cycles=unit_cycles,unit_consumption_cycles=consume,bell_cycles=bell_cycles,
                tiles=tiles,data_tiles_per_unit=data,correction_tiles_per_unit=correction,
                buffer_tiles_per_unit=buffers,copies=copies,runs_per_copy=runs,
                loss_patches=patches,loss_magic=magic,loss_synthesis=synthesis,
                passes=ok & (runs >= 1),T_layers=layers)


def evaluate_design(w, ensemble, hw, ch, d, f, nt, factories_per_unit, units, architecture=PBC_UNITS):
    _validate(w,ensemble,ch,architecture)
    if (d < 3 or d % 2 == 0 or factories_per_unit < 1 or units < 2
            or int(d) != d or int(factories_per_unit) != factories_per_unit or int(units) != units):
        raise ValueError('Invalid distance/factory/unit count')
    d,factories_per_unit,units=int(d),int(factories_per_unit),int(units)
    ix=np.flatnonzero(surface.NT == nt)
    if not len(ix):
        raise ValueError('T count outside menu')
    eps=float(surface.EPS[ix[0]])
    x=_costs(w,ensemble,hw,ch,architecture,d,f,nt,eps,factories_per_unit,units)
    row={k:float(v) for k,v in x.items() if k != 'passes'}
    for k in ('tiles','data_tiles_per_unit','correction_tiles_per_unit','buffer_tiles_per_unit',
              'copies','runs_per_copy','T_layers'):
        row[k]=int(row[k])
    row.update(passes=bool(x['passes']),distance=d,factory=f.name,
               factories=units*factories_per_unit,factories_per_unit=factories_per_unit,
               units=units,parallel_layers=units-1,T_per_rotation=nt,
               T_total=w*row['T_layers'],eps=eps,run_seconds=row['cycles']*ch.cycle_seconds,
               padded_width=ceil_sqrt(w)**2,cnot_cycles=0,final_measurement_cycles=w*d,
               schedule='pbc_units',layout='fast_columns',
               schedule_status='published unit protocol with declared factory interface and storage assumptions')
    return row


def min_qubits(w, ensemble, hw, ch, architecture=PBC_UNITS):
    _validate(w,ensemble,ch,architecture)
    nt,eps=surface.NT[:,None,None],surface.EPS[:,None,None]
    units=np.array(architecture.units)[None,None,:]
    best=None
    for f in menu(ch.factories):
        for d in range(3,101,2):
            minimum_tiles=unit_data_tiles(w)+2*w*architecture.correction_storage_factor+3*f.states_per_batch
            if best is not None and min(architecture.units)*(minimum_tiles*2*d*d+f.qubits) >= best['qubits']:
                break
            service=ch.t_per_injection_steps*d+architecture.delivery_cycles
            latency=f.batch_cycles+architecture.pair_preparation_steps*d
            saturation=max(1,ceil(latency/(f.states_per_batch*service)))
            counts=np.arange(1,saturation+1)[None,:,None]
            x=_costs(w,ensemble,hw,ch,architecture,d,f,nt,eps,counts,units)
            q=np.where(x['passes'],x['qubits'],np.inf)
            index=np.unravel_index(np.argmin(q),q.shape)
            if np.isfinite(q[index]) and (best is None or q[index] < best['qubits']):
                i,j,k=index
                best=evaluate_design(w,ensemble,hw,ch,d,f,int(surface.NT[i]),j+1,
                                     int(units[0,0,k]),architecture)
    return best
