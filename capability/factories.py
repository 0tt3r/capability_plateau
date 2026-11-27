"""Magic-state factories as black boxes (draft Sec. 3.5): output error, footprint, and output period.

Each factory is used unchanged for both noise models, as the draft specifies. Output errors come from the
source papers' own circuit-noise simulations at p = 1e-3, not from the draft's P_SC law.
Sequential PBC adds output buffers with patch-cycle error
and exposes extra delivery cycles; its default zero delivery time remains an uncalibrated interface.

  cultivation  Output error 2e-9 per T state: source estimate (Gidney, Shutty, and Jones, arXiv:2409.17595v1,
               Table 1 and Figs. 1, 15, 16, uniform depolarizing noise at p = 1e-3). The source state ends in
               a grafted matchable code, so transfer, storage, and consumption are outside the estimate.
               Footprint 1024 qubits and period 100 cycles per delivered state: our supply assumptions (the
               workshop values), a packing of 1.0e5 qubit-cycles per state against the source's expected
               volume of about 6.6e4 qubit-rounds (Fig. 15, with early rejection and sharing). No cultivation
               factory has been laid out or scheduled for this model.
  litinski     Litinski, Quantum 3, 205 (2019), arXiv:1905.06903v3, Table 1, its seven p = 1e-3 rows: output
               error per state, physical qubits, and code cycles per round as reported. Period = cycles per
               round / states per round (four for 20-to-4), which assumes rounds run back to back.
"""
from dataclasses import dataclass


@dataclass(frozen=True)
class Factory:
    family: str
    name: str
    error: float               # logical error per delivered T state
    qubits: int
    cycles_per_output: float   # average surface-code cycles per delivered T state
    states_per_batch: int = 1  # PBC supply schedule retains multi-output batches

    @property
    def batch_cycles(self):
        """Declared cold-start batch latency = published period times batch size.

        This deterministic service model omits random factory rejection stalls.
        """
        return self.cycles_per_output*self.states_per_batch


# Fold-transversal surface-code cultivation (Sahay, Tsai, Chang, arXiv:2509.05212v2, PRX Quantum 7, 033006), SD6 noise at
# p = 1e-3, which uses non-local operations and is therefore used only on all-to-all hardware. Expected spacetime volumes
# per accepted state are read from its Fig. 2b in qubit x gate steps: 3.6e4 for f = 5 at an error of 1e-9 (d_fin = 15) and
# 5.7e3 for f = 3 at 2e-6 (d_fin = 13). [choice] One gate step is one t_op, so a surface-code cycle is 6 steps, and the
# same packing factor as the 'cultivation' row (1.0e5 qubit-cycles against the source's 6.6e4, 1.55) is applied.
FOLD_PACKING = 1024*100/6.6e4
FOLD_F5 = dict(error=1e-9, qubits=2*15*15 - 1, volume=3.6e4*FOLD_PACKING/6)     # qubit-cycles per delivered state
FOLD_F3 = dict(error=2e-6, qubits=2*13*13 - 1, volume=5.7e3*FOLD_PACKING/6)     # an input state for distillation

FACTORIES = [
    Factory('cultivation', 'cultivation', 2e-9, 1024, 100.),
    Factory('fold_cultivation', 'fold-transversal cultivation', FOLD_F5['error'], int(FOLD_F5['qubits']),
            FOLD_F5['volume']/FOLD_F5['qubits']),
    Factory('litinski', '(15-to-1)_{17,7,7}', 4.5e-8, 4620, 42.6),
    Factory('litinski', '(15-to-1)_{9,5,5} x (15-to-1)_{21,9,11}', 6.1e-10, 7780, 469.),
    Factory('litinski', '(15-to-1)^6_{13,5,5} x (20-to-4)_{23,11,13}', 1.4e-10, 43300, 130/4, 4),
    Factory('litinski', '(15-to-1)^4_{13,5,5} x (20-to-4)_{27,13,15}', 2.6e-11, 46800, 157/4, 4),
    Factory('litinski', '(15-to-1)^6_{11,5,5} x (15-to-1)_{25,11,11}', 2.7e-12, 30700, 82.5),
    Factory('litinski', '(15-to-1)^6_{13,5,5} x (15-to-1)_{29,11,13}', 3.3e-14, 39100, 97.5),
    Factory('litinski', '(15-to-1)^6_{17,7,7} x (15-to-1)_{41,17,17}', 4.5e-20, 73400, 128.),
]


# Layouts of the arXiv:2610.06535 protocols that set a point on the extended surface-code curve (Sec. 4.10 appendix).
JWS_WINNERS = ('[[15,1,3]], cultivated d_a=7, d=(19,5,5)',
               '[[15,1,3]], cultivated d_a=7, d=(21,5,5)',
               '[[15,1,3]], cultivated d_a=9, d=(19,5,5)',
               '[[15,1,3]], cultivated d_a=9, d=(21,5,5)',
               '[[15,1,3]], cultivated d_a=9, d=(23,5,5)',
               '[[15,1,3]], faulty T, d=(11,5,5)',
               '[[15,1,3]], faulty T, d=(13,5,5)',
               '[[15,1,3]], faulty T, d=(15,5,5)',
               '[[15,1,3]], faulty T, d=(17,5,5)',
               '[[15,1,3]], faulty T, d=(9,5,5)',)


CATALOGUE_WINNERS = ('[[14,1,2]] T0, faulty T, d=(11,5,5)',
                     '[[14,1,2]] T0, faulty T, d=(7,5,5)',
                     '[[14,1,2]] T0, faulty T, d=(9,5,5)',
                     '[[8,3,2]] CCZ012, faulty T, d=(7,5,5)',
                     '[[8,3,2]] CCZ012, faulty T, d=(9,5,5)',)


# All-to-all hardware (fold-transversal f = 3 inputs allowed): the layouts that set a point on the exact optimized
# surface-code curves of H2 and H4 (scripts/find_winners.py), together with fold-transversal cultivation.
FOLD_JWS_WINNERS = ('[[15,1,3]], fold f=3, d=(19,5,5)',
                    '[[15,1,3]], fold f=3, d=(21,5,5)',
                    '[[15,1,3]], fold f=3, d=(23,5,5)',)
FOLD_CATALOGUE_WINNERS = ('[[14,1,2]] T0, fold f=3, d=(19,5,5)',
                          '[[14,1,2]] T0, fold f=3, d=(21,5,5)',
                          '[[14,2,2]] T0,T1, fold f=3, d=(19,5,5)',
                          '[[14,2,2]] T0,T1, fold f=3, d=(21,5,5)',
                          '[[14,2,2]] T0,T1, fold f=3, d=(23,5,5)',
                          '[[8,3,2]] CCZ012, fold f=3, d=(19,5,5)',)


def catalogue(names=None):
    """Surface-code layouts of the distance-2 and distance-3 catalogue factories (scripts/catalogue_factories.py),
    per delivered T state, with the T states of one run as its batch, as for the multi-output rows of Litinski.
    names, if given, keeps only those designs."""
    import json
    from pathlib import Path
    import os
    path = Path(os.environ.get('CATALOGUE_JSON', Path(__file__).resolve().parents[1]/'results/catalogue_factories.json'))
    if not path.exists():
        return []
    return [Factory('catalogue', r['name'], r['error'], int(r['qubits']), r['cycles_per_T'], int(r['T_per_run']))
            for r in json.loads(path.read_text())['designs'] if names is None or r['name'] in names]


def menu(families):
    """'jws' and 'catalogue' are the full sets; 'jws_winners' and 'catalogue_winners' keep the designs that set a point
    on the extended surface-code curve. Pinnacle and Frozen take the full sets and apply their own winner lists."""
    out = [f for f in FACTORIES if f.family in families]
    if 'jws' in families or 'jws_winners' in families:   # surface-code layouts of the arXiv:2610.06535 protocols
        from .surface_factories import factories
        out += [f for f in factories() if 'jws' in families or f.name in JWS_WINNERS]
    if 'jws_fold_winners' in families:                     # all-to-all hardware (fold-transversal f = 3 inputs allowed)
        from .surface_factories import factories
        out += [f for f in factories() if f.name in FOLD_JWS_WINNERS]
    if 'catalogue' in families or 'catalogue_winners' in families:
        out += catalogue(None if 'catalogue' in families else CATALOGUE_WINNERS)
    if 'catalogue_fold_winners' in families:
        out += catalogue(FOLD_CATALOGUE_WINNERS)
    return out
