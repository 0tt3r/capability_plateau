"""Small-circuit PBC audit with dense operators, independent of resource counts.

Used for exact Clifford+T equivalence checks, not large-width resource sweeps.
The retained identity is U = C R, with Pauli rotations in R and final readouts
C^dagger Z_q C. A common frame preserves commutation within an input T layer.
"""
import numpy as np

I = np.eye(2, dtype=complex)
X = np.array([[0,1],[1,0]], dtype=complex)
Y = np.array([[0,-1j],[1j,0]])
Z = np.diag([1,-1]).astype(complex)
H = (X+Z)/np.sqrt(2)
S = np.diag([1,1j])
T = np.diag([1,np.exp(1j*np.pi/4)])
GATES = dict(H=H,S=S,SDG=S.conj().T,X=X,T=T,W=np.exp(1j*np.pi/4)*I)


def operator(name, qubits, width):
    if name == 'CNOT':
        c,t=qubits
        order=np.arange(2**width)
        target=order ^ (((order >> (width-1-c)) & 1) << (width-1-t))
        return np.eye(2**width,dtype=complex)[:,target]
    out=np.ones((1,1),dtype=complex)
    for q in range(width):
        out=np.kron(out,GATES[name] if q == qubits[0] else I)
    return out


def pauli_label(p, width):
    x=int(np.argmax(abs(p[:,0])))
    phase=p[x,0]
    labels=[]
    for q in range(width):
        bit=1 << (width-1-q)
        z=np.real(p[x ^ bit,bit]/phase) < 0
        labels.append('Y' if x & bit and z else 'X' if x & bit else 'Z' if z else 'I')
    sign=int(round(float(np.real(phase/(1j**labels.count('Y'))))))
    return ('+' if sign == 1 else '-')+''.join(labels)


def compile_layers(layers, width):
    """Layers contain only Clifford gates or only T gates; return an audit."""
    dim=2**width
    frame=np.eye(dim,dtype=complex)
    rotations=np.eye(dim,dtype=complex)
    original=np.eye(dim,dtype=complex)
    labels, groups, axes=[],[],[]
    max_commutator=0.
    for layer in layers:
        t_layer=[g for g in layer if g[0] == 'T']
        if t_layer and len(t_layer) != len(layer):
            raise ValueError('Cliffords and T gates must be in separate layers')
        group=[]
        for name, qubits in layer:
            gate=operator(name,qubits,width)
            original=gate @ original
            if name == 'T':
                # Avoid using the gate-conjugation expression as the equivalence check.
                z=np.ones((1,1),dtype=complex)
                for q in range(width):
                    z=np.kron(z,Z if q == qubits[0] else I)
                p=frame.conj().T @ z @ frame
                r=np.cos(np.pi/8)*np.eye(dim)-1j*np.sin(np.pi/8)*p
                rotations=r @ rotations
                group.append(p)
                labels.append(pauli_label(p,width))
            else:
                frame=gate @ frame
        if group:
            for i,p in enumerate(group):
                for q in group[:i]:
                    max_commutator=max(max_commutator,float(np.linalg.norm(p@q-q@p)))
            groups.append(len(group))
            axes.extend(group)
    compiled=frame @ rotations
    phase=np.vdot(compiled,original)
    phase=phase/abs(phase)
    error=float(np.linalg.norm(original-phase*compiled,ord=2))
    finals=[]
    for q in range(width):
        z=np.ones((1,1),dtype=complex)
        for j in range(width):
            z=np.kron(z,Z if j == q else I)
        finals.append(frame.conj().T @ z @ frame)
    return dict(unitary_error=error,max_within_layer_commutator=max_commutator,
                T_total=len(labels),T_layers=len(groups),layer_widths=groups,
                axes=labels,final_axes=[pauli_label(p,width) for p in finals],
                original=original,rotations=rotations,final_operators=finals)


def synthesized_quops_layers(width, ensemble, seed, tolerance=.05):
    """Actual gridsynth words, preserving simultaneous T layers across wires."""
    import mpmath
    from pygridsynth import gridsynth_gates
    from .quops import layers as layer_count, random_pairing
    rng=np.random.default_rng(seed)
    mpmath.mp.dps=50
    result=[]
    for index in range(layer_count(width,ensemble)):
        if index % 2:
            result.append([('CNOT',(int(c),int(t))) for c,t in random_pairing(width,rng)])
            continue
        streams=[]
        for q in range(width):
            axis=rng.choice(list('XYZ'))
            angle=mpmath.mpf(float(rng.uniform(0,2*np.pi)))
            word=list(reversed(str(gridsynth_gates(angle,mpmath.mpf(tolerance)))))
            before={'Z':[],'X':['H'],'Y':['SDG','H']}[axis]
            after={'Z':[],'X':['H'],'Y':['H','S']}[axis]
            streams.append(before+word+after)
        positions=[0]*width
        while any(positions[q] < len(streams[q]) for q in range(width)):
            cliffords,ts=[],[]
            for q in range(width):
                while positions[q] < len(streams[q]) and streams[q][positions[q]] != 'T':
                    cliffords.append((streams[q][positions[q]],(q,)))
                    positions[q]+=1
                if positions[q] < len(streams[q]):
                    ts.append(('T',(q,)))
                    positions[q]+=1
            if cliffords:
                result.append(cliffords)
            if ts:
                result.append(ts)
    return result
