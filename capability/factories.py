"""Magic-state factories as black boxes (draft Sec. 3.5): output error, footprint, and output period.

Each factory is used unchanged for both noise models, as the draft specifies. Output errors come from the
source papers' own circuit-noise simulations at p = 1e-3, not from the draft's P_SC law. Moving a state from
the factory output to the consuming patch is not charged anywhere in this model.

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
    cycles_per_output: float   # surface-code cycles between delivered T states


FACTORIES = [
    Factory('cultivation', 'cultivation', 2e-9, 1024, 100.),
    Factory('litinski', '(15-to-1)_{17,7,7}', 4.5e-8, 4620, 42.6),
    Factory('litinski', '(15-to-1)_{9,5,5} x (15-to-1)_{21,9,11}', 6.1e-10, 7780, 469.),
    Factory('litinski', '(15-to-1)^6_{13,5,5} x (20-to-4)_{23,11,13}', 1.4e-10, 43300, 130/4),
    Factory('litinski', '(15-to-1)^4_{13,5,5} x (20-to-4)_{27,13,15}', 2.6e-11, 46800, 157/4),
    Factory('litinski', '(15-to-1)^6_{11,5,5} x (15-to-1)_{25,11,11}', 2.7e-12, 30700, 82.5),
    Factory('litinski', '(15-to-1)^6_{13,5,5} x (15-to-1)_{29,11,13}', 3.3e-14, 39100, 97.5),
    Factory('litinski', '(15-to-1)^6_{17,7,7} x (15-to-1)_{41,17,17}', 4.5e-20, 73400, 128.),
]


def menu(families):
    return [f for f in FACTORIES if f.family in families]
