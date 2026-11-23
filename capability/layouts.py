"""Tile-counted versions of Litinski's two-qubit-column fast block (Fig. 13b).

The certificate checks area, connected ancilla space and access to both sides
of every two-qubit patch. It is not a syndrome-circuit or decoder simulation.
"""
from dataclasses import dataclass
from functools import lru_cache
from math import isqrt


def ceil_sqrt(n):
    r = isqrt(n)
    return r + (r*r < n)


@dataclass(frozen=True)
class FastLayout:
    width: int
    heights: tuple
    data: frozenset
    ancilla: frozenset
    patches: tuple

    @property
    def area(self):
        return len(self.data) + len(self.ancilla)

    def certificate(self):
        reached = set()
        pending = [next(iter(self.ancilla))]
        while pending:
            x, y = pending.pop()
            if (x, y) in reached:
                continue
            reached.add((x, y))
            pending.extend((x+dx, y+dy) for dx,dy in ((1,0),(-1,0),(0,1),(0,-1))
                           if (x+dx, y+dy) in self.ancilla and (x+dx, y+dy) not in reached)
        access = all((x+side, y+dy) in self.ancilla
                     for x,y in self.patches for side in (-1,1) for dy in (0,1))
        return dict(width=self.width, pair_columns=list(self.heights), tiles=self.area,
                    data_tiles=len(self.data), ancilla_tiles=len(self.ancilla),
                    disjoint=not bool(self.data & self.ancilla),
                    connected_ancilla=reached == set(self.ancilla),
                    both_pauli_edges_accessible=access,
                    algorithm_qubits=2*len(self.patches))


@lru_cache(maxsize=2048)
def fast_layout(w):
    """Near-square columns, with the last column shortened by whole pairs.

    With m=w/2 pairs, r=ceil(sqrt(m)), c=ceil(m/r), the count is
    2w+2(r+c)+1. The source's continuous 2w+sqrt(8w)+1 is exact on
    square pair arrays; rounding that expression up can omit one tile.
    """
    if not isinstance(w, int) or w < 2 or w % 2:
        raise ValueError('Fast blocks require a positive even integer width')
    pairs = w//2
    rows = ceil_sqrt(pairs)
    columns = (pairs+rows-1)//rows
    heights = (rows,)*(columns-1) + (pairs-rows*(columns-1),)
    data, ancilla, patches = set(), set(), []
    for j, height in enumerate(heights):
        x = 2*j+1
        ancilla.add((x, 0))
        for k in range(height):
            y = 2*k+1
            patches.append((x,y))
            data.update(((x,y),(x,y+1)))
    for j in range(columns+1):
        adjacent = heights[max(0,j-1):min(columns,j+1)]
        ancilla.update((2*j,y) for y in range(2*max(adjacent)+1))
    return FastLayout(w, heights, frozenset(data), frozenset(ancilla), tuple(patches))


def unit_data_tiles(w):
    """Square padding for the published Bell-pair preparation schedule."""
    return (2*ceil_sqrt(w)+1)**2


def ladder_route_depth(pairing):
    """Disjoint-path depth on one routing row, for an old-layout example.

    Interval overlap is a routing conflict, even for disjoint algorithm pairs.
    This tests this ladder realization, not every imaginable 2w-tile layout.
    """
    intervals = [(min(a,b),max(a,b)) for a,b in pairing]
    return max((sum(a <= x <= b for a,b in intervals)
                for x in range(max(b for _,b in intervals)+1)), default=0)
