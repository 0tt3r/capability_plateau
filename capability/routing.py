"""SWAP-network routing of one QUOPS CNOT layer on a 2D grid (direct implementation, H1 and H3).

The w data qubits sit on an r x c grid (sites beyond w are empty). Before each CNOT layer we choose a
target placement in which every pair is adjacent, then reach it with rounds of disjoint nearest-neighbor
SWAPs. Two routers:

  snake  pairs occupy consecutive sites of a boustrophedon path; the permutation is realized by
         odd-even transposition sort along that path (at most r*c rounds). Simplest possible.
  rcr    pairs occupy horizontal dominoes chosen by minimum-distance assignment; the permutation is
         realized by row-column-row sorting (Alon, Chung, Graham 1994): a column phase, a row phase, and
         a column phase of odd-even transposition sorts, at most 2r + c rounds.

Placements map site -> token; tokens 0..w-1 are data qubits and tokens >= w are empty sites.
Each router returns the SWAP rounds (lists of swapped token pairs) and the new placement.
"""
import networkx as nx
import numpy as np
from scipy.optimize import linear_sum_assignment


def grid_shape(w):
    """Smallest rectangle with an even number of columns holding w qubits; ties go to the squarest."""
    best = None
    for c in range(2, w + 3, 2):
        r = -(-w//c)
        key = (r*c, abs(c - r))
        if best is None or key < best[0]:
            best = (key, (r, c))
    return best[1]


def snake_order(r, c):
    return [i*c + (j if i % 2 == 0 else c - 1 - j) for i in range(r) for j in range(c)]


def _odd_even_rounds(lines, key, placement):
    """Parallel odd-even transposition sort of tokens along each line (a list of sites) by key[token]."""
    rounds, parity, idle = [], 0, 0
    while idle < 2:
        swaps = []
        for line in lines:
            for i in range(parity, len(line) - 1, 2):
                a, b = line[i], line[i + 1]
                if key[placement[a]] > key[placement[b]]:
                    swaps.append((placement[a], placement[b]))
                    placement[a], placement[b] = placement[b], placement[a]
        if swaps:
            rounds.append(swaps)
            idle = 0
        else:
            idle += 1          # two consecutive quiet rounds (both parities) mean every line is sorted
        parity ^= 1
    return rounds


def route_snake(placement, pairs, r, c):
    placement = placement.copy()
    path = snake_order(r, c)
    where = np.empty(len(placement), int)
    for k, site in enumerate(path):
        where[placement[site]] = k
    w = 2*len(pairs)
    target = np.empty(len(placement), int)
    for slot, (a, b) in enumerate(sorted(pairs, key=lambda ab: where[ab[0]] + where[ab[1]])):
        first, second = sorted((a, b), key=lambda t: where[t])
        target[first], target[second] = 2*slot, 2*slot + 1
    holes = sorted(range(w, len(placement)), key=lambda t: where[t])
    for k, t in enumerate(holes):
        target[t] = w + k
    return _odd_even_rounds([path], target, placement), placement


def route_rcr(placement, pairs, r, c):
    placement = placement.copy()
    M, w = r*c, 2*len(pairs)
    row, col = np.empty(M, int), np.empty(M, int)
    for site, t in enumerate(placement):
        row[t], col[t] = divmod(site, c)
    # Target dominoes by minimum total Manhattan distance, each pair in its cheaper orientation.
    slots = [(i, 2*j) for i in range(r) for j in range(c//2)]
    cost = np.empty((len(pairs), len(slots)))
    orient = np.empty((len(pairs), len(slots)), bool)
    for k, (a, b) in enumerate(pairs):
        for m, (i, j) in enumerate(slots):
            ab = abs(row[a] - i) + abs(col[a] - j) + abs(row[b] - i) + abs(col[b] - j - 1)
            ba = abs(row[b] - i) + abs(col[b] - j) + abs(row[a] - i) + abs(col[a] - j - 1)
            cost[k, m], orient[k, m] = min(ab, ba), ba < ab
    dest = np.empty(M, int)
    used = set()
    for k, m in zip(*linear_sum_assignment(cost)):
        a, b = pairs[k]
        if orient[k, m]:
            a, b = b, a
        i, j = slots[m]
        dest[a], dest[b] = i*c + j, i*c + j + 1
        used.update((i*c + j, i*c + j + 1))
    free = [s for s in range(M) if s not in used]
    holes = list(range(w, M))
    if holes:
        hc = np.array([[abs(row[t] - s//c) + abs(col[t] - s % c) for s in free] for t in holes])
        for k, m in zip(*linear_sum_assignment(hc)):
            dest[holes[k]] = free[m]
    drow, dcol = dest//c, dest % c
    # Phase 1: route within columns so that each row then holds one token for every destination column.
    # The tokens form an r-regular bipartite multigraph (current column -> destination column), which
    # splits into r perfect matchings (Konig); matching k is sent to row k.
    bucket = {}
    for t in range(M):
        bucket.setdefault((col[t], dcol[t]), []).append(t)
    mid = np.empty(M, int)
    for k in range(r):
        g = nx.Graph()
        g.add_nodes_from((('L', j) for j in range(c)), bipartite=0)
        g.add_nodes_from((('R', j) for j in range(c)), bipartite=1)
        g.add_edges_from((('L', a), ('R', b)) for (a, b), ts in bucket.items() if ts)
        match = nx.bipartite.hopcroft_karp_matching(g, top_nodes=[('L', j) for j in range(c)])
        for j in range(c):
            t = bucket[j, match[('L', j)][1]].pop()
            mid[t] = k
    columns = [[i*c + j for i in range(r)] for j in range(c)]
    rows = [[i*c + j for j in range(c)] for i in range(r)]
    rounds = _odd_even_rounds(columns, mid, placement)
    rounds += _odd_even_rounds(rows, dcol, placement)
    rounds += _odd_even_rounds(columns, drow, placement)
    return rounds, placement


ROUTERS = dict(snake=route_snake, rcr=route_rcr)


def pairs_adjacent(placement, pairs, c):
    site = np.empty(len(placement), int)
    site[placement] = np.arange(len(placement))
    for a, b in pairs:
        (ia, ja), (ib, jb) = divmod(site[a], c), divmod(site[b], c)
        if abs(ia - ib) + abs(ja - jb) != 1:
            return False
    return True
