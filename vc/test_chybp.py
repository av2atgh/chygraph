"""Sanity: regions glued into a tree -> Bethe energy = exact minimum cover."""
import itertools, sys
import numpy as np
sys.path.insert(0, '.')
from chybp import ChyBP

def exact(n, edges):
    best = n
    for bits in range(1 << n):
        if all((bits >> a) & 1 or (bits >> b) & 1 for a, b in edges):
            best = min(best, bin(bits).count('1'))
    return best

rng = np.random.default_rng(1)
bad = 0
for trial in range(200):
    # regions of random size 1..5 with random internal edges, glued by a tree of external edges
    sizes = rng.integers(1, 6, size=rng.integers(2, 5))
    nodes, blocks, edges = 0, [], []
    for k, s in enumerate(sizes):
        vs = list(range(nodes, nodes + s)); nodes += s
        blocks += [k] * s
        for a, b in itertools.combinations(vs, 2):
            if rng.random() < 0.6: edges.append((a, b))
    starts = np.cumsum([0] + list(sizes[:-1]))
    for k in range(1, len(sizes)):          # tree between regions, random endpoints
        j = rng.integers(k)
        a = starts[k] + rng.integers(sizes[k]); b = starts[j] + rng.integers(sizes[j])
        edges.append((int(a), int(b)))
    if nodes > 16: continue
    ex = exact(nodes, edges)
    bp = ChyBP(nodes, np.array(edges), blocks, smax=10)
    it, diff = bp.iterate(iters=200, damping=0.3)
    E = bp.bethe_energy(); dec, _ = bp.decimate()
    ok = abs(E - ex) < 1e-6 and dec >= ex
    if not ok:
        bad += 1
        print('MISMATCH', sizes, 'exact', ex, 'bethe', round(E, 4), 'decimated', dec, 'iters', it, diff)
print('trials done, mismatches:', bad)
