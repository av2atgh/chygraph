"""Minimum vertex cover on the clique and merged chygraphs.

Chapters 14 to 16 measure one model, the Ising model, on three classes of
network.  This asks whether the meta-complex repair carries over to the model
Chapter 11 is about: minimum vertex cover, a hard constraint rather than a soft
coupling, which is where message passing is usually harder rather than easier.

The encoding.  A variable per atom, $x_v = 1$ when the atom is in the cover.
Every bond must be covered, so a complex carries the indicator that no bond
lying inside it has both ends free:

    f_A(x_A) = 1 if every bond {i,j} inside A has x_i + x_j >= 1, else 0,

which for a clique is the familiar ``at most one member may be left out'' of
Sec. 11.2 and for a merged meta-complex is the general statement -- a union of
cliques is not a clique, and only the bonds actually present may be demanded.
A chemical potential is attached to each atom as a complex of one,
f_v(x_v) = e^{-mu x_v}, so that large mu drives the cover to its minimum.

Everything then runs through ``ChygraphBP`` of ``cavity_clique.py``, the same
solver Chapters 14 and 16 use, and is compared against enumeration over all
2^n assignments: the exact ln Z at the same mu, and the exact minimum cover.

    python probe/cover_chygraph.py
"""

import json
import sys
from itertools import combinations
from pathlib import Path

import networkx as nx
import numpy as np
from scipy.special import logsumexp

sys.path.insert(0, str(Path.home() / 'av2atg' / 'chygraph' / 'src'))
sys.path.insert(0, str(Path.home() / 'av2atg' / 'chygraph_statmech'
                      / 'book' / 'figs'))
sys.path.insert(0, str(Path.home() / 'av2atg' / 'chygraph_statmech' / 'probe'))

from cavity_clique import ChygraphBP  # noqa: E402
from merge import merge_closure  # noqa: E402

RESULTS = Path(__file__).parent / 'results'
OUT = RESULTS / 'cover_chygraph.json'
NEG = -1e6            # a forbidden assignment, in logs
MU = (1.0, 2.0)       # chemical potentials


def cover_factors(complexes, edges, mu):
    """Log-tables for the constraint complexes and the unary weights."""
    bonds = {tuple(sorted(e)) for e in edges}
    fam, logf = [], []
    for c in complexes:
        c = tuple(sorted(c))
        inside = [(i, j) for i, j in combinations(range(len(c)), 2)
                  if (c[i], c[j]) in bonds]
        t = np.zeros((2,) * len(c))
        if inside:
            idx = np.indices(t.shape)
            bad = np.zeros(t.shape, bool)
            for i, j in inside:
                bad |= (idx[i] == 0) & (idx[j] == 0)
            t[bad] = NEG
        fam.append(c)
        logf.append(t)
    for v in sorted({v for c in complexes for v in c}):
        fam.append((v,))
        logf.append(np.array([0.0, -mu]))     # x_v = 0 free, x_v = 1 costs mu
    return fam, logf


def exact(edges, n, mu):
    """ln Z and the minimum cover, by enumeration over all 2^n assignments."""
    bits = np.arange(2 ** n)
    x = np.stack([(bits >> v) & 1 for v in range(n)])       # n x 2^n
    ok = np.ones(2 ** n, bool)
    for i, j in edges:
        ok &= (x[i] | x[j]).astype(bool)
    size = x.sum(axis=0)
    lz = float(logsumexp(-mu * size[ok]))
    return lz, int(size[ok].min())


def solve(G, n, mu):
    edges = sorted({tuple(sorted(e)) for e in G.edges()})
    cliques = [sorted(c) for c in nx.find_cliques(G)]
    merged = [sorted(c) for c in merge_closure(
        [frozenset(c) for c in cliques])[0]]
    lz, mincov = exact(edges, n, mu)
    out = {'n': n, 'mu': mu, 'exact_lnZ': lz, 'min_cover': mincov,
           'n_cliques': len(cliques), 'n_meta': len(merged)}
    for tag, fam in (('clique', cliques), ('merged', merged)):
        f, lf = cover_factors(fam, edges, mu)
        bp = ChygraphBP(f, damping=0.5, edges=edges, log_factors=lf).run()
        out[tag] = bp.log_Z() - lz
        out[tag + '_res'] = bp.residual
    return out


def validate():
    """Exact where the chygraph is treelike, as for the Ising model."""
    print('validation, vertex cover on treelike chygraphs:')
    cases = [('chain of 6', nx.path_graph(6)),
             ('star, 7 leaves', nx.star_graph(7)),
             ('two triangles at a node', nx.Graph(
                 [(0, 1), (1, 2), (2, 0), (2, 3), (3, 4), (4, 2)])),
             ('tree of triangles', nx.Graph(
                 [(0, 1), (1, 2), (2, 0), (2, 3), (3, 4), (4, 2),
                  (0, 5), (5, 6), (6, 0)]))]
    ok = True
    for name, G in cases:
        for mu in MU:
            r = solve(G, G.number_of_nodes(), mu)
            good = abs(r['clique']) < 1e-8
            ok &= good
            print(f"  {name:<24} mu={mu}  min cover={r['min_cover']}  "
                  f"error={r['clique']:+.2e}  {'ok' if good else 'FAIL'}")
    assert ok, 'the cover recursion is not exact on a treelike chygraph'
    print('  exact wherever complexes meet in at most one node\n')


def main():
    validate()
    from gbp_cliques import SIZES, instance
    from gbp_real import load
    from merge import karrer_graph

    rows, kbar = [], {n: k for n, k in SIZES}
    seen = set()
    for r in json.load(open(RESULTS / 'gbp_cliques.json')):
        if (r['n'], r['seed']) in seen:
            continue
        seen.add((r['n'], r['seed']))
        G = instance(r['n'], kbar[r['n']], r['seed'])
        for mu in MU:
            rows.append(dict(ensemble='hyperbolic', **solve(G, r['n'], mu)))

    seen = set()
    for r in json.load(open(RESULTS / 'gbp_karrer.json'))['runs']:
        if (r['n'], r['seed']) in seen:
            continue
        seen.add((r['n'], r['seed']))
        G, _ = karrer_graph(r['n'], r['s_mean'], {3: r['triangles']}, r['seed'])
        G.remove_edges_from(nx.selfloop_edges(G))
        for mu in MU:
            rows.append(dict(ensemble='karrer', **solve(G, r['n'], mu)))

    cache, seen = {}, set()
    for r in json.load(open(RESULTS / 'gbp_real.json')):
        if (r['network'], r['ego']) in seen:
            continue
        seen.add((r['network'], r['ego']))
        G = cache.setdefault(r['key'], load(r['key']))
        ego = sorted([r['ego']] + list(G[r['ego']]), key=str)
        H = nx.convert_node_labels_to_integers(G.subgraph(ego))
        for mu in MU:
            rows.append(dict(ensemble='real', network=r['network'],
                             **solve(H, H.number_of_nodes(), mu)))

    json.dump(rows, OUT.open('w'), indent=1)
    summarise(rows)
    print('\nwrote', OUT)


def summarise(rows):
    for ens in ('hyperbolic', 'karrer', 'real'):
        sel = [r for r in rows if r['ensemble'] == ens]
        c = [abs(r['clique']) for r in sel]
        m = [abs(r['merged']) for r in sel]
        ex = sum(1 for x in m if x < 1e-8)
        print(f'\n{ens}: {len(sel)} runs')
        print(f'  clique chygraph: {min(c):.2e} to {max(c):.2e}, '
              f'median {np.median(c):.3f}')
        print(f'  merged chygraph: {min(m):.2e} to {max(m):.2e}, '
              f'median {np.median(m):.3f}, exact on {ex}/{len(sel)}')


if __name__ == '__main__':
    main()
