"""Leaf removal generalised to a chygraph, for minimum vertex cover.

Chapters 14 to 16 measure the Ising model, where the question is how much of
ln Z the messages lose to overlap.  Vertex cover asks something different.  It
has an exact combinatorial algorithm -- leaf removal, Chapter 11's subject --
which either consumes the graph and returns a provably minimum cover, or stops
on a core and returns nothing.  So the question here is not how large the error
is but whether the algorithm finishes.

The generalisation.  Graph leaf removal takes the neighbour of a degree-one
vertex, because covering that edge costs one vertex either way and the
neighbour may cover others besides.  On a chygraph the same move is one level
up: if a complex A holds an atom v belonging to no other complex, then v is
needed nowhere else, so leaving v out and covering A minus v is optimal for A
and costs nothing anywhere else.  Cover it, delete every atom of A, repeat.
What survives, still constrained by a surviving complex, is the core.

Two things this buys over the hard-field cavity of Eq. (11.4).  It is exact
wherever it consumes the chygraph, overlap or no overlap -- two triangles
sharing an edge are consumed and answered correctly, which no message-passing
scheme in Part IV manages.  And on disjoint triangles it returns the true 2/3
where Sec. 10.4's hard-field map returns 1/2, because it never discards the
O(1) part of the field: there is no field to discard.

    python probe/cover_leafremoval.py
"""

import json
import sys
from itertools import combinations
from pathlib import Path

import networkx as nx
import numpy as np

sys.path.insert(0, str(Path.home() / 'av2atg' / 'chygraph' / 'src'))
sys.path.insert(0, str(Path.home() / 'av2atg' / 'chygraph_statmech'
                      / 'book' / 'figs'))
sys.path.insert(0, str(Path.home() / 'av2atg' / 'chygraph_statmech' / 'probe'))
sys.path.insert(0, str(Path.home() / 'av2atg' / 'computational_complexity'
                      / 'code'))

import leafremoval as lr  # noqa: E402  Ch. 11's own routine

from merge import merge_closure  # noqa: E402

RESULTS = Path(__file__).parent / 'results'
OUT = RESULTS / 'cover_leafremoval.json'


def _min_cover_induced(atoms, bonds, forbid=()):
    """Minimum cover of the subgraph induced on `atoms`, by enumeration."""
    atoms = sorted(atoms)
    idx = {v: i for i, v in enumerate(atoms)}
    inside = [(idx[u], idx[v]) for u, v in bonds
              if u in idx and v in idx]
    best, arg = len(atoms) + 1, None
    for m in range(1 << len(atoms)):
        if all((m >> i) & 1 or (m >> j) & 1 for i, j in inside):
            c = bin(m).count('1')
            if c < best:
                best, arg = c, m
    return {atoms[i] for i in range(len(atoms)) if (arg >> i) & 1}


def leaf_removal(complexes, n, bonds):
    """Cover what the chygraph forces, to a fixed point.

    Two moves, and only two, because only these two are safe.

    A *clique* complex A holding an atom v of chy-degree one: any cover leaves
    out at most one member of A, and if it leaves out some w other than v it
    contains v, so exchanging the two gives a cover of the same size leaving out
    v.  Leaving out v is therefore optimal, and A minus v goes into the cover.
    This is graph leaf removal one level up.

    An *isolated* complex, every atom of chy-degree one: it shares nothing with
    the rest, so its interior is enumerated and solved exactly.

    No move is made on a non-clique complex that is still attached to the rest.
    The exchange argument fails there --- a merged meta-complex is a union of
    cliques, and its minimum cover can be smaller than one less than its size,
    so covering all but one atom would overcount.  Getting this wrong is caught
    by the assertion in `run`.
    """
    cx = [set(c) for c in complexes if len(c) >= 2]
    bondset = {tuple(sorted(e)) for e in bonds}

    def is_clique(c):
        return all(tuple(sorted((u, v))) in bondset
                   for u, v in combinations(sorted(c), 2))

    cover = set()
    while True:
        deg = {}
        for c in cx:
            for v in c:
                deg[v] = deg.get(v, 0) + 1
        move = None
        for i, c in enumerate(cx):
            if all(deg[v] == 1 for v in c):
                move = (i, _min_cover_induced(c, bondset))
                break
        if move is None:
            for i, c in enumerate(cx):
                if not is_clique(c):
                    continue
                for v in c:
                    if deg[v] == 1:
                        move = (i, set(c) - {v})
                        break
                if move:
                    break
        if move is None:
            break
        i, take = move
        gone = set(cx[i])
        cover |= take
        cx = [c - gone for c in cx]
        cx = [c for c in cx if len(c) >= 2]
    core = set().union(*cx) if cx else set()
    return cover, core, cx


def exact_min_cover(G):
    """Minimum cover by enumeration; n <= 20 throughout Part IV."""
    n = G.number_of_nodes()
    edges = [(u, v) for u, v in G.edges()]
    bits = np.arange(2 ** n)
    x = [((bits >> v) & 1).astype(bool) for v in range(n)]
    ok = np.ones(2 ** n, bool)
    for u, v in edges:
        ok &= x[u] | x[v]
    size = np.zeros(2 ** n, np.int16)
    for v in range(n):
        size += x[v]
    return int(size[ok].min())


def is_cover(G, cover):
    return all(u in cover or v in cover for u, v in G.edges())


def ordinary_leaf_removal(G):
    """Chapter 11's leaf removal on the graph: core size, 0 when consumed."""
    n = G.number_of_nodes()
    A = nx.to_scipy_sparse_array(G, nodelist=range(n), format='csr')
    size, _ = lr.core(A.indptr.astype(np.int64), A.indices.astype(np.int64))
    return int(size)


def run(G, n):
    cl = [sorted(c) for c in nx.find_cliques(G) if len(c) >= 2]
    mg = [sorted(c) for c in merge_closure([frozenset(c) for c in cl])[0]]
    exact = exact_min_cover(G)
    gcore = ordinary_leaf_removal(G)
    out = {'n': n, 'min_cover': exact, 'n_cliques': len(cl), 'n_meta': len(mg),
           'graph_core': gcore / n, 'graph_solved': gcore == 0}
    for tag, fam in (('clique', cl), ('merged', mg)):
        cover, core, rest = leaf_removal(fam, n, G.edges())
        solved = not rest
        out[tag + '_core'] = len(core) / n
        out[tag + '_solved'] = bool(solved)
        out[tag + '_cover'] = len(cover) if solved else None
        if solved:
            assert is_cover(G, cover), 'leaf removal returned a non-cover'
            assert len(cover) == exact, (tag, len(cover), exact)
    return out


def validate():
    """Where it consumes the chygraph it must return the minimum cover."""
    print('validation of chygraph leaf removal:')
    cases = [('path of 6', nx.path_graph(6)),
             ('star, 5 leaves', nx.star_graph(5)),
             ('tree of triangles', nx.Graph(
                 [(0, 1), (1, 2), (2, 0), (2, 3), (3, 4), (4, 2),
                  (0, 5), (5, 6), (6, 0)])),
             ('three disjoint triangles',
              nx.disjoint_union_all([nx.complete_graph(3)] * 3)),
             ('two triangles sharing an edge',
              nx.Graph([(0, 1), (0, 2), (1, 2), (1, 3), (2, 3)])),
             ('5-cycle', nx.cycle_graph(5))]
    for name, G in cases:
        n = G.number_of_nodes()
        r = run(G, n)
        got = r['clique_cover']
        tag = (f"cover {got} = exact {r['min_cover']}" if r['clique_solved']
               else f"core {r['clique_core'] * n:.0f} of {n}, no answer")
        print(f'  {name:<31} n={n:>2}  {tag}')
    tri = run(nx.disjoint_union_all([nx.complete_graph(3)] * 3), 9)
    assert tri['clique_solved'] and tri['clique_cover'] == 6, tri
    print('  three disjoint triangles give 6 of 9, density 2/3, against the')
    print('  1/2 the hard-field map of Sec. 10.4 returns on the same object')
    cyc = run(nx.cycle_graph(5), 5)
    assert not cyc['clique_solved']
    print('  and a bare cycle is a core, as it must be\n')


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
        rows.append(dict(ensemble='hyperbolic',
                         **run(instance(r['n'], kbar[r['n']], r['seed']),
                               r['n'])))
    seen = set()
    for r in json.load(open(RESULTS / 'gbp_karrer.json'))['runs']:
        if (r['n'], r['seed']) in seen:
            continue
        seen.add((r['n'], r['seed']))
        G, _ = karrer_graph(r['n'], r['s_mean'], {3: r['triangles']}, r['seed'])
        G.remove_edges_from(nx.selfloop_edges(G))
        rows.append(dict(ensemble='karrer', **run(G, r['n'])))
    cache, seen = {}, set()
    for r in json.load(open(RESULTS / 'gbp_real.json')):
        if (r['network'], r['ego']) in seen:
            continue
        seen.add((r['network'], r['ego']))
        G = cache.setdefault(r['key'], load(r['key']))
        ego = sorted([r['ego']] + list(G[r['ego']]), key=str)
        H = nx.convert_node_labels_to_integers(G.subgraph(ego))
        rows.append(dict(ensemble='real', network=r['network'],
                         **run(H, H.number_of_nodes())))

    json.dump(rows, OUT.open('w'), indent=1)
    summarise(rows)
    print('\nwrote', OUT)


def summarise(rows):
    def tally(sel, key):
        return f"{sum(r[key] for r in sel)}/{len(sel)}"

    print(f'\n{"class":<12}{"inst":>5}{"ordinary":>11}{"clique chygraph":>18}'
          f'{"merged":>10}')
    for ens in ('hyperbolic', 'karrer', 'real'):
        sel = [r for r in rows if r['ensemble'] == ens]
        print(f'{ens:<12}{len(sel):>5}{tally(sel, "graph_solved"):>11}'
              f'{tally(sel, "clique_solved"):>18}'
              f'{tally(sel, "merged_solved"):>10}')
    print(f'{"total":<12}{len(rows):>5}{tally(rows, "graph_solved"):>11}'
          f'{tally(rows, "clique_solved"):>18}'
          f'{tally(rows, "merged_solved"):>10}')
    strict = sum(1 for r in rows if r['clique_solved'] and not r['graph_solved'])
    back = sum(1 for r in rows if r['graph_solved'] and not r['clique_solved'])
    print(f'  the chygraph rule solves {strict} that the graph rule cannot, '
          f'and {back} the other way')
    assert back == 0, 'the chygraph rule should dominate'


if __name__ == '__main__':
    main()
