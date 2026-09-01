"""Leaf-removal core of every Part IV instance.

Figures 14.2 and 16.3 plot the error in ln Z against the two structural
quantities that predict it -- the fraction of bonds lying inside two or more
complexes, and the number of meta-complexes the closure leaves.  This adds the
third quantity those figures are read against: the core-percolation core, the
fraction of nodes that leaf removal cannot strip.

Leaf removal is Bauer and Golinelli's: repeatedly delete a degree-one node
together with its neighbour, and whatever survives when no leaf remains is the
core.  A tree reduces to nothing; a cycle survives entire.  It is the structure
on which the algorithm has nothing to say, and Chapter 11 is about it.  The
routine is `leafremoval.core` from ~/av2atg/computational_complexity, the same
one Chapter 11 uses, so the two chapters count the same object.

The core fraction is a property of the graph and not of the coupling, so it is
recorded once per instance rather than once per run.

    python probe/core_fraction.py
"""

import json
import sys
from itertools import combinations
from pathlib import Path

import networkx as nx
import numpy as np

CC = Path.home() / 'av2atg' / 'computational_complexity' / 'code'
sys.path.insert(0, str(CC))
sys.path.insert(0, str(Path.home() / 'av2atg' / 'chygraph' / 'src'))
sys.path.insert(0, str(Path.home() / 'av2atg' / 'statmech'
                      / 'book' / 'figs'))
sys.path.insert(0, str(Path.home() / 'av2atg' / 'statmech' / 'probe'))

import leafremoval as lr  # noqa: E402

from merge import merge_closure  # noqa: E402

RESULTS = Path(__file__).parent / 'results'
OUT = RESULTS / 'core_fraction.json'


def graph_core(G):
    """Leaf-removal core of the graph itself, for reference."""
    n = G.number_of_nodes()
    A = nx.to_scipy_sparse_array(G, nodelist=range(n), format='csr')
    size, _ = lr.core(A.indptr.astype(np.int64), A.indices.astype(np.int64))
    return int(size), size / n


def chygraph_core(complexes, n):
    """Atoms the chygraph recursion cannot strip by propagation alone.

    Leaf removal on the incidence structure rather than on the graph, because
    that is the object the messages run on.  An atom lying in a single complex
    is determined by that complex and leaves; a complex holding one atom or none
    has nothing left to say and leaves with it.  Iterating to a fixed point is
    the 2-core of the bipartite incidence graph, and what survives is exactly
    what the propagation cannot resolve.

    It is zero precisely when the incidence structure is a forest, which is the
    condition Sec. 16.2 finds decides exactness -- so this is that criterion
    counted rather than tested.
    """
    B = nx.Graph()
    B.add_nodes_from(('v', v) for c in complexes for v in c)
    for i, c in enumerate(complexes):
        for v in c:
            B.add_edge(('c', i), ('v', v))
    core = nx.k_core(B, k=2)
    atoms = sum(1 for x in core if x[0] == 'v')
    return atoms, atoms / n


def profile(G, n):
    cl = [frozenset(c) for c in nx.find_cliques(G)]
    merged, _ = merge_closure(cl)
    cov = {}
    for c in cl:
        for e in combinations(sorted(c), 2):
            cov[e] = cov.get(e, 0) + 1
    m = G.number_of_edges()
    gsize, gfrac = graph_core(G)
    csize, cfrac = chygraph_core([sorted(c) for c in cl], n)
    msize, mfrac = chygraph_core([sorted(c) for c in merged], n)
    return {'n': n, 'n_bonds': m, 'n_meta': len(merged),
            'clustering': nx.average_clustering(G),
            'doubled_bonds': sum(1 for v in cov.values() if v > 1),
            'doubled_frac': sum(1 for v in cov.values() if v > 1) / max(m, 1),
            'graph_core': gsize, 'graph_core_frac': gfrac,
            'clique_core': csize, 'clique_core_frac': cfrac,
            'merged_core': msize, 'merged_core_frac': mfrac}


def karrer_core_examples():
    """The surviving core of two Karrer--Newman instances, for Fig. 16.5.

    The smallest core in the sample and the largest, recorded with the full
    membership of each surviving complex and with the atoms leaf removal
    stripped, so the figure can show what went and what stayed.
    """
    from merge import karrer_graph
    runs = json.load(open(RESULTS / 'gbp_karrer.json'))['runs']
    seen, found = set(), []
    for r in runs:
        k = (r['n'], r['seed'])
        if k in seen:
            continue
        seen.add(k)
        G, _ = karrer_graph(r['n'], r['s_mean'], {3: r['triangles']}, r['seed'])
        G.remove_edges_from(nx.selfloop_edges(G))
        merged, _ = merge_closure([frozenset(c) for c in nx.find_cliques(G)])
        B = nx.Graph()
        B.add_nodes_from(('v', v) for c in merged for v in c)
        for i, c in enumerate(merged):
            for v in c:
                B.add_edge(('c', i), ('v', v))
        K = nx.k_core(B, k=2)
        if K.number_of_nodes() == 0:
            continue
        atoms = sorted(x[1] for x in K if x[0] == 'v')
        cx = [sorted(merged[i]) for i in range(len(merged)) if ('c', i) in K]
        found.append({
            'n': r['n'], 'seed': r['seed'], 'atoms': atoms,
            'complexes': cx,
            'cycles': K.number_of_edges() - K.number_of_nodes()
                      + nx.number_connected_components(K),
            'components': nx.number_connected_components(K)})
    found.sort(key=lambda d: (len(d['atoms']), d['cycles']))
    return {'smallest': found[0], 'largest': found[-1]}


def main():
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
        rows.append(dict(ensemble='hyperbolic', **profile(G, r['n'])))

    seen = set()
    for r in json.load(open(RESULTS / 'gbp_karrer.json'))['runs']:
        if (r['n'], r['seed']) in seen:
            continue
        seen.add((r['n'], r['seed']))
        G, _ = karrer_graph(r['n'], r['s_mean'], {3: r['triangles']}, r['seed'])
        G.remove_edges_from(nx.selfloop_edges(G))
        rows.append(dict(ensemble='karrer', **profile(G, r['n'])))

    cache, seen = {}, set()
    for r in json.load(open(RESULTS / 'gbp_real.json')):
        if (r['network'], r['ego']) in seen:
            continue
        seen.add((r['network'], r['ego']))
        G = cache.setdefault(r['key'], load(r['key']))
        ego = sorted([r['ego']] + list(G[r['ego']]), key=str)
        H = nx.convert_node_labels_to_integers(G.subgraph(ego))
        rows.append(dict(ensemble='real', network=r['network'],
                         **profile(H, H.number_of_nodes())))

    json.dump({'instances': rows, 'karrer_examples': karrer_core_examples()},
              OUT.open('w'), indent=1)
    print(f'{"":12} {"":>3}   clique chygraph      merged chygraph')
    for ens in ('hyperbolic', 'karrer', 'real'):
        sel = [r for r in rows if r['ensemble'] == ens]
        c = [r['clique_core_frac'] for r in sel]
        m = [r['merged_core_frac'] for r in sel]
        print(f'{ens:<12} {len(sel):>3}   median {np.median(c):.2f}, '
              f'{sum(1 for x in c if x < 1e-9)} core-free    '
              f'median {np.median(m):.2f}, '
              f'{sum(1 for x in m if x < 1e-9)} core-free')
    print('wrote', OUT)


if __name__ == '__main__':
    main()
