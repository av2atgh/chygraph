"""Error in ln Z of the chygraph built from meta-complexes.

Chapter 16 prices the merge closure and says the answer is exact where the
closure terminates.  That claim deserves testing rather than asserting, because
merging removes *pairwise* violations and not every violation: a family in which
no two complexes share more than one atom can still close a cycle through three
or more of them, complex-node-complex-node-complex, and the chygraph recursion
is a Bethe calculation on that incidence structure.  Where the cycle survives,
so does an error.

For every instance the three GBP probes used, this merges the maximal cliques to
a fixed point and then solves the merged chygraph two ways -- the Bethe counting
on isolated regions, which is what Chapter 14 called the cavity answer, and the
messages run to a fixed point on the merged family -- against the exact ln Z.

    python probe/merge_lnz.py
"""

import json
import sys
from pathlib import Path

import networkx as nx
import numpy as np

sys.path.insert(0, str(Path.home() / 'av2atg' / 'chygraph' / 'src'))
sys.path.insert(0, str(Path.home() / 'av2atg' / 'chygraph_statmech'
                      / 'book' / 'figs'))
sys.path.insert(0, str(Path.home() / 'av2atg' / 'chygraph_statmech' / 'probe'))

from merge import merge_closure  # noqa: E402

from chygraph_statmech.gbp import (GBP, exact_log_Z, ising_factors,  # noqa: E402
                                   static_log_Z)
from chygraph_statmech.region import RegionGraph, overlap_profile  # noqa: E402

RESULTS = Path(__file__).parent / 'results'
OUT = RESULTS / 'merge_lnz.json'
COUPLINGS = (0.3, 0.8)


def incidence_acyclic(complexes):
    """Is the bipartite complex-node incidence structure a forest?

    This is the condition the chygraph recursion actually needs.  Pairwise
    intersections of at most one node are necessary for it and not sufficient.
    """
    B = nx.Graph()
    for i, c in enumerate(complexes):
        for v in c:
            B.add_edge(('c', i), ('v', v))
    return nx.is_forest(B)


def solve(G, n, bJ):
    # every maximal clique, singletons included: an isolated vertex is a
    # complex of one atom and still contributes ln 2 to Z.  Filtering it out
    # makes the chygraph answer short by ln 2 per isolated vertex, which looks
    # exactly like a real error and is not one.
    cl = [frozenset(c) for c in nx.find_cliques(G)]
    merged, rounds = merge_closure(cl)
    prof = overlap_profile([sorted(c) for c in merged])
    edges = sorted({tuple(sorted(e)) for e in G.edges()})
    f = ising_factors(edges, bJ)
    exact = exact_log_Z(f, range(n))
    rg = RegionGraph([sorted(c) for c in merged], max_rounds=6)
    out = {
        'n': n, 'beta_J': bJ, 'largest': max(len(x) for x in merged),
        'n_meta': len(merged), 'rounds': rounds,
        'treelike': bool(prof['treelike']),
        'acyclic': bool(incidence_acyclic(merged)),
        'exact': exact,
    }
    # the static chygraph answer on the merged family
    out['bethe'] = static_log_Z(rg.bethe_counting(), f) - exact
    # and with the messages run on it
    if rg.counting_is_valid():
        g = GBP(rg, f, damping=0.5).run(8000)
        out['gbp'] = g.log_Z() - exact
        out['residual'] = g.residual
    else:
        out['gbp'] = None
        out['residual'] = None
    return out


def hrg():
    from gbp_cliques import SIZES, instance
    kbar = {n: k for n, k in SIZES}
    seen, rows = set(), []
    for r in json.load(open(RESULTS / 'gbp_cliques.json')):
        if (r['n'], r['seed']) in seen:
            continue
        seen.add((r['n'], r['seed']))
        G = instance(r['n'], kbar[r['n']], r['seed'])
        for bJ in COUPLINGS:
            rows.append(dict(ensemble='hyperbolic', **solve(G, r['n'], bJ)))
    return rows


def karrer():
    from merge import karrer_graph
    p = RESULTS / 'gbp_karrer.json'
    seen, rows = set(), []
    for r in json.load(open(p))['runs']:
        if (r['n'], r['seed']) in seen:
            continue
        seen.add((r['n'], r['seed']))
        G, _ = karrer_graph(r['n'], r['s_mean'], {3: r['triangles']}, r['seed'])
        G.remove_edges_from(nx.selfloop_edges(G))
        for bJ in COUPLINGS:
            rows.append(dict(ensemble='karrer', **solve(G, r['n'], bJ)))
    return rows


def real():
    from gbp_real import load
    cache, seen, rows = {}, set(), []
    for r in json.load(open(RESULTS / 'gbp_real.json')):
        if (r['network'], r['ego']) in seen:
            continue
        seen.add((r['network'], r['ego']))
        G = cache.setdefault(r['key'], load(r['key']))
        ego = sorted([r['ego']] + list(G[r['ego']]), key=str)
        H = nx.convert_node_labels_to_integers(G.subgraph(ego))
        for bJ in COUPLINGS:
            rows.append(dict(ensemble='real', network=r['network'],
                             family=r['family'],
                             **solve(H, H.number_of_nodes(), bJ)))
    return rows


def summarise(rows):
    for ens in ('hyperbolic', 'karrer', 'real'):
        sel = [r for r in rows if r['ensemble'] == ens]
        if not sel:
            continue
        tl = sum(r['treelike'] for r in sel)
        ac = sum(r['acyclic'] for r in sel)
        b = [abs(r['bethe']) for r in sel]
        g = [abs(r['gbp']) for r in sel if r['gbp'] is not None]
        nz = [r for r in sel if abs(r['bethe']) > 1e-9]
        print(f'\n{ens}: {len(sel)} runs')
        print(f'  merged family treelike in {tl}/{len(sel)}, '
              f'incidence acyclic in {ac}/{len(sel)}')
        print(f'  |error in ln Z|, static chygraph: {min(b):.1e} to {max(b):.1e}'
              f'  (median {np.median(b):.1e})')
        if g:
            print(f'  |error in ln Z|, with messages: {min(g):.1e} to '
                  f'{max(g):.1e}  (median {np.median(g):.1e})')
        print(f'  runs where the merged chygraph is NOT exact: '
              f'{len(nz)}/{len(sel)}')
        if nz:
            w = max(nz, key=lambda r: abs(r['bethe']))
            print(f'    worst: n={w["n"]} bJ={w["beta_J"]} c={w["largest"]} '
                  f'meta={w["n_meta"]} acyclic={w["acyclic"]} '
                  f'error={w["bethe"]:+.3f}')


def main():
    rows = hrg() + karrer() + real()
    json.dump(rows, OUT.open('w'), indent=1)
    summarise(rows)
    print('\nwrote', OUT)


if __name__ == '__main__':
    main()
