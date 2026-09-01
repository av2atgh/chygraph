"""Merging against generalised belief propagation, on the same instances.

Sections 15.3 to 15.5 run GBP on maximal-clique region graphs; Sec. 16.1 runs
the merge closure on whole ensembles.  The two are never compared on the same
object, and they should be, because they are alternative repairs for one
problem.  This regenerates every instance the three GBP probes used and asks
what the merge closure would have cost there instead.

The comparison is asymmetric by construction and that is the point.  Where the
closure terminates at an enumerable size it returns the exact ln Z -- there is
no iteration, no fixed point to converge to, and no counting scheme to be
unbounded -- at a cost of 2^c on the largest meta-complex.  Where it does not,
it returns nothing at all and GBP is what remains.  So the useful number is the
distribution of c against the GBP error on the same instance.

    python probe/merge_vs_gbp.py
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

RESULTS = Path(__file__).parent / 'results'
OUT = RESULTS / 'merge_vs_gbp.json'
ENUMERABLE = 25          # 2^25 states, the limit Fig. 14.4 uses


def _largest(G):
    """Cost of each repair on one graph.

    Returns (largest meta-complex, number of cliques, largest maximal clique).
    The last is what GBP pays: its regions are the cliques and their
    intersections, so its largest table is 2^(largest clique).  Merging pays
    2^(largest meta-complex).  Comparing the two exponents is the whole of the
    cost question.
    """
    cl = [frozenset(c) for c in nx.find_cliques(G) if len(c) >= 2]
    merged, _ = merge_closure(cl)
    return max(len(x) for x in merged), len(cl), max(len(c) for c in cl)


def hrg_instances():
    from gbp_cliques import SIZES, instance
    kbar = {n: k for n, k in SIZES}
    rows = []
    for r in json.load(open(RESULTS / 'gbp_cliques.json')):
        G = instance(r['n'], kbar[r['n']], r['seed'])
        c, ncl, mx = _largest(G)
        rows.append(dict(ensemble='hyperbolic', n=r['n'], beta_J=r['beta_J'],
                         chordal=r['chordal'], gbp=r['gbp'],
                         kikuchi=r['kikuchi'], residual=r['residual'],
                         consistency=r['consistency'], largest=c, cliques=ncl, maxclique=mx))
    return rows


def real_instances():
    from gbp_real import load
    cache, rows = {}, []
    for r in json.load(open(RESULTS / 'gbp_real.json')):
        G = cache.setdefault(r['key'], load(r['key']))
        ego = sorted([r['ego']] + list(G[r['ego']]), key=str)
        H = nx.convert_node_labels_to_integers(G.subgraph(ego))
        c, ncl, mx = _largest(H)
        rows.append(dict(ensemble='real', network=r['network'],
                         family=r['family'], n=r['n'], beta_J=r['beta_J'],
                         chordal=r['chordal'], gbp=r['gbp'],
                         kikuchi=r['kikuchi'], residual=r['residual'],
                         consistency=r['consistency'], largest=c, cliques=ncl, maxclique=mx))
    return rows


def karrer_instances():
    path = RESULTS / 'gbp_karrer.json'
    if not path.exists():
        return []
    from merge import karrer_graph
    rows = []
    for r in json.load(open(path))['runs']:
        G, _ = karrer_graph(r['n'], r['s_mean'], {3: r['triangles']},
                            r['seed'])
        G.remove_edges_from(nx.selfloop_edges(G))
        c, ncl, mx = _largest(G)
        rows.append(dict(ensemble='karrer', n=r['n'], beta_J=r['beta_J'],
                         chordal=r['chordal'], gbp=r['gbp'],
                         kikuchi=r['kikuchi'], residual=r['residual'],
                         consistency=r['consistency'], largest=c, cliques=ncl, maxclique=mx))
    return rows


def sound(r):
    return (r['residual'] < 1e-9 and r['consistency'] is not None
            and r['consistency'] < 1e-6)


def summarise(rows):
    for ens in ('hyperbolic', 'karrer', 'real'):
        sel = [r for r in rows if r['ensemble'] == ens]
        if not sel:
            continue
        c = [r['largest'] for r in sel]
        ok = [r for r in sel if r['largest'] <= ENUMERABLE]
        print(f'\n{ens}: {len(sel)} runs, largest meta-complex '
              f'{min(c)} to {max(c)} atoms')
        print(f'  merge closure enumerable (c <= {ENUMERABLE}) on '
              f'{len(ok)}/{len(sel)} of them')
        bad = [r for r in sel if not r['chordal'] and not sound(r)]
        bad_ok = [r for r in bad if r['largest'] <= ENUMERABLE]
        print(f'  of the {len(bad)} runs where GBP gives no sound fixed point, '
              f'merging is exact on {len(bad_ok)}')
        both = [r for r in sel if not r['chordal'] and sound(r)
                and r['largest'] <= ENUMERABLE]
        if both:
            g = [abs(r['gbp']) for r in both]
            print(f'  on the {len(both)} where both work, GBP errs '
                  f'{min(g):.1e} to {max(g):.1e} and merging is exact')


def main():
    rows = hrg_instances() + karrer_instances() + real_instances()
    json.dump(rows, OUT.open('w'), indent=1)
    summarise(rows)
    print('\nwrote', OUT)


if __name__ == '__main__':
    main()
