"""Chapter 14: merging overlapping complexes, and when that is affordable.

Two complexes that share more than one atom break the treelike assumption.
The obvious repair is to merge them into one meta-complex and sum its interior
exactly -- which is legitimate, because a complex whose vertices are complexes
is what a chygraph is.  The cost is 2^c on the merged complex, so the question
is whether the merging terminates at an enumerable size.

It is a percolation problem.  Put a node on every complex and a link between two
complexes that share two or more atoms; the meta-complexes are the connected
components of that graph, and the repair is affordable exactly while it has no
giant one.  This script measures where that giant appears on the ensemble the
book cares about.

  fig-merge   the largest meta-complex against density and against n

Graphs come from `~/av2atg/computational_complexity/code/hrg.py`; cliques from
networkx.
"""

import sys
from itertools import combinations
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path.home() / 'av2atg' / 'computational_complexity' / 'code'))
from hrg import hrg_calibrated  # noqa: E402

OUT = Path(__file__).resolve().parent
DARK, MID, LIGHT = '0.10', '0.45', '0.70'


def _mpl():
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    return plt


def _tidy(ax):
    ax.tick_params(labelsize=8)
    for sp in ('top', 'right'):
        ax.spines[sp].set_visible(False)


def merge_closure(complexes, max_rounds=40):
    """Merge complexes sharing two or more atoms, iterated to a fixed point.

    One pass is not enough in principle: two meta-complexes built from
    components that were themselves disjoint can still share two atoms between
    their unions.  The loop runs until a pass changes nothing.
    """
    cx = [frozenset(c) for c in complexes]
    for rnd in range(max_rounds):
        by_pair = {}
        for i, a in enumerate(cx):
            for p in combinations(sorted(a), 2):
                by_pair.setdefault(p, []).append(i)
        parent = list(range(len(cx)))

        def find(x):
            while parent[x] != x:
                parent[x] = parent[parent[x]]
                x = parent[x]
            return x

        merged = False
        for ids in by_pair.values():
            if len(ids) > 1:
                merged = True
                r0 = find(ids[0])
                for j in ids[1:]:
                    r = find(j)
                    if r != r0:
                        parent[r] = r0
        if not merged:
            return cx, rnd
        groups = {}
        for i in range(len(cx)):
            groups.setdefault(find(i), set()).update(cx[i])
        cx = [frozenset(v) for v in groups.values()]
    raise RuntimeError('merge closure did not converge')


def cliques_of_hrg(n, tau, kbar, seed):
    import networkx as nx
    src, dst = hrg_calibrated(n, tau=tau, kbar=kbar,
                              rng=np.random.default_rng(seed))[:2]
    G = nx.Graph()
    G.add_nodes_from(range(n))
    G.add_edges_from(zip(src.tolist(), dst.tolist()))
    return [frozenset(c) for c in nx.find_cliques(G) if len(c) >= 2]


def largest(n, tau, kbar, seeds=3):
    out = []
    for s in range(1, seeds + 1):
        cl = cliques_of_hrg(n, tau, kbar, s)
        m, _ = merge_closure(cl)
        out.append(max(len(x) for x in m))
    return float(np.mean(out))


def check_two_triangles():
    """The repair is exact on Sec. 14.4's example, trivially."""
    m, rnd = merge_closure([(0, 1, 2), (1, 2, 3)])
    assert len(m) == 1 and m[0] == frozenset({0, 1, 2, 3}), m
    print(f'  two triangles sharing an edge merge into one complex of '
          f'{len(m[0])} atoms')
    print(f'    its interior is 2^{len(m[0])} = {2**len(m[0])} states, summed '
          f'exactly by Eq. (8.4);')
    print('    the overlap is gone because there is nothing left to overlap '
          'with')


def check_percolates():
    """On hyperbolic random graphs the merging runs away, and early."""
    print('  largest meta-complex, tau = 2.5, mean over 3 seeds')
    print('     <k>     n=1000   n=2000   n=4000    growth   verdict')
    for kbar in (0.3, 0.5, 0.8, 1.5, 3.0):
        r = [largest(n, 2.5, kbar) for n in (1000, 2000, 4000)]
        g = r[-1] / max(r[0], 1e-9)
        v = 'bounded' if g < 2.0 else 'grows with n'
        print(f'  {kbar:>6}   {r[0]:>8.1f} {r[1]:>8.1f} {r[2]:>8.1f}   '
              f'{g:>7.1f}   {v}')
    print('    a meta-complex that grows with n cannot be enumerated: 2^c with')
    print('    c extensive is the original problem back again')


def figure_merge():
    plt = _mpl()
    fig, axes = plt.subplots(1, 2, figsize=(4.6, 2.5))

    ax = axes[0]
    ks = (0.3, 0.5, 0.8, 1.2, 2.0, 3.0, 4.0)
    for n, col, mk in ((1000, LIGHT, '^'), (2000, MID, 's'), (4000, DARK, 'o')):
        ax.semilogy(ks, [largest(n, 2.5, k) for k in ks], mk + '-', ms=3.4,
                    lw=1.2, color=col, label=f'$n={n}$')
    ax.axhline(25, ls=':', lw=0.9, color='0.5')
    ax.set_xlabel(r'mean degree $\langle k\rangle$', fontsize=8.5)
    ax.set_ylabel('largest meta-complex', fontsize=8.5)
    ax.legend(frameon=False, fontsize=7, loc='lower right')
    _tidy(ax)

    ax = axes[1]
    ns = (1000, 2000, 4000, 8000)
    for kbar, col, mk in ((0.3, LIGHT, '^'), (0.8, MID, 's'), (2.0, DARK, 'o')):
        ax.loglog(ns, [largest(n, 2.5, kbar, seeds=2) for n in ns], mk + '-',
                  ms=3.4, lw=1.2, color=col, label=rf'$\langle k\rangle={kbar}$')
    ax.loglog(ns, [0.05 * n for n in ns], ':', lw=0.9, color='0.55')
    ax.annotate('slope 1', xy=(2600, 90), fontsize=6.4, color='0.45')
    ax.set_xticks(list(ns))
    ax.set_xticklabels(['1k', '2k', '4k', '8k'])
    ax.minorticks_off()
    ax.set_xlabel(r'$n$', fontsize=8.5)
    ax.legend(frameon=False, fontsize=7, loc='upper left')
    _tidy(ax)

    fig.tight_layout()
    fig.savefig(OUT / 'fig-merge.pdf')
    print(f'  wrote {OUT / "fig-merge.pdf"}')


if __name__ == '__main__':
    print('the repair, on the two-triangle example:')
    check_two_triangles()
    print('what it costs on a real ensemble:')
    check_percolates()
    print('figure:')
    figure_merge()
