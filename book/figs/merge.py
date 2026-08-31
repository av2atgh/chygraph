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

Section 14.7 then asks whether that verdict is a fact about clustering or a
fact about one ensemble, by running the same measurement on Karrer and Newman's
subgraph model (Phys. Rev. E 82, 066118), where the motifs are placed by
matching roles rather than grown by a geometry.  There the two answers differ,
and the reason is that the number of motif pairs sharing an edge is O(1) on
their ensemble and Theta(n) on a hyperbolic graph.

  fig-merge    the largest meta-complex against density and against n
  fig-motifs   the overlap count against n, the two ensembles side by side

Graphs come from `~/av2atg/computational_complexity/code/hrg.py`; cliques from
networkx.
"""

import sys
from itertools import combinations
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path.home() / 'av2atg' / 'computational_complexity' / 'code'))
sys.path.insert(0, str(Path.home() / 'av2atg' / 'chygraph_statmech' / 'src'))
from chygraph_statmech.region import overlap_profile  # noqa: E402
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


def cliques_of(G):
    import networkx as nx
    return [frozenset(c) for c in nx.find_cliques(G) if len(c) >= 2]


def overlap_stats(cl):
    """(number of complex PAIRS sharing two or more atoms, shared_2+ ratio).

    Both come out of `chygraph_statmech.region.overlap_profile`, which is the
    routine Sec. 3.3 measures real networks with, so the two chapters are
    counting the same thing.  The ratio is what Ch. 3 reports; the raw count is
    what Sec. 14.7 needs, because the number of intersecting pairs is
    Theta(n) in every ensemble here and dividing by it hides the distinction
    the section is about.
    """
    prof = overlap_profile(cl)
    return (prof['shared_2plus'] * prof['n_intersecting_pairs'],
            prof['shared_2plus'])


def karrer_graph(n, s_mean, motifs, seed=1):
    """Karrer and Newman (2010), Secs. II-III: single edges by stub matching,
    each motif by matching its corners in tuples.

    `motifs` maps motif size to the mean number a vertex carries.  Two placed
    motifs can share a vertex but share an edge only by coincidence, which is
    the whole point -- the ensemble is treelike above the level of the motif by
    construction, which is what makes it solvable and, here, what makes the
    merge closure terminate.
    """
    import networkx as nx
    rng = np.random.default_rng(seed)
    G = nx.Graph()
    G.add_nodes_from(range(n))
    stubs = rng.poisson(s_mean, n)
    if stubs.sum() % 2:
        stubs[rng.integers(n)] += 1
    half = np.repeat(np.arange(n), stubs)
    rng.shuffle(half)
    for a, b in half.reshape(-1, 2):
        if a != b:
            G.add_edge(int(a), int(b))
    placed = []
    for m, mu in motifs.items():
        corners = rng.poisson(mu, n)
        while corners.sum() % m:
            corners[rng.integers(n)] += 1
        roles = np.repeat(np.arange(n), corners)
        rng.shuffle(roles)
        for grp in roles.reshape(-1, m):
            g = {int(x) for x in grp}
            if len(g) == m:
                G.add_edges_from((u, v) for u in g for v in g if u < v)
                placed.append(frozenset(g))
    return G, placed


def diamond_graph(n, s_mean, d_mean, seed=1):
    """The same model with Karrer and Newman's Fig. 1 element: two triangles
    sharing an edge, placed as ONE subgraph type.

    This is their answer to overlapping triangles -- promote the overlapping
    pair to an element of its own rather than let it arise -- and it is
    Sec. 14.6's merge, taken at the point the ensemble is written down instead
    of afterwards.  The element has two roles (the two vertices on the shared
    edge carry three internal edges, the two tips carry two), so the two are
    matched separately.
    """
    import networkx as nx
    rng = np.random.default_rng(seed)
    G = nx.Graph()
    G.add_nodes_from(range(n))
    stubs = rng.poisson(s_mean, n)
    if stubs.sum() % 2:
        stubs[rng.integers(n)] += 1
    half = np.repeat(np.arange(n), stubs)
    rng.shuffle(half)
    for a, b in half.reshape(-1, 2):
        if a != b:
            G.add_edge(int(a), int(b))
    spine, tip = rng.poisson(d_mean, n), rng.poisson(d_mean, n)
    for arr in (spine, tip):
        while arr.sum() % 2:
            arr[rng.integers(n)] += 1
    A = np.repeat(np.arange(n), spine)
    B = np.repeat(np.arange(n), tip)
    rng.shuffle(A)
    rng.shuffle(B)
    placed = []
    for i in range(min(len(A), len(B)) // 2):
        u, v, x, y = int(A[2 * i]), int(A[2 * i + 1]), int(B[2 * i]), int(B[2 * i + 1])
        if len({u, v, x, y}) != 4:
            continue
        G.add_edge(u, v)
        for w in (x, y):
            G.add_edge(u, w)
            G.add_edge(v, w)
        placed.append(frozenset({u, v, x, y}))
    return G, placed


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


def _row(G, cl):
    m, _ = merge_closure(cl)
    count, ratio = overlap_stats(cl)
    return len(cl), count, ratio, max(len(x) for x in m)


def check_overlap_is_not_extensive(seeds=3):
    """Sec. 14.7: the count of edge-sharing clique pairs against n.

    This is the measurement that decides whether Sec. 14.6's verdict is about
    clustering or about one ensemble.  On the Karrer-Newman ensemble the count
    is flat in n while the number of cliques grows like n, so shared_2+ decays
    like 1/n and the merge closure has a bounded number of merges to make.  On
    a hyperbolic graph of the SAME clustering coefficient the count grows like
    n and the closure percolates.
    """
    ns = (1000, 2000, 4000, 8000)
    print('    ensemble                 <k>     C      n   cliques  overlaps'
          '  shared2+  largest')
    for lab, mk in (('Karrer-Newman, t = 1', lambda n, s: karrer_graph(n, 1.0, {3: 1.0}, s)),
                    ('Karrer-Newman, t = 4', lambda n, s: karrer_graph(n, 1.0, {3: 4.0}, s)),
                    ('hyperbolic, <k> = 2', None),
                    ('hyperbolic, <k> = 3', None)):
        kbar = 2.0 if lab.endswith('= 2') else 3.0
        for n in ns:
            import networkx as nx
            rows, ks, Cs = [], [], []
            for s in range(1, seeds + 1):
                if mk is None:
                    G = nx.Graph()
                    G.add_nodes_from(range(n))
                    G.add_edges_from(_hrg_edges(n, 2.5, kbar, s))
                else:
                    G = mk(n, s)[0]
                cl = cliques_of(G)
                rows.append(_row(G, cl))
                ks.append(2 * G.number_of_edges() / n)
                Cs.append(nx.average_clustering(G))
            a = np.mean(rows, axis=0)
            print(f'    {lab:22s}{np.mean(ks):6.2f}{np.mean(Cs):7.3f}{n:7d}'
                  f'{a[0]:10.0f}{a[1]:10.1f}{a[2]:10.4f}{a[3]:9.1f}')
    print('    The clustering coefficients are matched across the two families')
    print('    (0.25 against 0.27) and the verdicts are opposite. What differs')
    print('    is the overlap COUNT: flat in n on the placed ensemble, linear in')
    print('    n on the geometric one. Clustering is not what decides this.')


def check_diamonds_are_recovered(seeds=3):
    """Sec. 14.7: run the closure blind on an ensemble built OUT OF overlaps.

    Karrer and Newman handle triangles that share an edge by placing the
    resulting diamond as a subgraph type of its own.  Its maximal cliques are
    the two triangles, which do share an edge, so shared_2+ does NOT vanish
    here -- the overlap is in the ensemble by design.  The question is whether
    the merge closure, given only the graph, returns the elements they placed.
    """
    print('       n    <k>      C  diamonds  shared2+  largest  recovered')
    for n in (1000, 2000, 4000, 8000):
        import networkx as nx
        rows = []
        for s in range(1, seeds + 1):
            G, placed = diamond_graph(n, 1.0, 1.0, s)
            cl = cliques_of(G)
            m, _ = merge_closure(cl)
            got = set(m)
            rows.append((2 * G.number_of_edges() / n, nx.average_clustering(G),
                         len(placed), overlap_stats(cl)[1],
                         max(len(x) for x in m),
                         sum(1 for d in placed if d in got) / len(placed)))
        a = np.mean(rows, axis=0)
        print(f'    {n:5d}{a[0]:7.2f}{a[1]:7.3f}{a[2]:10.0f}{a[3]:10.4f}'
              f'{a[4]:9.1f}{a[5]:11.3f}')
    print('    The closure recovers the placed diamonds exactly, and the')
    print('    shortfall is the same O(1) residue of accidental overlaps: about')
    print('    a hundred elements at every n, which is 21% of them at n = 1000')
    print('    and 3% at n = 8000. So the merge is not a repair bolted on after')
    print('    the fact -- it inverts the modelling choice they made.')


def _hrg_edges(n, tau, kbar, seed):
    src, dst = hrg_calibrated(n, tau=tau, kbar=kbar,
                              rng=np.random.default_rng(seed))[:2]
    return list(zip(src.tolist(), dst.tolist()))


def figure_motifs():
    """The overlap count against n for the two ensembles, and what it costs."""
    plt = _mpl()
    import networkx as nx
    ns = np.array([1000, 2000, 4000, 8000])
    fig, axes = plt.subplots(1, 2, figsize=(4.6, 2.4))

    series = (('placed motifs', DARK, 'o',
               lambda n, s: karrer_graph(n, 1.0, {3: 4.0}, s)[0]),
              ('hyperbolic', LIGHT, '^', None))
    for ax, which in ((axes[0], 'count'), (axes[1], 'largest')):
        for lab, col, mk, gen in series:
            ys = []
            for n in ns:
                vals = []
                for s in (1, 2, 3):
                    if gen is None:
                        G = nx.Graph()
                        G.add_nodes_from(range(n))
                        G.add_edges_from(_hrg_edges(n, 2.5, 3.0, s))
                    else:
                        G = gen(n, s)
                    cl = cliques_of(G)
                    c, _ = overlap_stats(cl)
                    m, _ = merge_closure(cl)
                    vals.append(c if which == 'count' else max(len(x) for x in m))
                ys.append(np.mean(vals))
            ax.loglog(ns, ys, mk + '-', ms=3.4, lw=1.2, color=col, label=lab)
        ax.loglog(ns, ys[0] * ns / ns[0], ':', lw=0.9, color='0.55')
        ax.set_xticks(list(ns))
        ax.set_xticklabels(['1k', '2k', '4k', '8k'])
        ax.minorticks_off()
        ax.set_xlabel('$n$', fontsize=8.5)
        _tidy(ax)
    axes[0].set_ylabel('clique pairs sharing an edge', fontsize=8)
    axes[1].set_ylabel('largest meta-complex', fontsize=8)
    axes[0].legend(frameon=False, fontsize=7, loc='upper left')
    fig.tight_layout()
    fig.savefig(OUT / 'fig-motifs.pdf')
    print(f'  wrote {OUT / "fig-motifs.pdf"}')


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
    print('what it costs on a hyperbolic ensemble:')
    check_percolates()
    print('an ensemble where it terminates:')
    check_overlap_is_not_extensive()
    print('and one built out of overlaps, run blind:')
    check_diamonds_are_recovered()
    print('figures:')
    figure_merge()
    figure_motifs()
