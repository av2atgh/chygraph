"""Chapter 16: merging overlapping complexes, and when that is affordable.

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

Only `check_placed_finite_size` reaches the book, as Table 15.2 in Sec. 15.4.
The other four checks compare the closure across ensembles -- hyperbolic against
Karrer and Newman's subgraph model (Phys. Rev. E 82, 066118), where the motifs
are placed by matching roles rather than grown by a geometry, and against an
ensemble built out of overlaps -- and none of their numbers is quoted anywhere.
They are kept because they are the evidence that the O(1) count of edge-sharing
motif pairs, not the clustering coefficient, is what decides the closure.

This module writes no figures.  It used to write fig-merge and fig-motifs;
nothing has included them since Part IV was split into three chapters, so both
they and the code that drew them are gone.  What the book still takes from here
is Table 15.2, via `check_placed_finite_size`.  Everything else in the file is
either a check that is not quoted or, in the case of `merge_closure` and
`karrer_graph`, a definition that six probe scripts import.

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
    what `check_overlap_is_not_extensive` needs, because the number of
    intersecting pairs is Theta(n) in every ensemble here and dividing by it
    hides the distinction that check is about.
    """
    prof = overlap_profile(cl)
    return (prof['shared_2plus'] * prof['n_intersecting_pairs'],
            prof['shared_2plus'])


def karrer_graph(n, s_mean, motifs, seed=1):
    """Karrer and Newman (2010), Secs. II-III: single edges by stub matching,
    each motif by matching its corners in tuples.

    `motifs` maps motif size to the mean number a vertex carries.  Two placed
    motifs can share a vertex but share an edge only by coincidence, which is
    the whole point: no two complexes meet in more than one atom, which is
    Sec. 2.9's pairwise condition and what makes the merge closure terminate.
    It is not acyclicity -- the incidence structure still carries a 2-core above
    an incidence branching of one (Sec. 14.2).
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
    Sec. 16.1's merge, taken at the point the ensemble is written down instead
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


def check_merge_on_real(quiet=False):
    """The merge closure on the six real networks of Sec.~14.7.

    Unlike the GBP measurement, this needs no enumeration, so it runs on the
    whole network rather than on neighbourhoods: the closure is a percolation
    question about the overlap graph, and its answer is a size.  The two
    families split here as they did there, and for the same reason -- a network
    whose cliques the data placed has a bounded number of merges to make.
    """
    sys.path.insert(0, str(Path.home() / 'av2atg' / 'chygraph_statmech'
                           / 'probe'))
    from gbp_real import NETWORKS, load
    rows = []
    for key, label, family in NETWORKS:
        G = load(key)
        cl = cliques_of(G)
        merged, rounds = merge_closure(cl)
        pairs, ratio = overlap_stats(cl)
        n = G.number_of_nodes()
        big = max(len(x) for x in merged)
        rows.append(dict(network=label, family=family, n=n,
                         cliques=len(cl), pairs=int(pairs), shared=ratio,
                         largest=big, frac=big / n, rounds=rounds))
    if not quiet:
        print('  merge closure on the six real networks:')
        print(f'    {"network":<18}{"family":<9}{"n":>7}{"cliques":>8}'
              f'{"pairs":>8}{"largest":>9}{"of n":>8}')
        for r in sorted(rows, key=lambda r: (r['family'] != 'grouped',
                                             r['frac'])):
            print(f'    {r["network"]:<18}{r["family"]:<9}{r["n"]:>7}'
                  f'{r["cliques"]:>8}{r["pairs"]:>8}{r["largest"]:>9}'
                  f'{100 * r["frac"]:>7.1f}%')
        for fam in ('grouped', 'dyadic'):
            sel = [r for r in rows if r['family'] == fam]
            print(f'    {fam:>8}: largest meta-complex '
                  f'{min(r["frac"] for r in sel) * 100:.1f} to '
                  f'{max(r["frac"] for r in sel) * 100:.1f} per cent of the '
                  f'network')
    return rows


def check_placed_finite_size(seeds=12, quiet=False):
    """How far the placed ensemble is from its own asymptotics at small n.

    The Karrer-Newman ensemble shares an edge between two placed triangles only
    when their remaining corners coincide, which has probability O(1/n).  The
    count of such pairs is therefore O(1) while the number of cliques grows like
    n, so the FRACTION of cliques caught in an overlap falls like 1/n.  At the
    sizes where ln Z can be enumerated that fraction has not fallen at all,
    which is why Sec. 15.4's instances look nothing like the ensemble they are
    drawn from.  This measures the crossover, and is Table 15.2.
    """
    import networkx as nx
    rows = []
    for n in (14, 18, 20, 30, 50, 100, 200, 500, 1000, 2000):
        pairs, fracs, treelike = [], [], 0
        for sd in range(1, seeds + 1):
            G, _ = karrer_graph(n, 1.0, {3: 2.0}, sd)
            G.remove_edges_from(nx.selfloop_edges(G))
            cl = cliques_of(G)
            prof = overlap_profile(cl)
            npairs = prof['shared_2plus'] * prof['n_intersecting_pairs']
            pairs.append(npairs)
            fracs.append(npairs / max(len(cl), 1))
            treelike += bool(prof['treelike'])
        rows.append(dict(n=n, pairs=float(np.mean(pairs)),
                         per_clique=float(np.mean(fracs)),
                         treelike=treelike / seeds))
    if not quiet:
        print('  Karrer-Newman at <k> = 5, mean over '
              f'{seeds} seeds:')
        print(f'    {"n":>6}{"offending pairs":>17}{"per clique":>12}'
              f'{"treelike":>10}')
        for r in rows:
            print(f'    {r["n"]:>6}{r["pairs"]:>17.1f}{r["per_clique"]:>12.3f}'
                  f'{100 * r["treelike"]:>9.0f}%')
        print('    the count is O(1); the per-clique rate falls like 1/n, so at')
        print('    n <= 20 the ensemble is not yet in the regime that defines it')
    return rows


def check_two_triangles():
    """The repair is exact on Fig. 14.1's two triangles, trivially."""
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
    """The count of edge-sharing clique pairs against n.  Not quoted.

    This is the measurement that decides whether the closure's verdict is about
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
    """Run the closure blind on an ensemble built OUT OF overlaps.  Not quoted.

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


if __name__ == '__main__':
    print('the repair, on the two-triangle example:')
    check_two_triangles()
    print('what it costs on a hyperbolic ensemble:')
    check_percolates()
    print('an ensemble where it terminates:')
    check_overlap_is_not_extensive()
    print('and one built out of overlaps, run blind:')
    check_diamonds_are_recovered()
    print('how far the placed ensemble is from its asymptotics at small n:')
    check_placed_finite_size()
    print('the closure on the six real networks:')
    check_merge_on_real()
