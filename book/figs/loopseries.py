"""Chapter 24: the loop series on the pairwise and the promoted factor graph.

One figure, two tables, and the numbers Sec. 24.3 quotes.

  fig-loopseries       (a) partial sums of the pairwise loop series on a ring
                       of six triangles, by loop size, against the whole; the
                       promoted series is one term. (b) the error of belief
                       propagation on the pairwise factor graph against the
                       error on the promoted one, forty small clustered
                       instances at two couplings.
  tab-loopseries-ring  rings of k triangles: loops and errors, both graphs.
  tab-loopseries-dc    the chygraph recursion's error on two triangles
                       sharing an edge, split into its double-count part and
                       its loop part, at Chapter 14's couplings; and the same
                       split averaged over the hyperbolic instances of
                       Chapter 14 at n = 14.

Checks run first: the single-loop term is tanh^3 on a triangle; the identity
Z = Z_BP (1 + sum r_C) holds to 1e-12 on two triangles (three factor
graphs), on rings of four to six triangles (both graphs) and on a random
graph; and the double-count factor graph reproduces the recursion of
`statmech/probe/cavity_clique.py` to 1e-10.
"""

import json
import math
import sys
import time
from pathlib import Path

import networkx as nx
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'statmech' / 'src'))
sys.path.insert(0, str(ROOT / 'percolation' / 'src'))
sys.path.append(str(ROOT / 'statmech' / 'probe'))
from statmech.loopseries import (  # noqa: E402
    BinaryFactorGraph, loop_series, partial_sums,
)

OUT = Path(__file__).resolve().parent
DARK, MID, LIGHT = '0.10', '0.45', '0.70'
TWO_TRIANGLES = [(0, 1, 2), (1, 2, 3)]
TWO_EDGES = [(0, 1), (0, 2), (1, 2), (1, 3), (2, 3)]


def _mpl():
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    return plt


def _tidy(ax):
    ax.tick_params(labelsize=8)
    for sp in ('top', 'right'):
        ax.spines[sp].set_visible(False)


def ring_of_triangles(k):
    G = nx.Graph()
    hub = [3 * i for i in range(k)]
    for i in range(k):
        a, b, c = hub[i], 3 * i + 1, hub[(i + 1) % k]
        G.add_edges_from([(a, b), (b, c), (a, c)])
    return nx.convert_node_labels_to_integers(G)


def cliques_and_edges(G):
    cx = [tuple(sorted(c)) for c in nx.find_cliques(G)]
    e = sorted(tuple(sorted(x)) for x in G.edges())
    return cx, e


def series(fg):
    """BP, then the series; returns (exact - Bethe, terms)."""
    fg.bp(damping=0.5)
    assert fg.converged(), fg.residual
    return fg.exact_log_Z() - fg.log_Z_bethe(), loop_series(fg)


# ------------------------------------------------------------------ checks

def check_identity():
    fg = BinaryFactorGraph.pairwise([(0, 1), (1, 2), (0, 2)], 0.5)
    err, terms = series(fg)
    assert len(terms) == 1 and abs(terms[0][2] - math.tanh(0.5) ** 3) < 1e-14
    print(f'  triangle: r_C = tanh^3 = {terms[0][2]:.6f}')
    worst = 0.0
    cases = []
    for bJ in (0.5, 1.0):
        cases += [(f'two triangles, promoted, bJ={bJ}',
                   BinaryFactorGraph.promoted(TWO_TRIANGLES, TWO_EDGES, bJ)),
                  (f'two triangles, pairwise, bJ={bJ}',
                   BinaryFactorGraph.pairwise(TWO_EDGES, bJ)),
                  (f'two triangles, double count, bJ={bJ}',
                   BinaryFactorGraph.promoted(TWO_TRIANGLES, TWO_EDGES, bJ, 'all'))]
    for k in (4, 5, 6):
        cx, e = cliques_and_edges(ring_of_triangles(k))
        cases += [(f'ring k={k}, promoted', BinaryFactorGraph.promoted(cx, e, 0.8)),
                  (f'ring k={k}, pairwise', BinaryFactorGraph.pairwise(e, 0.8))]
    G = nx.gnm_random_graph(9, 14, seed=3)
    cases.append(('random n=9 m=14', BinaryFactorGraph.pairwise(
        sorted(tuple(sorted(x)) for x in G.edges()), 0.4)))
    for name, fg in cases:
        err, terms = series(fg)
        d = abs(err - math.log1p(sum(r for _, _, r in terms)))
        worst = max(worst, d)
        print(f'  {name:40s} loops {len(terms):5d}  identity to {d:.1e}')
    assert worst < 1e-12
    print(f'  Z = Z_BP (1 + sum r_C) holds to {worst:.1e} on every case')


def check_double_count_is_the_recursion():
    from cavity_clique import ChygraphBP
    for k in (4, 6):
        cx, e = cliques_and_edges(ring_of_triangles(k))
        for bJ in (0.3, 0.8):
            a = BinaryFactorGraph.promoted(cx, e, bJ, 'all').bp().log_Z_bethe()
            b = ChygraphBP(cx, bJ, damping=0.5, edges=e).run().log_Z()
            assert abs(a - b) < 1e-10, (k, bJ, a, b)
    for bJ in (0.2, 0.5, 1.0, 2.0):
        a = BinaryFactorGraph.promoted(TWO_TRIANGLES, TWO_EDGES, bJ, 'all').bp().log_Z_bethe()
        b = ChygraphBP(TWO_TRIANGLES, bJ, damping=0.5, edges=TWO_EDGES).run().log_Z()
        assert abs(a - b) < 1e-10, (bJ, a, b)
    print('  the double-count factor graph is the recursion of Chapter 14, to 1e-10')


# ------------------------------------------------------------- (1) rings

def table_rings():
    rows = []
    for k in (4, 5, 6):
        cx, e = cliques_and_edges(ring_of_triangles(k))
        for bJ in (0.3, 0.8):
            ep, tp = series(BinaryFactorGraph.promoted(cx, e, bJ))
            eq, tq = series(BinaryFactorGraph.pairwise(e, bJ))
            ps = partial_sums(tq)
            total = ps[-1][1]
            tri = ps[0][1] / total
            within = next(s for s, v in ps if abs(v - total) < 0.01 * total)
            rows.append((k, 2 * k, bJ, len(tq), eq, tri, within, ep))
            print(f'  k={k} bJ={bJ}: pairwise {len(tq):5d} loops, error {eq:.4f}, '
                  f'triangles {tri:.2f} of the sum, within 1% at size {within}; '
                  f'promoted 1 loop, error {ep:.4f}')
    with open(OUT / 'tab-loopseries-ring.tex', 'w') as f:
        f.write('% generated by figs/loopseries.py; do not edit\n')
        f.write('\\begin{tabular}{rrrrrrrr}\n\\hline\\hline\n')
        f.write('$k$ & $n$ & $\\beta J$ & \\multicolumn{4}{c}{pairwise} & promoted\\\\\n')
        f.write(' & & & loops & error & triangles & size at $1\\%$ & error\\\\\n\\hline\n')
        for k, n, bJ, nl, eq, tri, within, ep in rows:
            f.write(f'{k} & {n} & {bJ} & {nl} & {eq:.4f} & {tri:.2f} & {within} & {ep:.4f}\\\\\n')
        f.write('\\hline\\hline\n\\end{tabular}\n')
    print('  wrote tab-loopseries-ring.tex')
    return rows


# ---------------------------------------------------- (2) the decomposition

def decomposition(cx, e, bJ):
    """(recursion error, double-count part, loop part, assigned-once error)
    on one instance:  ln Z_BP^dc - ln Z = (ln Z^dc - ln Z) + (ln Z_BP^dc - ln Z^dc),
    and ln Z_BP^once - ln Z for the same factor graph with every bond in one
    clique, which is belief propagation on alpha + I for the true model."""
    dc = BinaryFactorGraph.promoted(cx, e, bJ, 'all').bp()
    assert dc.converged(), dc.residual
    once = BinaryFactorGraph.promoted(cx, e, bJ, 'once').bp()
    assert once.converged(), once.residual
    z = BinaryFactorGraph.pairwise(e, bJ).exact_log_Z()
    zdc = dc.exact_log_Z()
    zb = dc.log_Z_bethe()
    return zb - z, zdc - z, zb - zdc, once.log_Z_bethe() - z


def hyperbolic_graph(n, kbar, tau, rng, R=None):
    """A hyperbolic random graph in the disc model of Krioukov et al.: radii
    with density proportional to sinh(alpha r) on [0, R], alpha = (tau-1)/2,
    angles uniform, an edge where the hyperbolic distance is below R.  The
    generator behind Chapter 14's instances is no longer on disk, so these
    are fresh draws of the same ensemble at the same n, kbar and tau; the
    radius R is set by bisection on the mean degree averaged over draws."""
    alpha = (tau - 1) / 2

    def draw(R, rng):
        u = rng.random(n)
        r = np.arccosh(1 + u * (np.cosh(alpha * R) - 1)) / alpha
        th = rng.random(n) * 2 * np.pi
        dth = np.pi - np.abs(np.pi - np.abs(th[:, None] - th[None, :]))
        d = np.arccosh(np.maximum(1.0, np.cosh(r[:, None]) * np.cosh(r[None, :])
                                  - np.sinh(r[:, None]) * np.sinh(r[None, :]) * np.cos(dth)))
        A = (d < R) & ~np.eye(n, dtype=bool)
        return A

    if R is None:
        lo, hi = 0.5, 20.0
        probe = np.random.default_rng(12345)
        for _ in range(40):
            R = (lo + hi) / 2
            kk = np.mean([draw(R, probe).sum() / n for _ in range(200)])
            # a larger disc is a sparser graph, so too few edges means a smaller R
            lo, hi = (lo, R) if kk < kbar else (R, hi)
    A = draw(R, rng)
    G = nx.from_numpy_array(A.astype(int))
    return G, R


def hrg_instances(n=14, kbar=3.5, tau=2.5, seeds=range(1, 41), want=20):
    from statmech.region import overlap_profile
    out = []
    R = None
    for seed in seeds:
        G, R = hyperbolic_graph(n, kbar, tau, np.random.default_rng(seed), R)
        G.remove_nodes_from([v for v in list(G) if G.degree(v) == 0])
        G = nx.convert_node_labels_to_integers(G)
        if G.number_of_nodes() < 8:
            continue
        cx, e = cliques_and_edges(G)
        if overlap_profile([list(c) for c in cx])['treelike']:
            continue
        out.append((seed, cx, e))
        if len(out) == want:
            break
    return out


def table_decomposition():
    rows = []
    print('  two triangles sharing an edge: recursion error = double count + loop')
    for bJ in (0.2, 0.5, 1.0, 2.0):
        tot, dcp, lp, once = decomposition(TWO_TRIANGLES, TWO_EDGES, bJ)
        rows.append(('two triangles', bJ, tot, dcp, lp, once))
        print(f'    bJ={bJ}: {tot:+.4f} = {dcp:+.4f} {lp:+.4f};  assigned once {once:+.4f}')
    inst = hrg_instances()
    print(f'  hyperbolic instances at n = 14, {len(inst)} with overlap '
          f'(mean bonds {np.mean([len(e) for _, _, e in inst]):.1f}, '
          f'cliques {np.mean([len(cx) for _, cx, _ in inst]):.1f}):')
    for bJ in (0.3, 0.8):
        parts = np.array([decomposition(cx, e, bJ)[:] for _, cx, e in inst])
        m = parts.mean(axis=0)
        cancel = np.mean(np.sign(parts[:, 1]) != np.sign(parts[:, 2]))
        rows.append((f'hyperbolic $n=14$, mean of {len(inst)}', bJ, *m))
        print(f'    bJ={bJ}: mean {m[0]:+.4f} = {m[1]:+.4f} {m[2]:+.4f}; '
              f'assigned once {m[3]:+.4f}; opposite signs on {cancel:.0%}; '
              f'|dc| > |loop| on {np.mean(np.abs(parts[:, 1]) > np.abs(parts[:, 2])):.0%}; '
              f'|once| < |recursion| on {np.mean(np.abs(parts[:, 3]) < np.abs(parts[:, 0])):.0%}, '
              f'ratio of mean |errors| {np.abs(parts[:, 3]).mean() / np.abs(parts[:, 0]).mean():.2f}')
    with open(OUT / 'tab-loopseries-dc.tex', 'w') as f:
        f.write('% generated by figs/loopseries.py; do not edit\n')
        f.write('\\begin{tabular}{lrrrrr}\n\\hline\\hline\n')
        f.write('instance & $\\beta J$ & recursion & double count & loop & assigned once\\\\\n\\hline\n')
        for name, bJ, tot, dcp, lp, once in rows:
            f.write(f'{name} & {bJ} & ${tot:+.4f}$ & ${dcp:+.4f}$ & ${lp:+.4f}$ & ${once:+.4f}$\\\\\n')
        f.write('\\hline\\hline\n\\end{tabular}\n')
    print('  wrote tab-loopseries-dc.tex')
    return rows


# -------------------------------------------- (3) small clustered instances

def clustered_instance(seed, n=10, triangles=4, extra=3):
    rng = np.random.default_rng(seed)
    while True:
        G = nx.Graph()
        G.add_nodes_from(range(n))
        for _ in range(triangles):
            a, b, c = rng.choice(n, 3, replace=False)
            G.add_edges_from([(a, b), (b, c), (a, c)])
        while G.number_of_edges() < 3 * triangles + extra:
            u, v = rng.choice(n, 2, replace=False)
            G.add_edge(int(u), int(v))
        G.remove_nodes_from([v for v in list(G) if G.degree(v) == 0])
        G = nx.convert_node_labels_to_integers(G)
        if nx.is_connected(G) and G.number_of_edges() <= 16:
            return G


def clustered_panel(ax, seeds=range(40)):
    rows = []
    for bJ, mk, col in ((0.3, 'o', DARK), (0.8, 's', MID)):
        pts = []
        for seed in seeds:
            G = clustered_instance(seed)
            cx, e = cliques_and_edges(G)
            eq, tq = series(BinaryFactorGraph.pairwise(e, bJ))
            ep, tp = series(BinaryFactorGraph.promoted(cx, e, bJ))
            pts.append((eq, ep, len(tq), len(tp)))
        pts = np.array(pts)
        rows.append((bJ, pts))
        ax.plot(pts[:, 0], pts[:, 1], mk, color=col, ms=3.2, mfc='white',
                mew=0.8, label=f'$\\beta J = {bJ}$')
        print(f'  clustered instances, bJ={bJ}: pairwise error {pts[:, 0].mean():.4f} '
              f'({pts[:, 2].mean():.0f} loops), promoted error {pts[:, 1].mean():.4f} '
              f'({pts[:, 3].mean():.1f} loops); promoted smaller on '
              f'{np.mean(pts[:, 1] < pts[:, 0]):.0%}, ratio of means '
              f'{pts[:, 1].mean() / pts[:, 0].mean():.2f}')
    lim = max(r[1][:, 0].max() for r in rows) * 1.05
    ax.plot([0, lim], [0, lim], '-', color=LIGHT, lw=0.8)
    ax.set_xlabel(r'$\ln Z - \ln Z_{\rm BP}$, pairwise', fontsize=8)
    ax.set_ylabel(r'$\ln Z - \ln Z_{\rm BP}$, promoted', fontsize=8)
    ax.set_xlim(0, lim)
    ax.set_ylim(0, lim)
    ax.legend(fontsize=7, frameon=False, loc='upper left')
    _tidy(ax)
    return rows


def ring_panel(ax, k=6):
    cx, e = cliques_and_edges(ring_of_triangles(k))
    for bJ, col in ((0.3, DARK), (0.8, MID)):
        eq, tq = series(BinaryFactorGraph.pairwise(e, bJ))
        ps = partial_sums(tq)
        total = ps[-1][1]
        xs = [s for s, _ in ps]
        ys = [v / total for _, v in ps]
        ax.step(xs, ys, where='post', color=col, lw=1.0,
                label=f'$\\beta J = {bJ}$, pairwise')
        ep, tp = series(BinaryFactorGraph.promoted(cx, e, bJ))
        ax.plot([2 * k], [1.0], 'o' if bJ == 0.3 else 's', color=col, ms=4,
                mfc='white', mew=0.9)
    ax.axhline(1.0, color=LIGHT, lw=0.8)
    ax.set_xlabel('loop size (legs)', fontsize=8)
    ax.set_ylabel('partial sum / whole', fontsize=8)
    ax.set_ylim(0, 1.08)
    ax.legend(fontsize=7, frameon=False, loc='lower right')
    _tidy(ax)


def figure():
    plt = _mpl()
    fig, axes = plt.subplots(2, 1, figsize=(3.4, 4.2))
    ring_panel(axes[0])
    clustered_panel(axes[1])
    for ax, tag in zip(axes, 'ab'):
        ax.text(-0.2, 1.02, f'({tag})', transform=ax.transAxes, fontsize=9)
    fig.tight_layout()
    fig.savefig(OUT / 'fig-loopseries.pdf')
    print('  wrote fig-loopseries.pdf')


if __name__ == '__main__':
    t0 = time.time()
    print('identity:'); check_identity()
    print('recursion:'); check_double_count_is_the_recursion()
    print('rings:'); table_rings()
    print('decomposition:'); table_decomposition()
    print('figure:'); figure()
    print(f'done in {time.time() - t0:.0f} s')
