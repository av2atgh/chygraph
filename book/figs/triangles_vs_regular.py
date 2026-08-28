"""Chapter 4: do triangles help percolation or hurt it?

The published comparison (Vazquez 2024, Sec. 4.3) puts a Poisson graph with
triangles against a Poisson graph with the same number of links and no
triangles, and finds the triangles help: dtheta = 2 q^2 (1-q) <k>_tri > 0.
Chapter 9 finds the opposite sign for the Ising model, comparing two triangles
per vertex against the 4-regular graph.

The two comparisons are not the same comparison.  Putting links into triangles
raises a node's degree by two at a time, so at fixed link count the degree
distribution is broader -- and it is the excess degree, not the triangle, that
is doing the work.  Held at fixed degree the sign reverses.

This script establishes that three ways:

  * exact enumeration of the intra-complex excess component size of a triangle
    under bond percolation, checking sbar_tri(q) = 2q(1 + q - q^2);
  * the two thresholds from the chygraph branching condition, 1/3 for the
    4-regular graph against 0.403 for two triangles per vertex;
  * Monte Carlo bond percolation on finite instances of both.

Writes fig-triangles.pdf into the book directory.
"""

from itertools import product
from pathlib import Path

import networkx as nx
import numpy as np

OUT = Path(__file__).resolve().parent
RNG = np.random.default_rng(20260827)


# --------------------------------------------------------------- exact pieces
def sbar_triangle_exact(q):
    """Mean excess component size inside a bond-percolated triangle.

    Arrive at node A; B and C are the others; the three links are present
    independently with probability q.  Enumerate all eight configurations and
    count how many of {B, C} are reachable from A.
    """
    tot = 0.0
    for ab, ac, bc in product((0, 1), repeat=3):
        pr = q ** (ab + ac + bc) * (1 - q) ** (3 - ab - ac - bc)
        reach = set()
        if ab:
            reach.add('B')
        if ac:
            reach.add('C')
        if bc and ('B' in reach or 'C' in reach):
            reach |= {'B', 'C'}
        tot += pr * len(reach)
    return tot


def sbar_triangle_closed(q):
    return 2 * q * (1 + q - q * q)


def theta_regular(q):
    """4-regular graph: one link layer, kappa = 4, excess 3, sbar = q."""
    return 3 * q - 1


def theta_triangles(q):
    """Two triangles per vertex: one triangle layer, kappa = 2, excess 1."""
    return 1 * sbar_triangle_closed(q) - 1


def threshold(f, lo=1e-6, hi=1.0):
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if f(mid) < 0:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


# ------------------------------------------------------------- the two graphs
def two_triangles_per_vertex(n, rng):
    """Each node in exactly two triangles, triangles meeting in <= 1 node."""
    assert n % 3 == 0
    for _ in range(200):
        stubs = rng.permutation(np.repeat(np.arange(n), 2))
        tris = stubs.reshape(-1, 3)
        if any(len(set(t)) < 3 for t in tris):
            continue
        g, ok = nx.Graph(), True
        for t in tris:
            a, b, c = t
            if g.has_edge(a, b) or g.has_edge(a, c) or g.has_edge(b, c):
                ok = False   # two triangles sharing an edge: not treelike
                break
            g.add_edges_from([(a, b), (a, c), (b, c)])
        if ok and g.number_of_nodes() == n:
            return g
    raise RuntimeError('no clean triangle graph found')


def giant(g, q, rng):
    h = nx.Graph()
    h.add_nodes_from(g)
    e = np.array(g.edges())
    keep = rng.random(len(e)) < q
    h.add_edges_from(map(tuple, e[keep]))
    return len(max(nx.connected_components(h), key=len)) / g.number_of_nodes()


def main():
    print('--- exact check of sbar_tri(q) ---')
    for q in (0.0, 0.25, 0.5, 0.75, 1.0):
        a, b = sbar_triangle_exact(q), sbar_triangle_closed(q)
        print(f'  q={q:.2f}  enumerated={a:.6f}  closed form={b:.6f}  '
              f'{"ok" if abs(a - b) < 1e-12 else "MISMATCH"}')

    qc_reg = threshold(theta_regular)
    qc_tri = threshold(theta_triangles)
    print(f'\n--- thresholds ---\n  4-regular graph        q_c = {qc_reg:.4f}'
          f'\n  two triangles / vertex q_c = {qc_tri:.4f}'
          f'\n  triangles raise it by {100 * (qc_tri / qc_reg - 1):.1f} per cent')

    print('\n--- Monte Carlo ---')
    n, seeds = 30000, 4
    qs = np.linspace(0.20, 0.70, 18)
    mc = {}
    for name, build in (('4-regular', lambda r: nx.random_regular_graph(4, n, seed=int(r.integers(1 << 30)))),
                        ('triangles', lambda r: two_triangles_per_vertex(n, r))):
        curves = []
        for s in range(seeds):
            rng = np.random.default_rng(1000 + s)
            g = build(rng)
            curves.append([giant(g, q, rng) for q in qs])
        mc[name] = np.mean(curves, axis=0)
        print(f'  {name}: <k> = {2 * g.number_of_edges() / n:.2f}, '
              f'C = {nx.average_clustering(g):.3f}')
    figure(qs, mc, qc_reg, qc_tri)
    return qc_reg, qc_tri


def figure(qs, mc, qc_reg, qc_tri):
    """Two stacked panels: the measured order parameter, and the theta that
    predicts where it lifts off.  Kept apart because they are different
    quantities and sharing an axis invites reading one as the other."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(4.3, 4.3), sharex=True,
                                   gridspec_kw=dict(height_ratios=[1.7, 1]))
    fine = np.linspace(0.20, 0.70, 300)

    # ---- (a) measured giant component
    ax1.plot(qs, mc['4-regular'], 'o-', ms=4, lw=0.8, mfc='white', mec='0.55',
             color='0.55', mew=0.9, label='4-regular (no triangles)')
    ax1.plot(qs, mc['triangles'], 's-', ms=4, lw=0.8, color='0.15',
             label='two triangles per vertex')
    for qc, c in ((qc_reg, '0.55'), (qc_tri, '0.15')):
        ax1.axvline(qc, color=c, lw=0.8, ls=':')
    ax1.annotate(f'$q_c={qc_reg:.3f}$', (qc_reg, 0.42), xytext=(-4, 0),
                 textcoords='offset points', fontsize=7.5, ha='right', color='0.35')
    ax1.annotate(f'$q_c={qc_tri:.3f}$', (qc_tri, 0.20), xytext=(5, 0),
                 textcoords='offset points', fontsize=7.5, ha='left', color='0.15')
    ax1.set_ylabel('giant component $S$')
    ax1.set_ylim(-0.03, 1.03)
    ax1.legend(fontsize=7.5, frameon=False, loc='upper left')
    ax1.set_title('(a) measured, $n=3\\times10^{4}$', fontsize=8.5, loc='left')

    # ---- (b) the prediction
    ax2.plot(fine, [theta_regular(q) for q in fine], '-', color='0.55', lw=1.1)
    ax2.plot(fine, [theta_triangles(q) for q in fine], '-', color='0.15', lw=1.1)
    ax2.axhline(0, color='0.8', lw=0.7)
    for qc, c in ((qc_reg, '0.55'), (qc_tri, '0.15')):
        ax2.axvline(qc, color=c, lw=0.8, ls=':')
    ax2.set_xlabel('bond occupation probability $q$')
    ax2.set_ylabel(r'$\theta$')
    ax2.set_xlim(0.20, 0.70); ax2.set_ylim(-0.45, 0.75)
    ax2.set_title('(b) predicted', fontsize=8.5, loc='left')

    for ax in (ax1, ax2):
        ax.tick_params(labelsize=8)
        for sp in ('top', 'right'):
            ax.spines[sp].set_visible(False)
    fig.tight_layout()
    fig.savefig(OUT / 'fig-triangles.pdf')
    print('  wrote fig-triangles.pdf')


if __name__ == '__main__':
    main()
