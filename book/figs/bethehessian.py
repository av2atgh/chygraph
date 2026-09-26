"""Chapter 24: the Bethe Hessian of a chygraph (Sec. 24.6).

One figure and one table, and the numbers the section quotes.

  fig-bethehessian  (a) on one instance, the Perron root of the linearised
                    recursion less one and the smallest eigenvalue of the
                    vertex Bethe Hessian against the coupling: both cross
                    zero at the same point.  (b) the instance threshold
                    against the size of the instance, for three ensembles,
                    against the branching-matrix threshold of Eq. (8.x).
                    (c) how far the block of a complex is from the product
                    form at a random-field fixed point, by cardinality and
                    field spread; and at zero field with +-J couplings.
  tab-bethehessian  the three ensembles: the ensemble threshold, the
                    instance threshold at the largest size, the finite-size
                    gap.

Checks run first: the block at the trivial point is u'(c) for c = 2..5;
the weighted Ihara-Bass identity to 1e-12 at the trivial point, at a
polarised and at a random-field fixed point with complexes of at most three
members; its failure with a four-member complex; the Schur complement on a
graph is Saade's Bethe Hessian; the instance threshold by the Perron root
and by the Hessian gap agree to 1e-7.
"""

import json
import sys
import time
from itertools import combinations
from pathlib import Path

import networkx as nx
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'statmech' / 'src'))
sys.path.insert(0, str(ROOT / 'percolation' / 'src'))
from statmech.bethehessian import (  # noqa: E402
    BinaryFactorGraph, bethe_hessian, clique_factor_graph, edge_weights,
    factorise, instance_threshold, jacobian_block, linearised_operator,
    random_chygraph, smallest_eigenvalue, spectral_radius, trivial_blocks,
    trivial_vertex_hessian, vertex_hessian,
)
from statmech.ising import clique_derivative, critical_coupling  # noqa: E402

OUT = Path(__file__).resolve().parent
CACHE = ROOT / 'statmech' / 'probe' / 'results' / 'bethehessian.json'
DARK, MID, LIGHT = '0.10', '0.45', '0.70'
SPIN = np.array([1.0, -1.0])

ENSEMBLES = {
    r'Poisson graph, $\langle k\rangle=4$': dict(cards=[2], means=[4.0], regular=False),
    r'Poisson triangles, $\langle\kappa\rangle=2$': dict(cards=[3], means=[2.0], regular=False),
    'Poisson links, triangles, 4-cliques':
        dict(cards=[2, 3, 4], means=[1.5, 0.7, 0.3], regular=False),
}
REGULAR = {
    '4-regular graph': dict(cards=[2], means=[4.0], regular=True),
    'two triangles per atom': dict(cards=[3], means=[2.0], regular=True),
}
SIZES = [100, 200, 400, 800, 1600, 3200]
SEEDS = range(8)


def _mpl():
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    return plt


def _tidy(ax):
    ax.tick_params(labelsize=8)
    for sp in ('top', 'right'):
        ax.spines[sp].set_visible(False)


def ensemble_threshold(spec):
    exc = None
    if spec['regular']:
        exc = [m - 1 for m in spec['means']]
    return critical_coupling(spec['cards'], spec['means'], excess=exc)


# ---------------------------------------------------------------------------
# checks
# ---------------------------------------------------------------------------

def _identity(fg):
    T = linearised_operator(fg)
    alpha, beta, worst = edge_weights(fg)
    H, pref = bethe_hessian(fg, alpha, beta)
    lhs = np.linalg.det(np.eye(T.shape[0]) - T.toarray())
    return lhs, pref * np.linalg.det(H.toarray()), worst


def checks():
    print('checks')
    for c in range(2, 6):
        fg = clique_factor_graph([tuple(range(c))], 0.5)
        M = trivial_blocks(fg)[0]
        u = clique_derivative(c, 0.5)
        assert np.allclose(M[~np.eye(c, dtype=bool)], u, atol=1e-14)
        print(f"  c = {c}: block entry {M[0, 1]:.6f} = u'(c) {u:.6f}")
    rng = np.random.default_rng(1)
    mixed = [(0, 1, 2), (1, 2, 3), (3, 4), (2, 4, 5, 6), (5, 6, 7)]
    small = [(0, 1, 2), (1, 2, 3), (3, 4), (2, 4, 5), (5, 6, 7), (0, 7)]
    fg = clique_factor_graph(mixed, 0.4)
    blocks = trivial_blocks(fg)
    T = linearised_operator(fg, blocks)
    al, be, _ = edge_weights(fg, blocks)
    H, pref = bethe_hessian(fg, al, be)
    S, pref2 = trivial_vertex_hessian(fg, [M[0, 1] for M in blocks])
    lhs = np.linalg.det(np.eye(T.shape[0]) - T.toarray())
    print(f"  trivial point, c up to 4: det(I-T) = {lhs:.12f}, "
          f"identity {pref * np.linalg.det(H.toarray()):.12f}, "
          f"closed form {pref2 * np.linalg.det(S.toarray()):.12f}")
    fg = clique_factor_graph(small, 0.6, fields=rng.normal(0, 0.7, 8)).bp()
    lhs, rhs, worst = _identity(fg)
    print(f"  random field, c <= 3: {lhs:.12f} = {rhs:.12f}, "
          f"residual {worst:.1e}")
    fg = clique_factor_graph(small, 1.2)
    for k in fg.m_va:
        fg.m_va[k] = np.array([0.3, -0.3])
    fg.bp()
    lhs, rhs, worst = _identity(fg)
    print(f"  polarised, c <= 3, m = {fg.magnetisation(0):.3f}: "
          f"{lhs:.12f} = {rhs:.12f}, residual {worst:.1e}")
    fg = clique_factor_graph(mixed, 0.6, fields=rng.normal(0, 0.7, 8)).bp()
    lhs, rhs, worst = _identity(fg)
    print(f"  random field with a 4-clique: {lhs:.6f} vs {rhs:.6f}, "
          f"residual {worst:.3f}")
    G = nx.random_regular_graph(3, 20, seed=1)
    E = [tuple(e) for e in G.edges()]
    fg = clique_factor_graph(E, 0.4)
    t = np.tanh(0.4)
    S, _ = trivial_vertex_hessian(fg, [t] * len(E))
    A = nx.to_numpy_array(G, nodelist=fg.nodes)
    r = 1.0 / t
    Hs = (r * r - 1) * np.eye(20) - r * A + np.diag(A.sum(1))
    print(f"  Saade on a 3-regular graph: max deviation "
          f"{np.abs(S.toarray() - (t * t / (1 - t * t)) * Hs).max():.1e}")
    rng = np.random.default_rng(3)
    cx = random_chygraph(400, [2, 3], [1.5, 0.7], rng)
    bT, bH = instance_threshold(cx, via='T'), instance_threshold(cx, via='H')
    print(f"  n = 400 instance: beta* by Perron root {bT:.9f}, "
          f"by Hessian gap {bH:.9f}")
    for name, spec in REGULAR.items():
        rng = np.random.default_rng(7)
        cx = random_chygraph(102, spec["cards"], spec["means"], rng,
                             regular=True)
        bc = ensemble_threshold(spec)
        print(f"  {name}, n = 102: beta* = {instance_threshold(cx):.9f}, "
              f"ensemble {bc:.9f} (exact: the all-ones Perron vector)")


# ---------------------------------------------------------------------------
# panels
# ---------------------------------------------------------------------------

def panel_instance(ax, n=300, seed=5):
    spec = ENSEMBLES['Poisson links, triangles, 4-cliques']
    rng = np.random.default_rng(seed)
    cx = random_chygraph(n, spec['cards'], spec['means'], rng)
    fg = clique_factor_graph(cx, 1.0)
    cards = [len(sc) for sc, _ in fg.factors]
    bjs = np.linspace(0.05, 0.9, 35)
    rho, lam = [], []
    for bj in bjs:
        u = [clique_derivative(c, bj) for c in cards]
        blocks = [u[a] * (np.ones((c, c)) - np.eye(c))
                  for a, c in enumerate(cards)]
        rho.append(spectral_radius(linearised_operator(fg, blocks)) - 1)
        S, _ = trivial_vertex_hessian(fg, u)
        lam.append(smallest_eigenvalue(S))
    bstar = instance_threshold(cx)
    bc = ensemble_threshold(spec)
    ax.axhline(0, color=LIGHT, lw=0.6)
    ax.plot(bjs, rho, '-', color=DARK, lw=1.0,
            label=r'$\rho(T)-1$')
    ax.plot(bjs, lam, '--', color=MID, lw=1.0,
            label=r'$\lambda_{\min}(H_{V})$')
    ax.axvline(bstar, color=DARK, lw=0.5, ls=':')
    ax.axvline(bc, color=LIGHT, lw=0.8)
    ax.set_xlabel(r'$\beta J$', fontsize=8)
    ax.set_ylabel('stability', fontsize=8)
    ax.legend(fontsize=7, frameon=False, loc='center right')
    ax.set_title(f'(a) one instance, $n={n}$', fontsize=8, loc='left')
    _tidy(ax)
    return dict(n=n, seed=seed, beta_star=bstar, beta_c=bc,
                bjs=bjs.tolist(), rho=rho, lam=lam)


def thresholds(cache):
    key = 'thresholds'
    if key in cache:
        return cache[key]
    out = {}
    for name, spec in ENSEMBLES.items():
        bc = ensemble_threshold(spec)
        rows = {}
        for n in SIZES:
            t0 = time.time()
            vals = []
            for seed in SEEDS:
                rng = np.random.default_rng(1000 * n + seed)
                cx = random_chygraph(n, spec['cards'], spec['means'], rng,
                                     regular=spec['regular'])
                vals.append(instance_threshold(cx))
            rows[str(n)] = vals
            print(f"  {name}: n = {n}: beta* = {np.mean(vals):.5f} "
                  f"+- {np.std(vals):.5f} (ensemble {bc:.5f}), "
                  f"{time.time() - t0:.0f} s")
        out[name] = dict(beta_c=bc, sizes=rows)
    cache[key] = out
    return out


def panel_sizes(ax, res):
    markers = ['o', 's', '^']
    ax.axhline(0, color=LIGHT, lw=0.6)
    for (name, r), mk, shift in zip(res.items(), markers, (0.93, 1.0, 1.07)):
        bc = r['beta_c']
        ns = sorted(int(k) for k in r['sizes'])
        mean = [np.mean(r['sizes'][str(n)]) for n in ns]
        sd = [np.std(r['sizes'][str(n)]) for n in ns]
        ax.errorbar(np.array(ns) * shift, np.array(mean) - bc, yerr=sd,
                    fmt=mk + '-', color=DARK, mfc='white' if mk != 'o' else DARK,
                    ms=4, lw=0.8, capsize=2, label=name)
    ax.set_xscale('log')
    ax.set_xlabel('$n$', fontsize=8)
    ax.set_ylabel(r'$\beta^{*}J-\beta_{c}J$', fontsize=8)
    ax.legend(fontsize=6.5, frameon=False, loc='upper right')
    ax.set_title('(b) instance against ensemble', fontsize=8, loc='left')
    _tidy(ax)


def residual_field(c, sigma, bj, rng, reps=20):
    fg = clique_factor_graph([tuple(range(c))], bj)
    out = []
    for _ in range(reps):
        h = rng.normal(0, sigma, c)
        for v in range(c):
            fg.m_va[(v, 0)] = h[v] * SPIN
        out.append(factorise(jacobian_block(fg, 0))[2])
    return float(np.mean(out)), float(np.max(out))


def pm_clique(c, signs, bj):
    s = np.array([[1.0 if (i >> (c - 1 - b)) & 1 == 0 else -1.0
                   for b in range(c)] for i in range(2 ** c)])
    e = np.zeros(2 ** c)
    for k, (p, q) in enumerate(combinations(range(c), 2)):
        e += signs[k] * s[:, p] * s[:, q]
    return BinaryFactorGraph([(tuple(range(c)), (bj * e).reshape((2,) * c))])


def residual_disorder(c, bj, rng, reps=20):
    out = []
    for _ in range(reps):
        signs = rng.choice([-1.0, 1.0], c * (c - 1) // 2)
        fg = pm_clique(c, signs, bj)
        out.append(factorise(jacobian_block(fg, 0))[2])
    return float(np.mean(out)), float(np.max(out))


def panel_residual(ax, cache):
    key = 'residuals'
    if key not in cache:
        rng = np.random.default_rng(11)
        sig = np.array([0.05, 0.1, 0.2, 0.4, 0.8, 1.6])
        r = {'sigma': sig.tolist(), 'field': {}, 'disorder': {}}
        for c in (3, 4, 5):
            r['field'][str(c)] = [residual_field(c, s, 0.5, rng)[0] for s in sig]
        bjs = np.array([0.05, 0.1, 0.2, 0.4, 0.8, 1.6])
        r['bjs'] = bjs.tolist()
        for c in (3, 4, 5):
            r['disorder'][str(c)] = [residual_disorder(c, b, rng)[0]
                                     for b in bjs]
        cache[key] = r
    r = cache[key]
    sig = np.array(r['sigma'])
    for c, mk in zip((3, 4, 5), ('o', 's', '^')):
        y = np.maximum(np.array(r['field'][str(c)]), 1e-16)
        ax.plot(sig, y, mk + '-', color=DARK, ms=4, lw=0.8,
                mfc='white' if mk != 'o' else DARK, label=f'$c={c}$, field')
    for c, mk in zip((4, 5), ('s', '^')):
        y = np.maximum(np.array(r['disorder'][str(c)]), 1e-16)
        ax.plot(r['bjs'], y, mk + '--', color=MID, ms=4, lw=0.8, mfc='white',
                label=f'$c={c}$, $\\pm J$')
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_ylim(1e-17, 1)
    ax.set_xlabel(r'field spread $\sigma$ (solid), $\beta J$ (dashed)',
                  fontsize=8)
    ax.set_ylabel('residual of the product form', fontsize=8)
    ax.legend(fontsize=6.5, frameon=False, loc='center left')
    ax.set_title('(c) the block of one complex', fontsize=8, loc='left')
    _tidy(ax)
    for c in (3, 4, 5):
        print(f"  c = {c}: field residual at sigma = 0.4: "
              f"{r['field'][str(c)][3]:.2e}, at 1.6: {r['field'][str(c)][5]:.2e}"
              + (f"; +-J at beta J = 0.4: {r['disorder'][str(c)][3]:.2e}"
                 if str(c) in r['disorder'] else ''))


def table(res):
    lines = [r'\begin{tabular}{@{}lccc@{}}', r'\hline\hline',
             r'ensemble & $\beta_{c}J$, Eq.~\eqref{eq:branch} & '
             r'$\beta^{*}J$, $n=3200$ & gap\\', r'\hline']
    for name, r in res.items():
        bc = r['beta_c']
        v = r['sizes'][str(SIZES[-1])]
        lines.append(f"{name} & {bc:.5f} & ${np.mean(v):.5f}\\pm{np.std(v):.5f}$"
                     f" & {abs(np.mean(v) - bc) / bc * 100:.2f}\\%\\\\")
    lines += [r'\hline\hline', r'\end{tabular}']
    (OUT / 'tab-bethehessian.tex').write_text('\n'.join(lines) + '\n')
    print('\n'.join(lines))


def main():
    t0 = time.time()
    checks()
    cache = json.loads(CACHE.read_text()) if CACHE.exists() else {}
    plt = _mpl()
    fig, axes = plt.subplots(3, 1, figsize=(3.4, 7.2))
    print('panel (a)')
    inst = panel_instance(axes[0])
    print(f"  beta* = {inst['beta_star']:.5f}, ensemble {inst['beta_c']:.5f}")
    print('panel (b)')
    res = thresholds(cache)
    CACHE.write_text(json.dumps(cache))
    panel_sizes(axes[1], res)
    print('panel (c)')
    panel_residual(axes[2], cache)
    CACHE.write_text(json.dumps(cache))
    fig.tight_layout()
    fig.savefig(OUT / 'fig-bethehessian.pdf')
    table(res)
    print(f'done in {time.time() - t0:.0f} s')


if __name__ == '__main__':
    main()
