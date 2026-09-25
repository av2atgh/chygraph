"""Localization chapter: Anderson mobility edge and LLT percolation on
regular chygraphs with triangles.

  fig-anderson-mechanism  left: lambda(1/2) against E - E_min at W = 3 on the
                          degree-three tree and cactus, the mobility edge as a
                          crossing of one; right: the percolation growth factor
                          against E' on the tree at W = 1.5, with the (Sigma,
                          eta, p) triple drawn jointly and separately
  fig-anderson            the two lines measured from the bottom of the band,
                          against disorder, one panel per degree
  fig-anderson-gap        the gap E_c^perc - E_c^loc against disorder, five
                          ensembles in one panel
  fig-anderson-centre     lambda(1/2) at the band centre against disorder:
                          where each ensemble localises entirely
  tab-anderson-wc         band bottoms and band-centre critical disorder

Data: ../statmech/probe/results/anderson_{lines,wc,stability,growth,centre}.csv,
the cached outputs of ../statmech/probe/anderson.py scan and side.  The
checks recompute the small things: the Schur-complement cavity against exact
inversion on incidence trees (1e-15) and the pure-hopping band bottoms
-2 sqrt(2) and -2 sqrt(3) for (3,0) and (4,0).
"""
import csv
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
for _p in (ROOT / 'statmech' / 'src', ROOT / 'percolation' / 'src'):
    sys.path.insert(0, str(_p))
from statmech import Chygraph  # noqa: E402
from statmech import resolvent as R  # noqa: E402

OUT = Path(__file__).resolve().parent
PROBE = ROOT / 'statmech' / 'probe' / 'results'
DARK, MID, LIGHT = '0.10', '0.45', '0.70'
STYLE = {(3, 0): (DARK, 'o'), (1, 1): (LIGHT, '^'),
         (4, 0): (DARK, 'o'), (2, 1): (MID, 's'), (0, 2): (LIGHT, '^')}
GAPSTYLE = {(3, 0): (DARK, 'o', '-'), (1, 1): (DARK, '^', '--'),
            (4, 0): (LIGHT, 'o', '-'), (2, 1): (LIGHT, 's', '-.'), (0, 2): (LIGHT, '^', '--')}


def _mpl():
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    return plt


def _tidy(ax):
    ax.tick_params(labelsize=8)
    for sp in ('top', 'right'):
        ax.spines[sp].set_visible(False)


def _rows(name):
    return [{k: float(v) for k, v in r.items()} for r in csv.DictReader(open(PROBE / f'anderson_{name}.csv'))]


def load():
    by = defaultdict(list)
    for r in _rows('lines'):
        by[(int(r['s']), int(r['t']))].append(r)
    for k in by:
        by[k].sort(key=lambda r: r['W'])
        for r in by[k]:
            # the bisection runs up to the band centre; landing there means every
            # state is localised at this W (W above the band-centre W_c): no edge
            if r['Ec_loc'] > -0.01:
                r['Ec_loc'] = np.nan
                r['gap'] = np.nan
    return by


def _lab(k):
    return rf'$(s,t)=({k[0]},{k[1]})$'


def figure_mechanism():
    plt = _mpl()
    fig, (ax, bx) = plt.subplots(1, 2, figsize=(4.6, 2.5))
    st = defaultdict(list)
    for r in _rows('stability'):
        st[(int(r['s']), int(r['t']))].append(r)
    for k, rows in sorted(st.items()):
        rows.sort(key=lambda r: r['E'])
        Emin = Chygraph([2, 3], [k[0], k[1]], regular=True).spectral_bottom(3.0)
        col, mk = STYLE[k]
        ax.plot([r['E'] - Emin for r in rows], [r['lam'] for r in rows], '-', marker=mk, ms=3.5,
                color=col, mfc='white', label=_lab(k))
    ax.axhline(1, color=LIGHT, lw=0.8, ls=':')
    ax.set_xlabel(r'$E-E_{\min}$', fontsize=8)
    ax.set_ylabel(r'$\lambda(1/2)$ at $W=3$', fontsize=8)
    ax.legend(fontsize=6.5, frameon=False, loc='upper left')
    _tidy(ax)
    g = sorted(_rows('growth'), key=lambda r: r['Eprime'])
    bx.plot([r['Eprime'] for r in g], [r['growth_joint'] for r in g], 'o-', ms=3.5, color=DARK, mfc='white',
            label='triple drawn jointly')
    bx.plot([r['Eprime'] for r in g], [r['growth_indep'] for r in g], 's--', ms=3.5, color=LIGHT, mfc='white',
            label='$p$ drawn separately')
    bx.axhline(1, color=LIGHT, lw=0.8, ls=':')
    bx.set_xlabel(r"$E'$", fontsize=8)
    bx.set_ylabel(r'percolation growth, $(3,0)$, $W=1.5$', fontsize=8)
    bx.legend(fontsize=6.5, frameon=False, loc='upper left')
    _tidy(bx)
    fig.tight_layout()
    fig.savefig(OUT / 'fig-anderson-mechanism.pdf')
    plt.close(fig)


def figure_lines():
    plt = _mpl()
    by = load()
    fig, axes = plt.subplots(2, 1, figsize=(4.3, 5.0))
    for ax, degree in zip(axes, (3, 4)):
        for k, rows in sorted(by.items(), key=lambda kv: -kv[0][1]):
            if k[0] + 2 * k[1] != degree:
                continue
            col, mk = STYLE[k]
            W = [r['W'] for r in rows]
            ax.plot(W, [r['Ec_loc'] - r['Emin'] for r in rows], '-', marker=mk, ms=4, color=col, mfc='white',
                    label=_lab(k) + ' mobility edge')
            ax.plot(W, [r['Ec_perc'] - r['Emin'] for r in rows], '--', marker=mk, ms=4, color=col, mfc=col,
                    label=_lab(k) + ' LLT percolation')
        ax.set_title(f'degree {degree}', fontsize=9)
        ax.set_ylabel(r'$E_{c}-E_{\min}$, above the band bottom', fontsize=8)
        ax.legend(fontsize=6, frameon=False, loc='upper left')
        _tidy(ax)
    axes[1].set_xlabel('disorder $W$', fontsize=8)
    fig.tight_layout()
    fig.savefig(OUT / 'fig-anderson.pdf')
    plt.close(fig)


def figure_phase():
    """The (E, W) plane at the bottom of the band, per degree: bulk edge,
    isolated eigenvalue below W_min, LLT percolation line and mobility
    edge.  The bottom-edge half of Tonetti et al.'s Fig. S7."""
    plt = _mpl()
    by = load()
    fig, axes = plt.subplots(1, 2, figsize=(4.6, 2.6))
    for ax, degree in zip(axes, (3, 4)):
        Wmax = 12.6
        for k, rows in sorted(by.items(), key=lambda kv: -kv[0][1]):
            if k[0] + 2 * k[1] != degree:
                continue
            col, mk = STYLE[k]
            emin = rows[0]['band_bottom']
            iso = -degree
            Wmin = 2 * (degree + emin)
            Wg = np.linspace(Wmin, Wmax, 50)
            ax.plot(emin - Wg / 2, Wg, ':', color=col, lw=0.8)
            ax.plot([iso, iso], [0, Wmin], '-', color=col, lw=2.2, solid_capstyle='butt')
            ax.plot([r['Ec_perc'] for r in rows], [r['W'] for r in rows], '--', marker=mk, ms=4,
                    color=col, mfc=col, label=f'({k[0]},{k[1]}) percolation')
            ax.plot([r['Ec_loc'] for r in rows], [r['W'] for r in rows], '-', marker=mk, ms=4,
                    color=col, mfc='white', label=f'({k[0]},{k[1]}) mobility edge')
        ax.set_ylim(0, Wmax)
        ax.set_title(f'degree {degree}', fontsize=9)
        ax.set_xlabel('energy $E$', fontsize=8)
        ax.legend(fontsize=5.5, frameon=False, loc='lower left')
        _tidy(ax)
    axes[0].set_ylabel('disorder $W$', fontsize=8)
    fig.tight_layout()
    fig.savefig(OUT / 'fig-anderson-phase.pdf')
    plt.close(fig)


def figure_gap():
    plt = _mpl()
    by = load()
    fig, ax = plt.subplots(figsize=(4.3, 2.8))
    for k, rows in sorted(by.items(), key=lambda kv: (kv[0][0] + 2 * kv[0][1], -kv[0][1])):
        col, mk, ls = GAPSTYLE[k]
        pts = [(r['W'], r['gap']) for r in rows if np.isfinite(r['gap'])]
        ax.plot([w for w, _ in pts], [g for _, g in pts], ls, marker=mk, ms=4, color=col, mfc='white',
                label=_lab(k) + (', degree 3' if k[0] + 2 * k[1] == 3 else ', degree 4'))
        if len(pts) < len(rows):        # the band localised entirely beyond the last point
            w0 = pts[-1][0]
            ax.annotate('all localised', xy=(w0 + 0.3, pts[-1][1]), fontsize=6, color=col, va='center')
    ax.axhline(0, color=LIGHT, lw=0.8, ls=':')
    ax.set_xlabel('disorder $W$', fontsize=8)
    ax.set_ylabel(r'$E_{c}^{\mathrm{perc}}-E_{c}^{\mathrm{loc}}$', fontsize=8)
    ax.legend(fontsize=6, frameon=False, loc='upper left', ncol=2)
    _tidy(ax)
    fig.tight_layout()
    fig.savefig(OUT / 'fig-anderson-gap.pdf')
    plt.close(fig)


def figure_centre():
    plt = _mpl()
    c = defaultdict(list)
    for r in _rows('centre'):
        c[(int(r['s']), int(r['t']))].append(r)
    fig, ax = plt.subplots(figsize=(4.3, 2.8))
    for k, rows in sorted(c.items(), key=lambda kv: (kv[0][0] + 2 * kv[0][1], -kv[0][1])):
        rows.sort(key=lambda r: r['W'])
        col, mk, ls = GAPSTYLE[k]
        ax.plot([r['W'] for r in rows], [r['lam'] for r in rows], ls, marker=mk, ms=4, color=col, mfc='white',
                label=_lab(k))
    ax.axhline(1, color=LIGHT, lw=0.8, ls=':')
    ax.set_xlabel('disorder $W$', fontsize=8)
    ax.set_ylabel(r'$\lambda(1/2)$ at the band centre', fontsize=8)
    ax.legend(fontsize=6.5, frameon=False, ncol=2)
    _tidy(ax)
    fig.tight_layout()
    fig.savefig(OUT / 'fig-anderson-centre.pdf')
    plt.close(fig)


def table_wc():
    rows = list(csv.DictReader(open(PROBE / 'anderson_wc.csv')))
    rows.sort(key=lambda r: (int(r['degree']), -int(r['t'])))
    lines = [r'\begin{tabular}{ccrrr}', r'\hline\hline',
             r'$(s,t)$ & degree & band bottom & branching & $W_{c}$ at $E=0$\\', r'\hline']
    bb = {(int(r['s']), int(r['t'])): float(r['band_bottom']) for r in
          csv.DictReader(open(PROBE / 'anderson_lines.csv'))}
    for r in rows:
        k = (int(r['s']), int(r['t']))
        s_, t_ = k
        # branching of the incidence tree, Eq. branchingloc: the leading eigenvalue
        br = max(abs(np.linalg.eigvals(Chygraph([2, 3], [s_, t_], regular=True).incidence_branching())))
        lines.append(f'$({k[0]},{k[1]})$ & {int(r["degree"])} & ${bb[k]:.4f}$ & ${br:.2f}$ & ${float(r["Wc"]):.1f}$\\\\')
    lines += [r'\hline\hline', r'\end{tabular}']
    (OUT / 'tab-anderson-wc.tex').write_text('\n'.join(lines) + '\n')


def print_lines():
    """The numbers behind the figures, for the record."""
    by = load()
    print('(s,t)  W   E_min   E_c^loc   E_c^perc   gap')
    for k, rows in sorted(by.items(), key=lambda kv: (kv[0][0] + 2 * kv[0][1], -kv[0][1])):
        for r in rows:
            print(f"{k}  {r['W']:4.1f}  {r['Emin']:7.3f}  {r['Ec_loc']:7.3f}  {r['Ec_perc']:7.3f}  {r['gap']:+6.3f}")


def check_cavity():
    """The recursion against exact inversion on incidence trees, and the pure
    hopping band bottoms against their closed forms -- both from
    statmech.resolvent, whose own tests pin the same things."""
    rng = np.random.default_rng(1)
    worst = 0.0
    for (s, t, depth) in ((3, 0, 5), (2, 1, 4), (1, 2, 3)):
        n, cx = R.incidence_tree([2, 3], [s, t], depth)
        eps = rng.uniform(-1.5, 1.5, n)
        H = R.hamiltonian(n, cx, eps)
        z = 0.7 - 0.2j
        Gex = np.diag(np.linalg.inv(H - z * np.eye(n)))
        G, _, _ = R.cavity_instance(n, cx, eps, z, sweeps=2 * depth + 5)
        worst = max(worst, np.abs(G - Gex).max())
        print(f'incidence tree (s,t)=({s},{t}), n={n}: max |G_cavity - G_exact| = {np.abs(G - Gex).max():.1e}')
    assert worst < 1e-12
    for (s, t), want in (((3, 0), -2 * np.sqrt(2)), ((4, 0), -2 * np.sqrt(3)), ((0, 2), -(2 * np.sqrt(2) + 1))):
        bb = Chygraph([2, 3], [s, t], regular=True).band_bottom()
        print(f'band bottom ({s},{t}) = {bb:.6f} against {want:.6f}')
        assert abs(bb - want) < 1e-3


if __name__ == '__main__':
    check_cavity()
    print_lines()
    figure_lines()
    figure_gap()
    figure_phase()
    table_wc()
    if (PROBE / 'anderson_centre.csv').exists():
        figure_mechanism()
        figure_centre()
    print('wrote fig-anderson.pdf, fig-anderson-gap.pdf, fig-anderson-phase.pdf, fig-anderson-mechanism.pdf, fig-anderson-centre.pdf, tab-anderson-wc.tex')
