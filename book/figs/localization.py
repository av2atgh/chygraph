"""Localization chapter: Anderson mobility edge and LLT percolation on
regular chygraphs with triangles.

  fig-anderson   the two lines at the bottom of the band against disorder,
                 one panel per degree: the mobility edge E_c^loc (solid) and
                 the LLT percolation threshold E_c^perc (dashed) for each
                 (s, t), with the band bottom E_min dotted
  tab-anderson-wc  the band-centre critical disorder W_c per ensemble

Data: ../statmech/probe/results/anderson_lines.csv and anderson_wc.csv,
the cached outputs of ../statmech/probe/anderson.py scan.  The checks
recompute the small things: the Schur-complement cavity against exact
inversion on incidence trees (1e-15) and on random instances (error falling
with the imaginary part), and the pure-hopping band bottoms -2 sqrt(2) and
-2 sqrt(3) for (3,0) and (4,0).
"""
import csv
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'statmech' / 'probe'))
import anderson as A  # noqa: E402

OUT = Path(__file__).resolve().parent
PROBE = ROOT / 'statmech' / 'probe' / 'results'
DARK, MID, LIGHT = '0.10', '0.45', '0.70'
STYLE = {(3, 0): (DARK, 'o'), (1, 1): (LIGHT, '^'),
         (4, 0): (DARK, 'o'), (2, 1): (MID, 's'), (0, 2): (LIGHT, '^')}


def _mpl():
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    return plt


def _tidy(ax):
    ax.tick_params(labelsize=8)
    for sp in ('top', 'right'):
        ax.spines[sp].set_visible(False)


def load():
    rows = list(csv.DictReader(open(PROBE / 'anderson_lines.csv')))
    by = defaultdict(list)
    for r in rows:
        by[(int(r['s']), int(r['t']))].append({k: float(v) for k, v in r.items()})
    for k in by:
        by[k].sort(key=lambda r: r['W'])
        for r in by[k]:
            # the bisection runs up to the band centre; landing there means every
            # state is localised at this W (W above the band-centre W_c): no edge
            if r['Ec_loc'] > -0.01:
                r['Ec_loc'] = np.nan
                r['gap'] = np.nan
    return by


def figure_anderson():
    plt = _mpl()
    by = load()
    fig, axes = plt.subplots(2, 1, figsize=(4.3, 5.2))
    for ax, degree in zip(axes, (3, 4)):
        for (s, t), rows in sorted(by.items()):
            if s + 2 * t != degree:
                continue
            col, mk = STYLE[(s, t)]
            W = [r['W'] for r in rows]
            ax.plot(W, [r['Ec_loc'] for r in rows], '-', marker=mk, ms=4, color=col, mfc='white',
                    label=rf'$(s,t)=({s},{t})$ mobility edge')
            ax.plot(W, [r['Ec_perc'] for r in rows], '--', marker=mk, ms=4, color=col, mfc=col,
                    label=rf'$(s,t)=({s},{t})$ LLT percolation')
            ax.plot(W, [r['Emin'] for r in rows], ':', color=col, lw=0.8)
        ax.set_title(f'degree {degree}', fontsize=9)
        ax.set_ylabel('$E$ at the bottom of the band', fontsize=8)
        ax.legend(fontsize=6, frameon=False, ncol=1)
        _tidy(ax)
    axes[1].set_xlabel('disorder $W$', fontsize=8)
    fig.tight_layout()
    fig.savefig(OUT / 'fig-anderson.pdf')
    plt.close(fig)


def table_wc():
    rows = list(csv.DictReader(open(PROBE / 'anderson_wc.csv')))
    rows.sort(key=lambda r: (int(r['degree']), -int(r['t'])))
    lines = [r'\begin{tabular}{ccrr}', r'\hline\hline',
             r'$(s,t)$ & degree & band bottom & $W_{c}$ at $E=0$\\', r'\hline']
    bb = {(int(r['s']), int(r['t'])): float(r['band_bottom']) for r in
          csv.DictReader(open(PROBE / 'anderson_lines.csv'))}
    for r in rows:
        k = (int(r['s']), int(r['t']))
        lines.append(f'$({k[0]},{k[1]})$ & {int(r["degree"])} & ${bb[k]:.4f}$ & ${float(r["Wc"]):.2f}$\\\\')
    lines += [r'\hline\hline', r'\end{tabular}']
    (OUT / 'tab-anderson-wc.tex').write_text('\n'.join(lines) + '\n')


def table_lines():
    by = load()
    lines = [r'\begin{tabular}{ccrrrr}', r'\hline\hline',
             r'$(s,t)$ & $W$ & $E_{\min}$ & $E_{c}^{\mathrm{loc}}$ & $E_{c}^{\mathrm{perc}}$ & gap\\', r'\hline']
    for (s, t), rows in sorted(by.items(), key=lambda kv: (kv[0][0] + 2 * kv[0][1], -kv[0][1])):
        for r in rows:
            loc = f'${r["Ec_loc"]:.3f}$' if np.isfinite(r['Ec_loc']) else 'none'
            gap = f'${r["gap"]:+.3f}$' if np.isfinite(r['gap']) else '---'
            lines.append(f'$({s},{t})$ & {r["W"]:g} & ${r["Emin"]:.3f}$ & {loc} & ${r["Ec_perc"]:.3f}$ & {gap}\\\\')
    lines += [r'\hline\hline', r'\end{tabular}']
    (OUT / 'tab-anderson-lines.tex').write_text('\n'.join(lines) + '\n')


def check_cavity():
    rng = np.random.default_rng(1)

    def complex_tree(s, t, depth):
        complexes, n, frontier = [], 1, [(0, 0)]
        while frontier:
            a, g = frontier.pop()
            if g >= depth:
                continue
            for _ in range(s if g == 0 else s - 1):
                complexes.append((a, n)); frontier.append((n, g + 1)); n += 1
            for _ in range(t if g == 0 else t - 1):
                complexes.append((a, n, n + 1)); frontier += [(n, g + 1), (n + 1, g + 1)]; n += 2
        return n, complexes
    worst = 0.0
    for (s, t, depth) in ((3, 0, 5), (2, 1, 4), (1, 2, 3)):
        n, cx = complex_tree(s, t, depth)
        eps = rng.uniform(-1.5, 1.5, n)
        H = A.hamiltonian(n, cx, eps)
        z = 0.7 - 0.2j
        Gex = np.diag(np.linalg.inv(H - z * np.eye(n)))
        G, _, _ = A.cavity_instance(n, cx, eps, z, sweeps=2 * depth + 5)
        worst = max(worst, np.abs(G - Gex).max())
        print(f'incidence tree (s,t)=({s},{t}), n={n}: max |G_cavity - G_exact| = {np.abs(G - Gex).max():.1e}')
    assert worst < 1e-12
    for (s, t) in ((3, 0), (4, 0)):
        bb = A.Ensemble(s, t, 1.0, P=10).band_bottom()
        print(f'band bottom ({s},{t}) = {bb:.6f} against -2 sqrt(K) = {-2 * np.sqrt(s - 1):.6f}')
        assert abs(bb + 2 * np.sqrt(s - 1)) < 1e-4


if __name__ == '__main__':
    check_cavity()
    figure_anderson()
    table_wc()
    table_lines()
    print('wrote fig-anderson.pdf, tab-anderson-wc.tex, tab-anderson-lines.tex')
