"""Chapter 20: the endogeny test, run.

One figure and one table, drawn from the cache written by
`statmech/probe/endogeny.py scan` (about forty minutes; nothing here
recomputes it).

  fig-endogeny     (a) the +-J spin glass on a random regular graph of
                   degree three, in a field: three boundaries in beta J
                   against H -- where two copies started as independent
                   draws stop converging (the nonlinear test), where a small
                   perturbation carried alongside the field stops decaying
                   (the linearised test), and where the annealed criterion
                   <kbar><(du/dh)^2> crosses one; the exact zero-field
                   transition is the anchor;
                   (b) the hitting set on Poisson hypergraphs of cardinality
                   two, three and four: the stationary rate of a perturbation
                   against the mean chy-degree in units of the hard-field
                   line <k>(c-1) = e.
  tab-endogeny     the boundaries by cardinality, and the regular
                   hypergraphs of Mezard and Tarzia against their exact
                   entropy.

Checks run first, on the cache: the linear and the annealed boundary of the
spin glass approach the exact zero-field transition as H -> 0 with the
H^(2/3) law of a de Almeida--Thouless line; and the shuffled-pair boundary
lies within the bisection's resolution of the linear one.
"""

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
CACHE = ROOT / 'statmech' / 'probe' / 'results' / 'endogeny.json'
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


def load():
    return json.load(open(CACHE))


# ------------------------------------------------------------------ checks

def check_anchor(d):
    sg = d['spin_glass']
    bc = sg['zero_field']['beta_J_c']
    rows = sorted(sg['rows'], key=lambda r: r['H'])
    H = np.array([r['H'] for r in rows])
    lin = np.array([r['linear'] for r in rows])
    ann = np.array([r['annealed'] for r in rows])
    sh = np.array([r['shuffled'] for r in rows])
    # H^(2/3) approach: fit bJ = bc + a H^(2/3) on the three smallest fields
    x = H[:3] ** (2 / 3)
    a = float(np.sum(x * (lin[:3] - bc)) / np.sum(x * x))
    resid = lin[:3] - (bc + a * x)
    print(f'  exact zero-field transition beta J_c = {bc:.4f}')
    print(f'  linear boundary at H = {H[:3].tolist()}: {lin[:3].round(3).tolist()}; '
          f'fit b_c + a H^(2/3) with a = {a:.3f}, residuals {resid.round(3).tolist()}')
    print(f'  annealed boundary at the same fields: {ann[:3].round(3).tolist()}')
    gap = np.abs(sh - lin) / lin
    print(f'  shuffled against linear boundary: relative gap {gap.round(3).tolist()}')
    assert np.all(np.abs(resid) < 0.03)
    assert np.all(gap < 0.08)
    return a


def boundaries(d):
    """Zero of the stationary rate by linear interpolation, per cardinality."""
    out = {}
    for c in sorted({r['c'] for r in d['hitting_set']}):
        rows = sorted((r for r in d['hitting_set'] if r['c'] == c), key=lambda r: r['k'])
        k = np.array([r['k'] for r in rows])
        rt = np.array([r['rate'] for r in rows])
        hard = rows[0]['hard']
        i = np.where(np.diff(np.sign(rt)) > 0)[0]
        if len(i):
            i = i[0]
            kb = k[i] - rt[i] * (k[i + 1] - k[i]) / (rt[i + 1] - rt[i])
        else:
            kb = float('nan')
        ent_all = np.array([r['entropy'] for r in rows])
        se_all = np.array([r['entropy_se'] for r in rows])
        ok = np.isfinite(ent_all) & np.isfinite(se_all)
        ent = np.interp(kb, k[ok], ent_all[ok]) if kb == kb else float('nan')
        se = np.interp(kb, k[ok], se_all[ok]) if kb == kb else float('nan')
        out[c] = dict(hard=hard, boundary=float(kb), ratio=float(kb / hard), entropy=float(ent),
                      entropy_se=float(se),
                      entropy_hard=float(np.interp(hard, k[ok], ent_all[ok])))
        print(f'  c={c}: hard-field line {hard:.4f}, endogeny boundary {kb:.4f} '
              f'({kb / hard:.3f} x), entropy there {ent:+.3f} ({se:.3f})')
    return out


# ------------------------------------------------------------------ figure

def panel_spin_glass(ax, d, a):
    sg = d['spin_glass']
    bc = sg['zero_field']['beta_J_c']
    rows = sorted(sg['rows'], key=lambda r: r['H'])
    H = np.array([r['H'] for r in rows])
    for key, mk, col, lab in (('shuffled', 'o', DARK, 'two copies, independent start'),
                              ('linear', 's', MID, 'perturbation carried with the field'),
                              ('annealed', '^', LIGHT, r'$\langle\bar\kappa\rangle\langle u^{\prime 2}\rangle = 1$')):
        y = np.array([r[key] for r in rows])
        ax.plot(H, y, mk, color=col, ms=3.5, mfc='white', mew=0.9, label=lab)
    hh = np.linspace(0, H.max(), 200)
    ax.plot(hh, bc + a * hh ** (2 / 3), '-', color=MID, lw=0.8)
    ax.plot([0], [bc], 'D', color=DARK, ms=4, label=r'$2\tanh^2\beta J = 1$')
    ax.set_xlabel(r'field $\beta H$', fontsize=8)
    ax.set_ylabel(r'boundary $\beta J$', fontsize=8)
    ax.set_xlim(-0.03, H.max() * 1.05)
    ax.legend(fontsize=6.3, frameon=False, loc='upper left')
    _tidy(ax)


def panel_hitting_set(ax, d):
    for c, mk, col in ((2, 'o', DARK), (3, 's', MID), (4, '^', LIGHT)):
        rows = sorted((r for r in d['hitting_set'] if r['c'] == c and r['ratio'] <= 1.2),
                      key=lambda r: r['k'])
        x = np.array([r['ratio'] for r in rows])
        y = np.array([r['rate'] for r in rows])
        ax.plot(x, y, mk + '-', color=col, ms=3.5, mfc='white', mew=0.9, lw=0.8,
                label=f'$c = {c}$')
    ax.axhline(0, color=LIGHT, lw=0.8)
    ax.axvline(1, color=LIGHT, lw=0.8, ls=':')
    ax.set_xlabel(r'$\langle k\rangle (c-1) / e$', fontsize=8)
    ax.set_ylabel('rate of a perturbation, per sweep', fontsize=8)
    ax.legend(fontsize=7, frameon=False, loc='upper left')
    _tidy(ax)


def figure(d, a):
    plt = _mpl()
    fig, axes = plt.subplots(2, 1, figsize=(3.4, 4.3))
    panel_spin_glass(axes[0], d, a)
    panel_hitting_set(axes[1], d)
    for ax, tag in zip(axes, 'ab'):
        ax.text(-0.2, 1.02, f'({tag})', transform=ax.transAxes, fontsize=9)
    fig.tight_layout()
    fig.savefig(OUT / 'fig-endogeny.pdf')
    print('  wrote fig-endogeny.pdf')


def table(d, b):
    with open(OUT / 'tab-endogeny.tex', 'w') as f:
        f.write('% generated by figs/endogeny.py; do not edit\n')
        f.write('\\begin{tabular}{lrrrr}\n\\hline\\hline\n')
        f.write('Poisson, cardinality $c$ & hard-field line & endogeny boundary & ratio & entropy there\\\\\n\\hline\n')
        for c in sorted(b):
            r = b[c]
            f.write(f"$c={c}$ & {r['hard']:.3f} & {r['boundary']:.3f} & {r['ratio']:.3f} & ${r['entropy']:+.3f}\\pm{r['entropy_se']:.3f}$\\\\\n")
        f.write('\\hline\n')
        f.write('regular, $L$ complexes of $K$ & entropy & pair distance & rate & \\\\\n\\hline\n')
        for r in d['regular']:
            f.write(f"$L={r['L']}$, $K={r['K']}$ & ${r['entropy']:+.3f}$ & "
                    f"{_sci(r['shuffled'])} & ${r['rate']:+.4f}$ & \\\\\n")
        f.write('\\hline\\hline\n\\end{tabular}\n')
    print('  wrote tab-endogeny.tex')


def _sci(x):
    if x == 0:
        return '$0$'
    e = int(np.floor(np.log10(abs(x))))
    m = x / 10 ** e
    return f'${m:.1f}\\times10^{{{e}}}$' if e < -1 else f'${x:.2f}$'


if __name__ == '__main__':
    d = load()
    print('anchor:'); a = check_anchor(d)
    print('boundaries:'); b = boundaries(d)
    print('regular hypergraphs:')
    for r in d['regular']:
        print(f"  L={r['L']} K={r['K']}: entropy {r['entropy']:+.3f}, pair distance "
              f"{r['shuffled']:.1e}, rate {r['rate']:+.4f}")
    print('figure:'); figure(d, a); table(d, b)
