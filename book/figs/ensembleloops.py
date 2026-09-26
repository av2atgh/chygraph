"""Chapter 24: the loop series on the ensemble (Sec. 24.3).

One figure and the numbers the section quotes.

  fig-ensembleloops  E[ln Z - ln Z_BP] against the coupling for two random
                     chygraphs, from the closed form (solid), its first
                     order -1/2 ln det(I - B) - 1/2 tr B (dashed), the mean
                     of the cycle sum on instances of 300 atoms (filled
                     markers) and the mean of the exact ln Z - ln Z_BP by
                     enumeration on instances of 18 atoms (open markers);
                     the vertical line is the threshold of Eq. (8.det).

Checks run first: on six instances of sixteen atoms the cycle sum equals
the exact ln Z - ln Z_BP whenever the cycles are disjoint and the full
loop series equals it always; the truncated series approaches the closed
form as the cut-off grows.
"""

import json
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'statmech' / 'src'))
sys.path.insert(0, str(ROOT / 'percolation' / 'src'))
from statmech.bethehessian import clique_factor_graph, random_chygraph  # noqa: E402
from statmech.ensembleloops import (  # noqa: E402
    cycle_sum, ensemble_series, exact_minus_bethe,
)
from statmech.ising import critical_coupling  # noqa: E402
from statmech.loopseries import loop_series  # noqa: E402

OUT = Path(__file__).resolve().parent
CACHE = ROOT / 'statmech' / 'probe' / 'results' / 'ensembleloops.json'
DARK, MID, LIGHT = '0.10', '0.45', '0.70'

ENSEMBLES = {
    'links and triangles': dict(cards=[2, 3], means=[1.2, 0.5]),
    'links, triangles, 4-cliques': dict(cards=[2, 3, 4], means=[1.0, 0.4, 0.15]),
}
N_BIG, SEEDS_BIG, LMAX = 300, 400, 10
N_SMALL, SEEDS_SMALL = 18, 300


def _mpl():
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    return plt


def _tidy(ax):
    ax.tick_params(labelsize=8)
    for sp in ('top', 'right'):
        ax.spines[sp].set_visible(False)


def checks():
    print('checks')
    cards, means = [2, 3], [1.2, 0.5]
    rng = np.random.default_rng(1)
    for s in range(6):
        cx = random_chygraph(16, cards, means, rng)
        fg = clique_factor_graph(cx, 0.25)
        ex = exact_minus_bethe(fg)
        cyc, counts = cycle_sum(cx, 0.25, lmax=16)
        full = float('nan')
        if len(cx) <= 13:
            fg.bp()
            full = np.log1p(sum(r for _, _, r in loop_series(fg)))
        print(f"  n = 16, seed {s}: exact {ex:.6f}, cycles {cyc:.6f}, "
              f"full series {full:.6f}, cycles by length {counts}")
    for lmax in (6, 10, 14):
        print(f"  truncated at {lmax}: {ensemble_series(cards, means, 0.35, lmax=lmax)[0]:.5f}"
              f" against closed form {ensemble_series(cards, means, 0.35)[0]:.5f}")


def scan(cache, name, spec):
    if name in cache:
        return cache[name]
    cards, means = spec['cards'], spec['means']
    bc = critical_coupling(cards, means)
    bjs = np.linspace(0.05, 0.92 * bc, 8)
    out = dict(beta_c=bc, bjs=bjs.tolist(), big=[], big_se=[], small=[],
               small_se=[])
    for bj in bjs:
        t0 = time.time()
        rng = np.random.default_rng(int(1e4 * bj))
        vals = [cycle_sum(random_chygraph(N_BIG, cards, means, rng), bj,
                          lmax=LMAX)[0] for _ in range(SEEDS_BIG)]
        out['big'].append(float(np.mean(vals)))
        out['big_se'].append(float(np.std(vals) / np.sqrt(len(vals))))
        rng = np.random.default_rng(int(1e4 * bj) + 1)
        vals = [exact_minus_bethe(clique_factor_graph(
            random_chygraph(N_SMALL, cards, means, rng), bj))
            for _ in range(SEEDS_SMALL)]
        out['small'].append(float(np.mean(vals)))
        out['small_se'].append(float(np.std(vals) / np.sqrt(len(vals))))
        clos, first = ensemble_series(cards, means, bj)
        print(f"  {name}, beta J = {bj:.3f}: closed {clos:.4f} (first order "
              f"{first:.4f}), n = {N_BIG} cycles {out['big'][-1]:.4f} +- "
              f"{out['big_se'][-1]:.4f}, n = {N_SMALL} exact "
              f"{out['small'][-1]:.4f} +- {out['small_se'][-1]:.4f} "
              f"({time.time() - t0:.0f} s)")
    cache[name] = out
    return out


def small_cycles(cache, name, spec):
    """The cycle sum on the same eighteen-atom instances as the exact
    enumeration, added to a cached scan."""
    out = cache[name]
    if 'small_cycles' in out:
        return out
    cards, means = spec['cards'], spec['means']
    out['small_cycles'], out['small_cycles_se'] = [], []
    for bj in out['bjs']:
        rng = np.random.default_rng(int(1e4 * bj) + 1)
        vals = [cycle_sum(random_chygraph(N_SMALL, cards, means, rng), bj,
                          lmax=N_SMALL)[0] for _ in range(SEEDS_SMALL)]
        out['small_cycles'].append(float(np.mean(vals)))
        out['small_cycles_se'].append(float(np.std(vals) / np.sqrt(len(vals))))
    return out


def table(cache):
    lines = [r'\begin{tabular}{@{}lcccccc@{}}', r'\hline\hline',
             r'ensemble & $\beta J$ & closed form & $\ell\le10$ & cycles, $n=300$ & '
             r'exact, $n=18$ & cycles, $n=18$\\', r'\hline']
    for name, spec in ENSEMBLES.items():
        r = cache[name]
        for k in (1, 4, 6, 7):
            bj = r['bjs'][k]
            clos = ensemble_series(spec['cards'], spec['means'], bj)[0]
            trunc = ensemble_series(spec['cards'], spec['means'], bj, lmax=LMAX)[0]
            lines.append(
                f"{name if k == 1 else ''} & {bj:.3f} & {clos:.3f} & {trunc:.3f} & "
                f"${r['big'][k]:.3f}\\pm{r['big_se'][k]:.3f}$ & "
                f"${r['small'][k]:.4f}\\pm{r['small_se'][k]:.4f}$ & "
                f"${r['small_cycles'][k]:.4f}\\pm{r['small_cycles_se'][k]:.4f}$\\\\")
    lines += [r'\hline\hline', r'\end{tabular}']
    (OUT / 'tab-ensembleloops.tex').write_text('\n'.join(lines) + '\n')
    print('\n'.join(lines))


def panel(ax, name, spec, res, letter):
    cards, means = spec['cards'], spec['means']
    bc = res['beta_c']
    grid = np.linspace(0.02, 0.985 * bc, 120)
    clos = [ensemble_series(cards, means, b)[0] for b in grid]
    first = [ensemble_series(cards, means, b)[1] for b in grid]
    ax.plot(grid, clos, '-', color=DARK, lw=1.0, label='closed form')
    trunc = [ensemble_series(cards, means, b, lmax=LMAX)[0] for b in grid]
    ax.plot(grid, trunc, ':', color=DARK, lw=1.0, label=f'closed form, $\\ell \\leq {LMAX}$')
    ax.plot(grid, first, '--', color=MID, lw=1.0,
            label=r'$-\frac{1}{2}\ln\det(I-B)-\frac{1}{2}\,\mathrm{tr}\,B$')
    ax.errorbar(res['bjs'], res['big'], yerr=res['big_se'], fmt='o', color=DARK,
                ms=4, capsize=2, lw=0.8, label=f'cycles, $n={N_BIG}$')
    ax.errorbar(res['bjs'], res['small'], yerr=res['small_se'], fmt='s',
                color=DARK, mfc='white', ms=4, capsize=2, lw=0.8,
                label=f'exact, $n={N_SMALL}$')
    ax.axvline(bc, color=LIGHT, lw=0.8)
    ax.set_xlabel(r'$\beta J$', fontsize=8)
    ax.set_ylabel(r'$\mathbb{E}[\ln Z-\ln Z_{\mathrm{BP}}]$', fontsize=8)
    ax.set_ylim(0, 1.6)
    ax.legend(fontsize=6.5, frameon=False, loc='upper left')
    ax.set_title(f'({letter}) {name}', fontsize=8, loc='left')
    _tidy(ax)


def main():
    t0 = time.time()
    checks()
    cache = json.loads(CACHE.read_text()) if CACHE.exists() else {}
    plt = _mpl()
    fig, axes = plt.subplots(1, 2, figsize=(6.0, 2.6))
    for ax, (name, spec), letter in zip(axes, ENSEMBLES.items(), 'ab'):
        print(name)
        res = scan(cache, name, spec)
        CACHE.write_text(json.dumps(cache))
        small_cycles(cache, name, spec)
        CACHE.write_text(json.dumps(cache))
        panel(ax, name, spec, res, letter)
        for k, bj in enumerate(res['bjs']):
            print(f"  beta J = {bj:.3f}: n = 18 exact {res['small'][k]:.4f}, "
                  f"cycles on the same instances {res['small_cycles'][k]:.4f}")
    fig.tight_layout()
    fig.savefig(OUT / 'fig-ensembleloops.pdf')
    table(cache)
    print(f'done in {time.time() - t0:.0f} s')


if __name__ == '__main__':
    main()
