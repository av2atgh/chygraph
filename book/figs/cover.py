"""Chapter 11: vertex cover, leaf removal and core percolation.

  fig-core     the leaf-removal core against chy-degree, for complexes of
               cardinality two, three and four: a transition at <k> = e in the
               first case and no transition at all in the others
  fig-hrg      the empirical payoff. Core of hyperbolic random graphs against
               mean degree, with the chygraph prediction evaluated on the
               measured clique ensemble and the degree-matched control

Chapter 11's first figure, leaf removal frame by frame, is TikZ in `cover.tex`.

The hyperbolic-graph points are read from `../probe/results/prediction4.csv`,
the cached output of `../probe/prediction4.py`; that script builds the graphs,
runs leaf removal on them and measures the clique ensemble, none of which the
book recomputes.
"""

import sys
from pathlib import Path

import numpy as np
from scipy.special import lambertw
from sympy import exp as sexp

sys.path.insert(0, str(Path.home() / 'av2atg' / 'statmech' / 'src'))
import statmech.hittingset as hs  # noqa: E402
from statmech import core as co  # noqa: E402
from statmech import cover as cv  # noqa: E402

OUT = Path(__file__).resolve().parent
PROBE = Path.home() / 'av2atg' / 'statmech' / 'probe' / 'results'
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


# ------------------------------------------------------------ (1) one point
def check_threshold():
    """<k> = e, reached from core percolation and from Chapter 10's map."""
    t = co.core_threshold(lambda k: co.graph(mean=k))
    assert abs(t - np.e) < 1e-9, t
    khs = hs.rsb_point([2], [1.0])
    print(f'  core-percolation threshold   <k> = {t:.10f}')
    print(f'  hard-field instability       <k> = {khs:.10f}')
    print(f'  e                            <k> = {np.e:.10f}')
    print(f'  the three agree to {abs(t - khs):.1e}: one point, computed twice')


# ---------------------------------------------- (2) the induced-graph cover
def check_induced_cover():
    """Vertex cover of the induced graph: same family, tighter constraint."""
    print('     <k>   Weigt-Hartmann   induced cover   hitting set')
    for k in (1.0, 2.0, 2.5):
        W = float(np.real(lambertw(k)))
        wh = 1 - (2 * W + W**2) / (2 * k)
        a = cv.poisson([2], [k]).cover_size()
        b = hs.poisson([2], [k]).cover_size()
        assert abs(a - wh) < 1e-9 and abs(b - wh) < 1e-9
        print(f'  {k:>6.1f}   {wh:>14.6f}   {a:>13.6f}   {b:>11.6f}')
    print('    the two problems coincide at cardinality two, as they must')
    x = hs.layer_symbols(1)
    tri = cv.CliqueCover([3], x[0]).cover_size()
    assert abs(tri - 0.5) < 1e-12
    print(f'  isolated triangles: hard field {tri:.4f}, truth {2/3:.4f}')
    print('    the mirror of Sec. 10.3 -- same defect, different problem, and')
    print('    the same wrong answer 1/2')


# ------------------------------------------- (3) no core-free branch above c=2
def check_core_free():
    """A layer of cardinality three or more removes the core-free branch."""
    for c in (2, 3, 4, 5, 8):
        free = co.clique_network(c, 1.0).has_core_free_branch()
        assert free == (c == 2), c
        print(f'  c = {c}: core-free branch {"exists" if free else "does not exist"}')
    # and then the core is exactly the non-isolated fraction
    print('  for a pure clique network of cardinality >= 3 the core is'
          ' 1 - e^{-<k>}:')
    for c in (3, 4, 6):
        for k in (0.3, 1.0, 2.5):
            got = co.clique_network(c, k).core_fraction()
            assert abs(got - (1 - np.exp(-k))) < 1e-9, (c, k)
        print(f'    c = {c}: every non-isolated vertex is core   OK')


def check_factorisation():
    """The factorisation that kills the core-free branch, symbolically."""
    from sympy import Symbol, expand, factor, simplify, symbols
    d = Symbol('delta')
    for c in range(2, 10):
        f = (c - 1) * d**(c - 2) - (c - 2) * d**(c - 1) - 1
        g = -(1 - d)**2 * sum((j + 1) * d**j for j in range(c - 2))
        assert simplify(expand(f - g)) == 0, c
    print('  f(delta) = -(1-delta)^2 sum_{j<c-2} (j+1) delta^j, '
          'for c = 2 ... 9   OK')
    print('    at c = 2 the sum is empty and f == 0: the condition is an '
          'identity')
    print('    at c >= 3 every coefficient is positive, so the only root is '
          'delta = 1')


def check_against_leaf_removal():
    """Eq. (11.4) against pure leaf removal on graphs built to match it."""
    path = PROBE / 'validate_core.txt'
    rows = []
    for line in path.read_text().splitlines()[2:]:
        f = line.split()
        if len(f) == 6:
            rows.append((int(f[0]), float(f[4]), float(f[5])))
    err = [e for _, sim, e in rows if sim > 1e-3]
    print(f'  {len(rows)} ensembles at n = 4e5, cardinalities '
          f'{sorted({c for c, _, _ in rows})}:')
    print(f'    error against pure leaf removal runs {min(err):.1e} to '
          f'{max(err):.1e}')
    print(f'    (from {path.name}, the cached output of probe/validate_core.py)')


def check_core_can_be_tiny():
    """The strong reading is false: the core need not be large."""
    x = hs.layer_symbols(2)
    for kL, kT in ((1.0, 1e-2), (1.0, 1e-3), (1.0, 1e-4)):
        phi = sexp(kL * (x[0] - 1)) * sexp(kT * (x[1] - 1))
        f = co.CorePercolation([2, 3], phi).core_fraction()
        print(f'  links <k> = {kL:g}, triangles <k> = {kT:g}:  core = {f:.3e}')
    print('    the algebra forbids a core-free branch; it does not make the'
          ' core big')


def figure_core():
    plt = _mpl()
    fig, ax = plt.subplots(figsize=(4.3, 2.6))
    ks = np.linspace(0.05, 6.0, 120)
    # c = 3 and c = 4 fall exactly on top of each other, which is the point;
    # they are drawn in different styles so that both are visible.
    for c, col, style, kw in ((2, DARK, '-', {}), (3, MID, '-', {}),
                              (4, LIGHT, '--', {'dashes': (4, 2)})):
        ax.plot(ks, [co.clique_network(c, k).core_fraction() for k in ks],
                style, lw=1.6, color=col, label=f'$c={c}$', **kw)
    ax.axvline(np.e, ls=':', lw=0.9, color=DARK)
    ax.annotate(r'$\langle k\rangle=e$', xy=(np.e + 0.13, 0.06), fontsize=7.5,
                color='0.3')
    ax.annotate(r'$1-e^{-\langle k\rangle}$: $c=3$ and $4$ coincide',
                xy=(1.55, 0.86), fontsize=7, color='0.35')
    ax.set_xlabel(r'chy-degree $\langle k\rangle$', fontsize=8.5)
    ax.set_ylabel(r'$P_C(\mathrm{VC},\,\mathrm{BP}_\chi)$', fontsize=9)
    ax.set_ylim(-0.02, 1.02)
    ax.legend(frameon=False, fontsize=8, loc='center right')
    _tidy(ax)
    fig.tight_layout()
    fig.savefig(OUT / 'fig-core.pdf')
    print(f'  wrote {OUT / "fig-core.pdf"}')


# ------------------------------------------------ (4) hyperbolic random graphs
def _hrg():
    import csv
    rows = []
    with open(PROBE / 'prediction4.csv') as fh:
        for r in csv.DictReader(fh):
            rows.append({k: (float(v) if k != 'seed' else int(v))
                         for k, v in r.items()})
    return rows


def _fit(k, y, kmax=1.0):
    """Power-law exponent of y ~ k^theta, fitted on k <= kmax."""
    k, y = np.asarray(k), np.asarray(y)
    m = (k <= kmax) & (y > 0)
    return float(np.polyfit(np.log(k[m]), np.log(y[m]), 1)[0])


def check_hrg():
    rows = _hrg()
    print('    tau   theta measured   theta chygraph   control')
    for tau in (2.5, 2.9, 2.1):
        r = [x for x in rows if abs(x['tau'] - tau) < 1e-9]
        k = [x['kbar'] for x in r]
        tm = _fit(k, [x['hrg_measured'] for x in r])
        tc = _fit(k, [x['hrg_chygraph'] for x in r])
        ctrl = max(x['cfg_measured'] for x in r if x['kbar'] <= 1.0)
        print(f'  {tau:>5.1f}   {tm:>13.3f}   {tc:>14.3f}   '
              f'{"no core" if ctrl <= 0 else f"{ctrl:.2e}"}')
    print('    the predicted exponent barely moves with tau; the measured one'
          ' does')
    # the one place the control does have a core
    r = [x for x in rows if abs(x['tau'] - 2.9) < 1e-9 and x['kbar'] > 5]
    cp = np.mean([x['cfg_chygraph'] for x in r])
    cm = np.mean([x['cfg_measured'] for x in r])
    print(f'  the control at tau = 2.9, <k> = 6 does have a core: predicted '
          f'{cp:.3f},\n    measured {cm:.3f} -- agreement to '
          f'{100*abs(cp-cm)/cm:.0f} per cent, same map, same pipeline')
    # how much of the real core the chygraph accounts for
    fr = [x['hrg_chygraph'] / x['hrg_measured'] for x in rows
          if x['hrg_measured'] > 0 and abs(x['tau'] - 2.9) < 1e-9]
    print(f'  fraction of the measured core accounted for at tau = 2.9: '
          f'{min(fr):.2f} to {max(fr):.2f}')


def figure_hrg():
    plt = _mpl()
    rows = _hrg()
    fig, axes = plt.subplots(1, 2, figsize=(4.6, 2.5), sharey=True)
    for ax, tau in zip(axes, (2.5, 2.9)):
        r = sorted((x for x in rows if abs(x['tau'] - tau) < 1e-9),
                   key=lambda z: z['kbar'])
        k = np.array([x['kbar'] for x in r])
        ax.loglog(k, [x['hrg_measured'] for x in r], 'o', ms=3.4, color=DARK,
                  label='measured')
        ax.loglog(k, [x['hrg_chygraph'] for x in r], 's', ms=3.4, mfc='white',
                  mew=0.9, color=DARK, label='chygraph')
        ctrl = np.array([x['cfg_measured'] for x in r])
        floor = 3e-6
        ax.loglog(k, np.where(ctrl > 0, ctrl, floor), '^', ms=3.2,
                  color=LIGHT, label='degree-matched')
        th = _fit(k, [x['hrg_chygraph'] for x in r])
        kk = np.array([k.min(), 1.0])
        c0 = r[0]['hrg_chygraph'] / k.min()**th
        ax.loglog(kk, c0 * kk**th, '-', lw=0.9, color='0.55')
        ax.set_title(rf'$\tau={tau}$', fontsize=8.5)
        ax.set_xlabel(r'mean degree $\langle k\rangle$', fontsize=8.5)
        ax.set_ylim(1e-6, 2.0)
        _tidy(ax)
    axes[0].set_ylabel('core fraction', fontsize=8.5)
    axes[0].legend(frameon=False, fontsize=6.6, loc='upper left')
    axes[1].annotate('control has no core here', xy=(0.055, 6e-6),
                     fontsize=6.2, color='0.45')
    fig.tight_layout()
    fig.savefig(OUT / 'fig-hrg.pdf')
    print(f'  wrote {OUT / "fig-hrg.pdf"}')


if __name__ == '__main__':
    print('one point, three names:')
    check_threshold()
    print('vertex cover of the induced graph:')
    check_induced_cover()
    print('above cardinality two:')
    check_core_free()
    check_factorisation()
    check_core_can_be_tiny()
    print('against pure leaf removal:')
    check_against_leaf_removal()
    print('hyperbolic random graphs:')
    check_hrg()
    print('figures:')
    figure_core()
    figure_hrg()
