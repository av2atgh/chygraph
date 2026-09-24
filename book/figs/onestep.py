"""Sec. 13.10: one-step replica symmetry breaking at m != 0 for 3-SAT.

  fig-onestep   left: the complexity Sigma(m) against m at several clause
                densities; right: the Parisi parameter m*(alpha) at which
                Sigma vanishes, with the block density of Sec. 8.10 drawn
                across for the comparison the conjecture asks for

Data are read from ../statmech/probe/results/onestep_sat.csv, the cached
output of ../statmech/probe/onestep_sat.py scan (population of populations,
M = 1000 surveys of P = 200 messages, 300 sweeps, 200000 samples for the
potential; about five minutes per (alpha, m) point).  The checks recompute
the small anchors: BP's Bethe free entropy equals ln(#solutions) on a tree
instance, and the m = 1 potential on a replica-symmetric population equals
the replica-symmetric free entropy.
"""
import csv
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'statmech' / 'probe'))

OUT = Path(__file__).resolve().parent
PROBE = ROOT / 'statmech' / 'probe' / 'results'
DARK, MID, LIGHT = '0.10', '0.45', '0.70'
ALPHA_S = 4.267        # Mertens, Mezard, Zecchina 2006; Montanari et al. 2008
ALPHA_D = 3.86         # Krzakala et al. 2007


def _mpl():
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    return plt


def _tidy(ax):
    ax.tick_params(labelsize=8)
    for sp in ('top', 'right'):
        ax.spines[sp].set_visible(False)


def _rows():
    rows = list(csv.DictReader(open(PROBE / 'onestep_sat.csv')))
    extra = PROBE / 'onestep_sat_extra.csv'
    if extra.exists():
        rows += list(csv.DictReader(open(extra)))
    return rows


SPREAD_SOFT = 0.085   # a population with a smaller survey spread at alpha >= 4 has frozen


def is_soft(r):
    """A population that has not frozen: wide surveys, a finite potential
    and a complexity in the physical range.  Frozen populations have a
    spread below 0.085 and either a divergent potential or a complexity of
    order 0.1 to 0.8, a hundred times the physical one."""
    sg = float(r['Sigma'])
    return (np.isfinite(sg) and abs(sg) < FROZEN
            and (float(r['alpha']) < 4.0 or float(r['m']) < 0.5 or float(r['spread']) >= SPREAD_SOFT))


def load():
    rows = _rows()
    g = defaultdict(list)
    for r in rows:
        g[(float(r['alpha']), float(r['m']))].append(r)
    out = {}
    for k, rs in g.items():
        soft = [r for r in rs if is_soft(r)]
        out[k] = {c: (float(np.mean([float(r[c]) for r in soft])) if soft else np.nan)
                  for c in ('F', 's', 'Sigma', 'spread')}
        out[k]['n'] = len(soft)
        out[k]['frozen'] = len(rs) - len(soft)
        out[k]['Sigma_sd'] = float(np.std([float(r['Sigma']) for r in soft])) if len(soft) > 1 else np.nan
    return out


NOISE = 5e-4     # |Sigma| below this is indistinguishable from zero at this population size


def m_bracket(curve, frozen_at=None, noise=NOISE):
    """[m_low, m_high]: the largest m at which Sigma is clearly positive and
    the smallest later m at which it is clearly negative, or at which every
    population has frozen, whichever comes first.  (nan, nan) if Sigma is
    never clearly positive; (m_low, 1.0) if it is never clearly negative."""
    pts = [(m, sg) for m, sg in curve if np.isfinite(sg)]
    pos = [m for m, sg in pts if sg > noise]
    if not pos:
        return np.nan, np.nan
    lo = max(pos)
    neg = [m for m, sg in pts if m > lo and sg < -noise]
    hi = min(neg) if neg else 1.0
    if frozen_at is not None and lo < frozen_at < hi:
        hi = frozen_at
    return lo, hi
FROZEN = 1e-2    # |Sigma| above this only occurs where the population has frozen


def m_star(curve, noise=NOISE):
    """Zero of Sigma(m) by linear interpolation between the last m at which
    Sigma is clearly positive and the first later m at which it is clearly
    negative; nan if the curve never leaves the noise band on both sides.

    curve: sorted list of (m, Sigma), nan entries ignored.
    """
    pts = [(m, sg) for m, sg in curve if np.isfinite(sg) and abs(sg) < FROZEN]
    for i, (m0, s0) in enumerate(pts):
        if s0 <= noise:
            continue
        for m1, s1 in pts[i + 1:]:
            if s1 < -noise:
                return m0 + (m1 - m0) * s0 / (s0 - s1)
            if s1 > noise:
                m0, s0 = m1, s1
    return np.nan


def figure_onestep(block_density=None):
    plt = _mpl()
    data = load()
    alphas = sorted({k[0] for k in data})
    fig, (ax, bx) = plt.subplots(2, 1, figsize=(4.3, 4.8))
    shades = np.linspace(0.75, 0.05, len(alphas))
    ms_star = []
    for a, sh in zip(alphas, shades):
        curve = sorted((k[1], v['Sigma']) for k, v in data.items() if k[0] == a)
        ms_star.append((a, None, curve))
        if a < 3.86:
            continue
        # a frozen population -- surveys collapsed onto hard fields, potential
        # diverging -- shows up as |Sigma| far outside the physical range; those
        # points are not drawn, and the curve stops where they begin
        kept = [(m, sg) for m, sg in curve if np.isfinite(sg) and abs(sg) < FROZEN]
        ax.plot([m for m, _ in kept], [1e3 * sg for _, sg in kept], 'o-', ms=3.5, color=str(sh),
                mfc='white', label=rf'$\alpha={a:g}$')
    ax.axhline(0, color=LIGHT, lw=0.8, ls=':')
    ax.set_xlabel('$m$', fontsize=8)
    ax.set_ylabel(r'complexity $\Sigma(m)\times10^{3}$', fontsize=8)
    ax.set_ylim(-6, 8)
    ax.legend(fontsize=6.5, frameon=False, ncol=2, loc='upper left')
    _tidy(ax)
    # m*(alpha) as a bracket: 1 where Sigma(1) is within noise of zero, else
    # [last clearly positive m, first clearly negative or frozen m]
    first = True
    for a, ms, curve in ms_star:
        frozen_at = min((k[1] for k, v in data.items() if k[0] == a and v['n'] == 0), default=None)
        lo, hi = m_bracket(curve, frozen_at)
        s1 = dict(curve).get(1.0, np.nan)
        if np.isfinite(s1) and s1 >= -NOISE and not np.isfinite(lo):
            lo = hi = 1.0
        if not np.isfinite(lo):
            continue
        bx.plot([a, a], [lo, hi], '-', color=DARK, lw=1.6, solid_capstyle='butt',
                label='$m^{*}(\\alpha)$, bracketed' if first else None)
        bx.plot([a], [lo], 'o', ms=4, color=DARK, mfc='white')
        bx.plot([a], [hi], 'o', ms=4, color=DARK, mfc='white')
        first = False
    if block_density is not None:
        bx.axhline(block_density, color=MID, lw=1.0, ls='--',
                   label='blocks per free variable, Sec. 8.10')
    bx.axvline(ALPHA_S, color=LIGHT, lw=0.8, ls=':')
    bx.axvline(ALPHA_D, color=LIGHT, lw=0.8, ls=':')
    bx.set_xlabel(r'clause density $\alpha$', fontsize=8)
    bx.set_ylabel('$m^{*}$', fontsize=8)
    bx.set_ylim(-0.02, 1.05)
    bx.set_xlim(3.78, 4.3)
    bx.legend(fontsize=7, frameon=False, loc='lower left')
    _tidy(bx)
    fig.tight_layout()
    fig.savefig(OUT / 'fig-onestep.pdf')
    plt.close(fig)
    return ms_star


def load_by_seed():
    rows = _rows()
    out = defaultdict(dict)
    for r in rows:
        if is_soft(r):
            out[int(r['seed'])][(float(r['alpha']), float(r['m']))] = float(r['Sigma'])
    return out


def table_onestep():
    data = load()
    by_seed = load_by_seed()
    alphas = sorted({k[0] for k in data})
    ms = sorted({k[1] for k in data})
    print('Sigma(m) x 1000, populations averaged; columns m = ' + ' '.join(f'{m:.1f}' for m in ms))
    for a in alphas:
        print(f'{a:5.2f} ' + ' '.join(f'{1000 * data[(a, m)]["Sigma"]:+6.2f}' if (a, m) in data else '   .  ' for m in ms)
              + '   n=' + ','.join(str(data[(a, m)]['n']) for m in ms if (a, m) in data))
    print('frozen populations (soft/total) per alpha: ' + '; '.join(
        f'{a:g}: {sum(data[k]["n"] for k in data if k[0] == a)}/{sum(data[k]["n"] + data[k]["frozen"] for k in data if k[0] == a)}'
        for a in alphas))
    print('alpha  Sigma(1)  F(1)  s(1)  spread(1)  first frozen m   m* bracket')
    mstars = []
    for a in alphas:
        curve = sorted((k[1], v['Sigma']) for k, v in data.items() if k[0] == a)
        d = dict(curve)
        v1 = data.get((a, 1.0), {})
        frozen_at = min((k[1] for k, v in data.items() if k[0] == a and v['n'] == 0), default=None)
        lo, hi = m_bracket(curve, frozen_at)
        if np.isfinite(d.get(1.0, np.nan)) and d[1.0] >= -NOISE and not np.isfinite(lo):
            lo, hi = 1.0, 1.0          # uncondensed at this precision
        mstars.append((a, lo, hi))
        print(f'{a:5.2f}  {d.get(1.0, np.nan):+.4f}  {v1.get("F", np.nan):.4f}  {v1.get("s", np.nan):.4f}  '
              f'{v1.get("spread", np.nan):.4f}  {frozen_at}   [{lo:.2f}, {hi:.2f}]')
    # seed scatter of Sigma at fixed (alpha, m), over the condensed points
    if len(by_seed) >= 2:
        keys = [k for k in data if k[0] >= 3.9 and data[k]['n'] >= 2 and abs(data[k]['Sigma']) < FROZEN]
        sds = [np.nanstd([bs[k] for bs in by_seed.values() if k in bs and abs(bs[k]) < FROZEN]) for k in keys]
        ses = [sd / np.sqrt(data[k]['n']) for sd, k in zip(sds, keys)]
        print(f'Sigma scatter between populations (alpha >= 3.9): mean sd {np.mean(sds):.5f}, max {np.max(sds):.5f}; '
              f'mean standard error of the mean {np.mean(ses):.5f}')
    print('block density 0.55 inside the bracket at alpha = ' + ', '.join(f'{a:g}' for a, lo, hi in mstars if np.isfinite(lo) and lo <= 0.55 <= hi))
    return mstars


def check_anchors():
    import onestep_sat as o
    o.check_bp_exact_on_tree()
    alpha = 3.5
    eta = o.rs_population(alpha, size=100000, sweeps=300, seed=0)
    frs = o.rs_free_entropy(eta, alpha, nsamp=400000)
    sub = eta[np.random.default_rng(1).integers(0, len(eta), 1000)]
    one = o.OneStep(alpha, 1.0, M=1000, P=200, seed=0)
    one.pop = np.repeat(sub[:, None], 200, axis=1)
    r = one.potential(200000)
    print(f'alpha={alpha}: F_RS = {frs:.4f}; one-step potential on the same messages at m=1: {r["F"]:.4f}, Sigma = {r["Sigma"]:+.1e}')
    assert abs(frs - r['F']) < 0.02


if __name__ == '__main__':
    check_anchors()
    table_onestep()
    figure_onestep(block_density=0.55)
    print('wrote fig-onestep.pdf')
