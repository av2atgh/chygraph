"""Manuscript figures generated from probe results.

Fig. 2: leaf-removal core of hyperbolic random graphs against mean degree, on
log-log axes, with the chygraph prediction from each graph's own measured
clique ensemble and the degree-matched configuration model as control.

Palette validated for colour-vision deficiency (categorical, three slots);
identity is carried by marker shape as well as hue so the figure survives
greyscale printing.
"""

import csv
from collections import defaultdict
from pathlib import Path

import matplotlib as mpl
import numpy as np

mpl.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402

HERE = Path(__file__).parent
CSV = HERE / 'probe' / 'results' / 'prediction4.csv'
OUT = HERE / 'fig_prediction4.pdf'
OUT_SIM = HERE / 'fig_simplicial.pdf'

MEAS, CHY, CTRL = '#0072B2', '#D55E00', '#8C6BB1'
FLOOR = 3e-5                     # control is exactly zero; drawn on the floor
FIT_BELOW = 1.0

mpl.rcParams.update({
    'font.size': 8, 'axes.labelsize': 8, 'xtick.labelsize': 7,
    'ytick.labelsize': 7, 'legend.fontsize': 7, 'axes.linewidth': 0.6,
    'xtick.major.width': 0.6, 'ytick.major.width': 0.6,
    'font.family': 'serif', 'mathtext.fontset': 'dejavuserif',
})


def load():
    rows = []
    with CSV.open() as fh:
        for r in csv.DictReader(fh):
            for k, v in r.items():
                if k != 'seed':
                    try:
                        r[k] = float(v)
                    except ValueError:
                        pass
            rows.append(r)
    g = defaultdict(list)
    for r in rows:
        g[(r['tau'], r['kbar_target'])].append(r)
    out = defaultdict(lambda: defaultdict(list))
    for (tau, kt), v in g.items():
        m = lambda key: float(np.mean([x[key] for x in v]))
        out[tau]['kbar'].append(m('kbar'))
        for key in ('hrg_chygraph', 'hrg_measured', 'cfg_measured'):
            out[tau][key].append(m(key))
    for tau in out:
        order = np.argsort(out[tau]['kbar'])
        for key in out[tau]:
            out[tau][key] = np.asarray(out[tau][key])[order]
    return out


def slope(k, y):
    ok = (k > 0) & (y > 1e-9) & (k <= FIT_BELOW)
    if ok.sum() < 3:
        return None
    return np.polyfit(np.log(k[ok]), np.log(y[ok]), 1)


def main():
    d = load()
    taus = sorted(d)
    fig, axes = plt.subplots(1, 3, figsize=(7.05, 2.45), sharey=True)

    for ax, tau in zip(axes, taus):
        k = d[tau]['kbar']
        meas, chy, ctrl = (d[tau]['hrg_measured'], d[tau]['hrg_chygraph'],
                           d[tau]['cfg_measured'])

        for y, col in ((meas, MEAS), (chy, CHY)):
            f = slope(k, y)
            if f is not None:
                xs = np.array([k.min() * 0.85, 1.15])
                ax.plot(xs, np.exp(f[1]) * xs ** f[0], color=col,
                        lw=1.1, alpha=0.55, zorder=1)

        ax.plot(k, meas, 'o', ms=4.2, color=MEAS, mec='white', mew=0.5,
                zorder=3, label='measured')
        ax.plot(k, chy, 's', ms=4.0, mfc='none', mec=CHY, mew=1.1,
                zorder=3, label='chygraph')
        zero, pos = ctrl <= 0, ctrl > 0
        ax.plot(k[zero], np.full(zero.sum(), FLOOR), 'v', ms=3.6, color=CTRL,
                mec='white', mew=0.4, clip_on=False, zorder=3, label='control')
        if pos.any():
            ax.plot(k[pos], ctrl[pos], 'v', ms=4.2, color=CTRL, mec='white',
                    mew=0.5, zorder=3)

        bm, bc = slope(k, meas), slope(k, chy)
        ax.set_xscale('log'); ax.set_yscale('log')
        ax.set_xlim(0.035, 9); ax.set_ylim(FLOOR, 1.6)
        ax.set_xlabel(r'mean degree $\langle k\rangle$')
        ax.set_title(rf'$\tau={tau}$', pad=3)
        # direct labels rather than a colour key: identity must not rest on hue
        ax.text(0.045, 0.90,
                rf'measured  $\beta={bm[0]:.2f}$' + '\n'
                + rf'chygraph  $\beta={bc[0]:.2f}$',
                transform=ax.transAxes, fontsize=6.4, va='top',
                linespacing=1.6, color='0.30')
        for s in ('top', 'right'):
            ax.spines[s].set_visible(False)
        ax.grid(True, which='major', lw=0.35, color='0.88', zorder=0)
        ax.set_axisbelow(True)

    axes[0].set_ylabel('core fraction')
    axes[0].legend(loc='lower right', frameon=False, handletextpad=0.4,
                   borderaxespad=0.3, labelspacing=0.25)
    fig.tight_layout(pad=0.4, w_pad=0.9)
    fig.savefig(OUT, bbox_inches='tight')
    print('wrote', OUT)


def simplicial_panel():
    """Fig. 3: the q=16 double transition, with the free energies that locate it."""
    from chygraph_statmech.simplicial import SimplicialChygraph
    M = SimplicialChygraph([2, 16], [4, 4], [0.7, 0.3])
    lo, hi = M.coexistence(0.96, 0.99, n=200)
    Tc = M.transition(scan=(0.96, 0.99), n=200)

    fig, ax = plt.subplots(1, 2, figsize=(7.05, 2.45))

    Ts = np.linspace(0.955, 1.02, 140)
    for u0, col, mk, lab in ((8.0, MEAS, '-', 'from ordered'),
                             (1e-8, CHY, '--', 'from disordered')):
        m = [M.magnetisation(T, u0=u0)[0] for T in Ts]
        ax[0].plot(Ts, m, mk, color=col, lw=1.4, label=lab)
    ax[0].axvline(Tc, color='0.45', lw=0.8, ls=':')
    ax[0].axvspan(lo, hi, color='0.88', zorder=0)
    ax[0].set_xlabel('$T$'); ax[0].set_ylabel('$m$')
    ax[0].set_title('order parameter', pad=3)
    ax[0].legend(frameon=False, loc='center left', bbox_to_anchor=(0.30, 0.30))
    ax[0].annotate('$T_c$', xy=(Tc, 0.70), xytext=(0.9655, 0.70), fontsize=7,
                   color='0.35', va='center',
                   arrowprops=dict(arrowstyle='-', lw=0.6, color='0.55'))
    ax[0].text(1.0115, 0.05, 'continuous', fontsize=6.5, color='0.35', ha='right')

    Tf = np.linspace(lo + 1e-4, hi - 1e-4, 60)
    d = [M.minus_beta_f(T, u0=8.0) - M.minus_beta_f(T, u0=1e-8) for T in Tf]
    ax[1].axhline(0, color='0.75', lw=0.6)
    ax[1].plot(Tf, np.array(d) * 1e3, color=CTRL, lw=1.4)
    ax[1].axvline(Tc, color='0.45', lw=0.8, ls=':')
    ax[1].set_xlabel('$T$')
    ax[1].set_ylabel(r'$10^{3}\,\Delta(-\beta f)$')
    ax[1].set_title('free-energy difference', pad=3)
    ax[1].text(Tc + 0.0002, min(np.array(d) * 1e3) * 0.75,
               f'$T_c={Tc:.4f}$', fontsize=7, color='0.35')

    for a in ax:
        for sp in ('top', 'right'):
            a.spines[sp].set_visible(False)
        a.grid(True, lw=0.35, color='0.9')
        a.set_axisbelow(True)
    fig.tight_layout(pad=0.4, w_pad=1.2)
    fig.savefig(OUT_SIM, bbox_inches='tight')
    print('wrote', OUT_SIM, f'(T*={lo:.5f}, Tc={Tc:.5f}, T**={hi:.5f})')


if __name__ == '__main__':
    main()
    simplicial_panel()
