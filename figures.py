"""Manuscript figures.

Fig. 2  clustering_panel   T_c against clustering at fixed degree
Fig. 3  simplicial_panel   the q=16 double transition and the free energies
Fig. 4  hittingset_panel   hitting-set density: hard field, soft field, MT
Fig. 5  main               leaf-removal core of hyperbolic random graphs

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
OUT_HS = HERE / 'fig_hittingset.pdf'
HS_JSON = HERE / 'probe' / 'results' / 'hittingset.json'
OUT_SIM = HERE / 'fig_simplicial.pdf'
OUT_CLUST = HERE / 'fig_clustering.pdf'

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
    """Fig. 5: the leaf-removal core of hyperbolic random graphs."""
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
                rf'measured  $\theta={bm[0]:.2f}$' + '\n'
                + rf'chygraph  $\theta={bc[0]:.2f}$',
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


def clustering_panel():
    """Fig. 2: T_c against clustering fraction at genuinely fixed degree.

    A vertex carries k_T = f n / 2 triangles and k_L = (1-f) n links, so its
    degree is exactly n for every f and the whole family shares one p_d.  n is
    an ordered magnitude, so the curves take a sequential ramp rather than
    categorical hues.
    """
    from chygraph_statmech import Chygraph

    def Tc(f, n):
        kT, kL = f * n / 2.0, (1 - f) * n
        cs, ks, ex = [], [], []
        if kL > 1e-12:
            cs.append(2); ks.append(kL); ex.append(kL - 1)
        if kT > 1e-12:
            cs.append(3); ks.append(kT); ex.append(kT - 1)
        return 1.0 / Chygraph(cs, ks, excess=ex).critical_coupling()

    ns = (4, 6, 10, 20)
    ramp = ['#9ecae1', '#6baed6', '#3182bd', '#08519c']       # one hue, light->dark
    fs = np.linspace(0, 1, 41)

    fig, ax = plt.subplots(1, 2, figsize=(7.05, 2.45))

    for n, col in zip(ns, ramp):
        y = np.array([Tc(f, n) for f in fs])
        ax[0].plot(fs, y / y[0], color=col, lw=1.6, label=f'$n={n}$')
        ax[0].text(1.012, y[-1] / y[0], f'{n}', color=col, fontsize=6.5, va='center')
    ax[0].axhline(1.0, color='0.75', lw=0.6, zorder=0)
    ax[0].set_xlabel('fraction $f$ of neighbours inside triangles')
    ax[0].set_ylabel(r'$T_c(f)\,/\,T_c(0)$')
    ax[0].set_title('at fixed degree $d=n$', pad=3)
    ax[0].set_xlim(0, 1.06)
    ax[0].text(1.03, 1.0006, '$n$', fontsize=6.5, color='0.35', va='bottom')

    y4 = np.array([Tc(f, 4) for f in fs])
    ax[1].plot(fs, y4, color=ramp[2], lw=1.6, label='theory, Eq. (8)', zorder=2)
    ax[1].errorbar([0.0, 1.0], [2.894, 2.482], yerr=[0.006, 0.006], fmt='o',
                   ms=5, color=CHY, mec='white', mew=0.6, capsize=2, lw=1.0,
                   label='Monte Carlo', zorder=3)
    ax[1].set_xlabel('fraction $f$ of neighbours inside triangles')
    ax[1].set_ylabel('$T_c$')
    ax[1].set_title('$n=4$, with simulation', pad=3)
    ax[1].legend(frameon=False, loc='upper right')
    ax[1].set_xlim(-0.06, 1.06)

    for a in ax:
        for sp in ('top', 'right'):
            a.spines[sp].set_visible(False)
        a.grid(True, lw=0.35, color='0.9')
        a.set_axisbelow(True)
    fig.tight_layout(pad=0.4, w_pad=1.2)
    fig.savefig(OUT_CLUST, bbox_inches='tight')
    print('wrote', OUT_CLUST,
          f'(n=4: {Tc(0,4):.4f} -> {Tc(1,4):.4f}, {(Tc(1,4)/Tc(0,4)-1)*100:.1f}%)')


def hittingset_panel():
    """Fig. 4: minimum hitting-set density, hard field against soft field.

    Referee-requested: Table III made visible.  Three panels answer three
    separate questions.  Left, how the two treatments separate as cardinality
    grows at fixed chy-degree -- they agree at ``c = 2``, where both reproduce
    the Weigt-Hartmann cover size, and diverge above it.  Centre, the regular
    hypergraphs of Mezard & Tarzia, where the soft field returns their
    ``rho = 1/K`` and the hard field does not, and where the 1RSB answer sits
    between the two.  Right, that the size of the correction tracks the weight
    the ensemble puts on complexes of cardinality three or more.

    Data from ``probe/results/hittingset.json``; regenerate with
    ``python probe/hittingset_density.py``.
    """
    import json
    d = json.load(HS_JSON.open())
    HARD, SOFT, EXACT = CHY, MEAS, '0.35'

    fig, ax = plt.subplots(1, 3, figsize=(7.05, 2.35))

    # -- A: rho against <k> at cardinality 2, 3, 4 -------------------------
    A = d['A']
    marks = {'2': 'o', '3': 's', '4': '^'}
    for c in ('2', '3', '4'):
        k = np.asarray(A[c]['k'])
        ax[0].plot(k, A[c]['hard'], marks[c] + '--', color=HARD, ms=3.6,
                   lw=1.0, mfc='none', mew=1.0)
        ax[0].plot(k, A[c]['soft'], marks[c] + '-', color=SOFT, ms=3.6,
                   lw=1.3, mec='white', mew=0.4)
        ax[0].text(k[-1] * 1.03, A[c]['soft'][-1], f'$c={c}$', fontsize=6.5,
                   color='0.3', va='center')
    ax[0].plot(A['2']['k'], A['WH'], '-', color=EXACT,
               lw=2.6, alpha=0.25, zorder=0)
    ax[0].set_xlabel(r'chy-degree $\langle k\rangle$')
    ax[0].set_ylabel(r'hitting-set density $\rho$')
    ax[0].set_title('Poisson layers', pad=3)
    ax[0].set_xlim(0.35, 3.55)
    ax[0].plot([], [], '--', color=HARD, lw=1.0, label='hard field, Eq. (26)')
    ax[0].plot([], [], '-', color=SOFT, lw=1.3, label='soft field, Eq. (31)')
    ax[0].plot([], [], '-', color=EXACT, lw=2.6, alpha=0.25,
               label='Weigt--Hartmann')
    ax[0].legend(frameon=False, loc='lower right', handlelength=1.6,
                 borderaxespad=0.2, labelspacing=0.2)

    # -- B: regular hypergraphs against Mezard & Tarzia --------------------
    B = d['B']
    xs = np.arange(len(B))
    ax[1].plot(xs, [b['hard'] for b in B], 'o', color=HARD, ms=4.6,
               mfc='none', mew=1.2, label='hard field')
    settled = [b['spread'] < 1e-6 for b in B]
    ax[1].plot(xs[settled], [b['soft'] for b, ok in zip(B, settled) if ok], 's',
               color=SOFT, ms=4.0, mec='white', mew=0.5, label='soft field')
    if not all(settled):                     # drawn hollow: did not settle
        ax[1].plot(xs[[not v for v in settled]],
                   [b['soft'] for b, ok in zip(B, settled) if not ok], 's',
                   color=SOFT, ms=4.0, mfc='none', mew=1.1)
    ax[1].plot(xs, [b['mt'] for b in B], '_', color=EXACT, ms=11, mew=1.4,
               label=r'MT, $\rho=1/K$')
    ax[1].plot([4], [0.178], '*', color=CTRL, ms=8, mec='white', mew=0.4,
               label='MT 1RSB', zorder=4)
    for i, b in enumerate(B):                       # entropy sign: RS validity
        if b['s'] < 0:
            ax[1].plot([i], [b['soft']], 'x', color='0.55', ms=6, mew=1.0,
                       zorder=5)
    ax[1].set_xticks(xs)
    ax[1].set_xticklabels([f"{b['L']},{b['K']}" for b in B], fontsize=6.2)
    ax[1].set_xlabel(r'regular hypergraph $(L,K)$')
    ax[1].set_title('regular, against Ref. [MT07]', pad=3)
    ax[1].set_xlim(-0.6, len(B) - 0.4)
    ax[1].legend(frameon=False, loc='upper right', handletextpad=0.3,
                 borderaxespad=0.2, labelspacing=0.2)
    ax[1].text(0.02, 0.10, r'$\times$: $s<0$, RS invalid', fontsize=6.2,
               color='0.45', transform=ax[1].transAxes)
    ax[1].text(0.02, 0.02, 'hollow: iteration does not settle', fontsize=6.2,
               color='0.45', transform=ax[1].transAxes)

    # -- C: the correction tracks the weight on cardinality >= 3 -----------
    ramp = {'3': '#6baed6', '4': '#2171b5', '6': '#08306b'}
    for c in ('3', '4', '6'):
        x = np.asarray(d['C'][c]['x'])
        y = 100.0 * np.asarray(d['C'][c]['rel'])
        ax[2].plot(x, y, '-o', color=ramp[c], ms=3.2, lw=1.4, mec='white',
                   mew=0.4)
        ax[2].text(1.02, y[-1], f'$c={c}$', fontsize=6.5, color=ramp[c],
                   va='center')
    ax[2].axhline(0, color='0.8', lw=0.6, zorder=0)
    ax[2].set_xlabel('fraction of complexes with $c\\geq3$')
    ax[2].set_ylabel(r'$(\rho_{\rm hard}-\rho_{\rm soft})/\rho_{\rm soft}$  (%)')
    ax[2].set_title(r'at $\langle k\rangle=1$', pad=3)
    ax[2].set_xlim(-0.04, 1.14)

    for a in ax:
        for sp in ('top', 'right'):
            a.spines[sp].set_visible(False)
        a.grid(True, lw=0.35, color='0.9')
        a.set_axisbelow(True)
    fig.tight_layout(pad=0.4, w_pad=1.1)
    fig.savefig(OUT_HS, bbox_inches='tight')
    print('wrote', OUT_HS)


if __name__ == '__main__':
    main()
    clustering_panel()
    simplicial_panel()
    hittingset_panel()
