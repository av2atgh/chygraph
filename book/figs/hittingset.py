"""Chapter 10: minimum hitting set, hard fields against soft fields.

  fig-hitting   the repair: minimum hitting-set density from the hard-field
                limit and from belief propagation that keeps the O(1) fields,
                for Poisson layers and for regular hypergraphs

Chapter 10's other figure, the counterexample, is TikZ in `hittingset.tex`.

Every soft-field number here comes from population dynamics at mu = 60; the
populations are small enough to run in seconds and are checked against the
closed forms of Mezard & Tarzia wherever those exist.
"""

import sys
from pathlib import Path

import numpy as np
from scipy.special import lambertw

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'percolation' / 'src'))
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'statmech' / 'src'))
import statmech.hittingset as hs  # noqa: E402
from statmech.softfield import (HittingSetBP, regular_density,  # noqa: E402
                                         regular_entropy, regular_field)

OUT = Path(__file__).resolve().parent
DARK, MID, LIGHT = '0.10', '0.45', '0.70'
MU = 60.0


def _mpl():
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    return plt


def _tidy(ax):
    ax.tick_params(labelsize=8)
    for sp in ('top', 'right'):
        ax.spines[sp].set_visible(False)


def soft(cards, degs, regular=False, size=20000, sweeps=300, seed=1):
    return HittingSetBP(cards, degs, regular=regular, mu=MU, size=size,
                        seed=seed, damping=0.5).run(sweeps)


def soft_density(cards, degs, regular=False, size=60000, sweeps=600, seeds=3):
    """Mean density over several seeds, with the spread across them.

    A population of finite size leaves a scatter of a few parts in a thousand
    that does not shrink cleanly with the population, so a single run is quoted
    to three decimals at most and the spread is reported alongside.
    """
    d = [soft(cards, degs, regular, size, sweeps, seed=s).density()
         for s in range(1, seeds + 1)]
    return float(np.mean(d)), float(np.std(d))


def weigt_hartmann(k):
    """The Poisson-graph cover size, Weigt & Hartmann."""
    W = float(np.real(lambertw(k)))
    return 1 - (2 * W + W**2) / (2 * k)


# ------------------------------------------------- (1) the hard-field threshold
def check_threshold():
    """<k>(c-1) = e: counted in neighbours, cardinality does not move it."""
    for c in (2, 3, 4, 5, 10, 20):
        k = hs.rsb_point([c], [1.0])
        assert abs(k * (c - 1) - np.e) < 1e-11, (c, k)
    print('  <k>(c-1) = e to eleven digits for c = 2 ... 20   OK')
    print(f'  c = 2 gives <k> = {hs.rsb_point([2], [1.0]):.6f}, the '
          f'Bauer-Golinelli value (Ch. 11)')


# ------------------------------------------------------- (2) the counterexample
def check_counterexample():
    """Disjoint 3-hyperedges: replica symmetry is trivially valid and the
    hard-field rule is still wrong."""
    x = hs.layer_symbols(1)
    hard = hs.HittingSet([3], x[0])          # chy-degree exactly one
    ch = hard.cover_size()
    cs = soft([3], [1], regular=True).density()
    assert abs(ch - 0.5) < 1e-12 and abs(cs - 1 / 3) < 2e-3, (ch, cs)
    print(f'  disjoint 3-hyperedges: hard {ch:.4f}, soft {cs:.4f}, '
          f'truth {1/3:.4f}')
    print('    the bracket is (0, 1), so the hard field cannot even '
          'distinguish them')


def check_weigt_hartmann():
    """At cardinality two hard and soft must agree, and do."""
    print('     <k>   Weigt-Hartmann      hard      soft')
    for k in (0.5, 1.0, 2.0, 2.5):
        wh = weigt_hartmann(k)
        h = hs.poisson([2], [k]).cover_size()
        s, err = soft_density([2], [k])
        assert abs(h - wh) < 1e-6 and abs(s - wh) < 6e-3, (k, h, s)
        print(f'  {k:>6.1f}  {wh:>14.6f}  {h:>8.6f}  {s:>7.4f} +- {err:.4f}')


# ---------------------------------------------------- (3) the regular solution
def check_regular():
    """h_RS carries a fraction of mu; rho = 1/K; the entropy closed form."""
    print('    L   K    h_RS closed   h_RS population    1/K     rho     s')
    for L, K in ((1, 3), (2, 6), (4, 6), (3, 4), (2, 3)):
        m = soft([K], [L], regular=True)
        hclosed, hpop = regular_field(L, K, MU), float(m.P[0].mean())
        assert abs(hclosed - hpop) < 1e-3, (L, K)
        assert abs(m.density() - regular_density(K)) < 2e-3
        assert abs(m.entropy() - regular_entropy(L, K)) < 1e-3
        print(f'  {L:>3} {K:>3}   {hclosed:>11.5f}   {hpop:>13.5f}   '
              f'{1/K:>6.4f}  {m.density():>6.4f}  {regular_entropy(L, K):>7.4f}')
    print(f'  h_RS = -mu/L - (L-1)/L ln(K-1): a *fraction* of mu, which is '
          f'exactly\n  what an integer ansatz cannot carry')


# ------------------------------------------- (4) mixed cardinality, correlation
def check_mixed():
    """Spread in cardinality postpones the hard-field point; correlation
    across cardinalities brings it forward."""
    print('  spread, at fixed <cbar> = 2 (neighbours per complex):')
    for label, w, cs in (('all c = 3', [1.0], [3]),
                         ('half 2, half 4', [.5, .5], [2, 4]),
                         ('3/4 of 2, 1/4 of 6', [.75, .25], [2, 6])):
        cb = sum(wi * (c - 1) for wi, c in zip(w, cs))
        k = hs.rsb_point(cs, w)
        print(f'    {label:<20} <k> = {k:.5f},  neighbours per node '
              f'{k*cb:.4f}')
    # Two equal-weight node classes with identical marginals, correlated
    # across layers.  `isolated` is Phi(0,...,0): where it grows the shift is
    # dilution rather than correlation, so only the small-spread entries of
    # the positive family are a clean measurement.
    CS, M = [2, 6], [0.75, 0.25]
    cbar = sum(m * (c - 1) for m, c in zip(M, CS)) / sum(M)

    def matched(sign, spread):
        def f(t):
            if sign > 0:
                a_ = [t * m * (1 + spread) for m in M]
                b_ = [t * m * (1 - spread) for m in M]
            else:
                a_ = [t * M[0] * (1 + spread), t * M[1] * (1 - spread)]
                b_ = [t * M[0] * (1 - spread), t * M[1] * (1 + spread)]
            return hs.HittingSet(CS, hs.two_class_phi(a_, b_, 0.5))
        return f

    print('  correlation across layers, cardinalities [2, 6], marginals '
          '[3/4, 1/4]:')
    print('    spread   positive  isolated   negative  isolated')
    for spread in (0.0, 0.4, 0.8, 0.95):
        row = f'    {spread:>6.2f}'
        for sign in (+1, -1):
            f = matched(sign, spread)
            t = hs.rsb_scale(f)
            row += f'  {t*cbar:>9.4f}  {f(t).isolated_fraction():>7.3f}'
        print(row)
    print('    negative correlation moves the point forward monotonically, '
          '6.83 to 2.64;')
    print('    the positive family is clean only at small spread, since its '
          'low class')
    print('    empties into isolated vertices and that is dilution, not '
          'correlation')


# ----------------------------------------- (5) the five structures of Fig. 2.5
FIVE = [
    ('A graph',                 [2]),
    ('A hypergraph',            [3]),
    ('A multiplex network',     [2, 2]),
    ('Interacting graphs',      [2, 2, 2]),
    ('Interacting hypergraphs', [3, 3, 3]),
    ('A clustered graph',       [2, 3]),
]


def _split(cards, nu):
    """Layer chy-degrees splitting ``nu`` neighbours per node evenly.

    A layer-``l`` complex supplies ``c_l - 1`` neighbours per member, so
    ``nu = sum_l kappa_l (c_l - 1)``.  Holding that fixed is what makes the
    rows comparable: the same density of constraints, arranged into complexes
    of different sizes.
    """
    return [nu / len(cards) / (c - 1) for c in cards]


def _instability_closed_form(cards, direction):
    """Neighbours per node at the Eq. (10.4) instability, in closed form.

    With Poisson layers and ``y_l = kappa_l sigma^{c_l-1}``, the fixed point is
    ``sigma = exp(-Y)`` with ``Y = sum_l y_l``, and ``|F'| = 1`` reads
    ``sum_l (c_l-1) y_l = 1``.  Substituting ``sigma`` leaves the pair

        sum_l kappa_l exp(-(c_l-1) Y) = Y,
        sum_l (c_l-1) kappa_l exp(-(c_l-1) Y) = 1,

    in the scale of ``kappa`` and in ``Y``.  When every ``c_l`` is the same
    ``c`` the second sum is ``(c-1)`` times the first, so dividing gives
    ``(c-1) Y = 1`` and hence ``(sum_l kappa_l)(c-1) = e`` -- Eq. (10.4) with
    the *total* chy-degree, at any number of layers.  A spread in cardinality
    is what leaves nothing to divide.
    """
    from scipy.optimize import brentq
    c = np.asarray(cards, float)
    d = np.asarray(direction, float)

    def gap(Y):
        # dividing the two sums removes the scale of kappa and leaves one
        # equation in Y; at a single cardinality it is (c-1) - 1/Y.
        e = np.exp(-(c - 1) * Y)
        return float(((c - 1) * d * e).sum() / (d * e).sum()) - 1.0 / Y

    Y = brentq(gap, 1e-9, 50.0, xtol=1e-14)
    e = np.exp(-(c - 1) * Y)
    t = Y / float((d * e).sum())               # the scale kappa = t * direction
    return float((t * d * (c - 1)).sum())


def check_five_instabilities():
    """One cardinality puts the point at e; only a spread moves it."""
    print('  instability in neighbours per node, closed form against rsb_point:')
    for label, cards in FIVE:
        kap = _split(cards, 1.0)
        closed = _instability_closed_form(cards, kap)
        t = hs.rsb_point(cards, kap)
        num = t * sum(k * (c - 1) for k, c in zip(kap, cards))
        assert abs(closed - num) < 1e-8, (label, closed, num)
        flag = '' if len(set(cards)) > 1 else '  = e'
        print(f'    {label:24s} {closed:.10f}{flag}')
    print(f'    e = {np.e:.10f}; only the mixed-cardinality row moves')


def table_five_structures(nu=2.0):
    """Table 10.1: Chapter 10's two treatments on the structures of Fig. 2.5.

    Hard field against soft field at one density of constraints, with the
    instability point of each structure beside them.  The soft-field column is
    population dynamics and carries the scatter of Sec. 10.7's partial results,
    so it is quoted to three decimals.
    """
    out = [r'\begin{tabular}{lccrrc}', r'\hline\hline',
           r'structure & $c$ & instability & hard & soft & hard exact\\',
           r'\hline']
    print(f'  at {nu} neighbours per node')
    for label, cards in FIVE:
        kap = _split(cards, nu)
        m = hs.HittingSet(cards, hs.poisson_phi(kap))
        sigma = m.solve()
        hard = m.cover_size(sigma)
        lo, hi = m.cover_bracket(sigma)
        soft, spread = soft_density(cards, kap)
        ok = m.certified()
        nu_c = _instability_closed_form(cards, kap)
        out.append(f"{label} & ${','.join(map(str, cards))}$ & {nu_c:.3f} & "
                   f"{hard:.4f} & {soft:.3f} & {'yes' if ok else 'no'}\\\\")
        print(f'    {label:24s} nu_c={nu_c:.4f}  hard={hard:.4f} '
              f'[{lo:.4f},{hi:.4f}]  soft={soft:.3f}+-{spread:.3f}  '
              f"exact={'yes' if ok else 'no'}")
        if ok:
            assert abs(hard - soft) < 5e-3, (label, hard, soft)
    out += [r'\hline\hline', r'\end{tabular}']
    (OUT / 'tab-five-structures-hs.tex').write_text('\n'.join(out) + '\n')
    print('  the cardinality-two rows agree hard with soft; the others do not')
    print('  wrote tab-five-structures-hs.tex')


# ------------------------------------------------------------ (5) the figure
def figure_density():
    plt = _mpl()
    fig, axes = plt.subplots(1, 2, figsize=(4.6, 2.5))

    ax = axes[0]
    ks = np.linspace(0.3, 3.0, 14)
    ax.plot(ks, [weigt_hartmann(k) for k in ks], '-', lw=3.2, color='0.86')
    for c, col in ((2, LIGHT), (3, MID), (4, DARK)):
        ax.plot(ks, [hs.poisson([c], [k]).cover_size() for k in ks], '--',
                lw=1.2, color=col, dashes=(3, 2))
        ax.plot(ks, [soft_density([c], [k], size=30000, sweeps=400, seeds=2)[0]
                     for k in ks], '-', lw=1.4, color=col, label=f'$c={c}$')
    ax.text(1.5, 0.455, 'Weigt\u2013Hartmann', fontsize=6.4, color='0.55')
    ax.set_xlabel(r'chy-degree $\langle k\rangle$', fontsize=8.5)
    ax.set_ylabel('hitting set density', fontsize=8.5)
    ax.legend(frameon=False, fontsize=7.5, loc='lower right')
    _tidy(ax)

    ax = axes[1]
    Ks = np.arange(2, 9)
    L = 4
    x = hs.layer_symbols(1)
    hard = [hs.HittingSet([int(K)], x[0]**L).cover_size() for K in Ks]
    ax.plot(Ks, [1.0 / K for K in Ks], '-', lw=1.5, color=DARK,
            label=r'soft, $=1/K$')
    ax.plot(Ks, hard, '--', lw=1.3, color=MID, dashes=(3, 2), label='hard')
    ax.plot([6], [0.178], '*', ms=9, color=DARK, mfc='white', mew=0.9)
    ax.annotate('1RSB', xy=(6.2, 0.176), fontsize=7, color='0.2')
    ax.annotate(r'$s<0$ throughout,' '\n' r'so soft is a lower bound',
                xy=(4.6, 0.44), fontsize=6.4, color='0.35')
    ax.set_xlabel(r'cardinality $K$   ($L=4$)', fontsize=8.5)
    ax.set_ylim(0.08, 0.60)
    ax.legend(frameon=False, fontsize=7.5, loc='upper right')
    _tidy(ax)

    fig.tight_layout()
    fig.savefig(OUT / 'fig-hitting.pdf')
    print(f'  at L = 4, K = 6:  soft {1/6:.4f} < 1RSB 0.178 < hard '
          f'{hard[4]:.4f}')
    print(f'  wrote {OUT / "fig-hitting.pdf"}')


def check_entropy_is_not_necessary():
    """A positive entropy does not prove replica symmetry holds."""
    m = soft([2], [1.0], size=40000, sweeps=400)
    s, err = m.entropy_averaged()
    print(f'  Poisson graph at <k> = 1: s = {s:+.4f} +- {err:.4f}, positive, '
          f'and yet\n  vertex cover on this ensemble breaks replica symmetry '
          f'at <k> = e   (Ch. 11)')
    print(f'  regular L = 4, K = 6: s = {regular_entropy(4, 6):+.4f}, so the '
          f'RS value 1/6 is an\n  underestimate of the 0.178 of Mezard and '
          f'Tarzia -- and the hard field\'s\n  0.252 is an overestimate, '
          f'further away')
    # the case the chapter quotes for the third signal: the closed form is
    # negative and the iteration does not settle, so no population number is
    # quoted for it and check_regular deliberately excludes it
    s6 = regular_entropy(6, 12)
    assert s6 < 0
    print(f'  regular L = 6, K = 12: s = {s6:+.4f}, and there the iteration does'
          f'\n  not settle at all -- Sec. 1.6\'s third signal, in place of a '
          f'number')


if __name__ == '__main__':
    print('the hard-field threshold:')
    check_threshold()
    print('the counterexample:')
    check_counterexample()
    print('cardinality two:')
    check_weigt_hartmann()
    print('regular hypergraphs:')
    check_regular()
    print('mixed cardinalities and correlation:')
    check_mixed()
    print('the five structures of Fig. 2.5:')
    check_five_instabilities()
    table_five_structures()
    print('the entropy criterion:')
    check_entropy_is_not_necessary()
    print('figure:')
    figure_density()
