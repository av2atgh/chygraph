"""Chapter 9: the Ising model on a chygraph.

  fig-transmit    what a complex transmits per traversal, u'(t) for c = 2, 3, 4
  fig-clustering  the headline: at fixed degree, clustering lowers T_c, with
                  the Monte Carlo of `probe/results/ising_mc.log` on the same
                  axes
  fig-unanimity   the unanimity interaction: the spinodal against cardinality
                  at finite and infinite chy-degree, and the double transition
                  at q = 16 where the spinodal is not the transition

The Monte Carlo points are read from the cached log of `probe/ising_mc.py`
rather than recomputed here; that script uses no cavity equation and is the one
test in the chapter from outside the formalism.

`figure_unanimity` solves the q = 16 coexistence window and takes about a
minute; everything else is immediate.
"""

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path.home() / 'av2atg' / 'statmech' / 'src'))
from statmech import Chygraph, ising  # noqa: E402
from statmech.simplicial import SimplicialChygraph, uprime  # noqa: E402

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


# ------------------------------------------------------- (1) what a clique transmits
def check_uprime():
    """The three closed forms, against the enumeration behind Eq. (8.4)."""
    for bJ in (0.05, 0.2, 0.5, 1.0):
        t = np.tanh(bJ)
        closed = {2: t,
                  3: t / (1 - t + t**2),
                  4: (t + t**3) / (1 - 2 * t + 3 * t**2)}
        for c, want in closed.items():
            got = ising.clique_derivative(c, bJ)
            assert abs(got - want) < 1e-12, (c, bJ, got, want)
    print('  u\' = t,  t/(1-t+t^2),  (t+t^3)/(1-2t+3t^2)  for c = 2, 3, 4   OK')
    # a triangle transmits more per traversal than an edge, at every coupling
    for bJ in np.linspace(0.02, 2.0, 40):
        assert ising.clique_derivative(3, bJ) > ising.clique_derivative(2, bJ)
    print('  u\'(3) > u\'(2) at every coupling   OK')


def figure_transmission():
    plt = _mpl()
    fig, ax = plt.subplots(figsize=(4.3, 2.5))
    bJ = np.linspace(1e-4, 2.6, 300)
    t = np.tanh(bJ)
    for c, col in ((2, LIGHT), (3, MID), (4, DARK)):
        u = [ising.clique_derivative(c, b) for b in bJ]
        ax.plot(t, u, '-', lw=1.4, color=col, label=f'$c={c}$')
    ax.set_xlabel(r'$t=\tanh\beta J$', fontsize=8.5)
    ax.set_ylabel(r"transmission $u'$", fontsize=8.5)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.legend(frameon=False, fontsize=8, loc='lower right')
    _tidy(ax)
    fig.tight_layout()
    fig.savefig(OUT / 'fig-transmit.pdf')
    print(f'  wrote {OUT / "fig-transmit.pdf"}')


# ------------------------------------------------------- (2) clustering lowers T_c
def Tc_family(n, f, regular=True):
    """T_c for a vertex carrying f n/2 triangles and (1-f) n links.

    Degree is exactly n for every f, so the whole family shares one degree
    distribution and only the clustering coefficient, f/(n-1), varies.
    """
    kL, kT = (1 - f) * n, f * n / 2.0
    cards, degs, exc = [], [], []
    if kL > 1e-12:
        cards.append(2)
        degs.append(kL)
        exc.append(kL - 1 if regular else kL)
    if kT > 1e-12:
        cards.append(3)
        degs.append(kT)
        exc.append(kT - 1 if regular else kT)
    return 1.0 / Chygraph(cards, degs, excess=exc).critical_coupling()


def check_clustering():
    """The regular nulls of the chapter's table, and the mechanism."""
    print('   n    T_c links    T_c triangles    change')
    for n in (4, 6, 8, 10, 20):
        TL, TT = Tc_family(n, 0.0), Tc_family(n, 1.0)
        assert TT < TL, n
        print(f'  {n:>2}   {TL:9.4f}    {TT:9.4f}      {100*(TT/TL-1):+.1f}%')
    # the family is monotone in f at every degree
    for n in (4, 6, 10, 20):
        row = [Tc_family(n, f) for f in np.linspace(0, 1, 21)]
        assert all(b < a for a, b in zip(row, row[1:])), n
    print('  T_c falls monotonically with f at every degree   OK')
    # The break-even arithmetic at n = 4, evaluated at the *link* ensemble's
    # critical coupling: that is the comparison, not each ensemble at its own.
    bJ = np.arctanh(1.0 / 3.0)
    u3 = ising.clique_derivative(3, bJ)
    assert abs(u3 - 3 / 7) < 1e-12
    print(f'  at n = 4, both at t = 1/3: links need 3t = 1, triangles need '
          f'2u\' = 1, so u\' must reach 1/2 and reaches {u3:.4f} = 3/7')


def check_misleading_comparison():
    """Sec. 9.4's "comparison that misleads", and why the sign flips.

    A Poisson triangle layer matched to a Poisson link layer on *mean degree*
    is not a matched control.  The triangle construction gives degree d = 2X
    with X ~ Poisson(n/2), so Var[d] = 2n rather than n and the excess degree
    is <kbar> = <d^2>/<d> - 1 = n + 1, not n.

    Eq. (8.7) sees that through the multiplicity rather than through the degree:
    a Poisson layer is its own excess, so the triangle layer carries
    (c-1)<kappabar> = 2 (n/2) = n, exactly the link layer's <kappabar> = n.  The
    branch lost in the regular comparison is not lost here, only the
    transmission enhancement survives, and T_c comes out *higher* for the
    triangles.  Matching the control on the triangle ensemble's own excess
    degree n + 1 restores the sign the regular comparison found.
    """
    print('    n   links (mean-matched)   triangles   links (excess-matched)')
    for n in (4, 6, 8):
        TL = Tc_family(n, 0.0, regular=False)
        TT = Tc_family(n, 1.0, regular=False)
        # a link layer whose excess chy-degree is the triangle ensemble's own
        # graph excess degree, n + 1
        TC = 1.0 / Chygraph([2], [n + 1.0], excess=[n + 1.0]).critical_coupling()
        # the analytic moments of d = 2X, X ~ Poisson(n/2)
        lam = n / 2.0
        d1, d2 = 2 * lam, 4 * (lam + lam * lam)
        assert abs(d2 / d1 - 1 - (n + 1)) < 1e-12, n
        assert abs((d2 - d1 ** 2) - 2 * n) < 1e-12, n
        assert TT > TL and TT < TC, n          # sign reverses, then is restored
        print(f'   {n:>2}   {TL:18.4f}   {TT:9.4f}   {TC:20.4f}')
    print('  mean-degree matching reverses the sign; excess-degree matching '
          'does not   OK')


MC = {'links': (2.894, [(2.600, .6627, .6656), (2.750, .6466, .6610),
                        (2.850, .5604, .6052), (2.890, .4745, .4803),
                        (2.950, .3204, .2333), (3.100, .0947, .0428)]),
      'triangles': (2.482, [(2.200, .6643, .6661), (2.350, .6538, .6631),
                            (2.450, .5662, .6148), (2.490, .4720, .4597),
                            (2.550, .2746, .1346), (2.700, .0372, -.0144)])}


def figure_clustering():
    plt = _mpl()
    fig, axes = plt.subplots(1, 2, figsize=(4.6, 2.5))
    fs = np.linspace(0, 1, 41)

    ax = axes[0]
    for n, col in ((4, DARK), (6, '0.30'), (10, MID), (20, LIGHT)):
        row = np.array([Tc_family(n, f) for f in fs])
        ax.plot(fs, row / row[0], '-', lw=1.4, color=col, label=f'$n={n}$')
    ax.set_xlabel(r'fraction $f$ through triangles', fontsize=8)
    ax.set_ylabel(r'$T_c(f)/T_c(0)$', fontsize=8.5)
    ax.legend(frameon=False, fontsize=7, loc='lower left')
    _tidy(ax)

    ax = axes[1]
    row = np.array([Tc_family(4, f) for f in fs])
    ax.plot(fs, row, '-', lw=1.4, color=DARK)
    for (label, (Tmc, _)), f, mk in (((k, v), ff, m) for (k, v), ff, m
                                     in zip(MC.items(), (0.0, 1.0), ('o', 's'))):
        ax.plot([f], [Tmc], mk, ms=4.5, color=DARK, mfc='white', mew=1.0)
    ax.text(0.06, 2.86, 'Monte Carlo', fontsize=7, color='0.25')
    ax.set_xlabel(r'fraction $f$ through triangles', fontsize=8)
    ax.set_ylabel(r'$T_{c}$ at $n=4$', fontsize=8.5)
    _tidy(ax)

    fig.tight_layout()
    fig.savefig(OUT / 'fig-clustering.pdf')
    for name, (Tmc, _) in MC.items():
        Tth = Tc_family(4, 0.0 if name == 'links' else 1.0)
        print(f'  {name:<10} Monte Carlo {Tmc:.3f}  theory {Tth:.4f}  '
              f'({100*abs(Tmc-Tth)/Tth:.2f}%)')
    print(f'  shift: Monte Carlo {100*(MC["triangles"][0]/MC["links"][0]-1):+.1f}%, '
          f'theory {100*(Tc_family(4,1)/Tc_family(4,0)-1):+.1f}%')
    print(f'  wrote {OUT / "fig-clustering.pdf"}')


# --------------------------------------------------------------- (3) the AT line
def check_at_line():
    """Squaring u' gives the spin-glass boundary; it always needs more coupling."""
    print('   c    beta J_c    AT line')
    for c in range(2, 7):
        g = Chygraph([c], [2.0])
        fm, at = g.critical_coupling(), g.critical_coupling(squared=True)
        assert at > fm, c
        print(f'  {c:>2}    {fm:.4f}      {at:.4f}')
    print("  the spin-glass line always needs the stronger coupling, "
          "since u'^2 < u'   OK")


# ------------------------------------------------------- (4) the unanimity rule
def Tstar(q):
    return q * (q - 1) / 2.0**(q - 1)


def Cq(q):
    return q * (2 * q - 2.0**(q - 1)) / 2.0**q


def check_unanimity():
    """The closed form, the Bragg-Williams limit, and its residual."""
    from sympy import Rational, exp, simplify, symbols
    a, q = symbols('a q')
    for qq in (2, 3, 5):
        for bJ in (0.3, 1.1):
            want = float((exp(bJ) - 1) / (2**(qq - 1) + exp(bJ) - 1))
            assert abs(uprime(qq, bJ) - want) < 1e-12
    print("  u' = (e^a - 1)/(2^{q-1} + e^a - 1)   OK")
    print('   q     T*        C_q       T-T* at k=200    -C_q/200')
    for qq in (2, 3, 4, 5, 6, 8):
        m = SimplicialChygraph([qq], [200], [qq / 200.0])
        d = m.spinodal() - Tstar(qq)
        pred = -Cq(qq) / 200.0
        assert abs(d - pred) < 3e-4, (qq, d, pred)
        print(f'  {qq:>2}   {Tstar(qq):.5f}   {Cq(qq):+.4f}   {d:+.3e}      '
              f'{pred:+.3e}')
    print('  C_q vanishes at q = 4 and only there, so T* is approached as '
          '1/k^2 there')


def _window_edges(q, k, scan=1.5, n=60):
    """The coexistence window itself, for the one case the chapter quotes."""
    m = SimplicialChygraph([q], [k], [1.0])
    Ts = m.spinodal()
    return m.coexistence(Ts * 1.0001, Ts * scan, n=n)


def _window(q, k, scan=1.5, n=60):
    """(relative width of the coexistence window, first-order T_c or nan)."""
    m = SimplicialChygraph([q], [k], [1.0])
    Ts = m.spinodal()
    lo, hi = m.coexistence(Ts * 1.0001, Ts * scan, n=n)
    if not np.isfinite(lo):
        return 0.0, float('nan')
    return (hi - lo) / Ts, m.transition(lo=lo, hi=hi, iters=40)


def check_tricritical():
    """Where the transition stops being continuous.

    The naive test -- iterate from an ordered start just above the spinodal and
    see whether the magnetisation survives -- cannot tell a narrow first-order
    window from critical slowing down, and reports coexistence in both cases.
    The test that can: a coexistence window of non-negligible width AND a
    free-energy crossing strictly inside it.  At q = 4 the apparent window is
    5e-4 of T* and holds no crossing at any chy-degree; at q = 5 it is 1.5% to
    18% and holds one.
    """
    print('    k   q   window/T*    first-order T_c')
    for k in (3, 4, 6, 10, 40):
        for q in ((5, 6) if k == 3 else (4, 5)):
            w, Tc = _window(q, k)
            first = np.isfinite(Tc) and w > 1e-3
            assert first == (q == (6 if k == 3 else 5)), (k, q, w, Tc)
            lo_hi = _window_edges(q, k) if (q == 5 and k == 4) else None
            if lo_hi:
                print(f'       window {lo_hi[0]:.4f} to {lo_hi[1]:.4f}')
            print(f'  {k:>3}  {q:>2}   {w:.5f}      '
                  f'{"none" if not np.isfinite(Tc) else f"{Tc:.4f}"}'
                  f'   {"first order" if first else "continuous"}')
    print('  boundary: between q = 5 and 6 at chy-degree 3, between q = 4 and 5'
          '\n            from chy-degree 4 upward')


def figure_unanimity():
    plt = _mpl()
    fig, axes = plt.subplots(1, 2, figsize=(4.6, 2.5))

    ax = axes[0]
    qs = np.arange(2, 10)
    ax.plot(qs, [Tstar(q) for q in qs], 'o-', ms=3, lw=1.2, color=DARK,
            label=r'$k\to\infty$')
    for k, col in ((6, MID), (4, LIGHT)):
        Ts = [SimplicialChygraph([int(q)], [k], [q / k]).spinodal() for q in qs]
        ax.plot(qs, Ts, 's--', ms=2.6, lw=1.0, color=col, dashes=(3, 2),
                label=f'$k={k}$')
    ax.set_xlabel(r'cardinality $q$', fontsize=8.5)
    ax.set_ylabel(r'spinodal $T^{*}$', fontsize=8.5)
    ax.legend(frameon=False, fontsize=7.5, loc='upper right')
    _tidy(ax)

    ax = axes[1]
    m = SimplicialChygraph([16, 2], [4, 4], [0.3, 0.7])
    lo, hi = m.coexistence(0.95, 1.02)
    Tc = m.transition(lo=lo, hi=hi)
    Tcont = m.spinodal()
    # dense inside the window, where the two branches differ, coarse outside
    Ts = np.unique(np.concatenate([np.linspace(0.970, 1.016, 70),
                                   np.linspace(lo - 5e-4, hi + 5e-4, 60)]))
    ord_ = [m.magnetisation(T, u0=8.0)[0] for T in Ts]
    dis = [m.magnetisation(T, u0=1e-8)[0] for T in Ts]
    ax.axvspan(lo, hi, color='0.86', lw=0)
    ax.plot(Ts, ord_, '-', lw=1.5, color=DARK, label='from ordered')
    ax.plot(Ts, dis, '--', lw=1.2, color=LIGHT, dashes=(3, 2),
            label='from disordered')
    ax.axvline(Tc, ls=':', lw=0.9, color=DARK)
    ax.axvline(Tcont, ls=':', lw=0.9, color='0.55')
    ax.annotate(r'$T_{c}$', xy=(Tc - 0.0015, 0.93), fontsize=7.5, ha='right',
                color=DARK)
    ax.annotate('continuous', xy=(Tcont - 0.0015, 0.30), fontsize=7,
                ha='right', color='0.35')
    ax.set_xlim(0.970, 1.016)
    ax.set_ylim(0, 1.0)
    ax.set_xlabel(r'$T$   ($q=16$)', fontsize=8.5)
    ax.set_ylabel(r'magnetisation $m$', fontsize=8.5)
    ax.legend(frameon=False, fontsize=6.8, loc='lower left')
    _tidy(ax)

    fig.tight_layout()
    fig.savefig(OUT / 'fig-unanimity.pdf')
    mo = m.magnetisation(Tc - 1e-5, u0=8.0)[0]
    md = m.magnetisation(Tc - 1e-5, u0=1e-8)[0]
    jump = mo - md
    print(f'  q = 16 double transition: continuous at T = {m.spinodal():.4f}, '
          f'coexistence {lo:.4f}-{hi:.4f}, first order at T_c = {Tc:.4f}')
    print(f'    m jumps from {md:.3f} to {mo:.3f}, by {jump:.3f}; quoting the spinodal instead of T_c '
          f'would be wrong by {100*abs(lo-Tc)/Tc:.1f}%')
    print(f'  wrote {OUT / "fig-unanimity.pdf"}')


if __name__ == '__main__':
    print('what a clique transmits:')
    check_uprime()
    print('clustering lowers T_c:')
    check_clustering()
    print('the comparison that misleads:')
    check_misleading_comparison()
    print('the de Almeida-Thouless line:')
    check_at_line()
    print('the unanimity interaction:')
    check_unanimity()
    print('where it stops being continuous:')
    check_tricritical()
    print('figures:')
    figure_transmission()
    figure_clustering()
    figure_unanimity()
