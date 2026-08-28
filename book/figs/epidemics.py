"""Chapter 6: epidemics and higher-order spreading.

Two figures, and the numbers quoted in the chapter.

  fig-households   final size of an SIR epidemic with two levels of mixing,
                   households of size three on a Poisson global contact
                   network, computed with `chygraph.household_epidemic`; the
                   marked thresholds are the household reproduction number
                   R* = 1;
  fig-contagion    onset of a contagion in which a group transmits only when
                   all its other members are infected, against the pairwise
                   coupling; one dial, beta_2, carries it from continuous
                   through tricritical to discontinuous.

The second figure is deliberately not a percolation calculation. Sec. 6.4
argues that a threshold rule leaves the class of Chapter 4's map, and the
equation solved here is the simplest self-consistent equation that shows what
happens when it does.
"""

import sys
from pathlib import Path

import numpy as np
from sympy import Rational, exp, nsolve, simplify, symbols

sys.path.insert(0, str(Path.home() / 'av2atg' / 'chygraph' / 'src'))
from chygraph import (clique_excess_pgf, household_epidemic,  # noqa: E402
                      size_biased)

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


# ------------------------------------------------------- (1) two levels of mixing
def figure_households(kbar=2):
    """S against the global transmissibility T, households of size three."""
    plt = _mpl()
    T = symbols('T')
    M = household_epidemic({3: 1})
    print('  theta =', simplify(M.theta()))

    fig, ax = plt.subplots(figsize=(4.3, 2.9))
    for pH, col in ((0.0, LIGHT), (0.5, MID), (1.0, DARK)):
        sub = {'k': kbar, 'p_H': pH}
        Tc = float(nsolve(M.theta().subs(sub), T, 0.3))
        Ts = np.linspace(0.02, 1.0, 99)
        S = [M.node_fraction({**sub, 'T': v}) for v in Ts]
        ax.plot(Ts, S, '-', lw=1.4, color=col,
                label=rf'$p_H={pH:g}$   ($T_c={Tc:.4f}$)')
        ax.axvline(Tc, ls=':', lw=0.8, color=col)
        # mu_H, the mean number of additional household members infected
        gbar = list(size_biased({3: 1}).items())
        y = symbols('y')
        mu = sum(w * clique_excess_pgf(n, pH)(y).diff(y).subs(y, 1)
                 for n, w in gbar)
        B = M.amplitude_numeric({**sub, 'T': Tc}, 0)
        print(f'  p_H = {pH:g}:  mu_H = {float(mu):.4f},  T_c = {Tc:.4f},  '
              f'R*(T=1) = {float(1 + M.theta().subs({**sub, "T": 1})):.4f},  '
              f'B = {B:.4f},  S(T=1) = {S[-1]:.4f}')
    ax.set_xlabel(r'global transmissibility $T$')
    ax.set_ylabel(r'final size $S$')
    ax.legend(frameon=False, fontsize=7.5, loc='lower right')
    _tidy(ax)
    fig.tight_layout()
    fig.savefig(OUT / 'fig-households.pdf')


def check_triangle_identity():
    """The size-three household final size IS Ch. 4's triangle enumeration."""
    q, y = symbols('q y')
    gbar = clique_excess_pgf(3, q)(y).expand()
    tri = ((3 * q**2 - 2 * q**3) * y**2 + 2 * q * (1 - q)**2 * y
           + (1 - q)**3 + q * (1 - q)**2).expand()
    assert simplify(gbar - tri) == 0, (gbar, tri)
    mu = simplify(gbar.diff(y).subs(y, 1))
    assert simplify(mu - 2 * q * (1 + q - q**2)) == 0, mu
    print(f'  Gbar^1_0 for n = 3 equals the bond-percolated triangle; '
          f'mu_H = {mu}')


def check_Rstar():
    """theta + 1 = T [<kbar> + mu_H <k>], symbolically, for several sizes."""
    T, k, pH, y = symbols('T k p_H y')
    for sizes in ({3: 1}, {2: Rational(1, 2), 4: Rational(1, 2)},
                  {1: Rational(1, 3), 5: Rational(2, 3)}):
        M = household_epidemic(sizes)
        mu = sum(w * clique_excess_pgf(n, pH)(y).diff(y).subs(y, 1)
                 for n, w in size_biased(sizes).items())
        assert simplify(M.theta() + 1 - T * (k + mu * k)) == 0, sizes
        print(f'  R* verified for household sizes {dict(sizes)}')


def check_fixed_budget():
    """Households at a *fixed contact budget*, four contacts per person.

    Sec. 6.5's "say what is held fixed": Sec. 6.3's table adds household
    transmission on top of an unchanged global network, and households then
    help the epidemic.  Taking the household contacts out of the same budget
    instead -- size-three households (two contacts each) plus a Poisson global
    network of mean two, against no households and a Poisson network of mean
    four -- reverses the sign.
    """
    from scipy.optimize import brentq
    T = symbols('T')
    M = household_epidemic({3: 1})
    th = M.theta().subs({'k': 2, 'p_H': T})
    Tc_hh = float(nsolve(th, T, 0.2))
    print(f'  fixed budget of four contacts:  T_c = {Tc_hh:.4f} with '
          f'households, {1 / 4:.4f} without')
    for t in (0.3, 0.5):
        S_hh = M.node_fraction({'k': 2, 'p_H': t, 'T': t})
        S_pl = (brentq(lambda x: 1 - np.exp(-4 * t * x) - x, 1e-12, 1 - 1e-12)
                if 4 * t > 1 else 0.0)
        print(f'    T = {t:g}:  S = {S_hh:.4f} with households, '
              f'{S_pl:.4f} without')


# ---------------------------------------------------- (2) contagion inside a group
def branch(beta2, rho):
    """beta_1 along the branch rho = 1 - exp(-beta_1 rho - beta_2 rho^2)."""
    return (-np.log1p(-rho) - beta2 * rho**2) / rho


def saddle_node(beta2):
    """The fold: the minimum over rho of beta_1(rho)."""
    r = np.linspace(1e-6, 1 - 1e-9, 4000001)
    b = branch(beta2, r)
    i = int(np.argmin(b))
    return float(b[i]), float(r[i])


def figure_contagion():
    plt = _mpl()
    fig, ax = plt.subplots(figsize=(4.3, 3.0))
    rho = np.linspace(1e-6, 1 - 1e-9, 20000)
    for beta2, col, lab in ((0.0, LIGHT, r'$\beta_2=0$  (pairwise only)'),
                            (0.5, MID, r'$\beta_2=1/2$  (tricritical)'),
                            (2.0, DARK, r'$\beta_2=2$')):
        b1 = branch(beta2, rho)
        bmin, rmin = saddle_node(beta2)
        stable = rho >= rmin
        ax.plot(b1[stable], rho[stable], '-', lw=1.4, color=col, label=lab)
        if bmin < b1[-1] and rmin > 1e-4:
            ax.plot(b1[~stable], rho[~stable], '--', lw=1.0, color=col,
                    dashes=(3, 2))
            ax.plot([bmin], [rmin], 'o', ms=3.5, color=col)
        print(f'  beta_2 = {beta2:g}:  fold at beta_1 = {bmin:.4f}, '
              f'rho = {rmin:.4f};  invasion threshold beta_1 = 1')
    # the jump at the invasion threshold, for the discontinuous case
    from scipy.optimize import brentq
    up = brentq(lambda x: 1 - np.exp(-x - 2.0 * x * x) - x, 0.5, 1 - 1e-9)
    ax.annotate('', xy=(1.0, up), xytext=(1.0, 0.0),
                arrowprops=dict(arrowstyle='-|>', lw=0.9, color=DARK,
                                shrinkA=0, shrinkB=0))
    ax.text(0.93, up / 2, r'jump', fontsize=7.5, color=DARK, rotation=90,
            va='center', ha='right')
    print(f'  beta_2 = 2: at beta_1 = 1 the outbreak branch is rho = {up:.4f}')
    ax.axvline(1.0, ls=':', lw=0.8, color='0.6')
    ax.text(1.04, 0.03, r'$\beta_1=1$', fontsize=7.5, color='0.35')
    ax.plot([0.0, 1.0], [0.0, 0.0], '-', lw=1.4, color=DARK)
    ax.set_xlim(0.0, 2.6)
    ax.set_ylim(0, 1)
    ax.set_xlabel(r'pairwise coupling $\beta_1=T\langle k\rangle_{|}$')
    ax.set_ylabel(r'final size $\rho$')
    ax.legend(frameon=False, fontsize=7.5, loc='lower right')
    _tidy(ax)
    fig.tight_layout()
    fig.savefig(OUT / 'fig-contagion.pdf')


def check_tricritical():
    """The quadratic coefficient at beta_1 = 1 changes sign at beta_2 = 1/2."""
    r, b1, b2 = symbols('rho beta_1 beta_2')
    f = (1 - exp(-(b1 * r + b2 * r**2))).series(r, 0, 4).removeO()
    c2 = simplify(f.coeff(r, 2).subs(b1, 1))
    c3 = simplify(f.coeff(r, 3).subs(b1, 1))
    print(f'  at beta_1 = 1:  quadratic coefficient {c2}, cubic {c3}')
    assert simplify(c2.subs(b2, Rational(1, 2))) == 0
    # and the fold is at beta_1 < 1 exactly when beta_2 > 1/2
    for b, expect in ((0.4, False), (0.5, False), (0.6, True), (2.0, True)):
        bmin, _ = saddle_node(b)
        assert (bmin < 1 - 1e-6) == expect, (b, bmin)
    print('  fold below beta_1 = 1 exactly for beta_2 > 1/2')


def check_pure_group():
    """With no pairwise route there is no invasion threshold at all."""
    from scipy.optimize import brentq
    # tangency of rho = 1 - exp(-beta_2 rho^2): eliminating the exponential
    # leaves -2 ln(1-rho)(1-rho) = rho, and then beta_2 = 1/[2 rho (1-rho)]
    r = brentq(lambda x: -2 * np.log1p(-x) * (1 - x) - x, 1e-9, 1 - 1e-9)
    b2 = 1 / (2 * r * (1 - r))
    assert abs(b2 + np.log1p(-r) / r**2) < 1e-9
    print(f'  beta_1 = 0: no outbreak for any beta_2 < {b2:.4f}; at that '
          f'point rho jumps to {r:.4f}')
    return b2, r


def check_tricritical_exponent():
    """At beta_2 = 1/2 the branch leaves as sqrt(3 (beta_1 - 1))."""
    from scipy.optimize import brentq
    for d in (1e-2, 1e-3, 1e-4):
        g = lambda x: 1 - np.exp(-(1 + d) * x - 0.5 * x * x) - x  # noqa: E731
        rr = brentq(g, 1e-12, 1 - 1e-9)
        print(f'  beta_1 - 1 = {d:g}:  rho = {rr:.6f},  '
              f'rho / sqrt(3 (beta_1-1)) = {rr / np.sqrt(3 * d):.4f}')


if __name__ == '__main__':
    print('checks:')
    check_triangle_identity()
    check_Rstar()
    check_fixed_budget()
    check_tricritical()
    check_pure_group()
    check_tricritical_exponent()
    print('households:')
    figure_households()
    print('contagion:')
    figure_contagion()
    print('wrote two PDFs')
