"""Chapter 13: satisfiability on a chygraph.

Like Chapter 12 this has no manuscript behind it, so the calculations are here.

  fig-sat    the warning-propagation map. At k = 2 the trivial fixed point
             loses stability at alpha = 1, the exact 2-SAT threshold; at k >= 3
             the linearisation vanishes identically and the non-trivial branch
             arrives by a fold instead.
  fig-nested three levels: variables, clauses of variables, and clauses whose
             members are variables *and* clauses. What flattening one to CNF
             costs.

Two of the three results are exact and convention-free, being linearisations:
alpha = 1 at k = 2, and the vanishing of the derivative for k >= 3. The fold is
not: it depends on how a variable receiving contradictory warnings is counted,
and the chapter says so.

Benchmarks, from `~/Downloads/chygraph_references/`:
  Mezard, Parisi & Zecchina, Science 297, 812 (2002)
  Mezard & Zecchina, Phys. Rev. E 66, 056126 (2002)
  Montanari, Ricci-Tersenghi & Semerjian, Table I  (alpha_d, alpha_c, alpha_s)
"""

import itertools
import sys
from fractions import Fraction as F
from pathlib import Path

import numpy as np

OUT = Path(__file__).resolve().parent
DARK, MID, LIGHT = '0.10', '0.45', '0.70'

# Montanari, Ricci-Tersenghi & Semerjian, Table I.
PUBLISHED = {3: (3.86, 3.86, 4.267), 4: (9.38, 9.547, 9.931),
             5: (19.16, 20.80, 21.117), 6: (36.53, 43.08, 43.37)}


def _mpl():
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    return plt


def _tidy(ax):
    ax.tick_params(labelsize=8)
    for sp in ('top', 'right'):
        ax.spines[sp].set_visible(False)


# ---------------------------------------------------------------------------
# the interior sum
# ---------------------------------------------------------------------------

def clause_Z(k, eta, falsify):
    """``Z(sigma_0)`` for one clause, by enumeration.

    ``falsify[j]`` is the value of variable ``j`` that falsifies its literal;
    the clause forbids the single assignment in which every literal is false.
    ``eta[j][s]`` are the incoming messages, ``j = 1 .. k-1``.
    """
    Z = []
    for s0 in (0, 1):
        tot = F(0)
        for rest in itertools.product((0, 1), repeat=k - 1):
            if s0 == falsify[0] and all(r == falsify[j + 1]
                                        for j, r in enumerate(rest)):
                continue                      # the one forbidden assignment
            w = F(1)
            for j, sj in enumerate(rest):
                w *= eta[j][sj]
            tot += w
        Z.append(tot)
    return Z


def check_interior():
    """The interior sum, and the one thing that separates SAT from colouring."""
    for k in (2, 3, 4):
        half = [[F(1, 2), F(1, 2)] for _ in range(k - 1)]
        Z = clause_Z(k, half, [0] * k)
        want = [1 - F(1, 2**(k - 1)), F(1)]
        assert Z == want, (k, Z, want)
        print(f'  k = {k}: at eta = 1/2 the clause emits '
              f'({Z[0]}, {Z[1]}) -- not uniform')
    print('    A clause forbids ONE assignment, and which one depends on its')
    print('    signs, so the uniform message is not a fixed point. Chapter 12\'s')
    print('    constraint was symmetric under permuting colours; this one is')
    print('    not, and that is the whole reason SAT has no closed-form RS')
    print('    threshold where colouring does.')
    # and the shape is Ch. 12's hypergraph rule at q = 2, one configuration
    # forbidden instead of two
    half = [[F(1, 2), F(1, 2)] for _ in range(2)]
    Zs = clause_Z(3, half, [0, 0, 0])
    print(f'    Its shape is Eq. (12.2) at q = 2 with one forbidden '
          f'configuration\n    rather than two: 1 - prod eta_j(u_j), here '
          f'{Zs[0]} against 1')


# ---------------------------------------------------------------------------
# warning propagation
# ---------------------------------------------------------------------------

def pi_forced(x):
    """``P(a cavity variable is forced to one given value)``.

    ``x = <kbar> eta / 2``: the mean number of warnings arriving in one
    direction, the 1/2 being the chance that a warning points the way that
    falsifies the clause being computed.  Strict convention: forced means at
    least one warning that way and none the other.  A variable with warnings
    both ways is a contradiction, not a forcing.
    """
    z = np.exp(-x)
    return (1.0 - z) * z


def wp_map(eta, k, alpha):
    """``eta -> [pi(<kbar> eta / 2)]^(k-1)``, Poisson chy-degree."""
    return pi_forced(k * alpha * np.asarray(eta) / 2.0) ** (k - 1)


def check_two_sat():
    """k = 2: the linearisation gives alpha = 1, the exact 2-SAT threshold."""
    # d(eta')/d(eta) at 0 is <kbar>/2 = k alpha / 2, so k = 2 gives alpha
    for a in (0.5, 0.9, 1.0, 1.1, 2.0):
        slope = 2 * a / 2.0
        num = (wp_map(1e-7, 2, a) / 1e-7)
        assert abs(slope - num) < 1e-5, (a, slope, num)
    print('  k = 2: d(eta\')/d(eta)|_0 = <kbar>/2 = alpha, so the trivial fixed')
    print('    point loses stability at alpha = 1 -- the exact 2-SAT threshold')
    e = 1e-9
    for a in (0.90, 0.99, 1.01, 1.10, 1.50):
        v = branch(2, a)
        print(f'    alpha = {a:.2f}: non-trivial branch eta* = {v:.6f}')
    assert branch(2, 0.99) == 0.0 and branch(2, 1.01) > 0
    print('    and it grows continuously out of alpha = 1   OK')


def check_higher_k():
    """k >= 3: the derivative vanishes identically, at every density.

    Tested the way a vanishing derivative has to be: the finite difference at
    ``eps`` must itself go to zero with ``eps``, and at the rate ``eps^(k-2)``
    that ``(lambda eta / 2)^(k-1)`` predicts.
    """
    for k in (3, 4, 5, 6):
        for a in (1.0, 10.0, 100.0):
            r = [wp_map(e, k, a) / e for e in (1e-6, 5e-7)]
            # the finite difference halves as eps^(k-2): the ratio is the test,
            # an absolute bound would only say eps was chosen small enough
            assert abs(r[0] / r[1] / 2.0**(k - 2) - 1) < 1e-2, (k, a, r)
        print(f'  k = {k}: the finite difference halves as eps^{k-2} at '
              f'alpha = 1, 10, 100')
    print('    so the derivative at the trivial fixed point is exactly zero,')
    print('    and that point is stable at EVERY clause density. A single')
    print('    warning cannot start anything, because a clause needs k-1 of')
    print('    them at once -- Sec. 6.4\'s exclusion, exactly.')


def pi_net(x, nmax=140):
    """The other convention: forced means the *net* warning count points that
    way, so contradictory warnings partly cancel instead of being a
    contradiction.  ``P(N_+ > N_-)`` with both Poisson(x)."""
    from scipy.stats import poisson
    x = np.atleast_1d(np.asarray(x, float))
    n = np.arange(nmax)
    p = poisson.pmf(n[None, :], x[:, None])
    c = np.cumsum(p, axis=1)
    return np.sum(p[:, 1:] * c[:, :-1], axis=1)


def _fold_with(k, pi, lo=1e-3, hi=400.0, iters=45, n=20001):
    e = np.linspace(1e-10, 1.0, n)
    for _ in range(iters):
        mid = 0.5 * (lo + hi)
        g = pi(k * mid * e / 2.0) ** (k - 1) - e
        if np.any(np.diff(np.sign(g)) != 0):
            hi = mid
        else:
            lo = mid
    return hi


def check_conventions():
    """Which results survive the convention and which do not."""
    print('     k   strict   net field   alpha_s')
    for k in (2, 3, 4, 5):
        a = _fold_with(k, pi_forced)
        b = _fold_with(k, pi_net)
        a_s = PUBLISHED[k][2] if k in PUBLISHED else 1.0
        print(f'  {k:>4}   {a:>6.3f}   {b:>9.3f}   {a_s:>7.3f}')
        if k == 2:
            assert abs(a - 1) < 2e-3 and abs(b - 1) < 2e-3
        else:
            assert a > a_s and b > a_s, (k, a, b, a_s)
    print('    At k = 2 the two agree, because the answer is a linearisation.')
    print('    Above it they do not -- but under either convention the fold')
    print('    lies above alpha_s, so the direction is robust and the number')
    print('    is not.')


def branch(k, alpha, n=200001):
    """The largest non-trivial fixed point of the warning map, or 0."""
    e = np.linspace(1e-12, 1.0, n)
    g = wp_map(e, k, alpha) - e
    s = np.where(np.diff(np.sign(g)) != 0)[0]
    return float(e[s[-1]]) if len(s) else 0.0


def fold(k, lo=1e-3, hi=400.0, iters=60):
    """Smallest clause density carrying a non-trivial warning branch."""
    for _ in range(iters):
        mid = 0.5 * (lo + hi)
        if branch(k, mid, 40001) > 1e-8:
            hi = mid
        else:
            lo = mid
    return hi


def check_fold_against_alpha_s():
    """Where the hard field first notices anything, against the truth."""
    print('     k   WP fold   alpha_s (published)   ratio')
    for k in (3, 4, 5):
        f, a_s = fold(k), PUBLISHED[k][2]
        assert f > a_s, (k, f, a_s)
        print(f'  {k:>4}   {f:>7.3f}   {a_s:>19.3f}   {f/a_s:>5.2f}')
    print('    The warning map acquires a branch only well ABOVE the density at')
    print('    which the formula stops being satisfiable, and the gap widens')
    print('    with k. Hard fields are blind here, exactly as in Sec. 10.3 --')
    print('    and this fold is convention-dependent, unlike the two results')
    print('    above; see the docstring of pi_forced.')


def figure_sat():
    plt = _mpl()
    fig, axes = plt.subplots(1, 2, figsize=(4.6, 2.5))

    ax = axes[0]
    e = np.linspace(0, 0.35, 300)
    ax.plot(e, e, ':', lw=0.9, color='0.6')
    for a, col, ls in ((0.8, LIGHT, '-'), (1.0, MID, '-'), (1.3, DARK, '-')):
        ax.plot(e, wp_map(e, 2, a), ls, lw=1.4, color=col,
                label=rf'$\alpha={a:g}$')
    ax.set_title(r'$k=2$', fontsize=8.5)
    ax.annotate('slope $\\alpha$\nat the origin', xy=(0.135, 0.045),
                fontsize=6.4, color='0.35')
    ax.set_xlabel(r'$\eta$', fontsize=8.5)
    ax.set_ylabel(r'warning map', fontsize=8.5)
    ax.legend(frameon=False, fontsize=7, loc='upper left')
    ax.set_xlim(0, 0.35)
    ax.set_ylim(0, 0.35)
    _tidy(ax)

    ax = axes[1]
    e = np.linspace(0, 0.12, 300)
    ax.plot(e, e, ':', lw=0.9, color='0.6')
    f3 = fold(3)
    for a, col, lab in ((4.267, LIGHT, r'$\alpha_{s}=4.27$'),
                        (f3, MID, rf'$\alpha={f3:.2f}$, the fold'),
                        (8.0, DARK, r'$\alpha=8$')):
        ax.plot(e, wp_map(e, 3, a), '-', lw=1.4, color=col, label=lab)
    ax.set_title(r'$k=3$', fontsize=8.5)
    ax.set_xlabel(r'$\eta$', fontsize=8.5)
    ax.legend(frameon=False, fontsize=7, loc='upper left')
    ax.set_xlim(0, 0.12)
    ax.set_ylim(0, 0.12)
    _tidy(ax)

    fig.tight_layout()
    fig.savefig(OUT / 'fig-sat.pdf')
    print(f'  k = 3 fold at alpha = {f3:.3f}')
    print(f'  wrote {OUT / "fig-sat.pdf"}')


# ---------------------------------------------------------------------------
# three levels: a clause whose members are clauses
# ---------------------------------------------------------------------------
# Layer 0  variables.
# Layer 1  C = OR of k1 literals -- a sub-expression, not a constraint.
# Layer 2  D = (NOT C) OR (k2 variable literals), i.e. "if C holds then one of
#          these k2 literals must".  The formula is the conjunction of the
#          layer-2 clauses.  Each C is the antecedent of exactly one D, so
#          alpha_1 = alpha_2 = alpha and a sub-clause is forced false exactly
#          when its single parent says so.
#
# Warning propagation carries two kinds of warning where the flat problem has
# one: eta, a layer-2 clause forcing a plain variable member, and mu, a layer-2
# clause forcing its *sub-clause* member to be false -- which forces every one
# of that sub-clause's k1 variables at once.  The system closes on pi.


def nested_G(pi, k1, k2, alpha):
    """One step of the three-level warning map, reduced to a scalar in ``pi``."""
    pi = np.asarray(pi, float)
    mu = delta = pi**k2                  # D warns its sub-clause; one parent
    gamma = 1.0 - (1.0 - pi)**k1         # some member of C forced to satisfy it
    eta = pi**(k2 - 1) * gamma           # D warns a plain variable member
    lam = k2 * alpha * eta + k1 * alpha * delta
    z = np.exp(-lam / 2.0)
    return (1.0 - z) * z


def _branch(f, n=400001, **kw):
    p = np.linspace(1e-12, 1.0, n)
    g = f(p, **kw) - p
    s = np.where(np.diff(np.sign(g)) != 0)[0]
    return float(p[s[-1]]) if len(s) else 0.0


def _flat_G(pi, k, alpha):
    """Flat k-SAT, in the same variable, for comparison."""
    pi = np.asarray(pi, float)
    lam = k * alpha * pi**(k - 1)
    z = np.exp(-lam / 2.0)
    return (1.0 - z) * z


def _fold(f, key, lo=1e-4, hi=800.0, iters=60, n=60001, **kw):
    for _ in range(iters):
        m = 0.5 * (lo + hi)
        if _branch(f, n=n, **{**kw, key: m}) > 1e-8:
            hi = m
        else:
            lo = m
    return hi


def check_nested_reduces():
    """k1 = 1 makes the middle layer a relay: the three-level system must be
    the two-level one at k = k2 + 1, exactly."""
    for k2 in (2, 3, 4):
        for a in (4.0, 10.0, 30.0, 80.0):
            n = _branch(nested_G, k1=1, k2=k2, alpha=a)
            f = _branch(_flat_G, k=k2 + 1, alpha=a)
            assert n == f, (k2, a, n, f)
        print(f'  k2 = {k2}: identical to flat {k2+1}-SAT at every density tried')
    print('    so the middle layer costs nothing when it holds nothing')


def check_nested_flatten():
    """What the CNF flattening costs, and where it starts costing it.

    Distributing NOT C over the k2 plain literals gives k1 clauses of length
    k2 + 1 that pairwise share those k2 literals.  At k2 = 1 they share one
    variable and the factor graph is still a tree; at k2 >= 2 they share two or
    more, which is exactly Sec. 14.1's condition for a chygraph to stop being
    treelike.  The thresholds agree in the first case and not in the second.
    """
    print('  k2 = 1: the distributed clauses share one variable, still treelike')
    for k1 in (2, 3, 4):
        # both are linear at the origin; the condition is k1 alpha = 1
        for d, want in ((0.99, 0.0), (1.01, None)):
            a = d / k1
            n = _branch(nested_G, k1=k1, k2=1, alpha=a)
            f = _branch(_flat_G, k=2, alpha=k1 * a)
            assert (n > 0) == (f > 0), (k1, d, n, f)
        print(f'    k1 = {k1}: both thresholds at alpha = 1/k1 = {1/k1:.4f}')
    print('  k2 >= 2: they share k2 >= 2 variables, and the folds part company')
    print('     k1  k2   nested    flattened    flat/nested')
    for k1 in (2, 3, 4):
        for k2 in (2, 3):
            fn = _fold(nested_G, 'alpha', k1=k1, k2=k2)
            ff = _fold(_flat_G, 'alpha', k=k2 + 1) / k1
            assert ff < fn, (k1, k2, fn, ff)
            print(f'    {k1:>3} {k2:>3}   {fn:>7.4f}   {ff:>9.4f}   {ff/fn:>11.3f}')
    print('    flattening moves the fold down by 7 to 22 per cent, and the gap')
    print('    grows with k1 -- the size of the sub-clause, which is how much')
    print('    correlation the flattening throws away')


def figure_nested():
    plt = _mpl()
    fig, axes = plt.subplots(1, 2, figsize=(4.6, 2.5))

    ax = axes[0]
    al = np.linspace(0.05, 6.0, 200)
    for k1, col in ((2, LIGHT), (3, MID), (4, DARK)):
        ax.plot(al, [_branch(nested_G, n=40001, k1=k1, k2=2, alpha=a)
                     for a in al], '-', lw=1.4, color=col, label=f'$k_1={k1}$')
        ax.plot(al, [_branch(_flat_G, n=40001, k=3, alpha=k1 * a) for a in al],
                '--', lw=1.0, color=col, dashes=(3, 2))
    ax.set_xlabel(r'$\alpha$   ($k_2=2$)', fontsize=8.5)
    ax.set_ylabel(r'forced fraction $\pi$', fontsize=8.5)
    ax.legend(frameon=False, fontsize=7, loc='center right')
    ax.annotate('dashed:\nflattened to CNF', xy=(0.15, 0.185), fontsize=6.4,
                color='0.35')
    _tidy(ax)

    ax = axes[1]
    ks = (2, 3, 4, 5)
    for k2, col, mk in ((1, LIGHT, '^'), (2, MID, 's'), (3, DARK, 'o')):
        r = []
        for k1 in ks:
            if k2 == 1:
                r.append(1.0)          # thresholds identical, ratio exactly one
            else:
                fn = _fold(nested_G, 'alpha', k1=k1, k2=k2)
                r.append(_fold(_flat_G, 'alpha', k=k2 + 1) / k1 / fn)
        ax.plot(ks, r, mk + '-', ms=3.6, lw=1.3, color=col, label=f'$k_2={k2}$')
    ax.axhline(1.0, ls=':', lw=0.8, color='0.6')
    ax.set_xticks(ks)
    ax.set_xlabel(r'sub-clause size $k_1$', fontsize=8.5)
    ax.set_ylabel('flattened / nested', fontsize=8.5)
    ax.legend(frameon=False, fontsize=7, loc='lower left')
    _tidy(ax)

    fig.tight_layout()
    fig.savefig(OUT / 'fig-nested.pdf')
    print(f'  wrote {OUT / "fig-nested.pdf"}')


if __name__ == '__main__':
    print('the interior sum:')
    check_interior()
    print('two-SAT:')
    check_two_sat()
    print('k >= 3:')
    check_higher_k()
    print('the fold, against the published satisfiability threshold:')
    check_fold_against_alpha_s()
    print('how much of that depends on the convention:')
    check_conventions()
    print('three levels, reduction check:')
    check_nested_reduces()
    print('three levels, what flattening costs:')
    check_nested_flatten()
    print('figures:')
    figure_sat()
    figure_nested()
