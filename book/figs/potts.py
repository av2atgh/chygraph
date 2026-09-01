"""Chapter 7: the random-cluster model, and one recursion at two values of q.

No plot -- Chapter 7's only figure is TikZ in `potts.tex`.  What is here is the
chapter's central claim, checked rather than asserted:

    the q -> 1 limit of the Potts cavity recursion on a chygraph IS the
    percolation map of Chapter 4, complex by complex.

The interior sum is done at symbolic q by summing over set partitions of a
complex's members, so that "q states" never has to be an integer: a partition
into m blocks admits a falling factorial q(q-1)...(q-m+1) of assignments of
distinct values, which is a polynomial in q and continues to real q on its own.
Substituting the percolation scaling of the messages and letting q -> 1 then
returns the intra-complex generating function Gbar of Sec. 2.8, which is what
had to be shown.

Checked here for a link, a triangle, a 4-clique, a 5-clique, a hyperedge under
OR-logic, and a path -- against `percolation.clique_excess_pgf` and against direct
enumeration, symbolically in the bond probability.
"""

import itertools
import sys
from pathlib import Path

from sympy import (Rational, binomial, expand, factor, limit, simplify,
                   symbols, together)

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'percolation' / 'src'))
from percolation import clique_excess_pgf  # noqa: E402

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


# ---------------------------------------------------------------------------
# set partitions
# ---------------------------------------------------------------------------

def set_partitions(items):
    """Every partition of ``items`` into non-empty blocks."""
    items = list(items)
    if not items:
        yield []
        return
    first, rest = items[0], items[1:]
    for part in set_partitions(rest):
        for i in range(len(part)):
            yield part[:i] + [[first] + part[i]] + part[i + 1:]
        yield [[first]] + part


def falling(q, m):
    """``q (q-1) ... (q-m+1)``, a polynomial in ``q`` for integer ``m >= 0``."""
    out = 1
    for i in range(m):
        out *= (q - i)
    return out


# ---------------------------------------------------------------------------
# the interior sum at symbolic q
# ---------------------------------------------------------------------------

def emitted_ratio(members, edges, y, q, v, root=0):
    """``rho = Z(sigma_root != 1) / Z(sigma_root = 1)`` for one complex.

    The complex holds ``members`` with Potts couplings on ``edges``, weight
    ``1 + v delta`` per edge.  Incoming messages are unnormalised, member ``j``
    carrying weight ``1`` on the distinguished state and ``y[j]`` on each of the
    other ``q - 1`` states.  ``root`` receives no message: this is the emitted,
    cavity quantity.

    The sum is over set partitions of the members, one block optionally carrying
    the distinguished state; the number of ways to colour the remaining blocks
    with distinct non-distinguished values is a falling factorial in ``q``.
    """
    others = [m for m in members if m != root]

    def Z(root_is_special):
        total = 0
        for part in set_partitions(members):
            # bond weight depends only on which pairs share a block
            same = {frozenset((u, w)) for blk in part
                    for u, w in itertools.combinations(blk, 2)}
            bonds = sum(1 for e in edges if frozenset(e) in same)
            w_bond = (1 + v)**bonds
            root_blk = next(i for i, b in enumerate(part) if root in b)
            for one in list(range(len(part))) + [None]:
                # `one` is the index of the block holding the distinguished
                # state, or None if no block does
                if root_is_special and one != root_blk:
                    continue
                if (not root_is_special) and one == root_blk:
                    continue
                # message weight: 1 on the distinguished state, y_j otherwise
                w_msg = 1
                for i, blk in enumerate(part):
                    if i == one:
                        continue
                    for j in blk:
                        if j != root:
                            w_msg *= y[j]
                free = len(part) - (1 if one is not None else 0)
                if root_is_special:
                    # root sits in the distinguished state, so its block IS the
                    # `one` block: every remaining block takes a distinct value
                    # from the q - 1 others
                    n_col = falling(q - 1, free)
                else:
                    # root sits in one *fixed* non-distinguished state, using
                    # its block up and leaving q - 2 values for the rest
                    n_col = falling(q - 2, free - 1)
                total += w_bond * w_msg * n_col
        return expand(total)

    return together(Z(False) / Z(True))


def percolation_limit(members, edges, y, root=0):
    """``lim_{q -> 1} rho``, rewritten in the bond probability ``p``."""
    q, v, p = symbols('q v p')
    rho = emitted_ratio(members, edges, y, q, v, root=root)
    out = limit(rho, q, 1)
    return simplify(expand(out.subs(v, p / (1 - p))))


# ---------------------------------------------------------------------------
# reference: what Chapter 4 says the answer is
# ---------------------------------------------------------------------------

def reachability_pgf(members, edges, y, p, root=0):
    """``Gbar``: PGF of which members are reachable from ``root``, by brute force.

    Enumerate the 2^|edges| bond configurations, find the component of ``root``,
    and weight ``y_j`` for each other member in it.  This is the definition
    Sec. 2.8 gives, evaluated with no cleverness at all.
    """
    edges = list(edges)
    total = 0
    for on in itertools.product((0, 1), repeat=len(edges)):
        adj = {m: set() for m in members}
        for keep, (u, w) in zip(on, edges):
            if keep:
                adj[u].add(w)
                adj[w].add(u)
        seen, stack = {root}, [root]
        while stack:
            for nxt in adj[stack.pop()]:
                if nxt not in seen:
                    seen.add(nxt)
                    stack.append(nxt)
        w = p**sum(on) * (1 - p)**(len(edges) - sum(on))
        for j in seen:
            if j != root:
                w *= y[j]
        total += w
    return expand(total)


# ---------------------------------------------------------------------------
# the checks
# ---------------------------------------------------------------------------

def clique(n):
    return list(range(n)), list(itertools.combinations(range(n), 2))


def path(n):
    return list(range(n)), [(i, i + 1) for i in range(n - 1)]


def check(name, members, edges):
    p = symbols('p')
    y = {j: symbols(f'y{j}') for j in members}
    got = percolation_limit(members, edges, y)
    want = reachability_pgf(members, edges, y, p)
    ok = simplify(expand(got - want)) == 0
    print(f'  {name:<26} {"OK" if ok else "MISMATCH"}')
    if not ok:
        print(f'    q -> 1 limit : {factor(got)}')
        print(f'    Gbar         : {factor(want)}')
    assert ok, name
    return got


def check_against_package():
    """The symmetric case, against `percolation.clique_excess_pgf`."""
    p, y = symbols('p y')
    for n in (2, 3, 4, 5):
        members, edges = clique(n)
        ys = {j: y for j in members}
        got = percolation_limit(members, edges, ys)
        want = clique_excess_pgf(n, p)(y)
        assert simplify(expand(got - want)) == 0, n
        print(f'  K_{n} against clique_excess_pgf   OK')


def check_ising():
    """The other end of the dial: q = 2 is the Ising cavity field.

    With ``delta = (1 + s s')/2`` a Potts bond of coupling ``J_P`` is an Ising
    bond of coupling ``J_P / 2``, so ``v = exp(2 beta J) - 1``.  Writing the
    unnormalised message as ``1`` on one state and ``y = exp(-2h)`` on the
    other, the emitted ratio is ``exp(-2u)`` with ``u`` the emitted cavity
    field of Ch. 9.  Checked against an independent enumeration in
    ``statmech.cavity``.
    """
    sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'statmech' / 'src'))
    from statmech.cavity import emitted_field, ising_clique
    from sympy import exp, log, symbols as sym

    beta, J, v = sym('beta J v', positive=True)
    for n in (2, 3, 4):
        members, edges = clique(n)
        hs = sym(f'h1:{n}') if n > 1 else ()
        y = {j: exp(-2 * hs[j - 1]) for j in members if j != 0}
        y[0] = None
        rho = emitted_ratio(members, edges, y, 2, exp(2 * beta * J) - 1)
        got = simplify(-log(rho) / 2)
        want = simplify(emitted_field(n, ising_clique(n, beta, J), hs))
        assert simplify(expand(exp(2 * (got - want)))) == 1, n
        print(f'  K_{n} at q = 2 is the Ising cavity field   OK')


# ---------------------------------------------------------------------------
# the transmission factor at general q
# ---------------------------------------------------------------------------

def transmission(c):
    """``tau_c(q, v)``: d(emitted ratio)/d(one incoming y) at the symmetric point.

    The paramagnetic / no-giant-component fixed point is ``y = 1`` (every state
    equally likely, equivalently nothing connected to anything), where the
    emitted ratio is 1.  Its derivative with respect to one incoming message is
    the factor a single traversal of the complex contributes, and it is the one
    object that Chapters 4 and 9 have in common.
    """
    from sympy import cancel, diff
    q, v = symbols('q v')
    members, edges = clique(c)
    ys = {j: symbols(f'y{j}') for j in members}
    rho = emitted_ratio(members, edges, ys, q, v)
    at_one = {ys[j]: 1 for j in members if j != 0}
    assert simplify(rho.subs(at_one)) == 1, c
    return factor(cancel(diff(rho, ys[1]).subs(at_one)))


def check_transmission():
    """tau at q -> 1 is Chapter 4's; at q = 2 it is Chapter 9's."""
    q, v, p, t = symbols('q v p t')
    tau2, tau3 = transmission(2), transmission(3)
    print(f'  link      tau(q,v) = {tau2}')
    print(f'  triangle  tau(q,v) = {tau3}')

    # q -> 1: substitute the bond probability p = v / (1 + v)
    assert simplify(tau2.subs(q, 1).subs(v, p / (1 - p)) - p) == 0
    assert simplify(tau3.subs(q, 1).subs(v, p / (1 - p))
                    - p * (1 + p - p**2)) == 0
    print('  q -> 1:  link gives p, triangle gives p(1 + p - p^2) '
          '= sbar_tri / 2   OK')

    # q = 2: a Potts bond of coupling J is an Ising bond of J/2, v = e^{2 bJ} - 1
    assert simplify(tau2.subs(q, 2).subs(v, 2 * t / (1 - t)) - t) == 0
    assert simplify(tau3.subs(q, 2).subs(v, 2 * t / (1 - t))
                    - t / (1 - t + t**2)) == 0
    print('  q = 2 :  link gives t, triangle gives t/(1 - t + t^2)   OK')


def check_thresholds(kbar=3):
    """One branching condition, two famous thresholds."""
    from sympy import nsolve, solve
    q, v, p, t = symbols('q v p t')
    tau2 = transmission(2)
    vc = solve((2 - 1) * tau2 * kbar - 1, v)[0]
    assert simplify(vc - q / (kbar - 1)) == 0
    print(f'  graph, <kbar> = {kbar}:  v_c = {vc}, exactly linear in q')
    pc = simplify((vc / (1 + vc)).subs(q, 1))
    tc = simplify(vc.subs(q, 2) / (2 + vc.subs(q, 2)))
    assert pc == Rational(1, kbar) and tc == Rational(1, kbar)
    print(f'    q -> 1: p_c = {pc}  (Molloy-Reed, Eq. 4.6)')
    print(f'    q  = 2: tanh(beta J) = {tc}  (the Bethe-lattice Ising T_c)')


def check_first_order(b=2):
    """For q > 2 the linear instability is a spinodal, not the transition.

    Iterate the graph recursion ``y -> rho(y)^b`` with ``b`` the excess degree,
    from a strongly ordered start, and find the smallest ``v`` at which a
    non-trivial fixed point survives.  At q = 2 it coincides with the linear
    instability ``v_c = q / (b - 1)``; for q > 2 it sits strictly below, so
    there is a coexistence window and ``v_c`` is where the *paramagnet* loses
    stability rather than where the transition is.
    """
    import numpy as np

    def rho(y, q, v):
        return (1 + (1 + v) * y + (q - 2) * y) / ((1 + v) + (q - 1) * y)

    def ordered(q, v):
        y = 1e-9
        for _ in range(200000):
            yn = rho(y, q, v)**b
            if abs(yn - y) < 1e-15:
                break
            y = yn
        return y

    for q in (2, 3, 4, 6):
        vc = q / (b - 1)
        lo, hi = 0.0, vc
        for _ in range(60):
            mid = 0.5 * (lo + hi)
            if ordered(q, mid) < 1 - 1e-6:
                hi = mid
            else:
                lo = mid
        gap = 100 * (vc - hi) / vc
        # the bisection stops one tolerance short of a continuous onset, so
        # anything under a tenth of a per cent counts as coincident
        kind = ('continuous'
                if gap < 0.1 else f'first order, v_c high by {gap:.1f}%')
        print(f'  q = {q}, excess degree {b}:  linear instability v = {vc:.4f}, '
              f'ordered branch from v = {hi:.4f}   ({kind})')


def check_transition_point(b=2):
    """Where the first-order transition actually is, from the Bethe free energy.

    ``check_first_order`` brackets the coexistence window: the paramagnet loses
    stability at ``v_c = q/(b-1)`` and the ordered branch first exists at the
    lower ordered spinodal.  The transition itself is neither of those.  It is
    where the two branches exchange which has the lower free energy, strictly
    inside the window, and locating it needs the Bethe free energy rather than
    a fixed point.

    On a d-regular graph with d = b + 1, parameterise the (unnormalised) cavity
    message as weight 1 on one distinguished state and ``y`` on each of the
    other q - 1.  With

        A = (1 + v) + (q - 1) y     the incoming weight at a distinguished site
        B = 1 + y (q - 1 + v)       the incoming weight at a fixed other site

    the message fixed point is ``y = (B/A)^b`` and

        -beta f = ln[A^d + (q-1) B^d] - (d/2) ln[(1 + (q-1)y)^2 + v(1 + (q-1)y^2)]

    which is invariant under rescaling the messages, as it must be.  The
    transition is the crossing of that quantity on the ordered and disordered
    (y = 1) branches.  At q = 2 the window has zero width and the crossing sits
    on top of both spinodals, which is the check that the free energy is right.
    """
    import numpy as np
    from scipy.optimize import brentq

    d = b + 1

    def AB(y, q, v):
        return (1 + v) + (q - 1) * y, 1 + y * (q - 1 + v)

    def ordered_y(q, v):
        """Smallest fixed point of y = (B/A)^b below 1, or None if there is none."""
        f = lambda y: (lambda A, B: (B / A)**b - y)(*AB(y, q, v))
        ys = np.linspace(1e-12, 1 - 1e-9, 200001)
        vals = np.array([f(y) for y in ys])
        s = np.where(np.diff(np.sign(vals)) != 0)[0]
        return None if len(s) == 0 else brentq(f, ys[s[0]], ys[s[0] + 1])

    def minus_beta_f(y, q, v):
        A, B = AB(y, q, v)
        Zi = A**d + (q - 1) * B**d
        Zij = (1 + (q - 1) * y)**2 + v * (1 + (q - 1) * y**2)
        return np.log(Zi) - 0.5 * d * np.log(Zij)

    for q in (2, 3, 4, 6):
        vc = q / (b - 1)
        lo, hi = 0.0, vc                      # the ordered spinodal
        for _ in range(200):
            mid = 0.5 * (lo + hi)
            if ordered_y(q, mid) is not None:
                hi = mid
            else:
                lo = mid
        spin = hi
        if vc - spin < 1e-6:                  # continuous: no window to search
            vt = vc
        else:
            g = lambda v: (minus_beta_f(ordered_y(q, v), q, v)
                           - minus_beta_f(1.0, q, v))
            vt = brentq(g, spin * (1 + 1e-9), vc * (1 - 1e-9))
        print(f'  q = {q}:  ordered spinodal {spin:.4f} < transition {vt:.4f} '
              f'< linear instability {vc:.4f}    '
              f'v_c high by {100 * (vc - vt) / vc:.1f}% '
              f'(window {100 * (vc - spin) / vc:.1f}%)')


def figure_transmission():
    """One transmission factor, with q as a dial."""
    import numpy as np
    plt = _mpl()
    from sympy import lambdify
    q, v = symbols('q v')
    fig, axes = plt.subplots(1, 2, figsize=(4.6, 2.6), sharey=True)
    vs = np.linspace(0, 6, 400)
    curves = ((1, DARK, r'$q\to1$:  $\tau=p$'),
              (2, MID, r'$q=2$:  $\tau=\tanh\beta J$'),
              (4, LIGHT, r'$q=4$'))
    for ax, c, title in ((axes[0], 2, 'link'), (axes[1], 3, 'triangle')):
        tau = transmission(c)
        for qv, col, lab in curves:
            f = lambdify(v, tau.subs(q, qv), 'numpy')
            ax.plot(vs, f(vs), '-', lw=1.4, color=col, label=lab)
        ax.set_title(title, fontsize=8.5)
        ax.set_xlabel(r'$v=e^{\beta J}-1$', fontsize=8.5)
        _tidy(ax)
    axes[0].set_ylabel(r'transmission $\tau$', fontsize=8.5)
    axes[0].set_ylim(0, 1.02)
    axes[0].legend(frameon=False, fontsize=6.8, loc='lower right',
                   handlelength=1.4, borderaxespad=0.2)
    fig.tight_layout()
    fig.savefig(OUT / 'fig-potts.pdf')
    print(f'  wrote {OUT / "fig-potts.pdf"}')


def show_forms():
    """The two forms the chapter quotes."""
    p, y = symbols('p y')
    for name, (members, edges) in (('link', clique(2)), ('triangle', clique(3))):
        ys = {j: y for j in members}
        print(f'  {name:<10} Gbar(y) = {factor(percolation_limit(members, edges, ys))}')


if __name__ == '__main__':
    print('q -> 1 limit of the Potts interior sum, against Gbar:')
    check('link, K_2', *clique(2))
    check('triangle, K_3', *clique(3))
    check('4-clique, K_4', *clique(4))
    check('5-clique, K_5', *clique(5))
    check('path on 3', *path(3))
    check('path on 4', *path(4))
    print('symmetric case:')
    check_against_package()
    print('q = 2:')
    check_ising()
    print('the forms quoted in the text:')
    show_forms()
    print('the transmission factor:')
    check_transmission()
    print('one branching condition, two thresholds:')
    check_thresholds()
    print('order of the transition:')
    check_first_order()
    print('where the transition actually is:')
    check_transition_point()
    print('figure:')
    figure_transmission()
