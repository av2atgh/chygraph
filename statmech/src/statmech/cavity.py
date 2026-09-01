"""Cavity field updates emitted by a complex.

In the percolation map of :mod:`chygraph` a message is a scalar and the field
update is the identity: the occupation probability is folded into the moment
tables by thinning, so nothing is left to differentiate.  For a general
Hamiltonian the message is an effective field ``h`` and each traversal carries a
factor ``u'``, the derivative of the emitted field with respect to one incoming
field, evaluated at the trivial fixed point.

For a complex holding ``c`` members with internal energy ``-beta H(sigma)``, the
field emitted onto member 0 given cavity fields ``h_1..h_{c-1}`` on the others is

    u = (1/2) ln [ Z(sigma_0=+1) / Z(sigma_0=-1) ],
    Z(sigma_0) = sum_{sigma_1..sigma_{c-1}}
                 exp( -beta H(sigma) + sum_j h_j sigma_j ).

This module computes that sum *exactly inside the complex*.  That is the step
percolation does not need and a general model cannot avoid: the internal
structure of a complex is not summarised by a component-size generating
function once the interaction is not reachability.  It is also why complexes
are regions in the Kikuchi sense (WP5 of the README) — this module is that idea
at the smallest scale, exact enumeration for complexes small enough to enumerate.

Reference:
    Vazquez & Weigt, Phys. Rev. E 67, 027101 (2003), Eq. (5).
"""

import itertools

from sympy import Rational, Symbol, diff, exp, log, simplify, symbols, tanh


def emitted_field(c, minus_betaH, hs):
    """Field emitted onto member 0 of a ``c``-member complex.

    Args:
        c: number of members.
        minus_betaH: callable taking a length-``c`` tuple of ``+-1`` spins and
            returning ``-beta H``.
        hs: length ``c-1`` sequence of incoming cavity fields on members 1..c-1.
    """
    Z = {}
    for s0 in (1, -1):
        total = 0
        for rest in itertools.product((1, -1), repeat=c - 1):
            sigma = (s0,) + rest
            total += exp(minus_betaH(sigma) + sum(h * s for h, s in zip(hs, rest)))
        Z[s0] = total
    return Rational(1, 2) * log(Z[1] / Z[-1])


def cavity_derivative(c, minus_betaH, at=0):
    """``u'``: derivative of the emitted field w.r.t. one incoming field.

    Evaluated at ``h = at`` on every incoming edge, which for a spin-flip
    symmetric complex is the paramagnetic fixed point ``h = 0``.  Returns a
    single number per (complex type, member role) because at the trivial fixed
    point every incoming field sits at the same value; away from it the average
    ``<u'>`` must be taken over the field distribution and this shortcut fails.
    """
    hs = symbols(f'h1:{c}') if c > 1 else ()
    u = emitted_field(c, minus_betaH, hs)
    d = diff(u, hs[0])
    return simplify(d.subs({h: at for h in hs}))


# ---------------------------------------------------------------------------
# Named complexes
# ---------------------------------------------------------------------------

def ising_clique(n, beta=None, J=None):
    """``-beta H`` for ``n`` spins with ferromagnetic coupling on every pair."""
    beta = Symbol('beta', positive=True) if beta is None else beta
    J = Symbol('J', positive=True) if J is None else J

    def energy(sigma):
        return beta * J * sum(sigma[i] * sigma[j]
                              for i in range(n) for j in range(i + 1, n))
    return energy


def ising_edge_derivative(beta=None, J=None):
    """``u' = tanh(beta J)`` — the textbook Ising cavity factor."""
    return cavity_derivative(2, ising_clique(2, beta, J))


def ising_triangle_derivative(beta=None, J=None):
    """``u'`` emitted by a 3-clique.

    Equals ``t / (1 - t + t^2)`` with ``t = tanh(beta J)``, which exceeds ``t``
    for ``0 < t < 1``: a triangle transmits more per neighbour than an
    independent edge does.  This is the clustering effect that the
    ``{p_d, e_dd'}`` ensemble of Vazquez & Weigt cannot represent.
    """
    return cavity_derivative(3, ising_clique(3, beta, J))


def in_tanh(expr, beta=None, J=None, t=None):
    """Rewrite an expression in ``t = tanh(beta J)``.

    Substitutes ``exp(2 beta J) -> (1 + t) / (1 - t)`` by way of the algebraic
    variable ``X = exp(2 beta J)``, which is what the enumeration in
    :func:`emitted_field` naturally produces.
    """
    beta = Symbol('beta', positive=True) if beta is None else beta
    J = Symbol('J', positive=True) if J is None else J
    t = Symbol('t', positive=True) if t is None else t
    X = Symbol('_X', positive=True)
    e = expr.rewrite(exp).subs(exp(beta * J), X**Rational(1, 2))
    e = simplify(e.subs(X, (1 + t) / (1 - t)))
    return simplify(e)


def tanh_of(beta=None, J=None):
    beta = Symbol('beta', positive=True) if beta is None else beta
    J = Symbol('J', positive=True) if J is None else J
    return tanh(beta * J)
