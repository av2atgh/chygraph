"""The three calculations are three orders of one object.

A chygraph is specified once, by its generating functions, and then answers
where the transition is (theta, Lambda), how steeply the order parameter rises
out of it (B), and what the order parameter is away from it (S).  These tests
pin that surface and check the delegated route against the underlying one.
"""

import math

from sympy import symbols, nsolve, simplify, sqrt

from percolation import (Chygraph, GiantComponent, JointChygraph,
                      JointGiantComponent, CriticalAmplitude,
                      hypergraph_giant, graph_with_triangles_giant,
                      stc_percolation)

p, q, k, c = symbols('p q k c')


def test_one_object_answers_all_three_questions():
    M = hypergraph_giant()
    assert simplify(M.theta() - (c * k * p * q - 1)) == 0          # threshold
    assert simplify(M.Lambda() - (sqrt(c * k * p * q) - 1)) == 0   # order param of it
    assert simplify(M.amplitude() - 4 * p / (c * p + sqrt(c * k * p * q))) == 0
    assert abs(M.node_fraction({'k': 3, 'c': 3, 'p': 0.6, 'q': 0.8}) - 0.5082196) < 1e-6


def test_delegation_agrees_with_the_underlying_class():
    for M in (hypergraph_giant(), graph_with_triangles_giant(), stc_percolation()):
        C = CriticalAmplitude(M)
        assert M.core() == C.core
        assert simplify(M.perron_root() - C.perron_root()) == 0
        assert simplify(M.Lambda() - C.Lambda()) == 0
        assert simplify(M.amplitude() - C.amplitude()) == 0
        assert simplify(M.curvature() - C.curvature()) == 0
        assert M.is_continuous() == C.is_continuous()
        assert M.verify() == C.verify()


def test_symbolic_and_numeric_amplitudes_agree():
    """amplitude() works on the reduced core symbolically; amplitude_numeric()
    uses numpy eigenvectors on the full index set.  They are independent."""
    for M, subs in ((hypergraph_giant(), {'k': 2, 'c': 2, 'p': 1, 'q': 0.25}),
                    (hypergraph_giant(graph=True), {'k': 4, 'p': 0.5, 'q': 0.5})):
        sym = float(M.amplitude().subs({symbols(a): v for a, v in subs.items()}))
        assert abs(sym - M.amplitude_numeric(subs)) < 1e-9


def test_amplitude_is_the_slope_of_the_order_parameter():
    """The three answers are consistent: S/Lambda -> B as Lambda -> 0."""
    T = graph_with_triangles_giant()
    kL, kT = symbols('k_L k_T')
    base = {kL: 1, kT: 0.5}
    qc = float(nsolve(T.Lambda().subs(base), q, 0.4))
    B = float(T.amplitude().subs({**base, q: qc}))
    S = T.node_fraction({'k_L': 1, 'k_T': 0.5, 'q': qc + 1e-3},
                        tol=1e-16, maxiter=4_000_000)
    lam = float(T.perron_root().subs({**base, q: qc + 1e-3})) - 1
    assert abs(S / lam - B) < 5e-3 * B


def test_joint_input_shares_the_whole_surface():
    """JointChygraph differs only in how the chygraph is specified."""
    M = stc_percolation()
    assert isinstance(M, JointChygraph) and isinstance(M, Chygraph)
    for name in ('theta', 'Lambda', 'amplitude', 'curvature', 'node_fraction',
                 'perron_root', 'verify', 'core'):
        assert hasattr(M, name), name


def test_old_names_still_resolve():
    assert GiantComponent is Chygraph
    assert JointGiantComponent is JointChygraph
