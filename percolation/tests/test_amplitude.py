"""Tests for the symbolic critical amplitude.

Every closed form here is cross-checked against the independent numeric path in
``GiantComponent.amplitude`` (numpy eigenvectors, full index space), which is in
turn anchored to Monte Carlo through ``test_giant.py``.
"""

import math

from sympy import symbols, nsolve, simplify, expand, Rational

from percolation.giant import (
    hypergraph_giant, multiplex_hypergraph_giant, graph_with_triangles_giant,
)
from percolation.amplitude import CriticalAmplitude


# ---------------------------------------------------------------------------
# Core reduction
# ---------------------------------------------------------------------------

def test_core_is_the_physical_unknowns():
    """8 -> 2 for a hypergraph, 18 -> 4 for a graph with triangles."""
    assert CriticalAmplitude(hypergraph_giant()).labels() == [('-', 0, 1), ('+', 1, 0)]
    assert CriticalAmplitude(graph_with_triangles_giant()).labels() == [
        ('-', 0, 1), ('-', 0, 2), ('+', 1, 0), ('+', 2, 0)]


def test_reduction_is_exact_in_the_full_space():
    for G in (hypergraph_giant(), hypergraph_giant(graph=True),
              multiplex_hypergraph_giant(2), graph_with_triangles_giant()):
        assert all(CriticalAmplitude(G).verify().values())


# ---------------------------------------------------------------------------
# Perron root reproduces the published threshold
# ---------------------------------------------------------------------------

def test_perron_root_is_sqrt_theta_plus_one():
    """Lambda = sqrt(theta + 1) - 1 for these models."""
    for G in (hypergraph_giant(), multiplex_hypergraph_giant(2),
              graph_with_triangles_giant()):
        C = CriticalAmplitude(G)
        assert simplify(C.perron_root()**2 - (G.theta() + 1)) == 0


# ---------------------------------------------------------------------------
# Closed forms on the critical manifold
# ---------------------------------------------------------------------------

def test_graph_amplitude_is_four_p():
    p, q = symbols('p q')
    B = CriticalAmplitude(hypergraph_giant(graph=True)).amplitude_at_threshold(0, q)
    assert simplify(B - 4 * p) == 0


def test_hypergraph_amplitude():
    p, q, c = symbols('p q c')
    B = CriticalAmplitude(hypergraph_giant()).amplitude_at_threshold(0, q)
    assert simplify(B - 4 * p / (1 + c * p)) == 0


def test_multiplex_amplitude_is_one_over_theta_weights():
    """B = 4 / (1 + sum_l <c>_l0^2 <k>_0l q_l) on the critical manifold."""
    k01, k02, c10, c20, q1, q2 = symbols('k_01 k_02 c_10 c_20 q_1 q_2')
    B = CriticalAmplitude(multiplex_hypergraph_giant(2)).amplitude_at_threshold(0, q1)
    want = 4 / (1 + c10**2 * k01 * q1 + c20**2 * k02 * q2)
    # equal only on the manifold, so eliminate q1 from the target as well
    sub = {q1: (1 - c20 * k02 * q2) / (c10 * k01)}
    assert simplify(B.subs(sub) - want.subs(sub)) == 0


def test_triangle_amplitude_does_not_depend_on_the_link_layer():
    """B = 4 / (1 + 2 k_T q^2 (3 - 2q)): links contribute no curvature because a
    link's excess component size is Bernoulli, whose PGF is affine."""
    q, kL, kT = symbols('q k_L k_T')
    B = CriticalAmplitude(graph_with_triangles_giant()).amplitude_at_threshold(0, kL)
    assert simplify(B - 4 / (1 + 2 * kT * q**2 * (3 - 2 * q))) == 0
    assert kL not in B.free_symbols


# ---------------------------------------------------------------------------
# Agreement with the independent numeric path, and with the map itself
# ---------------------------------------------------------------------------

def test_symbolic_matches_numeric_amplitude():
    cases = [
        (hypergraph_giant(), {'k': 2, 'c': 2, 'p': 1, 'q': 0.25}),
        (hypergraph_giant(graph=True), {'k': 1, 'p': 1, 'q': 1}),
        (hypergraph_giant(graph=True), {'k': 4, 'p': 0.5, 'q': 0.5}),
    ]
    for G, subs in cases:
        sym = float(CriticalAmplitude(G).amplitude(0).subs(
            {symbols(k): v for k, v in subs.items()}))
        assert abs(sym - G.amplitude_numeric(subs)) < 1e-9, (subs, sym)


def test_amplitude_predicts_the_order_parameter_near_threshold():
    """S / Lambda -> B as the threshold is approached from above."""
    q, kL, kT = symbols('q k_L k_T')
    T = graph_with_triangles_giant()
    C = CriticalAmplitude(T)
    base = {kL: 1, kT: 0.5}
    qc = float(nsolve(C.Lambda().subs(base), q, 0.4))
    B = float(C.amplitude(0).subs({**base, q: qc}))
    ratios = []
    for dq in (0.002, 0.0005):
        S = T.node_fraction({'k_L': 1, 'k_T': 0.5, 'q': qc + dq},
                            tol=1e-16, maxiter=4_000_000)
        lam = float(C.perron_root().subs({**base, q: qc + dq}))
        ratios.append(S / (lam - 1))
    assert abs(ratios[-1] - B) < 5e-3
    assert abs(ratios[-1] - B) < abs(ratios[0] - B)   # converging


def test_curvature_is_positive_for_a_continuous_transition():
    C = CriticalAmplitude(hypergraph_giant())
    assert C.is_continuous()
    val = float(C.curvature().subs({symbols(k): v for k, v in
                                    {'k': 2, 'c': 2, 'p': 1, 'q': 0.25}.items()}))
    assert val > 0


# ---------------------------------------------------------------------------
# Arbitrary distributions: only the first two factorial moments enter B
# ---------------------------------------------------------------------------

def _hypergraph_moments(graph=False):
    """Site-bond percolation with the degree and cardinality distributions
    represented by their first two factorial moments only."""
    from percolation.giant import GiantComponent, thin, moment_pgf as mp
    p, q = symbols('p q')
    k, kb, kb1 = symbols('k kbar kbar1')
    c, cb, cb1 = symbols('c cbar cbar1')
    card, cardbar = ((mp(2, 2), mp(1, 0)) if graph else (mp(c, c * cb), mp(cb, cb1)))
    return GiantComponent(
        [[None, thin(mp(k, k * kb), p)], [None, None]],
        [[None, thin(mp(kb, kb1), p)], [None, None]],
        [[None, None], [thin(card, q), None]],
        [[None, None], [thin(cardbar, q), None]],
    )


_PROBE = {symbols(s): v for s, v in
          {'p': .4, 'q': .5, 'k': 1.3, 'kbar': 1.7, 'kbar1': 2.1,
           'c': 1.1, 'cbar': 1.9, 'cbar1': 2.3}.items()}


def test_general_graph_amplitude():
    """B = 4 p <k> <kbar> / <kbar(kbar-1)> on p q <kbar> = 1."""
    p, q, k, kb, kb1 = symbols('p q k kbar kbar1')
    B = CriticalAmplitude(_hypergraph_moments(graph=True),
                          probe=_PROBE).amplitude_at_threshold(0, q)
    assert simplify(B - 4 * p * k * kb / kb1) == 0


def test_general_hypergraph_amplitude():
    """B = 4 p <k><kbar><cbar> / (<cbar><kbar(kbar-1)> + p <kbar>^2 <cbar(cbar-1)>)."""
    p, q, k, kb, kb1, cb, cb1 = symbols('p q k kbar kbar1 cbar cbar1')
    B = CriticalAmplitude(_hypergraph_moments(),
                          probe=_PROBE).amplitude_at_threshold(0, q)
    want = 4 * p * k * kb * cb / (cb * kb1 + p * kb**2 * cb1)
    assert simplify(B - want) == 0


def test_general_formula_reproduces_the_poisson_case():
    """Poisson: <kbar> = <k>, <kbar(kbar-1)> = <k>^2, and likewise for c."""
    p, q, k, kb, kb1, c, cb, cb1 = symbols('p q k kbar kbar1 c cbar cbar1')
    B = CriticalAmplitude(_hypergraph_moments(),
                          probe=_PROBE).amplitude_at_threshold(0, q)
    poisson = {kb: k, kb1: k**2, cb: c, cb1: c**2}
    assert simplify(B.subs(poisson) - 4 * p / (1 + c * p)) == 0


def test_general_formula_reproduces_the_two_point_distributions_of_fig2():
    """The two degree distributions with identical <k>, <kbar> of Fig. 2."""
    p, q, k, kb, kb1 = symbols('p q k kbar kbar1')
    B = CriticalAmplitude(_hypergraph_moments(graph=True),
                          probe=_PROBE).amplitude_at_threshold(0, q)
    # A: P(k) = {0: 1/2, 4: 1/2}   -> kbar = 3 always, <kbar(kbar-1)> = 6
    # B: P(k) = {0: .2, 1: .5, 5: .3} -> kbar in {0, 4}, <kbar> = 3, second = 9
    at = lambda m: float(B.subs({k: 2, kb: 3, kb1: m, p: 1 / 3}))
    assert abs(at(6) - 4 / 3) < 1e-12
    assert abs(at(9) - 8 / 9) < 1e-12


def test_poisson_multiplex_amplitude_is_additive_like_theta():
    """B = 4 / (1 + sum_l <c>_l0^2 <k>_0l q_l), mirroring theta = sum_l q_l <k>_l <c>_l - 1."""
    for L in (1, 2):
        C = CriticalAmplitude(multiplex_hypergraph_giant(number_of_types=L))
        q1 = symbols('q_1')
        B = C.amplitude_at_threshold(0, q1)
        want = 4 / (1 + sum(symbols(f'c_{l}0')**2 * symbols(f'k_0{l}') * symbols(f'q_{l}')
                            for l in range(1, L + 1)))
        sub = {q1: (1 - sum(symbols(f'c_{l}0') * symbols(f'k_0{l}') * symbols(f'q_{l}')
                            for l in range(2, L + 1)))
                   / (symbols('c_10') * symbols('k_01'))}
        assert simplify(B.subs(sub) - want.subs(sub)) == 0


# ---------------------------------------------------------------------------
# Can B diverge?  Only in the degenerate one-dimensional limit
# ---------------------------------------------------------------------------

def test_curvature_vanishes_only_for_the_two_regular_graph():
    """C = l . M[r, r] is a sum of non-negative terms.  It vanishes only when
    every step along the critical direction is deterministic and unique, which
    for a graph means k = 2: a union of cycles, where lambda = sqrt(pq) is
    critical only at the boundary p = q = 1 and every Q is a fixed point."""
    from percolation.giant import constant_pgf
    p, q = symbols('p q')
    for kk, want_zero in ((2, True), (3, False), (4, False)):
        G = hypergraph_giant(degree=constant_pgf(kk),
                             excess_degree=constant_pgf(kk - 1),
                             graph=True, poisson=False)
        C = CriticalAmplitude(G)
        assert (simplify(C.curvature()) == 0) is want_zero, kk
        assert C.is_continuous() is not want_zero
    G2 = hypergraph_giant(degree=constant_pgf(2), excess_degree=constant_pgf(1),
                          graph=True, poisson=False)
    # lambda = sqrt(pq): the transition sits at the edge of the parameter range
    assert simplify(CriticalAmplitude(G2).perron_root() - (p * q)**Rational(1, 2)) == 0


def test_amplitude_diverges_on_approach_to_the_two_regular_limit():
    """P(k) = (1-eps) delta_{k2} + eps delta_{k4} tends to a union of cycles as
    eps -> 0.  Then B = 4<k>^2 / [m(m-1)(m-2) eps] diverges like 1/eps: an
    amplitude can be made arbitrarily large, but by degenerating towards one
    dimension, not by any explosive mechanism."""
    from percolation.giant import finite_pgf
    p, q = symbols('p q')
    m = 4
    for eps in (Rational(1, 4), Rational(1, 10), Rational(1, 100)):
        P = {2: 1 - eps, m: eps}
        kmean = sum(kk * v for kk, v in P.items())
        G = hypergraph_giant(
            degree=finite_pgf(P),
            excess_degree=finite_pgf({kk - 1: kk * v / kmean for kk, v in P.items()}),
            graph=True, poisson=False)
        C = CriticalAmplitude(G)
        pc = float(nsolve(C.Lambda().subs(q, 1), p, 0.5))
        B = float(C.amplitude(0).subs({q: 1, p: pc}))
        want = 4 * float(kmean)**2 / (m * (m - 1) * (m - 2) * float(eps))
        assert abs(B - want) < 1e-6 * want, (eps, B, want)
