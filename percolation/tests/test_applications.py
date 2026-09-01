"""Tests for literature constructions mapped to chygraphs.

Reference values marked "MC" were confirmed against Monte Carlo simulation with
n = 2-4 x 10^5 and agree within one to two standard errors; the assertions here
pin the theory so a regression in the mapping fails loudly.
"""

import math

from sympy import Rational, expand, nsolve, simplify, symbols, diff

from percolation.applications import (
    clique_cluster_distribution, clique_excess_pgf, size_biased,
    and_or_hypergraph, correlated_cardinality_hypergraph,
    two_class_joint_degree, household_epidemic, clique_network, _mixture,
    stc_percolation,
)
from percolation.amplitude import CriticalAmplitude

q, p, pH, T, k, y = symbols('q p p_H T k y')


# ---------------------------------------------------------------------------
# Bond percolation inside K_n
# ---------------------------------------------------------------------------

def test_clique_cluster_distribution_normalises():
    for n in (2, 3, 4, 5):
        assert simplify(sum(clique_cluster_distribution(n, q).values()) - 1) == 0


def test_K3_reproduces_the_published_triangle_enumeration():
    """cnae047 Fig. 3: excess size 2, 1, 0 with the stated probabilities."""
    d = clique_cluster_distribution(3, q)
    want = {3: 3 * q**2 - 2 * q**3,
            2: 2 * q * (1 - q)**2,
            1: (1 - q)**3 + q * (1 - q)**2}
    for j in (1, 2, 3):
        assert simplify(d[j] - want[j]) == 0
    mean_excess = sum((j - 1) * pr for j, pr in d.items())
    assert simplify(mean_excess - 2 * q * (1 + q - q**2)) == 0


# ---------------------------------------------------------------------------
# 1. AND-logic versus OR-logic hypergraph percolation
# ---------------------------------------------------------------------------

def test_or_logic_recovers_the_factor_graph_threshold():
    """OR: theta = p <kbar><cbar> - 1, the standard factor graph result."""
    c = symbols('c')
    assert simplify(and_or_hypergraph('or').theta() - (c * k * p - 1)) == 0


def test_and_logic_threshold():
    """AND: <sbar> = p Gbar_c'(p), giving theta = p<kbar> c e^{c(p-1)} - 1
    for Poisson cardinalities."""
    c = symbols('c')
    from sympy import exp
    got = and_or_hypergraph('and').theta()
    assert simplify(got - (c * k * p * exp(c * (p - 1)) - 1)) == 0


def test_and_threshold_exceeds_or_threshold():
    """Bianconi & Dorogovtsev: node percolation on hypergraphs is harder than
    on the corresponding factor graph."""
    OR, AND = and_or_hypergraph('or'), and_or_hypergraph('and')
    sub = {k: 2, symbols('c'): 3}
    p_or = float(nsolve(OR.theta().subs(sub), p, 0.2))
    p_and = float(nsolve(AND.theta().subs(sub), p, 0.6))
    assert p_and > p_or
    assert abs(p_or - 1 / 6) < 1e-9


def test_and_and_or_coincide_at_full_occupation():
    OR, AND = and_or_hypergraph('or'), and_or_hypergraph('and')
    sub = {'k': 2, 'c': 3, 'p': 1}
    assert abs(OR.node_fraction(sub) - AND.node_fraction(sub)) < 1e-12


def test_and_or_order_parameters_match_simulation():
    ref = {('or', 0.4): 0.26723, ('or', 0.6): 0.46712, ('or', 0.8): 0.65693,
           ('or', 1.0): 0.84112, ('and', 0.8): 0.44308, ('and', 1.0): 0.84112}
    for (logic, pv), want in ref.items():
        got = and_or_hypergraph(logic).node_fraction({'k': 2, 'c': 3, 'p': pv})
        assert abs(got - want) < 1e-4, (logic, pv, got)


# ---------------------------------------------------------------------------
# 2. Hyperdegree-cardinality correlations
# ---------------------------------------------------------------------------

CARD = [2, 5]
MA, MB, SP = Rational(1), Rational(3, 5), Rational(4, 5)


def _corr(w):
    return correlated_cardinality_hypergraph(
        CARD, two_class_joint_degree(MA, MB, SP, Rational(w).limit_denominator(100)))


def test_correlation_family_has_fixed_marginals():
    """Every member has the same <kappa>_0l and the same within-layer excess,
    so the published threshold tensor cannot tell them apart."""
    base = _corr(Rational(1, 2))
    for w in (0, Rational(1, 4), 1):
        M = _corr(w)
        assert M.kappa(0, 1) == base.kappa(0, 1) == 1
        assert M.kappa(0, 2) == base.kappa(0, 2) == Rational(3, 5)
        assert simplify(M.kappa_bar(1, 0, 1) - base.kappa_bar(1, 0, 1)) == 0
        assert M.layers_independent() == (w == Rational(1, 2))


def test_correlation_shifts_the_threshold_monotonically():
    pc = [float(nsolve(_corr(w).theta(), p, 0.3)) for w in (0, 0.5, 1.0)]
    assert abs(pc[0] - 0.24604) < 1e-4
    assert abs(pc[1] - 0.21204) < 1e-4
    assert abs(pc[2] - 0.17934) < 1e-4
    assert pc[0] > pc[1] > pc[2]     # positive correlation lowers the threshold


def test_correlation_order_parameter_matches_simulation():
    """The ordering of S in the correlation reverses between p = 0.5 and p = 1:
    the curves cross, so no single statement 'correlation helps or hurts' holds."""
    ref = {(0, 0.5): 0.22043, (0, 1.0): 0.74272,
           (0.5, 0.5): 0.24226, (0.5, 1.0): 0.66980,
           (1, 0.5): 0.23690, (1, 1.0): 0.59052}
    for (w, pv), want in ref.items():
        got = _corr(w).node_fraction({'p': pv})
        assert abs(got - want) < 1e-4, (w, pv, got)
    # non-monotonic at p = 0.5, monotonic at p = 1
    assert ref[(0.5, 0.5)] > ref[(0, 0.5)] and ref[(0.5, 0.5)] > ref[(1, 0.5)]
    assert ref[(0, 1.0)] > ref[(0.5, 1.0)] > ref[(1, 1.0)]


# ---------------------------------------------------------------------------
# 3. Epidemics with two levels of mixing
# ---------------------------------------------------------------------------

def _mu_H(sizes):
    gb = _mixture(size_biased(sizes), lambda n: clique_excess_pgf(n, pH))
    return simplify(diff(gb(y), y).subs(y, 1))


def test_household_threshold_is_the_ball_mollison_reproduction_number():
    """theta + 1 = R* = T [<kbar> + mu_H <k>], with Poisson global degree."""
    for sizes in ({3: 1}, {2: Rational(1, 2), 4: Rational(1, 2)},
                  {1: Rational(1, 4), 5: Rational(3, 4)}):
        M = household_epidemic(sizes)
        Rstar = T * (k + _mu_H(sizes) * k)
        assert simplify(M.theta() - (Rstar - 1)) == 0, sizes


def test_household_with_no_within_household_transmission_is_a_plain_graph():
    """p_H = 0 leaves only the global network: theta = T<kbar> - 1."""
    M = household_epidemic({3: 1})
    assert simplify(M.theta().subs(pH, 0) - (T * k - 1)) == 0


def test_household_final_size_matches_simulation():
    ref = {(0.5, 0.3, 2.0): 0.41961, (0.8, 0.2, 2.0): 0.23292,
           (0.3, 0.5, 1.5): 0.35394, (0.0, 0.6, 2.0): 0.31370}
    H = household_epidemic({3: 1})
    for (a, b, c), want in ref.items():
        got = H.node_fraction({'p_H': a, 'T': b, 'k': c})
        assert abs(got - want) < 1e-4, (a, b, c, got)


def test_household_subcritical_below_R_star_one():
    H = household_epidemic({3: 1})
    sub = {'p_H': 0.9, 'T': 0.1, 'k': 3.0}
    assert float(H.theta().subs({pH: 0.9, T: 0.1, k: 3.0})) < 0
    assert H.node_fraction(sub) < 1e-6


def test_households_alone_never_percolate():
    """With T = 0 an individual belongs to exactly one household, so there is no
    route from one household to another and theta = -1 for any p_H."""
    assert simplify(household_epidemic({3: 1}, T=0).theta() + 1) == 0
    assert simplify(household_epidemic({2: Rational(1, 2), 6: Rational(1, 2)},
                                       T=0).theta() + 1) == 0


def test_clique_network_of_pairs_is_ordinary_bond_percolation():
    """Cliques of size two: theta = q<kbar> - 1."""
    assert simplify(clique_network({2: 1}, p_bond=q).theta() - (q * k - 1)) == 0


def test_clique_network_of_triangles_matches_the_triangle_model():
    """Cliques of size three must agree with graph_with_triangles_giant at
    <k>_| = 0, which is an independent code path."""
    from percolation.giant import graph_with_triangles_giant
    got = clique_network({3: 1}, p_bond=q).theta()
    want = graph_with_triangles_giant().theta().subs(
        {symbols('k_L'): 0, symbols('k_T'): k})
    assert expand(simplify(got - want)) == 0


# ---------------------------------------------------------------------------
# The amplitude machinery composes with all of them
# ---------------------------------------------------------------------------

def test_amplitude_composes_with_the_application_models():
    for M, sub, guess, var in (
        (and_or_hypergraph('and'), {k: 2, symbols('c'): 3}, 0.6, p),
        (household_epidemic({3: 1}), {T: 0.3, k: 2.0}, 0.2, pH),
        (_corr(1), {}, 0.2, p),
    ):
        C = CriticalAmplitude(M)
        assert all(C.verify().values())
        xc = float(nsolve(C.Lambda().subs(sub), var, guess))
        B = float(C.amplitude(0).subs({**sub, var: xc}))
        assert B > 0
        S = M.node_fraction({**{str(s): v for s, v in sub.items()},
                             str(var): xc + 1e-3},
                            tol=1e-16, maxiter=4_000_000)
        lam = float(C.perron_root().subs({**sub, var: xc + 1e-3}))
        assert abs(S / (lam - 1) - B) < 0.15 * B


# ---------------------------------------------------------------------------
# 4. Static Triadic Closure / extended-range R = 2 percolation
# ---------------------------------------------------------------------------

phi = symbols('phi')


def _cirigliano_S(p, b):
    """Independent implementation of Eqs. (9)-(11) of arXiv:2506.17175 for a
    Poisson backbone, iterated directly rather than through the library."""
    g = lambda z: math.exp(b * (z - 1))
    u1 = u2 = 0.0
    for _ in range(400000):
        n1, n2 = g(p * u1 + (1 - p) * u2), g(p * u1 + 1 - p)
        if abs(n1 - u1) + abs(n2 - u2) < 1e-15:
            break
        u1, u2 = n1, n2
    return p * (1 - g(p * u1 + (1 - p) * u2))


def test_stc_threshold_is_cirigliano_eq12():
    """theta = -b^2 phi^2 + b(1+b) phi - 1, so
    phi_c = [(1+b) - sqrt((1+b)^2 - 4)] / (2b)."""
    from sympy import sqrt, solve
    b = symbols('k')                      # Poisson backbone: b = <k(k-1)>/<k> = <k>
    M = stc_percolation()
    assert expand(simplify(M.theta() - (-b**2 * phi**2 + b * (1 + b) * phi - 1))) == 0
    want = (1 + b - sqrt((1 + b)**2 - 4)) / (2 * b)
    assert any(simplify(s - want) == 0 for s in solve(M.theta(), phi))


def test_stc_order_parameter_matches_the_exact_solution():
    """Machine-precision agreement with Eqs. (9)-(11), themselves validated
    against simulation of STC graphs (n = 1.2e5)."""
    M = stc_percolation()
    for p in (0.10, 0.15, 0.25, 0.40, 0.60, 0.90):
        got = M.node_fraction({'phi': p, 'k': 3.0})
        assert abs(got - _cirigliano_S(p, 3.0)) < 1e-11, (p, got)


def test_stc_matches_simulation():
    """Reference values confirmed against Monte Carlo on STC graphs built from
    a Poisson(3) configuration-model backbone, n = 1.2e5."""
    M = stc_percolation()
    for p, want in ((0.10, 0.016155), (0.15, 0.083846), (0.25, 0.200429),
                    (0.40, 0.358276), (0.60, 0.557487)):
        assert abs(M.node_fraction({'phi': p, 'k': 3.0}) - want) < 1e-5


def test_the_naive_stc_mapping_is_wrong():
    """Taking G1 itself as a hypergraph of closed neighbourhoods gives an
    incidence graph with four-cycles, and the prediction fails.  Pinned so that
    nobody 'simplifies' stc_percolation onto G1."""
    from sympy import exp, nsolve
    k = symbols('k')
    Phi = lambda z: z * exp(3 * (z - 1))                       # kappa = 1 + Poisson(3)
    Phibar = lambda z: exp(3 * (z - 1)) * (1 + 3 * z) / 4      # its excess
    naive = and_or_hypergraph('or', degree=Phi, excess_degree=Phibar,
                              cardinality=Phi, excess_cardinality=Phibar,
                              poisson=False)
    p = symbols('p')
    pc_naive = float(nsolve(naive.theta(), p, 0.1))
    pc_exact = float(nsolve(stc_percolation().theta().subs(symbols('k'), 3), phi, 0.1))
    assert abs(pc_exact - 0.089316) < 1e-5
    assert abs(pc_naive - 0.071111) < 1e-5
    assert pc_naive < pc_exact            # and it is wrong by 20%
