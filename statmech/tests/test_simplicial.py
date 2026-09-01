"""The simplicial Ising model of Son, Lee & Goh, as a chygraph.

arXiv:2411.19080; Commun. Phys. (2026), doi:10.1038/s42005-026-02724-2.
"""

import numpy as np
import pytest

from statmech.simplicial import (SimplicialChygraph, emitted,
                                          uprime)


# ---------------------------------------------------------------------------
# The Hamiltonian drops into the existing intra-complex step
# ---------------------------------------------------------------------------

@pytest.mark.parametrize('bJ', [0.3, 0.91, 2.0])
@pytest.mark.parametrize('h', [-1.5, -0.2, 0.0, 0.3, 2.0])
def test_pairwise_simplicial_is_a_half_coupling_bond(bJ, h):
    """delta_{S0 S1} = (1 + S0 S1)/2, so a simplicial pair is an ordinary Ising
    bond of half the coupling and the emitted field must be the textbook one."""
    assert float(emitted(2, bJ, h)) == pytest.approx(
        np.arctanh(np.tanh(bJ / 2) * np.tanh(h)), abs=1e-12)


def _emit_independent(q, a, h):
    """Emitted field with q-1 *independent* member fields, by enumeration."""
    import itertools
    Z = {}
    for s0 in (1, -1):
        tot = 0.0
        for rest in itertools.product((1, -1), repeat=q - 1):
            e = a if all(s == s0 for s in rest) else 0.0
            tot += np.exp(e + sum(hi * si for hi, si in zip(h, rest)))
        Z[s0] = tot
    return 0.5 * np.log(Z[1] / Z[-1])


@pytest.mark.parametrize('q', [2, 3, 5])
@pytest.mark.parametrize('bJ', [0.2, 0.8, 3.0])
def test_uprime_is_per_neighbour(q, bJ):
    """u' = (e^a - 1)/(2^{q-1} + e^a - 1), the derivative with respect to ONE
    other member's field, matching ising.clique_derivative's convention."""
    eps = 1e-6
    hp = [eps] + [0.0] * (q - 2)
    hm = [-eps] + [0.0] * (q - 2)
    num = (_emit_independent(q, bJ, hp) - _emit_independent(q, bJ, hm)) / (2 * eps)
    assert uprime(q, bJ) == pytest.approx(num, abs=1e-7)


@pytest.mark.parametrize('q', [2, 3, 5, 16])
@pytest.mark.parametrize('bJ', [0.2, 0.8, 3.0])
def test_common_field_derivative_is_q_minus_one_times_uprime(q, bJ):
    """The multiplicity trap: differentiating with respect to a field common to
    all q-1 members gives (q-1) u'.  Using that in the branching matrix, which
    supplies c-1 itself, counts the multiplicity twice."""
    eps = 1e-6
    common = float((emitted(q, bJ, eps) - emitted(q, bJ, -eps)) / (2 * eps))
    assert common == pytest.approx((q - 1) * uprime(q, bJ), rel=1e-5)


@pytest.mark.parametrize('bJ', [1e-6, 0.8, 50.0, 900.0])
def test_uprime_survives_low_temperature(bJ):
    """e^{beta J} - 1 overflows past beta J ~ 700, which any spinodal search
    reaches; the log-space form does not."""
    u = uprime(5, bJ)
    assert np.isfinite(u) and 0.0 <= u <= 4.0


def test_uprime_saturates_at_one():
    """Per neighbour, u' -> 1 at zero temperature for every cardinality."""
    for q in (2, 5, 16):
        assert uprime(q, 900.0) == pytest.approx(1.0, rel=1e-12)


def test_both_routes_to_the_spinodal_agree():
    """SimplicialChygraph supplies (q-1) itself; Chygraph.branching_matrix does
    too.  Their Perron roots must both be 1 at the spinodal."""
    from statmech import Chygraph
    for q, k in ((2, 4), (3, 4), (6, 4)):
        M = SimplicialChygraph([q], [k], [1.0])
        Ts = M.spinodal(lo=1e-2, hi=60.0)
        g = Chygraph([q], [float(k)], excess=[float(k - 1)])
        rho = float(np.max(np.abs(np.linalg.eigvals(
            g.branching_matrix(1.0 / Ts, 'simplicial')))))
        assert rho == pytest.approx(1.0, abs=1e-6)


# ---------------------------------------------------------------------------
# q-uniform: the tricritical cardinality and the ambivalent T_c
# ---------------------------------------------------------------------------

def _kind(q, k=4, J=1.0):
    M = SimplicialChygraph([q], [k], [J])
    Ts = M.spinodal(lo=1e-2, hi=60.0)
    above, _ = M.magnetisation(Ts * 1.02, u0=10.0)
    return Ts, ('discontinuous' if above > 1e-4 else 'continuous')


@pytest.mark.parametrize('q,want', [(2, 'continuous'), (3, 'continuous'),
                                    (4, 'continuous'), (5, 'continuous'),
                                    (6, 'discontinuous'), (8, 'discontinuous'),
                                    (16, 'discontinuous')])
def test_tricritical_cardinality_sits_between_five_and_six(q, want):
    """Son-Lee-Goh put the tricritical point at q = 4 in Bragg-Williams.  On a
    Bethe hyperlattice of chy-degree 4 it moves up: continuous through q = 5,
    discontinuous from q = 6.  Finite connectivity delays the discontinuity."""
    assert _kind(q)[1] == want


def test_bragg_williams_is_recovered_at_large_connectivity():
    """The tricritical cardinality is chy-degree dependent, and the shift to
    5--6 is confined to k = 3, 4.  By k = 6 it is back to the Bragg-Williams
    value of 4 and stays there, so the two treatments agree in the limit where
    the infinite-range one is exact."""
    def first_discontinuous(k):
        for q in range(2, 12):
            M = SimplicialChygraph([q], [k], [1.0])
            Ts = M.spinodal(lo=1e-3, hi=200.0)
            above, _ = M.magnetisation(Ts * 1.02, u0=10.0)
            if above > 1e-4:
                return q
        return None

    assert first_discontinuous(4) == 6
    for k in (6, 20, 100):
        assert first_discontinuous(k) == 5      # continuous through q = 4


def test_spinodal_is_nonmonotonic_in_cardinality():
    """Their 'ambivalent effect of group size', in the Bethe setting: the
    spinodal rises from q = 2 to q = 3 and falls thereafter."""
    T = [_kind(q)[0] for q in (2, 3, 4, 5, 6)]
    assert T[1] > T[0] and T[1] > T[2] > T[3] > T[4]


def test_spinodal_closed_forms():
    """With k = 4 the condition is 3 u' = 1, solvable in closed form.
    q = 2 and q = 4 both give e^{beta J} = 2, hence the same T = 1/ln 2."""
    assert _kind(2)[0] == pytest.approx(1 / np.log(2), rel=1e-6)
    assert _kind(4)[0] == pytest.approx(1 / np.log(2), rel=1e-6)
    assert _kind(3)[0] == pytest.approx(1 / np.log(1.8), rel=1e-6)


# ---------------------------------------------------------------------------
# The double transition, Fig. 4 of the paper
# ---------------------------------------------------------------------------

@pytest.fixture(scope='module')
def fig4():
    return SimplicialChygraph([2, 16], [4, 4], [0.7, 0.3])


def test_continuous_transition_is_driven_by_the_pairwise_layer(fig4):
    """The q = 16 layer contributes essentially nothing to the linear
    instability, so the spinodal is the pairwise one, 3 tanh(J_2/2T) = 1."""
    assert fig4.spinodal() == pytest.approx(0.35 / np.arctanh(1 / 3), rel=2e-3)


def test_intermediate_phase_is_pairwise_dominated(fig4):
    """Their m_2 >> m_q between the two transitions."""
    for T in (1.008, 1.000, 0.990):
        c, _ = fig4.components(T, u0=8.0)
        assert c[0] > 100 * c[1]


def test_hoi_becomes_comparable_below_the_discontinuous_transition(fig4):
    """'Only at temperatures below the discontinuous transition temperature
    does the contribution of HOIs become comparable.'"""
    for T in (0.95, 0.80, 0.60):
        c, _ = fig4.components(T, u0=8.0)
        assert 1.0 < c[0] / c[1] < 4.0


def test_the_transition_is_discontinuous_and_hysteretic(fig4):
    """Cooling and heating give different magnetisations in a window, which is
    the second, discontinuous transition."""
    cool, _ = fig4.magnetisation(0.980, u0=8.0)
    heat, _ = fig4.magnetisation(0.980, u0=1e-8)
    assert cool - heat > 0.2
    for T in (1.005, 0.950):                       # outside the window
        a, _ = fig4.magnetisation(T, u0=8.0)
        b, _ = fig4.magnetisation(T, u0=1e-8)
        assert abs(a - b) < 1e-3


def test_both_components_jump_at_the_discontinuous_transition(fig4):
    """'Both m_2 and m_q undergo both transitions at the same temperatures.'"""
    hi, _ = fig4.components(0.990, u0=8.0)
    lo, _ = fig4.components(0.975, u0=8.0)
    assert lo[0] > 1.5 * hi[0]                     # m_2 jumps
    assert lo[1] > 10 * hi[1]                      # m_q jumps much harder


# ---------------------------------------------------------------------------
# The discontinuous transition itself, not merely its signature
# ---------------------------------------------------------------------------

def test_continuous_case_has_no_coexistence_window():
    """q = 3 at chy-degree 4: one branch everywhere, so nothing to compare."""
    M = SimplicialChygraph([3], [4], [1.0])
    Ts = M.spinodal(lo=1e-2, hi=60.0)
    lo, hi = M.coexistence(Ts * 0.9, Ts * 2.2, n=100)
    assert np.isnan(lo) and np.isnan(hi)
    assert np.isnan(M.transition(scan=(Ts * 0.9, Ts * 2.2), n=100))


@pytest.mark.parametrize('q', [6, 8, 12])
def test_first_order_transition_lies_inside_its_metastability_limits(q):
    """T* < T_c < T**, with T_c from the free energies rather than from
    hysteresis.  The spinodal alone gives T*, not the transition."""
    M = SimplicialChygraph([q], [4], [1.0])
    Ts = M.spinodal(lo=1e-2, hi=60.0)
    lo, hi = M.coexistence(Ts * 0.9, Ts * 2.2, n=100)
    Tc = M.transition(scan=(Ts * 0.9, Ts * 2.2), n=100)
    assert lo < Tc < hi
    assert lo == pytest.approx(Ts, rel=1e-2)         # lower limit IS the spinodal
    a, _ = M.magnetisation(Tc + 1e-6, u0=1e-8)
    b, _ = M.magnetisation(Tc - 1e-6, u0=8.0)
    assert b - a > 0.5                                # a real jump


def test_fig4_double_transition_is_located_thermodynamically(fig4):
    """Both transitions of Fig. 4, with the second one placed by free energy:
    continuous at 1.0108, first-order at 0.9790 with Delta m = 0.286, inside
    metastability limits 0.9770 and 0.9810."""
    assert fig4.spinodal() == pytest.approx(1.0108, abs=1e-3)
    lo, hi = fig4.coexistence(0.96, 0.99, n=60)
    Tc = fig4.transition(scan=(0.96, 0.99), n=60)
    assert lo == pytest.approx(0.9770, abs=1e-3)
    assert hi == pytest.approx(0.9810, abs=1e-3)
    assert Tc == pytest.approx(0.9790, abs=1e-3)
    assert lo < Tc < hi
    a, _ = fig4.magnetisation(Tc + 1e-6, u0=1e-8)
    b, _ = fig4.magnetisation(Tc - 1e-6, u0=8.0)
    assert b - a == pytest.approx(0.286, abs=0.01)


def test_free_energy_paramagnetic_closed_form():
    """A simplicial pair is an Ising bond at half coupling, so the
    paramagnetic Bethe free energy is ln 2 + (k/2) ln[(1 + e^{beta J})/2]."""
    M = SimplicialChygraph([2], [4], [1.0])
    for T in (5.0, 3.0, 1.6):
        want = np.log(2) + 2 * np.log((1 + np.exp(1 / T)) / 2)
        assert M.minus_beta_f(T, u0=0.0) == pytest.approx(want, abs=1e-12)


def test_branch_gap_is_zero_outside_the_window(fig4):
    """Which is why a free-energy bracket has to sit inside it: outside, the
    difference is rounding noise of random sign.

    Sampled below the window only.  Just below the *continuous* transition at
    1.0108 the near-zero start converges slowly and a finite budget leaves a
    spurious gap -- critical slowing down, not coexistence.
    """
    for T in (0.9840, 0.9760, 0.9700):
        assert abs(fig4.branch_gap(T)) < 1e-9


# ---------------------------------------------------------------------------
# Against the published mean-field result
# ---------------------------------------------------------------------------

@pytest.mark.parametrize('q', [2, 3, 4, 5, 6, 8])
def test_bragg_williams_limit_is_recovered(q):
    """Son, Lee & Goh Eq. (8): T* = q(q-1)/2^(q-1), under their normalisation
    rho_q J_q = 1, i.e. J_q = q/k for k complexes of size q per vertex.

    Tested as a *rate* rather than a tolerance: the residual halves with each
    doubling of the chy-degree, so the limit is demonstrated rather than
    asserted at one k.  At fixed J the two treatments are not comparable at all.
    """
    want = q * (q - 1) / 2.0 ** (q - 1)
    res = []
    for k in (200.0, 400.0, 800.0):
        got = SimplicialChygraph([q], [k], [q / k]).spinodal(lo=1e-4, hi=200.0)
        res.append(abs(got - want))
    # 1/k everywhere except q = 4, where the leading correction vanishes and
    # the approach is 1/k^2 -- the mean-field tricritical cardinality.
    rate = 4.0 if q == 4 else 2.0
    assert res[0] / res[1] == pytest.approx(rate, abs=0.05)
    assert res[1] / res[2] == pytest.approx(rate, abs=0.05)
    assert res[-1] < 5e-3


def test_the_maximum_is_shared_at_q_three_and_four_in_the_limit():
    """Their Eq. (8) attains 3/2 at both q = 3 and q = 4.  At finite chy-degree
    the degeneracy is lifted, which is a connectivity effect and not a
    disagreement."""
    T = {q: SimplicialChygraph([q], [1600.0], [q / 1600.0]).spinodal(
        lo=1e-4, hi=200.0) for q in (2, 3, 4, 5)}
    assert T[3] == pytest.approx(1.5, abs=1e-3)
    assert T[4] == pytest.approx(1.5, abs=1e-3)
    assert T[2] < T[3] and T[5] < T[4]
