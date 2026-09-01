"""Every derivation written out in the Supplemental Material, checked.

The referee's closing request was that each quantitative claim be independently
checked, and the paper's own record is that the claims which turned out wrong
were the ones reported as a number with no derivation behind them.  So each
step restored in the book's calculation boxes gets a test here: the closed form against
the code that the manuscript actually runs, or against an independent
enumeration, whichever is the stronger check.

Section numbers below are those of the Supplemental Material.
"""

import numpy as np
import pytest
import sympy as sp

from chygraph_statmech import Chygraph
from chygraph_statmech import freeenergy, ising
from chygraph_statmech import simplicial as simp
from chygraph_statmech.hittingset import layer_symbols
from chygraph_statmech.softfield import (regular_density, regular_entropy,
                                         regular_field)


# ---------------------------------------------------------------------------
# S2: the branching matrix
# ---------------------------------------------------------------------------

def test_size_biased_excess_cardinality():
    """S2 A.  ``<sbar> = <c^2>/<c> - 1``, and it is not ``<c> - 1``.

    Following an inclusion into a complex samples it in proportion to its
    cardinality, so ``<sbar> = sum_c (c p_c / <c>)(c - 1)``, which telescopes.
    """
    for dist in ({3: 0.5, 5: 0.5}, {2: 0.7, 4: 0.2, 9: 0.1}, {4: 1.0}):
        c = np.array(sorted(dist), dtype=float)
        p = np.array([dist[int(x)] for x in c])
        g = Chygraph([dist], [1.0])
        want = (p * c * c).sum() / (p * c).sum() - 1.0
        assert g.excess_cardinality_layer(0) == pytest.approx(want, rel=1e-12)
    # the trap the referee named: the two differ
    g = Chygraph([{3: 0.5, 5: 0.5}], [1.0])
    assert g.excess_cardinality_layer(0) == pytest.approx(3.25)
    assert g.c[0] - 1 == pytest.approx(3.0)


def test_links_and_triangles_determinant():
    """S2 B.  ``det(I - B) = 0`` for two layers is Eq. (15), expanded by hand.

    The cross term carries ``<kappa>_L<kappa>_T - <kbar>_L<kbar>_T`` and so
    vanishes for Poisson layers, where the excess equals the mean.
    """
    for kL, kT, eL, eT, bJ in ((2.0, 1.0, 2.0, 1.0, 0.4),
                               (3.0, 2.0, 2.5, 1.5, 0.25),
                               (1.5, 0.5, 0.9, 0.2, 0.6)):
        g = Chygraph([2, 3], [kL, kT], excess=[eL, eT])
        det = float(np.linalg.det(np.eye(2) - g.branching_matrix(bJ)))
        uL = ising.clique_derivative(2, bJ)
        uT = ising.clique_derivative(3, bJ)
        eq15 = 1.0 - (eL * uL + 2 * eT * uT
                      + 2 * uL * uT * (kL * kT - eL * eT))
        assert det == pytest.approx(eq15, abs=1e-13)
    # Poisson: the cross term is gone and the condition is a plain sum
    kL, kT, bJ = 2.0, 1.0, 0.4
    g = Chygraph([2, 3], [kL, kT])
    det = float(np.linalg.det(np.eye(2) - g.branching_matrix(bJ)))
    uL = ising.clique_derivative(2, bJ)
    uT = ising.clique_derivative(3, bJ)
    assert det == pytest.approx(1.0 - kL * uL - 2 * kT * uT, abs=1e-13)


# ---------------------------------------------------------------------------
# S3: the Ising model
# ---------------------------------------------------------------------------

def test_clique_transmission_closed_forms():
    """S3 A.  Eq. (17): ``t``, ``t/(1-t+t^2)``, ``(t+t^3)/(1-2t+3t^2)``.

    The enumeration of Eq. (5) against the closed forms, and the statement that
    a triangle transmits more per neighbour than an edge.
    """
    for bJ in (0.05, 0.3, 0.9, 2.0):
        t = np.tanh(bJ)
        closed = {2: t,
                  3: t / (1 - t + t ** 2),
                  4: (t + t ** 3) / (1 - 2 * t + 3 * t ** 2)}
        for c, want in closed.items():
            assert ising.clique_derivative(c, bJ) == pytest.approx(want, rel=1e-12)
        assert closed[3] > closed[2]              # more per traversal


def test_triangle_transmission_symbolically():
    """S3 A.  ``u'_3 = t/(1-t+t^2)`` derived, not fitted.

    ``u' = (1/2)[<sigma_1>_+ - <sigma_1>_-]`` inside the triangle at zero field
    is ``(e^{4a}-1)/(e^{4a}+3)``; substituting ``e^{2a} = (1+t)/(1-t)`` gives
    the closed form.
    """
    a, t = sp.symbols('a t', positive=True)
    up = (sp.exp(4 * a) - 1) / (sp.exp(4 * a) + 3)
    sub = sp.simplify(up.subs(sp.exp(a), sp.sqrt((1 + t) / (1 - t))))
    assert sp.simplify(sp.together(sub - t / (1 - t + t ** 2))) == 0


def test_bethe_free_energy_collapse_at_zero_field():
    """S3 B.  Eq. (14) is Eq. (13) with every cavity field zero.

    ``Z_i = 2`` and ``ln Z_a = ln Z_c``, so the node term contributes
    ``(1 - sum_l <k_l>) ln 2`` and the complex terms ``<k_l>/c_l ln Z_c``,
    which regroups into Eq. (14).
    """
    for cs, ks, bJ in (([2], [3.0], 0.2), ([2, 3], [1.5, 0.8], 0.35),
                       ([4], [2.0], 0.15), ([3, 5], [1.0, 0.4], 0.28)):
        want = np.log(2.0)
        for c, k in zip(cs, ks):
            m = int(c)
            bits = ((np.arange(2 ** m)[:, None] >> np.arange(m)) & 1)
            cfg = 1.0 - 2.0 * bits
            S = cfg.sum(axis=1)
            Zc = np.exp(bJ * (S ** 2 - m) / 2.0).sum()
            want += (k / m) * np.log(Zc / 2 ** m)
        assert freeenergy.paramagnetic(cs, ks, bJ) == pytest.approx(want, rel=1e-13)
    # the textbook graph limit
    assert freeenergy.paramagnetic([2], [3.0], 0.4) == pytest.approx(
        freeenergy.graph_paramagnetic(3.0, 0.4), rel=1e-13)


# ---------------------------------------------------------------------------
# S4: the simplicial interaction
# ---------------------------------------------------------------------------

def test_simplicial_emitted_field_by_enumeration():
    """S4 A.  Eq. (20) against a brute-force sum over the complex interior.

    ``Z(sigma_0) = (2 cosh h)^{q-1} + (e^a - 1) e^{+-h(q-1)}`` says the
    unanimity rule adds one term to the free sum; the enumeration confirms it.
    """
    import itertools
    for q in (2, 3, 4, 5):
        for a, h in ((0.4, 0.3), (1.3, -0.7), (2.0, 0.0)):
            for s0 in (1, -1):
                brute = 0.0
                for rest in itertools.product((1, -1), repeat=q - 1):
                    agree = all(s == s0 for s in rest)
                    brute += np.exp(a * agree + h * sum(rest))
                closed = (2 * np.cosh(h)) ** (q - 1) + \
                    (np.exp(a) - 1) * np.exp(s0 * h * (q - 1))
                assert brute == pytest.approx(closed, rel=1e-12)
            u = 0.5 * np.log(
                ((2 * np.cosh(h)) ** (q - 1) + (np.exp(a) - 1) * np.exp(h * (q - 1)))
                / ((2 * np.cosh(h)) ** (q - 1) + (np.exp(a) - 1) * np.exp(-h * (q - 1))))
            assert simp.emitted(q, a, h) == pytest.approx(u, rel=1e-12)


def test_simplicial_uprime_and_its_two_conventions():
    """S4 A.  Eq. (21), per neighbour, and the ``q-1`` that must not be doubled.

    Differentiating with respect to one member's field gives
    ``(e^a-1)/(2^{q-1}+e^a-1)``; with respect to a field common to all ``q-1``
    others it gives ``q-1`` times that.  Using the second in Eq. (8) counts the
    multiplicity twice.
    """
    eps = 1e-6
    for q in (2, 3, 5, 8):
        for a in (0.4, 1.3, 3.0):
            want = (np.exp(a) - 1) / (2 ** (q - 1) + np.exp(a) - 1)
            assert simp.uprime(q, a) == pytest.approx(want, rel=1e-12)
            # common-field derivative, numerically
            common = (simp.emitted(q, a, eps) - simp.emitted(q, a, -eps)) / (2 * eps)
            assert common == pytest.approx((q - 1) * want, rel=1e-6)
    # q = 2 is an ordinary bond of half the coupling
    for a in (0.3, 1.1, 2.5):
        assert simp.uprime(2, a) == pytest.approx(np.tanh(a / 2), rel=1e-12)


def test_simplicial_bragg_williams_limit_and_its_residual():
    """S4 B.  Eq. (22), and *why* ``q = 4`` approaches it as ``1/k^2``.

    At the normalisation ``J_q = q/k`` of Ref. [SLG26], expanding
    ``(k-1)(q-1)u' = 1`` in ``1/k`` gives

        T = T* - C/k + O(1/k^2),   T* = q(q-1)/2^{q-1},
        C  = q(2q - 2^{q-1}) / 2^q,

    and ``C`` vanishes at ``q = 4`` because ``2q = 2^{q-1}`` there --- the one
    cardinality at which the leading finite-connectivity correction is absent.
    The manuscript quotes the residual at ``q = 3`` without this coefficient;
    it is ``3/4``, which reproduces the quoted numbers.
    """
    for q in (2, 3, 4, 5, 6, 8):
        Tstar = q * (q - 1) / 2.0 ** (q - 1)
        C = q * (2 * q - 2.0 ** (q - 1)) / 2.0 ** q
        res = {k: simp.SimplicialChygraph([q], [k], [q / k]).spinodal() - Tstar
               for k in (200, 800)}
        if q == 4:
            assert C == pytest.approx(0.0, abs=1e-12)
            # no 1/k term left: the residual falls four times faster per doubling
            assert abs(res[800]) < abs(res[200]) / 8
        else:
            for k in (200, 800):
                assert res[k] == pytest.approx(-C / k, rel=0.02)
    # the coefficient behind the numbers quoted for q = 3
    assert 3 * (6 - 4) / 8 == pytest.approx(0.75)
    for k, quoted in ((50, 1.5e-2), (200, 3.8e-3), (800, 9.4e-4)):
        assert 0.75 / k == pytest.approx(quoted, rel=0.02)


# ---------------------------------------------------------------------------
# S5: hitting set
# ---------------------------------------------------------------------------

def test_regular_hitting_set_closed_form():
    """S5 A.  Eq. (32) from Eq. (31), and ``rho = 1/K`` from the full field.

    As ``mu -> inf`` the complex step becomes ``v = -ln(K-1) - h``; substituting
    into ``h = -mu + (L-1)v`` gives ``L h = -mu - (L-1) ln(K-1)``.  The full
    field is ``H = h + v = -ln(K-1)``, independent of ``mu`` and of ``L``, so
    ``rho = 1/(1 + e^{-H}) = 1/K``.
    """
    for L, K in ((1, 3), (2, 3), (3, 4), (2, 6), (4, 6), (6, 12)):
        for mu in (40.0, 60.0, 120.0):
            h = regular_field(L, K, mu)
            v = -np.log(K - 1) - h                  # the mu -> inf complex step
            assert h == pytest.approx(-mu + (L - 1) * v, rel=1e-12)
            H = h + v
            assert H == pytest.approx(-np.log(K - 1), rel=1e-12)
            assert 1.0 / (1.0 + np.exp(-H)) == pytest.approx(regular_density(K),
                                                             rel=1e-12)
        # h carries a *fraction* of mu, which no integer ansatz can represent
        assert (regular_field(L, K, 60.0) + 60.0 / L) == pytest.approx(
            -(L - 1) / L * np.log(K - 1), rel=1e-12)


def test_regular_hitting_set_entropy_reduction():
    """S5 B.  Eq. (34) is Eq. (33) evaluated on Eq. (32).

    With ``ln Z_a -> ln K + h``, ``Z_i = K/(K-1)`` and ``rho = 1/K``, the two
    ``O(mu)`` pieces cancel and what is left regroups with
    ``a = (L-1)(K-1)`` into Eq. (34).
    """
    for L, K in ((1, 3), (2, 3), (3, 4), (2, 6), (4, 6), (6, 12)):
        mu = 60.0
        h = regular_field(L, K, mu)
        lnZa = np.log(K) + h                        # (1+e^h)^K - 1 ~ K e^h
        lnZi = np.log(K / (K - 1.0))
        s = (L / K) * lnZa + (1 - L) * lnZi + mu / K
        assert s == pytest.approx(regular_entropy(L, K), abs=1e-9)
        a = (L - 1) * (K - 1)
        assert regular_entropy(L, K) == pytest.approx(
            (a * np.log(K - 1) - (a - 1) * np.log(K)) / K, rel=1e-13)
    # where it turns negative the replica-symmetric answer is an underestimate
    assert regular_entropy(1, 3) > 0
    assert regular_entropy(4, 6) < 0 and regular_entropy(6, 12) < 0


def test_hard_field_threshold_elimination():
    """S5 C.  ``<k>(c-1) = e``, by the elimination the manuscript states.

    ``sigma = e^{-<k> sigma^{c-1}}``; with ``y = <k> sigma^{c-1}`` the
    instability ``|F'| = <k>(c-1)sigma^{c-1} = 1`` is ``(c-1)y = 1``, and
    ``<k> = y e^{(c-1)y}`` then gives ``<k> = e/(c-1)``.
    """
    from scipy.optimize import brentq
    for c in range(2, 21):
        k = np.e / (c - 1)
        sigma = brentq(lambda s: s - np.exp(-k * s ** (c - 1)), 1e-12, 1.0)
        y = k * sigma ** (c - 1)
        assert (c - 1) * y == pytest.approx(1.0, rel=1e-10)      # instability
        assert k * (c - 1) * sigma ** (c - 1) == pytest.approx(1.0, rel=1e-10)
        assert k * (c - 1) == pytest.approx(np.e, rel=1e-12)


# ---------------------------------------------------------------------------
# S6: core percolation
# ---------------------------------------------------------------------------

def test_no_core_free_branch_above_cardinality_two():
    """S6 B.  ``gamma = 0`` forces ``lambda = 0`` for ``c >= 3``, in closed form.

    ``gamma = 0`` needs ``zeta = 1 - kappa`` with ``lambda = 1 - delta``, i.e.

        f(delta) = (c-1) delta^{c-2} - (c-2) delta^{c-1} - 1 = 0,

    and ``f`` factors as ``-(1-delta)^2 sum_{j=0}^{c-3} (j+1) delta^j``.  Every
    coefficient of that sum is positive, so ``f < 0`` on ``[0,1)`` and
    ``delta = 1`` is the only root: no core-free branch survives.  At ``c = 2``
    the sum is empty, ``f`` vanishes identically, and the branch exists.
    """
    d = sp.Symbol('d')
    for c in range(2, 15):
        f = sp.expand((c - 1) * d ** (c - 2) - (c - 2) * d ** (c - 1) - 1)
        poly = sum((j + 1) * d ** j for j in range(c - 2))
        assert sp.simplify(f + (1 - d) ** 2 * poly) == 0
        if c == 2:
            assert f == 0                        # identity: the branch exists
        else:
            roots = [r for r in sp.solve(sp.Eq(f, 0), d)
                     if r.is_real and 0 <= r <= 1]
            assert roots == [1]                  # delta = 1, so lambda = 0


def test_core_map_reduces_to_the_graph_at_cardinality_two():
    """S6 A.  ``zeta = delta^{c-1}`` and ``kappa = (c-1)delta^{c-2}lambda`` are
    the graph recursion at ``c = 2`` and nothing else."""
    from chygraph_statmech.core import CorePercolation
    from chygraph_statmech.hittingset import poisson_phi
    m = CorePercolation([2], poisson_phi([1.0]))
    for lam, delta in ((0.3, 0.4), (0.0, 1.0), (0.62, 0.11)):
        assert m.zeta(np.array([delta]))[0] == pytest.approx(delta)
        assert m.kappa(np.array([lam]), np.array([delta]))[0] == pytest.approx(lam)
