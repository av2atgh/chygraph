"""Tests for chygraphs with statistically dependent layers.

The published threshold tensor uses the unconditional <kappa>_{nk} in the slots
with k != m.  That is the independent-layer special case; in general those slots
carry the inclusion-biased moment <kappa_m kappa_k>/<kappa_m>.  These tests pin
both the reduction to the published tensor and the correlated generalisation,
the latter against Monte Carlo reference values.
"""

from sympy import symbols, exp, simplify, zeros, Rational, nsolve

from chygraph.giant import hypergraph_giant, graph_with_triangles_giant
from chygraph.joint import JointGiantComponent
from chygraph.amplitude import CriticalAmplitude


q, kL, kT, k, c, p = symbols('q k_L k_T k c p')
HALF = Rational(1, 2)


def _triangle_parts():
    """Intra-complex generating functions of the link and triangle layers."""
    tri_bar = lambda y: ((3 * q**2 - 2 * q**3) * y[0]**2
                         + 2 * q * (1 - q)**2 * y[0]
                         + (1 - q)**3 + q * (1 - q)**2)
    Gbar = [[None] * 3 for _ in range(3)]
    Gbar[1][0] = lambda y: (1 - q) + q * y[0]
    Gbar[2][0] = tri_bar
    G = [None,
         lambda y: (1 - q) * y[0] + q * y[0]**2,
         lambda y: y[0] * tri_bar(y)]
    return G, Gbar


# joint chy-degree generating functions with identical marginals:
# k_| in {1, 3} and k_T in {0, 2}, each with probability 1/2
PHI = {
    'correlated':  lambda x: HALF * x[1] + HALF * x[1]**3 * x[2]**2,
    'anticorr':    lambda x: HALF * x[1] * x[2]**2 + HALF * x[1]**3,
    'independent': lambda x: (HALF * x[1] + HALF * x[1]**3) * (HALF + HALF * x[2]**2),
}


def _model(name):
    G, Gbar = _triangle_parts()
    return JointGiantComponent(Phi=[PHI[name], None, None], G=G, Gbar=Gbar)


# ---------------------------------------------------------------------------
# Reduction to the factorised construction
# ---------------------------------------------------------------------------

def test_hypergraph_reduces_to_the_published_tensor():
    J = JointGiantComponent(Phi=[lambda x: exp(k * (x[1] - 1)), None],
                            G=[None, lambda y: exp(c * (y[0] - 1))],
                            occupation=[p, q])
    assert simplify(J.A() - hypergraph_giant().A()) == zeros(8, 8)
    assert J.layers_independent()
    for s in ({'k': 3, 'c': 3, 'p': 0.6, 'q': 0.8},
              {'k': 2, 'c': 4, 'p': 0.9, 'q': 0.5}):
        assert abs(J.node_fraction(s) - hypergraph_giant().node_fraction(s)) < 1e-12


def test_triangles_reduce_to_the_published_tensor():
    G, Gbar = _triangle_parts()
    J = JointGiantComponent(
        Phi=[lambda x: exp(kL * (x[1] - 1)) * exp(kT * (x[2] - 1)), None, None],
        G=G, Gbar=Gbar)
    assert simplify(J.A() - graph_with_triangles_giant().A()) == zeros(18, 18)
    for s in ({'k_L': 1, 'k_T': 0.5, 'q': 0.7}, {'k_L': 3, 'k_T': 2, 'q': 0.35}):
        assert abs(J.node_fraction(s) - graph_with_triangles_giant().node_fraction(s)) < 1e-12


def test_supplied_excess_is_required_for_motifs():
    """A triangle's excess component size is not the size-biased cardinality, so
    deriving it by differentiation gives the wrong map."""
    G, Gbar = _triangle_parts()
    Phi = [lambda x: exp(kL * (x[1] - 1)) * exp(kT * (x[2] - 1)), None, None]
    derived = JointGiantComponent(Phi=Phi, G=G)              # no Gbar
    supplied = JointGiantComponent(Phi=Phi, G=G, Gbar=Gbar)
    s = {'k_L': 1, 'k_T': 0.5, 'q': 0.7}
    assert abs(supplied.node_fraction(s) - graph_with_triangles_giant().node_fraction(s)) < 1e-12
    assert abs(derived.node_fraction(s) - graph_with_triangles_giant().node_fraction(s)) > 1e-3


# ---------------------------------------------------------------------------
# The correlated generalisation
# ---------------------------------------------------------------------------

def test_inclusion_biased_moment():
    """<kappabar^(1)>_{0,2} = <kappa_1 kappa_2>/<kappa_1>, which equals
    <kappa>_{0,2} only when the layers are independent."""
    for name, want in (('correlated', Rational(3, 2)),
                       ('anticorr', Rational(1, 2)),
                       ('independent', 1)):
        M = _model(name)
        assert M.kappa(0, 1) == 2 and M.kappa(0, 2) == 1     # marginals shared
        assert simplify(M.kappa_bar(1, 0, 2) - want) == 0
        assert M.layers_independent() == (name == 'independent')


def test_identical_marginals_give_different_thresholds():
    qc = {}
    for name in PHI:
        qc[name] = float(nsolve(_model(name).theta(), q, 0.3))
    assert abs(qc['correlated'] - 0.193884) < 1e-5
    assert abs(qc['anticorr'] - 0.316110) < 1e-5
    assert abs(qc['independent'] - 0.240915) < 1e-5


def test_correlated_order_parameter_matches_simulation():
    """Reference values from configuration model simulation with n = 3e5 nodes
    and the prescribed joint degree pairs; agreement is within one standard
    deviation.  The marginals-only prediction is the 'independent' column and
    differs from the other two."""
    ref = {
        ('correlated', 0.3): 0.46601, ('correlated', 0.5): 0.67262,
        ('correlated', 1.0): 0.87500,
        ('anticorr', 0.5): 0.79267, ('anticorr', 1.0): 1.00000,
        ('independent', 0.3): 0.35811, ('independent', 0.5): 0.74554,
        ('independent', 1.0): 0.96656,
    }
    for (name, qv), want in ref.items():
        got = _model(name).node_fraction({'q': qv})
        assert abs(got - want) < 1e-4, (name, qv, got)
    # the correlated models are genuinely distinguishable from the marginals
    for name in ('correlated', 'anticorr'):
        assert abs(ref[(name, 1.0)] - ref[('independent', 1.0)]) > 0.03


def test_amplitude_composes_with_joint_models():
    M = _model('correlated')
    C = CriticalAmplitude(M)
    assert all(C.verify().values())
    qc = float(nsolve(C.Lambda(), q, 0.2))
    B = float(C.amplitude(0).subs(q, qc))
    S = M.node_fraction({'q': qc + 0.001}, tol=1e-16, maxiter=4_000_000)
    lam = float(C.perron_root().subs(q, qc + 0.001))
    assert abs(S / (lam - 1) - B) < 5e-2
