"""Regression tests pinning the supplementary equations from

    Alexei Vazquez, "Percolation in higher order networks via mapping to chygraphs"
    https://doi.org/10.1093/comnet/cnae047

Each ``theta`` / ``Lambda`` here is compared against its published closed form
via ``simplify(... ) == 0``, so a change in the symbolic engine (or an
accidental edit to the model definitions) fails loudly. The ``Lambda`` tests
also pin the eigenvalue *index*, guarding against sympy changing the order in
which ``Matrix.eigenvals()`` returns its keys.
"""

from sympy import symbols, sqrt, simplify, Float

from chygraph import (
    HypergraphPercolation,
    MultiplexHypergraph,
    GraphWithTriangles,
)


def _eq(a, b):
    """True iff two symbolic expressions are equal after simplification."""
    return simplify(a - b) == 0


# ---------------------------------------------------------------------------
# S.VI Graphs and hypergraphs
# ---------------------------------------------------------------------------

def test_hypergraph_theta_and_lambda():
    p, q, k, K, c, C = symbols('p q k K c C')  # noqa: F841 (documenting the symbols)
    H = HypergraphPercolation()
    assert _eq(H.A.theta(), C * K * p * q - 1)
    assert _eq(H.A.eigenvals()[1], sqrt(C * K * p * q) - 1)


def test_graph_theta_and_lambda():
    p, q, K = symbols('p q K')
    G = HypergraphPercolation(graph=True)
    assert _eq(G.A.theta(), K * p * q - 1)
    assert _eq(G.A.eigenvals()[1], sqrt(K * p * q) - 1)


# ---------------------------------------------------------------------------
# S.V Multiplex hypergraphs
# ---------------------------------------------------------------------------

def test_multiplex_two_types_theta():
    K01, K02, S10, S20, k01, k02 = symbols('K_01 K_02 S_10 S_20 k_01 k_02')
    MH = MultiplexHypergraph(number_of_types=2)
    expected = K01 * S10 + K02 * S20 + S10 * S20 * (k01 * k02 - K01 * K02) - 1
    assert _eq(MH.A.theta(), expected)


# ---------------------------------------------------------------------------
# S.V.B Network motifs: graph with links and triangles
# ---------------------------------------------------------------------------

def test_graph_with_triangles_theta():
    q, k_L, k_T, K_L, K_T = symbols('q k_L k_T K_L K_T')
    S10 = q
    S20 = 2 * q * (1 + q - q**2)
    expected = K_L * S10 + K_T * S20 - 1 + S10 * S20 * (k_L * k_T - K_L * K_T)
    NM = GraphWithTriangles()
    assert _eq(NM.A.theta(), expected)


def test_graph_with_triangles_is_exact():
    """The (3/2) coefficient must stay a Rational, never a Python float."""
    NM = GraphWithTriangles()
    assert not NM.A.theta().atoms(Float)
    assert not NM.A.eigenvals()[3].atoms(Float)


# ---------------------------------------------------------------------------
# Caching behaviour
# ---------------------------------------------------------------------------

def test_theta_and_eigenvals_are_cached():
    H = HypergraphPercolation()
    assert H.A.theta() is H.A.theta()
    assert H.A.eigenvals() is H.A.eigenvals()
