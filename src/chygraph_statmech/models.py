"""Named models: a chygraph plus a Hamiltonian.

Each constructor returns a :class:`~chygraph_statmech.stability.StabilityMatrix`
built from the same four moment tables the percolation calculation uses, with
the cavity weights of :mod:`chygraph_statmech.cavity` on the intra-complex
channel.
"""

from sympy import Symbol, eye, ones, symbols, zeros

from chygraph_statmech.cavity import (ising_edge_derivative,
                                      ising_triangle_derivative)
from chygraph_statmech.stability import StabilityMatrix


def _bare_graph_tables():
    """Layer 0 nodes, layer 1 links.  Cardinality 2, excess 1."""
    k, K, s, S = (zeros(2, 2) for _ in range(4))
    k[0, 1] = Symbol('k')          # mean degree
    K[0, 1] = Symbol('K')          # mean excess degree <kbar>
    s[1, 0] = 2                    # a link holds two nodes
    S[1, 0] = 1                    # ... one of them other than the entry node
    return k, K, s, S


def graph_percolation(p=None):
    """Bond percolation on a configuration-model graph.

    The weights are 1: occupation is folded into the moment tables by thinning,
    which is the percolation convention.  Reproduces
    ``chygraph.HypergraphPercolation(graph=True)`` at ``q = p``.
    """
    p = Symbol('p') if p is None else p
    k, K, s, S = _bare_graph_tables()
    s[1, 0] = 1 + p
    S[1, 0] = p
    return StabilityMatrix(k, K, s, S)


def graph_ising(beta=None, J=None, squared=False):
    """Ferromagnetic Ising on a configuration-model graph.

    ``squared=False`` gives the ferromagnetic (RS) instability, ``True`` the
    de Almeida-Thouless line.  Expected instability condition:

        <kbar> tanh(beta J) = 1        (ferro)
        <kbar> tanh^2(beta J) = 1      (AT)
    """
    k, K, s, S = _bare_graph_tables()
    u = ising_edge_derivative(beta, J)
    ws = zeros(2, 2)
    ws[1, 0] = u**2 if squared else u
    return StabilityMatrix(k, K, s, S, ws=_fill(ws))


def graph_with_triangles_ising(beta=None, J=None, squared=False):
    """Ferromagnetic Ising on a graph with links and triangles.

    Layer 0 nodes, layer 1 links, layer 2 triangles, following
    ``chygraph.percolation.GraphWithTriangles``.  The triangle is solved exactly
    inside the complex, so its cavity factor is ``t / (1 - t + t^2)`` rather than
    ``t``: the clustering that the ``{p_d, e_dd'}`` ensemble of Vazquez & Weigt
    cannot encode enters here as a different weight on layer 2.
    """
    k, K, s, S = (zeros(3, 3) for _ in range(4))
    k[0, 1], K[0, 1] = symbols('k_L K_L')
    k[0, 2], K[0, 2] = symbols('k_T K_T')
    s[1, 0], S[1, 0] = 2, 1        # a link holds 2 nodes, 1 other than entry
    s[2, 0], S[2, 0] = 3, 2        # a triangle holds 3 nodes, 2 other than entry
    uL = ising_edge_derivative(beta, J)
    uT = ising_triangle_derivative(beta, J)
    ws = zeros(3, 3)
    ws[1, 0] = uL**2 if squared else uL
    ws[2, 0] = uT**2 if squared else uT
    return StabilityMatrix(k, K, s, S, ws=_fill(ws))


def _fill(ws):
    """Zero entries of a weight table mean 'unused'; make them 1 so the
    elementwise multiply leaves the (already zero) moment untouched."""
    L = ws.shape[0]
    out = ones(L, L)
    for i in range(L):
        for j in range(L):
            if ws[i, j] != 0:
                out[i, j] = ws[i, j]
    return out
