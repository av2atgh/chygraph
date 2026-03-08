"""Core symbolic percolation classes for Chygraphs.

References:
    Alexei Vazquez, "Percolation in higher order networks via mapping to chygraphs"
    https://doi.org/10.1093/comnet/cnae047
    https://arxiv.org/abs/2308.00987
"""

from sympy import zeros, symbols


# ---------------------------------------------------------------------------
# Core matrix class
# ---------------------------------------------------------------------------

class PercolationMatrix:
    """Symbolic matrix for percolation calculations in Chygraphs."""

    def __init__(self, k, K, s, S):
        L = k.shape[0]
        Lv = L ** 2
        A = zeros(2 * Lv, 2 * Lv)

        def delta(i, j):
            return int(i == j)

        for i in range(2):
            for j in range(2):
                for m in range(L):
                    for l in range(L):
                        for n in range(L):
                            for o in range(L):
                                row = i * Lv + m * L + l
                                col = j * Lv + n * L + o
                                if i == 0 and j == 0:
                                    A[row, col] = (
                                        delta(m, n) * delta(l, o) - k[n, o] * delta(n, l)
                                    )
                                elif i == 0 and j == 1:
                                    A[row, col] = (
                                        -(S[n, o] * delta(m, o) + (1 - delta(m, o)) * s[n, o]) * delta(n, l)
                                    )
                                elif i == 1 and j == 0:
                                    A[row, col] = (
                                        -(K[n, o] * delta(m, o) + (1 - delta(m, o)) * k[n, o]) * delta(n, l)
                                    )
                                elif i == 1 and j == 1:
                                    A[row, col] = (
                                        delta(m, n) * delta(l, o) - s[n, o] * delta(n, l)
                                    )
        self.A = A

    def theta(self):
        """Return the percolation threshold expression det(-A)."""
        return -self.A.det()

    def eigenvals(self):
        """Return the eigenvalues of -A as a list."""
        return list((-self.A).eigenvals().keys())


# Backward-compatibility alias
vec2A = PercolationMatrix


# ---------------------------------------------------------------------------
# Network model classes
# ---------------------------------------------------------------------------

class HypergraphPercolation:
    """Symbolic percolation calculations for Hypergraphs."""

    def __init__(self, graph=False):
        L = 2
        k, K, s, S = zeros(L, L), zeros(L, L), zeros(L, L), zeros(L, L)
        p, q = symbols('p q')
        k[0, 1] = p * symbols('k')
        K[0, 1] = p * symbols('K')
        s[1, 0] = 1 + q if graph else q * symbols('c')
        S[1, 0] = q if graph else q * symbols('C')
        self.A = PercolationMatrix(k, K, s, S)


class MultiplexHypergraph:
    """Symbolic percolation calculations for Multiplex Hypergraphs."""

    def __init__(self, number_of_types=2, graph=False, poisson=False):
        L = number_of_types + 1
        k, K, s, S = zeros(L, L), zeros(L, L), zeros(L, L), zeros(L, L)
        for l in range(1, L):
            k[0, l] = symbols(f'k_0{l}')
            K[0, l] = symbols(f'K_0{l}')
            s[l, 0] = 2 if graph else symbols(f's_{l}0')
            S[l, 0] = 1 if graph else symbols(f'S_{l}0')
        if poisson:
            K = k
        if not graph and poisson:
            S = s
        self.A = PercolationMatrix(k, K, s, S)


class InteractingHypergraphs:
    """Symbolic percolation calculations for Interacting Hypergraphs.

    Args:
        g:       Number of graphs/hypergraphs.
        graph:   True for graphs, False for hypergraphs.
        poisson: True for Poisson degree/cardinality distributions.
    """

    def __init__(self, g=2, graph=False, poisson=False):
        L = int(g + (g * (g + 1)) / 2)
        k, K, s, S = zeros(L, L), zeros(L, L), zeros(L, L), zeros(L, L)
        for l in range(g):
            for m in range(l, g):
                i = g + l * g - int((l * (l - 1)) / 2) + m - l
                if l == m:
                    k[l, i] = symbols(f'k_{l}{i}')
                    K[l, i] = symbols(f'K_{l}{i}')
                    s[i, l] = 2 if graph else symbols(f'c_{i}{l}')
                    S[i, l] = 1 if graph else symbols(f'C_{i}{l}')
                else:
                    k[l, i] = symbols(f'k_{l}{i}')
                    k[m, i] = symbols(f'k_{m}{i}')
                    K[l, i] = symbols(f'K_{l}{i}')
                    K[m, i] = symbols(f'K_{m}{i}')
                    s[i, l] = 1 if graph else symbols(f'c_{i}{l}')
                    s[i, m] = 1 if graph else symbols(f'c_{i}{m}')
                    S[i, l] = 0 if graph else symbols(f'C_{i}{l}')
                    S[i, m] = 0 if graph else symbols(f'C_{i}{m}')
        if poisson:
            K = k
        if poisson and not graph:
            S = s
        self.A = PercolationMatrix(k, K, s, S)


class GraphWithTriangles:
    """Percolation in a graph with links and triangles.

    The graph is mapped to a multiplex hypergraph with two layers:
    one for links (l=1) and one for triangles (l=2).

    Args:
        poisson: True for Poisson degree distribution (sets K = k).
    """

    def __init__(self, poisson=False):
        L = 3
        k, K, s, S = zeros(L, L), zeros(L, L), zeros(L, L), zeros(L, L)
        q = symbols('q')
        k[0, 1] = symbols('k_L')
        k[0, 2] = symbols('k_T')
        s[1, 0] = 1 + q
        s[2, 0] = 3 * (q**3 + 3*q**2*(1 - q)) + (3/2)*3*q*(1 - q)**2 + (1 - q)**3
        S[1, 0] = q
        S[2, 0] = 2 * q * (1 + q - q**2)
        if poisson:
            K = k
        else:
            K[0, 1] = symbols('K_L')
            K[0, 2] = symbols('K_T')
        self.A = PercolationMatrix(k, K, s, S)


