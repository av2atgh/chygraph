"""The reweighted stability tensor (WP1).

``chygraph.percolation.PercolationMatrix`` builds ``A = I - J`` from four
``L x L`` moment tables ``k``, ``K``, ``s``, ``S``.  Its entries are *first
moments* because in percolation the message derivative is identically the
occupation probability, and thinning has already folded that into
``<kappa>_{0l} = p<k>``.  Nothing else multiplies a branch.

For a general Hamiltonian the linearised cavity recursion is

    delta h^(i|j) = sum_{k != j} u'(h^(k|i)) delta h^(k|i)

so every branch carries a factor ``u'``.  At the trivial fixed point that factor
is one number per traversal type, and the reweight is therefore a **pure
elementwise multiply on the moment tables** — the ``2L^2`` index structure, the
block layout, ``theta = -det(A)`` and ``Lambda = max eig(-A)`` are all unchanged.
That claim is the content of prediction 1 in the README, and it is what makes
this module fourteen lines of arithmetic rather than a rewrite.

Two instabilities, same tensor, different weight:

    <kappa> -> <kappa> <u'>       ferromagnetic / RS instability
    <kappa> -> <kappa> <u'^2>     de Almeida-Thouless / spin-glass line

Where the weights go.  A chygraph has two kinds of step.  Going *up* through a
chy-degree, a complex reports its own cavity field to a complex that includes
it; that is a pure relay and carries weight 1.  Going *down* through the
intra-complex hypergraph, the interaction inside the complex transforms the
field, and the weight is the ``u'`` of :mod:`chygraph_statmech.cavity`.  So
``wkappa = 1`` and ``ws`` carries the physics.  For percolation both are 1 and
``StabilityMatrix`` reduces to ``PercolationMatrix`` identically.

Caveat on factorisation.  The correct entry is ``<sbar * u'>``, the joint
average over the complex's internal structure.  Writing it as
``<sbar> * <u'>`` assumes ``u'`` is the same for every complex in the layer.
That is not a restriction so much as a definition of what a layer is: a class of
complexes with the same statistical properties.  When ``u'`` depends on
cardinality, split the layer by cardinality — the same device
``chygraph.applications.correlated_cardinality_hypergraph`` already uses.
"""

from sympy import Matrix, ones, sympify

from chygraph.percolation import PercolationMatrix


def _elementwise(table, weights):
    """Hadamard product, tolerating ``None`` for 'all ones'."""
    if weights is None:
        return table
    L = table.shape[0]
    return Matrix(L, L, lambda i, j: sympify(table[i, j]) * sympify(weights[i, j]))


def reweight(k, K, s, S, wkappa=None, ws=None):
    """Apply ``wkappa`` to the inclusion channel and ``ws`` to the intra-complex one."""
    return (_elementwise(k, wkappa), _elementwise(K, wkappa),
            _elementwise(s, ws), _elementwise(S, ws))


class StabilityMatrix(PercolationMatrix):
    """``A = I - J`` for a statistical-mechanics model on a chygraph.

    Args:
        k, K, s, S: the four moment tables, exactly as in
            ``chygraph.percolation.PercolationMatrix``.
        wkappa: ``L x L`` weights on the inclusion channel.  ``None`` means 1.
        ws: ``L x L`` weights on the intra-complex channel.  ``None`` means 1.

    With both weights ``None`` this *is* ``PercolationMatrix``, entry for entry.
    """

    def __init__(self, k, K, s, S, wkappa=None, ws=None):
        super().__init__(*reweight(k, K, s, S, wkappa, ws))

    # -- diagnostics --------------------------------------------------------

    def perron_root(self, probe=None):
        """Largest eigenvalue of ``J = I - A``; the instability is at 1."""
        eigs = [e for e in (-self.A).eigenvals() if e != 0]
        vals = [sympify(e) + 1 for e in eigs]
        if probe is None:
            return max(vals, key=lambda e: (e.free_symbols == set(), str(e)))
        num = [(complex(v.subs(probe)).real, v) for v in vals]
        return max(num)[1]

    def instability(self, probe=None):
        """``Lambda = lambda - 1``; the transition is at ``Lambda = 0``."""
        return self.perron_root(probe) - 1


def uniform_weights(L, value):
    """``L x L`` table with every entry equal to ``value``."""
    return ones(L, L) * sympify(value)
