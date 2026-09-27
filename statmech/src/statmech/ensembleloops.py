"""The loop series averaged over a random chygraph (Sec. 24.3).

On a sparse random chygraph the generalised loops of the promoted factor
graph that survive as ``n -> infinity`` are the unicyclic ones: a connected
loop with a vertex of degree three has cyclomatic number two and there are
``O(1/n)`` of them.  So at leading order the loops are the simple cycles
of the incidence graph -- the loops over complexes of Chapter 17 -- and
their disjoint unions, and since the term of a disjoint union is the
product of the terms,

    ln Z - ln Z_BP  ->  sum_{cycles C} ln(1 + r_C).

At the trivial fixed point every atom on a cycle has ``mu_i = 1`` (degree
two, zero magnetisation) and every complex contributes
``mu_a = <sigma_v sigma_u>_a``, the correlation of its two members on the
cycle under the clique belief at zero field, which is ``u'(c_a)`` of
Eq. (8.uprime).  So ``r_C = prod_{a in C} u'(c_a)``.

The number of cycles through ``l`` complexes with a given sequence of
layers is asymptotically Poisson, and the weighted mean is
``tr(B^l) / 2l`` with ``B`` the branching matrix of Eq. (8.branch): its
entries count the expected continuations from a complex of one layer to a
complex of the next, size-biased and with the excess bracket, each
carrying ``u'``.  Expanding the logarithm,

    E[ln Z - ln Z_BP]  ->  sum_{l >= 2} (1/2l) sum_{j >= 1} (-1)^{j+1}/j tr(B_j^l),

``B_j`` the branching matrix with ``u'`` replaced by ``u'^j``; summing over
``l`` first,

    E[ln Z - ln Z_BP]  ->  -1/2 sum_j (-1)^{j+1}/j [ ln det(I - B_j) + tr B_j ].

The ``j = 1`` term is ``-1/2 ln det(I - B)`` less the one-complex term: the
Gaussian correction around the Bethe saddle, which diverges at Eq. (8.det),
the ferromagnetic threshold.  The ``j = 2`` term carries the
de Almeida--Thouless matrix.  The instance-level counterpart of the
``j = 1`` term is ``-1/2 ln det(I - T)``, the Bethe Hessian of
:mod:`statmech.bethehessian`; the series sharpens ``-ln(1 - r)`` per cycle
to ``ln(1 + r)``, which is exact for a ring.
"""

import math

import networkx as nx
import numpy as np

from statmech.ising import branching_matrix, clique_derivative
from statmech.loopseries import BinaryFactorGraph


def branching_power(cardinalities, means, beta_J, j, excess=None):
    """``B_j``: the branching matrix with ``u'`` raised to the power ``j``."""
    B = branching_matrix(cardinalities, means, beta_J, excess)
    u = np.array([clique_derivative(int(c), beta_J) for c in cardinalities])
    return B * (u ** (j - 1))[None, :]


def ensemble_series(cardinalities, means, beta_J, excess=None, lmax=None,
                    jmax=12):
    """``E[ln Z - ln Z_BP]`` on the ensemble: the closed form when ``lmax``
    is ``None``, else the sum truncated at cycles of ``lmax`` complexes.
    Returns ``(total, first_order)`` with the ``j = 1`` term separately."""
    total, first = 0.0, 0.0
    for j in range(1, jmax + 1):
        Bj = branching_power(cardinalities, means, beta_J, j, excess)
        if lmax is None:
            sign, ld = np.linalg.slogdet(np.eye(len(Bj)) - Bj)
            if sign <= 0:
                return math.inf, math.inf
            s = -0.5 * (ld + np.trace(Bj))
        else:
            P = Bj.copy()
            s = 0.0
            for l in range(2, lmax + 1):
                P = P @ Bj
                s += np.trace(P) / (2 * l)
        term = (-1) ** (j + 1) / j * s
        total += term
        if j == 1:
            first = term
    return float(total), float(first)


def incidence_graph(complexes):
    G = nx.Graph()
    for a, c in enumerate(complexes):
        for v in c:
            G.add_edge(('v', int(v)), ('a', a))
    return G


def cycle_sum(complexes, beta_J, lmax=8):
    """``sum_C ln(1 + r_C)`` over the simple cycles of the incidence graph
    through at most ``lmax`` complexes, with ``r_C = prod u'(c_a)``; also
    the number of cycles by length."""
    G = incidence_graph(complexes)
    u = {a: clique_derivative(len(c), beta_J) for a, c in enumerate(complexes)}
    total, counts = 0.0, {}
    for cyc in nx.simple_cycles(G, length_bound=2 * lmax):
        cs = [x[1] for x in cyc if x[0] == 'a']
        r = float(np.prod([u[a] for a in cs]))
        total += math.log1p(r)
        counts[len(cs)] = counts.get(len(cs), 0) + 1
    return total, counts


def exact_minus_bethe(fg):
    """``ln Z - ln Z_BP`` at the trivial fixed point, ``Z`` by enumeration
    (small instances only)."""
    for k in fg.m_va:
        fg.m_va[k] = np.zeros(2)
    for k in fg.m_av:
        fg.m_av[k] = np.zeros(2)
    return fg.exact_log_Z() - fg.log_Z_bethe()


# ---------------------------------------------------------------------------
# the hard-core gas of cycles, and the O(1/n) term
# ---------------------------------------------------------------------------
#
# At the trivial fixed point every atom on a cycle has mu_i = 1 when its
# degree in the loop is even and mu_i = 0 when it is odd, and a complex
# carries the correlation of the members the loop uses, which vanishes for
# an odd number of them.  On a family of edges and triangles no complex can
# hold two cycles without an odd vertex, so the generalised loops that
# survive are the collections of cycles pairwise sharing no complex --
# sharing atoms is free -- and the series is the partition function of a
# hard-core gas of cycles,
#
#     Z / Z_BP = sum_{S: pairwise complex-disjoint} prod_{C in S} r_C,
#
# exactly.  Its Mayer expansion gives ln(Z / Z_BP) = sum_C ln(1 + r_C)
# - sum_{incompatible pairs} r_C r_C' + ..., the pair sum being O(1/n) on
# a sparse random chygraph.  Cliques of four or more members admit two
# cycles through disjoint member pairs, with the four-point correlation
# in place of the product; that case is not covered here.

def cycle_gas(complexes, beta_J, lmax=10):
    """The cycles through at most ``lmax`` complexes, as ``(complexes, r)``."""
    G = incidence_graph(complexes)
    u = {a: clique_derivative(len(c), beta_J) for a, c in enumerate(complexes)}
    out = []
    for cyc in nx.simple_cycles(G, length_bound=2 * lmax):
        cs = frozenset(x[1] for x in cyc if x[0] == 'a')
        out.append((cs, float(np.prod([u[a] for a in cs]))))
    return out


def incompatible_pairs(cycles):
    by = {}
    for i, (cs, _) in enumerate(cycles):
        for a in cs:
            by.setdefault(a, []).append(i)
    pairs = set()
    for lst in by.values():
        for i in range(len(lst)):
            for j in range(i + 1, len(lst)):
                pairs.add((lst[i], lst[j]))
    return pairs


def hardcore_log_sum(cycles):
    """``ln sum_S prod r`` over pairwise complex-disjoint collections, by
    the independence-polynomial recursion (fine while conflicts are sparse)."""
    pairs = incompatible_pairs(cycles)
    adj = {i: set() for i in range(len(cycles))}
    for i, j in pairs:
        adj[i].add(j)
        adj[j].add(i)

    def Z(alive):
        if not alive:
            return 1.0
        v = max(alive, key=lambda i: len(adj[i] & alive))
        rest = alive - {v}
        return Z(rest) + cycles[v][1] * Z(rest - adj[v])

    return math.log(Z(frozenset(range(len(cycles)))))


def pair_correction(cycles):
    """``sum_{incompatible pairs} r_C r_C'``, the first Mayer term."""
    return sum(cycles[i][1] * cycles[j][1] for i, j in incompatible_pairs(cycles))


def pair_correction_ensemble(cardinalities, means, n, beta_J, lmax=10):
    """The annealed estimate: cycles through a given layer-``m`` complex
    number ``sum_l (B^l)_{mm} / 2 M_m`` on average, ``M_m = n <kappa>_m / c_m``
    complexes in the layer, so unordered pairs at a complex sum to
    ``sum_m [sum_l (B^l)_{mm}]^2 / 8 M_m``.  Pairs sharing a path are
    counted once per shared complex and correlated pairs are not, so the
    estimate is a leading-order one."""
    B = branching_power(cardinalities, means, beta_J, 1)
    L = len(cardinalities)
    M = [n * means[m] / cardinalities[m] for m in range(L)]
    S = np.zeros((L, L))
    P = np.linalg.matrix_power(B, 2)
    for _ in range(2, lmax + 1):
        S += P
        P = P @ B
    return float(sum(S[m, m] ** 2 / (8 * M[m]) for m in range(L)))
