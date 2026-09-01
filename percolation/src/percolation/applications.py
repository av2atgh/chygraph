"""Percolation and spreading problems from the literature, mapped to chygraphs.

Every construction here is stated as an explicit chygraph mapping: which layer
holds which objects, what the chy-degree generating functions ``Phi^l`` are, and
what the intra-complex generating functions ``G^l`` are.  Once those are written
down the threshold, the giant component fraction and the critical amplitude all
follow from the generic machinery in :mod:`percolation.percolation`,
:mod:`percolation.giant`, :mod:`percolation.joint` and :mod:`percolation.amplitude`.

The three families are

1. AND-logic versus OR-logic hypergraph percolation
   (Bianconi & Dorogovtsev, Phys. Rev. E 109, 014306 (2024); Ha, Neri &
   Annibale, J. Phys. Complexity 6, 045011 (2025)).  Same chygraph, two
   different ``Gbar^l``.

2. Hypergraphs with hyperdegree-cardinality correlations
   (Ha, Neri & Annibale; Valdez & La Rocca, arXiv:2604.14065; Fujiki &
   Mizutaka, arXiv:2403.05162).  Hyperedges are split into layers by
   cardinality class, and the correlation becomes a joint chy-degree.

3. Epidemics with two levels of mixing: households plus a global contact
   network (Ball, Mollison & Scalia-Tomba, Ann. Appl. Probab. 7, 46 (1997);
   Ball, Sirl & Trapman, Math. Biosci. 224, 53 (2010)).  The within-household
   final size distribution is exactly ``Gbar^l``.  The same machinery covers
   networks of cliques and bipartite projections with clique size
   fluctuations.

4. Site percolation on Static Triadic Closure graphs, equivalently
   extended-range percolation with ``R = 2`` (Cirigliano, EPL 152, 31002
   (2025), arXiv:2506.17175).  The chygraph lives on the *backbone*, not on the
   clustered graph, and the range-2 rule is carried by an entry-dependent
   ``Gbar^{l,(m)}``.
"""

from sympy import Rational, binomial, expand, symbols, sympify

from percolation.giant import GiantComponent, poisson_pgf
from percolation.joint import JointGiantComponent


# ---------------------------------------------------------------------------
# Bond percolation inside a complete graph: the Reed-Frost final size
# ---------------------------------------------------------------------------

def _connected_prob(n, p):
    """Probability that ``G(j, p)`` is connected, for ``j = 1..n``."""
    conn = [None] * (n + 1)
    conn[1] = sympify(1)
    for j in range(2, n + 1):
        conn[j] = expand(1 - sum(binomial(j - 1, i - 1) * conn[i]
                                 * (1 - p)**(i * (j - i)) for i in range(1, j)))
    return conn


def clique_cluster_distribution(n, p):
    """Size of the bond-percolation cluster of a given vertex of ``K_n``.

    Returns ``{size: probability}``.  For a Reed-Frost epidemic in a closed
    group of ``n`` individuals with per-pair transmission probability ``p``,
    this is the final size distribution including the initial case, which is
    the classical result of Ball (Adv. Appl. Probab. 18, 289 (1986)) in the
    form obtained by identifying the epidemic with bond percolation on ``K_n``.
    """
    conn = _connected_prob(n, p)
    return {j: expand(binomial(n - 1, j - 1) * conn[j] * (1 - p)**(j * (n - j)))
            for j in range(1, n + 1)}


def clique_excess_pgf(n, p):
    """``Gbar``: PGF of the number of *other* members of ``K_n`` reached.

    For ``n = 3`` this reproduces the triangle enumeration of Vazquez,
    J. Complex Netw. 12, cnae047 (2024), Fig. 3.
    """
    d = clique_cluster_distribution(n, p)
    return lambda y: sum(pr * y**(j - 1) for j, pr in d.items())


def clique_pgf(n, p):
    """``G``: PGF of the component size of a uniformly chosen member of ``K_n``."""
    d = clique_cluster_distribution(n, p)
    return lambda y: sum(pr * y**j for j, pr in d.items())


def _mixture(weights, factory):
    """PGF of a mixture: ``sum_n w_n f_n(y)``."""
    parts = [(w, factory(n)) for n, w in weights.items()]
    return lambda y: sum(w * f(y) for w, f in parts)


def size_biased(weights):
    """Size-bias a distribution over sizes: ``n P(n) / <n>``."""
    mean = sum(n * w for n, w in weights.items())
    return {n: Rational(n, 1) * w / mean for n, w in weights.items()}


# ---------------------------------------------------------------------------
# 1. AND-logic versus OR-logic hypergraph percolation
# ---------------------------------------------------------------------------

def and_or_hypergraph(logic='or', p=None, degree=None, excess_degree=None,
                      cardinality=None, excess_cardinality=None, poisson=True):
    """Site percolation on a hypergraph under AND- or OR-logic.

    Chygraph mapping
    ----------------
    ``layer 0``  nodes, as atoms.
    ``layer 1``  hyperedges, as complexes whose intra-complex hypergraph is the
                 single hyperedge holding their member nodes.

    ``Phi^0_1(x)``  hyperdegree generating function of a node; ``Phibar^0_1``
                    its excess.  Neither is thinned: under both logics a node
                    reached along a *functioning* hyperedge is present.  Node
                    occupation enters through ``root_occupation = [p, 1]``, so
                    ``S_0`` is the fraction of all nodes.
    ``Phi^1_k``     absent; hyperedges are included in nothing.
    ``Gbar^1_0(y)`` is where the two logics differ, and is the only difference:

        OR  (factor graph percolation; the hyperedge keeps connecting whichever
             members survive)          Gbar(y) = Gbar_c(1 - p + p y)

        AND (hypergraph percolation; the hyperedge fails if any one of its
             members is removed)       Gbar(y) = 1 - Gbar_c(p) + Gbar_c(p y)

    where ``Gbar_c`` is the excess cardinality generating function.  Both reduce
    to ``Gbar_c(y)`` at ``p = 1``, so the two problems share a chygraph and
    differ in one generating function.

    Args:
        logic: ``'or'`` or ``'and'``.
        p: node occupation probability (symbol ``p`` by default).
        poisson: build the degree and cardinality generating functions from
            their means ``k`` and ``c``.
    """
    if logic not in ('or', 'and'):
        raise ValueError("logic must be 'or' or 'and'")
    k, c = symbols('k c')
    p = symbols('p') if p is None else p

    if poisson:
        degree = degree or poisson_pgf(k)
        excess_degree = excess_degree or degree
        cardinality = cardinality or poisson_pgf(c)
        excess_cardinality = excess_cardinality or cardinality
    if excess_cardinality is None:
        raise ValueError("supply excess_cardinality for non-Poisson input")
    Gc = excess_cardinality

    if logic == 'or':
        gbar = lambda y: Gc(1 - p + p * y)
    else:
        gbar = lambda y: 1 - Gc(p) + Gc(p * y)

    phi = [[None, degree], [None, None]]
    phibar = [[None, excess_degree], [None, None]]
    g = [[None, None], [lambda y: cardinality(y), None]]
    gb = [[None, None], [gbar, None]]
    return GiantComponent(phi, phibar, g, gb, root_occupation=[p, 1])


# ---------------------------------------------------------------------------
# 2. Hyperdegree-cardinality correlations
# ---------------------------------------------------------------------------

def correlated_cardinality_hypergraph(cardinalities, joint_degree_pgf, p=None):
    """Site percolation on a hypergraph correlating hyperdegree and cardinality.

    Chygraph mapping
    ----------------
    ``layer 0``       nodes, as atoms.
    ``layer 1..L``    hyperedges, one layer per cardinality class ``c_l``.

    ``Phi^0(x_1, ..., x_L)`` is the *joint* generating function of the vector
    ``(kappa_1, ..., kappa_L)`` counting how many hyperedges of each cardinality
    class contain the node.  Correlation between a node's hyperdegree and the
    cardinality of the hyperedges it joins is exactly a correlation in this
    joint distribution, and enters the threshold tensor through the
    inclusion-biased moments ``<kappabar^(m)>_{0k} = <kappa_m kappa_k>/<kappa_m>``
    of :mod:`percolation.joint`.  Marginal hyperdegree distributions alone do not
    determine it.

    ``Gbar^l_0(y) = (1 - p + p y)^{c_l - 1}``: a hyperedge of cardinality
    ``c_l`` reaches each of its other ``c_l - 1`` members when that member is
    present (OR-logic).  Node occupation is reinstated at the root.

    Args:
        cardinalities: list of cardinalities, one per hyperedge layer.
        joint_degree_pgf: callable ``f(x)`` on a length-``(L+1)`` sequence,
            using ``x[1] ... x[L]``.
        p: node occupation probability.
    """
    p = symbols('p') if p is None else p
    L = len(cardinalities) + 1
    G = [None] * L
    Gbar = [[None] * L for _ in range(L)]
    for i, cl in enumerate(cardinalities, start=1):
        G[i] = (lambda cl: lambda y: y[0] * (1 - p + p * y[0])**(cl - 1))(cl)
        Gbar[i][0] = (lambda cl: lambda y: (1 - p + p * y[0])**(cl - 1))(cl)
    return JointGiantComponent(Phi=[joint_degree_pgf] + [None] * (L - 1),
                               G=G, Gbar=Gbar,
                               root_occupation=[p] + [1] * (L - 1))


def two_class_joint_degree(mean_a, mean_b, spread, mixing):
    """A one-parameter family of joint hyperdegrees with *fixed* marginals.

    A node is one of two equally likely types.  With probability ``w`` the types
    pair the high mean in class 1 with the high mean in class 2, and with
    probability ``1 - w`` they pair high with low:

        Phi(x) = w   [ (1/2) e^{a+(x1-1)+b+(x2-1)} + (1/2) e^{a-(x1-1)+b-(x2-1)} ]
             + (1-w) [ (1/2) e^{a+(x1-1)+b-(x2-1)} + (1/2) e^{a-(x1-1)+b+(x2-1)} ]

    with ``a± = mean_a (1 ± spread)`` and ``b± = mean_b (1 ± spread)``.  Both
    marginals are the same mixture of Poissons for every ``w``, so the published
    threshold tensor, which sees only marginals, cannot distinguish the family;
    ``Cov(kappa_1, kappa_2) = (2w - 1) spread^2 mean_a mean_b`` runs from
    maximally negative at ``w = 0`` through independent at ``w = 1/2`` to
    maximally positive at ``w = 1``.
    """
    ap, am = mean_a * (1 + spread), mean_a * (1 - spread)
    bp, bm = mean_b * (1 + spread), mean_b * (1 - spread)
    w = mixing
    from sympy import exp

    def pgf(x):
        e = lambda u, v: exp(u * (x[1] - 1) + v * (x[2] - 1))
        return (w * (e(ap, bp) + e(am, bm)) / 2
                + (1 - w) * (e(ap, bm) + e(am, bp)) / 2)
    return pgf


# ---------------------------------------------------------------------------
# 3. Epidemics with two levels of mixing
# ---------------------------------------------------------------------------

def household_epidemic(household_sizes, p_H=None, T=None,
                       degree=None, excess_degree=None, poisson=True):
    """SIR with within-household and global transmission, as a chygraph.

    Chygraph mapping
    ----------------
    ``layer 0``  individuals, as atoms.
    ``layer 1``  households, as complexes whose intra-complex graph is the
                 complete graph on their members.
    ``layer 2``  global contacts, as complexes holding two individuals.

    ``Phi^0_1(x) = x``      every individual belongs to exactly one household,
                            so ``kappa_1 = 1`` and ``Phibar^0_1 = 1``: an
                            individual reached from their household has no
                            second household to pass through.  This is what
                            makes the household layer a sink rather than a route.
    ``Phi^0_2(x)``          global degree generating function; ``Phibar^0_2`` its
                            excess.
    ``Gbar^1_0(y)``         the within-household final size distribution, minus
                            the index case.  A Reed-Frost epidemic in a closed
                            group is bond percolation on ``K_n``, so this is
                            ``clique_excess_pgf(n, p_H)``, mixed over the
                            *size-biased* household size distribution because a
                            household is reached through one of its members.
    ``Gbar^2_0(y) = 1 - T + T y``   a global contact transmits with probability
                            ``T``.

    The threshold then reproduces the household reproduction number of
    Ball, Mollison & Scalia-Tomba,

        theta + 1 = R* = T [ <kbar> + mu_H <k> ],

    with ``mu_H = Gbar^{1'}_0(1)`` the mean number of additional household
    members infected, and ``S_0`` is the final size of a major outbreak.

    Setting ``T = 0`` (or omitting layer 2) leaves a network of cliques, which
    is also the bipartite projection studied by Fujiki & Mizutaka.

    Args:
        household_sizes: ``{n: probability}`` over households (not individuals).
        p_H: within-household per-pair transmission probability.
        T: per-contact global transmission probability.
        degree: global degree PGF; built from the mean ``k`` if ``poisson``.
    """
    p_H = symbols('p_H') if p_H is None else p_H
    T = symbols('T') if T is None else T
    k = symbols('k')
    if poisson:
        degree = degree or poisson_pgf(k)
        excess_degree = excess_degree or degree
    if degree is None or excess_degree is None:
        raise ValueError("supply degree and excess_degree for non-Poisson input")

    biased = size_biased(household_sizes)
    gbar_H = _mixture(biased, lambda n: clique_excess_pgf(n, p_H))
    g_H = _mixture(biased, lambda n: clique_pgf(n, p_H))

    phi = [[None, lambda x: x, degree], [None] * 3, [None] * 3]
    phibar = [[None, lambda x: sympify(1), excess_degree], [None] * 3, [None] * 3]
    g = [[None] * 3,
         [g_H, None, None],
         [lambda y: (1 - T) * y + T * y**2, None, None]]
    gb = [[None] * 3,
          [gbar_H, None, None],
          [lambda y: 1 - T + T * y, None, None]]
    return GiantComponent(phi, phibar, g, gb)


def stc_percolation(backbone_degree=None, phi=None, poisson=True):
    """Site percolation on a Static Triadic Closure (STC) random graph.

    An STC graph ``G1`` is built from a treelike backbone ``G0`` by closing every
    triad, so each node's closed neighbourhood becomes a clique.  ``G1`` has
    scale-free degrees, overlapping short loops and assortative correlations,
    and it is *not* locally treelike.  Site percolation on ``G1`` is nonetheless
    exactly equivalent to extended-range percolation with range ``R = 2`` on the
    backbone: two occupied nodes are connected when at most one unoccupied node
    lies between them.

    Chygraph mapping
    ----------------
    Not on ``G1`` -- on the backbone, with the occupation state promoted to a
    layer:

    ``layer 0``  occupied backbone nodes, as atoms.
    ``layer 1``  unoccupied backbone nodes, as atoms.
    ``layer 2``  backbone edges, as complexes of cardinality two.

    ``Phi^0(x) = Phi^1(x) = g_0(x_2)``   a node sits in ``k`` backbone edges;
                            ``Phibar^{0,(2)} = g_1`` is derived, since
                            ``g_1 = g_0' / g_0'(1)`` by definition.
    ``Gbar^{2,(0)}(y) = phi y_0 + (1 - phi) y_1``   entered from an occupied
                            node: the far end may be either.
    ``Gbar^{2,(1)}(y) = phi y_0 + (1 - phi)``       entered from an unoccupied
                            node: a second unoccupied node in a row kills the
                            path, which is the ``R = 2`` constraint.

    The whole range-2 rule is carried by letting ``Gbar`` depend on the layer the
    complex was entered from -- the same freedom needed for motifs in
    :class:`percolation.joint.JointGiantComponent`.  Occupation is reinstated at
    the root, so ``S_0`` is the fraction of *all* backbone nodes.

    For a Poisson backbone this gives

        theta = -b^2 phi^2 + b(1 + b) phi - 1,
        phi_c = [(1 + b) - sqrt((1 + b)^2 - 4)] / (2b),     b = <k(k-1)>/<k>,

    which is Eq. (12) of Cirigliano, and the order parameter reproduces his
    Eqs. (9)-(11) to machine precision.

    Why the obvious mapping does not work.  Taking ``G1`` itself as a hypergraph
    -- nodes as atoms, closed neighbourhoods as complexes -- fails, because every
    backbone edge ``u ~ v`` creates the four-cycle ``u - C(v) - v - C(u) - u`` in
    the node/complex incidence graph.  Chygraphs are exact when that incidence
    structure is locally treelike, and there it is not: for a Poisson backbone
    with mean degree 3 the naive mapping gives ``phi_c = 0.0711`` against the
    exact ``0.0893``, no better than heterogeneous mean field's ``0.0635``.

    Args:
        backbone_degree: PGF ``g_0`` of the backbone degree distribution.
            Built from the mean ``k`` when ``poisson`` is set.
        phi: occupation probability (symbol ``phi`` by default).
    """
    from sympy import exp
    phi = symbols('phi') if phi is None else phi
    if poisson:
        k = symbols('k')
        backbone_degree = backbone_degree or (lambda z: exp(k * (z - 1)))
    if backbone_degree is None:
        raise ValueError("supply backbone_degree for a non-Poisson backbone")
    g0 = backbone_degree

    Phi = [lambda x: g0(x[2]), lambda x: g0(x[2]), None]
    G = [None, None, lambda y: (phi * y[0] + (1 - phi) * y[1])**2]
    Gbar = [[None] * 3 for _ in range(3)]
    Gbar[2][0] = lambda y: phi * y[0] + (1 - phi) * y[1]
    Gbar[2][1] = lambda y: phi * y[0] + (1 - phi)
    return JointGiantComponent(Phi=Phi, G=G, Gbar=Gbar,
                               root_occupation=[phi, 1, 1])


def clique_network(clique_sizes, p_bond=None, degree=None, excess_degree=None,
                   poisson=True):
    """Bond percolation on a network of cliques.

    Chygraph mapping
    ----------------
    ``layer 0``  individuals, as atoms.
    ``layer 1``  cliques, as complexes whose intra-complex graph is complete and
                 carries bond percolation with probability ``p_bond``.

    This is :func:`household_epidemic` with the global layer dropped *and*
    ``Phi^0_1`` a genuine degree distribution rather than the constant 1: a node
    belongs to many cliques, so the cliques themselves form a connected
    structure.  Households alone never percolate precisely because
    ``Phi^0_1(x) = x`` there.

    Covers SIR on clique random networks with clique-type dependent
    transmission, and the connected components of bipartite projections
    (Fujiki & Mizutaka) when ``p_bond = 1``.  With all cliques of size two it
    reduces to ordinary bond percolation; with size three it reproduces the
    triangle layer of ``graph_with_triangles_giant``.
    """
    p_bond = symbols('q') if p_bond is None else p_bond
    k = symbols('k')
    if poisson:
        degree = degree or poisson_pgf(k)
        excess_degree = excess_degree or degree
    biased = size_biased(clique_sizes)
    phi = [[None, degree], [None, None]]
    phibar = [[None, excess_degree], [None, None]]
    g = [[None, None], [_mixture(biased, lambda n: clique_pgf(n, p_bond)), None]]
    gb = [[None, None],
          [_mixture(biased, lambda n: clique_excess_pgf(n, p_bond)), None]]
    return GiantComponent(phi, phibar, g, gb)


__all__ = [
    "clique_cluster_distribution", "clique_excess_pgf", "clique_pgf",
    "size_biased",
    "and_or_hypergraph",
    "correlated_cardinality_hypergraph", "two_class_joint_degree",
    "household_epidemic", "clique_network",
    "stc_percolation",
]
