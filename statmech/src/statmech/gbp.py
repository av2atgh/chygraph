"""Generalised belief propagation on the region graph (WP5, the standing item).

:mod:`chygraph_statmech.region` builds the region graph and assigns Mobius
counting numbers, and it stops there: it *measures* how far a chygraph is from
treelike without repairing anything.  The repair is generalised belief
propagation [Yedidia, Freeman & Weiss, IEEE Trans. Inf. Theory 51, 2282
(2005)], and this module is it.

Why it is needed.  Everything before Sec. VII treats a node's complexes as
independent, which is exact when complexes meet in at most one node.  Two
cliques sharing an edge are a loop the chygraph cannot see.  Closing the
complex family under intersection and counting by Mobius inversion says *where*
the Bethe counting is wrong; it does not fix the messages, and evaluating the
counting on isolated regions --- what Table IV of the manuscript does --- is a
static estimate that never lets the regions talk to each other.

The parent-to-child algorithm.  Messages live on directed edges of the region
graph, from a region ``P`` to a direct child ``R``, and are functions of the
child's variables.  With ``U_R = {R} + descendants(R)``, the belief on a region
is

    b_R(x_R) ~ f_R(x_R) prod_{(I->J): J in U_R, I not in U_R} m_{I->J}(x_J)

with ``f_R`` the product of every factor whose scope lies inside ``R``.
Imposing ``sum_{x_P \\ x_R} b_P = b_R``, the messages common to both sides
cancel and what is left is the update

    m_{P->R}(x_R) = [ sum_{x_P \\ x_R} f_{P\\R}(x_P) prod_{N(P,R)} m ]
                    / prod_{D(P,R) \\ (P,R)} m,

    N(P,R) = {(I,J): J in U_P \\ U_R, I not in U_P},
    D(P,R) = {(I,J): J in U_R,       I in U_P \\ U_R}.

On the two-layer region graph of a factor model --- one region per factor plus
one per variable --- ``N(P,R)`` is the incoming messages to the factor's other
variables, ``D(P,R)`` is empty, and the update is ordinary belief propagation.
So GBP contains BP the way the Mobius counting contains the Bethe counting, and
:class:`GBP` runs both through one code path.

**What it buys, on the manuscript's own example.**  Two triangles sharing an
edge, with all pairs coupled: the closed region family is
``{0,1,2}, {1,2,3}, {1,2}`` once the zero-counting singletons are pruned, which
is a junction tree, so GBP is *exact* there.  The Bethe counting of Eq. (13)
errs by up to ``1.3`` in ``ln Z`` and the static Mobius counting by up to
``7 x 10^-2``; GBP recovers ``ln Z`` to machine precision.  That is the whole of
what Table IV said the present treatment leaves on the table, recovered.

**What it does not do.**  This is an *instance-level* calculation: it takes an
explicit list of complexes and an explicit interaction and returns numbers for
that structure.  The rest of this package works at the level of an *ensemble*,
where a message is a distribution over a chy-degree class rather than a function
on one region.  Lifting the parent-to-child update to ensembles --- messages
indexed by region *type*, with the intersection structure itself random --- is
what would close the overlap deficit on hyperbolic random graphs, and it is not
attempted here.
"""

from itertools import combinations

import numpy as np
from scipy.special import logsumexp


# ---------------------------------------------------------------------------
# lifting and projecting log-arrays between regions
# ---------------------------------------------------------------------------

def _lift(arr, sub, sup):
    """Reshape a log-array on variables ``sub`` for broadcasting onto ``sup``.

    Both tuples are sorted, so the axes of ``sub`` occur in ``sup`` in the same
    order and a reshape with singleton axes is all that is required; NumPy
    broadcasting does the rest when the result is added.
    """
    if sub == sup:
        return arr
    shape = [2 if v in sub else 1 for v in sup]
    return arr.reshape(shape)


def _project(arr, sup, sub):
    """``logsumexp`` a log-array on ``sup`` down to its sub-tuple ``sub``."""
    if sub == sup:
        return arr
    axes = tuple(i for i, v in enumerate(sup) if v not in sub)
    return logsumexp(arr, axis=axes)


LOG_CLIP = 100.0


def _normalise(arr, clip=LOG_CLIP):
    """Shift a log-array so that it sums to one, then clip it to a finite range.

    The parent-to-child update *divides* by the messages of ``D(P,R)``, so on a
    region graph where the iteration does not converge the log-messages can run
    away to ``+-inf`` and then to ``nan`` within a few sweeps.  Clipping bounds
    them without touching any run that is converging: a message within
    ``e^-100`` of a delta function is already numerically a delta function, and
    a run that reaches the clip is reported as unconverged by its residual
    rather than as ``nan``.
    """
    out = arr - logsumexp(arr)
    return np.clip(out, -clip, clip)


# ---------------------------------------------------------------------------
# factors
# ---------------------------------------------------------------------------

def ising_factors(edges, beta_J, field=0.0):
    """Log-factors for ``-beta H = beta J sum_<ij> s_i s_j + B sum_i s_i``.

    ``edges`` is an iterable of pairs, each counted **once** however many
    complexes contain it: the interaction is a property of the graph, not of
    the cover by complexes.  State ``0`` is spin ``+1`` and state ``1`` is spin
    ``-1``.
    """
    s = np.array([1.0, -1.0])
    pair = beta_J * np.outer(s, s)
    out = [((int(u), int(v)), pair) for u, v in
           sorted({tuple(sorted(map(int, e))) for e in edges})]
    if field:
        nodes = sorted({v for e in out for v in e[0]})
        out += [((v,), field * s) for v in nodes]
    return out


def clique_edges(complexes):
    """Every pair inside every complex, deduplicated."""
    return sorted({tuple(sorted((int(u), int(v))))
                   for a in complexes for u, v in combinations(sorted(a), 2)})


# ---------------------------------------------------------------------------
# the algorithm
# ---------------------------------------------------------------------------

class GBP:
    """Parent-to-child generalised belief propagation on a region graph.

    Args:
        region_graph: a :class:`~chygraph_statmech.region.RegionGraph`, or any
            object exposing ``counting`` as ``{frozenset(vars): counting
            number}``.
        log_factors: sequence of ``(scope, array)``.  ``scope`` is a tuple of
            variable labels in sorted order and ``array`` has shape
            ``(2,) * len(scope)`` holding ``ln f`` --- so for the Ising model
            ``+beta J s_i s_j`` and not its negative.  Every factor's scope must
            lie inside at least one region.
        damping: fraction of the *old* log-message kept at each update.  Zero is
            the undamped map.
        prune: drop regions of counting number zero.  They contribute nothing to
            the free energy, and removing them turns the closed family of two
            edge-sharing triangles into the junction tree it really is.

    Raises:
        ValueError: if the region graph does not count every variable and every
            factor exactly once, which is the condition under which the Kikuchi
            free energy is exact on a tree and a valid approximation otherwise.
    """

    def __init__(self, region_graph, log_factors, damping=0.5, prune=True):
        counting = dict(region_graph.counting)
        if prune:
            counting = {r: c for r, c in counting.items() if c != 0}
        if not counting:
            raise ValueError("no regions with a non-zero counting number")

        self.counting = counting
        self.regions = sorted(counting, key=lambda r: (-len(r), sorted(r)))
        self.vars = {r: tuple(sorted(r)) for r in self.regions}
        self.factors = [(tuple(sorted(sc)), np.asarray(a, dtype=float))
                        for sc, a in log_factors]
        self.damping = float(damping)

        self._check_counting()
        self._build_edges()
        self._region_factor()
        self._message_sets()
        self.reset()

    # -- validity -----------------------------------------------------------

    def _check_counting(self):
        """Every variable and every factor covered with total weight one."""
        per_var = {}
        for r, c in self.counting.items():
            for v in r:
                per_var[v] = per_var.get(v, 0) + c
        bad = {v: t for v, t in per_var.items() if t != 1}
        if bad:
            raise ValueError(f"variables not counted once: {bad}")
        for sc, _ in self.factors:
            s = frozenset(sc)
            tot = sum(c for r, c in self.counting.items() if s <= r)
            if tot != 1:
                raise ValueError(
                    f"factor on {tuple(sc)} counted {tot} times, not once; "
                    "the region family does not cover the interaction")

    # -- region-graph structure ---------------------------------------------

    def _build_edges(self):
        """Direct parents and children, plus the descendant sets ``U_R``."""
        rs = self.regions
        strict = {r: {s for s in rs if r < s} for r in rs}        # ancestors
        self.parents = {}
        for r in rs:
            anc = strict[r]
            self.parents[r] = [p for p in anc
                               if not any(p < q for q in anc)]
        self.children = {r: [] for r in rs}
        for r, ps in self.parents.items():
            for p in ps:
                self.children[p].append(r)
        self.edges = [(p, r) for r in rs for p in self.parents[r]]
        # U_R = {R} + every region strictly inside R
        self.U = {r: frozenset({r} | {s for s in rs if s < r}) for r in rs}

    def _region_factor(self):
        """``ln f_R``: every factor whose scope lies inside ``R``."""
        self.logf, self.inside = {}, {}
        for r in self.regions:
            vs = self.vars[r]
            acc = np.zeros((2,) * len(vs))
            own = []
            for n, (sc, a) in enumerate(self.factors):
                if frozenset(sc) <= r:
                    acc = acc + _lift(a, sc, vs)
                    own.append(n)
            self.logf[r] = acc
            self.inside[r] = own

    def _message_sets(self):
        """``N(P,R)`` and ``D(P,R)`` minus the edge being updated."""
        self.N, self.D, self.into = {}, {}, {}
        for r in self.regions:
            self.into[r] = [(i, j) for (i, j) in self.edges
                            if j in self.U[r] and i not in self.U[r]]
        self._only = {}
        for (p, r) in self.edges:
            self._only[(p, r)] = [n for n in self.inside[p]
                                  if n not in set(self.inside[r])]
            Up, Ur = self.U[p], self.U[r]
            outer = Up - Ur
            self.N[(p, r)] = [(i, j) for (i, j) in self.edges
                              if j in outer and i not in Up]
            self.D[(p, r)] = [(i, j) for (i, j) in self.edges
                              if j in Ur and i in outer and (i, j) != (p, r)]

    # -- iteration ----------------------------------------------------------

    def reset(self):
        """Uniform messages."""
        self.m = {(p, r): np.zeros((2,) * len(self.vars[r]))
                  for (p, r) in self.edges}
        return self

    def _update(self, p, r):
        vp, vr = self.vars[p], self.vars[r]
        # factors inside P but not inside R
        acc = np.zeros((2,) * len(vp))
        for n in self._only[(p, r)]:
            sc, a = self.factors[n]
            acc = acc + _lift(a, sc, vp)
        for (i, j) in self.N[(p, r)]:
            acc = acc + _lift(self.m[(i, j)], self.vars[j], vp)
        out = _project(acc, vp, vr)
        for (i, j) in self.D[(p, r)]:
            out = out - _lift(self.m[(i, j)], self.vars[j], vr)
        return _normalise(out)

    def sweep(self):
        """One pass over every edge, parents before children.

        Returns the largest absolute change in any log-message, which is the
        convergence diagnostic.
        """
        lam, delta = self.damping, 0.0
        for (p, r) in self.edges:
            new = self._update(p, r)
            old = self.m[(p, r)]
            mix = _normalise(lam * old + (1.0 - lam) * new)
            if not np.all(np.isfinite(mix)):
                self.m[(p, r)] = np.zeros_like(old)
                return np.inf
            delta = max(delta, float(np.max(np.abs(mix - old))))
            self.m[(p, r)] = mix
        return delta

    def run(self, sweeps=2000, tol=1e-12):
        """Iterate to convergence.  Sets ``self.sweeps`` and ``self.residual``."""
        self.residual, self.sweeps = np.inf, 0
        for n in range(1, sweeps + 1):
            self.residual, self.sweeps = self.sweep(), n
            if self.residual < tol:
                break
        return self

    def converged(self, tol=1e-9):
        return getattr(self, 'residual', np.inf) < tol

    # -- observables --------------------------------------------------------

    def belief(self, region):
        """Normalised belief on one region, as probabilities."""
        r = frozenset(region)
        vs = self.vars[r]
        acc = self.logf[r].copy()
        for (i, j) in self.into[r]:
            acc = acc + _lift(self.m[(i, j)], self.vars[j], vs)
        return np.exp(_normalise(acc))

    def log_Z(self):
        """``ln Z`` from the region-based free energy, ``-F = sum_R c_R (U+H)``.

        With ``b_R`` the region beliefs and ``ln f_R`` the sum of every factor
        inside ``R``,

            ln Z ~ sum_R c_R sum_{x_R} b_R(x_R) [ ln f_R(x_R) - ln b_R(x_R) ].

        Exact when the region graph is a junction tree, which is what makes the
        two-triangle example of Table IV come out to machine precision.
        """
        tot = 0.0
        for r in self.regions:
            b = self.belief(r)
            lb = np.where(b > 0, np.log(np.where(b > 0, b, 1.0)), 0.0)
            tot += self.counting[r] * float(np.sum(b * (self.logf[r] - lb)))
        return tot

    def consistency(self):
        """Largest violation of ``sum_{x_P \\ x_R} b_P = b_R`` over the edges.

        The condition the parent-to-child update is derived from.  A converged
        run that leaves this at zero is a genuine fixed point, so an error in
        :meth:`log_Z` there is the approximation\'s and not the solver\'s.
        """
        worst = 0.0
        for (p, r) in self.edges:
            down = self.belief(p)
            vs, vr = self.vars[p], self.vars[r]
            axes = tuple(i for i, v in enumerate(vs) if v not in r)
            marg = down.sum(axis=axes) if axes else down
            worst = max(worst, float(np.max(np.abs(marg - self.belief(r)))))
        return worst

    def magnetisation(self):
        """``<sigma_v>`` per variable, from the smallest region holding it."""
        out = {}
        for v in sorted({v for r in self.regions for v in r}):
            r = min((r for r in self.regions if v in r), key=len)
            vs = self.vars[r]
            b = self.belief(r)
            axes = tuple(i for i, w in enumerate(vs) if w != v)
            p = b.sum(axis=axes) if axes else b
            out[v] = float(p[0] - p[1])
        return out


# ---------------------------------------------------------------------------
# the static comparisons the manuscript quotes
# ---------------------------------------------------------------------------

def static_log_Z(counting, log_factors):
    """``sum_R c_R ln Z_R`` with each region's partition function isolated.

    This is the counting applied without messages --- Table IV's ``Kikuchi``
    and ``Bethe`` columns, depending on which counting is passed.  Regions never
    see each other, so it is an estimate rather than a fixed point.
    """
    factors = [(tuple(sorted(sc)), np.asarray(a, float)) for sc, a in log_factors]
    tot = 0.0
    for r, c in counting.items():
        if not c:
            continue
        vs = tuple(sorted(r))
        acc = np.zeros((2,) * len(vs))
        for sc, a in factors:
            if frozenset(sc) <= r:
                acc = acc + _lift(a, sc, vs)
        tot += c * float(logsumexp(acc))
    return tot


def exact_log_Z(log_factors, nodes=None):
    """``ln Z`` by enumeration over every spin configuration."""
    factors = [(tuple(sorted(sc)), np.asarray(a, float)) for sc, a in log_factors]
    vs = tuple(sorted(nodes)) if nodes is not None else \
        tuple(sorted({v for sc, _ in factors for v in sc}))
    acc = np.zeros((2,) * len(vs))
    for sc, a in factors:
        acc = acc + _lift(a, sc, vs)
    return float(logsumexp(acc))
