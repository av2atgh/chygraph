"""The loop series of Chertkov and Chernyak on a binary factor graph.

For a model of binary variables with factors ``f_a(sigma_a)``, the partition
function expands around any belief-propagation fixed point as

    Z = Z_BP * (1 + sum_C r_C),

the sum running over the *generalised loops* ``C`` of the factor graph --
subgraphs in which every vertex present has degree at least two -- with

    r_C = prod_{i in C} mu_i * prod_{a in C} mu_a,
    mu_i = [ (1 - m_i)^{q_i - 1} + (-1)^{q_i} (1 + m_i)^{q_i - 1} ]
           / [ 2 (1 - m_i^2)^{q_i - 1} ],
    mu_a = sum_{sigma_a} b_a(sigma_a) prod_{i in a cap C} (sigma_i - m_i),

``m_i`` the belief-propagation magnetisation of variable ``i``, ``b_a`` the
factor belief, and ``q_i`` the degree of ``i`` in ``C`` (Chertkov and
Chernyak, J. Stat. Mech. P06009, 2006; Phys. Rev. E 73, 065102, 2006).  The
identity is exact on any finite factor graph; the point of it here is which
factor graph.

Two factor graphs represent the same Ising model on a clustered graph.  The
*pairwise* one has a factor per edge; its generalised loops include every
triangle, every union of triangles, and every longer cycle.  The *promoted*
one is the factor graph of alpha + I of Chapter 19 with the maximal cliques
as factors, every bond assigned to exactly one of the cliques that contain
it, so that the model is the same; its generalised loops are the loops over
complexes.  Both expansions are exact; they run around different fixed
points, and Sec. 24.2 of the book asks what one buys over the other.  A
third factor graph, ``assign='all'``, gives every clique every bond inside
it -- the chygraph recursion's own convention, Chapter 14's double count --
and represents a different model; its loop series expands *that* model's
partition function, which is what separates the double-count part of
Chapter 14's error from the loop part.
"""

from itertools import combinations

import numpy as np
from scipy.special import logsumexp

from statmech.gbp import exact_log_Z, ising_factors, _lift

SPIN = np.array([1.0, -1.0])   # state 0 is +1, state 1 is -1, as in gbp


def _norm(v):
    return v - logsumexp(v)


class BinaryFactorGraph:
    """A factor graph over binary variables, its belief propagation, and its
    exact partition function.

    Args:
        factors: iterable of ``(scope, log_table)`` with ``scope`` a tuple of
            variable labels and ``log_table`` an array of shape ``(2,) * len``.
    """

    def __init__(self, factors):
        self.factors = [(tuple(int(v) for v in sc), np.asarray(t, float))
                        for sc, t in factors]
        self.nodes = sorted({v for sc, _ in self.factors for v in sc})
        self.of = {v: [a for a, (sc, _) in enumerate(self.factors) if v in sc]
                   for v in self.nodes}
        self.m_va = {(v, a): np.zeros(2) for a, (sc, _) in enumerate(self.factors)
                     for v in sc}
        self.m_av = {(a, v): np.zeros(2) for a, (sc, _) in enumerate(self.factors)
                     for v in sc}
        self.residual = np.inf
        self.sweeps = 0

    # -- constructors -------------------------------------------------------

    @classmethod
    def pairwise(cls, edges, beta_J, field=0.0):
        """One factor per edge: the ordinary factor graph of the graph."""
        return cls(ising_factors(edges, beta_J, field))

    @classmethod
    def promoted(cls, complexes, edges, beta_J, assign='once', field=0.0):
        """One factor per complex, the factor graph of alpha + I.

        ``assign='once'`` gives each bond to the largest complex containing
        it (ties to the first), so that the model is the pairwise model;
        ``'all'`` gives every complex every bond inside it, the chygraph
        recursion's double count.
        """
        cx = [tuple(sorted(int(v) for v in c)) for c in complexes]
        bonds = {tuple(sorted(map(int, e))) for e in edges}
        owner = {}
        if assign == 'once':
            order = sorted(range(len(cx)), key=lambda a: (-len(cx[a]), a))
            for a in order:
                for p, q in combinations(cx[a], 2):
                    if (p, q) in bonds and (p, q) not in owner:
                        owner[(p, q)] = a
            missing = bonds - set(owner)
            if missing:
                raise ValueError(f'bonds inside no complex: {sorted(missing)}')
        factors = []
        for a, c in enumerate(cx):
            n = len(c)
            # axis p of the C-ordered table is bit n-1-p of the flat index,
            # so read the spin of position p from that bit; getting this
            # backwards reverses the scope and misplaces every bond of a
            # complex whose bond set is not reversal-symmetric
            s = np.array([[1.0 if (i >> (n - 1 - b)) & 1 == 0 else -1.0
                           for b in range(n)] for i in range(2 ** n)])
            e = np.zeros(2 ** n)
            for p, q in combinations(range(n), 2):
                bond = (c[p], c[q])
                if bond not in bonds:
                    continue
                if assign == 'all' or owner.get(bond) == a:
                    e += s[:, p] * s[:, q]
            table = (beta_J * e).reshape((2,) * n)
            if field:
                for p in range(n):
                    sh = [1] * n
                    sh[p] = 2
                    table = table + (field * SPIN).reshape(sh) / max(1, len(
                        [b for b in cx if c[p] in b]))
            factors.append((c, table))
        return cls(factors)

    # -- belief propagation -------------------------------------------------

    def _factor_to_var(self, a, v):
        sc, t = self.factors[a]
        acc = t.copy()
        for i, u in enumerate(sc):
            if u == v:
                continue
            sh = [1] * len(sc)
            sh[i] = 2
            acc = acc + self.m_va[(u, a)].reshape(sh)
        ax = tuple(i for i, u in enumerate(sc) if u != v)
        return _norm(logsumexp(acc, axis=ax) if ax else acc)

    def sweep(self, damping=0.5):
        d, delta = damping, 0.0
        for a, (sc, _) in enumerate(self.factors):
            for v in sc:
                new = self._factor_to_var(a, v)
                delta = max(delta, float(np.abs(new - self.m_av[(a, v)]).max()))
                self.m_av[(a, v)] = d * self.m_av[(a, v)] + (1 - d) * new
        for v in self.nodes:
            for a in self.of[v]:
                acc = np.zeros(2)
                for b in self.of[v]:
                    if b != a:
                        acc = acc + self.m_av[(b, v)]
                new = _norm(acc)
                delta = max(delta, float(np.abs(new - self.m_va[(v, a)]).max()))
                self.m_va[(v, a)] = d * self.m_va[(v, a)] + (1 - d) * new
        return delta

    def bp(self, damping=0.5, sweeps=20000, tol=1e-13):
        """Run belief propagation to a fixed point; returns ``self``."""
        for k in range(sweeps):
            self.residual = self.sweep(damping)
            self.sweeps = k + 1
            if self.residual < tol:
                break
        return self

    def converged(self, tol=1e-9):
        return self.residual < tol

    # -- beliefs and free energies ------------------------------------------

    def node_belief(self, v):
        acc = np.zeros(2)
        for a in self.of[v]:
            acc = acc + self.m_av[(a, v)]
        return np.exp(_norm(acc))

    def magnetisation(self, v):
        return float(self.node_belief(v) @ SPIN)

    def factor_belief(self, a):
        sc, t = self.factors[a]
        acc = t.copy()
        for i, u in enumerate(sc):
            sh = [1] * len(sc)
            sh[i] = 2
            acc = acc + self.m_va[(u, a)].reshape(sh)
        return np.exp(acc - logsumexp(acc))

    def log_Z_bethe(self):
        """``ln Z`` at the fixed point: factors, nodes, minus the inclusions."""
        tot = 0.0
        for a, (sc, t) in enumerate(self.factors):
            acc = t.copy()
            for i, u in enumerate(sc):
                sh = [1] * len(sc)
                sh[i] = 2
                acc = acc + self.m_va[(u, a)].reshape(sh)
            tot += float(logsumexp(acc))
        for v in self.nodes:
            acc = np.zeros(2)
            for a in self.of[v]:
                acc = acc + self.m_av[(a, v)]
            tot += float(logsumexp(acc))
        for (a, v), m in self.m_av.items():
            tot -= float(logsumexp(m + self.m_va[(v, a)]))
        return tot

    def exact_log_Z(self):
        return exact_log_Z(self.factors, self.nodes)


# ---------------------------------------------------------------------------
# generalised loops
# ---------------------------------------------------------------------------

def generalised_loops(fg, max_edges=None):
    """Every generalised loop of ``fg``: a choice, for each factor, of none
    or at least two of its legs, such that every variable touched is touched
    at least twice.  Yields tuples ``((a, legs), ...)`` over the factors
    present.  ``max_edges`` bounds the number of legs in a loop.

    Enumeration is by backtracking over the factors, pruning a variable's
    degree as soon as its last factor has been decided.
    """
    F = len(fg.factors)
    last = {}
    for a, (sc, _) in enumerate(fg.factors):
        for v in sc:
            last[v] = a
    finish = {}
    for v, a in last.items():
        finish.setdefault(a, []).append(v)
    options = []
    for sc, _ in fg.factors:
        opts = [()]
        for k in range(2, len(sc) + 1):
            opts.extend(combinations(sc, k))
        options.append(opts)
    deg = {v: 0 for v in fg.nodes}
    chosen = []

    def rec(a, edges):
        if a == F:
            if edges:
                yield tuple(chosen)
            return
        for legs in options[a]:
            if max_edges is not None and edges + len(legs) > max_edges:
                continue
            for v in legs:
                deg[v] += 1
            ok = all(deg[v] != 1 for v in finish.get(a, ()))
            if ok:
                if legs:
                    chosen.append((a, legs))
                yield from rec(a + 1, edges + len(legs))
                if legs:
                    chosen.pop()
            for v in legs:
                deg[v] -= 1

    yield from rec(0, 0)


def loop_term(fg, loop):
    """``r_C`` for one generalised loop, from the current fixed point."""
    q = {}
    for a, legs in loop:
        for v in legs:
            q[v] = q.get(v, 0) + 1
    r = 1.0
    for v, qv in q.items():
        m = fg.magnetisation(v)
        r *= ((1 - m) ** (qv - 1) + (-1) ** qv * (1 + m) ** (qv - 1)) \
            / (2 * (1 - m * m) ** (qv - 1))
    for a, legs in loop:
        sc, _ = fg.factors[a]
        b = fg.factor_belief(a)
        n = len(sc)
        acc = b
        for i, u in enumerate(sc):
            if u in legs:
                sh = [1] * n
                sh[i] = 2
                acc = acc * (SPIN - fg.magnetisation(u)).reshape(sh)
        r *= float(acc.sum())
    return r


def loop_series(fg, max_edges=None):
    """All terms ``(loop, number_of_legs, r_C)`` of the series, in order of
    the number of legs."""
    out = []
    for loop in generalised_loops(fg, max_edges):
        size = sum(len(legs) for _, legs in loop)
        out.append((loop, size, loop_term(fg, loop)))
    out.sort(key=lambda t: (t[1], t[0]))
    return out


def partial_sums(terms):
    """``(size, cumulative sum of r_C over loops of at most that size)``."""
    sizes = sorted({s for _, s, _ in terms})
    acc, out = 0.0, []
    for s in sizes:
        acc += sum(r for _, ss, r in terms if ss == s)
        out.append((s, acc))
    return out


__all__ = ['BinaryFactorGraph', 'generalised_loops', 'loop_term',
           'loop_series', 'partial_sums']
