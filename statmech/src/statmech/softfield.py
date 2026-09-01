"""Hitting set with the O(1) cavity fields kept (soft fields).

:mod:`statmech.hittingset` takes ``mu -> inf`` with *hard* fields --
warning propagation, the ansatz of Vazquez & Weigt.  That is exact at
cardinality two and wrong above it, because the limit discards an ``O(1)``
contribution that survives it.  Two symptoms: on disjoint 3-hyperedges, where
replica symmetry is trivially valid, the hard-field cover size is ``1/2`` where
the truth is ``1/3``; and on random regular hypergraphs the answer misses
``rho = 1/K`` [Mezard & Tarzia, Phys. Rev. E 76, 041124 (2007)].

This module keeps the fields.  Writing ``x_i = 1`` for a vertex *in* the hitting
set and weighting configurations by ``exp(-mu sum_i x_i)``, belief propagation on
the chygraph's factor graph reads

    v_{a->i} = -ln[ 1 - prod_{j in a\\i} 1/(1 + e^{h_{j->a}}) ],
    h_{i->a} = -mu + sum_{b ni i, b != a} v_{b->i},

with ``rho = sigma(-mu + sum_{b ni i} v_{b->i})``.  The chy-degree step is again
the convolution of Sec. II; the complex step is again a sum over the interior,
here of the single forbidden configuration in which no member is taken.  Nothing
is scaled by ``mu``, so the ``O(1)`` parts survive and ``rho`` comes out as a
genuine ``O(1)`` number rather than a count of integer fields.

Two checks fix the implementation.  For a regular hypergraph with ``L``
complexes per node and cardinality ``K``, symmetry gives a closed form,

    h_RS = -mu/L - ((L-1)/L) ln(K-1),      rho = 1/K,

which is Eq. (11) of Mezard & Tarzia; and for a graph at ``K = 2`` it must
return the Weigt-Hartmann cover size, since there the hard and soft treatments
agree.  Both hold.

What this buys beyond Ref. [MT] is the chygraph structure: arbitrary chy-degree
distributions, mixed cardinalities and correlated layers, none of which the
regular ansatz reaches.

One caveat on the validity criterion.  A negative :meth:`entropy` proves the
replica-symmetric answer wrong, but a positive one does not prove it right: the
solution can be unstable while its entropy is still positive.  Vertex cover on
Erdos-Renyi graphs is the example -- replica symmetry breaks at mean degree
``e``, yet the entropy measured here is ``+0.09`` at mean degree 1 and still
positive at 3.  Ref. [MT] pairs the entropy with a separate stability criterion
for exactly this reason; only the entropy is implemented here.
"""

import numpy as np


class HittingSetBP:
    """Minimum hitting set on a chygraph, with the ``O(1)`` fields kept.

    Args:
        cardinalities: ``c_l`` per complex layer.
        degrees: chy-degree per layer -- the Poisson mean, or the exact degree
            when ``regular`` is set.
        regular: draw exactly ``degrees[l]`` complexes per node instead of
            Poisson.
        mu: chemical potential.  The minimum hitting set is ``mu -> inf``;
            ``rho`` converges and 40--80 is ample.
        size: population size.
        damping: fraction of the population replaced per sweep.  The undamped
            map flips between ``h = -mu`` and ``h = +mu(k-1)`` and never
            settles; see :meth:`_mix` for why this is a replacement and not an
            average.
    """

    def __init__(self, cardinalities, degrees, regular=False, mu=60.0,
                 size=100_000, seed=0, damping=0.5):
        self.c = np.asarray(cardinalities, dtype=int)
        self.L = len(self.c)
        self.deg = np.asarray(degrees, dtype=float)
        self.regular = regular
        self.mu = float(mu)
        self.size = size
        self.rng = np.random.default_rng(seed)
        self.damping = float(damping)
        # Start at the right scale.  The fixed point has h ~ -mu / <k>, since
        # the -mu of the chy-degree step is shared among the node's complexes;
        # starting at -mu instead makes the first sweep overshoot to +mu(k-1)
        # and the iteration then flips between the two extremes forever.
        ktot = max(float(self.deg.sum()), 1.0)
        self.P = [np.full(size, -mu / ktot) for _ in range(self.L)]
        self.Q = [np.full(size, mu / ktot) for _ in range(self.L)]

    # -- the complex step ---------------------------------------------------

    @staticmethod
    def _emit(h):
        """``-ln[1 - prod_j sigma(-h_j)]`` for ``h`` of shape ``(n, c-1)``.

        Computed as ``S = sum_j softplus(h_j)`` so the product is ``e^{-S}``;
        for the deeply negative fields the limit produces, ``S`` is tiny and
        ``v -> -ln S`` without cancellation.
        """
        S = np.logaddexp(0.0, h).sum(axis=1)
        return -np.log(-np.expm1(-S))

    # -- the chy-degree step ------------------------------------------------

    def _draw(self, n, l, excess):
        """Number of layer-``l`` complexes at a node, excess or full."""
        if self.regular:
            k = int(round(self.deg[l])) - (1 if excess else 0)
            return np.full(n, max(k, 0))
        return self.rng.poisson(self.deg[l], n)      # Poisson is its own excess

    def _sum_v(self, n, exclude=None):
        out = np.zeros(n)
        for l in range(self.L):
            d = self._draw(n, l, excess=(l == exclude))
            total = int(d.sum())
            if not total:
                continue
            draws = self.Q[l][self.rng.integers(0, self.size, total)]
            ends = np.cumsum(d)
            csum = np.concatenate(([0.0], np.cumsum(draws)))
            out += csum[ends] - csum[ends - d]
        return out

    # -- iteration ----------------------------------------------------------

    def sweep(self):
        newQ = []
        for l in range(self.L):
            idx = self.rng.integers(0, self.size,
                                    (self.size, int(self.c[l]) - 1))
            newQ.append(self._emit(self.P[l][idx]))
        self.Q = [self._mix(o, q) for o, q in zip(self.Q, newQ)]
        newP = [-self.mu + self._sum_v(self.size, exclude=l)
                for l in range(self.L)]
        self.P = [self._mix(o, q) for o, q in zip(self.P, newP)]

    def _mix(self, old, new):
        """Damp by replacing a fraction of the *population*, not by averaging.

        Averaging field values elementwise -- ``(1-lam) h_i + lam h_i^new`` --
        looks like damping but is not: the two entries are different messages,
        not two estimates of one, so averaging contracts the width of the
        distribution and moves the fixed point.  On Erdos-Renyi graphs at
        cardinality two it returns ``0.172`` where the Weigt-Hartmann cover
        size, confirmed by exact leaf removal, is ``0.272``.  Replacing a random
        subset leaves the fixed point alone and still kills the oscillation.
        """
        keep = self.rng.random(self.size) >= self.damping
        return np.where(keep, old, new)

    def run(self, sweeps=300):
        for _ in range(sweeps):
            self.sweep()
        return self

    # -- observables --------------------------------------------------------

    def density(self, samples=None):
        """``rho``, the fraction of vertices in the minimum hitting set."""
        n = samples or self.size
        H = -self.mu + self._sum_v(n)
        return float(np.mean(1.0 / (1.0 + np.exp(-H))))

    def log_Z(self, samples=None):
        """``(ln Z)/N`` from the Bethe free energy.

        With ``nu_{i->a}(x) ~ e^{h x}`` normalised to ``nu(0) = 1``, the overlap
        term collapses as in Sec. II -- ``h_{i->a} + v_{a->i}`` is the same full
        field for every complex containing ``i`` -- and

            Z_i = 1 + e^{H_i},
            Z_a = prod_{i in a} (1 + e^{h_{i->a}}) - 1,

        the complex term being the sum over its interior with the one forbidden
        configuration, all members untaken, removed.  Then

            ln Z / N = sum_l n_l <ln Z_a> + <(1 - k) ln Z_i>,  n_l = <k_l>/c_l.
        """
        n = samples or self.size
        total = 0.0
        for l in range(self.L):
            idx = self.rng.integers(0, self.size, (n, int(self.c[l])))
            S = np.logaddexp(0.0, self.P[l][idx]).sum(axis=1)
            lnZa = np.log(np.expm1(S))          # ln(e^S - 1), accurate small S
            total += (self.deg[l] / self.c[l]) * float(np.mean(lnZa))
        k = np.zeros(n)
        H = np.full(n, -self.mu)
        for l in range(self.L):
            d = self._draw(n, l, excess=False)
            k += d
            tot = int(d.sum())
            if tot:
                draws = self.Q[l][self.rng.integers(0, self.size, tot)]
                ends = np.cumsum(d)
                cs = np.concatenate(([0.0], np.cumsum(draws)))
                H += cs[ends] - cs[ends - d]
        total += float(np.mean((1.0 - k) * np.logaddexp(0.0, H)))
        return total

    def entropy(self, samples=None):
        """``s = ln Z/N + mu rho``, the ground-state entropy density.

        Positive means an extensive number of minimum hitting sets and a
        replica-symmetric answer that can be believed; negative means the ansatz
        has failed and one-step breaking is required.  This is the criterion
        that replaces the hard-field instability, and unlike it, it holds for an
        arbitrary chygraph rather than only for regular hypergraphs.
        """
        return self.log_Z(samples) + self.mu * self.density(samples)

    def entropy_averaged(self, keep=200, samples=None):
        """Time-average :meth:`entropy` over ``keep`` further sweeps.

        Returns ``(mean, standard error)``.  Necessary for heterogeneous
        ensembles: ``sum_l n_l <ln Z_a>`` and ``mu rho`` are each ``O(mu)``
        while ``s`` is ``O(1)``, so a single snapshot carries the cancellation
        error of both.  Symmetry removes the sampling entirely in the regular
        case, which is why :meth:`entropy` is exact there.
        """
        vals = []
        for _ in range(keep):
            self.sweep()
            vals.append(self.entropy(samples))
        v = np.asarray(vals)
        return float(v.mean()), float(v.std() / np.sqrt(len(v)))

    def cover_size(self, samples=None):
        """Alias for :meth:`density`, matching the hard-field module's name."""
        return self.density(samples)


# ---------------------------------------------------------------------------
# Closed forms for the regular case (Mezard & Tarzia)
# ---------------------------------------------------------------------------

def regular_field(L, K, mu):
    """``h_RS = -mu/L - ((L-1)/L) ln(K-1)``, Eq. (11) of Ref. [MT]."""
    return -mu / L - (L - 1) / L * np.log(K - 1)


def regular_density(K):
    """``rho = 1/K``: every complex holds exactly one member of the set."""
    return 1.0 / K


def regular_entropy(L, K):
    """``s(mu -> inf)``, Eq. (13) of Ref. [MT].

    Positive means an extensive number of minimum hitting sets and a
    replica-symmetric solution that can be believed; negative means the ansatz
    has failed and one-step replica symmetry breaking is needed.
    """
    a = (L - 1) * (K - 1)
    return (a * np.log(K - 1) - (a - 1) * np.log(K)) / K
