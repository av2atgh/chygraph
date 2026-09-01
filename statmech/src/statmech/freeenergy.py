"""The Bethe free energy on a chygraph (WP6).

Eq.~(12) of Vazquez & Weigt is a site-plus-link sum: one term per vertex with
weight ``d - 1`` and one per edge.  It has no counterpart in ``chygraph``,
because a percolation threshold never needs it.  Anything variational does.

On a chygraph the sum runs over complexes and nodes.  Parameterising the
messages as ``nu_{i->a}(sigma) ~ exp(h_{i->a} sigma)`` and
``nu_{a->i}(sigma) ~ exp(u_{a->i} sigma)``, the overlap term of the general
Bethe free energy collapses --- because ``h_{i->a} + u_{a->i}`` is the *same*
full field ``h_i`` for every complex ``a`` containing ``i`` --- and what is left
is one term per complex and ``1 - k_i`` per node:

    -beta f = sum_l n_l <ln Z_a^(l)> + <(1 - k) ln Z_i>,     n_l = <k_l> / c_l

with ``Z_a`` the exact partition function *inside* a complex given the cavity
fields of its members, and ``Z_i = 2 cosh(h_i)``.

Those weights are not a coincidence.  They are exactly the Mobius counting
numbers of :mod:`statmech.region` in the treelike case: ``1`` on every
complex and ``1 - k_v`` on every node.  WP5 computes the counting, WP6
evaluates the free energy with it, and where complexes overlap both are wrong
in the same way.

At zero cavity field this reduces to a closed form that serves as the test:

    -beta f_para = ln 2 + sum_l (<k_l> / c_l) ln( Z_c^(l) / 2^{c_l} )

with ``Z_c`` the isolated-complex partition function.  For a graph that is the
textbook ``ln 2 + (c/2) ln cosh(beta J)``.
"""

import numpy as np
from scipy.special import logsumexp

from statmech.population import CavityPopulation


class BetheFreeEnergy:
    """Free energy per spin of an Ising chygraph, from cavity populations.

    Args:
        population: a converged :class:`~statmech.population.CavityPopulation`.
    """

    def __init__(self, population):
        self.pop = population
        self._full = []
        for c in population.c:
            m = int(c)
            bits = ((np.arange(2 ** m)[:, None] >> np.arange(m)) & 1)
            self._full.append(1.0 - 2.0 * bits)        # (2^c, c)

    # -- the two terms ------------------------------------------------------

    def log_Z_complex(self, l, n=None):
        """``<ln Z_a>`` for a layer-``l`` complex, over its members' cavity fields."""
        p = self.pop
        n = n or p.size
        c, bJ, cfg = int(p.c[l]), p.beta_J[l], self._full[l]
        S = cfg.sum(axis=1)
        pair = bJ * (S ** 2 - c) / 2.0                  # (2^c,)
        h = p.P[l][p.rng.integers(0, p.size, (n, c))]   # (n, c)
        return float(np.mean(logsumexp(h @ cfg.T + pair, axis=1)))

    def log_Z_node(self, n=None):
        """``<(1 - k) ln Z_i>`` over the node's full field and complex count.

        Split as ``ln Z_i = ln 2 + ln cosh h`` and take the first piece
        exactly: ``<(1-k) ln 2> = (1 - sum_l <k_l>) ln 2`` needs no sampling,
        and sampling it costs a bias of order ``sqrt(<k>/n) ln 2`` -- about
        ``5e-3`` at ``<k> = 6`` and ``n = 8e4``, which is larger than the
        free-energy differences this module exists to resolve.  Only the second
        piece is sampled, and it vanishes identically at zero field, so the
        paramagnetic value comes out exact.
        """
        p = self.pop
        n = n or p.size
        h = np.zeros(n)
        k = np.zeros(n)
        for l in range(p.L):
            d = p.rng.poisson(p.means[l], n)
            k += d
            total = int(d.sum())
            if total:
                draws = p.Q[l][p.rng.integers(0, p.size, total)]
                ends = np.cumsum(d)
                csum = np.concatenate(([0.0], np.cumsum(draws)))
                h += csum[ends] - csum[ends - d]
        exact = (1.0 - float(p.means.sum())) * np.log(2.0)
        sampled = float(np.mean((1.0 - k) * np.log(np.cosh(h))))
        return exact + sampled

    def minus_beta_f(self, n=None):
        """``-beta f`` per node."""
        p = self.pop
        out = self.log_Z_node(n)
        for l in range(p.L):
            out += (p.means[l] / float(p.c[l])) * self.log_Z_complex(l, n)
        return out


# ---------------------------------------------------------------------------
# Closed forms
# ---------------------------------------------------------------------------

def paramagnetic(cardinalities, means, beta_J):
    """``-beta f`` with every cavity field zero, in closed form.

    ``ln 2 + sum_l (<k_l>/c_l) ln(Z_c / 2^{c_l})`` with ``Z_c`` the partition
    function of one isolated complex.  Exact at any temperature in the
    paramagnetic phase, and the high-temperature limit of the ordered branch.
    """
    c = np.asarray(cardinalities, dtype=int)
    means = np.asarray(means, dtype=float)
    bj = np.asarray(beta_J, dtype=float)
    bj = np.full(len(c), float(bj)) if bj.ndim == 0 else bj
    out = np.log(2.0)
    for l, cl in enumerate(c):
        m = int(cl)
        bits = ((np.arange(2 ** m)[:, None] >> np.arange(m)) & 1)
        cfg = 1.0 - 2.0 * bits
        S = cfg.sum(axis=1)
        lnZ = logsumexp(bj[l] * (S ** 2 - m) / 2.0)
        out += (means[l] / m) * (lnZ - m * np.log(2.0))
    return float(out)


def graph_paramagnetic(mean, beta_J):
    """``ln 2 + (c/2) ln cosh(beta J)``, the textbook Bethe result."""
    return float(np.log(2.0) + 0.5 * mean * np.log(np.cosh(beta_J)))


# ---------------------------------------------------------------------------
# Locating the transition thermodynamically
# ---------------------------------------------------------------------------

def free_energy_gap(cardinalities, means, beta_J, sweeps=250, size=60_000,
                    seed=0):
    """``(-beta f)_ordered - (-beta f)_paramagnetic``.

    Zero above ``T_c`` where the only solution is the paramagnet, positive
    below it since the ordered branch has the lower free energy.  Locating the
    transition this way uses neither the linearisation of WP1 nor the order
    parameter of WP4.
    """
    pop = CavityPopulation(cardinalities, means, beta_J, size=size,
                           seed=seed).run(sweeps)
    ordered = BetheFreeEnergy(pop).minus_beta_f()
    return ordered - paramagnetic(cardinalities, means, beta_J)
