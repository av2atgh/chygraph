"""SP on the three-level chygraph of Sec. 13.6, against its CNF flattening.

Messages reaching a variable are of two kinds and are kept as two populations:
  eta  a layer-2 clause D warns a plain variable member
  del  a layer-1 clause C, forced false by its parent D, passes that down
Both arrive with a random sign, so each splits Poisson(.../2) by direction.

  pi     = Pi^u / (Pi^u + Pi^s + Pi^0)     normalised, contradiction discarded
  gamma  = 1 - prod_{k1} (1 - pi)          some member of C is forced to satisfy C
  eta'   = [prod_{k2-1} pi] * gamma        D needs its other members violating
  del'   = prod_{k2} pi                    D forces C false
"""
import numpy as np

def _p1m(pop, n, rng):
    out = np.ones(n.size)
    tot = int(n.sum())
    if tot:
        d = np.minimum(pop[rng.integers(0, pop.size, tot)], 1.0 - 1e-12)
        c = np.concatenate(([0.0], np.cumsum(np.log1p(-d))))
        ends = np.cumsum(n)
        out = np.exp(c[ends] - c[ends - n])
    return out

def _pi(eta, del_, k1, k2, alpha, rng, size, full=False):
    """pi^u for a cavity variable (full=False) or the site weights (full=True)."""
    lD, lC = k2 * alpha / 2.0, k1 * alpha / 2.0
    Pu = _p1m(eta, rng.poisson(lD, size), rng) * _p1m(del_, rng.poisson(lC, size), rng)
    Ps = _p1m(eta, rng.poisson(lD, size), rng) * _p1m(del_, rng.poisson(lC, size), rng)
    wu, ws, w0 = (1 - Ps) * Pu, (1 - Pu) * Ps, Pu * Ps
    Z = np.maximum(wu + ws + w0, 1e-300)
    return Z if full else np.clip(wu / Z, 0.0, 1.0)

def nested_converge(k1, k2, alpha, size=20000, sweeps=200, seed=0):
    rng = np.random.default_rng(seed)
    eta = rng.uniform(0.3, 0.5, size)
    dl  = rng.uniform(0.3, 0.5, size)
    for _ in range(sweeps):
        pis = [_pi(eta, dl, k1, k2, alpha, rng, size) for _ in range(k2)]
        gam = 1.0 - np.prod([1.0 - _pi(eta, dl, k1, k2, alpha, rng, size)
                             for _ in range(k1)], axis=0)
        eta = np.prod(pis[:k2 - 1], axis=0) * gam if k2 > 1 else gam
        dl  = np.prod(pis, axis=0)
    return eta, dl, rng

def flat_converge(k, alpha, size=20000, sweeps=200, seed=0):
    """Flat k-SAT at clause density alpha, same conventions."""
    rng = np.random.default_rng(seed)
    eta = rng.uniform(0.3, 0.5, size)
    lam = k * alpha / 2.0
    for _ in range(sweeps):
        new = np.ones(size)
        for _ in range(k - 1):
            Pu = _p1m(eta, rng.poisson(lam, size), rng)
            Ps = _p1m(eta, rng.poisson(lam, size), rng)
            wu, ws, w0 = (1 - Ps) * Pu, (1 - Pu) * Ps, Pu * Ps
            new *= np.clip(wu / np.maximum(wu + ws + w0, 1e-300), 0.0, 1.0)
        eta = new
    return eta, rng

def onset(f, lo, hi, iters=10, **kw):
    for _ in range(iters):
        m = 0.5 * (lo + hi)
        if f(m, **kw) > 1e-4:
            hi = m
        else:
            lo = m
    return 0.5 * (lo + hi)

if __name__ == '__main__':
    print('k1 = 1: the relay. nested must equal flat (k2+1)-SAT.')
    for k2 in (2, 3):
        for a in (3.0, 6.0, 12.0):
            n = nested_converge(1, k2, a)[0].mean()
            f = flat_converge(k2 + 1, a)[0].mean()
            print(f'  k2={k2} alpha={a:5.1f}   nested <eta>={n:.5f}   '
                  f'flat {k2+1}-SAT <eta>={f:.5f}   diff={abs(n-f):.2e}')


# --- complexity -------------------------------------------------------------
#
# Hypothesis: a layer-1 complex is a sub-expression, not a constraint, so it
# contributes no factor of its own -- its interior is already summed into the
# message D receives. Then, per variable,
#
#   Sigma = <ln Z_i> + alpha <ln Z_D> - k2 alpha <ln Z_iD> - k1 alpha <ln Z_iC>
#
# At k1 = 1 this must reduce to flat (k2+1)-SAT identically: the two inclusion
# channels merge to (k2+1) alpha, and Z_D = 1 - [prod_k2 pi] gamma collapses to
# 1 - prod_{k2+1} pi because gamma = pi for a one-member sub-clause.

def nested_sigma(k1, k2, alpha, size=20000, sweeps=200, seed=0, nsamp=300000):
    eta, dl, rng = nested_converge(k1, k2, alpha, size, sweeps, seed)
    if eta.mean() < 1e-9 and dl.mean() < 1e-9:
        return 0.0
    kw = dict(k1=k1, k2=k2, alpha=alpha, rng=rng, size=nsamp)
    site = np.log(_pi(eta, dl, full=True, **kw)).mean()
    pis = [_pi(eta, dl, **kw) for _ in range(k2)]
    gam = 1.0 - np.prod([1.0 - _pi(eta, dl, **kw) for _ in range(k1)], axis=0)
    ZD = np.log(np.maximum(1.0 - np.prod(pis, axis=0) * gam, 1e-300)).mean()
    e = eta[rng.integers(0, eta.size, nsamp)]
    d = dl[rng.integers(0, dl.size, nsamp)]
    ZiD = np.log(np.maximum(1.0 - e * _pi(eta, dl, **kw), 1e-300)).mean()
    ZiC = np.log(np.maximum(1.0 - d * _pi(eta, dl, **kw), 1e-300)).mean()
    return site + alpha * ZD - k2 * alpha * ZiD - k1 * alpha * ZiC


def flat_sigma(k, alpha, size=20000, sweeps=200, seed=0, nsamp=300000):
    eta, rng = flat_converge(k, alpha, size, sweeps, seed)
    if eta.mean() < 1e-9:
        return 0.0
    lam = k * alpha / 2.0
    def piu(n):
        Pu = _p1m(eta, rng.poisson(lam, n), rng)
        Ps = _p1m(eta, rng.poisson(lam, n), rng)
        wu, ws, w0 = (1 - Ps) * Pu, (1 - Pu) * Ps, Pu * Ps
        return np.clip(wu / np.maximum(wu + ws + w0, 1e-300), 0, 1), wu + ws + w0
    _, Zi = piu(nsamp)
    site = np.log(np.maximum(Zi, 1e-300)).mean()
    prod = np.ones(nsamp)
    for _ in range(k):
        prod *= piu(nsamp)[0]
    ZD = np.log(np.maximum(1.0 - prod, 1e-300)).mean()
    e = eta[rng.integers(0, eta.size, nsamp)]
    Zed = np.log(np.maximum(1.0 - e * piu(nsamp)[0], 1e-300)).mean()
    return site + alpha * ZD - k * alpha * Zed
