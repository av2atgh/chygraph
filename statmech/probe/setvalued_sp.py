"""Set-valued SP for proper colouring, carried as a POPULATION of surveys.

setsp.py collapsed both indices: colour symmetry legitimately averages over
colours, but averaging over variables too is the scalar-closure error one level
up -- the map is nonlinear, so <f(v)> != f(<v>). It failed its c = 2 gate by
4-10%. Here each population member is its own (q+1)-vector.

The interior is precomputed once as a tensor, so the per-sample cost is a small
contraction rather than a double sum over subsets:

    T2[sj, so]      = P(emitted size so | |A_j| = sj)            c = 2
    T3[sj, sk, so]  = P(emitted size so | |A_j| = sj, |A_k| = sk) c = 3
"""
from itertools import combinations
from math import comb
import numpy as np

def _subsets(q):
    return [frozenset(T) for s in range(q + 1) for T in combinations(range(q), s)]

def _avail_for_0(others, q):
    ok = set()
    for g in range(q):
        rest = [A - {g} for A in others]
        if any(len(A) == 0 for A in rest):
            continue
        if len(rest) == 2:
            a, b = rest
            if len(a) == 1 and len(b) == 1 and a == b:
                continue
        ok.add(g)
    return ok

def interior_tensor(q, c):
    subs = _subsets(q)
    if c == 2:
        T = np.zeros((q + 1, q + 1))
        for A in subs:
            T[len(A), len(_avail_for_0([A], q))] += 1.0 / comb(q, len(A))
        return T
    T = np.zeros((q + 1, q + 1, q + 1))
    for A in subs:
        for B in subs:
            w = 1.0 / (comb(q, len(A)) * comb(q, len(B)))
            T[len(A), len(B), len(_avail_for_0([A, B], q))] += w
    return T

def _g_rows(V, q):
    """g[n, t] = P(member n's set contains a given t-set)."""
    C = np.array([[comb(q - t, s - t) / comb(q, s) if s >= t else 0.0
                   for s in range(q + 1)] for t in range(q + 1)])
    return V @ C.T

def sweep(V, q, c, kappa, T, rng):
    n = V.shape[0]
    if c == 2:
        M = V @ T
    else:
        a = V[rng.integers(0, n, n)]
        b = V[rng.integers(0, n, n)]
        M = np.einsum('ni,nj,ijk->nk', a, b, T)
    # Product of g over d ~ Poisson(kappa) ACTUAL draws, per sample. Replacing
    # this by exp(kappa (<g> - 1)) uses only the mean and destroys the
    # variable-to-variable fluctuation the population exists to carry.
    G = _g_rows(M, q)                                # (n, q+1), one per complex
    d = rng.poisson(kappa, n)
    tot = int(d.sum())
    gt = np.ones((n, q + 1))
    if tot:
        draws = np.log(np.maximum(G[rng.integers(0, n, tot)], 1e-300))
        cs = np.concatenate([np.zeros((1, q + 1)), np.cumsum(draws, axis=0)])
        ends = np.cumsum(d)
        gt = np.exp(cs[ends] - cs[ends - d])
    out = np.zeros_like(V)
    for u in range(q + 1):
        acc = np.zeros(n)
        for j in range(q - u + 1):
            acc += (-1) ** j * comb(q - u, j) * gt[:, u + j]
        out[:, u] = comb(q, u) * acc
    out = np.clip(out, 0, None)
    out[:, 0] = 0.0
    ssum = out.sum(axis=1, keepdims=True)
    return np.where(ssum > 0, out / np.maximum(ssum, 1e-300), 0.0)

def converge(q, c, kappa, size=4000, sweeps=300, seed=0):
    rng = np.random.default_rng(seed)
    T = interior_tensor(q, c)
    V = np.zeros((size, q + 1)); V[:, 1] = 0.9; V[:, 2] = 0.1
    for _ in range(sweeps):
        V = sweep(V, q, c, kappa, T, rng)
    return V

def onset(q, c, lo, hi, iters=12, **kw):
    for _ in range(iters):
        mid = 0.5 * (lo + hi)
        if converge(q, c, mid, **kw)[:, 1].mean() > 1e-4:
            hi = mid
        else:
            lo = mid
    return 0.5 * (lo + hi)

if __name__ == '__main__':
    print('GATE: c = 2, Poisson, branch onset against Mulet c_d')
    print('  q    onset    c_d    error')
    for q, cd in ((3, 4.42), (4, 8.27), (5, 12.67)):
        r = onset(q, 2, 2.0, 20.0)
        print(f'  {q}   {r:6.3f}   {cd:5.2f}   {100*(r-cd)/cd:+6.2f}%')


# --- complexity -------------------------------------------------------------
# Mirrors the form figs/colouring.py validates to 0.1% against Mulet's c_q:
#
#     Sigma = <ln Z_i>  -  (kappa/c) <ln Z_a>
#
# Z_i is the site normalisation using ALL kappa complexes (Poisson, not excess);
# Z_a is the weight that a complex is satisfiable given its members' cavity
# surveys. At c = 2, Z_a = 1 - v_i[1] v_j[1] / q, which is the scalar code's
# edge term with e = v[1], so the two must agree there. That is the gate.

def _sat_c2(vi, vj, q):
    """P(an edge is properly colourable) = 1 - P(both ends forced alike)."""
    return 1.0 - vi[:, 1] * vj[:, 1] / q


def _sat_tensor(q, c):
    """S[s1,..,sc] = P(c members with those available-set sizes admit an SDR)."""
    subs = _subsets(q)
    idx = {}
    for A in subs:
        idx.setdefault(len(A), []).append(A)
    S = np.zeros((q + 1,) * c)
    from itertools import product as iproduct
    for sizes in iproduct(range(q + 1), repeat=c):
        if any(s == 0 for s in sizes):
            continue
        tot = ok = 0
        for combo in iproduct(*[idx[s] for s in sizes]):
            tot += 1
            # SDR exists iff distinct colours can be picked, Hall's condition
            found = False
            for pick in iproduct(*[sorted(A) for A in combo]):
                if len(set(pick)) == c:
                    found = True
                    break
            ok += found
        S[sizes] = ok / tot if tot else 0.0
    return S


def site_and_sigma(q, c, kappa, size=4000, sweeps=300, seed=0):
    rng = np.random.default_rng(seed)
    T = interior_tensor(q, c)
    V = converge(q, c, kappa, size=size, sweeps=sweeps, seed=seed)
    if V[:, 1].mean() < 1e-6:
        return 0.0, 0.0
    n = V.shape[0]
    # site: intersect the FULL degree, Poisson(kappa)
    M = (V @ T) if c == 2 else np.einsum(
        'ni,nj,ijk->nk', V[rng.integers(0, n, n)], V[rng.integers(0, n, n)], T)
    G = _g_rows(M, q)
    d = rng.poisson(kappa, n)
    tot = int(d.sum())
    gt = np.ones((n, q + 1))
    if tot:
        draws = np.log(np.maximum(G[rng.integers(0, n, tot)], 1e-300))
        cs = np.concatenate([np.zeros((1, q + 1)), np.cumsum(draws, axis=0)])
        ends = np.cumsum(d)
        gt = np.exp(cs[ends] - cs[ends - d])
    P0 = sum((-1) ** j * comb(q, j) * gt[:, j] for j in range(q + 1))
    site = np.log(np.maximum(1.0 - P0, 1e-300)).mean()
    if c == 2:
        Za = _sat_c2(V[rng.integers(0, n, n)], V[rng.integers(0, n, n)], q)
    else:
        Sm = _sat_tensor(q, c)
        parts = [V[rng.integers(0, n, n)] for _ in range(c)]
        Za = np.einsum('ni,nj,nk,ijk->n', *parts, Sm)
    # Three-term Bethe count, as in the hypergraph case. The two-term form in
    # figs/colouring.py is a c = 2 collapse and does not generalise: mirroring
    # it here made Sigma large and negative with no crossing.
    #   Z_ia = 1 - P(complex forbids gamma AND member forced to gamma) summed
    #        = 1 - (1 - g(1)) v[1]
    g1 = _g_rows(M, q)[:, 1]
    Zia = 1.0 - (1.0 - g1[rng.integers(0, n, n)]) * V[rng.integers(0, n, n), 1]
    lnZa = np.log(np.maximum(Za, 1e-300)).mean()
    lnZia = np.log(np.maximum(Zia, 1e-300)).mean()
    return site + (kappa / c) * lnZa - kappa * lnZia, V[:, 1].mean()


def threshold(q, c, lo, hi, iters=10, **kw):
    assert site_and_sigma(q, c, lo, **kw)[0] > 0, (q, c, lo)
    for _ in range(iters):
        mid = 0.5 * (lo + hi)
        if site_and_sigma(q, c, mid, **kw)[0] > 0:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)
