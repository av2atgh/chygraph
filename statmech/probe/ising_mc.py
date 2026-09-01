"""Monte Carlo test of the clustering/T_c claim (referee 4.8).

The theory says: at a genuinely matched degree distribution, clustering
*lowers* T_c.  The cleanest comparison is a 4-regular random graph against a
network of two triangles per node -- both have every vertex at degree exactly 4,
both have 2 edges per node, and they differ only in whether those edges close
into triangles.  Eqs. (6) and (12) give

    links      3 t = 1                 -> t_c = 1/3,        T_c = 2.8854
    triangles  2 t/(1-t+t^2) = 1       -> t_c = 0.381966,   T_c = 2.4853

so the prediction is a 13.9% *reduction*.  Nothing in this file uses the cavity
equations; it is a direct simulation.

Wolff cluster updates, Binder cumulant U = 1 - <m^4>/(3<m^2>^2) at several sizes.
U is size-independent at T_c, so the crossing of the curves locates it.
"""

import sys
from pathlib import Path

import numpy as np
from numba import njit

SEEDS = 4


def regular_graph(n, deg, rng):
    """Simple random regular graph by stub pairing with rejection."""
    for _ in range(200):
        stubs = np.repeat(np.arange(n), deg)
        rng.shuffle(stubs)
        a, b = stubs[0::2], stubs[1::2]
        if (a == b).any():
            continue
        key = np.minimum(a, b) * n + np.maximum(a, b)
        if len(np.unique(key)) != len(key):
            continue
        return a, b
    raise RuntimeError('could not build a simple regular graph')


def triangle_network(n, per_node, rng):
    """`per_node` triangles at every node; degree = 2 * per_node."""
    for _ in range(200):
        stubs = np.repeat(np.arange(n), per_node)
        rng.shuffle(stubs)
        tri = stubs.reshape(-1, 3)
        if any((tri[:, i] == tri[:, j]).any()
               for i, j in ((0, 1), (0, 2), (1, 2))):
            continue
        a = np.concatenate([tri[:, 0], tri[:, 1], tri[:, 0]])
        b = np.concatenate([tri[:, 1], tri[:, 2], tri[:, 2]])
        key = np.minimum(a, b) * n + np.maximum(a, b)
        if len(np.unique(key)) != len(key):
            continue
        return a, b
    raise RuntimeError('could not build a simple triangle network')


def csr(n, a, b):
    src = np.concatenate((a, b))
    dst = np.concatenate((b, a))
    o = np.argsort(src, kind='stable')
    src, dst = src[o], dst[o]
    indptr = np.zeros(n + 1, np.int64)
    np.add.at(indptr, src + 1, 1)
    np.cumsum(indptr, out=indptr)
    return indptr, dst.astype(np.int64)


@njit(cache=True)
def _wolff(indptr, indices, spin, padd, nsweep, seed):
    np.random.seed(seed)
    n = spin.size
    stack = np.empty(n, np.int64)
    m2 = 0.0
    m4 = 0.0
    cnt = 0
    for s in range(nsweep):
        flipped = 0
        while flipped < n:                      # one sweep ~ n flipped spins
            i = np.random.randint(n)
            old = spin[i]
            spin[i] = -old
            top = 0
            stack[0] = i
            top = 1
            flipped += 1
            while top > 0:
                top -= 1
                v = stack[top]
                for e in range(indptr[v], indptr[v + 1]):
                    w = indices[e]
                    if spin[w] == old and np.random.random() < padd:
                        spin[w] = -old
                        stack[top] = w
                        top += 1
                        flipped += 1
        if s >= nsweep // 3:                    # discard the first third
            m = 0.0
            for k in range(n):
                m += spin[k]
            m = abs(m) / n
            m2 += m * m
            m4 += m * m * m * m
            cnt += 1
    return m2 / cnt, m4 / cnt


def binder(indptr, indices, n, T, nsweep, seed):
    padd = 1.0 - np.exp(-2.0 / T)
    spin = np.ones(n, np.int64)
    m2, m4 = _wolff(indptr, indices, spin, padd, nsweep, seed)
    return 1.0 - m4 / (3.0 * m2 * m2), np.sqrt(m2)


def run(kind, label, sizes, temps, nsweep=1600):
    print(f"\n{label}")
    print(f"{'T':>7}" + "".join(f"{('U n=%d' % n):>12}" for n in sizes))
    table = {}
    for T in temps:
        row = f"{T:>7.3f}"
        for n in sizes:
            us = []
            for s in range(SEEDS):
                rng = np.random.default_rng(1000 * s + n)
                a, b = (regular_graph(n, 4, rng) if kind == 'links'
                        else triangle_network(n, 2, rng))
                ip, idx = csr(n, a, b)
                u, _ = binder(ip, idx, n, T, nsweep, 7 * s + 13)
                us.append(u)
            table[(T, n)] = (float(np.mean(us)), float(np.std(us) / np.sqrt(SEEDS)))
            row += f"{np.mean(us):>9.4f}±{np.std(us)/np.sqrt(SEEDS):.3f}"
        print(row, flush=True)
    return table


if __name__ == '__main__':
    # n must be divisible by 3 so that 2n stubs make whole triangles
    sizes = (3000, 12000)
    print("Binder cumulant; T_c is where the curves for different n cross.")
    print("theory: links T_c = 2.8854, triangles T_c = 2.4853 (13.9% lower)")
    run('links', 'links (4-regular graph)', sizes,
        (2.60, 2.75, 2.85, 2.89, 2.95, 3.10))
    run('tri', 'triangles (two per node, degree 4)', sizes,
        (2.20, 2.35, 2.45, 2.49, 2.55, 2.70))
