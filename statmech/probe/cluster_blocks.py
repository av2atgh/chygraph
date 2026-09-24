"""Is a solution cluster itself a chygraph?  Exhaustive 3-SAT, small N.

A cluster of solutions is a set of assignments connected by single flips.
The one-step ansatz treats it as one pure state with one field per variable.
This probe asks what the cluster looks like from inside: which variables are
frozen (the same in every solution of the cluster), and whether the free
ones fall into blocks that vary *independently* -- the cluster being exactly
the Cartesian product of its projections onto the blocks.  If it is, the
cluster is a chygraph: a frozen core plus complexes (the blocks) whose
interiors are enumerated and whose product is the whole.

Method.  Random 3-SAT at N = 16..24 variables, clause density alpha, every
assignment tried (2^N x M, vectorised).  Solutions are clustered at Hamming
distance one.  For each cluster with at least two solutions: the frozen
set; on the free set F, the dependency graph -- i and j are dependent when
the projection of the cluster onto (i, j) is not the product of the two
marginals -- and its connected components are the candidate blocks; the
cluster is a product iff its size equals the product of the block
projection sizes.  Control: a uniformly random subset of {0,1}^|F| of the
same size, put through the same test; a random set of that size is
essentially never a product of two or more blocks.

    python probe/cluster_blocks.py [seeds]

Writes probe/results/cluster_blocks.csv, one row per cluster.
"""

import csv
import sys
import time
from collections import defaultdict
from pathlib import Path

import numpy as np

OUT = Path(__file__).parent / 'results' / 'cluster_blocks.csv'
SIZES = (16, 18, 20, 22, 24)


def random_3sat(n, alpha, rng):
    m = int(round(alpha * n))
    clauses = np.empty((m, 3), dtype=np.int64)
    for c in range(m):
        clauses[c] = rng.choice(n, 3, replace=False)
    signs = rng.integers(0, 2, (m, 3)).astype(bool)
    return clauses, signs


def solutions(n, clauses, signs):
    """All satisfying assignments as an (S, n) boolean array."""
    total = 1 << n
    sols = []
    chunk = 1 << 20
    bits = 1 << np.arange(n, dtype=np.int64)
    for start in range(0, total, chunk):
        idx = np.arange(start, min(total, start + chunk), dtype=np.int64)
        x = (idx[:, None] & bits[None, :]) != 0         # (chunk, n)
        lit = x[:, clauses] == signs[None, :, :]        # literal true?
        sat = lit.any(axis=2).all(axis=1)
        sols.append(x[sat])
    return np.vstack(sols) if sols else np.zeros((0, n), bool)


def clusters(sols):
    """Connected components at Hamming distance one, by union-find on a hash."""
    s = sols.shape[0]
    n = sols.shape[1]
    weights = 1 << np.arange(n, dtype=np.int64)
    keys = (sols.astype(np.int64) * weights).sum(1)
    index = {int(k): i for i, k in enumerate(keys)}
    parent = np.arange(s)

    def find(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    for i, k in enumerate(keys):
        for b in range(n):
            j = index.get(int(k) ^ (1 << b))
            if j is not None and j > i:
                ri, rj = find(i), find(j)
                if ri != rj:
                    parent[rj] = ri
    groups = defaultdict(list)
    for i in range(s):
        groups[find(i)].append(i)
    return [sols[g] for g in groups.values()]


def block_decomposition(S):
    """Frozen set, blocks of free variables, and whether S is their product.

    S is a (size, n) boolean array of assignments.  Returns (frozen, blocks,
    is_product, block_sizes) with block_sizes the number of distinct
    projections onto each block.
    """
    size, n = S.shape
    frozen = [i for i in range(n) if S[:, i].all() or (~S[:, i]).all()]
    free = [i for i in range(n) if i not in frozen]
    f = len(free)
    # pairwise dependency on the free variables
    dep = np.zeros((f, f), bool)
    for a in range(f):
        xa = S[:, free[a]]
        pa = xa.mean()
        for b in range(a + 1, f):
            xb = S[:, free[b]]
            joint = np.array([[(~xa & ~xb).any(), (~xa & xb).any()],
                              [(xa & ~xb).any(), (xa & xb).any()]])
            # product iff every combination of seen marginal values occurs
            if not joint.all():
                dep[a, b] = dep[b, a] = True
    # components
    seen = [False] * f
    blocks = []
    for a in range(f):
        if seen[a]:
            continue
        stack, comp = [a], []
        seen[a] = True
        while stack:
            u = stack.pop()
            comp.append(free[u])
            for v in np.nonzero(dep[u])[0]:
                if not seen[v]:
                    seen[v] = True
                    stack.append(v)
        blocks.append(sorted(comp))
    sizes = []
    for B in blocks:
        proj = np.unique(S[:, B], axis=0)
        sizes.append(len(proj))
    is_product = int(np.prod(sizes)) == size
    return frozen, blocks, is_product, sizes


def finest_product(S, blocks):
    """The finest partition of the free variables over which S is a product.

    Product decompositions of a set of tuples are closed under common
    refinement (if S = A x B over one cut and C x D over another, it is a
    product over the four intersections), so a finest one exists and is
    unique.  Every product partition is coarser than the pairwise one, so it
    is a merging of the candidate `blocks`: a union U of candidates is
    *separable* when |S| = |proj_U(S)| |proj_{F\\U}(S)|, and the
    inclusion-minimal separable unions are the finest blocks.  Enumerates the
    2^k unions for k candidates, k <= 12; above that, returns one block.
    """
    k = len(blocks)
    size = S.shape[0]
    if k <= 1:
        return [list(b) for b in blocks]
    if k > 12:
        return [sorted(v for b in blocks for v in b)]
    free = sorted(v for b in blocks for v in b)
    cache = {}

    def nproj(cols):
        key = tuple(cols)
        if key not in cache:
            cache[key] = len(np.unique(S[:, list(cols)], axis=0)) if cols else 1
        return cache[key]

    separable = []
    for mask in range(1, 1 << k):
        if mask == (1 << k) - 1:
            continue
        U = sorted(v for i, b in enumerate(blocks) if mask >> i & 1 for v in b)
        Uc = [v for v in free if v not in set(U)]
        if nproj(U) * nproj(Uc) == size:
            separable.append(mask)
    minimal = [m for m in separable if not any(o != m and (o & m) == o for o in separable)]
    if not minimal:
        return [free]
    covered = 0
    out = []
    for m in sorted(minimal, key=lambda m: bin(m).count('1')):
        if covered & m:
            continue
        out.append(sorted(v for i, b in enumerate(blocks) if m >> i & 1 for v in b))
        covered |= m
    rest = [v for v in free if v not in {x for b in out for x in b}]
    if rest:
        out.append(rest)
    return out


def random_control(size, f, rng):
    """A uniformly random subset of {0,1}^f with `size` elements."""
    if f == 0:
        return np.zeros((size, 0), bool)
    total = 1 << f
    if size >= total:
        idx = np.arange(total)
    else:
        idx = rng.choice(total, size, replace=False)
    bits = 1 << np.arange(f, dtype=np.int64)
    return (idx[:, None] & bits[None, :]) != 0


def main(alphas=(3.5, 3.8, 4.0), seeds=range(30)):
    rows = []
    for n in SIZES:
        for alpha in alphas:
            for seed in seeds:
                rng = np.random.default_rng(1000 * n + seed)
                clauses, signs = random_3sat(n, alpha, rng)
                t0 = time.time()
                sols = solutions(n, clauses, signs)
                if sols.shape[0] == 0:
                    continue
                cl = clusters(sols)
                for ci, S in enumerate(cl):
                    if S.shape[0] < 2:
                        continue
                    frozen, blocks, prod, sizes = block_decomposition(S)
                    f = S.shape[1] - len(frozen)
                    fine = finest_product(S, blocks)
                    C = random_control(S.shape[0], f, rng)
                    _, cblocks, cprod, csizes = block_decomposition(C) if f else ([], [], True, [])
                    rows.append(dict(
                        n=n, seed=seed, alpha=alpha, nsol=sols.shape[0], nclusters=len(cl),
                        cluster=ci, size=S.shape[0], frozen=len(frozen), free=f,
                        nblocks=len(blocks), product=int(prod),
                        largest_block=max(len(b) for b in blocks) if blocks else 0,
                        entropy=float(np.log(S.shape[0])),
                        block_entropy=float(sum(np.log(s) for s in sizes)),
                        ctrl_nblocks=len(cblocks), ctrl_product=int(cprod),
                        finest_nblocks=len(fine), finest_largest=max(len(b) for b in fine) if fine else 0,
                        ctrl_finest_nblocks=len(finest_product(C, cblocks)) if f else 0,
                    ))
                print(f'n={n} alpha={alpha} seed={seed} solutions={sols.shape[0]} '
                      f'clusters={len(cl)} ({time.time() - t0:.1f}s)', flush=True)
            with open(OUT, 'w', newline='') as fh:
                w = csv.DictWriter(fh, fieldnames=list(rows[0]))
                w.writeheader()
                w.writerows(rows)


if __name__ == '__main__':
    seeds = range(int(sys.argv[1])) if len(sys.argv) > 1 else range(30)
    main(seeds=seeds)
