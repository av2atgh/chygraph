"""Question 1: community sizes of the degree-correlated ensemble under the
nested degree-corrected SBM, against n and r.

For each (tau, r, n): fit minimize_nested_blockmodel_dl(deg_corr=True), take
the bottom level, and record the block-size distribution, the internal density
of the blocks, the fraction of edges between blocks, and how tree-like the
block graph is (edges of the block multigraph vs blocks-1; its 2-core).
Cached as results/sbm_sizes.json.
"""
import json
import sys
import time
from collections import Counter
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from ensemble import graph  # noqa: E402

OUT = Path(__file__).resolve().parent / 'results' / 'sbm_sizes.json'


def fit(n, edges, seed):
    import graph_tool.all as gt
    gt.seed_rng(seed)
    np.random.seed(seed)
    g = gt.Graph(directed=False)
    g.add_vertex(n)
    g.add_edge_list(edges)
    st = gt.minimize_nested_blockmodel_dl(g, state_args=dict(deg_corr=True))
    b = np.asarray(st.get_bs()[0], dtype=int)
    b = np.unique(b, return_inverse=True)[1]      # contiguous labels
    return b, [len(set(np.asarray(x).tolist())) for x in st.get_bs()], float(st.entropy())


def block_graph_stats(b, edges):
    B = b.max() + 1
    bi, bj = b[edges[:, 0]], b[edges[:, 1]]
    inter = bi != bj
    m_inter = int(inter.sum())
    # simple block graph
    pairs = np.unique(np.sort(np.stack([bi[inter], bj[inter]], 1), axis=1), axis=0)
    # 2-core of the simple block graph
    deg = Counter()
    adj = [set() for _ in range(B)]
    for a, c in pairs.tolist():
        adj[a].add(c); adj[c].add(a)
    alive = np.ones(B, bool)
    stack = [v for v in range(B) if len(adj[v]) <= 1]
    d = np.array([len(a) for a in adj])
    while stack:
        v = stack.pop()
        if not alive[v] or d[v] > 1:
            continue
        alive[v] = False
        for u in adj[v]:
            if alive[u]:
                d[u] -= 1
                if d[u] <= 1:
                    stack.append(u)
    core2 = int((alive & (d > 0)).sum())
    return dict(blocks=int(B), m_inter=m_inter, m_inter_frac=m_inter / len(edges),
                block_edges=int(len(pairs)), block_cycle_excess=int(len(pairs) - (B - 1)),
                block_2core_frac=core2 / B)


def one(tau, r, n, seed):
    e = graph(tau, r, n, seed)
    t0 = time.time()
    b, levels, dl = fit(n, e, seed)
    t = time.time() - t0
    sizes = np.bincount(b)
    sizes = sizes[sizes > 0]
    # internal density per block
    same = b[e[:, 0]] == b[e[:, 1]]
    m_in = np.bincount(b[e[same, 0]], minlength=len(sizes))
    deg = np.bincount(e.ravel(), minlength=n)
    kmean = np.bincount(b, weights=deg) / sizes
    kstd = np.sqrt(np.maximum(np.bincount(b, weights=deg ** 2) / sizes - kmean ** 2, 0))
    per_block = sorted(zip(sizes.tolist(), kmean.round(2).tolist(), kstd.round(2).tolist(),
                           (m_in / np.maximum(sizes * (sizes - 1) / 2, 1)).round(4).tolist()),
                       reverse=True)[:12]
    dens = np.where(sizes > 1, m_in / np.maximum(sizes * (sizes - 1) / 2, 1), 0)
    w = sizes / sizes.sum()                       # node-weighted
    row = dict(tau=tau, r=r, n=n, seed=seed, m=int(len(e)), kbar=2 * len(e) / n,
               levels=levels, dl=dl, secs=t, per_block=per_block,
               size_mean=float(sizes.mean()), size_median=float(np.median(sizes)),
               size_max=int(sizes.max()),
               size_node_mean=float((w * sizes).sum()),
               size_node_median=float(sizes[np.searchsorted(np.cumsum(np.sort(sizes) / sizes.sum()), 0.5)]),
               frac_nodes_le50=float(w[sizes <= 50].sum()),
               frac_nodes_le100=float(w[sizes <= 100].sum()),
               dens_node_mean=float((w * dens).sum()),
               frac_edges_internal=float(same.mean()),
               **block_graph_stats(b, e))
    np.save(OUT.parent / f'blocks_tau{tau}_r{r}_n{n}_s{seed}.npy', b)
    return row


def main():
    taus = [2.5]
    rs = [0.0, 0.5, 0.8, 1.0]
    ns = [int(x) for x in sys.argv[1:]] or [2000, 5000, 10000, 20000, 50000]
    rows = json.loads(OUT.read_text()) if OUT.exists() else []
    done = {(x['tau'], x['r'], x['n'], x['seed']) for x in rows}
    for n in ns:
        for tau in taus:
            for r in rs:
                if (tau, r, n, 0) in done:
                    continue
                row = one(tau, r, n, 0)
                rows.append(row)
                OUT.write_text(json.dumps(rows, indent=1))
                print(f"tau={tau} r={r} n={n:6d} B={row['blocks']:5d} "
                      f"size mean={row['size_mean']:6.1f} node-median={row['size_node_median']:5.0f} "
                      f"max={row['size_max']:5d} <=50:{row['frac_nodes_le50']:.2f} "
                      f"int.edges={row['frac_edges_internal']:.2f} dens={row['dens_node_mean']:.3f} "
                      f"blockgraph excess={row['block_cycle_excess']} 2core={row['block_2core_frac']:.2f} "
                      f"levels={row['levels']} [{row['secs']:.0f}s]", flush=True)
                print('   blocks (size, <k>, sd k, density):', row['per_block'][:8], flush=True)


if __name__ == '__main__':
    main()
