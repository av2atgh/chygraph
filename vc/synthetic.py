"""The regime the bracket is about, planted: dense blocks (near-cliques of
5..40 vertices, internal density rho) joined by a sparse random graph of
inter-block links (mean external degree c per vertex).  Does the nested
DC-SBM recover the blocks, and does the complex-level cavity get the exact
cover where the node-level one does not?"""
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from chybp import ChyBP                     # noqa: E402
from compare import exact_vc, run_bp        # noqa: E402
from ensemble import leaf_removal           # noqa: E402
from sbm_sizes import fit                   # noqa: E402


def planted(n_blocks, rho, c, seed, smin=5, smax=40):
    rng = np.random.default_rng(seed)
    sizes = rng.integers(smin, smax + 1, size=n_blocks)
    blocks = np.repeat(np.arange(n_blocks), sizes)
    n = len(blocks)
    starts = np.concatenate([[0], np.cumsum(sizes)[:-1]])
    edges = []
    for b, (s0, s) in enumerate(zip(starts, sizes)):
        for i in range(s):
            for j in range(i):
                if rng.random() < rho:
                    edges.append((s0 + i, s0 + j))
    m_ext = int(c * n / 2)
    a = rng.integers(n, size=m_ext)
    b = rng.integers(n, size=m_ext)
    for x, y in zip(a, b):
        if blocks[x] != blocks[y]:
            edges.append((min(x, y), max(x, y)))
    e = np.unique(np.array(edges), axis=0)
    return n, e, blocks


def main():
    rows = []
    cases = [(float(sys.argv[1]), float(sys.argv[2]))] if len(sys.argv) > 2 else \
        [(rho, c) for rho in (1.0, 0.8) for c in (0.5, 1.0, 2.0)]
    for rho, c in cases:
        if True:
            n, e, planted_b = planted(150, rho, c, seed=1)
            t = time.time()
            res = exact_vc(n, e, 900)
            ex, exs = res.fun, time.time() - t
            lr, core = leaf_removal(n, e)
            sbm_b, levels, _ = fit(n, e, 0)
            # agreement of SBM blocks with planted: fraction of vertices whose SBM block is exactly a planted block
            same = 0
            for k in np.unique(sbm_b):
                vs = np.flatnonzero(sbm_b == k)
                pb = np.unique(planted_b[vs])
                if len(pb) == 1 and (planted_b == pb[0]).sum() == len(vs):
                    same += len(vs)
            print(f'rho={rho} c={c} n={n} m={len(e)} exact={ex:.0f} ({exs:.0f}s) leaf={lr} core={len(core)} '
                  f'SBM blocks={len(np.unique(sbm_b))} planted=150 exact-recovered={same / n:.2f}', flush=True)
            for label, bl in (('planted', planted_b), ('sbm', sbm_b)):
                for smax in (0, 64):
                    if smax == 0 and label == 'sbm':
                        continue
                    r = run_bp(n, e, bl, smax, iters=400, damping=0.5)
                    print(f"   {label:8s} smax={smax:3d} regions={r['regions']:5d} ext={r['ext_edges']:5d} "
                          f"iters={r['iters']:4d} resid={r['resid']:.1e} bethe={r['bethe']:.2f} "
                          f"cover={r['cover']} (exact {ex:.0f}) [{r['secs']:.0f}s]", flush=True)
                    rows.append(dict(rho=rho, c=c, n=n, exact=ex, leaf=lr, core=len(core), label=label, **r))
    import json
    tag = f'_{cases[0][0]}_{cases[0][1]}' if len(cases) == 1 else ''
    Path(f'results/synthetic{tag}.json').write_text(json.dumps(rows, indent=1))


if __name__ == '__main__':
    main()
