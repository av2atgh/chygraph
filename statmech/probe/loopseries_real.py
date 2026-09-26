"""The double-count / loop split of Sec. 24.3 on Chapter 14's real neighbourhoods.

``gbp_real.py`` measured the chygraph recursion on 120 ego-networks of six
real networks (ten per network, 8 to 20 vertices, maximal cliques not
treelike), at two couplings, against exact enumeration.  This regenerates the
same neighbourhoods -- same cache, same seeds -- and evaluates on each

    recursion   ln Z_BP^dc - ln Z    (every clique given every bond inside it)
    double count   ln Z^dc - ln Z
    loop           ln Z_BP^dc - ln Z^dc
    assigned once  ln Z_BP^once - ln Z  (each bond in its largest clique)

with the factor graphs of ``statmech.loopseries``.  Results are cached in
``results/loopseries_real.json``; ``book/figs/loopseries.py`` tabulates them.

    python probe/loopseries_real.py
"""

import json
import sys
import time
from pathlib import Path

import networkx as nx
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'src'))
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'percolation' / 'src'))
import statmech.gbp  # noqa: E402,F401  (before the probe path, see cavity_clique)
sys.path.append(str(Path(__file__).resolve().parent))
from statmech.loopseries import BinaryFactorGraph  # noqa: E402
from gbp_real import NETWORKS, SEED, egonets, load  # noqa: E402

OUT = Path(__file__).parent / 'results' / 'loopseries_real.json'
COUPLINGS = (0.3, 0.8)


def split(cx, e, bJ):
    dc = BinaryFactorGraph.promoted(cx, e, bJ, 'all').bp(damping=0.5)
    once = BinaryFactorGraph.promoted(cx, e, bJ, 'once').bp(damping=0.5)
    z = BinaryFactorGraph.pairwise(e, bJ).exact_log_Z()
    zdc = dc.exact_log_Z()
    zb = dc.log_Z_bethe()
    return dict(recursion=zb - z, double_count=zdc - z, loop=zb - zdc,
                once=once.log_Z_bethe() - z,
                converged=bool(dc.converged() and once.converged()))


def main():
    t0 = time.time()
    rows = []
    for key, label, family in NETWORKS:
        G = load(key)
        rng = np.random.default_rng(SEED + [k for k, _, _ in NETWORKS].index(key))
        eg = egonets(G, rng)
        print(f'{label}: {len(eg)} ego-networks', flush=True)
        for v, H in eg:
            cx = [tuple(sorted(c)) for c in nx.find_cliques(H)]
            e = sorted(tuple(sorted(x)) for x in H.edges())
            for bJ in COUPLINGS:
                t = time.time()
                r = split(cx, e, bJ)
                r.update(network=label, family=family, ego=str(v), n=H.number_of_nodes(),
                         m=len(e), beta_J=bJ, n_cliques=len(cx))
                rows.append(r)
                print(f"  ego={str(v)[:12]:>12s} n={r['n']:>2d} bJ={bJ}: recursion {r['recursion']:+.3f} "
                      f"= dc {r['double_count']:+.3f} + loop {r['loop']:+.3f}; once {r['once']:+.3f} "
                      f"conv={r['converged']} ({time.time() - t:.0f}s)", flush=True)
    json.dump(rows, OUT.open('w'), indent=1)
    print(f'wrote {OUT} in {time.time() - t0:.0f} s')


if __name__ == '__main__':
    main()
