"""Probe for TODO item 1: can an HRG be given a chygraph complex ensemble?

A chygraph maps a clustered graph by making its dense motifs into *complexes*.
The natural choice for a hyperbolic random graph is the maximal cliques.  For
the threshold tensor to exist at all, the excess cardinality seen from a member
must have a finite mean, and

    <sbar> = <c^2>/<c> - 1

so the maximal-clique size distribution needs a **finite second moment**.  The
critical amplitude needs the third.  If those grow with ``n`` instead of
converging, there is no ``n``-independent complex ensemble, prediction 4's
clique route is dead, and TODO item 1 is answered without building anything.

The control is the degree-matched erased configuration model: same ``P(k)``,
same assortativity up to erasure, no geometry.  If the HRG's moments grow and
the control's do not, the growth is clustering, which is the whole question.

Writes probe/results/clique_moments.csv.
"""

import sys
import time
from collections import Counter
from pathlib import Path

import networkx as nx
import numpy as np

HRG_CODE = Path.home() / 'av2atg' / 'computational_complexity' / 'code'
sys.path.insert(0, str(HRG_CODE))
from hrg import erased_configuration_model, hrg_calibrated  # noqa: E402

TAUS = (2.1, 2.5, 2.9)
KBARS = (2.0, 4.0, 8.0)
NS = (1_000, 3_000, 10_000, 30_000, 100_000, 300_000)
SEEDS = (1, 2, 3)


def clique_moments(n, edges):
    """Moments of the maximal-clique size distribution, plus the largest."""
    G = nx.Graph()
    G.add_nodes_from(range(n))
    G.add_edges_from(zip(edges[0].tolist(), edges[1].tolist()))
    cnt = Counter(len(c) for c in nx.find_cliques(G))
    c = np.array(sorted(cnt), dtype=float)
    w = np.array([cnt[int(i)] for i in c], dtype=float)
    w /= w.sum()
    return {
        'n_cliques': int(sum(cnt.values())),
        'c_max': int(c.max()),
        'm1': float((w * c).sum()),
        'm2': float((w * c**2).sum()),
        'm3': float((w * c**3).sum()),
        # what the chygraph actually needs: excess cardinality from a member
        'sbar': float((w * c**2).sum() / (w * c).sum() - 1.0),
        'degeneracy': int(max(nx.core_number(G).values())),
    }


def run(out):
    rows = []
    for tau in TAUS:
        for kbar in KBARS:
            for n in NS:
                for seed in SEEDS:
                    t0 = time.time()
                    src, dst, r, th, R, C = hrg_calibrated(
                        n, tau=tau, kbar=kbar, rng=seed, tol=0.01, max_iter=25)
                    deg = np.bincount(np.concatenate((src, dst)), minlength=n)
                    for family, edges in (('hrg', (src, dst)),
                                          ('config', erased_configuration_model(
                                              deg, rng=seed))):
                        m = clique_moments(n, edges)
                        m.update(family=family, tau=tau, kbar_target=kbar, n=n,
                                 seed=seed,
                                 kbar=2.0 * edges[0].size / n,
                                 secs=round(time.time() - t0, 2))
                        rows.append(m)
                        print(f"{family:>6} tau={tau} kbar={m['kbar']:.2f} "
                              f"n={n:>7} seed={seed} c_max={m['c_max']:>4} "
                              f"m2={m['m2']:>9.3f} sbar={m['sbar']:>8.3f}",
                              flush=True)
    import csv
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open('w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    print(f"\nwrote {out} ({len(rows)} rows)")


if __name__ == '__main__':
    run(Path(__file__).parent / 'results' / 'clique_moments.csv')
