"""Prediction 4: does a chygraph of the HRG's cliques reproduce its core?

`~/av2atg/computational_complexity` measured that hyperbolic random graphs keep
an extensive leaf-removal core at *every* mean degree, with no transition, and
that at small `kbar` the core is a pure power law

    core ~ kbar^beta,   beta = 1.22 (tau=2.1), 1.50 (2.5), 1.63 (2.9),

down to kbar = 0.05, while the degree-matched configuration model is core-free
until kbar ~ 4.5-10.  Sec. VI writes core percolation as a chygraph fixed point;
Sec. VIII showed the clique-ensemble chygraph accounts for 65-83% of the core at
kbar ~ 2-4.  The open question is the *shape*: a chygraph of independent
complexes predicts a core of `1 - Phi(0)`, the chance of belonging to any
complex, which is not obviously a power law in kbar at all.

For each (tau, kbar) this measures three things:

  measured   pure leaf removal on the HRG itself
  chygraph   Sec. VI on the maximal-clique ensemble measured from that graph
  control    the same, on the degree-matched configuration model

and fits `core ~ kbar^beta` to each over the small-kbar range.

Writes probe/results/prediction4.csv.
"""

import csv
import sys
import time
from collections import Counter
from pathlib import Path

import networkx as nx
import numpy as np

CC = Path.home() / 'av2atg' / 'computational_complexity' / 'code'
sys.path.insert(0, str(CC))
import leafremoval as lr  # noqa: E402
from hrg import erased_configuration_model, hrg_calibrated, to_csr  # noqa: E402

from chygraph_statmech import Chygraph  # noqa: E402

TAUS = (2.1, 2.5, 2.9)
KBARS = (0.05, 0.1, 0.2, 0.4, 0.8, 1.5, 3.0, 6.0)
N = 200_000
SEEDS = (1, 2)
FIT_BELOW = 1.0          # fit the power law on kbar <= this


def clique_chygraph(n, src, dst):
    """Maximal cliques of the graph as the complex ensemble."""
    G = nx.Graph()
    G.add_nodes_from(range(n))
    G.add_edges_from(zip(src.tolist(), dst.tolist()))
    cliques = [c for c in nx.find_cliques(G) if len(c) >= 2]
    if not cliques:
        return None
    cards = sorted({len(c) for c in cliques})
    idx = {c: i for i, c in enumerate(cards)}
    K = np.zeros((n, len(cards)))
    for a in cliques:
        j = idx[len(a)]
        for v in a:
            K[v, j] += 1
    keep = [i for i in range(len(cards)) if K[:, i].mean() > 0]
    return Chygraph.from_samples([cards[i] for i in keep], K[:, keep])


def measure(n, src, dst):
    g = clique_chygraph(n, src, dst)
    pred = g.core_from_samples().core_fraction() if g is not None else 0.0
    meas = lr.core(*to_csr(n, src, dst))[0] / n
    return float(pred), float(meas)


def fit(kb, y):
    kb, y = np.asarray(kb, float), np.asarray(y, float)
    ok = (kb > 0) & (y > 1e-9) & (kb <= FIT_BELOW)
    if ok.sum() < 3:
        return float('nan')
    return float(np.polyfit(np.log(kb[ok]), np.log(y[ok]), 1)[0])


def run(out):
    rows = []
    for tau in TAUS:
        for kbar in KBARS:
            for seed in SEEDS:
                t0 = time.time()
                src, dst, r, th, R, C = hrg_calibrated(
                    N, tau=tau, kbar=kbar, rng=seed, tol=0.02, max_iter=25)
                deg = np.bincount(np.concatenate((src, dst)), minlength=N)
                kb = 2.0 * src.size / N
                hp, hm = measure(N, src, dst)
                cs, cd = erased_configuration_model(deg, rng=seed)
                cp, cm = measure(N, cs, cd)
                rows.append(dict(tau=tau, kbar_target=kbar, kbar=kb, seed=seed,
                                 hrg_chygraph=hp, hrg_measured=hm,
                                 cfg_chygraph=cp, cfg_measured=cm,
                                 secs=round(time.time() - t0, 1)))
                print(f"tau={tau} kbar={kb:.3f} seed={seed}  "
                      f"HRG chy={hp:.5f} meas={hm:.5f}   "
                      f"cfg chy={cp:.5f} meas={cm:.5f}", flush=True)
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open('w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)

    print(f"\npower-law exponents on kbar <= {FIT_BELOW}: core ~ kbar^beta")
    print("computational_complexity measured 1.22 / 1.50 / 1.63 "
          "at tau = 2.1 / 2.5 / 2.9")
    print(f"{'tau':>5}{'beta measured':>15}{'beta chygraph':>15}{'beta control':>14}")
    for tau in TAUS:
        sel = [x for x in rows if x['tau'] == tau]
        ks = sorted({x['kbar_target'] for x in sel})
        avg = lambda key: [float(np.mean([x[key] for x in sel
                                          if x['kbar_target'] == k])) for k in ks]
        kb = [float(np.mean([x['kbar'] for x in sel
                             if x['kbar_target'] == k])) for k in ks]
        print(f"{tau:>5.1f}{fit(kb, avg('hrg_measured')):>15.3f}"
              f"{fit(kb, avg('hrg_chygraph')):>15.3f}"
              f"{fit(kb, avg('cfg_measured')):>14.3f}")
    print(f"\nwrote {out} ({len(rows)} rows)")


if __name__ == '__main__':
    run(Path(__file__).parent / 'results' / 'prediction4.csv')
