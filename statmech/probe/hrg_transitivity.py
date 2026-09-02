"""Transitivity of the graphs behind prediction4.csv, so Fig. 11.4 can be
replotted against structure instead of against density.

Figure 11.4 puts the leaf-removal core of hyperbolic random graphs against mean
degree, log-log, with the chygraph prediction and the degree-matched control.
Section 11.8's Fig. 11.6 does the same for the Vazquez-Weigt ensemble, and the
two are worth putting in the same coordinate: there the knob is r and it makes
triangles as a side effect, here the knob is density and the clustering it
produces is the geometry's.

This regenerates exactly the graphs of probe/prediction4.py -- same taus, same
target mean degrees, same seeds, same calibration -- and measures the
average local clustering of the HRG and of its degree-matched configuration
model, with the transitivity alongside -- on a heavy tail the two disagree by
orders of magnitude, and transitivity turns out not to be usable here at all,
being set by whichever few hubs the seed drew rather than by the density being
swept.  No
clique enumeration and no chygraph, which is where prediction4's time went, so
the run is minutes rather than an hour.  The core fractions are not recomputed:
they are joined from prediction4.csv on (tau, kbar_target, seed).

Writes probe/results/hrg_transitivity.csv.
"""

import csv
import sys
import time
from pathlib import Path

import networkx as nx
import numpy as np

CC = Path.home() / 'av2atg' / 'computational_complexity' / 'code'
sys.path.insert(0, str(CC))

from hrg import erased_configuration_model, hrg_calibrated            # noqa: E402

TAUS = (2.5, 2.9)                 # the two panels of Fig. 11.4
KBARS = (0.05, 0.1, 0.2, 0.4, 0.8, 1.5, 3.0, 6.0)
N = 200_000
SEEDS = (1, 2)
OUT = Path(__file__).resolve().parent / 'results' / 'hrg_transitivity.csv'


def structure(n, src, dst):
    """Transitivity and average local clustering, the latter two ways.

    Transitivity is 3 x triangles / triples over the whole graph, so a hub
    contributes d^2 triples and dilutes it; on a tau = 2.5 tail a handful of
    hubs set the denominator and T measures which hubs the seed drew, not how
    clustered the graph is.  Average local clustering weights every vertex
    equally instead.

    Which vertices, though.  These graphs are mostly isolated vertices at the
    sparse end, and a vertex of degree < 2 cannot be in a triangle, so counting
    it as c = 0 makes the average track the *density* -- the fraction of
    vertices that have two neighbours at all -- rather than the clustering.

      C_all  the mean over all n vertices, networkx's default and the
             convention of Chapter 3's table, where isolated vertices are rare
      C      the mean over vertices of degree >= 2, which is the quantity that
             says how clustered the part of the graph that can be clustered is

    C is what Fig. 11.7 plots, and the two are reported together because they
    differ by three orders of magnitude at kbar = 0.05.
    """
    G = nx.Graph()
    G.add_nodes_from(range(n))
    G.add_edges_from(zip(src.tolist(), dst.tolist()))
    c = nx.clustering(G)
    deg = dict(G.degree())
    eligible = [c[v] for v in range(n) if deg[v] >= 2]
    return (float(nx.transitivity(G)),
            float(np.mean(list(c.values()))),
            float(np.mean(eligible)) if eligible else 0.0,
            len(eligible) / n)


def run():
    rows = []
    for tau in TAUS:
        for kbar in KBARS:
            for seed in SEEDS:
                t0 = time.time()
                src, dst, *_ = hrg_calibrated(N, tau=tau, kbar=kbar, rng=seed,
                                              tol=0.02, max_iter=25)
                deg = np.bincount(np.concatenate((src, dst)), minlength=N)
                kb = 2.0 * src.size / N
                ht, hca, hc, hf = structure(N, src, dst)
                cs, cd = erased_configuration_model(deg, rng=seed)
                ct, cca, cc, cf = structure(N, cs, cd)
                rows.append(dict(tau=tau, kbar_target=kbar, kbar=kb, seed=seed,
                                 hrg_trans=ht, hrg_C_all=hca, hrg_C=hc,
                                 hrg_frac2=hf,
                                 cfg_trans=ct, cfg_C_all=cca, cfg_C=cc,
                                 cfg_frac2=cf,
                                 secs=round(time.time() - t0, 1)))
                print(f'tau={tau} kbar={kb:.3f} seed={seed}  '
                      f'HRG T={ht:.4f} C={hc:.4f} C_all={hca:.4f} '
                      f'f2={hf:.3f}   cfg C={cc:.5f}  '
                      f'[{rows[-1]["secs"]}s]', flush=True)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    with OUT.open('w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    print(f'\nwrote {OUT} ({len(rows)} rows)')


if __name__ == '__main__':
    run()
