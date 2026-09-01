"""Generalised belief propagation on maximal-clique region graphs.

Section VII of the manuscript measures how far a chygraph is from treelike and
attributes the resulting deficit to overlap; the repair it names is generalised
belief propagation on the region graph.  ``statmech.gbp`` implements
it, and this probe asks how far it actually gets on the structures that
motivated the paper.

Method.  Sample a hyperbolic random graph small enough that ``ln Z`` can be
enumerated exactly, take its maximal cliques as complexes, close the family
under intersection, and compare four numbers: the exact ``ln Z``, the Bethe
counting of Eq. (13), the Mobius counting of Eq. (37) evaluated on isolated
regions, and generalised belief propagation run to a fixed point.  Each pair
inside a clique is an interaction counted once, however many cliques contain it.

What it finds, and the honest form of it: GBP is exact exactly when the graph is
chordal, since then the maximal cliques and their intersections are a junction
tree.  On non-chordal instances the parent-to-child iteration converges on some
and not others even at damping 0.999; where it converges the error is orders of
magnitude below either static counting, and where it does not the residual says
so and no number should be quoted.  That is what makes the ensemble-level lift
the open problem it is.

    python probe/gbp_cliques.py
"""

import json
import sys
from pathlib import Path

import networkx as nx
import numpy as np

sys.path.insert(0, str(Path.home() / 'av2atg' / 'computational_complexity' / 'code'))
from hrg import hrg_calibrated  # noqa: E402

from statmech.gbp import (GBP, exact_log_Z, ising_factors,  # noqa: E402
                                   static_log_Z)
from statmech.region import RegionGraph, overlap_profile  # noqa: E402

OUT = Path(__file__).parent / 'results' / 'gbp_cliques.json'
SIZES = ((14, 3.5), (18, 5.0), (20, 7.0))
SEEDS = range(1, 11)
COUPLINGS = (0.3, 0.8)
DAMPING = (0.5, 0.9, 0.97, 0.995, 0.999)
CONVERGED = 1e-9


def instance(n, kbar, seed, tau=2.5):
    src, dst, *_ = hrg_calibrated(n=n, kbar=kbar, tau=tau, rng=seed)
    G = nx.Graph()
    G.add_nodes_from(range(n))
    G.add_edges_from(zip(src.tolist(), dst.tolist()))
    G.remove_edges_from(nx.selfloop_edges(G))
    return G


def run_one(G, n, bJ):
    cx = [sorted(c) for c in nx.find_cliques(G)]
    prof = overlap_profile(cx)
    rg = RegionGraph(cx, max_rounds=6)
    if not rg.counting_is_valid():
        return None
    edges = sorted({tuple(sorted(e)) for e in G.edges()})
    f = ising_factors(edges, bJ)
    exact = exact_log_Z(f, range(n))
    best = None
    for d in DAMPING:
        g = GBP(rg, f, damping=d).run(8000)
        if best is None or g.residual < best[0].residual:
            best = (g, d)
        if g.residual < CONVERGED:
            break
    g, d = best
    return {
        'n': n, 'beta_J': bJ, 'n_cliques': len(cx),
        'shared_2plus': prof['shared_2plus'], 'chordal': bool(nx.is_chordal(G)),
        'exact': exact,
        'bethe': static_log_Z(rg.bethe_counting(), f) - exact,
        'kikuchi': static_log_Z(rg.counting, f) - exact,
        'gbp': g.log_Z() - exact, 'residual': g.residual, 'damping': d,
        'consistency': g.consistency() if np.isfinite(g.residual) else None,
    }


def main():
    rows = []
    for n, kbar in SIZES:
        for seed in SEEDS:
            G = instance(n, kbar, seed)
            if overlap_profile([sorted(c) for c in nx.find_cliques(G)])['treelike']:
                continue
            for bJ in COUPLINGS:
                r = run_one(G, n, bJ)
                if r is None:
                    continue
                r['seed'] = seed
                rows.append(r)
                print(f"n={n} seed={seed} bJ={bJ} chordal={r['chordal']} "
                      f"gbp={r['gbp']:+.2e} res={r['residual']:.0e} "
                      f"kik={r['kikuchi']:+.2e} bethe={r['bethe']:+.2e}",
                      flush=True)

    OUT.parent.mkdir(parents=True, exist_ok=True)
    json.dump(rows, OUT.open('w'), indent=1)
    summarise(rows)
    print('wrote', OUT)


def summarise(rows):
    def worst(sel, key):
        v = [abs(r[key]) for r in sel]
        return (min(v), max(v)) if v else (float('nan'),) * 2

    ch = [r for r in rows if r['chordal']]
    nc = [r for r in rows if not r['chordal']]
    ncc = [r for r in nc if r['residual'] < CONVERGED]
    print(f"\nchordal        {len(ch):>3} runs, GBP error "
          f"{worst(ch,'gbp')[1]:.1e} at worst")
    print(f"non-chordal    {len(nc):>3} runs, {len(ncc)} converged")
    if ncc:
        lo, hi = worst(ncc, 'gbp')
        klo, khi = worst(ncc, 'kikuchi')
        print(f"  converged:   GBP {lo:.1e}-{hi:.1e}, "
              f"static Kikuchi {klo:.1e}-{khi:.1e}")
    print(f"static Bethe   {worst(rows,'bethe')[0]:.1e}-"
          f"{worst(rows,'bethe')[1]:.1e} over all runs")


if __name__ == '__main__':
    main()
