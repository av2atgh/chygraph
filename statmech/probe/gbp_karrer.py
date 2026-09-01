"""Generalised belief propagation on Karrer-Newman placed-motif graphs.

The third of the chapter's three ensembles, run under exactly the protocol of
``gbp_cliques.py`` (hyperbolic random graphs) and ``gbp_real.py`` (real
networks), so that the three are comparable line for line: exact ``ln Z`` by
enumeration, the Bethe counting, the Mobius counting on isolated regions, and
GBP to a fixed point over the same damping ladder.

The ensemble is \\citet{karrer2010}'s: single edges by stub matching, each motif
by matching its corners in tuples.  Two placed triangles may share a vertex
freely; they share an *edge* only when their remaining corners happen to
coincide, which is a coincidence of the matching rather than a feature of the
model.  So the ensemble is treelike above the level of the motif by
construction, and the question this probe asks is what that buys the repair.

The sizes and mean degrees match ``gbp_cliques.py`` --- $n=14,18,20$ at
$\\ave{k}=3.5,5,7$ --- with $\\ave{k}\\simeq s+2t$ for $s$ single-edge stubs and
$t$ triangles per vertex.  Because most samples come out treelike and are
filtered, seeds are drawn until ``PER_SIZE`` non-treelike instances are found;
how many had to be drawn is itself reported, being the ensemble's own answer to
how often placed motifs overlap at all.

    python probe/gbp_karrer.py
"""

import json
import sys
from pathlib import Path

import networkx as nx
import numpy as np

sys.path.insert(0, str(Path.home() / 'av2atg' / 'chygraph' / 'src'))
sys.path.insert(0, str(Path.home() / 'av2atg' / 'chygraph_statmech'
                      / 'book' / 'figs'))

from merge import karrer_graph  # noqa: E402  the book's one definition of it

from chygraph_statmech.gbp import (GBP, exact_log_Z, ising_factors,  # noqa: E402
                                   static_log_Z)
from chygraph_statmech.region import RegionGraph, overlap_profile  # noqa: E402

OUT = Path(__file__).parent / 'results' / 'gbp_karrer.json'

# (n, single-edge stubs per vertex, triangles per vertex) -> <k> = s + 2t,
# matching gbp_cliques.py's (14, 3.5), (18, 5.0), (20, 7.0)
SIZES = ((14, 1.5, 1.0), (18, 1.0, 2.0), (20, 1.0, 3.0))
PER_SIZE = 10          # non-treelike instances kept per size
MAX_SEEDS = 400        # give up after this many draws
COUPLINGS = (0.3, 0.8)
DAMPING = (0.5, 0.9, 0.97, 0.995, 0.999)
CONVERGED = 1e-9


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
        'n': n, 'm': len(edges), 'beta_J': bJ,
        'n_cliques': len(cx), 'n_regions': len(rg.regions),
        'shared_2plus': prof['shared_2plus'],
        'chordal': bool(nx.is_chordal(G)),
        'exact': exact,
        'bethe': static_log_Z(rg.bethe_counting(), f) - exact,
        'kikuchi': static_log_Z(rg.counting, f) - exact,
        'gbp': g.log_Z() - exact, 'residual': g.residual, 'damping': d,
        'consistency': g.consistency() if np.isfinite(g.residual) else None,
        'mean_abs_m': (float(np.mean(np.abs(list(g.magnetisation().values()))))
                       if np.isfinite(g.residual) else None),
    }


def main():
    rows, census = [], []
    for n, s_mean, t in SIZES:
        kept = drawn = 0
        for seed in range(1, MAX_SEEDS + 1):
            if kept >= PER_SIZE:
                break
            drawn = seed
            G, _ = karrer_graph(n, s_mean, {3: t}, seed)
            G.remove_edges_from(nx.selfloop_edges(G))
            cx = [sorted(c) for c in nx.find_cliques(G)]
            if overlap_profile(cx)['treelike']:
                continue
            kept += 1
            for bJ in COUPLINGS:
                r = run_one(G, n, bJ)
                if r is None:
                    continue
                r.update(seed=seed, s_mean=s_mean, triangles=t,
                         kbar=2 * G.number_of_edges() / n)
                rows.append(r)
                print(f"n={n} seed={seed} bJ={bJ} chordal={r['chordal']:d} "
                      f"gbp={r['gbp']:+.2e} res={r['residual']:.0e} "
                      f"kik={r['kikuchi']:+.2e} bethe={r['bethe']:+.2e}",
                      flush=True)
        census.append({'n': n, 's_mean': s_mean, 'triangles': t,
                       'kept': kept, 'drawn': drawn})
        print(f'  n={n}: {kept} non-treelike out of {drawn} drawn', flush=True)

    OUT.parent.mkdir(parents=True, exist_ok=True)
    json.dump({'runs': rows, 'census': census}, OUT.open('w'), indent=1)
    summarise(rows, census)
    print('wrote', OUT)


def summarise(rows, census):
    tot_kept = sum(c['kept'] for c in census)
    tot_drawn = sum(c['drawn'] for c in census)
    print(f'\n{tot_kept} non-treelike instances out of {tot_drawn} graphs '
          f'drawn ({100 * tot_kept / tot_drawn:.0f} per cent)')
    ch = [r for r in rows if r['chordal']]
    nc = [r for r in rows if not r['chordal']]
    print(f'chordal      {len(ch):>3} runs, worst GBP error '
          f'{max((abs(r["gbp"]) for r in ch), default=float("nan")):.1e}')
    snd = [r for r in nc if r['residual'] < CONVERGED
           and r['consistency'] is not None and r['consistency'] < 1e-6]
    print(f'non-chordal  {len(nc):>3} runs, {len(snd)} with a sound fixed point')
    if snd:
        g = [abs(r['gbp']) for r in snd]
        k = [abs(r['kikuchi']) for r in snd]
        print(f'  sound: GBP {min(g):.1e} to {max(g):.1e}, '
              f'Mobius {min(k):.1e} to {max(k):.1e}')
        worse = [r for r in snd if abs(r['gbp']) >= abs(r['kikuchi'])]
        print(f'  no better than static: {len(worse)}')


if __name__ == '__main__':
    main()
