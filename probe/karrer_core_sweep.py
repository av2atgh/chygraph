"""Does the merged chygraph's core appear at a threshold in the stub density?

Sections 16.2 and 16.3 report the core as an outcome -- so many meta-complexes,
so much core -- which cannot say what produces it.  The Karrer-Newman ensemble
has actual control parameters: n, the single edges per vertex s, and the
triangles per vertex t.  The mechanism of Sec. 16.2 says the cycles are carried
by the single-edge stubs, which form a sparse random graph threading several
complexes at single atoms, so s should be the knob.

This sweeps it.  For each (n, s) the merge closure is run on the maximal
cliques, the atom-complex incidence structure is built, and its 2-core is taken:
the atoms propagation cannot strip, which Sec. 16.2 finds is zero exactly when
the recursion is exact.  Averaging over realisations gives an order parameter,
and increasing n says whether it sharpens into a transition or stays a crossover.

    python probe/karrer_core_sweep.py
"""

import json
import sys
from pathlib import Path

import networkx as nx
import numpy as np

sys.path.insert(0, str(Path.home() / 'av2atg' / 'chygraph' / 'src'))
sys.path.insert(0, str(Path.home() / 'av2atg' / 'chygraph_statmech'
                      / 'book' / 'figs'))

from merge import karrer_graph, merge_closure  # noqa: E402

OUT = Path(__file__).parent / 'results' / 'karrer_core_sweep.json'
SIZES = (200, 400, 800, 1600)
SEEDS = 24

# Three lines through the (s, t) plane, all crossing the predicted threshold.
# Following an atom into a complex of cardinality c and out through its c-1
# other atoms, the incidence branching is s.1 + t.2, so the 2-core should
# appear at b = s + 2t = 1 whatever mixture of stubs and triangles produces it.
# At t = 0 the complexes are bare edges and the statement reduces to the
# Erdos-Renyi 2-core at mean degree one, which is the anchor.
FAMILIES = (('edges only, $t=0$', lambda b: (b, 0.0)),
            ('triangles only, $s=0$', lambda b: (0.0, b / 2)),
            ('equal mixture', lambda b: (b / 2, b / 4)))
B_GRID = [round(0.1 * i, 2) for i in range(0, 26)]


def core_fraction(n, s, t, seed):
    G, _ = karrer_graph(n, s, {3: t}, seed)
    G.remove_edges_from(nx.selfloop_edges(G))
    merged, _ = merge_closure([frozenset(c) for c in nx.find_cliques(G)])
    B = nx.Graph()
    B.add_nodes_from(('v', v) for c in merged for v in c)
    for i, c in enumerate(merged):
        for v in c:
            B.add_edge(('c', i), ('v', v))
    K = nx.k_core(B, k=2)
    return sum(1 for x in K if x[0] == 'v') / n


def main():
    rows = []
    for label, fam in FAMILIES:
        for n in SIZES:
            for bb in B_GRID:
                sm, tt = fam(bb)
                f = [core_fraction(n, sm, tt, sd) for sd in range(1, SEEDS + 1)]
                rows.append({'family': label, 'n': n, 'b': bb, 's': sm, 't': tt,
                             'core_mean': float(np.mean(f)),
                             'core_sem': float(np.std(f) / np.sqrt(len(f)))})
            print(f'  {label:<22} n={n:<5} done', flush=True)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    json.dump(rows, OUT.open('w'), indent=1)
    print('\nwrote', OUT)


if __name__ == '__main__':
    main()
