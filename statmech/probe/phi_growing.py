"""Does phi on the interactomes vanish in a growing ensemble with the same
degree tail?  (Ch. 17's item (e).)

For each real network of Ch. 17's table, the degree sequence is replicated
k = 1, 2, 4, 8 times and a configuration-model graph drawn from it; on its
giant component the fraction of vertices on a short chordless cycle,
phi_short, and the giant of that set are measured as in ``acyclicity.py``.
k = 1 is the chapter's rewired control.  If phi falls with k the loops the
degree sequence makes are a finite-size effect; if it does not they are a
property of the tail.  Results in ``results/phi_growing.csv``.
"""

import csv
import sys
import time
from pathlib import Path

import networkx as nx
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'src'))
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'percolation' / 'src'))
import statmech.gbp  # noqa: E402,F401
sys.path.append(str(Path(__file__).resolve().parents[2] / 'book' / 'figs'))
sys.path.append(str(Path(__file__).resolve().parent))
from acyclicity import on_short_chordless_cycle  # noqa: E402
from real_acyclicity import NET14, NETS3, load14, load_stem  # noqa: E402

OUT = Path(__file__).parent / 'results' / 'phi_growing.csv'
REPS = (1, 2, 4, 8)
SEEDS = (0, 1, 2)


def main():
    rows = []
    jobs = [(lbl, stem, load_stem) for (lbl, stem, fam) in NETS3]
    jobs += [(lbl, key, load14) for (key, lbl, fam) in NET14]
    for lbl, key, loader in jobs:
        G = loader(key)
        G = G.subgraph(max(nx.connected_components(G), key=len))
        deg = [d for _, d in G.degree()]
        n0 = G.number_of_nodes()
        real = on_short_chordless_cycle(nx.convert_node_labels_to_integers(G))
        rows.append(dict(network=lbl, k=0, n=n0, phi=len(real) / n0, giant=float('nan'), sec=0.0))
        for k in REPS:
            t0 = time.time()
            vals = []
            for seed in SEEDS:
                H = nx.Graph(nx.configuration_model(deg * k, seed=seed))
                H.remove_edges_from(nx.selfloop_edges(H))
                H = H.subgraph(max(nx.connected_components(H), key=len)).copy()
                bad = on_short_chordless_cycle(H)
                comp = [len(c) for c in nx.connected_components(H.subgraph(bad))] if bad else []
                vals.append((len(bad) / H.number_of_nodes(), (max(comp) / H.number_of_nodes()) if comp else 0.0))
            vals = np.array(vals)
            rows.append(dict(network=lbl, k=k, n=n0 * k, phi=vals[:, 0].mean(), giant=vals[:, 1].mean(),
                             sec=round(time.time() - t0, 1)))
            print(f"{lbl:22s} k={k} n={n0 * k:6d} phi={vals[:, 0].mean():.4f} giant={vals[:, 1].mean():.4f} "
                  f"(real phi {len(real) / n0:.4f}) {time.time() - t0:.0f}s", flush=True)
        with open(OUT, 'w', newline='') as f:
            w = csv.DictWriter(f, fieldnames=list(rows[0]))
            w.writeheader()
            w.writerows(rows)


if __name__ == '__main__':
    main()
