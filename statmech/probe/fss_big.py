"""Ch. 17's item (d): a fourth size for the two-dimensional finite-size
scaling of the acyclicity threshold, n = 320000 at the couplings that
bracket it.  Appends to ``results/acyclicity_fss.csv``."""
import csv
import sys
import time
from pathlib import Path

import numpy as np
import networkx as nx

sys.path.append(str(Path(__file__).resolve().parent))
from acyclicity import OUT, _done, _existing, measure, rgg_torus  # noqa: E402

path = OUT / 'acyclicity_fss.csv'
rows = _existing(path)
for kbar in (6.0, 7.0, 8.0, 9.0, 5.0, 10.0):
    for seed in (1, 2):
        n = 320000
        if _done(rows, param=2, kbar=kbar, seed=seed, n=n):
            continue
        rng = np.random.default_rng(seed)
        t0 = time.time()
        G = rgg_torus(n, 2, kbar, rng)
        row = dict(ensemble='rgg', param=2, kbar=kbar, seed=seed)
        row['kbar_meas'] = 2 * G.number_of_edges() / n
        row['gc'] = max(len(c) for c in nx.connected_components(G)) / n
        row.update(measure(G))
        row['sec'] = round(time.time() - t0, 1)
        rows.append(row)
        print(' '.join(f'{k}={v:.4g}' if isinstance(v, float) else f'{k}={v}' for k, v in row.items()), flush=True)
        with open(path, 'w', newline='') as f:
            w = csv.DictWriter(f, fieldnames=list(rows[0]))
            w.writeheader()
            w.writerows(rows)
