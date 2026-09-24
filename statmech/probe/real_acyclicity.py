"""How alpha-acyclic are real networks?  GYO residue on the sixteen networks of
Chs. 3 and 14, with the merge closure of Ch. 16 alongside.

    PYTHONPATH=../src:../../percolation/src python probe/real_acyclicity.py

Writes probe/results/real_acyclicity.csv.
"""

import csv
import sys
import time
from pathlib import Path

import networkx as nx
import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
from acyclicity import measure, on_short_chordless_cycle  # noqa: E402
from gbp_real import load as load14  # noqa: E402
from gbp_real import NETWORKS as NET14  # noqa: E402

OUT = Path(__file__).parent / 'results' / 'real_acyclicity.csv'
NREWIRE = 5
DATA = (Path.home() / 'av2atg' / 'LocalNetworkGrowth' / 'figs' / 'data' / 'netzschleuder')

# Chapter 3's ten, copied from real_core.py (which imports modules that no longer exist).
NETS3 = [
    ('Collins yeast', 'collins_yeast__collins_yeast', 'grouped'),
    ('Yeast (Y2H)', 'interactome_yeast__interactome_yeast', 'dyadic'),
    ('Human (Vidal)', 'interactome_vidal__interactome_vidal', 'dyadic'),
    ('Human (Stelzl)', 'interactome_stelzl__interactome_stelzl', 'dyadic'),
    ('Human (Figeys)', 'interactome_figeys__interactome_figeys', 'grouped'),
    ('PDZ domains', 'interactome_pdz__interactome_pdz', 'dyadic'),
    ('C. elegans WI8', 'celegans_interactomes__WI8', 'dyadic'),
    ('C. elegans 2007', 'celegans_interactomes__wi2007', 'dyadic'),
    ('Power grid', 'power__power', 'dyadic'),
    ('Euroroad', 'euroroad__euroroad', 'dyadic'),
]


def load_stem(stem):
    import io as _io, zipfile
    with zipfile.ZipFile(DATA / f'{stem}.csv.zip') as z:
        raw = z.read('edges.csv').decode()
    G = nx.Graph()
    for line in _io.StringIO(raw):
        if line.startswith('#') or not line.strip():
            continue
        u, v = line.split(',')[:2]
        u, v = u.strip(), v.strip()
        if u != v:
            G.add_edge(u, v)
    return G



def main():
    rows = []
    jobs = [(lbl, stem, fam, load_stem) for (lbl, stem, fam) in NETS3]
    jobs += [(lbl, key, fam, load14) for (key, lbl, fam) in NET14]
    for lbl, key, fam, loader in jobs:
        t0 = time.time()
        G = loader(key)
        G = nx.convert_node_labels_to_integers(G.subgraph(max(nx.connected_components(G), key=len)))
        row = dict(network=lbl, key=key, family=fam, kbar=2 * G.number_of_edges() / G.number_of_nodes(),
                   transitivity=nx.transitivity(G))
        row.update(measure(G))
        # degree-matched control: the same degree sequence, rewired.  Short
        # chordless cycles that survive the rewiring are a hub effect, present
        # in any graph with this degree sequence at this n; the ones that do
        # not are structure.
        ctrl = []
        for seed in range(NREWIRE):
            H = nx.configuration_model([d for _, d in G.degree()], seed=seed)
            H = nx.Graph(H)
            H.remove_edges_from(nx.selfloop_edges(H))
            H = H.subgraph(max(nx.connected_components(H), key=len)).copy()
            bad = on_short_chordless_cycle(H)
            comp = [len(c) for c in nx.connected_components(H.subgraph(bad))] if bad else []
            ctrl.append((len(bad) / H.number_of_nodes(), (max(comp) / H.number_of_nodes()) if comp else 0.0,
                         nx.transitivity(H)))
        ctrl = np.array(ctrl)
        row['phi_short_ctrl'], row['giant_short_ctrl'], row['transitivity_ctrl'] = ctrl.mean(0)
        row['sec'] = round(time.time() - t0, 1)
        rows.append(row)
        print(' '.join(f'{k}={v:.4g}' if isinstance(v, float) else f'{k}={v}' for k, v in row.items()), flush=True)
        with open(OUT, 'w', newline='') as f:
            w = csv.DictWriter(f, fieldnames=list(rows[0]))
            w.writeheader()
            w.writerows(rows)


if __name__ == '__main__':
    main()
