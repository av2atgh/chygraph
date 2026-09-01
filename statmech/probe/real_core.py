"""Leaf removal on real networks, against the chygraph of their cliques.

Chapter 11 tests Eq. (11.4) on hyperbolic random graphs: an ensemble built to be
clustered, where the degree-matched control has no core and the real object has
one everywhere.  That test is clean because the ensemble is synthetic -- the
clustering can be turned, the control is exact, and n can be taken to 2e5.

This runs the same three-way comparison on networks nobody designed:

  measured   pure leaf removal on the graph itself
  chygraph   Eq. (11.4) on the maximal-clique ensemble measured from that same
             graph, with nothing fitted
  control    leaf removal on the degree-matched configuration model, which keeps
             P(k) and destroys the clustering

The six networks are Chapter 14's, in its two families -- group-generated, where
an edge exists because the data recorded a group and every group is a clique of
its members, and dyadic, where edges were recorded one pair at a time.  Chapter
14 uses induced neighbourhoods of at most twenty vertices because it enumerates
ln Z; leaf removal is linear, so whole networks are used here.

Caching: probe/results/real_core.json.
"""

import json
import sys
from pathlib import Path

import networkx as nx
import numpy as np

CC = Path.home() / 'av2atg' / 'computational_complexity' / 'code'
sys.path.insert(0, str(CC))
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'percolation' / 'src'))
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'statmech' / 'src'))
sys.path.insert(0, str(Path(__file__).resolve().parent))

import leafremoval as lr                                        # noqa: E402
from hrg import erased_configuration_model, to_csr              # noqa: E402
from statmech import Chygraph                                   # noqa: E402
from gbp_real import NETWORKS as NETS14, load                   # noqa: E402

# Chapter 3's ten, from the same netzschleuder cache.  They matter because they
# reach below the leaf-removal threshold, where the control has no core and the
# comparison can say something; Chapter 14's six are all far above it.
DATA = (Path.home() / 'av2atg' / 'LocalNetworkGrowth' / 'figs' / 'data'
        / 'netzschleuder')
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

OUT = Path(__file__).resolve().parent / 'results' / 'real_core.json'
CONTROL_SEEDS = 20


def as_arrays(G):
    """(n, src, dst) over 0..n-1 for the largest connected component."""
    G = G.subgraph(max(nx.connected_components(G), key=len))
    idx = {v: i for i, v in enumerate(G.nodes())}
    src = np.array([idx[u] for u, v in G.edges()], dtype=np.int64)
    dst = np.array([idx[v] for u, v in G.edges()], dtype=np.int64)
    return len(idx), src, dst, G


def clique_chygraph(n, G):
    """Maximal cliques of the graph as the complex ensemble.  Nothing fitted."""
    cliques = [c for c in nx.find_cliques(G) if len(c) >= 2]
    if not cliques:
        return None
    cards = sorted({len(c) for c in cliques})
    idx = {v: i for i, v in enumerate(G.nodes())}
    col = {c: j for j, c in enumerate(cards)}
    K = np.zeros((n, len(cards)))
    for a in cliques:
        j = col[len(a)]
        for v in a:
            K[idx[v], j] += 1
    keep = [i for i in range(len(cards)) if K[:, i].mean() > 0]
    return Chygraph.from_samples([cards[i] for i in keep], K[:, keep])


def measured_core(n, src, dst):
    return float(lr.core(*to_csr(n, src, dst))[0]) / n


def run():
    todo = ([(lab, fam, 'ch14', lambda k=k: load(k)) for k, lab, fam in NETS14]
            + [(lab, fam, 'ch3', lambda st=st: load_stem(st)) for lab, st, fam in NETS3])
    rows = []
    for label, family, source, getter in todo:
        n, src, dst, G = as_arrays(getter())
        meas = measured_core(n, src, dst)

        g = clique_chygraph(n, G)
        pred = float(g.core_from_samples().core_fraction()) if g else 0.0

        deg = np.array([d for _, d in G.degree()], dtype=np.int64)
        ctrl = []
        for s in range(CONTROL_SEEDS):
            cs, cd = erased_configuration_model(deg, rng=s)
            ctrl.append(measured_core(n, cs, cd))
        ctrl = np.array(ctrl)

        cliques = [c for c in nx.find_cliques(G) if len(c) >= 2]
        sizes = np.array([len(c) for c in cliques], float)
        sbar = float((sizes ** 2).sum() / sizes.sum() - 1)

        row = dict(label=label, family=family, source=source, n=n,
                   m=int(G.number_of_edges()),
                   kbar=2 * G.number_of_edges() / n,
                   C=float(nx.average_clustering(G)),
                   sbar=sbar,
                   measured=meas, chygraph=pred,
                   control=float(ctrl.mean()), control_sd=float(ctrl.std()),
                   recovered=pred / meas if meas > 0 else float('nan'))
        rows.append(row)
        print(f"  {label:20s} n={n:4d} <k>={row['kbar']:5.2f}  "
              f"measured={meas:.4f}  chygraph={pred:.4f} ({row['recovered']*100:5.1f}%)  "
              f"control={ctrl.mean():.4f}({ctrl.std():.4f})", flush=True)

    OUT.parent.mkdir(exist_ok=True)
    OUT.write_text(json.dumps(rows, indent=1))
    print(f"\n  wrote {OUT}")
    return rows


if __name__ == '__main__':
    run()
