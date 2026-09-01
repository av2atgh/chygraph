"""Chapter 3: what a chygraph built from a real network looks like.

A clustered graph becomes a chygraph by promoting its dense motifs to
complexes, and for a graph with no group data of its own the natural choice is
the maximal cliques.  This script does that for ten real networks and asks the
two questions Chapter 3 needs answered:

  * what ensemble comes out -- the cardinality distribution of the complexes,
    its excess average sbar = <c^2>/<c> - 1, which is what the threshold tensor
    needs to be finite, and the chy-degree of a node;

  * how far from treelike the result is -- `shared_2plus`, the fraction of
    intersecting complex pairs that share two or more nodes, which are exactly
    the pairs the chygraph mapping cannot represent (Fig. 2.7b).

Each network is carried alongside its degree-matched configuration model, which
holds P(k) fixed and destroys the clustering, so the comparison is paired.

Data: netzschleuder CSV dumps cached under LocalNetworkGrowth/figs/data.
Writes tab-real-chygraphs.tex and fig-real-overlap.pdf into the book directory.
"""

import io
import sys
import zipfile
from collections import Counter
from pathlib import Path

import networkx as nx
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'percolation' / 'src'))
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'statmech' / 'src'))
from statmech.region import overlap_profile  # noqa: E402

DATA = Path.home() / 'av2atg' / 'LocalNetworkGrowth' / 'figs' / 'data' / 'netzschleuder'
OUT = Path(__file__).resolve().parent

# name in the book -> cached archive; the last two are the non-biological
# controls, networks with no group-level interaction behind them at all
NETWORKS = [
    ('Collins yeast',   'collins_yeast__collins_yeast'),
    ('Yeast (Y2H)',     'interactome_yeast__interactome_yeast'),
    ('Human (Vidal)',   'interactome_vidal__interactome_vidal'),
    ('Human (Stelzl)',  'interactome_stelzl__interactome_stelzl'),
    ('Human (Figeys)',  'interactome_figeys__interactome_figeys'),
    ('PDZ domains',     'interactome_pdz__interactome_pdz'),
    ('C. elegans WI8',  'celegans_interactomes__WI8'),
    ('C. elegans 2007', 'celegans_interactomes__wi2007'),
    ('Power grid',      'power__power'),
    ('Euroroad',        'euroroad__euroroad'),
]


def load(stem):
    """Read one cached netzschleuder dump as a simple undirected graph."""
    with zipfile.ZipFile(DATA / f'{stem}.csv.zip') as z:
        raw = z.read('edges.csv').decode()
    g = nx.Graph()
    for line in io.StringIO(raw):
        if line.startswith('#') or not line.strip():
            continue
        u, v = line.split(',')[:2]
        u, v = u.strip(), v.strip()
        if u != v:
            g.add_edge(u, v)
    return g


def ensemble(g):
    """The chygraph a graph has when its maximal cliques are its complexes."""
    cliques = [c for c in nx.find_cliques(g) if len(c) >= 2]
    card = np.array([len(c) for c in cliques], dtype=float)
    # excess cardinality seen from a member: the size-biased average, minus the
    # vertex arrived by.  This is the moment the threshold tensor needs finite.
    sbar = (card ** 2).sum() / card.sum() - 1.0
    kappa = Counter()
    for c in cliques:
        for v in c:
            kappa[v] += 1
    kap = np.array([kappa.get(v, 0) for v in g], dtype=float)
    # degree heterogeneity: what the degree-matched control has to work with
    # when it manufactures cliques at hubs (Sec. 3.3)
    d = np.array([x for _, x in g.degree()], dtype=float)
    kexc = (d ** 2).mean() / d.mean()
    prof = overlap_profile(cliques)
    return dict(
        n=g.number_of_nodes(), m=g.number_of_edges(),
        kmean=2 * g.number_of_edges() / g.number_of_nodes(),
        ncx=len(cliques), cmean=card.mean(), cmax=int(card.max()),
        sbar=sbar, kappa=kap.mean(),
        cc=nx.average_clustering(g), kexc=kexc,
        shared=prof['shared_2plus'], cover=prof['edge_cover_mean'],
    )


def control(g, seed):
    """Degree-matched configuration model: same P(k), clustering destroyed."""
    deg = [d for _, d in g.degree()]
    h = nx.configuration_model(deg, seed=seed)
    h = nx.Graph(h)
    h.remove_edges_from(nx.selfloop_edges(h))
    return h


def cardinality_hist(g):
    """Distribution of maximal-clique cardinalities, as a normalised count."""
    card = Counter(len(c) for c in nx.find_cliques(g) if len(c) >= 2)
    tot = sum(card.values())
    xs = sorted(card)
    return np.array(xs), np.array([card[x] / tot for x in xs])


def figure(rows, hists):
    """Two panels: how far from treelike, and what the complexes look like."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(4.3, 4.9))

    # ---- (a) overlap against clustering, real and degree-matched control
    for name, real, ctrl in rows:
        ax1.plot(ctrl['cc'], max(ctrl['shared'], 3e-4), 'o', ms=5,
                 mfc='white', mec='0.55', mew=0.9, zorder=2)
        ax1.plot(real['cc'], max(real['shared'], 3e-4), 'o', ms=5,
                 color='0.15', zorder=3)
    ax1.set_xscale('log'); ax1.set_yscale('log')
    ax1.set_xlabel('clustering coefficient $C$')
    ax1.set_ylabel(r'shared$_{2+}$')
    coll = rows[0][1]
    ax1.annotate('Collins yeast\n(AP/MS)', (coll['cc'], coll['shared']),
                 xytext=(-8, -34), textcoords='offset points',
                 fontsize=7.5, ha='right',
                 arrowprops=dict(arrowstyle='-', lw=0.7, color='0.4'))
    ax1.annotate('binary interactomes,\nroads, power grid',
                 (0.022, 0.014), xytext=(-14, 34), textcoords='offset points',
                 fontsize=7.5, ha='right',
                 arrowprops=dict(arrowstyle='-', lw=0.7, color='0.4'))
    ax1.plot([], [], 'o', ms=5, color='0.15', label='real')
    ax1.plot([], [], 'o', ms=5, mfc='white', mec='0.55', label='degree-matched')
    ax1.legend(fontsize=7.5, frameon=False, loc='lower right')
    ax1.set_ylim(2e-4, 4.0)
    ax1.set_title('(a) how far from treelike', fontsize=8.5, loc='left')

    # ---- (b) what the complexes are
    marks = {'Collins yeast': ('s', '0.15'), 'Yeast (Y2H)': ('o', '0.45'),
             'Power grid': ('^', '0.65')}
    for name, (m, c) in marks.items():
        xs, ps = hists[name]
        ax2.plot(xs, ps, m + '-', ms=4.5, lw=1.0, color=c, mfc='white',
                 mec=c, label=name)
    ax2.set_yscale('log')
    ax2.set_xlabel('complex cardinality $c$')
    ax2.set_ylabel('fraction of complexes')
    ax2.legend(fontsize=7.5, frameon=False)
    ax2.set_title('(b) what the complexes are', fontsize=8.5, loc='left')

    for ax in (ax1, ax2):
        ax.tick_params(labelsize=8)
        for sp in ('top', 'right'):
            ax.spines[sp].set_visible(False)
    fig.tight_layout()
    fig.savefig(OUT / 'fig-real-overlap.pdf')
    print('wrote fig-real-overlap.pdf')


def main():
    rows, hists = [], {}
    for name, stem in NETWORKS:
        g = load(stem)
        g = g.subgraph(max(nx.connected_components(g), key=len)).copy()
        real = ensemble(g)
        ctrls = [ensemble(control(g, s)) for s in (1, 2, 3)]
        ctrl = {k: float(np.mean([c[k] for c in ctrls])) for k in real}
        rows.append((name, real, ctrl))
        hists[name] = cardinality_hist(g)
        print(f"{name:16s} n={real['n']:5d} kbar={real['kmean']:5.2f} "
              f"C={real['cc']:.3f} sbar={real['sbar']:6.2f} "
              f"(ctrl {ctrl['sbar']:5.2f})  shared2+={real['shared']:.3f} "
              f"(ctrl {ctrl['shared']:.3f})  cover={real['cover']:.2f} "
              f"cmax={real['cmax']} <k^2>/<k>={real['kexc']:.1f}")
    write_table(rows)
    figure(rows, hists)
    return rows


def write_table(rows):
    """Table 3.1, sized for the 5x8 trim: seven columns, real against control."""
    lines = [
        r'\begin{tabular}{lrrrrrr}',
        r'\hline\hline',
        r' & & & \multicolumn{2}{c}{$\bar s$} & \multicolumn{2}{c}{shared$_{2+}$}\\',
        r'network & $\ave{k}$ & $C$ & real & ctrl & real & ctrl\\',
        r'\hline',
    ]
    for name, real, ctrl in rows:
        lines.append(
            f"{name} & {real['kmean']:.2f} & {real['cc']:.3f} & "
            f"{real['sbar']:.2f} & {ctrl['sbar']:.2f} & "
            f"{real['shared']:.3f} & {ctrl['shared']:.3f}\\\\")
    lines += [r'\hline\hline', r'\end{tabular}']
    (OUT / 'tab-real-chygraphs.tex').write_text('\n'.join(lines) + '\n')


if __name__ == '__main__':
    main()
