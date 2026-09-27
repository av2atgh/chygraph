"""Chapter 21: the finite-component distribution by layer.

One figure and the numbers Sec. 21.4 quotes.

  fig-components   (a) the outbreak-size distribution of an SIR epidemic
                   with two levels of mixing (households on a Poisson global
                   network, Ch. 6's construction), from the marked map of
                   Eq. (21.6) against a bond-percolation simulation, below
                   and above R* = 1;
                   (b) the finite components of two interactomes, real
                   against the prediction from the degree distribution and
                   from the clique chygraph of Ch. 3.

Checks run first: the Borel distribution on an Erdos-Renyi graph, the exact
series against the transform, Good's inversion against the series, the
joint (individuals, households) distribution against the simulation, and
the mass of the finite components against 1 - S.  Everything comes from
`percolation.components`; the simulation is bond percolation on households
plus global contacts, which is exact for SIR final sizes (Sec. 6.1).
"""

import sys
import time
from collections import Counter
from pathlib import Path

import numpy as np
from sympy import Rational, nsolve, symbols

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'percolation' / 'src'))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from percolation import household_epidemic, hypergraph_giant  # noqa: E402
from percolation.giant import Chygraph, _tables, finite_pgf  # noqa: E402
from percolation.components import (  # noqa: E402
    ComponentDistribution, atom_distribution, atom_finite_fraction, borel,
    default_bins, good_by_layer, joint_clique_model, symbolic_series,
)

OUT = Path(__file__).resolve().parent
DARK, MID, LIGHT = '0.10', '0.45', '0.70'


def _mpl():
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    return plt


def _tidy(ax):
    ax.tick_params(labelsize=8)
    for sp in ('top', 'right'):
        ax.spines[sp].set_visible(False)


# ------------------------------------------------------------------ checks

def check_borel():
    """Erdos-Renyi: the marked map returns the Borel distribution."""
    G = hypergraph_giant(graph=True)
    for c in (0.5, 1.5, 2.5):
        D = ComponentDistribution(G, {'k': c, 'p': 1, 'q': 1})
        P = D.distribution(300, radius=0.99)
        err = np.max(np.abs(P - borel(c, 300)))
        print(f'  Borel at <k> = {c}: max error {err:.1e}, mass {P.sum():.6f} '
              f'against 1 - S = {D.finite_fraction():.6f}')
        assert err < 1e-11


def check_exact_routes():
    """Exact series = transform = Good's inversion on a finite two-layer case."""
    phi, phibar, g, gbar = _tables(2)
    deg = {1: Rational(1, 2), 2: Rational(1, 4), 3: Rational(1, 4)}
    mean = sum(k * v for k, v in deg.items())
    exc = {k - 1: k * v / mean for k, v in deg.items()}
    card = {2: Rational(1, 2), 3: Rational(1, 2)}
    cm = sum(k * v for k, v in card.items())
    cexc = {k - 1: k * v / cm for k, v in card.items()}
    phi[0][1], phibar[0][1] = finite_pgf(deg), finite_pgf(exc)
    g[1][0], gbar[1][0] = finite_pgf(card), finite_pgf(cexc)
    M = Chygraph(phi, phibar, g, gbar)
    S = symbolic_series(M, 6)
    J = ComponentDistribution(M).joint(8)
    worst = max(abs(float(c) - J[mon]) for mon, c in S.items())
    print(f'  exact series against the transform, {len(S)} terms: {worst:.1e}')
    for mon in ((2, 1), (3, 1), (3, 2)):
        gc = good_by_layer(M, (mon[0] - 1, mon[1]))
        print(f'  Good {mon}: {gc} = {float(gc):.6f}, series {float(S[mon]):.6f}')
        assert gc == S[mon]


# ------------------------------------------------- (a) household epidemics

SIZES = {1: Rational(1, 5), 3: Rational(1, 2), 5: Rational(3, 10)}
P_H, K = 0.5, 2.0


def threshold_T():
    T = symbols('T')
    M = household_epidemic(SIZES)
    return float(nsolve(M.theta().subs({'k': K, 'p_H': P_H}), T, 0.2))


def simulate(T, n_households=100_000, seed=0):
    """Bond percolation on households (p_H) plus a Poisson global network
    kept with probability T.  Returns component size per node and the
    number of distinct households per component."""
    from scipy.sparse import coo_matrix
    from scipy.sparse.csgraph import connected_components
    rng = np.random.default_rng(seed)
    sizes = np.array(list(SIZES)); probs = np.array([float(v) for v in SIZES.values()])
    hs = rng.choice(sizes, size=n_households, p=probs)
    n = int(hs.sum())
    household = np.repeat(np.arange(n_households), hs)
    start = np.concatenate([[0], np.cumsum(hs)[:-1]])
    rows, cols = [], []
    for h in np.unique(hs):
        idx = np.nonzero(hs == h)[0]
        if h < 2:
            continue
        pairs = [(a, b) for a in range(h) for b in range(a + 1, h)]
        keep = rng.random((len(idx), len(pairs))) < P_H
        for (a, b), col in zip(pairs, keep.T):
            base = start[idx[col]]
            rows.append(base + a); cols.append(base + b)
    m = rng.poisson(K * T * n / 2)
    u = rng.integers(0, n, m); v = rng.integers(0, n, m)
    ok = u != v
    rows.append(u[ok]); cols.append(v[ok])
    rows = np.concatenate(rows); cols = np.concatenate(cols)
    A = coo_matrix((np.ones(len(rows)), (rows, cols)), shape=(n, n))
    ncomp, label = connected_components(A, directed=False)
    csize = np.bincount(label, minlength=ncomp)
    # households per component
    pairs = np.unique(np.stack([label, household], axis=1), axis=0)
    hcount = np.bincount(pairs[:, 0], minlength=ncomp)
    return csize[label], csize, hcount


def outbreak_panel(ax, seeds=12):
    M = household_epidemic(SIZES)
    Tc = threshold_T()
    print(f'  household threshold T_c = {Tc:.4f} (R* = 1)')
    nmax = 60
    results = {}
    for frac, style in ((0.6, ('o', DARK)), (1.5, ('s', MID))):
        T = frac * Tc
        subs = {'k': K, 'p_H': P_H, 'T': T}
        D = ComponentDistribution(M, subs)
        P = D.distribution(nmax, weights=[1, 0, 0], radius=0.95)
        S = M.node_fraction(subs)
        # simulation: component size of a random node, giant removed
        hist = np.zeros(nmax + 1); tot = 0; giant = []
        for seed in range(seeds):
            per_node, csize, _ = simulate(T, seed=seed)
            n = len(per_node)
            big = csize.max()
            is_giant = (per_node == big) if S > 0.01 else np.zeros(n, bool)
            giant.append(is_giant.mean())
            s = per_node[~is_giant]
            hist += np.bincount(np.minimum(s, nmax), minlength=nmax + 1)[:nmax + 1]
            tot += n
        hist /= tot
        results[frac] = (T, S, np.mean(giant), P, hist)
        ss = np.arange(1, nmax + 1)
        ax.plot(ss, P[1:], '-', color=style[1], lw=1.0,
                label=f'$T = {frac:g}\\,T_c$, theory')
        ax.plot(ss[::2], hist[1:nmax + 1:2], style[0], color=style[1], ms=3,
                mfc='white', mew=0.8, label=f'$T = {frac:g}\\,T_c$, simulation')
        print(f'  T = {frac:g} T_c = {T:.4f}: S = {S:.4f} (theory), giant '
              f'{np.mean(giant):.4f} (simulation); finite mass {P.sum():.4f}')
        for s in (1, 2, 5, 10, 20, 40):
            print(f'    P({s:2d}) theory {P[s]:.3e}  simulation {hist[s]:.3e}')
    ax.set_yscale('log')
    ax.set_xlabel('outbreak size $s$ (individuals)', fontsize=8)
    ax.set_ylabel('$P(s)$', fontsize=8)
    ax.set_ylim(1e-5, 1)
    ax.legend(fontsize=6.5, frameon=False, loc='upper right')
    _tidy(ax)
    return results


def check_joint(seeds=4):
    """Mean number of households touched by an outbreak of s individuals."""
    M = household_epidemic(SIZES)
    T = 0.6 * threshold_T()
    D = ComponentDistribution(M, {'k': K, 'p_H': P_H, 'T': T})
    J = D.joint((30, 30), layers=(0, 1), radius=0.9)
    ss = np.arange(31)
    mean_h = (J * ss[None, :]).sum(axis=1) / np.maximum(J.sum(axis=1), 1e-300)
    sim = Counter(); cnt = Counter()
    for seed in range(seeds):
        _, csize, hcount = simulate(T, seed=seed)
        for s, h in zip(csize, hcount):
            if s <= 30:
                sim[s] += h; cnt[s] += 1
    print('  households touched by an outbreak of s individuals '
          '(theory / simulation):')
    for s in (1, 2, 3, 5, 8, 12, 20):
        print(f'    s = {s:2d}: {mean_h[s]:.3f} / {sim[s] / cnt[s]:.3f}  '
              f'({cnt[s]} outbreaks)')


# ------------------------------------------------ (b) interactome components

INTERACTOMES = [('Yeast (Y2H)', 'interactome_yeast__interactome_yeast'),
                ('C. elegans (WI-2007)', 'celegans_interactomes__wi2007')]

# every cached network with at least five components, for the table
TABLE = [('Yeast (Y2H)', 'interactome_yeast__interactome_yeast'),
         ('C. elegans WI-2007', 'celegans_interactomes__wi2007'),
         ('C. elegans WI8', 'celegans_interactomes__WI8'),
         ('PDZ', 'interactome_pdz__interactome_pdz'),
         ('Stelzl', 'interactome_stelzl__interactome_stelzl'),
         ('Vidal', 'interactome_vidal__interactome_vidal'),
         ('Figeys', 'interactome_figeys__interactome_figeys'),
         ('Collins yeast', 'collins_yeast__collins_yeast'),
         ('Euroroad', 'euroroad__euroroad')]


def graph_model(g):
    import networkx as nx
    deg = Counter(d for _, d in g.degree())
    n = g.number_of_nodes()
    pk = {k: Rational(c, n) for k, c in deg.items()}
    mean = sum(k * v for k, v in pk.items())
    exc = {k - 1: k * v / mean for k, v in pk.items() if k >= 1}
    return hypergraph_giant(degree=finite_pgf(pk), excess_degree=finite_pgf(exc),
                            poisson=False, graph=True)


def clique_model(g):
    import networkx as nx
    cliques = [c for c in nx.find_cliques(g) if len(c) >= 2]
    kappa = Counter()
    for c in cliques:
        for v in c:
            kappa[v] += 1
    n = g.number_of_nodes()
    kap = Counter(kappa.get(v, 0) for v in g)
    pk = {k: Rational(c, n) for k, c in kap.items()}
    mean = sum(k * v for k, v in pk.items())
    exc = {k - 1: k * v / mean for k, v in pk.items() if k >= 1}
    card = Counter(len(c) for c in cliques)
    m = len(cliques)
    pc = {c: Rational(v, m) for c, v in card.items()}
    cmean = sum(c * v for c, v in pc.items())
    biased = {c: c * v / cmean for c, v in pc.items()}
    phi, phibar, g_, gbar = _tables(2)
    phi[0][1], phibar[0][1] = finite_pgf(pk), finite_pgf(exc)
    g_[1][0] = finite_pgf(biased)
    gbar[1][0] = finite_pgf({c - 1: v for c, v in biased.items()})
    return Chygraph(phi, phibar, g_, gbar)


def interactome_panel(ax, nmax=16):
    import networkx as nx
    from real_chygraphs import load
    out = {}
    for (name, stem), style in zip(INTERACTOMES, (('o', DARK), ('s', MID))):
        g = load(stem)
        n = g.number_of_nodes()
        comps = sorted((len(c) for c in nx.connected_components(g)), reverse=True)
        lcc = comps[0] / n
        real = np.zeros(nmax + 1)
        for s in comps[1:]:
            if s <= nmax:
                real[s] += s / n
        Pg = ComponentDistribution(graph_model(g), {'p': 1, 'q': 1}).distribution(nmax, radius=0.9)
        Sg = 1 - ComponentDistribution(graph_model(g), {'p': 1, 'q': 1}).finite_fraction()
        Dc = ComponentDistribution(clique_model(g))
        Pc = Dc.distribution(nmax, radius=0.9)
        Sc = 1 - Dc.finite_fraction()
        Dj = ComponentDistribution(joint_clique_model(g, complexes=merged_family(g)))
        Pj = atom_distribution(Dj, nmax)
        out[name] = dict(n=n, lcc=lcc, Sg=Sg, Sc=Sc, real=real, Pg=Pg, Pc=Pc,
                         Pj=Pj, Sj=1 - atom_finite_fraction(Dj), small=sum(comps[1:]))
        print(f'  {name}: n = {n}, largest component {lcc:.3f} of the nodes; '
              f'graph model S = {Sg:.3f}, clique chygraph S = {Sc:.3f}; '
              f'{sum(comps[1:])} nodes in {len(comps) - 1} finite components')
        for s in (1, 2, 3, 4, 6, 8):
            print(f'    P({s}): real {real[s]:.4f}  graph {Pg[s]:.4f}  clique {Pc[s]:.4f}')
        ss = np.arange(1, nmax + 1)
        ax.plot(ss, Pg[1:], ':', color=style[1], lw=1.0)
        ax.plot(ss, Pc[1:], '-', color=style[1], lw=1.0)
        ax.plot(ss, Pj[1:], '--', color=style[1], lw=1.0)
        mask = real[1:] > 0
        ax.plot(ss[mask], real[1:][mask], style[0], color=style[1], ms=3.2,
                mfc='white', mew=0.8, label=name)
    ax.plot([], [], ':', color=DARK, lw=1.0, label='degree distribution')
    ax.plot([], [], '-', color=DARK, lw=1.0, label='clique chygraph')
    ax.plot([], [], '--', color=DARK, lw=1.0, label='merged, joint laws')
    ax.set_yscale('log')
    ax.set_xlabel('finite component size $s$', fontsize=8)
    ax.set_ylabel('fraction of nodes', fontsize=8)
    ax.set_ylim(1e-3, 1)
    ax.set_xlim(0.5, 16)
    ax.set_xticks(range(2, 17, 2))
    ax.legend(fontsize=6.5, frameon=False)
    _tidy(ax)
    return out


def table(nmax=30):
    """Table 21.1: real against predicted components, nine networks."""
    import networkx as nx
    from real_chygraphs import load
    rows = []
    for name, stem in TABLE:
        g = load(stem)
        n = g.number_of_nodes()
        comps = sorted((len(c) for c in nx.connected_components(g)), reverse=True)
        real = np.zeros(nmax + 1)
        for s in comps[1:]:
            if s <= nmax:
                real[s] += s / n
        Dg = ComponentDistribution(graph_model(g), {'p': 1, 'q': 1})
        Pg = Dg.distribution(nmax, radius=0.9)
        Dc = ComponentDistribution(clique_model(g))
        Pc = Dc.distribution(nmax, radius=0.9)
        rows.append((name, n, comps[0] / n, 1 - Dg.finite_fraction(),
                     1 - Dc.finite_fraction(), real[2], Pg[2], Pc[2],
                     real[3], Pg[3], Pc[3]))
        print(f'  {name:20s} n={n:5d} largest {comps[0] / n:.3f} | S graph '
              f'{1 - Dg.finite_fraction():.3f} clique {1 - Dc.finite_fraction():.3f} | '
              f'P2 {real[2]:.3f} {Pg[2]:.3f} {Pc[2]:.3f} | P3 {real[3]:.3f} {Pg[3]:.3f} {Pc[3]:.3f}')
    with open(OUT / 'tab-components.tex', 'w') as f:
        f.write('% generated by figs/components.py; do not edit\n')
        f.write('\\begin{tabular}{lrrrrrrrr}\n\\hline\\hline\n')
        f.write('network & $n$ & largest & \\multicolumn{2}{c}{$S$} '
                '& \\multicolumn{2}{c}{$P(2)$} & \\multicolumn{2}{c}{$P(3)$}\\\\\n')
        f.write(' & & real & graph & clique & real & clique & real & clique\\\\\n\\hline\n')
        for r in rows:
            name, n, lcc, Sg, Sc, r2, g2, c2, r3, g3, c3 = r
            vals = ' & '.join(f'{v:.3f}' for v in (lcc, Sg, Sc, r2, c2, r3, c3))
            f.write(f'{name} & {n} & {vals}\\\\\n')
        f.write('\\hline\\hline\n\\end{tabular}\n')
    print('  wrote tab-components.tex')


def check_joint_reduces(stem='interactome_yeast__interactome_yeast', nmax=16):
    """The joint model with independent class draws is the clique model."""
    from real_chygraphs import load
    g = load(stem)
    Dc = ComponentDistribution(clique_model(g))
    Dt = ComponentDistribution(joint_clique_model(g, thinned=True))
    Pc, Pt = Dc.distribution(nmax, radius=0.9), atom_distribution(Dt, nmax)
    dS = abs(Dc.finite_fraction() - atom_finite_fraction(Dt))
    dP = float(np.abs(Pc - Pt).max())
    print(f'  thinned joint model against the clique model: S differs by {dS:.1e}, '
          f'P(s) by {dP:.1e}')
    assert dS < 1e-8 and dP < 1e-8


def merged_family(g):
    """Chapter 16's merge closure of the maximal cliques of ``g``."""
    import sys
    import types
    import networkx as nx
    sys.modules.setdefault('hrg', types.SimpleNamespace(hrg_calibrated=None))
    from merge import merge_closure
    merged, _ = merge_closure([frozenset(c) for c in nx.find_cliques(g) if len(c) >= 2])
    return [tuple(sorted(c)) for c in merged]


def table_joint(nmax=30):
    """Table 21.2: the measured joint law against the clique model, on the
    maximal cliques and on their merge closure."""
    import networkx as nx
    from real_chygraphs import load
    rows = []
    for name, stem in TABLE:
        g = load(stem)
        n = g.number_of_nodes()
        comps = sorted((len(c) for c in nx.connected_components(g)), reverse=True)
        real = np.zeros(nmax + 1)
        for s in comps[1:]:
            if s <= nmax:
                real[s] += s / n
        Dc = ComponentDistribution(clique_model(g))
        Pc = Dc.distribution(nmax, radius=0.9)
        Dn = ComponentDistribution(joint_clique_model(g, types=[(1, 10 ** 9)]))
        Pn = atom_distribution(Dn, nmax)
        Dj = ComponentDistribution(joint_clique_model(g))
        Pj = atom_distribution(Dj, nmax)
        merged = merged_family(g)
        Dm = ComponentDistribution(joint_clique_model(g, complexes=merged))
        Pm = atom_distribution(Dm, nmax)
        S = [comps[0] / n, 1 - Dc.finite_fraction(), 1 - atom_finite_fraction(Dn),
             1 - atom_finite_fraction(Dj), 1 - atom_finite_fraction(Dm)]
        P2 = [real[2], Pc[2], Pn[2], Pj[2], Pm[2]]
        tail = [real[4:].sum(), Pc[4:].sum(), Pn[4:].sum(), Pj[4:].sum(), Pm[4:].sum()]
        rows.append((name, len(merged), max(len(c) for c in merged), S, P2, tail))
        print(f'  {name:20s} ({len(merged)} meta-complexes, largest {max(len(c) for c in merged)}): '
              f'S ' + ' '.join(f'{v:.3f}' for v in S) + ' | P2 ' + ' '.join(f'{v:.3f}' for v in P2)
              + f' | tail(4..{nmax}) ' + ' '.join(f'{v:.3f}' for v in tail)
              + '   [real, independent, node law, both laws, merged + both laws]')
    with open(OUT / 'tab-components-joint.tex', 'w') as f:
        f.write('% generated by figs/components.py; do not edit\n')
        f.write('\\begin{tabular}{lrrrrrrrrrrrr}\n\\hline\\hline\n')
        f.write('network & \\multicolumn{4}{c}{$S$} & \\multicolumn{4}{c}{$P(2)$} '
                '& \\multicolumn{4}{c}{$\\sum_{s\\ge4}P(s)$}\\\\\n')
        f.write(' & real & ind. & joint & merged & real & ind. & joint & merged '
                '& real & ind. & joint & merged\\\\\n\\hline\n')
        for name, nm, top, S, P2, tail in rows:
            vals = [S[0], S[1], S[3], S[4], P2[0], P2[1], P2[3], P2[4], tail[0], tail[1], tail[3], tail[4]]
            f.write(f'{name} & ' + ' & '.join(f'{v:.3f}' for v in vals) + '\\\\\n')
        f.write('\\hline\\hline\n\\end{tabular}\n')
    print('  wrote tab-components-joint.tex')
    return rows


def table_periphery(nmax=30):
    """Table 21.3: the remaining tail deficit is the periphery.  The tail
    of the merged joint model with two and with five atom types, and the
    neighbour statistics of chy-degree-one nodes inside the finite
    components against inside the giant."""
    import networkx as nx
    from real_chygraphs import load
    five = [(1, 1), (2, 2), (3, 3), (4, 5), (6, 10 ** 9)]
    rows = []
    for name, stem in TABLE:
        if name in ('Collins yeast', 'Euroroad'):
            continue
        g = load(stem)
        n = g.number_of_nodes()
        merged = merged_family(g)
        comps = sorted((len(c) for c in nx.connected_components(g)), reverse=True)
        real = sum(s for s in comps[1:] if 4 <= s <= nmax) / n
        tails = []
        for ty in ([(1, 1), (2, 10 ** 9)], five):
            D = ComponentDistribution(joint_clique_model(g, complexes=merged, types=ty))
            tails.append(float(atom_distribution(D, nmax)[4:].sum()))
        kappa = Counter(v for c in merged for v in c)
        giant = max(nx.connected_components(g), key=len)
        nb = {v: set() for v in g}
        for c in merged:
            for u in c:
                nb[u].update(x for x in c if x != u)

        def low(nodes):
            vs = [v for v in nodes if kappa[v] == 1]
            return float(np.mean([np.mean([kappa[x] <= 2 for x in nb[v]]) for v in vs]))

        fin = [v for v in g if v not in giant]
        share = np.mean([v not in giant for v in g if kappa[v] <= 2])
        rows.append((name, real, tails[0], tails[1], share, low(fin), low(giant)))
        print(f'  {name:20s} tail real {real:.3f}, two types {tails[0]:.3f}, five types {tails[1]:.3f}; '
              f'share of low-chy-degree nodes in finite components {share:.2f}; '
              f'P(neighbour has chy-degree <= 2 | chy-degree 1): finite {low(fin):.2f}, giant {low(giant):.2f}')
    with open(OUT / 'tab-components-periphery.tex', 'w') as f:
        f.write('% generated by figs/components.py; do not edit\n')
        f.write('\\begin{tabular}{lrrrrrr}\n\\hline\\hline\n')
        f.write('network & \\multicolumn{3}{c}{$\\sum_{s\\ge4}P(s)$} & low-degree share & '
                '\\multicolumn{2}{c}{$P(\\kappa_{\\mathrm{nb}}\\le2\\mid\\kappa=1)$}\\\\\n')
        f.write(' & real & 2 types & 5 types & in finite & finite & giant\\\\\n\\hline\n')
        for name, *vals in rows:
            f.write(f'{name} & ' + ' & '.join(f'{v:.3f}' if i < 3 else f'{v:.2f}'
                                              for i, v in enumerate(vals)) + '\\\\\n')
        f.write('\\hline\\hline\n\\end{tabular}\n')
    print('  wrote tab-components-periphery.tex')


def figure():
    plt = _mpl()
    fig, axes = plt.subplots(2, 1, figsize=(3.4, 4.0))
    outbreak_panel(axes[0])
    interactome_panel(axes[1])
    for ax, tag in zip(axes, 'ab'):
        ax.text(-0.18, 1.02, f'({tag})', transform=ax.transAxes, fontsize=9)
    fig.tight_layout()
    fig.savefig(OUT / 'fig-components.pdf')
    print('  wrote fig-components.pdf')


if __name__ == '__main__':
    t0 = time.time()
    print('Borel anchor:'); check_borel()
    print('exact routes:'); check_exact_routes()
    print('joint distribution:'); check_joint()
    print('table:'); table()
    print('the joint laws reduce:'); check_joint_reduces()
    print('table, joint laws:'); table_joint()
    print('table, the periphery:'); table_periphery()
    print('figure:'); figure()
    print(f'done in {time.time() - t0:.0f} s')
