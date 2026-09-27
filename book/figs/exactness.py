"""Chapter on exactness: where the clique family stops being locally acyclic.

  fig-acyclicity      phi, the fraction of vertices on a chordless cycle of
                      length four or five, and the largest connected piece of
                      that set, against mean degree, for random geometric
                      graphs in d = 1, 2, 3 and hyperbolic random graphs at
                      three tau (n = 20000)
  fig-acyclicity-fss  the largest piece against mean degree at three sizes,
                      d = 2 and d = 3: acyclicity percolation
  tab-real-acyclicity phi on sixteen real networks against a degree-matched
                      rewired control, with the GYO residue alongside

Data are read from ../statmech/probe/results/acyclicity_{rgg,hrg,fss}.csv and
real_acyclicity.csv, the cached outputs of ../statmech/probe/acyclicity.py and
real_acyclicity.py.  The checks recompute the small things: GYO reduction
against chordality, the line against the ring, the short-cycle finder against
networkx's enumeration, and the join tree against enumeration on 1D and 2D
geometric graphs of sixteen vertices.
"""
import csv
import sys
from collections import defaultdict
from pathlib import Path

import networkx as nx
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'percolation' / 'src'))
sys.path.insert(0, str(ROOT / 'statmech' / 'src'))
sys.path.insert(0, str(ROOT / 'statmech' / 'probe'))
import acyclicity as A  # noqa: E402

OUT = Path(__file__).resolve().parent
PROBE = ROOT / 'statmech' / 'probe' / 'results'
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


def _read(name):
    path = PROBE / f'acyclicity_{name}.csv'
    if not path.exists():
        return []
    return list(csv.DictReader(open(path)))


def _avg(rows, keys, ycols):
    """Average ycols over seeds, grouped by keys; returns {key: {ycol: mean}}."""
    acc = defaultdict(lambda: defaultdict(list))
    for r in rows:
        k = tuple(float(r[c]) for c in keys)
        for y in ycols:
            acc[k][y].append(float(r[y]))
    return {k: {y: float(np.mean(v[y])) for y in ycols} for k, v in acc.items()}


# ------------------------------------------------------------- figures

STYLES = {
    ('rgg', 1.0): dict(marker='s', color=DARK, ls='-', label='$d=1$'),
    ('rgg', 2.0): dict(marker='o', color=DARK, ls='-', label='$d=2$'),
    ('rgg', 3.0): dict(marker='^', color=DARK, ls='-', label='$d=3$'),
    ('hrg', 2.5): dict(marker='o', color=LIGHT, ls='--', label=r'HRG $\tau=2.5$'),
    ('hrg', 3.0): dict(marker='^', color=LIGHT, ls='--', label=r'HRG $\tau=3$'),
    ('hrg', 4.0): dict(marker='s', color=LIGHT, ls='--', label=r'HRG $\tau=4$'),
}


def figure_acyclicity():
    plt = _mpl()
    rows = _read('rgg') + _read('hrg')
    av = _avg(rows, ('kbar',), ('phi_short', 'giant_short')) if False else None
    fig, (ax, bx) = plt.subplots(2, 1, figsize=(4.3, 4.6), sharex=True)
    for (ens, p), st in STYLES.items():
        sub = [r for r in rows if r['ensemble'] == ens and float(r['param']) == p]
        if not sub:
            continue
        g = _avg(sub, ('kbar',), ('phi_short', 'giant_short'))
        ks = sorted(g)
        ax.plot([k[0] for k in ks], [g[k]['phi_short'] for k in ks], ms=4, mfc='white', **st)
        bx.plot([k[0] for k in ks], [g[k]['giant_short'] for k in ks], ms=4, mfc='white', **st)
    ax.set_ylabel(r'$\phi$  (on a short chordless cycle)', fontsize=8)
    bx.set_ylabel('largest piece of that set', fontsize=8)
    bx.set_xlabel(r'mean degree $\langle k\rangle$', fontsize=8)
    ax.set_xscale('log')
    ax.legend(fontsize=7, frameon=False, ncol=2)
    for a in (ax, bx):
        _tidy(a)
        a.set_ylim(-0.02, 1.02)
    fig.tight_layout()
    fig.savefig(OUT / 'fig-acyclicity.pdf')
    plt.close(fig)


def figure_fss():
    plt = _mpl()
    rows = _read('fss')
    # stacked, not side by side: two panels across the 4.05in measure are
    # illegible at the 5x8 trim
    fig, axes = plt.subplots(2, 1, figsize=(4.3, 4.6), sharex=False)
    shades = {5000.0: '0.80', 20000.0: LIGHT, 80000.0: MID, 320000.0: DARK}
    for ax, d in zip(axes, (2.0, 3.0)):
        sub = [r for r in rows if r['ensemble'] == 'rgg' and float(r['param']) == d]
        g = _avg(sub, ('n', 'kbar'), ('giant_short', 'gc'))
        for n in sorted({k[0] for k in g}):
            ks = sorted(k for k in g if k[0] == n)
            ax.plot([k[1] for k in ks], [g[k]['giant_short'] for k in ks], 'o-', ms=4,
                    color=shades.get(n, DARK), mfc='white', label=f'$n={int(n)}$')
        kc = {2.0: 4.512, 3.0: 2.736}[d]
        ax.axvline(kc, color=LIGHT, lw=0.8, ls=':')
        ax.set_title(f'$d={int(d)}$', fontsize=9)
        ax.set_ylabel('largest piece, fraction of $n$', fontsize=8)
        _tidy(ax)
    axes[1].set_xlabel(r'mean degree $\langle k\rangle$', fontsize=8)
    axes[0].legend(fontsize=7, frameon=False)
    fig.tight_layout()
    fig.savefig(OUT / 'fig-acyclicity-fss.pdf')
    plt.close(fig)


# --------------------------------------------------------------- table

def table_real():
    rows = list(csv.DictReader(open(PROBE / 'real_acyclicity.csv')))
    rows.sort(key=lambda r: float(r['kbar']))
    lines = [r'\begin{tabular}{lrrrrrr}', r'\hline\hline',
             r' & & & \multicolumn{2}{c}{$\phi$} & & \\',
             r'network & $\ave{k}$ & $C$ & real & rewired & piece & residue\\',
             r'\hline']
    for r in rows:
        lines.append('%s & %.2f & %.3f & %.3f & %.3f & %.3f & %.3f\\\\' % (
            r['network'], float(r['kbar']), float(r['transitivity']), float(r['phi_short']),
            float(r['phi_short_ctrl']), float(r['giant_short']), float(r['residue'])))
    lines += [r'\hline\hline', r'\end{tabular}']
    (OUT / 'tab-real-acyclicity.tex').write_text('\n'.join(lines) + '\n')


def table_growing():
    """Table 17.x: phi on rewirings of the real degree sequences replicated
    k times (probe/phi_growing.py)."""
    rows = list(csv.DictReader(open(PROBE / 'phi_growing.csv')))
    nets = []
    for r in rows:
        if r['network'] not in nets:
            nets.append(r['network'])
    lines = [r'\begin{tabular}{lrrrrr}', r'\hline\hline',
             r'network & real & $k=1$ & $k=2$ & $k=4$ & $k=8$\\', r'\hline']
    for net in nets:
        d = {int(r['k']): float(r['phi']) for r in rows if r['network'] == net}
        lines.append('%s & %.2f & %.2f & %.2f & %.2f & %.2f\\\\' % (
            net, d[0], d[1], d[2], d[4], d[8]))
    lines += [r'\hline\hline', r'\end{tabular}']
    (OUT / 'tab-phi-growing.tex').write_text('\n'.join(lines) + '\n')
    print('  wrote tab-phi-growing.tex')


def fss_collapse():
    """The largest piece at four sizes in d = 2, scaled with the exponents of
    two-dimensional percolation: the crossing locates the threshold."""
    rows = _read('fss')
    sub = [r for r in rows if r['ensemble'] == 'rgg' and float(r['param']) == 2.0]
    g = _avg(sub, ('n', 'kbar'), ('giant_short',))
    ns = sorted({k[0] for k in g})
    beta_nu = 5.0 / 48.0
    print('  d = 2, largest piece x n^{beta/nu} (beta/nu = 5/48):')
    for kb in sorted({k[1] for k in g}):
        vals = [(n, g[(n, kb)]['giant_short'] * n ** beta_nu) for n in ns if (n, kb) in g]
        print(f'    kbar = {kb:4.1f}: ' + '  '.join(f'n={int(n)}: {v:.3f}' for n, v in vals))
    # crossing between consecutive sizes, by linear interpolation in kbar
    out = {}
    for n1, n2 in zip(ns[:-1], ns[1:]):
        ks = sorted(k for (n, k) in g if n == n1 and (n2, k) in g)
        diff = [g[(n1, k)]['giant_short'] * n1 ** beta_nu - g[(n2, k)]['giant_short'] * n2 ** beta_nu for k in ks]
        for i in range(len(ks) - 1):
            if diff[i] * diff[i + 1] < 0:
                kc = ks[i] + (ks[i + 1] - ks[i]) * diff[i] / (diff[i] - diff[i + 1])
                out[(n1, n2)] = kc
                print(f'    crossing of n = {int(n1)} and {int(n2)}: kbar_a = {kc:.2f}')
                break
    return out


# -------------------------------------------------------------- checks

def check_gyo_is_chordality(trials=200):
    ok = 0
    for s in range(trials):
        g = nx.gnp_random_graph(12, 0.35, seed=s)
        cl = [c for c in nx.find_cliques(g) if len(c) >= 2]
        ok += (len(A.gyo_residue(cl, 12)) == 0) == nx.is_chordal(g)
    print(f'GYO residue empty iff chordal: {ok}/{trials}')
    assert ok == trials


def check_short_cycles(trials=150):
    ok = 0
    for s in range(trials):
        g = nx.gnp_random_graph(14, 0.3, seed=s)
        ref = set()
        for c in nx.chordless_cycles(g, length_bound=5):
            if len(c) >= 4:
                ref.update(c)
        ok += ref == A.on_short_chordless_cycle(g)
    print(f'short-cycle finder agrees with nx.chordless_cycles(<=5): {ok}/{trials}')
    assert ok == trials


def check_line_vs_ring():
    rng = np.random.default_rng(0)
    line = A.rgg_torus(4000, 1, 30.0, rng)
    ring = A.rgg_torus(4000, 1, 30.0, rng, wrap=True)
    print(f'1D line chordal: {nx.is_chordal(line)}; 1D ring chordal: {nx.is_chordal(ring)} '
          f'(both connected: {nx.is_connected(line)}, {nx.is_connected(ring)})')
    assert nx.is_chordal(line) and not nx.is_chordal(ring)


def check_join_tree_exact(n=16, kbar=5.0, per=6):
    """Node Bethe, chygraph Bethe and GBP over the clique region graph against
    enumeration, on 1D (chordal) and 2D geometric graphs of n vertices."""
    from statmech.gbp import GBP, exact_log_Z, ising_factors, static_log_Z
    from statmech.region import RegionGraph
    rng = np.random.default_rng(5)
    for d in (1, 2):
        for bJ in (0.3, 0.8):
            errs = defaultdict(list)
            chordal = []
            for _ in range(per):
                while True:
                    G = A.rgg_torus(n, d, kbar, rng)
                    if nx.is_connected(G):
                        break
                f = ising_factors(list(G.edges()), bJ)
                ex = exact_log_Z(f, range(n))
                cl = [frozenset(c) for c in nx.find_cliques(G)]
                rg = RegionGraph(cl)
                node_rg = RegionGraph([frozenset(e) for e in G.edges()] + [frozenset([v]) for v in G])
                errs['node Bethe'].append(static_log_Z(node_rg.bethe_counting(), f) - ex)
                errs['chygraph Bethe'].append(static_log_Z(rg.bethe_counting(), f) - ex)
                errs['GBP (join tree)'].append(GBP(rg, f, damping=0.7).run(6000).log_Z() - ex)
                chordal.append(nx.is_chordal(G))
            print(f'd={d} bJ={bJ}: chordal {sum(chordal)}/{per}; max |error| in ln Z: ' +
                  ', '.join(f'{k} {max(abs(x) for x in v):.4f}' for k, v in errs.items()))
            if d == 1:
                gbp_on_chordal = [abs(e) for e, c in zip(errs['GBP (join tree)'], chordal) if c]
                assert max(gbp_on_chordal) < 1e-6


def check_karrer_newman(s=1.0, t=1.0, sizes=(5000, 20000, 80000)):
    """The placed ensemble: the number of vertices on short chordless cycles
    is a constant, so phi falls as 1/n and the ensemble is exact at every b."""
    for n in sizes:
        rng = np.random.default_rng(1)
        ne, nt = rng.poisson(s * n / 2), rng.poisson(t * n / 3)
        G = nx.Graph()
        G.add_nodes_from(range(n))
        E = rng.integers(0, n, (ne, 2))
        G.add_edges_from((int(a), int(b)) for a, b in E if a != b)
        for a, b, c in rng.integers(0, n, (nt, 3)):
            if len({a, b, c}) == 3:
                G.add_edges_from([(int(a), int(b)), (int(b), int(c)), (int(a), int(c))])
        bad = A.on_short_chordless_cycle(G)
        print(f'Karrer-Newman s={s} t={t}: n={n} kbar={2 * G.number_of_edges() / n:.3f} '
              f'vertices on short chordless cycles {len(bad)} (phi={len(bad) / n:.5f})')


if __name__ == '__main__':
    check_gyo_is_chordality()
    check_karrer_newman()
    check_short_cycles()
    check_line_vs_ring()
    check_join_tree_exact()
    figure_acyclicity()
    figure_fss()
    fss_collapse()
    table_growing()
    table_real()
    print('wrote fig-acyclicity.pdf, fig-acyclicity-fss.pdf, tab-real-acyclicity.tex')
