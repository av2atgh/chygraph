"""Chapter 13: what overlapping complexes cost, and what recovers it.

  fig-gbp        the same three schemes over 60 maximal-clique region graphs of
                 hyperbolic random graphs, split by whether the clique
                 structure is chordal

  fig-gbp-real   the same experiment on the clustered neighbourhoods of six
                 real networks, and the chordal fraction network by network

Chapter 13's first figure, two triangles sharing an edge, is TikZ in
`overlap.tex`; its table of two-triangle errors is printed by
`check_two_triangles` below.

The 60-instance data is read from `../probe/results/gbp_cliques.json`, the
cached output of `../probe/gbp_cliques.py`, and the real-network data from
`../probe/results/gbp_real.json`, the cached output of `../probe/gbp_real.py`.
The two-triangle numbers are recomputed here, being a four-spin enumeration.
"""

import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path.home() / 'av2atg' / 'chygraph_statmech' / 'src'))
from chygraph_statmech import gbp, region  # noqa: E402

OUT = Path(__file__).resolve().parent
PROBE = Path.home() / 'av2atg' / 'chygraph_statmech' / 'probe' / 'results'
DARK, MID, LIGHT = '0.10', '0.45', '0.70'

# Two triangles sharing the edge {1,2}.
TWO_TRIANGLES = [(0, 1, 2), (1, 2, 3)]
COUPLINGS = (0.2, 0.5, 1.0, 2.0)

# The residual distribution over the non-chordal runs has a clear two-order gap
# between 3.2e-6 and 2.0e-2; anything below this counts as converged.
CONVERGED = 1e-4


def _mpl():
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    return plt


def _tidy(ax):
    ax.tick_params(labelsize=8)
    for sp in ('top', 'right'):
        ax.spines[sp].set_visible(False)


# --------------------------------------------------- (1) the two countings
def check_counting():
    """Mobius against Bethe on two triangles sharing an edge."""
    rg = region.RegionGraph(TWO_TRIANGLES)
    mob = {frozenset(k): v for k, v in rg.counting.items() if v}
    bet = {frozenset(k): v for k, v in rg.bethe_counting().items() if v}
    print(f'  Mobius: {sorted((sorted(k), v) for k, v in mob.items())}')
    print(f'  Bethe : {sorted((sorted(k), v) for k, v in bet.items())}')
    print('    both count every node once; only one subtracts the shared edge')
    # the factor-coverage test: the shared bond {1,2}
    shared = frozenset({1, 2})
    for name, cnt in (('Mobius', rg.counting), ('Bethe', rg.bethe_counting())):
        cov = sum(v for k, v in cnt.items() if shared <= frozenset(k))
        print(f'    {name}: the shared bond is covered {cov} time'
              f'{"" if cov == 1 else "s"}')
        assert cov == (1 if name == 'Mobius' else 2)


def check_two_triangles():
    """The table: three schemes against an exact four-spin enumeration."""
    rg = region.RegionGraph(TWO_TRIANGLES)
    edges = gbp.clique_edges(TWO_TRIANGLES)
    print('   beta J     exact        Bethe      Kikuchi          GBP')
    rows = []
    for bJ in COUPLINGS:
        lf = gbp.ising_factors(edges, bJ)
        ex = gbp.exact_log_Z(lf, nodes=[0, 1, 2, 3])
        be = gbp.static_log_Z(rg.bethe_counting(), lf) - ex
        ki = gbp.static_log_Z(rg.counting, lf) - ex
        g = gbp.GBP(rg, lf, damping=0.5)
        g.run()
        gb = g.log_Z() - ex
        assert be > 0 and ki < 0 and abs(gb) < 1e-12, (bJ, be, ki, gb)
        rows.append((bJ, ex, be, ki, gb))
        print(f'  {bJ:>6.1f}  {ex:>9.4f}  {be:>+11.2e}  {ki:>+11.2e}  '
              f'{gb:>+11.2e}')
    print('    the two countings err in opposite directions, and the messages')
    print('    remove the whole of it: the region family is a junction tree')
    return rows


# ------------------------------------------------- (2) sixty clique graphs
def _runs():
    return json.load(open(PROBE / 'gbp_cliques.json'))


def check_gbp_instances():
    d = _runs()
    ch = [r for r in d if r['chordal']]
    nc = [r for r in d if not r['chordal']]
    conv = [r for r in nc if r['residual'] < CONVERGED]
    bad = [r for r in nc if r['residual'] >= CONVERGED]
    worst = max(abs(r['gbp']) for r in ch)
    assert worst < 1e-11, worst
    print(f'  {len(d)} maximal-clique region graphs of hyperbolic random '
          f'graphs, n = 14, 18, 20')
    print(f'    chordal      {len(ch):>3}:  exact, worst error {worst:.1e}')
    print(f'    non-chordal  {len(nc):>3}:  {len(conv)} converge, '
          f'{len(bad)} do not')
    g = [abs(r['gbp']) for r in conv]
    k = [abs(r['kikuchi']) for r in conv]
    b = [abs(r['bethe']) for r in conv]
    print(f'      where they converge:  GBP {min(g):.1e} to {max(g):.1e}')
    print(f'                            M\'obius {min(k):.2f} to {max(k):.2f}')
    print(f'                            Bethe {min(b):.2f} to {max(b):.1f}')
    orders = [np.log10(abs(r['kikuchi']) / abs(r['gbp'])) for r in conv]
    print(f'      GBP better by {min(orders):.1f} to {max(orders):.1f} orders')
    print(f'      consistency sum_P b_P = b_R held to '
          f'{max(r["consistency"] for r in conv):.1e}')
    print(f'    the {len(bad)} that do not converge do so at damping 0.999; '
          f'their residual\n    says so, and no number from them is quoted')


def figure_gbp():
    plt = _mpl()
    d = _runs()
    fig, ax = plt.subplots(figsize=(4.3, 2.9))
    groups = (('chordal', [r for r in d if r['chordal']], 'o', DARK),
              ('non-chordal,\nconverged',
               [r for r in d if not r['chordal']
                and r['residual'] < CONVERGED], 's', MID),
              ('non-chordal,\nno fixed point',
               [r for r in d if not r['chordal']
                and r['residual'] >= CONVERGED], 'x', LIGHT))
    floor = 1e-14
    rng = np.random.default_rng(0)
    for i, (name, rows, mk, col) in enumerate(groups):
        for j, key in enumerate(('bethe', 'kikuchi', 'gbp')):
            x = i * 3.4 + j + rng.uniform(-0.16, 0.16, len(rows))
            y = np.maximum([abs(r[key]) for r in rows], floor)
            ax.semilogy(x, y, mk, ms=3.2, mew=0.9,
                        mfc='white' if mk != 'x' else col, color=col)
            # the four sound fixed points that came out no better than the
            # static estimate they were supposed to improve: filled, because
            # they are the ones that bear on Sec. 14.2's claim
            if key == 'gbp':
                w = [k for k, r in enumerate(rows)
                     if abs(r['gbp']) >= abs(r['kikuchi'])]
                if w and mk == 's':
                    ax.semilogy(x[w], y[w], mk, ms=3.6, mew=0.9,
                                mfc=DARK, color=DARK)
        ax.annotate(name, xy=(i * 3.4 + 1, 3e2), fontsize=6.8, ha='center',
                    va='bottom', color='0.3', annotation_clip=False)
    ax.set_xticks([i * 3.4 + j for i in range(3) for j in range(3)])
    ax.set_xticklabels(['B', 'M', 'G'] * 3, fontsize=7)
    ax.set_xlim(-0.8, 9.2)
    ax.set_ylim(3e-15, 1e4)
    ax.set_ylabel(r'error in $\ln Z$', fontsize=8.5)
    ax.set_xlabel('B: Bethe counting    M: M\u00f6bius    G: GBP',
                  fontsize=7.5)
    _tidy(ax)
    fig.tight_layout()
    fig.savefig(OUT / 'fig-gbp.pdf')
    print(f'  wrote {OUT / "fig-gbp.pdf"}')


# ------------------------------------------- (2b) real clustered neighbourhoods
# A fixed point counts only when the messages have stopped moving AND the
# beliefs agree where they overlap.  Fig. 14.2's hyperbolic instances never
# forced the distinction -- every residual-converged run there was consistent to
# 2e-6 -- but on real neighbourhoods the two part company, so both are checked.
SOUND_RES, SOUND_CONS = 1e-9, 1e-6


def _sound(r):
    return (r['residual'] < SOUND_RES and r['consistency'] is not None
            and r['consistency'] < SOUND_CONS)


def _real():
    return json.load(open(PROBE / 'gbp_real.json'))


def _by_network(d):
    """Per network: family, clustering, mean neighbourhood density, and the
    chordal fraction.  Counted over ego-networks rather than runs, the two
    couplings on one neighbourhood being one structure asked twice."""
    out = {}
    for r in d:
        e = out.setdefault(r['network'], {'family': r['family'],
                                          'C': r['clustering'], 'ego': {}})
        e['ego'][r['ego']] = (r['chordal'],
                              2 * r['m'] / (r['n'] * (r['n'] - 1)))
    for e in out.values():
        v = list(e['ego'].values())
        e['n_ego'] = len(v)
        e['frac'] = sum(c for c, _ in v) / len(v)
        e['density'] = float(np.mean([rho for _, rho in v]))
    return out


def check_gbp_real():
    d = _real()
    nets = _by_network(d)
    print(f'  {len(d)} runs over {sum(e["n_ego"] for e in nets.values())} '
          f'ego-networks of {len(nets)} real networks, 8 <= n <= 20')
    for fam in ('grouped', 'dyadic'):
        ego = {(r['network'], r['ego']): r['chordal']
               for r in d if r['family'] == fam}
        print(f'    {fam:>8}: {sum(ego.values())}/{len(ego)} ego-networks '
              f'chordal ({100 * sum(ego.values()) / len(ego):.0f}%)')
        for name, e in sorted(nets.items(), key=lambda kv: -kv[1]['frac']):
            if e['family'] == fam:
                print(f'        {name:<18} C = {e["C"]:.2f}  '
                      f'density {e["density"]:.2f}  '
                      f'{int(round(e["frac"] * e["n_ego"])):>2d}/{e["n_ego"]:<2d} '
                      f'chordal')
    # the confound: sparser neighbourhoods are chordal more easily, so the split
    # has to be shown not to be a density split, or a clustering one, in disguise
    gr = [e for e in nets.values() if e['family'] == 'grouped']
    dy = [e for e in nets.values() if e['family'] == 'dyadic']
    print(f'    the two families overlap in both controls: C '
          f'{min(e["C"] for e in gr):.2f}-{max(e["C"] for e in gr):.2f} against '
          f'{min(e["C"] for e in dy):.2f}-{max(e["C"] for e in dy):.2f}, density '
          f'{min(e["density"] for e in gr):.2f}-{max(e["density"] for e in gr):.2f} '
          f'against {min(e["density"] for e in dy):.2f}-'
          f'{max(e["density"] for e in dy):.2f}')

    ch = [r for r in d if r['chordal']]
    worst = max(abs(r['gbp']) for r in ch)
    print(f'    chordal runs {len(ch):>3}:  exact, worst error {worst:.1e}')
    assert worst < 1e-10, worst

    nc = [r for r in d if not r['chordal']]
    lowres = [r for r in nc if r['residual'] < SOUND_RES]
    snd = [r for r in nc if _sound(r)]
    print(f'    non-chordal  {len(nc):>3}:  {len(lowres)} reach a residual '
          f'below {SOUND_RES:g}, of which {len(snd)} are also consistent')
    print(f'      {len(lowres) - len(snd)} stop moving without the beliefs '
          f'agreeing: on real neighbourhoods a small residual\n'
          f'      no longer certifies a usable fixed point, which it did on '
          f'every hyperbolic instance')
    for fam in ('grouped', 'dyadic'):
        sel = [r for r in snd if r['family'] == fam]
        g = [abs(r['gbp']) for r in sel]
        orders = [np.log10(abs(r['kikuchi']) / max(abs(r['gbp']), 1e-16))
                  for r in sel]
        worse = [r for r in sel if abs(r['gbp']) >= abs(r['kikuchi'])]
        print(f'      {fam:>8}: {len(sel):>2} sound of '
              f'{sum(1 for r in nc if r["family"] == fam):>2};  GBP '
              f'{min(g):.1e} to {max(g):.1e};  gain {min(orders):+.1f} to '
              f'{max(orders):+.1f} orders;  no better than static: {len(worse)}')
    bad = [r for r in snd if abs(r['gbp']) >= abs(r['kikuchi'])]
    print(f'    {len(bad)} sound fixed points are no better than the static '
          f'M\u00f6bius estimate they were to improve:')
    for r in sorted(bad, key=lambda r: -abs(r['gbp'])):
        print(f'      {r["network"]:<18} ego {r["ego"]:>4}, n = {r["n"]:>2}, '
              f'{r["n_cliques"]:>2} maximal cliques, bJ = {r["beta_J"]}: '
              f'GBP {r["gbp"]:+.3f} against M\u00f6bius {r["kikuchi"]:+.3f} '
              f'(residual {r["residual"]:.0e}, consistency '
              f'{r["consistency"]:.0e})')
    print('      maximal cliques never nest, so nesting is not what causes it')


def figure_gbp_real():
    plt = _mpl()
    d = _real()
    fig, (ax, bx) = plt.subplots(1, 2, figsize=(6.4, 2.9),
                                 gridspec_kw={'width_ratios': [1.5, 1]})

    # left: Fig. 14.2's axes, on real neighbourhoods
    nc = [r for r in d if not r['chordal']]
    groups = (('chordal', [r for r in d if r['chordal']], 'o', DARK),
              ('non-chordal,\nsound fixed point',
               [r for r in nc if _sound(r)], 's', MID),
              ('non-chordal,\nnone',
               [r for r in nc if not _sound(r)], 'x', LIGHT))
    floor = 1e-14
    rng = np.random.default_rng(0)
    for i, (name, rows, mk, col) in enumerate(groups):
        for j, key in enumerate(('bethe', 'kikuchi', 'gbp')):
            x = i * 3.4 + j + rng.uniform(-0.16, 0.16, len(rows))
            y = np.maximum([abs(r[key]) for r in rows], floor)
            ax.semilogy(x, y, mk, ms=3.2, mew=0.9,
                        mfc='white' if mk != 'x' else col, color=col)
            # the four sound fixed points that came out no better than the
            # static estimate they were supposed to improve: filled, because
            # they are the ones that bear on Sec. 14.2's claim
            if key == 'gbp':
                w = [k for k, r in enumerate(rows)
                     if abs(r['gbp']) >= abs(r['kikuchi'])]
                if w and mk == 's':
                    ax.semilogy(x[w], y[w], mk, ms=3.6, mew=0.9,
                                mfc=DARK, color=DARK)
        ax.annotate(name, xy=(i * 3.4 + 1, 3e2), fontsize=6.8, ha='center',
                    va='bottom', color='0.3', annotation_clip=False)
    ax.set_xticks([i * 3.4 + j for i in range(3) for j in range(3)])
    ax.set_xticklabels(['B', 'M', 'G'] * 3, fontsize=7)
    ax.set_xlim(-0.8, 9.2)
    ax.set_ylim(3e-15, 1e4)
    ax.set_ylabel(r'error in $\ln Z$', fontsize=8.5)
    ax.set_xlabel('B: Bethe counting    M: M\u00f6bius    G: GBP', fontsize=7.5)
    _tidy(ax)

    # right: which networks land in the exact group, and why
    nets = _by_network(d)
    order = sorted(nets.items(), key=lambda kv: (kv[1]['family'] != 'grouped',
                                                 -kv[1]['frac']))
    y = np.arange(len(order))[::-1]
    for yi, (name, e) in zip(y, order):
        col = DARK if e['family'] == 'grouped' else LIGHT
        bx.barh(yi, e['frac'], height=0.62, color=col, edgecolor='none')
    bx.set_yticks(y)
    bx.set_yticklabels([f'{n}  ($C\\!=\\!{e["C"]:.2f}$)' for n, e in order],
                       fontsize=6.6)
    bx.set_xlim(0, 1)
    bx.set_xticks([0, 0.5, 1])
    bx.set_xlabel('neighbourhoods that are chordal', fontsize=7.5)
    bx.tick_params(labelsize=7)
    for sp in ('top', 'right'):
        bx.spines[sp].set_visible(False)
    bx.set_ylim(-0.7, len(order) - 0.1)
    bx.annotate('group-generated', xy=(0.98, y[0] + 0.52), fontsize=6.5,
                color=DARK, ha='right', va='center')
    n_grouped = sum(1 for _, e in order if e['family'] == 'grouped')
    if n_grouped < len(order):
        bx.axhline(y[n_grouped] + 0.5, color='0.85', lw=0.6)
        bx.annotate('dyadic', xy=(0.98, y[n_grouped] + 0.52), fontsize=6.5,
                    color='0.45', ha='right', va='center')

    fig.tight_layout()
    fig.savefig(OUT / 'fig-gbp-real.pdf')
    print(f'  wrote {OUT / "fig-gbp-real.pdf"}')


# ------------------------------------------- (2c) the cavity failure, Ch. 14
def _karrer():
    f = PROBE / 'gbp_karrer.json'
    return json.load(open(f))['runs'] if f.exists() else []


def _cavity():
    f = PROBE / 'cavity_clique.json'
    return json.load(open(f)) if f.exists() else []


def check_cavity():
    """The cavity method itself, on the clique chygraph of each class."""
    d = _cavity()
    if not d:
        print('  (probe/results/cavity_clique.json not present)')
        return
    print('  error in ln Z of the cavity method on the clique chygraph:')
    for ens in ('hyperbolic', 'karrer', 'real'):
        sel = [r for r in d if r['ensemble'] == ens]
        conv = [r for r in sel if r['residual'] < 1e-9]
        e = [abs(r['cavity']) for r in conv]
        fr = [r['doubled_bonds'] / r['n_bonds'] for r in sel]
        rho = np.corrcoef(np.log10(np.maximum(e, 1e-16)),
                          [r['doubled_bonds'] / r['n_bonds'] for r in conv])[0, 1]
        print(f'    {ens:<12} {len(sel):>3} runs, {len(conv)} converged   '
              f'{min(e):.2e} to {max(e):.2e}   median {np.median(e):.2f}   '
              f'corr {rho:+.2f}')
        print(f'      bonds inside two or more complexes: '
              f'{100 * min(fr):.0f} to {100 * max(fr):.0f} per cent')
    assert all(r['residual'] < 1e-9 for r in d), 'a cavity run did not converge'
    assert all(r['doubled_bonds'] > 0 for r in d), 'an instance had no doubled bond'
    print('    every run reaches a fixed point, and every instance has at least'
          ' one\n    bond lying inside two complexes --- which is the defect'
          ' being priced')


def figure_cavity():
    """Ch. 14: the cavity error against how much of the graph is double-counted."""
    plt = _mpl()
    d = _cavity()
    if not d:
        return
    fig, ax = plt.subplots(figsize=(4.3, 2.9))
    for name, ens, mk, col in (('hyperbolic', 'hyperbolic', 'o', DARK),
                               ('Karrer--Newman', 'karrer', '^', MID),
                               ('real', 'real', 's', LIGHT)):
        sel = [r for r in d if r['ensemble'] == ens]
        if not sel:
            continue
        x = [r['doubled_bonds'] / r['n_bonds'] for r in sel]
        y = np.maximum([abs(r['cavity']) for r in sel], 1e-3)
        ax.semilogy(x, y, mk, ms=3.4, mew=0.9, mfc='white', color=col,
                    label=f'{name} ({len(sel)})')
    ax.axhline(1.0, color='0.75', lw=0.7, ls=':')   # labelled in the caption
    ax.set_xlabel('fraction of bonds lying inside two or more complexes',
                  fontsize=8)
    ax.set_ylabel(r'error in $\ln Z$, cavity method', fontsize=8.5)
    ax.set_xlim(-0.03, 1.03)
    ax.legend(fontsize=6.6, frameon=False, loc='upper left')
    _tidy(ax)
    fig.tight_layout()
    fig.savefig(OUT / 'fig-cavity.pdf')
    print(f'  wrote {OUT / "fig-cavity.pdf"}')


# ------------------------------- (2c-bis) the pairwise condition is not enough
def ring_of_triangles(k):
    """k triangles glued in a ring, each meeting the next at a single vertex."""
    import networkx as nx
    G = nx.Graph()
    hub = [3 * i for i in range(k)]
    for i in range(k):
        a, b, c = hub[i], 3 * i + 1, hub[(i + 1) % k]
        G.add_edges_from([(a, b), (b, c), (a, c)])
    return nx.convert_node_labels_to_integers(G)


def check_ring():
    """Pairwise-treelike, cyclic incidence, and not exact.

    Section 2.9 states the condition as: no two complexes share two atoms.  That
    is necessary and not sufficient.  A loop can run through three or more
    complexes each meeting the next in a single node, and the recursion is then
    wrong while the pairwise test passes.
    """
    import sys
    from itertools import combinations
    import networkx as nx
    sys.path.insert(0, str(Path.home() / 'av2atg' / 'chygraph_statmech'
                           / 'probe'))
    from cavity_clique import ChygraphBP
    from chygraph_statmech.gbp import exact_log_Z, ising_factors

    print('  rings of triangles glued at single vertices:')
    print(f'    {"k":>3}{"n":>4}{"pairwise":>10}{"forest":>8}'
          f'{"bJ=0.3":>11}{"bJ=0.8":>11}')
    rows = []
    for k in (4, 5, 6, 8):
        G = ring_of_triangles(k)
        n = G.number_of_nodes()
        cx = [sorted(c) for c in nx.find_cliques(G)]
        pw = all(len(set(a) & set(b)) <= 1 for a, b in combinations(cx, 2))
        B = nx.Graph()
        for i, c in enumerate(cx):
            for v in c:
                B.add_edge(('c', i), ('v', v))
        forest = nx.is_forest(B)
        edges = sorted({tuple(sorted(e)) for e in G.edges()})
        errs = []
        for bJ in (0.3, 0.8):
            ex = exact_log_Z(ising_factors(edges, bJ), range(n))
            errs.append(ChygraphBP(cx, bJ, damping=0.5,
                                   edges=edges).run().log_Z() - ex)
        rows.append((k, n, pw, forest, errs))
        print(f'    {k:>3}{n:>4}{str(pw):>10}{str(forest):>8}'
              f'{errs[0]:>+11.4f}{errs[1]:>+11.4f}')
        assert pw and not forest, k
        assert abs(errs[0]) > 1e-4 and abs(errs[1]) > 1e-2, k
    print('    every one passes the pairwise test of Sec. 2.9 and none is a')
    print('    forest, and the recursion is wrong on all of them: the pairwise')
    print('    condition is necessary and not sufficient')
    print('    the error decays as the ring lengthens -- it is a short-loop'
          ' effect')
    return rows


# ------------------------------------------- (2c-ter) the leaf-removal core
def _core():
    f = PROBE / 'core_fraction.json'
    if not f.exists():
        return []
    d = json.load(open(f))
    return d['instances'] if isinstance(d, dict) else d


def _core_examples():
    f = PROBE / 'core_fraction.json'
    if not f.exists():
        return {}
    d = json.load(open(f))
    return d.get('karrer_examples', {}) if isinstance(d, dict) else {}


def check_core_fraction():
    """Atoms the propagation cannot strip, on the two chygraphs it runs on."""
    d = _core()
    if not d:
        print('  (probe/results/core_fraction.json not present)')
        return
    print('  fraction of atoms the messages cannot remove by propagation:')
    print(f'    {"":12}{"":>5}{"clique chygraph":>22}{"merged chygraph":>22}')
    for ens in ('hyperbolic', 'karrer', 'real'):
        sel = [r for r in d if r['ensemble'] == ens]
        c = [r['clique_core_frac'] for r in sel]
        m = [r['merged_core_frac'] for r in sel]
        print(f'    {ens:<12}{len(sel):>5}   median {np.median(c):.2f}, '
              f'{sum(1 for x in c if x < 1e-9):>2} core-free'
              f'   median {np.median(m):.2f}, '
              f'{sum(1 for x in m if x < 1e-9):>2} core-free')
    assert all(r['clique_core_frac'] > 0 for r in d), \
        'some clique chygraph was already core-free'
    print('    not one clique chygraph is core-free: on every instance the')
    print('    propagation is left with a core it cannot resolve, which is why')
    print('    Ch. 14 finds an error everywhere')
    z = [r for r in d if r['merged_core_frac'] < 1e-9]
    print(f'    merging drives the core to zero on {len(z)} of {len(d)} '
          f'instances, and those are')
    print('    exactly the ones Sec. 16.2 finds exact')


CLASSES = (('hyperbolic', 'hyperbolic', 'o', DARK),
           ('Karrer--Newman', 'karrer', '^', MID),
           ('real', 'real', 's', LIGHT))


def _by_class(ax, rows, key, floor=None):
    """One column per ensemble, jittered; the classes are the comparison."""
    rng = np.random.default_rng(0)
    for i, (name, ens, mk, col) in enumerate(CLASSES):
        sel = [r for r in rows if r['ensemble'] == ens]
        if not sel:
            continue
        x = i + rng.uniform(-0.22, 0.22, len(sel))
        y = [abs(r[key]) for r in sel]
        if floor is not None:
            y = np.maximum(y, floor)
        ax.plot(x, y, mk, ms=3.2, mew=0.9, mfc='white', color=col)
    ax.set_xticks(range(len(CLASSES)))
    ax.set_xticklabels([n for n, _, _, _ in CLASSES], fontsize=7.5)
    ax.set_xlim(-0.6, len(CLASSES) - 0.4)
    _tidy(ax)


def figure_core_meta():
    """Ch. 16: the core the merged chygraph leaves, by ensemble."""
    plt = _mpl()
    d = _core()
    if not d:
        return
    fig, ax = plt.subplots(figsize=(4.3, 2.3))
    _by_class(ax, d, 'merged_core_frac')
    ax.set_ylabel('core of the\nmerged chygraph', fontsize=8.5)
    ax.set_ylim(-0.05, 1.05)
    ax.set_yticks([0, 0.5, 1])
    fig.tight_layout()
    fig.savefig(OUT / 'fig-core-meta.pdf')
    print(f'  wrote {OUT / "fig-core-meta.pdf"}')


def figure_core_bonds():
    """Ch. 14: the core the clique chygraph leaves, by ensemble."""
    plt = _mpl()
    d = _core()
    if not d:
        return
    fig, ax = plt.subplots(figsize=(4.3, 2.3))
    _by_class(ax, d, 'clique_core_frac')
    ax.set_ylabel('core of the\nclique chygraph', fontsize=8.5)
    ax.set_ylim(-0.05, 1.05)
    ax.set_yticks([0, 0.5, 1])
    fig.tight_layout()
    fig.savefig(OUT / 'fig-core-bonds.pdf')
    print(f'  wrote {OUT / "fig-core-bonds.pdf"}')


def figure_core_example():
    """Ch. 16: what the meta-complex messages are actually left with.

    Two measured cores, drawn as inclusion diagrams in Ch. 2's vocabulary:
    filled circles are atoms, open squares complexes, a line is membership.
    Atoms leaf removal stripped are shown faint, so the reader sees what went.
    """
    plt = _mpl()
    ex = _core_examples()
    if not ex:
        return
    import networkx as nx
    fig, axes = plt.subplots(1, 2, figsize=(6.4, 2.9),
                             gridspec_kw={'width_ratios': [1, 1.35]})
    for ax, key, title in ((axes[0], 'smallest', 'the modal case'),
                           (axes[1], 'largest', 'the largest in the sample')):
        d = ex[key]
        core = set(d['atoms'])
        B = nx.Graph()
        for i, c in enumerate(d['complexes']):
            for v in c:
                B.add_edge(('c', i), ('v', v))
        K = B.subgraph([n for n in B if n[0] == 'c' or n[1] in core])
        pos = nx.kamada_kawai_layout(K)
        pos = nx.spring_layout(B, pos=dict(pos), fixed=list(K), seed=3,
                               k=0.55) if B.number_of_nodes() > K.number_of_nodes() else dict(pos)
        for u, v in B.edges():
            inside = all(w[0] == 'c' or w[1] in core for w in (u, v))
            ax.plot(*zip(pos[u], pos[v]), '-', lw=1.0 if inside else 0.6,
                    color='0.35' if inside else '0.85', zorder=1)
        for n in B:
            x, y = pos[n]
            if n[0] == 'c':
                ax.plot(x, y, 's', ms=6.5, mfc='white', mec='0.25', mew=1.0,
                        zorder=3)
            elif n[1] in core:
                ax.plot(x, y, 'o', ms=5.5, color='0.15', zorder=3)
            else:
                ax.plot(x, y, 'o', ms=4.0, mfc='white', mec='0.75', mew=0.8,
                        zorder=2)
        ax.set_title(f'{title}: {len(core)} atoms, '
                     f'{len(d["complexes"])} complexes,\n'
                     f'{d["cycles"]} independent '
                     f'{"cycle" if d["cycles"] == 1 else "cycles"}, '
                     f'{d["components"]} component',
                     fontsize=7.2, color='0.25', pad=6)
        ax.set_aspect('equal')
        ax.axis('off')
    fig.tight_layout()
    fig.savefig(OUT / 'fig-core-example.pdf')
    print(f'  wrote {OUT / "fig-core-example.pdf"}')


# ------------------------------------- (2d) the meta-complex error, Ch. 16
def _mergelnz():
    f = PROBE / 'merge_lnz.json'
    return json.load(open(f)) if f.exists() else []


def check_merge_error():
    """Is the merged chygraph exact?  Only when its incidence structure is."""
    d = _mergelnz()
    if not d:
        print('  (probe/results/merge_lnz.json not present)')
        return
    print(f'  {len(d)} runs; exactness follows the incidence structure, not the'
          ' pairwise condition:')
    for ens in ('hyperbolic', 'karrer', 'real'):
        sel = [r for r in d if r['ensemble'] == ens]
        ac = sum(r['acyclic'] for r in sel)
        e = [abs(r['gbp']) for r in sel]
        ex = sum(1 for x in e if x < 1e-9)
        print(f'    {ens:<12} {len(sel):>3} runs, treelike '
              f'{sum(r["treelike"] for r in sel)}, acyclic {ac}, exact {ex}, '
              f'worst error {max(e):.2e}')
    ac = [r for r in d if r['acyclic']]
    cy = [r for r in d if not r['acyclic']]
    assert all(abs(r['gbp']) < 1e-9 for r in ac), 'an acyclic run was not exact'
    assert all(abs(r['gbp']) > 1e-9 for r in cy), 'a cyclic run came out exact'
    print(f'    every one of the {len(ac)} acyclic runs is exact; every one of '
          f'the {len(cy)} cyclic runs is not')
    if cy:
        e = [abs(r['gbp']) for r in cy]
        print(f'      cyclic error {min(e):.1e} to {max(e):.1e}, '
              f'median {np.median(e):.2f}')
    thr = max(r['n_meta'] for r in d if r['acyclic'] and
              all(x['acyclic'] for x in d if x['n_meta'] <= r['n_meta']))
    print(f'    every family of {thr} meta-complexes or fewer is acyclic; '
          f'cycles start at {thr + 1}')


def figure_merge_error():
    """Ch. 16: the error the merged chygraph leaves, by ensemble."""
    plt = _mpl()
    d = _mergelnz()
    if not d:
        return
    fig, ax = plt.subplots(figsize=(4.3, 2.6))
    _by_class(ax, d, 'gbp', floor=3e-15)
    ax.set_yscale('log')
    ax.set_ylabel(r'error in $\ln Z$', fontsize=8.5)
    ax.set_ylim(1e-15, 1e1)
    for i, (_, ens, _, _) in enumerate(CLASSES):
        sel = [r for r in d if r['ensemble'] == ens]
        ex = sum(1 for r in sel if abs(r['gbp']) < 1e-9)
        ax.annotate(f'{ex}/{len(sel)}\nexact', xy=(i, 3e1), fontsize=6.6,
                    ha='center', va='bottom', color='0.35',
                    annotation_clip=False)
    fig.tight_layout()
    fig.savefig(OUT / 'fig-merge-error.pdf')
    print(f'  wrote {OUT / "fig-merge-error.pdf"}')


# ------------------------------- (2e) where the core appears, Ch. 16
def _sweep():
    f = PROBE / 'karrer_core_sweep.json'
    return json.load(open(f)) if f.exists() else []


def check_core_transition():
    """The core appears at incidence branching one, whatever produces it."""
    d = _sweep()
    if not d:
        print('  (probe/results/karrer_core_sweep.json not present)')
        return
    fams = sorted({r['family'] for r in d})
    nmax = max(r['n'] for r in d)
    print(f'  merged-chygraph core against b = s + 2t, at n = {nmax}:')
    print(f'    {"b":>5}' + ''.join(f'{f.split(",")[0]:>22}' for f in fams))
    for b in (0.4, 0.8, 1.0, 1.2, 1.6, 2.0):
        row = f'    {b:>5.1f}'
        for f in fams:
            v = [r for r in d if r['family'] == f and r['n'] == nmax
                 and abs(r['b'] - b) < 1e-9]
            row += f'{(v[0]["core_mean"] if v else float("nan")):>22.4f}'
        print(row)
    for f in fams:
        below = [r['core_mean'] for r in d
                 if r['family'] == f and r['n'] == nmax and r['b'] <= 0.8]
        assert max(below) < 0.01, (f, max(below))
    print('    every family is core-free below b = 1 and orders above it;')
    print('    at t = 0 the complexes are bare edges and this is the')
    print('    Erdos-Renyi 2-core at mean degree one, which fixes the scale')


def figure_core_transition():
    """Ch. 16: the core as an order parameter in the branching ratio."""
    plt = _mpl()
    d = _sweep()
    if not d:
        return
    fams = [('edges only, $t=0$', 'o', DARK),
            ('triangles only, $s=0$', '^', MID),
            ('equal mixture', 's', LIGHT)]
    sizes = sorted({r['n'] for r in d})
    fig, (ax, bx) = plt.subplots(1, 2, figsize=(6.4, 2.5))
    nmax = sizes[-1]
    for label, mk, col in fams:
        sel = sorted([r for r in d if r['family'] == label and r['n'] == nmax],
                     key=lambda r: r['b'])
        ax.plot([r['b'] for r in sel], [r['core_mean'] for r in sel],
                mk + '-', ms=3.0, lw=1.0, mfc='white', color=col,
                label=label.replace('$', '$'))
    ax.axvline(1.0, color='0.75', lw=0.7, ls=':')
    ax.annotate('$b=1$', xy=(1.06, 0.02), fontsize=6.8, color='0.45')
    ax.set_xlabel('incidence branching $b=s+2t$', fontsize=8)
    ax.set_ylabel('core of the\nmerged chygraph', fontsize=8.5)
    ax.set_ylim(-0.03, 1.0)
    ax.legend(fontsize=6.4, frameon=False, loc='upper left')
    _tidy(ax)

    lab = 'triangles only, $s=0$'
    for n, mk in zip(sizes, ('o', 's', '^', 'D')):
        sel = sorted([r for r in d if r['family'] == lab and r['n'] == n],
                     key=lambda r: r['b'])
        bx.plot([r['b'] for r in sel], [r['core_mean'] for r in sel],
                mk + '-', ms=2.8, lw=0.9, mfc='white',
                color=str(0.10 + 0.18 * sizes.index(n)), label=f'$n={n}$')
    bx.axvline(1.0, color='0.75', lw=0.7, ls=':')
    bx.set_xlabel('incidence branching $b=s+2t$', fontsize=8)
    bx.set_ylabel('core, triangles only', fontsize=8.5)
    bx.set_xlim(0.3, 1.9)
    bx.set_ylim(-0.02, 0.55)
    bx.legend(fontsize=6.4, frameon=False, loc='upper left')
    _tidy(bx)
    fig.tight_layout()
    fig.savefig(OUT / 'fig-core-transition.pdf')
    print(f'  wrote {OUT / "fig-core-transition.pdf"}')


# ------------------------------------------- (3) does the clique ensemble exist
def check_ensemble():
    """The paired ratio of excess cardinality, from the cached scan."""
    txt = (PROBE / 'analysis.txt').read_text().splitlines()
    i = next(j for j, L in enumerate(txt) if 'Paired ratio' in L)
    print('  paired ratio sbar(HRG)/sbar(control), identical degree sequence:')
    for line in txt[i + 2:i + 11]:
        f = line.split()
        if len(f) == 9:
            print(f'    tau = {f[0]}, <k> = {f[1]}:  {f[2]} -> {f[7]} '
                  f'over n = 10^3 to 3x10^5,  growth {f[8]}')
    print('    for tau >= 2.5 the ratio is a converged constant between 1.8 '
          'and 4.8:')
    print('    the clustering signal is a multiplicative factor, not a '
          'finite-size artefact.')
    print('    At tau = 2.1 both diverge and the control diverges faster, so '
          'nothing')
    print('    measured there separates clustering from the tail.')


if __name__ == '__main__':
    print('the two countings:')
    check_counting()
    print('two triangles, exactly:')
    check_two_triangles()
    print('the cavity failure, across the three ensembles:')
    check_cavity()
    print('sixty clique region graphs:')
    check_gbp_instances()
    print('real clustered neighbourhoods:')
    check_gbp_real()
    print('the pairwise condition is not sufficient:')
    check_ring()
    print('the leaf-removal core of each instance:')
    check_core_fraction()
    print('where the core appears:')
    check_core_transition()
    print('the merged chygraph, and what it gets wrong:')
    check_merge_error()
    print('when the clique ensemble exists:')
    check_ensemble()
    print('figures:')
    figure_cavity()
    figure_gbp()
    figure_gbp_real()
    figure_merge_error()
    figure_core_bonds()
    figure_core_meta()
    figure_core_example()
    figure_core_transition()
