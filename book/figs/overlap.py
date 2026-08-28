"""Chapter 13: what overlapping complexes cost, and what recovers it.

  fig-gbp        the same three schemes over 60 maximal-clique region graphs of
                 hyperbolic random graphs, split by whether the clique
                 structure is chordal

Chapter 13's first figure, two triangles sharing an edge, is TikZ in
`overlap.tex`; its table of two-triangle errors is printed by
`check_two_triangles` below.

The 60-instance data is read from `../probe/results/gbp_cliques.json`, the
cached output of `../probe/gbp_cliques.py`. The two-triangle numbers are
recomputed here, being a four-spin enumeration.
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
    print('sixty clique region graphs:')
    check_gbp_instances()
    print('when the clique ensemble exists:')
    check_ensemble()
    print('figure:')
    figure_gbp()
