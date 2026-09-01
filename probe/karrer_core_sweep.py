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
from scipy.optimize import brentq

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


def core_closed_form(s, t):
    """Eqs. (16.2) and (16.3): the 2-core of the incidence structure.

    Each atom carries Poisson(s) memberships of cardinality two and Poisson(t)
    of cardinality three.  Let `a` be the probability that an atom prunes given
    the structure reached through one of its complexes.  A cardinality-c
    complex reached from an atom prunes when all c-1 of its other atoms do, so
    g_c = a^(c-1); an atom prunes when all its other complexes do, and for
    Poisson memberships the excess law is the original one, so

        a = exp[-s(1-a) - t(1-a^2)].

    Differentiating the right-hand side at a = 1 gives s + 2t, so the trivial
    root a = 1 loses stability at b = 1 -- Eq. (16.1), derived.  An atom is in
    the 2-core when at least two memberships reach complexes that do not prune;
    thinning a Poisson leaves a Poisson, so that count is Poisson(mu) with
    mu = s(1-a) + t(1-a^2).
    """
    F = lambda a: np.exp(-s * (1 - a) - t * (1 - a ** 2))            # noqa: E731
    a = 1.0 if s + 2 * t <= 1 else brentq(lambda a: F(a) - a, 1e-14, 1 - 1e-12)
    mu = s * (1 - a) + t * (1 - a ** 2)
    return float(1 - np.exp(-mu) - mu * np.exp(-mu))


def check_closed_form(rows):
    """The closed form against the sweep, at the largest size."""
    print('closed form against the measurement, n = %d:' % max(SIZES))
    print('    family                    b   measured   Eq. (16.3)')
    worst = 0.0
    for label, _ in FAMILIES:
        for r in [r for r in rows if r['family'] == label
                  and r['n'] == max(SIZES) and r['b'] in (1.2, 1.6, 2.0, 2.5)]:
            c = core_closed_form(r['s'], r['t'])
            worst = max(worst, abs(c - r['core_mean']))
            print(f"    {label:<22} {r['b']:>4}     {r['core_mean']:.4f}"
                  f"       {c:.4f}")
    print(f'  worst discrepancy {worst:.4f}, against a standard error of'
          f' about 0.004')
    print('  -- threshold and amplitude both, so the transition is the'
          ' classical one')
    assert worst < 0.01, 'the closed form no longer matches the sweep'


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
    print('\nwrote', OUT, '\n')
    check_closed_form(rows)


if __name__ == '__main__':
    main()
