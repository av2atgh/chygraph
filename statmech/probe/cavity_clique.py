"""The cavity method itself, run on the clique chygraph.

Chapter 14 needs the error of the calculation the book actually performs, which
is message passing on the chygraph's incidence structure -- not the counting
evaluated on isolated regions, which is a cruder object and wrong even where the
chygraph is treelike.  This runs the real thing.

The chygraph is bipartite: atoms on one side, complexes on the other.  A complex
sums over its interior exactly, so its factor is

    f_A(sigma_A) = exp( beta J sum_{(i,j) in A} sigma_i sigma_j ),

every pair inside A included.  That is the definition, and it is also where the
damage comes from: a bond lying inside two complexes is summed by both, so its
interaction is counted twice.  The messages are the standard factor-graph pair

    m_{v->A}(s) ~ prod_{B ni v, B != A} m_{B->v}(s),
    m_{A->v}(s) ~ sum_{sigma_{A\\v}} f_A(sigma_A) prod_{u in A, u != v} m_{u->A},

and ln Z is the Bethe free energy read off the fixed point,

    ln Z = sum_A ln Z_A + sum_v (1 - k_v) ln Z_v.

Where the incidence structure is a forest and no bond lies in two complexes this
is exact, which is the check `validate` performs before any measurement is
quoted.

    python probe/cavity_clique.py
"""

import json
import sys
from itertools import combinations
from pathlib import Path

import networkx as nx
import numpy as np
from scipy.special import logsumexp

sys.path.insert(0, str(Path.home() / 'av2atg' / 'chygraph' / 'src'))
sys.path.insert(0, str(Path.home() / 'av2atg' / 'statmech'
                      / 'book' / 'figs'))
sys.path.insert(0, str(Path.home() / 'av2atg' / 'statmech' / 'probe'))

from statmech.gbp import exact_log_Z, ising_factors  # noqa: E402
from statmech.region import overlap_profile  # noqa: E402

RESULTS = Path(__file__).parent / 'results'
OUT = RESULTS / 'cavity_clique.json'
COUPLINGS = (0.3, 0.8)
DAMPING = (0.0, 0.5, 0.9, 0.97, 0.995)
TOL = 1e-12


def _norm(v):
    return v - logsumexp(v)


class ChygraphBP:
    """Belief propagation on the atom--complex incidence structure."""

    def __init__(self, complexes, bJ, damping=0.5, edges=None):
        """`edges` is the bond set; without it every pair inside a complex is
        taken to be one, which holds for maximal cliques and fails for anything
        built out of them -- a merged meta-complex is a union of cliques, not a
        clique, and assuming otherwise invents bonds that are not there."""
        self.A = [tuple(sorted(c)) for c in complexes]
        self.bonds = None if edges is None else {
            tuple(sorted(e)) for e in edges}
        self.bJ = float(bJ)
        self.damping = float(damping)
        self.nodes = sorted({v for c in self.A for v in c})
        self.k = {v: sum(1 for c in self.A if v in c) for v in self.nodes}
        self.logf = [self._factor(c) for c in self.A]
        self.m_va = {(v, a): np.zeros(2) for a, c in enumerate(self.A)
                     for v in c}
        self.m_av = {(a, v): np.zeros(2) for a, c in enumerate(self.A)
                     for v in c}
        self.residual = np.inf

    def _factor(self, c):
        """exp(bJ sum_{bonds inside c} s_i s_j) as a log-table over 2^|c|.

        A bond lying inside two complexes is summed by both.  That is the
        double count the chapter is about, and it is deliberate here.
        """
        n = len(c)
        s = np.array([[1 if (i >> b) & 1 == 0 else -1 for b in range(n)]
                      for i in range(2 ** n)], dtype=float)
        e = np.zeros(2 ** n)
        for p, q in combinations(range(n), 2):
            if self.bonds is None or (c[p], c[q]) in self.bonds:
                e += s[:, p] * s[:, q]
        return (self.bJ * e).reshape((2,) * n)

    def _sweep(self):
        d, delta = self.damping, 0.0
        for a, c in enumerate(self.A):
            for j, v in enumerate(c):
                acc = self.logf[a].copy()
                for i, u in enumerate(c):
                    if u == v:
                        continue
                    sh = [1] * len(c)
                    sh[i] = 2
                    acc = acc + self.m_va[(u, a)].reshape(sh)
                ax = tuple(i for i in range(len(c)) if i != j)
                new = _norm(logsumexp(acc, axis=ax) if ax else acc)
                delta = max(delta, float(np.abs(new - self.m_av[(a, v)]).max()))
                self.m_av[(a, v)] = d * self.m_av[(a, v)] + (1 - d) * new
        for a, c in enumerate(self.A):
            for v in c:
                acc = np.zeros(2)
                for b, cb in enumerate(self.A):
                    if b != a and v in cb:
                        acc = acc + self.m_av[(b, v)]
                new = _norm(acc)
                delta = max(delta, float(np.abs(new - self.m_va[(v, a)]).max()))
                self.m_va[(v, a)] = d * self.m_va[(v, a)] + (1 - d) * new
        return delta

    def run(self, sweeps=5000):
        for _ in range(sweeps):
            self.residual = self._sweep()
            if self.residual < TOL:
                break
        return self

    def log_Z(self):
        """Bethe free energy of the incidence structure at the fixed point.

        ln Z = sum_A ln Z_A + sum_v ln Z_v - sum_{(A,v)} ln Z_{Av}.

        The edge term is not optional.  Writing the node contribution as
        (1 - k_v) ln Z_v instead absorbs it, and that form is only valid when
        every Z_{Av} equals Z_v, which normalised messages do not satisfy.
        Dropping it costs whole multiples of ln 2.
        """
        tot = 0.0
        for a, c in enumerate(self.A):
            acc = self.logf[a].copy()
            for i, u in enumerate(c):
                sh = [1] * len(c)
                sh[i] = 2
                acc = acc + self.m_va[(u, a)].reshape(sh)
            tot += float(logsumexp(acc))
        for v in self.nodes:
            acc = np.zeros(2)
            for a, c in enumerate(self.A):
                if v in c:
                    acc = acc + self.m_av[(a, v)]
            tot += float(logsumexp(acc))
        for (a, v), m in self.m_av.items():
            tot -= float(logsumexp(m + self.m_va[(v, a)]))
        return tot


def solve(G, n, bJ):
    """Best fixed point over the damping ladder, plus the exact answer."""
    cx = [sorted(c) for c in nx.find_cliques(G)]
    edges = sorted({tuple(sorted(e)) for e in G.edges()})
    exact = exact_log_Z(ising_factors(edges, bJ), range(n))
    # a bond inside two complexes is summed by both: the double count
    cov = {}
    for c in cx:
        for e in combinations(sorted(c), 2):
            cov[e] = cov.get(e, 0) + 1
    best = None
    for d in DAMPING:
        bp = ChygraphBP(cx, bJ, damping=d, edges=edges).run()
        if best is None or bp.residual < best.residual:
            best = bp
        if bp.residual < TOL:
            break
    return {
        'n': n, 'beta_J': bJ, 'n_cliques': len(cx),
        'doubled_bonds': sum(1 for k in cov.values() if k > 1),
        'n_bonds': len(edges),
        # The two candidate measures of overlap, so Ch. 14 can report the
        # correlation against both from one file: `doubled_bonds` counts the
        # defect itself, a bond summed by two complexes, and `shared_2plus` is
        # the coarser proxy, the fraction of INTERSECTING complex pairs that
        # meet in two or more atoms -- the same quantity, and the same
        # definition, that gbp_cliques/gbp_karrer/gbp_real record.
        'shared_2plus': overlap_profile(cx)['shared_2plus'],
        'exact': exact, 'cavity': best.log_Z() - exact,
        'residual': best.residual,
    }


def validate():
    """The recursion must be exact where the chygraph is treelike."""
    print('validation, on chygraphs the assumption holds for:')
    cases = [('chain of 6', nx.path_graph(6)),
             ('star, 7 leaves', nx.star_graph(7)),
             ('two triangles at a node', nx.Graph(
                 [(0, 1), (1, 2), (2, 0), (2, 3), (3, 4), (4, 2)])),
             ('tree of triangles', nx.Graph(
                 [(0, 1), (1, 2), (2, 0), (2, 3), (3, 4), (4, 2),
                  (0, 5), (5, 6), (6, 0)]))]
    ok = True
    for name, G in cases:
        for bJ in (0.3, 0.8, 1.5):
            r = solve(G, G.number_of_nodes(), bJ)
            good = abs(r['cavity']) < 1e-9 and r['doubled_bonds'] == 0
            ok &= good
            print(f"  {name:<24} bJ={bJ}  doubled={r['doubled_bonds']}  "
                  f"error={r['cavity']:+.2e}  {'ok' if good else 'FAIL'}")
    assert ok, 'the recursion is not exact on a treelike chygraph'
    print('  exact wherever complexes meet in at most one node, as required\n')


def main():
    validate()
    from gbp_cliques import SIZES, instance
    from gbp_real import load
    from merge import karrer_graph

    rows = []
    kbar = {n: k for n, k in SIZES}
    seen = set()
    for r in json.load(open(RESULTS / 'gbp_cliques.json')):
        if (r['n'], r['seed']) in seen:
            continue
        seen.add((r['n'], r['seed']))
        G = instance(r['n'], kbar[r['n']], r['seed'])
        for bJ in COUPLINGS:
            rows.append(dict(ensemble='hyperbolic', **solve(G, r['n'], bJ)))

    seen = set()
    for r in json.load(open(RESULTS / 'gbp_karrer.json'))['runs']:
        if (r['n'], r['seed']) in seen:
            continue
        seen.add((r['n'], r['seed']))
        G, _ = karrer_graph(r['n'], r['s_mean'], {3: r['triangles']}, r['seed'])
        G.remove_edges_from(nx.selfloop_edges(G))
        for bJ in COUPLINGS:
            rows.append(dict(ensemble='karrer', **solve(G, r['n'], bJ)))

    cache, seen = {}, set()
    for r in json.load(open(RESULTS / 'gbp_real.json')):
        if (r['network'], r['ego']) in seen:
            continue
        seen.add((r['network'], r['ego']))
        G = cache.setdefault(r['key'], load(r['key']))
        ego = sorted([r['ego']] + list(G[r['ego']]), key=str)
        H = nx.convert_node_labels_to_integers(G.subgraph(ego))
        for bJ in COUPLINGS:
            rows.append(dict(ensemble='real', network=r['network'],
                             **solve(H, H.number_of_nodes(), bJ)))

    json.dump(rows, OUT.open('w'), indent=1)
    summarise(rows)
    print('\nwrote', OUT)


def summarise(rows):
    for ens in ('hyperbolic', 'karrer', 'real'):
        sel = [r for r in rows if r['ensemble'] == ens]
        conv = [r for r in sel if r['residual'] < 1e-9]
        e = [abs(r['cavity']) for r in conv]
        db = [r['doubled_bonds'] / r['n_bonds'] for r in sel]
        print(f'\n{ens}: {len(sel)} runs, {len(conv)} reach a fixed point')
        if e:
            print(f'  |error in ln Z|: {min(e):.2e} to {max(e):.2e}, '
                  f'median {np.median(e):.2f}')
        print(f'  bonds lying in two or more complexes: '
              f'{100 * min(db):.0f} to {100 * max(db):.0f} per cent')


if __name__ == '__main__':
    main()
