"""Chapter 16's merged-family recursion at the stable fixed point.

``merge_lnz.py`` runs the Ch. 14 solver on the merged family from symmetric
messages, which keeps the paramagnetic fixed point whether or not it is
stable (Sec. 14.3).  On the 200 runs whose merged incidence structure is a
forest the fixed point is unique and exact, so nothing changes; on the 40
cyclic Karrer--Newman runs the reported error was the paramagnetic point's.
This reruns the 60 Karrer--Newman runs (same generator, same seeds) with the
stability of the paramagnetic point from the linearised recursion and the
polarised fixed point beside it, patches the stable-point error into
``results/merge_lnz.json`` (the original kept as
``merge_lnz_paramagnetic.json``), and sets beside it the unmerged recursion
at its stable fixed point, with the double count and with each bond in one
clique, from ``cavity_assigned.json``.

    python probe/merge_stable.py
"""

import json
import shutil
import sys
import types
from pathlib import Path

import networkx as nx
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'src'))
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'percolation' / 'src'))
from statmech.bethehessian import (linearised_operator, spectral_radius,  # noqa: E402
                                   trivial_blocks)
from statmech.gbp import exact_log_Z, ising_factors  # noqa: E402
from statmech.loopseries import BinaryFactorGraph  # noqa: E402
sys.modules.setdefault('hrg', types.SimpleNamespace(hrg_calibrated=None))
sys.path.append(str(Path(__file__).resolve().parents[2] / 'book' / 'figs'))
sys.path.append(str(Path(__file__).resolve().parent))
from merge import karrer_graph, merge_closure  # noqa: E402
from cavity_clique import ChygraphBP  # noqa: E402
from cavity_assigned import magnetisation, polarised  # noqa: E402
from merge_lnz import incidence_acyclic  # noqa: E402

RESULTS = Path(__file__).parent / 'results'
COUPLINGS = (0.3, 0.8)


def solve(G, bJ):
    n = G.number_of_nodes()
    merged, _ = merge_closure([frozenset(c) for c in nx.find_cliques(G)])
    cx = [sorted(c) for c in merged]
    edges = sorted({tuple(sorted(e)) for e in G.edges()})
    exact = exact_log_Z(ising_factors(edges, bJ), range(n))
    para = ChygraphBP(cx, bJ, damping=0.0, edges=edges).run()
    fg = BinaryFactorGraph.promoted(cx, edges, bJ, 'all')
    rho = spectral_radius(linearised_operator(fg, trivial_blocks(fg)))
    pol = polarised(cx, edges, bJ, 'all')
    out = dict(gbp_para=para.log_Z() - exact, gbp_rho=rho,
               gbp_pol=pol.log_Z() - exact, gbp_pol_m=magnetisation(pol),
               gbp_pol_residual=pol.residual, acyclic=bool(incidence_acyclic(cx)))
    out['gbp'] = out['gbp_para'] if rho < 1 else out['gbp_pol']
    out['gbp_stable_para'] = bool(rho < 1)
    return out


def main():
    src = RESULTS / 'merge_lnz.json'
    backup = RESULTS / 'merge_lnz_paramagnetic.json'
    if not backup.exists():
        shutil.copy(src, backup)
    rows = json.load(open(backup))
    assigned = {(r['n'], r['seed'], r['beta_J']): r
                for r in json.load(open(RESULTS / 'cavity_assigned.json'))
                if r['ensemble'] == 'karrer'}
    runs = json.load(open(RESULTS / 'gbp_karrer.json'))['runs']
    spec = {(r['n'], r['seed']): r for r in runs}
    k = 0
    for row in rows:
        if row['ensemble'] != 'karrer':
            continue
        # merge_lnz wrote the karrer rows in the order of gbp_karrer.json,
        # two couplings per instance
        pass
    # rebuild the karrer rows from the generator, in the same order
    out = [r for r in rows if r['ensemble'] != 'karrer']
    seen = []
    for r in runs:
        if (r['n'], r['seed']) in seen:
            continue
        seen.append((r['n'], r['seed']))
    old = [r for r in rows if r['ensemble'] == 'karrer']
    assert len(old) == 2 * len(seen)
    for i, (n, seed) in enumerate(seen):
        r = spec[(n, seed)]
        G, _ = karrer_graph(n, r['s_mean'], {3: r['triangles']}, seed)
        G.remove_edges_from(nx.selfloop_edges(G))
        for j, bJ in enumerate(COUPLINGS):
            row = dict(old[2 * i + j])
            assert row['n'] == n and row['beta_J'] == bJ
            new = solve(G, bJ)
            assert new['acyclic'] == row['acyclic']
            assert abs(new['gbp_para'] - row['gbp']) < 1e-8, (new['gbp_para'], row['gbp'])
            row.update(new)
            row['seed'] = seed
            a = assigned[(n, seed, bJ)]
            assert abs(a['exact'] - row['exact']) < 1e-8
            row['unmerged'] = a['cavity']
            row['unmerged_once'] = a['cavity_once']
            out.append(row)
            k += 1
    json.dump(out, src.open('w'), indent=1)
    print(f'rewrote {k} Karrer--Newman rows of {src.name}; original in {backup.name}')
    summarise([r for r in out if r['ensemble'] == 'karrer'])


def summarise(sel):
    cy = [r for r in sel if not r['acyclic']]
    ac = [r for r in sel if r['acyclic']]
    print(f'{len(sel)} runs, {len(ac)} acyclic (all exact: '
          f'{all(abs(r["gbp"]) < 1e-9 for r in ac)}), {len(cy)} cyclic')
    print(f'  paramagnetic point unstable on {sum(not r["gbp_stable_para"] for r in sel)}/{len(sel)} '
          f'({sum(not r["gbp_stable_para"] for r in cy)} of the cyclic)')
    for lab, key in (('merged, paramagnetic', 'gbp_para'), ('merged, stable', 'gbp'),
                     ('unmerged double count, stable', 'unmerged'),
                     ('unmerged assigned, stable', 'unmerged_once')):
        e = np.array([abs(r[key]) for r in cy])
        print(f'  cyclic 40, {lab:<32} {e.min():.1e} to {e.max():.2f}, median {e.median() if hasattr(e, "median") else np.median(e):.3f}')
    for bJ in COUPLINGS:
        s = [r for r in sel if r['beta_J'] == bJ]
        for lab, key in (('merged', 'gbp'), ('unmerged dc', 'unmerged'), ('assigned', 'unmerged_once')):
            e = np.array([abs(r[key]) for r in s])
            print(f'  bJ={bJ} all 30, {lab:<12} median {np.median(e):.3f} max {e.max():.2f}')
    m = np.array([abs(r['gbp']) for r in sel]); u = np.array([abs(r['unmerged']) for r in sel])
    o = np.array([abs(r['unmerged_once']) for r in sel])
    print(f'  all 60: medians merged {np.median(m):.3f}, unmerged {np.median(u):.3f}, assigned {np.median(o):.3f}; '
          f'merged smaller than assigned on {int(np.sum(m < o))}/60')


if __name__ == '__main__':
    main()
