"""Which clique should a shared bond go to?  Rules against an oracle.

Sec. 14.3 gives each bond lying in several cliques to the largest of them,
ties to the first, and Sec. 24.3 reports the spread over random choices.
This asks whether a better rule exists that uses nothing but the structure.
For every run of ``cavity_assigned.py`` (same instances) it evaluates, on
the promoted factor graph with the bonds assigned by the rule, the error of
the recursion at its stable fixed point, the protocol of Sec. 14.3 (the
symmetric point while the linearised recursion finds it stable, else the
one reached from a polarised start):

    largest    the largest clique, ties to the first (the book's rule)
    smallest   the smallest clique, ties to the first
    strongest  the clique in which the pair is most correlated at zero
               field through the clique's other bonds
    weakest    the least correlated
    balance    greedily, the clique holding the smallest share of its
               own bonds so far
    random     ten uniform draws: the median and the best
    oracle     local search from `largest`, flipping one bond's owner at a
               time while the error falls, at most forty evaluations --
               needs the exact answer, so it is a floor and not a rule

Results in ``results/assignment_rules.json``.

    python probe/assignment_rules.py [--subset]
"""

import json
import sys
import time
import types
from itertools import combinations
from pathlib import Path

import networkx as nx
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'src'))
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'percolation' / 'src'))
from statmech.gbp import exact_log_Z, ising_factors  # noqa: E402
from statmech.loopseries import BinaryFactorGraph  # noqa: E402
sys.modules.setdefault('hrg', types.SimpleNamespace(hrg_calibrated=None))
sys.path.append(str(Path(__file__).resolve().parents[2] / 'book' / 'figs'))
sys.path.append(str(Path(__file__).resolve().parent))
from cavity_assigned import COUPLINGS, SEEDS, SIZES  # noqa: E402
from loopseries import hyperbolic_graph  # noqa: E402

RESULTS = Path(__file__).parent / 'results'
OUT = RESULTS / 'assignment_rules.json'
RANDOM = 10


def shared_bonds(cx, edges):
    holders = {}
    for a, c in enumerate(cx):
        for e in combinations(c, 2):
            if e in edges:
                holders.setdefault(e, []).append(a)
    return {e: hs for e, hs in holders.items() if len(hs) > 1}


def indirect_correlation(c, bond, edges, bJ):
    """<s_u s_v> at zero field in clique `c` with all its bonds but `bond`."""
    n = len(c)
    pos = {v: i for i, v in enumerate(c)}
    s = np.array([[1.0 if (i >> (n - 1 - b)) & 1 == 0 else -1.0
                   for b in range(n)] for i in range(2 ** n)])
    e = np.zeros(2 ** n)
    for p, q in combinations(c, 2):
        if (p, q) in edges and (p, q) != bond:
            e += s[:, pos[p]] * s[:, pos[q]]
    w = np.exp(bJ * e - (bJ * e).max())
    u, v = bond
    return float((w * s[:, pos[u]] * s[:, pos[v]]).sum() / w.sum())


def rule_owner(rule, cx, edges, bJ, shared, rng=None):
    owner = {}
    if rule in ('largest', 'smallest'):
        sign = -1 if rule == 'largest' else 1
        for e, hs in shared.items():
            owner[e] = min(hs, key=lambda a: (sign * len(cx[a]), a))
    elif rule in ('strongest', 'weakest'):
        for e, hs in shared.items():
            corr = {a: indirect_correlation(cx[a], e, edges, bJ) for a in hs}
            owner[e] = (max if rule == 'strongest' else min)(hs, key=lambda a: (corr[a], -a))
    elif rule == 'balance':
        held = {a: 0 for a in range(len(cx))}
        total = {a: sum(1 for x in combinations(c, 2) if x in edges) for a, c in enumerate(cx)}
        for e in sorted(shared):
            hs = shared[e]
            a = min(hs, key=lambda a: (held[a] / max(total[a], 1), a))
            owner[e] = a
            held[a] += 1
        # bonds in one clique count towards its share too
        return owner
    elif rule == 'random':
        for e, hs in shared.items():
            owner[e] = int(rng.choice(hs))
    else:
        raise ValueError(rule)
    return owner


def error(cx, edges, bJ, owner, exact):
    """The error at the stable fixed point, as in Sec. 14.3: the symmetric
    one while the linearised recursion says it is stable, otherwise the one
    reached from a polarised start.  A polarised point misses the ln 2 of
    the other sector, which is reported as part of the error."""
    from statmech.bethehessian import (linearised_operator, spectral_radius,
                                       trivial_blocks)
    from cavity_assigned import polarised
    from cavity_clique import ChygraphBP
    fg = BinaryFactorGraph.promoted(cx, edges, bJ, owner)
    rho = spectral_radius(linearised_operator(fg, trivial_blocks(fg)))
    cxl = [list(c) for c in cx]
    if rho < 1:
        bp = ChygraphBP(cxl, bJ, damping=0.0, edges=sorted(edges), assign=owner).run()
    else:
        bp = polarised(cxl, sorted(edges), bJ, owner)
    return bp.log_Z() - exact, bool(bp.residual < 1e-9)


def oracle(cx, edges, bJ, shared, exact, start, passes=6, max_evals=40):
    owner = dict(start)
    best, _ = error(cx, edges, bJ, owner, exact)
    evals = 1
    for _ in range(passes):
        improved = False
        for e, hs in shared.items():
            for a in hs:
                if a == owner[e] or evals >= max_evals:
                    continue
                trial = dict(owner)
                trial[e] = a
                err, _ = error(cx, edges, bJ, trial, exact)
                evals += 1
                if abs(err) < abs(best) - 1e-12:
                    best, owner, improved = err, trial, True
        if not improved:
            break
    return best, owner, evals


def study(G, bJ, rng):
    cx = [tuple(sorted(c)) for c in nx.find_cliques(G)]
    edges = {tuple(sorted(e)) for e in G.edges()}
    exact = exact_log_Z(ising_factors(sorted(edges), bJ), range(G.number_of_nodes()))
    shared = shared_bonds(cx, edges)
    out = dict(n=G.number_of_nodes(), beta_J=bJ, n_cliques=len(cx),
               shared=len(shared), exact=exact)
    owners = {}
    for rule in ('largest', 'smallest', 'strongest', 'weakest', 'balance'):
        owners[rule] = rule_owner(rule, cx, edges, bJ, shared)
        err, conv = error(cx, edges, bJ, owners[rule], exact)
        out[rule] = err
        out[rule + '_conv'] = conv
    rnd = []
    for _ in range(RANDOM):
        err, conv = error(cx, edges, bJ, rule_owner('random', cx, edges, bJ, shared, rng), exact)
        rnd.append(abs(err))
    out['random_median'] = float(np.median(rnd))
    out['random_best'] = float(min(rnd))
    best, owner, evals = oracle(cx, edges, bJ, shared, exact, owners['largest'])
    out['oracle'] = best
    out['oracle_evals'] = evals
    # how the oracle's assignment relates to the rules
    for rule in ('largest', 'smallest', 'strongest'):
        out['oracle_agrees_' + rule] = sum(owner[e] == owners[rule][e] for e in shared) / max(len(shared), 1)
    return out


def main(subset=False):
    from gbp_real import load
    from merge import karrer_graph
    rows, t0 = [], time.time()
    rng = np.random.default_rng(7)
    R = {}
    for n, kbar in SIZES:
        for seed in (SEEDS if not subset else list(SEEDS)[:3]):
            G, R[n] = hyperbolic_graph(n, kbar, 2.5, np.random.default_rng(seed), R.get(n))
            G.remove_nodes_from([v for v in list(G) if G.degree(v) == 0])
            G = nx.convert_node_labels_to_integers(G)
            for bJ in COUPLINGS:
                rows.append(dict(ensemble='hyperbolic', seed=seed, **study(G, bJ, rng)))
    print(f'hyperbolic done, {time.time() - t0:.0f} s', flush=True)
    seen = set()
    for r in json.load(open(RESULTS / 'gbp_karrer.json'))['runs']:
        if (r['n'], r['seed']) in seen or (subset and len(seen) >= 9):
            continue
        seen.add((r['n'], r['seed']))
        G, _ = karrer_graph(r['n'], r['s_mean'], {3: r['triangles']}, r['seed'])
        G.remove_edges_from(nx.selfloop_edges(G))
        for bJ in COUPLINGS:
            rows.append(dict(ensemble='karrer', seed=r['seed'], **study(G, bJ, rng)))
    print(f'karrer done, {time.time() - t0:.0f} s', flush=True)
    if not subset:
        json.dump(rows, OUT.open('w'), indent=1)
    cache, seen = {}, set()
    for r in json.load(open(RESULTS / 'gbp_real.json')):
        if (r['network'], r['ego']) in seen or (subset and len(seen) >= 18):
            continue
        seen.add((r['network'], r['ego']))
        G = cache.setdefault(r['key'], load(r['key']))
        ego = sorted([r['ego']] + list(G[r['ego']]), key=str)
        H = nx.convert_node_labels_to_integers(G.subgraph(ego))
        for bJ in COUPLINGS:
            rows.append(dict(ensemble='real', network=r['network'], ego=str(r['ego']),
                             **study(H, bJ, rng)))
        if not subset and len(rows) % 10 == 0:
            json.dump(rows, OUT.open('w'), indent=1)
            print(f'  {len(rows)} rows, {time.time() - t0:.0f} s', flush=True)
    print(f'real done, {time.time() - t0:.0f} s', flush=True)
    if not subset:
        json.dump(rows, OUT.open('w'), indent=1)
        print('wrote', OUT)
    summarise(rows)


def summarise(rows):
    keys = ('largest', 'smallest', 'strongest', 'weakest', 'balance', 'random_median', 'random_best', 'oracle')
    print(f'{len(rows)} runs; median |error| by rule, and how often each rule ties the oracle within 0.01:')
    for ens in ('hyperbolic', 'karrer', 'real', 'all'):
        for bJ in COUPLINGS + (None,):
            sel = [r for r in rows if (ens == 'all' or r['ensemble'] == ens) and (bJ is None or r['beta_J'] == bJ)]
            if not sel or (ens != 'all' and bJ is None):
                continue
            med = {k: np.median([abs(r[k]) for r in sel]) for k in keys}
            ties = {k: int(sum(abs(r[k]) < abs(r['oracle']) + 0.01 for r in sel)) for k in keys[:5]}
            print(f"  {ens:<10} bJ={bJ}: n={len(sel):3d}  " + '  '.join(f'{k} {med[k]:.3f}' for k in keys) + f'  | ties {ties}')
    ag = {k: np.mean([r['oracle_agrees_' + k] for r in rows if r['shared']]) for k in ('largest', 'smallest', 'strongest')}
    print('  oracle agrees with rule on fraction of shared bonds:', {k: round(v, 2) for k, v in ag.items()})
    print('  oracle evaluations per run, median', np.median([r['oracle_evals'] for r in rows]))


if __name__ == '__main__':
    main(subset='--subset' in sys.argv)
