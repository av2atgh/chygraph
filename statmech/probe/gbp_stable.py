"""Chapter 15's GBP fixed points: which one, and is it stable?

The three GBP probes start from uniform messages, and a symmetric iteration
keeps the symmetric fixed point (Sec. 14.3): the magnetisations they record
are zero on most runs.  On a junction tree that point is the unique fixed
point and the exact answer.  Off it, this file asks whether the symmetric
fixed point the probes reached is stable -- the spectral radius of the
undamped sweep map, by finite differences at the fixed point -- and runs
the same iteration from a polarised start to get the other fixed point.  The
fixed point to report is the symmetric one while stable and the polarised
one otherwise, as for the chygraph recursion in ``cavity_assigned.py``,
whose stable-point errors are set beside it for the same instances.
Karrer--Newman and real instances are the probes' own (same caches, same
seeds); hyperbolic ones are the fresh draws of ``cavity_assigned.py``.
Results in ``results/gbp_stable.json``; ``book/figs/overlap.py`` pairs them.

    python probe/gbp_stable.py
"""

import json
import sys
import time
import types
from pathlib import Path

import networkx as nx
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'src'))
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'percolation' / 'src'))
from statmech.gbp import GBP, exact_log_Z, ising_factors, _normalise  # noqa: E402
from statmech.region import RegionGraph  # noqa: E402
sys.modules.setdefault('hrg', types.SimpleNamespace(hrg_calibrated=None))
sys.path.append(str(Path(__file__).resolve().parents[2] / 'book' / 'figs'))
sys.path.append(str(Path(__file__).resolve().parent))
from cavity_assigned import COUPLINGS, SEEDS, SIZES  # noqa: E402
from loopseries import hyperbolic_graph  # noqa: E402

RESULTS = Path(__file__).parent / 'results'
OUT = RESULTS / 'gbp_stable.json'
DAMPING = (0.5, 0.9, 0.97, 0.995, 0.999)
SWEEPS = 8000
CONVERGED = 1e-9
CONSISTENT = 1e-6
DMAX = 3000        # message coordinates above which the Jacobian is not formed


def ladder(rg, f, init=None, quick=False):
    best = None
    for d in (DAMPING[:3] if quick else DAMPING):
        g = GBP(rg, f, damping=d)
        if init is not None:
            init(g)
        g.run(3000 if quick else SWEEPS)
        if best is None or g.residual < best[0].residual:
            best = (g, d)
        if g.residual < CONVERGED:
            break
    return best


def polarised_init(h=0.5):
    def init(g):
        for (p, r) in g.edges:
            vs = g.vars[r]
            acc = np.zeros((2,) * len(vs))
            for i in range(len(vs)):
                sh = [1] * len(vs)
                sh[i] = 2
                acc = acc + (h * np.array([1.0, -1.0])).reshape(sh)
            g.m[(p, r)] = _normalise(acc)
    return init


def sweep_radius(g, eps=1e-6):
    """Largest real part among the eigenvalues of the undamped sweep map,
    linearised at the current messages by finite differences."""
    keys = list(g.edges)
    shapes = [g.m[k].shape for k in keys]
    sizes = [int(np.prod(s)) for s in shapes]
    off = np.concatenate([[0], np.cumsum(sizes)[:-1]]).astype(int)
    base = np.concatenate([g.m[k].ravel() for k in keys])

    def load(x):
        for k, o, n, s in zip(keys, off, sizes, shapes):
            g.m[k] = x[o:o + n].reshape(s)

    def sweep_from(x):
        load(x)
        g.sweep()
        return np.concatenate([g.m[k].ravel() for k in keys])

    lam = g.damping
    g.damping = 0.0
    f0 = sweep_from(base)
    D = len(base)
    J = np.empty((D, D))
    for k in range(D):
        x = base.copy()
        x[k] += eps
        J[:, k] = (sweep_from(x) - f0) / eps
    ev = np.linalg.eigvals(J)
    gauge = np.abs(ev - 1.0) < 1e-4
    rho = float(ev.real[~gauge].max()) if (~gauge).any() else 0.0
    # the growth rate of successive differences from a perturbed start,
    # which the gauge directions (eigenvalue one) do not feed
    rng = np.random.default_rng(0)
    x = base + 1e-7 * rng.normal(size=D)
    d = []
    for _ in range(60):
        y = sweep_from(x)
        d.append(float(np.abs(y - x).max()))
        x = y
        if d[-1] < 1e-15 or d[-1] > 1e-1:
            break
    d = np.array(d)
    ok = d > 1e-15
    if ok.sum() >= 10:
        t = np.arange(len(d))[ok][-min(30, ok.sum()):]
        slope = np.polyfit(t, np.log(d[ok][-len(t):]), 1)[0]
        growth = float(np.exp(slope))
    else:
        growth = 0.0
    load(base)
    g.damping = lam
    return rho, growth, int(gauge.sum()), float(np.abs(f0 - base).max())


def arnoldi_radius(g, k=8, eps=1e-6, tol=1e-8):
    """Largest real part among the eigenvalues of the undamped sweep map at
    the current messages, by Arnoldi on finite-difference Jacobian-vector
    products (one sweep each), the gauge eigenvalues at one excluded."""
    import scipy.sparse.linalg as spla
    keys = list(g.edges)
    shapes = [g.m[k].shape for k in keys]
    sizes = [int(np.prod(sh)) for sh in shapes]
    off = np.concatenate([[0], np.cumsum(sizes)[:-1]]).astype(int)
    base = np.concatenate([g.m[k].ravel() for k in keys])

    def load(x):
        for kk, o, n, sh in zip(keys, off, sizes, shapes):
            g.m[kk] = x[o:o + n].reshape(sh)

    def sweep_from(x):
        load(x)
        g.sweep()
        return np.concatenate([g.m[kk].ravel() for kk in keys])

    lam = g.damping
    g.damping = 0.0
    f0 = sweep_from(base)
    D = len(base)

    def matvec(v):
        v = np.asarray(v).ravel()
        nv = np.linalg.norm(v)
        if nv == 0:
            return np.zeros(D)
        return (sweep_from(base + eps * v / nv) - f0) * (nv / eps)

    A = spla.LinearOperator((D, D), matvec=matvec, dtype=float)
    try:
        ev = spla.eigs(A, k=min(k, D - 2), which='LM', tol=tol, maxiter=2000,
                       return_eigenvectors=False)
    except spla.ArpackNoConvergence as e:
        ev = e.eigenvalues
    load(base)
    g.damping = lam
    gauge = np.abs(ev - 1.0) < 1e-3
    rest = ev[~gauge]
    return (float(rest.real.max()) if len(rest) else 0.0), int(gauge.sum()), len(ev)


def belief_growth(g, sweeps=120, eps=1e-6):
    """Stability of the current fixed point under the map at the damping it
    was reached with: perturb the messages, iterate, and fit the growth of
    the deviation of the *beliefs*, which the message gauge leaves alone.
    Returns the growth factor per sweep (below one: stable)."""
    keys = list(g.edges)
    base = {k: g.m[k].copy() for k in keys}
    b0 = np.concatenate([g.belief(r).ravel() for r in g.regions])
    rng = np.random.default_rng(0)
    for k in keys:
        g.m[k] = _normalise(base[k] + eps * rng.normal(size=base[k].shape))
    d = []
    for _ in range(sweeps):
        g.sweep()
        b = np.concatenate([g.belief(r).ravel() for r in g.regions])
        d.append(float(np.abs(b - b0).max()))
        if d[-1] > 1e-2 or d[-1] < 1e-15:
            break
    for k in keys:
        g.m[k] = base[k]
    d = np.array(d)
    if len(d) < 20:
        return 2.0 if d[-1] > 1e-2 else 0.0
    t = np.arange(len(d))[len(d) // 2:]
    slope = np.polyfit(t, np.log(np.maximum(d[len(d) // 2:], 1e-300)), 1)[0]
    return float(np.exp(slope))


def study(G, bJ, cached=None):
    n = G.number_of_nodes()
    cx = [sorted(c) for c in nx.find_cliques(G)]
    rg = RegionGraph(cx, max_rounds=6)
    edges = sorted({tuple(sorted(e)) for e in G.edges()})
    f = ising_factors(edges, bJ)
    exact = exact_log_Z(f, range(n))
    out = dict(n=n, beta_J=bJ, chordal=bool(nx.is_chordal(G)), exact=exact,
               n_regions=len(rg.regions))
    if not rg.counting_is_valid():
        out['valid'] = False
        return out
    out['valid'] = True
    # the symmetric fixed point: reproduce the probe's, or skip the rerun where
    # the probe already found no fixed point
    if cached is not None and cached['residual'] >= CONVERGED:
        out.update(sym_settled=False, gbp_sym=cached['gbp'],
                   consistency_sym=cached['consistency'])
    else:
        g, d = ladder(rg, f)
        cons = g.consistency() if np.isfinite(g.residual) else np.inf
        out.update(gbp_sym=g.log_Z() - exact, consistency_sym=cons,
                   sym_settled=bool(g.residual < CONVERGED and cons < CONSISTENT),
                   damping_sym=d)
        if cached is not None and cached['residual'] < CONVERGED:
            out['reproduced'] = abs(out['gbp_sym'] - cached['gbp'])
        out['dim'] = sum(int(np.prod(g.m[k].shape)) for k in g.edges)
        if out['sym_settled']:
            out['m_sym'] = float(np.mean(np.abs(list(g.magnetisation().values()))))
            if out['dim'] <= DMAX:
                rho, growth, ngauge, drift = sweep_radius(g)
                out.update(rho=rho, n_gauge=ngauge, sweep_drift=drift)
            else:
                out['rho'] = None
    quick = cached is not None and cached['residual'] >= CONVERGED
    g, d = ladder(rg, f, polarised_init(), quick=quick)
    cons = g.consistency() if np.isfinite(g.residual) else np.inf
    out.update(gbp_pol=g.log_Z() - exact, consistency_pol=cons,
               pol_settled=bool(g.residual < CONVERGED and cons < CONSISTENT),
               m_pol=(float(np.mean(np.abs(list(g.magnetisation().values()))))
                      if np.isfinite(g.residual) else None), damping_pol=d)
    if out['pol_settled'] and sum(int(np.prod(g.m[k].shape)) for k in g.edges) <= DMAX:
        rho, growth, ngauge, drift = sweep_radius(g)
        out['rho_pol'] = rho
    elif out['pol_settled']:
        out['rho_pol'] = None
    out['pol_returns'] = bool(out['pol_settled'] and out['m_pol'] is not None
                              and out['m_pol'] < 1e-6)
    # the fixed point to report: the symmetric one while stable (by the
    # Jacobian, or, where it is not formed, by a polarised start returning
    # to it), the polarised one otherwise
    sym_stable = out['sym_settled'] and (out['rho'] < 1 if out.get('rho') is not None
                                         else out['pol_returns'])
    if sym_stable:
        out['gbp'], out['point'] = out['gbp_sym'], 'symmetric'
    elif out['pol_settled'] and out.get('rho_pol') is None:
        out['gbp'], out['point'] = out['gbp_pol'], 'polarised, unassessed'
    elif out['pol_settled'] and out['rho_pol'] < 1:
        out['gbp'], out['point'] = out['gbp_pol'], 'polarised'
    elif out['sym_settled']:
        out['gbp'], out['point'] = out['gbp_sym'], 'symmetric, unstable'
    elif out['pol_settled']:
        out['gbp'], out['point'] = out['gbp_pol'], 'polarised, unstable'
    else:
        out['gbp'], out['point'] = None, 'none'
    return out


def main():
    from gbp_real import load
    from merge import karrer_graph
    assigned = json.load(open(RESULTS / 'cavity_assigned.json'))
    rows, t0 = [], time.time()
    hyp = [r for r in assigned if r['ensemble'] == 'hyperbolic']
    R, k = {}, 0
    for n, kbar in SIZES:
        for seed in SEEDS:
            G, R[n] = hyperbolic_graph(n, kbar, 2.5, np.random.default_rng(seed), R.get(n))
            G.remove_nodes_from([v for v in list(G) if G.degree(v) == 0])
            G = nx.convert_node_labels_to_integers(G)
            for bJ in COUPLINGS:
                a = hyp[k]
                k += 1
                assert a['n'] == G.number_of_nodes() and a['beta_J'] == bJ and a['seed'] == seed
                row = dict(ensemble='hyperbolic', seed=seed, **study(G, bJ))
                row.update(bp=a['cavity'], bp_once=a['cavity_once'],
                           doubled=a['doubled_bonds'] / a['n_bonds'])
                rows.append(row)
        print(f'hyperbolic n = {n} done, {time.time() - t0:.0f} s', flush=True)
        json.dump(rows, OUT.open('w'), indent=1)
    kar = {(r['n'], r['seed'], r['beta_J']): r for r in assigned if r['ensemble'] == 'karrer'}
    cached = {(r['n'], r['seed'], r['beta_J']): r
              for r in json.load(open(RESULTS / 'gbp_karrer.json'))['runs']}
    seen = set()
    for r in json.load(open(RESULTS / 'gbp_karrer.json'))['runs']:
        if (r['n'], r['seed']) in seen:
            continue
        seen.add((r['n'], r['seed']))
        G, _ = karrer_graph(r['n'], r['s_mean'], {3: r['triangles']}, r['seed'])
        G.remove_edges_from(nx.selfloop_edges(G))
        for bJ in COUPLINGS:
            key = (r['n'], r['seed'], bJ)
            row = dict(ensemble='karrer', seed=r['seed'], **study(G, bJ, cached[key]))
            a = kar[key]
            row.update(bp=a['cavity'], bp_once=a['cavity_once'],
                       doubled=a['doubled_bonds'] / a['n_bonds'])
            rows.append(row)
    print(f'karrer done, {time.time() - t0:.0f} s', flush=True)
    json.dump(rows, OUT.open('w'), indent=1)
    real = {(r['network'], r['ego'], r['beta_J']): r for r in assigned if r['ensemble'] == 'real'}
    cached = {(r['network'], str(r['ego']), r['beta_J']): r
              for r in json.load(open(RESULTS / 'gbp_real.json'))}
    cache, seen = {}, set()
    for r in json.load(open(RESULTS / 'gbp_real.json')):
        if (r['network'], r['ego']) in seen:
            continue
        seen.add((r['network'], r['ego']))
        G = cache.setdefault(r['key'], load(r['key']))
        ego = sorted([r['ego']] + list(G[r['ego']]), key=str)
        H = nx.convert_node_labels_to_integers(G.subgraph(ego))
        for bJ in COUPLINGS:
            key = (r['network'], str(r['ego']), bJ)
            row = dict(ensemble='real', network=r['network'], ego=str(r['ego']),
                       **study(H, bJ, cached[key]))
            a = real[key]
            row.update(bp=a['cavity'], bp_once=a['cavity_once'],
                       doubled=a['doubled_bonds'] / a['n_bonds'])
            rows.append(row)
        if len(rows) % 20 == 0:
            json.dump(rows, OUT.open('w'), indent=1)
            print(f'  {len(rows)} rows, {time.time() - t0:.0f} s', flush=True)
    print(f'real done, {time.time() - t0:.0f} s', flush=True)
    json.dump(rows, OUT.open('w'), indent=1)
    summarise(rows)
    print('wrote', OUT)


def summarise(rows):
    rep = [r['reproduced'] for r in rows if 'reproduced' in r]
    print(f'{len(rows)} runs; cached settled fixed points reproduced to {max(rep):.1e} on {len(rep)}')
    for ens in ('hyperbolic', 'karrer', 'real'):
        sel = [r for r in rows if r['ensemble'] == ens and r['valid']]
        nc = [r for r in sel if not r['chordal']]
        ch = [r for r in sel if r['chordal']]
        print(f'{ens}: {len(sel)} valid, {len(ch)} chordal (all symmetric settled: '
              f'{all(r["sym_settled"] for r in ch)}, stable where formed: '
              f'{all(r["rho"] < 1 for r in ch if r.get("rho") is not None)}, '
              f'Jacobian formed on {sum(r.get("rho") is not None for r in ch)}), {len(nc)} non-chordal:')
        s = [r for r in nc if r['sym_settled']]
        print(f'   symmetric point settled on {len(s)}, Jacobian formed on {sum(r["rho"] is not None for r in s)}, '
              f'unstable on {sum(r["rho"] is not None and r["rho"] >= 1 for r in s)} '
              f'(polarised start returns to it on {sum(r["pol_returns"] for r in s)}); '
              f'polarised settled on {sum(r["pol_settled"] for r in nc)}, of which stable '
              f'{sum(r["pol_settled"] and r.get("rho_pol") is not None and r["rho_pol"] < 1 for r in nc)}')
        pts = {}
        for r in nc:
            pts[r['point']] = pts.get(r['point'], 0) + 1
        print(f'   reported point: {pts}')
        have = [r for r in nc if r['gbp'] is not None]
        g = np.array([abs(r['gbp']) for r in have]); b = np.array([abs(r['bp']) for r in have])
        o = np.array([abs(r['bp_once']) for r in have])
        gs = np.array([abs(r['gbp_sym']) for r in nc if r['sym_settled']])
        if len(have):
            print(f'   non-chordal with a fixed point ({len(have)}): GBP median {np.median(g):.3f} max {g.max():.2f}; '
                  f'symmetric-point median {np.median(gs):.3f}; BP stable median {np.median(b):.2f}; '
                  f'assigned median {np.median(o):.3f}; GBP worse than BP on {int(np.sum(g >= b))}, '
                  f'worse than assigned on {int(np.sum(g >= o))}')


if __name__ == '__main__':
    main()
