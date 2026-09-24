"""Acyclicity percolation: where does the clique family stop being alpha-acyclic?

The chygraph cavity method is exact on an instance iff the family of maximal
cliques has a join tree, i.e. is alpha-acyclic (Beeri, Fagin, Maier, Yannakakis
1983), which for a graph is the same as chordality (Gavril 1974).  The
certificate is GYO reduction -- delete a vertex that lies in one clique only,
delete a clique contained in another, repeat -- which is leaf removal on the
incidence structure.  What GYO leaves is the *residue*: the part of the graph
that no join tree covers and that every tree-based calculation gets wrong.

This probe measures, on ensembles with extensive overlap (geometric ones),

  residue      fraction of vertices GYO does not remove
  giant        fraction of vertices in the largest connected residue component
  islands      mean size of a residue component (over vertices), i.e. the
               cardinality of the meta-complex one would have to sum exactly
  merged       fraction of vertices in the largest meta-complex of Ch. 16's
               closure (cliques sharing >= 2 atoms merged, to convergence)

for random geometric graphs on the torus in d = 1, 2, 3 and for hyperbolic
random graphs at several tau, against mean degree.  d = 1 is an interval graph
and chordal by construction, so the residue is zero identically: the anchor.

    python probe/acyclicity.py [rgg|hrg|all]

Writes probe/results/acyclicity_<ensemble>.csv.
"""

import csv
import sys
import time
from collections import defaultdict
from pathlib import Path

import networkx as nx
import numpy as np
from scipy.special import gamma as Gamma

OUT = Path(__file__).parent / 'results'


# ------------------------------------------------------------------ graphs

def rgg_torus(n, d, kbar, rng, wrap=None):
    """Random geometric graph on the unit torus [0,1)^d with mean degree kbar.

    In d = 1 the default is the *line* (no wrap): the ring closes one long
    chordless cycle through every clique, which GYO keeps although the cavity
    method does not mind it.  Interval graphs are chordal; circular-arc graphs
    are not.
    """
    if wrap is None:
        wrap = d > 1
    vol = np.pi ** (d / 2) / Gamma(d / 2 + 1)
    r = (kbar / (n * vol)) ** (1.0 / d)
    x = rng.random((n, d))
    # cell list
    m = max(1, int(1.0 / r))
    cell = (x * m).astype(int) % m
    key = np.ravel_multi_index(cell.T, (m,) * d)
    order = np.argsort(key)
    key_s = key[order]
    starts = np.searchsorted(key_s, np.arange(m ** d))
    ends = np.searchsorted(key_s, np.arange(m ** d), side='right')
    offsets = np.array(np.meshgrid(*([[-1, 0, 1]] * d), indexing='ij')).reshape(d, -1).T
    edges = []
    cells = np.array(np.unravel_index(np.arange(m ** d), (m,) * d)).T
    for c in range(m ** d):
        i = order[starts[c]:ends[c]]
        if len(i) == 0:
            continue
        for off in offsets:
            nb = cells[c] + off
            if not wrap and (np.any(nb < 0) or np.any(nb >= m)):
                continue
            nb = tuple(nb % m)
            c2 = np.ravel_multi_index(nb, (m,) * d)
            if c2 < c:
                continue
            j = order[starts[c2]:ends[c2]]
            if len(j) == 0:
                continue
            dx = x[i][:, None, :] - x[j][None, :, :]
            if wrap:
                dx = dx - np.round(dx)
            dist2 = (dx ** 2).sum(-1)
            a, b = np.nonzero(dist2 < r * r)
            if c2 == c:
                keep = a < b
                a, b = a[keep], b[keep]
            edges.append(np.column_stack([i[a], j[b]]))
    E = np.vstack(edges) if edges else np.zeros((0, 2), int)
    G = nx.Graph()
    G.add_nodes_from(range(n))
    G.add_edges_from(map(tuple, E))
    return G


def hrg(n, tau, kbar, rng, nbands=24):
    """Hyperbolic random graph (Krioukov et al. 2010), curvature 1, alpha=(tau-1)/2.

    R is tuned by bisection on the mean degree measured on a subsample; the
    graph is then built band by band in radius, each band searched only over
    the angular window its smallest radius allows.
    """
    alpha = (tau - 1) / 2
    theta = rng.random(n) * 2 * np.pi
    u = rng.random(n)

    def radii(R):
        return np.arccosh(1 + u * (np.cosh(alpha * R) - 1)) / alpha

    def mean_degree(R, m=3000):
        r = radii(R)
        idx = rng.choice(n, size=min(m, n), replace=False)
        deg = 0.0
        for i in idx:
            dth = np.pi - np.abs(np.pi - np.abs(theta - theta[i]))
            ch = np.cosh(r[i]) * np.cosh(r) - np.sinh(r[i]) * np.sinh(r) * np.cos(dth)
            deg += (ch < np.cosh(R)).sum() - 1
        return deg / len(idx)

    lo, hi = 2.0, 4 * np.log(n)
    for _ in range(25):
        mid = 0.5 * (lo + hi)
        if mean_degree(mid) > kbar:
            lo = mid
        else:
            hi = mid
    R = 0.5 * (lo + hi)
    r = radii(R)
    chR = np.cosh(R)
    # bands in radius; within a band, points sorted by angle
    edges_r = np.quantile(r, np.linspace(0, 1, nbands + 1))
    band = np.clip(np.searchsorted(edges_r, r, side='right') - 1, 0, nbands - 1)
    members = [np.nonzero(band == b)[0] for b in range(nbands)]
    members = [m[np.argsort(theta[m])] for m in members]
    thetas = [theta[m] for m in members]
    rmins = [r[m].min() if len(m) else np.inf for m in members]
    src_l, dst_l = [], []
    for i in range(n):
        ri, ti = r[i], theta[i]
        for b in range(nbands):
            m = members[b]
            if len(m) == 0:
                continue
            cosw = (np.cosh(ri) * np.cosh(rmins[b]) - chR) / (np.sinh(ri) * np.sinh(rmins[b]))
            if cosw >= 1:
                continue
            w = np.pi if cosw <= -1 else np.arccos(cosw)
            th = thetas[b]
            if w >= np.pi:
                js = np.arange(len(m))
            else:
                lo_t, hi_t = ti - w, ti + w
                a, c = np.searchsorted(th, lo_t), np.searchsorted(th, hi_t, side='right')
                js = np.arange(a, c)
                if lo_t < 0:
                    js = np.concatenate([js, np.arange(np.searchsorted(th, lo_t + 2 * np.pi), len(m))])
                if hi_t > 2 * np.pi:
                    js = np.concatenate([js, np.arange(0, np.searchsorted(th, hi_t - 2 * np.pi, side='right'))])
            j = m[js]
            j = j[j > i]
            if len(j) == 0:
                continue
            dth = np.pi - np.abs(np.pi - np.abs(theta[j] - ti))
            ch = np.cosh(ri) * np.cosh(r[j]) - np.sinh(ri) * np.sinh(r[j]) * np.cos(dth)
            hit = j[ch < chR]
            src_l.append(np.full(len(hit), i))
            dst_l.append(hit)
    G = nx.Graph()
    G.add_nodes_from(range(n))
    if src_l:
        G.add_edges_from(zip(np.concatenate(src_l).tolist(), np.concatenate(dst_l).tolist()))
    return G


# ------------------------------------------------------------------ GYO

def gyo_residue(cliques, n):
    """GYO reduction of the hypergraph of maximal cliques.

    Returns the residue hypergraph (list of frozensets) after removing, to a
    fixed point, vertices in exactly one hyperedge and hyperedges contained in
    another.  Empty iff the family is alpha-acyclic.
    """
    H = {i: set(c) for i, c in enumerate(cliques)}
    inc = defaultdict(set)
    for i, c in H.items():
        for v in c:
            inc[v].add(i)
    changed = True
    while changed:
        changed = False
        # ear vertices
        for v in [v for v, s in inc.items() if len(s) == 1]:
            (i,) = inc[v]
            H[i].discard(v)
            del inc[v]
            changed = True
        # contained hyperedges
        for i in list(H):
            c = H[i]
            if not c:
                del H[i]
                changed = True
                continue
            v = next(iter(c))
            for j in inc[v]:
                if j != i and c <= H[j]:
                    for w in c:
                        inc[w].discard(i)
                    del H[i]
                    changed = True
                    break
    return [frozenset(c) for c in H.values()]


def components(hyperedges):
    G = nx.Graph()
    for e in hyperedges:
        e = list(e)
        G.add_nodes_from(e)
        G.add_edges_from(zip(e[:-1], e[1:]))
    return [len(c) for c in nx.connected_components(G)]


def merge_closure(cliques):
    """Ch. 16: merge complexes sharing >= 2 atoms, to convergence; sizes."""
    parent = list(range(len(cliques)))

    def find(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    # pairs sharing an edge: index cliques by their edges
    by_edge = defaultdict(list)
    for i, c in enumerate(cliques):
        c = sorted(c)
        for a in range(len(c)):
            for b in range(a + 1, len(c)):
                by_edge[(c[a], c[b])].append(i)
    for lst in by_edge.values():
        for i in lst[1:]:
            ra, rb = find(lst[0]), find(i)
            if ra != rb:
                parent[rb] = ra
    # after one round the merged complexes may share >= 2 atoms with others:
    # iterate on the merged families until stable
    while True:
        groups = defaultdict(set)
        for i, c in enumerate(cliques):
            groups[find(i)] |= set(c)
        keys = list(groups)
        by_edge = defaultdict(list)
        merged_any = False
        # a cheap superset: two groups sharing >= 2 atoms
        atom_to_groups = defaultdict(list)
        for k in keys:
            for v in groups[k]:
                atom_to_groups[v].append(k)
        pair_count = defaultdict(int)
        for v, ks in atom_to_groups.items():
            for a in range(len(ks)):
                for b in range(a + 1, len(ks)):
                    pair_count[(ks[a], ks[b])] += 1
        for (a, b), cnt in pair_count.items():
            if cnt >= 2:
                ra, rb = find(a), find(b)
                if ra != rb:
                    parent[rb] = ra
                    merged_any = True
        if not merged_any:
            return sorted((len(s) for s in groups.values()), reverse=True)


def local_chordality(G, radius):
    """Vertices whose ball of the given radius is not chordal.

    At radius one a cycle through the centre always has a chord, so the ball
    is non-chordal iff the neighbourhood N(v) alone holds a chordless cycle:
    v is the hub of a wheel with four or more spokes, a loop of complexes
    around one atom.  This is not Ch. 14's criterion (two complexes sharing
    two atoms), which is stronger and is what the atom-message recursion
    needs; a wheel is what the join tree cannot absorb.
    """
    bad = []
    for v in G:
        ball = nx.single_source_shortest_path_length(G, v, cutoff=radius)
        if len(ball) > 2 and not nx.is_chordal(G.subgraph(ball)):
            bad.append(v)
    return list(G), bad


def on_short_chordless_cycle(G, adj=None):
    """Vertices lying on a chordless cycle of length four or five.

    alpha-acyclicity is an instance property and GYO keeps every cycle,
    long ones included, which the cavity method does not mind.  The
    ensemble notion is local: a *short* chordless cycle is a loop between
    complexes that no complex contains, the chygraph analogue of the short
    loop the tree assumption forbids.  Four and five are the two shortest
    lengths and are found from each vertex in time quadratic in its degree.
    """
    adj = adj or {v: set(G[v]) for v in G}
    bad = set()
    for v in G:
        if v in bad:
            pass  # still scan: a cycle through v marks its other vertices too
        nv = adj[v]
        closed = nv | {v}
        nbrs = sorted(nv)
        for i in range(len(nbrs)):
            a = nbrs[i]
            na = adj[a]
            for j in range(i + 1, len(nbrs)):
                b = nbrs[j]
                if b in na:
                    continue
                nb = adj[b]
                common = (na & nb) - closed
                if common:
                    bad.update((v, a, b))
                    bad.update(common)
                    continue
                X = na - closed - nb
                Y = nb - closed - na
                if not X or not Y:
                    continue
                for x in X:
                    hit = adj[x] & Y
                    if hit:
                        bad.update((v, a, b, x))
                        bad.update(hit)
                        break
    return bad


def measure(G, rng=None):
    rng = rng or np.random.default_rng(0)
    n = G.number_of_nodes()
    cliques = [c for c in nx.find_cliques(G) if len(c) >= 2]
    res = gyo_residue(cliques, n)
    comp = components(res)
    resid = sum(comp)
    out = dict(
        n=n, m=G.number_of_edges(),
        chordal=int(nx.is_chordal(G)),
        ncliques=len(cliques),
        maxclique=max((len(c) for c in cliques), default=0),
        residue=resid / n,
        giant=(max(comp) / n) if comp else 0.0,
        island=(sum(c * c for c in comp) / resid) if resid else 0.0,
        ncomp=len(comp),
    )
    for name, bad in (('wheel', local_chordality(G, 1)[1]), ('short', on_short_chordless_cycle(G))):
        out[f'phi_{name}'] = len(bad) / n
        comp = [len(c) for c in nx.connected_components(G.subgraph(bad))] if bad else []
        out[f'giant_{name}'] = (max(comp) / n) if comp else 0.0
        out[f'island_{name}'] = (sum(c * c for c in comp) / sum(comp)) if comp else 0.0
    sizes = merge_closure(cliques)
    out['merged'] = sizes[0] / n if sizes else 0.0
    out['merged_mean'] = (sum(s * s for s in sizes) / sum(sizes)) if sizes else 0.0
    return out


def _existing(path):
    """Rows already in the CSV, so an interrupted sweep resumes."""
    if not path.exists():
        return []
    return list(csv.DictReader(open(path)))


def _done(rows, **key):
    return any(all(float(r[k]) == float(v) for k, v in key.items()) for r in rows)


def run(ensemble, n, seeds, params, kbars):
    path = OUT / f'acyclicity_{ensemble}.csv'
    rows = _existing(path)
    for p in params:
        for kbar in kbars:
            for seed in seeds:
                if _done(rows, param=p, kbar=kbar, seed=seed, n=n):
                    continue
                rng = np.random.default_rng(seed)
                t0 = time.time()
                G = rgg_torus(n, p, kbar, rng) if ensemble == 'rgg' else hrg(n, p, kbar, rng)
                row = dict(ensemble=ensemble, param=p, kbar=kbar, seed=seed)
                row['kbar_meas'] = 2 * G.number_of_edges() / n
                row.update(measure(G))
                row['sec'] = round(time.time() - t0, 1)
                rows.append(row)
                print(' '.join(f'{k}={v:.4g}' if isinstance(v, float) else f'{k}={v}'
                               for k, v in row.items()), flush=True)
                with open(path, 'w', newline='') as f:
                    w = csv.DictWriter(f, fieldnames=list(rows[0]))
                    w.writeheader()
                    w.writerows(rows)


if __name__ == '__main__':
    which = sys.argv[1] if len(sys.argv) > 1 else "all"
    n = int(sys.argv[2]) if len(sys.argv) > 2 else 20000
    kbars = (0.5, 1.0, 1.5, 2.0, 3.0, 4.0, 6.0, 8.0, 12.0)
    if which in ('rgg', 'all'):
        run('rgg', n, (1, 2), (1,), (0.5, 1.0, 2.0, 4.0, 12.0))   # chordal: phi = 0 identically
        run('rgg', n, (1, 2), (2, 3), kbars)
    if which in ('hrg', 'all'):
        run('hrg', n, (1,), (2.5, 3.0, 4.0), (1.0, 2.0, 4.0, 8.0))


def fss():
    """Finite-size sweep: does the residue's giant component survive n -> inf?"""
    path = OUT / 'acyclicity_fss.csv'
    rows = _existing(path)
    jobs = [('rgg', 2, k, n) for k in (4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 12.0) for n in (5000, 20000, 80000)]
    jobs += [('rgg', 3, k, n) for k in (2.7, 3.5, 4.5, 6.0, 8.0, 10.0) for n in (5000, 20000, 80000)]
    jobs += [('hrg', 2.5, 4.0, n) for n in (5000, 20000)]
    for ens, p, kbar, n in jobs:
        for seed in (1, 2):
            if _done(rows, param=p, kbar=kbar, seed=seed, n=n):
                continue
            rng = np.random.default_rng(seed)
            t0 = time.time()
            G = rgg_torus(n, p, kbar, rng) if ens == 'rgg' else hrg(n, p, kbar, rng)
            row = dict(ensemble=ens, param=p, kbar=kbar, seed=seed)
            row['kbar_meas'] = 2 * G.number_of_edges() / n
            row['gc'] = max(len(c) for c in nx.connected_components(G)) / n
            row.update(measure(G))
            row['sec'] = round(time.time() - t0, 1)
            rows.append(row)
            print(' '.join(f'{k}={v:.4g}' if isinstance(v, float) else f'{k}={v}'
                           for k, v in row.items()), flush=True)
            with open(path, 'w', newline='') as f:
                w = csv.DictWriter(f, fieldnames=list(rows[0]))
                w.writeheader()
                w.writerows(rows)


if __name__ == '__main__' and len(sys.argv) > 1 and sys.argv[1] == 'fss':
    fss()
