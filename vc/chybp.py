"""Question 2: complexes solved exactly, tree-like propagation between them.

Regions are a partition of the vertices (here: the bottom-level blocks of the
nested degree-corrected SBM, kept as a region when the block has at most
``smax`` vertices, split into singletons otherwise -- so ``smax=0`` is the
ordinary node-level cavity and ``smax=inf`` is "every community a complex").
Edges inside a region are solved exactly (a weighted vertex cover of the
region, by branch and bound on the complement's maximum-weight independent
set); edges between regions carry min-sum (zero-temperature) messages, which
is the tree-like propagation.  The message from region R to one of its
vertices i is the exact minimum cost of R with x_i fixed, incoming fields on
the other vertices included, so the interior is never approximated.

Factor-graph bookkeeping.  Factors are the regions (cost |cover ∩ R| plus the
internal covering constraints) and the external edges (constraint
x_i + x_j >= 1).  Factor-to-variable messages are normalised to min 0:

    m_{e->i}(1) = 0,  m_{e->i}(0) = h_{j->i} = max(0, n_{j->e}(1) - n_{j->e}(0))

so h is the push "cover i" that the other end exerts, and

    m_{R->i}(x_i) = min_{interior, x_i fixed} [ |cover| + sum_{u in R, u != i} [x_u = 0] H_u ]

with H_u the total push on u from all its external edges.  Variable-to-factor
messages are sums of the others, and the min-sum Bethe energy is

    E = sum_R eps_R + sum_e eps_e - sum_i (d_i - 1) eps_i .

Besides E, a valid cover is read off by decimation (each region's optimal
interior given its fields, then repair of any violated external edge), whose
size is an upper bound that can be checked against the exact optimum.
"""
from collections import defaultdict

import numpy as np


class BudgetExceeded(Exception):
    """The exact interior is too expensive: treat the block as singletons."""


# ------------------------------------------------------------ exact interior
class Region:
    """A vertex set with its internal edges, as bitmasks over local indices."""

    def __init__(self, nodes, internal_edges):
        self.nodes = list(nodes)
        self.index = {v: k for k, v in enumerate(self.nodes)}
        m = len(self.nodes)
        self.nbr = [0] * m
        for a, b in internal_edges:
            ia, ib = self.index[a], self.index[b]
            self.nbr[ia] |= 1 << ib
            self.nbr[ib] |= 1 << ia
        self.full = (1 << m) - 1

    budget = 200000          # max branch-and-bound calls per solve

    def mwis(self, w):
        """Max-weight independent set over vertices with w > 0.

        Returns (value, mask).  Vertices with w <= 0 are excluded up front,
        which for the cover means they are taken.  Raises BudgetExceeded if
        the search tree passes ``budget`` nodes (sparse, large, loopy
        interiors -- exactly the ones that are not dense blocks).
        """
        alive = 0
        for k, x in enumerate(w):
            if x > 0:
                alive |= 1 << k
        memo = {}
        nbr = self.nbr
        calls = [0]

        def solve(mask):
            if mask == 0:
                return 0.0, 0
            got = memo.get(mask)
            if got is not None:
                return got
            calls[0] += 1
            if calls[0] > self.budget:
                raise BudgetExceeded
            # split into connected components
            comps = []
            rest = mask
            while rest:
                seed = rest & -rest
                comp, frontier = 0, seed
                while frontier:
                    comp |= frontier
                    nxt = 0
                    f = frontier
                    while f:
                        b = f & -f
                        nxt |= nbr[b.bit_length() - 1]
                        f ^= b
                    frontier = nxt & mask & ~comp
                comps.append(comp)
                rest &= ~comp
            if len(comps) > 1:
                val, sel = 0.0, 0
                for c in comps:
                    v, s = solve(c)
                    val += v
                    sel |= s
                memo[mask] = (val, sel)
                return val, sel
            # one component: reductions, then branch on the max-degree vertex
            best_v, best_d = -1, -1
            f = mask
            while f:
                b = f & -f
                f ^= b
                v = b.bit_length() - 1
                d = bin(nbr[v] & mask).count('1')
                if d == 0:
                    val, sel = solve(mask ^ b)
                    memo[mask] = (val + w[v], sel | b)
                    return memo[mask]
                if d == 1:
                    u = (nbr[v] & mask).bit_length() - 1
                    if w[v] >= w[u]:
                        val, sel = solve(mask & ~(b | (1 << u)))
                        memo[mask] = (val + w[v], sel | b)
                        return memo[mask]
                if d > best_d:
                    best_v, best_d = v, d
            v = best_v
            b = 1 << v
            v_in, s_in = solve(mask & ~(nbr[v] | b))
            v_in += w[v]
            s_in |= b
            v_out, s_out = solve(mask ^ b)
            memo[mask] = (v_in, s_in) if v_in >= v_out else (v_out, s_out)
            return memo[mask]

        return solve(alive)

    def min_cover(self, H, fixed=None):
        """min over covers C of R of  |C| + sum_{u not in C} H[u],  with x fixed.

        ``fixed`` maps local index -> 0/1.  Returns (cost, cover mask).
        """
        m = len(self.nodes)
        w = [1.0 - H[k] for k in range(m)]          # cost of covering, relative
        forced = 0                                   # covered
        excluded = 0                                 # uncovered (fixed to 0)
        if fixed:
            for k, x in fixed.items():
                if x == 1:
                    forced |= 1 << k
                else:
                    excluded |= 1 << k
                    forced |= self.nbr[k]            # its neighbours must cover
            if forced & excluded:
                return np.inf, 0
        for k in range(m):
            if (forced >> k) & 1:
                w[k] = -1.0                          # take it, whatever it costs
            elif (excluded >> k) & 1:
                w[k] = np.inf                        # never take it
        val, sel = self.mwis([x if x != np.inf else 1e9 for x in w])
        # cover = everything not in the independent set; forced ones are in
        cover = self.full & ~sel
        base = sum(H)                                # all uncovered
        cost = base
        for k in range(m):
            if (cover >> k) & 1:
                cost += 1.0 - H[k]
        return cost, cover


BIG = 1e6


# --------------------------------------------------------------- the solver
class ChyBP:
    def __init__(self, n, edges, blocks, smax=64):
        self.n = n
        edges = [tuple(e) for e in np.asarray(edges).tolist()]
        # regions: blocks of size <= smax, singletons otherwise
        members = defaultdict(list)
        for v, b in enumerate(blocks):
            members[int(b)].append(v)
        region_of = np.empty(n, dtype=int)
        groups = []
        for b, vs in members.items():
            if len(vs) <= smax:
                region_of[vs] = len(groups)
                groups.append(vs)
            else:
                for v in vs:
                    region_of[v] = len(groups)
                    groups.append([v])
        self.region_of = region_of
        internal = defaultdict(list)
        self.ext = []                                  # external edges (i, j)
        for a, b in edges:
            if region_of[a] == region_of[b]:
                internal[region_of[a]].append((a, b))
            else:
                self.ext.append((a, b))
        # trial-solve every multi-vertex region; too expensive -> singletons
        regions, new_groups = [], []
        region_of = np.empty(n, dtype=int)
        self.n_split = 0
        for k, vs in enumerate(groups):
            R = Region(vs, internal[k])
            ok = True
            if len(vs) > 1:
                try:
                    R.min_cover([0.0] * len(vs))
                except BudgetExceeded:
                    ok = False
            if ok:
                region_of[vs] = len(regions)
                regions.append(R)
                new_groups.append(vs)
            else:
                self.n_split += 1
                for v in vs:
                    region_of[v] = len(regions)
                    regions.append(Region([v], []))
                    new_groups.append([v])
        self.region_of = region_of
        groups = new_groups
        self.ext = []
        for a, b in edges:
            if region_of[a] != region_of[b]:
                self.ext.append((a, b))
        self.regions = regions
        self.n_regions = len(groups)
        self.n_singletons = sum(len(g) == 1 for g in groups)
        # per node: list of external edge ids
        self.ext_of = defaultdict(list)
        for eid, (a, b) in enumerate(self.ext):
            self.ext_of[a].append(eid)
            self.ext_of[b].append(eid)
        # message h[eid][side]: push exerted on endpoint `side` by the other
        self.ext_of = dict(self.ext_of)               # plain: no phantom keys
        self.clamp = {}                                # node -> 0/1, decimation
        E = len(self.ext)
        self.h = np.zeros((E, 2))
        self.ext_arr = np.array(self.ext, dtype=int).reshape(E, 2)

    # total push on node v from all its external edges
    def _H(self, v):
        tot = 0.0
        for eid in self.ext_of.get(v, ()):
            side = 0 if self.ext[eid][0] == v else 1
            tot += self.h[eid, side]
        return tot

    def _single_dm(self, v):
        """m_{R->v}(0) - m_{R->v}(1) for a singleton region: -1, or a clamp."""
        c = self.clamp.get(v)
        return -1.0 if c is None else (BIG if c == 1 else -BIG)

    def _region_messages(self, k, Hloc):
        """m_{R->i}(0) - m_{R->i}(1) for every i in region k with external edges."""
        R = self.regions[k]
        out = {}
        if len(R.nodes) == 1:
            out[R.nodes[0]] = self._single_dm(R.nodes[0])
            return out
        fixed = {R.index[v]: x for v, x in self.clamp.items() if v in R.index}
        for v in R.nodes:
            if v not in self.ext_of:
                continue
            i = R.index[v]
            if v in self.clamp:
                out[v] = BIG if self.clamp[v] == 1 else -BIG
                continue
            c0, _ = R.min_cover(Hloc, {**fixed, i: 0})
            c1, _ = R.min_cover(Hloc, {**fixed, i: 1})
            # m_{R->i}(x_i) must exclude i's own field H_i: min_cover put H_i in
            # the x_i = 0 branch, take it back out
            out[v] = (c0 - Hloc[i]) - c1
        return out

    def iterate(self, iters=300, damping=0.5, tol=1e-7, verbose=False):
        for it in range(iters):
            # H_u for all nodes with external edges
            Htot = {v: self._H(v) for v in self.ext_of}
            # region messages dm[v] = m_{R->v}(0) - m_{R->v}(1)
            dm = {}
            for k, R in enumerate(self.regions):
                if len(R.nodes) == 1:
                    dm[R.nodes[0]] = self._single_dm(R.nodes[0])
                    continue
                Hloc = [Htot.get(v, 0.0) for v in R.nodes]
                dm.update(self._region_messages(k, Hloc))
            # new pushes: for edge e=(a,b), n_{a->e}(0)-n_{a->e}(1) = dm[a] + (H_a - h_{b->a})
            new = np.empty_like(self.h)
            for eid, (a, b) in enumerate(self.ext):
                da = dm[a] + Htot[a] - self.h[eid, 0]
                db = dm[b] + Htot[b] - self.h[eid, 1]
                new[eid, 1] = max(0.0, -da)            # push on b from a
                new[eid, 0] = max(0.0, -db)            # push on a from b
            diff = np.abs(new - self.h).max() if len(self.h) else 0.0
            self.h = damping * self.h + (1 - damping) * new
            if verbose and it % 10 == 0:
                print(f'  it {it:4d}  max change {diff:.2e}', flush=True)
            if diff < tol:
                return it + 1, diff
        return iters, diff

    def bethe_energy(self):
        Htot = {v: self._H(v) for v in self.ext_of}
        E = 0.0
        dm = {}
        for k, R in enumerate(self.regions):
            if len(R.nodes) == 1:
                v = R.nodes[0]
                E += min(Htot.get(v, 0.0), 1.0)         # eps_R
                dm[v] = self._single_dm(v)
                continue
            Hloc = [Htot.get(v, 0.0) for v in R.nodes]
            c, _ = R.min_cover(Hloc)
            E += c
            dm.update(self._region_messages(k, Hloc))
        for eid, (a, b) in enumerate(self.ext):
            # n_{a->e}(x): normalised so that min over x is 0
            # n_{a->e}(x) = m_{R->a}(x) + [x=0](H_a - h_{b->a}), m_R normalised
            na = (max(dm[a], 0.0) + Htot[a] - self.h[eid, 0], max(-dm[a], 0.0))
            nb = (max(dm[b], 0.0) + Htot[b] - self.h[eid, 1], max(-dm[b], 0.0))
            E += min(na[0] + nb[1], na[1] + nb[0], na[1] + nb[1])
        for v, eids in self.ext_of.items():
            d = len(eids) + 1
            b0 = max(dm[v], 0.0) + Htot[v]            # m_R(0) + H
            b1 = max(-dm[v], 0.0)
            E -= (d - 1) * min(b0, b1)
        return E

    def beliefs(self):
        """b_i(0) - b_i(1): >0 means i prefers to be covered.  Nodes without
        external edges get the preference from their region's optimum."""
        Htot = {v: self._H(v) for v in self.ext_of}
        pref = np.zeros(self.n)
        for k, R in enumerate(self.regions):
            if len(R.nodes) == 1:
                v = R.nodes[0]
                pref[v] = Htot.get(v, 0.0) - 1.0
                continue
            Hloc = [Htot.get(v, 0.0) for v in R.nodes]
            fixed = {R.index[v]: x for v, x in self.clamp.items() if v in R.index}
            c, mask = R.min_cover(Hloc, fixed or None)
            for i, v in enumerate(R.nodes):
                if v in self.clamp:
                    pref[v] = BIG if self.clamp[v] == 1 else -BIG
                    continue
                other = 0 if (mask >> i) & 1 else 1
                c_other, _ = R.min_cover(Hloc, {**fixed, i: other})
                pref[v] = (c_other - c) if other == 0 else -(c_other - c)
        return pref

    def _adjacency(self):
        adj = defaultdict(list)
        for a, b in self.ext:
            adj[a].append(b)
            adj[b].append(a)
        for R in self.regions:
            for x in range(len(R.nodes)):
                for y in range(x):
                    if (R.nbr[x] >> y) & 1:
                        adj[R.nodes[x]].append(R.nodes[y])
                        adj[R.nodes[y]].append(R.nodes[x])
        return adj

    def decimate(self, frac=0.1, iters=25, damping=0.5, tie=0.5):
        """Sequential decimation.  Decided nodes (|b_0 - b_1| >= tie) are
        clamped in batches; when only ties are left, one node per connected
        component of the tie subgraph is clamped to `covered`, which breaks
        the degeneracy of a matching edge, a cycle or a clique in one round.  Then BP is re-run.  Ends with a repair and a prune."""
        adj = self._adjacency()
        free = set(range(self.n))
        while free:
            pref = self.beliefs()
            strong = [v for v in free if abs(pref[v]) >= tie]
            if strong:
                strong.sort(key=lambda v: -abs(pref[v]))
                batch = strong[:max(1, int(frac * len(strong)))]
            else:
                ties = [v for v in free]
                tie_set = set(ties)
                seen, batch = set(), []
                for v in ties:                     # one per tie component
                    if v in seen:
                        continue
                    batch.append(v)
                    stack = [v]
                    seen.add(v)
                    while stack:
                        u = stack.pop()
                        for w in adj[u]:
                            if w in tie_set and w not in seen:
                                seen.add(w)
                                stack.append(w)
            for v in batch:
                # ties go to `uncovered`: that forces every neighbour and so
                # settles a clique, a cycle or a matching edge in one round
                x = (1 if pref[v] >= 0 else 0) if strong else 0
                if x == 0 and any(self.clamp.get(u) == 0 for u in adj[v]):
                    continue                       # never two uncovered neighbours
                self.clamp[v] = x
                free.discard(v)
            if not any(v in self.clamp for v in batch):
                v = batch[0]                       # all skipped: cover the first
                self.clamp[v] = 1
                free.discard(v)
            if free:
                self.iterate(iters=iters, damping=damping)
        cover = np.array([self.clamp[v] == 1 for v in range(self.n)])
        self.clamp = {}
        # repair every edge with both ends uncovered
        all_edges = list(self.ext) + [(R.nodes[x], R.nodes[y]) for R in self.regions
                                      for x in range(len(R.nodes)) for y in range(x)
                                      if (R.nbr[x] >> y) & 1]
        viol = [(a, b) for a, b in all_edges if not cover[a] and not cover[b]]
        while viol:
            cnt = defaultdict(int)
            for a, b in viol:
                cnt[a] += 1
                cnt[b] += 1
            v = max(cnt, key=cnt.get)
            cover[v] = True
            viol = [(a, b) for a, b in viol if not cover[a] and not cover[b]]
        # prune: uncover any covered vertex whose neighbours are all covered
        changed = True
        while changed:
            changed = False
            for v in np.flatnonzero(cover):
                if all(cover[u] for u in adj[v]):
                    cover[v] = False
                    changed = True
        return int(cover.sum()), cover
