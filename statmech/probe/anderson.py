"""Anderson localization and the Localization Landscape on a chygraph.

Tonetti, Cugliandolo and Tarzia (arXiv:2512.04037) solve two problems on
the K+1 = 3 Bethe lattice with one cavity recursion for the resolvent: the
Anderson mobility edge, and the percolation of the Localization Landscape
u = (H - E_min)^{-1} 1, whose threshold the LLT proposes as the mobility
edge.  Their Eqs. (3)-(9) are the book's two steps at cardinality two with
a Green's function for a message.  This probe carries the same recursion
to a chygraph whose complexes are cliques -- edges and triangles here --
where the interior sum of a complex is a Schur complement.  The message
type itself -- the down step `transmit`, the instance solvers, the clean
band and the incidence tree's branching -- is `statmech.resolvent`; this
probe holds only what runs it on a random ensemble: the population
dynamics of the self-energies, the landscape and the LLT percolation
message, and the scans behind the chapter's figures.

Ensembles: regular chy-degrees, every site in s single edges and t
triangles ("Bethe with triangles"); (s, t) = (3, 0) is the paper's
lattice.  Disorder eps ~ U[-W/2, W/2], hopping t = 1.

    python probe/anderson.py check     # instance checks against exact inversion
"""

import sys
import time
from pathlib import Path

import numpy as np

_ROOT = Path(__file__).resolve().parents[2]
for _p in (_ROOT / 'statmech' / 'src', _ROOT / 'percolation' / 'src'):
    if str(_p) not in sys.path:
        sys.path.insert(0, str(_p))
from statmech import resolvent as R  # noqa: E402

OUT = Path(__file__).parent / 'results'


# --------------------------------------------------------------- checks

def check_instances(n=600, W=3.0, seed=0):
    rng = np.random.default_rng(seed)
    print('instance checks: cavity with Schur-complement interior against exact inversion, n =', n)
    for (s, t_) in ((3, 0), (1, 1), (0, 2), (2, 1)):
        complexes = R.regular_instance(n, [2, 3], [s, t_], rng)
        eps = rng.uniform(-W / 2, W / 2, n)
        H = R.hamiltonian(n, complexes, eps)
        # (1) resolvent at complex z
        z = 0.7 - 0.2j
        Gex = np.diag(np.linalg.inv(H - z * np.eye(n)))
        G, _, _ = R.cavity_instance(n, complexes, eps, z, sweeps=150)
        err = np.abs(G - Gex)
        # (2) landscape at E_min below the spectrum
        emin = np.linalg.eigvalsh(H.real).min() - 0.05
        uex = np.linalg.solve(H.real - emin * np.eye(n), np.ones(n))
        u = R.landscape_instance(n, complexes, eps, emin, sweeps=150)
        rel = np.abs(u - uex) / np.abs(uex)
        print(f'  (s,t)=({s},{t_}): {len(complexes)} complexes; |G_cav - G_exact| median {np.median(err):.2e} '
              f'max {err.max():.2e} (|G| ~ {np.median(np.abs(Gex)):.2f}); landscape rel err median {np.median(rel):.2e} max {rel.max():.2e}')


if __name__ == '__main__':
    what = sys.argv[1] if len(sys.argv) > 1 else 'check'
    if what == 'check':
        check_instances()


# ------------------------------------------------------------- ensemble

class Ensemble:
    """Regular chygraph: every atom in s single edges and t triangles.

    Two populations of complex-to-atom self-energies, one per cardinality,
    each element carrying (Sigma, eta) and, for percolation, p.  Hopping
    t = 1 (the letter t is the triangle count here; hopping is `hop`).
    """

    def __init__(self, s, t, W, P=200000, seed=0, hop=1.0):
        self.s, self.t, self.W, self.P, self.hop = s, t, W, P, hop
        self.rng = np.random.default_rng(seed)

    # -- a cavity atom: its self-energy sum from (ne edges, nt triangles)
    def _cavity(self, E_re, T_re, ne, nt, n, extra=None):
        rng = self.rng
        tot = np.zeros(n, dtype=E_re.dtype)
        if ne:
            tot += E_re[rng.integers(0, self.P, (n, ne))].sum(1)
        if nt:
            tot += T_re[rng.integers(0, self.P, (n, nt))].sum(1)
        return tot

    def _eps(self, n):
        return self.rng.uniform(-self.W / 2, self.W / 2, n)

    # -- one synchronous generation of the real self-energy populations at energy E
    def real_step(self, E_re, T_re, E):
        s, t, hop, P = self.s, self.t, self.hop, self.P
        newE = newT = None
        if s:
            # edge (i,k): k's cavity has s-1 edges and t triangles
            cav = self._cavity(E_re, T_re, s - 1, t, P)
            d = self._eps(P) - E - cav
            newE = R.transmit(d[:, None], hop).sigma
        if t:
            # triangle (i,j1,j2): each j has s edges and t-1 triangles
            m11 = self._eps(P) - E - self._cavity(E_re, T_re, s, t - 1, P)
            m22 = self._eps(P) - E - self._cavity(E_re, T_re, s, t - 1, P)
            newT = R.transmit(np.stack([m11, m22], 1), hop).sigma
        return (newE if s else E_re), (newT if t else T_re)

    def real_population(self, E, sweeps=300, init=None):
        P = self.P
        E_re = np.zeros(P) if init is None else init[0]
        T_re = np.zeros(P) if init is None else init[1]
        for _ in range(sweeps):
            E_re, T_re = self.real_step(E_re, T_re, E)
        return E_re, T_re

    # -- band bottom of the pure hopping chygraph, from the uniform fixed point
    def band_bottom(self, lo=-12.0, hi=0.0, iters=60):
        """Largest E below which the uniform real fixed point exists (W = 0)."""
        return R.band_bottom([2, 3], [self.s, self.t], self.hop, lo, hi, iters)

    def E_min(self):
        """Bottom of the spectrum at disorder W: band bottom minus W/2, or the
        isolated uniform eigenvalue -hop*degree when that lies below it."""
        return R.spectral_bottom([2, 3], [self.s, self.t], self.W, self.hop)

    # -- landscape: joint (Sigma, eta) populations at E_min
    def landscape_step(self, E_re, T_re, E_eta, T_eta, Emin):
        s, t, hop, P, rng = self.s, self.t, self.hop, self.P, self.rng
        newE = newT = None
        newEe = newTe = None
        if s:
            ie = rng.integers(0, P, (P, s - 1)); it = rng.integers(0, P, (P, t))
            cav = E_re[ie].sum(1) + T_re[it].sum(1)
            eta_k = 1 + E_eta[ie].sum(1) + T_eta[it].sum(1)
            out = R.transmit((self._eps(P) - Emin - cav)[:, None], hop, eta_k[:, None])
            newE, newEe = out.sigma, out.eta
        if t:
            ie1 = rng.integers(0, P, (P, s)); it1 = rng.integers(0, P, (P, t - 1))
            ie2 = rng.integers(0, P, (P, s)); it2 = rng.integers(0, P, (P, t - 1))
            m11 = self._eps(P) - Emin - E_re[ie1].sum(1) - T_re[it1].sum(1)
            m22 = self._eps(P) - Emin - E_re[ie2].sum(1) - T_re[it2].sum(1)
            e1 = 1 + E_eta[ie1].sum(1) + T_eta[it1].sum(1)
            e2 = 1 + E_eta[ie2].sum(1) + T_eta[it2].sum(1)
            out = R.transmit(np.stack([m11, m22], 1), hop, np.stack([e1, e2], 1))
            newT, newTe = out.sigma, out.eta
        return (newE if s else E_re), (newT if t else T_re), (newEe if s else E_eta), (newTe if t else T_eta)

    def landscape_population(self, Emin, sweeps=300):
        P = self.P
        E_re, T_re = self.real_population(Emin, sweeps)
        E_eta, T_eta = np.zeros(P), np.zeros(P)
        for _ in range(sweeps):
            E_re, T_re, E_eta, T_eta = self.landscape_step(E_re, T_re, E_eta, T_eta, Emin)
        return E_re, T_re, E_eta, T_eta

    def landscape_sample(self, pops, Emin, n):
        """Full-site landscape u_i = G_ii eta_i for n atoms."""
        E_re, T_re, E_eta, T_eta = pops
        rng, s, t, P = self.rng, self.s, self.t, self.P
        ie = rng.integers(0, P, (n, s)); it = rng.integers(0, P, (n, t))
        g = 1.0 / (self._eps(n) - Emin - E_re[ie].sum(1) - T_re[it].sum(1))
        return g * (1 + E_eta[ie].sum(1) + T_eta[it].sum(1))

    # -- LLT percolation: linearised growth of the open-cluster message at E' = E - E_min
    def percolation_growth(self, pops, Emin, Eprime, sweeps=200, transient=100, joint=True):
        """Growth factor per generation of the percolation message p at
        landscape threshold u >= 1/E', with p kept small (linearised).
        Above 1: a giant cluster.  Every population element is a triple
        (Sigma, eta, p) generated from one cavity draw, so the correlation
        between a site's landscape and the messages it receives -- the
        smoothness of u along the graph -- is carried exactly."""
        E_re, T_re, E_eta, T_eta = [np.array(x) for x in pops]
        s, t, hop, P, rng = self.s, self.t, self.hop, self.P, self.rng
        thr = 1.0 / Eprime
        E_p = np.full(P, 1e-3) if s else np.zeros(0)
        T_p = np.full(P, 1e-3) if t else np.zeros(0)
        logs = []
        for it_ in range(sweeps):
            old = (E_p.sum() if s else 0.0) + (T_p.sum() if t else 0.0)
            if s:
                # edge (i, k): k's cavity has s-1 edges and t triangles; the edge's own
                # message towards k is an independent draw standing for i's side
                ie = rng.integers(0, P, (P, s - 1)); itr = rng.integers(0, P, (P, t))
                cav = E_re[ie].sum(1) + T_re[itr].sum(1)
                eta_cav = E_eta[ie].sum(1) + T_eta[itr].sum(1)
                eps_k = self._eps(P)
                out = R.transmit((eps_k - Emin - cav)[:, None], hop, (1 + eta_cav)[:, None])
                nE_re, nE_eta = out.sigma, out.eta
                back = rng.integers(0, P, P)
                u_k = (1 + eta_cav + E_eta[back]) / (eps_k - Emin - cav - E_re[back])
                if not joint:      # p from draws unrelated to the ones that made u_k
                    ie = rng.integers(0, P, (P, s - 1)); itr = rng.integers(0, P, (P, t))
                reach = 1 - np.prod(1 - E_p[ie], axis=1) * np.prod(1 - T_p[itr], axis=1)
                nE_p = (u_k >= thr) * reach
            if t:
                idx = [(rng.integers(0, P, (P, s)), rng.integers(0, P, (P, t - 1))) for _ in range(3)]
                eps3 = [self._eps(P) for _ in range(3)]
                cav = [E_re[ie].sum(1) + T_re[itr].sum(1) for ie, itr in idx]
                eta = [E_eta[ie].sum(1) + T_eta[itr].sum(1) for ie, itr in idx]
                m = [eps3[q] - Emin - cav[q] for q in range(3)]          # atoms i, j1, j2

                def msg(o1, o2):
                    out = R.transmit(np.stack([m[o1], m[o2]], 1), hop,
                                     np.stack([1 + eta[o1], 1 + eta[o2]], 1))
                    return out.sigma, out.eta
                nT_re, nT_eta = msg(1, 2)                                # the message to i
                p_j = []
                for j, o1, o2 in ((1, 0, 2), (2, 0, 1)):
                    sig, et = msg(o1, o2)                                # the message to j
                    u_j = (1 + eta[j] + et) / (m[j] - sig)
                    ie, itr = idx[j]
                    if not joint:
                        ie = rng.integers(0, P, (P, s)); itr = rng.integers(0, P, (P, t - 1))
                    reach = 1 - np.prod(1 - E_p[ie], axis=1) * np.prod(1 - T_p[itr], axis=1)
                    p_j.append((u_j >= thr) * reach)
                nT_p = 1 - (1 - p_j[0]) * (1 - p_j[1])
            new = (nE_p.sum() if s else 0.0) + (nT_p.sum() if t else 0.0)
            if not new > 0:
                # extinction in a finite population: subcritical; report what was measured
                logs.append(np.log(1e-3))
                break
            if it_ >= transient:
                logs.append(np.log(new / old))
            scale = 1e-3 * P * ((1 if s else 0) + (1 if t else 0)) / new
            if s:
                E_re, E_eta, E_p = nE_re, nE_eta, np.minimum(nE_p * scale, 1.0)
            if t:
                T_re, T_eta, T_p = nT_re, nT_eta, np.minimum(nT_p * scale, 1.0)
        return float(np.exp(np.mean(logs))) if logs else 0.0

    def percolation_growth_short(self, pops, Emin, Eprime, sweeps=8, transient=3, joint=True):
        """The growth factor from a few sweeps starting at a uniform p.

        Renormalised power iteration goes extinct in the subcritical phase:
        the support of p shrinks geometrically in a finite population.  A
        few sweeps from a uniform start cannot, and give the same crossing;
        this is what the figure of growth against E' uses on both sides of
        the threshold."""
        return self.percolation_growth(pops, Emin, Eprime, sweeps=sweeps, transient=transient, joint=joint)

    def percolation_threshold(self, pops, Emin, lo, hi, iters=12, **kw):
        """E'_c where the growth factor crosses one, by bisection."""
        for _ in range(iters):
            mid = 0.5 * (lo + hi)
            if self.percolation_growth(pops, Emin, mid, **kw) > 1:
                hi = mid
            else:
                lo = mid
        return 0.5 * (lo + hi)

    # -- Anderson: growth rate of the beta-moment of the linearised imaginary parts
    def imaginary_growth(self, E, betas=(0.3, 0.4, 0.5, 0.6, 0.7), sweeps=400, transient=200, init=None):
        """lambda(beta) for the linear map Im Sigma_{a->i} = hop^2 sum_j v_j^2 Im Sigma_{J->a},
        v = M^{-1} 1, the real parts sampled from their alpha = 0 population at energy E.
        The localised phase is stable iff min_beta lambda(beta) < 1."""
        s, t, hop, P, rng = self.s, self.t, self.hop, self.P, self.rng
        E_re, T_re = self.real_population(E, sweeps=transient, init=init)
        E_im = rng.random(P) if s else np.zeros(0)
        T_im = rng.random(P) if t else np.zeros(0)
        acc = {b: [] for b in betas}
        for it_ in range(sweeps):
            # real parts keep evolving (stationary) so the weights are sampled from the fixed point
            newE = newT = None
            if s:
                ie = rng.integers(0, P, (P, s - 1)); itr = rng.integers(0, P, (P, t))
                out = R.transmit((self._eps(P) - E - E_re[ie].sum(1) - T_re[itr].sum(1))[:, None], hop)
                newE_re = out.sigma
                newE = R.linearised_imaginary(out.v, (E_im[ie].sum(1) + T_im[itr].sum(1))[:, None], hop)
            if t:
                ie1 = rng.integers(0, P, (P, s)); it1 = rng.integers(0, P, (P, t - 1))
                ie2 = rng.integers(0, P, (P, s)); it2 = rng.integers(0, P, (P, t - 1))
                m11 = self._eps(P) - E - E_re[ie1].sum(1) - T_re[it1].sum(1)
                m22 = self._eps(P) - E - E_re[ie2].sum(1) - T_re[it2].sum(1)
                out = R.transmit(np.stack([m11, m22], 1), hop)
                newT_re = out.sigma
                newT = R.linearised_imaginary(
                    out.v, np.stack([E_im[ie1].sum(1) + T_im[it1].sum(1),
                                     E_im[ie2].sum(1) + T_im[it2].sum(1)], 1), hop)
            oldv = np.concatenate([E_im, T_im])
            newv = np.concatenate([newE if s else np.zeros(0), newT if t else np.zeros(0)])
            if it_ >= transient:
                for b in betas:
                    acc[b].append(np.log(np.mean(newv ** b) / np.mean(oldv ** b)))
            scale = 1.0 / np.exp(np.mean(np.log(newv)))
            if s:
                E_re, E_im = newE_re, newE * scale
            if t:
                T_re, T_im = newT_re, newT * scale
        return {b: float(np.exp(np.mean(v))) for b, v in acc.items()}

    def mobility_edge(self, lo, hi, iters=10, **kw):
        """E_c at fixed W: bisection on min_beta lambda(beta) = 1, between an energy
        lo in the localised phase and hi in the delocalised one (or the reverse)."""
        f_lo = min(self.imaginary_growth(lo, **kw).values()) - 1
        for _ in range(iters):
            mid = 0.5 * (lo + hi)
            f = min(self.imaginary_growth(mid, **kw).values()) - 1
            if (f > 0) == (f_lo > 0):
                lo, f_lo = mid, f
            else:
                hi = mid
        return 0.5 * (lo + hi)


# ------------------------------------------------------------------ scan

ENSEMBLES = ((3, 0), (1, 1), (4, 0), (2, 1), (0, 2))
WS = (1.5, 3.0, 4.5, 6.0, 9.0, 12.0)


def _lines(args):
    """E_c^perc and E_c^loc at the bottom of the band for one (s, t, W)."""
    import csv
    s, t, W, P, seed = args
    t0 = time.time()
    # percolation: a classical problem, a smaller population suffices; the
    # threshold E'_c grows with W (u shrinks as 1/W), so the bracket scales with |E_min|
    e = Ensemble(s, t, W, P=min(P, 300000), seed=seed)
    Emin = e.E_min()
    pops = e.landscape_population(Emin, sweeps=300)
    Epc = e.percolation_threshold(pops, Emin, 0.02, 1.5 * abs(Emin), iters=11, sweeps=250, transient=100)
    # mobility edge: localised at the very bottom, extended at the band centre for W < W_c
    e = Ensemble(s, t, W, P=P, seed=seed)
    Ecl = e.mobility_edge(Emin + 0.02, 0.0, iters=11, betas=(0.4, 0.5, 0.6), sweeps=400, transient=200)
    row = dict(s=s, t=t, degree=s + 2 * t, W=W, P=P, seed=seed, Emin=Emin, band_bottom=e.band_bottom(),
               Eprime_c=Epc, Ec_perc=Emin + Epc, Ec_loc=Ecl, gap=(Emin + Epc) - Ecl, sec=round(time.time() - t0))
    print(' '.join(f'{k}={v:.5g}' if isinstance(v, float) else f'{k}={v}' for k, v in row.items()), flush=True)
    return row


def _wc(args):
    """Band-centre critical disorder for one (s, t)."""
    s, t, P, seed = args
    t0 = time.time()
    e = Ensemble(s, t, 10.0, P=P, seed=seed)

    def f(W):
        e.W = W
        return min(e.imaginary_growth(0.0, betas=(0.4, 0.5, 0.6), sweeps=400, transient=200).values()) - 1
    lo, hi = 5.0, 40.0
    for _ in range(10):
        mid = 0.5 * (lo + hi)
        if f(mid) > 0:
            lo = mid
        else:
            hi = mid
    row = dict(s=s, t=t, degree=s + 2 * t, P=P, seed=seed, Wc=0.5 * (lo + hi), sec=round(time.time() - t0))
    print(' '.join(f'{k}={v:.5g}' if isinstance(v, float) else f'{k}={v}' for k, v in row.items()), flush=True)
    return row


def scan(P=1000000, seed=0, procs=8):
    import csv
    from multiprocessing import Pool
    jobs = [(s, t, W, P, seed) for (s, t) in ENSEMBLES for W in WS]
    with Pool(procs) as pool:
        rows = pool.map(_lines, jobs)
    with open(OUT / 'anderson_lines.csv', 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    jobs = [(s, t, P, seed) for (s, t) in ENSEMBLES]
    with Pool(min(procs, len(jobs))) as pool:
        rows = pool.map(_wc, jobs)
    with open(OUT / 'anderson_wc.csv', 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)


if __name__ == '__main__' and len(sys.argv) > 1 and sys.argv[1] == 'scan':
    scan()


# --------------------------------------------------- figures' side data

def side_data(P=300000, seed=0):
    """Three small scans behind the chapter's mechanism figures.

    stability   lambda(1/2) against E at W = 3 for (3,0) and (1,1): the
                mobility edge as a crossing of one
    growth      the percolation growth factor against E' at W = 1.5 on
                (3,0), with the (Sigma, eta, p) triple drawn jointly and
                drawn separately: what the smoothness of u is worth
    centre      lambda(1/2) at E = 0 against W for the five ensembles:
                where each band centre localises
    """
    import csv
    rows = []
    for (s, t) in ((3, 0), (1, 1)):
        e = Ensemble(s, t, 3.0, P=P, seed=seed)
        Emin = e.E_min()
        for E in np.linspace(Emin + 0.1, Emin + 2.5, 13):
            lam = e.imaginary_growth(E, betas=(0.5,), sweeps=300, transient=150)[0.5]
            rows.append(dict(s=s, t=t, W=3.0, E=E, lam=lam))
            print(f'stability ({s},{t}) E={E:.3f} lambda={lam:.4f}', flush=True)
    with open(OUT / 'anderson_stability.csv', 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)

    rows = []
    e = Ensemble(3, 0, 1.5, P=P, seed=seed)
    Emin = e.E_min()
    pops = e.landscape_population(Emin, sweeps=300)
    for Ep in np.linspace(0.35, 0.65, 13):
        g_joint = e.percolation_growth_short(pops, Emin, Ep)
        g_indep = e.percolation_growth_short(pops, Emin, Ep, joint=False)
        rows.append(dict(s=3, t=0, W=1.5, Eprime=Ep, growth_joint=g_joint, growth_indep=g_indep))
        print(f"growth E'={Ep:.3f} joint={g_joint:.4f} indep={g_indep:.4f}", flush=True)
    with open(OUT / 'anderson_growth.csv', 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)

    rows = []
    for (s, t) in ENSEMBLES:
        e = Ensemble(s, t, 10.0, P=P, seed=seed)
        for W in (4.0, 6.0, 8.0, 10.0, 12.0, 15.0, 18.0, 22.0, 26.0, 30.0, 35.0):
            e.W = W
            lam = e.imaginary_growth(0.0, betas=(0.5,), sweeps=300, transient=150)[0.5]
            rows.append(dict(s=s, t=t, W=W, lam=lam))
            print(f'centre ({s},{t}) W={W} lambda={lam:.4f}', flush=True)
    with open(OUT / 'anderson_centre.csv', 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)


def side_growth(P=300000, seed=0):
    import csv
    rows = []
    e = Ensemble(3, 0, 1.5, P=P, seed=seed)
    Emin = e.E_min()
    pops = e.landscape_population(Emin, sweeps=300)
    for Ep in np.linspace(0.35, 0.65, 13):
        g_joint = e.percolation_growth_short(pops, Emin, Ep)
        g_indep = e.percolation_growth_short(pops, Emin, Ep, joint=False)
        rows.append(dict(s=3, t=0, W=1.5, Eprime=Ep, growth_joint=g_joint, growth_indep=g_indep))
        print(f"growth E'={Ep:.3f} joint={g_joint:.4f} indep={g_indep:.4f}", flush=True)
    with open(OUT / 'anderson_growth.csv', 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)


if __name__ == '__main__' and len(sys.argv) > 1 and sys.argv[1] == 'side':
    side_data()
if __name__ == '__main__' and len(sys.argv) > 1 and sys.argv[1] == 'growth':
    side_growth()
