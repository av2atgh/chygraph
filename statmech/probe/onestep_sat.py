"""One-step replica symmetry breaking at m != 0 for random 3-SAT: the
entropic (zero-energy) calculation, as a population of populations.

Sec. 13.9 does the one-step calculation at m = 0, where every cluster is
counted once and the message is a warning.  This probe carries the same
recursion at general m, where the message is a real number -- the BP
message of the solution-counting problem -- and a survey is a population
of them, reweighted by z^m at every update.  It returns the one-step
potential F(m), the internal entropy s(m) = dF/dm and the complexity
Sigma(m) = F - m s, from which: the condensation point alpha_c, where
Sigma(1) reaches zero; the Parisi parameter m*(alpha) in the condensed
phase, where Sigma(m*) = 0; and the entropy of solutions.

Conventions (Mezard & Montanari, Chs. 14, 19, 20).  Clause-to-variable
message eta_{a->i} in [0,1]: the probability that the other variables of
a all take their unsatisfying values; as a distribution over x_i it is
[1 - eta d(x_i, unsat)] / (2 - eta).  Variable-to-clause message u_{j->a}
in [0,1]: the probability that x_j takes the value unsatisfying a,
u = Pi_s / (Pi_s + Pi_o) with Pi_s the product of (1 - eta_b) over the
other clauses b in which j appears with the same sign as in a, Pi_o the
same over the opposite sign.  Then eta_{a->i} = prod_{j != i} u_{j->a}.
The normalisations that carry the reweighting are

    z_{j->a} = (Pi_s + Pi_o) / prod_b (2 - eta_b)        variable update
    z_{a->i} = 2 - eta_{a->i}                              clause update

and the free entropy is F = <ln z_i> + alpha <ln z_a> - k alpha <ln z_ia>,
with z_i = (Pi_+ + Pi_-)/prod_b(2 - eta_b) over all of a variable's
clauses, z_a = 1 - prod_{j in a} u_j and z_ia = (1 - u eta)/(2 - eta).
At m = 1 the one-step potential equals the replica-symmetric free entropy;
at m -> 0 the complexity would count clusters, but with soft messages every
combination is non-contradictory and Sigma(0) = 0 identically: the m = 0
count of Sec. 13.9 needs the hard fields that exist only above the rigidity
point (about 4.25 for 3-SAT), so the anchors used are the m = 1 identity and
the published condensation point 3.86.

    python probe/onestep_sat.py check            # RS vs enumeration, m=1 vs RS, m=0 vs alpha_s
    python probe/onestep_sat.py scan             # F, s, Sigma on a grid of (alpha, m)

Writes probe/results/onestep_sat.csv.
"""

import csv
import sys
import time
from multiprocessing import Pool
from pathlib import Path

import numpy as np

OUT = Path(__file__).parent / 'results' / 'onestep_sat.csv'
K = 3


# ------------------------------------------------------------------ RS

def rs_population(alpha, size=100000, sweeps=200, seed=0, k=K):
    """Replica-symmetric population of clause-to-variable messages eta."""
    rng = np.random.default_rng(seed)
    eta = rng.random(size)
    lam = k * alpha / 2
    for _ in range(sweeps):
        # each of the k-1 variables: d_s, d_o other clauses
        u = np.ones(size)
        for _j in range(k - 1):
            ds = rng.poisson(lam, size)
            do = rng.poisson(lam, size)
            u *= _u_from(eta, ds, do, rng)
        eta = u
    return eta


def _u_from(eta, ds, do, rng):
    """u = Pi_s/(Pi_s+Pi_o) for arrays of same/opposite-sign degrees."""
    n = len(ds)
    lps = _log_prod_1m(eta, ds, rng)
    lpo = _log_prod_1m(eta, do, rng)
    return 1.0 / (1.0 + np.exp(lpo - lps))


def _log_prod_1m(eta, d, rng):
    """log prod of (1 - eta) over d draws from the population, per row."""
    n = len(d)
    out = np.zeros(n)
    dmax = int(d.max()) if n else 0
    for r in range(dmax):
        rows = np.nonzero(d > r)[0]
        out[rows] += np.log1p(-eta[rng.integers(0, len(eta), len(rows))])
    return out


def rs_free_entropy(eta, alpha, nsamp=200000, seed=1, k=K):
    """F_RS = <ln z_i> + alpha <ln z_a> - k alpha <ln z_ia>."""
    rng = np.random.default_rng(seed)
    lam = k * alpha / 2
    # variable term
    dp = rng.poisson(lam, nsamp)
    dm = rng.poisson(lam, nsamp)
    lp = _log_prod_1m(eta, dp, rng)
    lm = _log_prod_1m(eta, dm, rng)
    # log prod (2 - eta) over all dp + dm clauses
    l2 = np.zeros(nsamp)
    for arr in (dp, dm):
        for r in range(int(arr.max())):
            rows = np.nonzero(arr > r)[0]
            l2[rows] += np.log(2 - eta[rng.integers(0, len(eta), len(rows))])
    fi = np.logaddexp(lp, lm) - l2
    # clause term: k independent u's
    logprod = np.zeros(nsamp)
    for _ in range(k):
        u = _u_from(eta, rng.poisson(lam, nsamp), rng.poisson(lam, nsamp), rng)
        logprod += np.log(u)
    fa = np.log1p(-np.exp(logprod))
    # edge term
    u = _u_from(eta, rng.poisson(lam, nsamp), rng.poisson(lam, nsamp), rng)
    e = eta[rng.integers(0, len(eta), nsamp)]
    fia = np.log1p(-u * e) - np.log(2 - e)
    return fi.mean() + alpha * fa.mean() - k * alpha * fia.mean()


# ------------------------------------------------------------------ 1RSB

class OneStep:
    """Population of M surveys, each a population of P messages eta."""

    def __init__(self, alpha, m, M=400, P=400, seed=0, k=K):
        self.alpha, self.m, self.M, self.P, self.k = alpha, m, M, P, k
        self.rng = np.random.default_rng(seed)
        self.lam = k * alpha / 2
        self.pop = self.rng.random((M, P))

    # -- one variable-to-clause survey, reweighted: returns (u[P], log z[P])
    def _u_survey(self):
        rng = self.rng
        ds, do = rng.poisson(self.lam), rng.poisson(self.lam)
        P = self.P
        lps = np.zeros(P)
        lpo = np.zeros(P)
        l2 = np.zeros(P)
        for _ in range(ds):
            e = self.pop[rng.integers(self.M)][rng.permutation(P)]
            lps += np.log1p(-e)
            l2 += np.log(2 - e)
        for _ in range(do):
            e = self.pop[rng.integers(self.M)][rng.permutation(P)]
            lpo += np.log1p(-e)
            l2 += np.log(2 - e)
        u = 1.0 / (1.0 + np.exp(lpo - lps))
        logz = np.logaddexp(lps, lpo) - l2
        return u, logz

    def _resample(self, values, logw):
        """Importance resampling of P values with weights exp(m logw).

        A combination with z = 0 -- a contradiction -- gets weight zero at
        every m, m = 0 included: 0^0 is 0 here, which is what discards the
        contradictory state in survey propagation.
        """
        finite = np.isfinite(logw)
        w = np.where(finite, self.m * np.where(finite, logw, 0.0), -np.inf)
        w = np.exp(w - w[finite].max())
        w /= w.sum()
        idx = self.rng.choice(self.P, self.P, p=w)
        return values[idx]

    def _eta_survey(self):
        """One clause-to-variable survey: product of k-1 reweighted u's, then z_a^m."""
        logeta = np.zeros(self.P)
        for _ in range(self.k - 1):
            u, logz = self._u_survey()
            u = self._resample(u, logz)
            logeta += np.log(u)
        eta = np.minimum(np.exp(logeta), 1 - 1e-12)
        return self._resample(eta, np.log(2 - eta))

    def sweep(self):
        for _ in range(self.M):
            self.pop[self.rng.integers(self.M)] = self._eta_survey()

    def run(self, sweeps=300):
        for _ in range(sweeps):
            self.sweep()
        return self

    # -- the potential and its m-derivative
    def _term(self, logz):
        """(ln <z^m>, <z^m ln z>/<z^m>) over an inner population of log z."""
        w = self.m * logz
        c = w.max()
        zm = np.exp(w - c)
        lnavg = c + np.log(zm.mean())
        s = (zm * logz).sum() / zm.sum()
        return lnavg, s

    # -- vectorised draws: rows x P messages from random surveys
    def _draw(self, rows):
        rng = self.rng
        surv = rng.integers(0, self.M, rows)
        cols = rng.integers(0, self.P, (rows, self.R))
        return self.pop[surv[:, None], cols]

    def _u_block(self, C):
        """C candidate u-surveys as (u, logz), each (C, P), *unweighted*:
        the reweighting z^m is applied as an importance weight downstream."""
        rng = self.rng
        ds, do = rng.poisson(self.lam, C), rng.poisson(self.lam, C)
        lps, lpo, l2 = np.zeros((C, self.R)), np.zeros((C, self.R)), np.zeros((C, self.R))
        for d, acc in ((ds, lps), (do, lpo)):
            for r in range(int(d.max()) if C else 0):
                rows = np.nonzero(d > r)[0]
                e = self._draw(len(rows))
                acc[rows] += np.log1p(-e)
                l2[rows] += np.log(2 - e)
        u = 1.0 / (1.0 + np.exp(lpo - lps))
        return u, np.logaddexp(lps, lpo) - l2

    @staticmethod
    def _wterm(m, logz, logw=None):
        """Per row: ln <z^m>_w and <z^m ln z>_w / <z^m>_w, weights exp(logw)."""
        lw = m * logz if logw is None else m * logz + logw
        bad = ~np.isfinite(lw)
        lw = np.where(bad, -np.inf, lw)
        c = lw.max(axis=1, keepdims=True)
        dead = ~np.isfinite(c[:, 0])            # every combination contradictory
        c = np.where(np.isfinite(c), c, 0.0)
        w = np.exp(lw - c)
        tot = w.sum(axis=1)
        lnavg = c[:, 0] + np.log(np.where(dead, 1.0, tot) / logz.shape[1])
        lnavg = np.where(dead, np.nan, lnavg)
        if logw is not None:                       # normalise the importance weights
            lw0 = np.where(np.isfinite(logw), logw, -np.inf)
            c0 = lw0.max(axis=1, keepdims=True)
            lnavg -= c0[:, 0] + np.log(np.exp(lw0 - c0).sum(axis=1) / logz.shape[1])
        s = (w * np.where(bad, 0.0, logz)).sum(axis=1) / np.where(dead, 1.0, tot)
        return lnavg, np.where(dead, np.nan, s)

    def potential(self, nsamp=200000, chunk=1000, reps=5):
        """F(m), s(m) = dF/dm and Sigma = F - m s, by importance weighting.

        Each sample combines R = reps * P draws (with replacement) from the
        surveys involved: ln <z^m> over a survey combination is estimated
        from R combinations, and the Jensen bias of the logarithm of a noisy
        average falls with R.  reps=1 uses P combinations.
        """
        rng = self.rng
        self.R = reps * self.P
        m, k, alpha = self.m, self.k, self.alpha
        F = np.zeros(3)
        S = np.zeros(3)
        cnt = np.zeros(3)
        dead = np.zeros(3)
        done = 0

        def _acc(i, a, b):
            ok = np.isfinite(a)
            F[i] += a[ok].sum()
            S[i] += b[ok].sum()
            cnt[i] += ok.sum()
            dead[i] += (~ok).sum()
        while done < nsamp:
            C = min(chunk, nsamp - done)
            # variable term
            dp, dm = rng.poisson(self.lam, C), rng.poisson(self.lam, C)
            lp, lm, l2 = np.zeros((C, self.R)), np.zeros((C, self.R)), np.zeros((C, self.R))
            for d, acc in ((dp, lp), (dm, lm)):
                for r in range(int(d.max())):
                    rows = np.nonzero(d > r)[0]
                    e = self._draw(len(rows))
                    acc[rows] += np.log1p(-e)
                    l2[rows] += np.log(2 - e)
            a, b = self._wterm(m, np.logaddexp(lp, lm) - l2)
            _acc(0, a, b)
            # clause term: k candidate u-surveys, weights prod z_j^m
            logprod = np.zeros((C, self.R))
            logw = np.zeros((C, self.R))
            for _ in range(k):
                u, logz = self._u_block(C)
                logprod += np.log(u)
                logw += m * logz
            a, b = self._wterm(m, np.log1p(-np.exp(logprod)), logw)
            _acc(1, a, b)
            # edge term
            u, logz = self._u_block(C)
            e = self._draw(C)
            a, b = self._wterm(m, np.log1p(-u * e) - np.log(2 - e), m * logz)
            _acc(2, a, b)
            done += C
        F /= np.maximum(cnt, 1)
        S /= np.maximum(cnt, 1)
        Fm = F[0] + alpha * F[1] - k * alpha * F[2]
        sm = S[0] + alpha * S[1] - k * alpha * S[2]
        # dead: fraction of draws in which every combination was contradictory,
        # by term; the potential is only meaningful where this is negligible
        return dict(F=Fm, s=sm, Sigma=Fm - m * sm,
                    spread=float(np.mean(self.pop.std(axis=1))),
                    dead=float(dead.sum() / (dead.sum() + cnt.sum())))


# ------------------------------------------------------------------ checks

def bp_on_instance(n, clauses, signs, sweeps=500, damping=0.3, seed=0):
    """BP on one instance; returns the Bethe free entropy ln Z_Bethe."""
    rng = np.random.default_rng(seed)
    m = len(clauses)
    # clause a, position p: variable clauses[a,p], sign signs[a,p] (True: satisfied by x=1)
    eta = rng.random((m, K)) * 0.5
    # incidence lists per variable
    inc = [[] for _ in range(n)]
    for a in range(m):
        for p in range(K):
            inc[clauses[a, p]].append((a, p))
    for _ in range(sweeps):
        new = np.empty_like(eta)
        for a in range(m):
            us = []
            for p in range(K):
                j, sg = clauses[a, p], signs[a, p]
                ls = lo = 0.0
                for (b, q) in inc[j]:
                    if b == a:
                        continue
                    if signs[b, q] == sg:
                        ls += np.log1p(-eta[b, q])
                    else:
                        lo += np.log1p(-eta[b, q])
                us.append(1.0 / (1.0 + np.exp(lo - ls)))
            for p in range(K):
                new[a, p] = np.prod([us[q] for q in range(K) if q != p])
        eta = damping * eta + (1 - damping) * new
    # free entropy
    F = 0.0
    for i in range(n):
        lp = lm = l2 = 0.0
        for (b, q) in inc[i]:
            if signs[b, q]:      # satisfied by x=1: unsat value is 0 -> contributes to x=0 branch
                lp += np.log1p(-eta[b, q])
            else:
                lm += np.log1p(-eta[b, q])
            l2 += np.log(2 - eta[b, q])
        F += np.logaddexp(lp, lm) - l2
    for a in range(m):
        us = []
        for p in range(K):
            j, sg = clauses[a, p], signs[a, p]
            ls = lo = 0.0
            for (b, q) in inc[j]:
                if b == a:
                    continue
                if signs[b, q] == sg:
                    ls += np.log1p(-eta[b, q])
                else:
                    lo += np.log1p(-eta[b, q])
            us.append(1.0 / (1.0 + np.exp(lo - ls)))
        F += np.log1p(-np.prod(us))
        for p in range(K):
            F -= np.log1p(-us[p] * eta[a, p]) - np.log(2 - eta[a, p])
    return F


def check_bp_exact_on_tree():
    """A tree-shaped 3-SAT instance: BP's Bethe free entropy equals ln(#solutions)."""
    from cluster_blocks import solutions
    rng = np.random.default_rng(3)
    # chain of clauses sharing one variable each: clause a uses variables 2a, 2a+1, 2a+2
    m = 7
    n = 2 * m + 1
    clauses = np.array([[2 * a, 2 * a + 1, 2 * a + 2] for a in range(m)])
    signs = rng.integers(0, 2, (m, K)).astype(bool)
    exact = np.log(solutions(n, clauses, signs).shape[0])
    bethe = bp_on_instance(n, clauses, signs)
    print(f'tree instance: ln #solutions = {exact:.6f}, Bethe = {bethe:.6f}')
    assert abs(exact - bethe) < 1e-6


def check_rs_against_enumeration(alpha=2.0, n=20, trials=10):
    """Loopy small instances: Bethe is close to, not equal to, the exact count."""
    from cluster_blocks import random_3sat, solutions
    errs = []
    for t in range(trials):
        rng = np.random.default_rng(t)
        clauses, signs = random_3sat(n, alpha, rng)
        S = solutions(n, clauses, signs).shape[0]
        if S == 0:
            continue
        errs.append((bp_on_instance(n, clauses, signs) - np.log(S)) / n)
    print(f'random N={n} alpha={alpha}: (Bethe - exact)/N over {len(errs)} instances: '
          f'mean {np.mean(errs):+.4f}, max |.| {np.max(np.abs(errs)):.4f}')


def check_m1_equals_rs(alpha=3.5, seed=0):
    eta = rs_population(alpha, seed=seed)
    frs = rs_free_entropy(eta, alpha)
    one = OneStep(alpha, 1.0, M=300, P=300, seed=seed).run(150).potential(800)
    print(f'alpha={alpha}: F_RS = {frs:.5f}; one-step at m=1: F = {one["F"]:.5f}, '
          f'Sigma = {one["Sigma"]:+.5f}, survey spread {one["spread"]:.4f}')
    return frs, one


def check_m0_threshold(alphas=(4.15, 4.25, 4.35), seed=0):
    out = []
    for a in alphas:
        r = OneStep(a, 0.0, M=300, P=300, seed=seed).run(200).potential(800)
        print(f'alpha={a}: Sigma(m=0) = {r["Sigma"]:+.5f} (F = {r["F"]:.5f}, s = {r["s"]:.5f})')
        out.append((a, r['Sigma']))
    return out


# ------------------------------------------------------------------ scan

def _point(args):
    alpha, m, seed, M, P, sweeps, nsamp = args
    t0 = time.time()
    r = OneStep(alpha, m, M=M, P=P, seed=seed).run(sweeps).potential(nsamp, reps=2)
    r.update(alpha=alpha, m=m, seed=seed, sec=round(time.time() - t0, 1))
    print(' '.join(f'{k}={v:.5g}' if isinstance(v, float) else f'{k}={v}' for k, v in r.items()), flush=True)
    return r


def scan(alphas=(3.7, 3.8, 3.86, 3.9, 3.95, 4.0, 4.05, 4.1, 4.15, 4.2, 4.25),
         ms=(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0),
         seeds=(0, 1), M=1000, P=200, sweeps=300, nsamp=400000, procs=8, out=OUT):
    # resume: skip the points already in the log of a previous, interrupted run
    done = set()
    log = out.with_name(out.stem + '_scan.log') if out == OUT else out.with_suffix('.log')
    if log.exists():
        for line in open(log):
            if line.startswith('F='):
                r = dict(tok.split('=') for tok in line.split())
                done.add((float(r['alpha']), float(r['m']), int(r['seed'])))
    jobs = [(a, m, s, M, P, sweeps, nsamp) for a in alphas for m in ms for s in seeds
            if (a, m, s) not in done]
    with Pool(procs) as pool:
        rows = pool.map(_point, jobs) if jobs else []
    rows += [dict((k, float(v) if k != 'seed' else int(v)) for k, v in
                  (tok.split('=') for tok in line.split()))
             for line in open(log) if line.startswith('F=')] if log.exists() and done else []
    with open(out, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)


def log_to_csv(log=OUT.parent / 'onestep_sat_scan.log'):
    """Recover the CSV from the per-point log lines if a scan was interrupted."""
    rows = []
    for line in open(log):
        if not line.startswith('F='):
            continue
        r = {}
        for tok in line.split():
            k, v = tok.split('=')
            r[k] = float(v) if k not in ('seed',) else int(v)
        rows.append(r)
    if rows:
        with open(OUT, 'w', newline='') as f:
            w = csv.DictWriter(f, fieldnames=list(rows[0]))
            w.writeheader()
            w.writerows(rows)
    return len(rows)


if __name__ == '__main__':
    what = sys.argv[1] if len(sys.argv) > 1 else 'check'
    if what == 'log2csv':
        print(log_to_csv(), 'rows')
    if what == 'check':
        check_bp_exact_on_tree()
        check_rs_against_enumeration()
        check_m1_equals_rs()
        check_m0_threshold()
    elif what == 'scan':
        scan()
    elif what == 'extra':
        # more populations at the points that decide m*(alpha): seeds 2..5
        scan(alphas=(4.0, 4.05, 4.1, 4.15, 4.2), ms=(0.4, 0.6, 0.8), seeds=(2, 3, 4),
             out=OUT.parent / 'onestep_sat_extra.csv')
