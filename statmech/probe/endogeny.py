"""The endogeny test of Chapter 20, run on the spin glass and on the hitting set.

Everything the chapter quotes is measured here and cached in
``results/endogeny.json``; ``book/figs/endogeny.py`` draws it.

    python probe/endogeny.py scan      # the whole thing, about forty minutes

Three measurements.

spin glass
    Ising with couplings +-J on a random regular graph of degree three, in a
    field H.  For each H three boundaries in beta J are bisected: where the
    jointly carried linear perturbation stops decaying (the linearised
    endogeny test, which is the de Almeida--Thouless line computed with the
    perturbation carried alongside the field); where the annealed criterion
    <kbar> <(du/dh)^2> crosses one (Sec. 9.7's substitution evaluated on the
    non-trivial fixed point, which is not the same thing); and where two
    copies started as independent draws from the fixed-point law stop
    converging to each other (the nonlinear test, Theorem 11(c)).  At H = 0
    the exact transition is at 2 tanh^2(beta J) = 1.

hitting set, Poisson
    For cardinality c = 2, 3, 4 and a grid of mean chy-degrees around the
    hard-field line <k>(c-1) = e, the stationary rate of a small perturbation
    (long windows, since the transient is slow near the boundary) and the
    final distance of a shuffled pair.  The zero of the rate is the endogeny
    boundary; the entropy of Sec. 10.8 is evaluated at the same points.

hitting set, regular
    Mezard--Tarzia's regular hypergraphs, where the entropy is exact: the
    pair test against the sign of the entropy.
"""

import json
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'src'))
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'percolation' / 'src'))
from statmech.endogeny import BivariateHittingSet, BivariateSpinGlass  # noqa: E402
from statmech.hittingset import rsb_point  # noqa: E402
from statmech.softfield import regular_entropy  # noqa: E402

OUT = Path(__file__).parent / 'results' / 'endogeny.json'
SIZE = 100_000


def rate(hist, a, b):
    h = np.array(hist[a:b])
    h = h[h > 0]
    if len(h) < 10:
        return float('nan')
    return float(np.polyfit(np.arange(len(h)), np.log(h), 1)[0])


def bisect(f, lo, hi, n):
    for _ in range(n):
        mid = 0.5 * (lo + hi)
        lo, hi = (mid, hi) if f(mid) else (lo, mid)
    return 0.5 * (lo + hi)


# ------------------------------------------------------------- spin glass

def sg_linear(bJ, H, seed=1):
    m = BivariateSpinGlass(3, bJ, field=H, regular=True, size=SIZE, seed=seed)
    m.run_single(300)
    g = m.at_growth()
    m.run_pair(80, perturb=1e-4)
    return rate(m.history, 5, 60), g, m.overlap()


def sg_shuffled(bJ, H, seed=2):
    m = BivariateSpinGlass(3, bJ, field=H, regular=True, size=SIZE, seed=seed)
    m.run_single(300)
    m.run_pair(600)
    return m.history[-1], rate(m.history, 400, 600)


def scan_spin_glass(fields=(0.02, 0.05, 0.1, 0.2, 0.3, 0.6, 1.0)):
    rows = []
    for H in fields:
        t = time.time()
        b_lin = bisect(lambda b: sg_linear(b, H)[0] < 0, 0.6, 3.0, 10)
        b_ann = bisect(lambda b: sg_linear(b, H)[1] < 1, 0.6, 3.0, 10)
        b_nl = bisect(lambda b: sg_shuffled(b, H)[0] < 1e-2, 0.6, 3.0, 8)
        rows.append({'H': H, 'linear': b_lin, 'annealed': b_ann, 'shuffled': b_nl})
        print(f'  H={H}: linear {b_lin:.3f}  annealed {b_ann:.3f}  shuffled {b_nl:.3f}  '
              f'({time.time() - t:.0f}s)', flush=True)
    # the zero-field transition, for the anchor
    rows_zero = {'beta_J_c': float(np.arctanh(1 / np.sqrt(2)))}
    return {'rows': rows, 'zero_field': rows_zero}


# ------------------------------------------------------ hitting set, Poisson

def hs_point(c, k, seed=3, sweeps=600, damping=0.5):
    m = BivariateHittingSet([c], [k], mu=60.0, size=SIZE, seed=seed, damping=damping)
    m.run(400)
    m.run_pair(sweeps, perturb=0.06)
    h = m.history
    r1, r2 = rate(h, int(0.4 * sweeps), int(0.7 * sweeps)), rate(h, int(0.7 * sweeps), sweeps)
    s, se = m.entropy_averaged(keep=100)
    m2 = BivariateHittingSet([c], [k], mu=60.0, size=SIZE, seed=seed + 1, damping=damping)
    m2.run(400)
    m2.run_pair(600)
    return {'c': c, 'k': k, 'rate_a': r1, 'rate_b': r2, 'rate': 0.5 * (r1 + r2),
            'entropy': s, 'entropy_se': se, 'shuffled': m2.history[-1],
            'shuffled_rate': rate(m2.history, 400, 600), 'density': m.density()}


def scan_hitting_set(cards=(2, 3, 4), ratios=(0.85, 0.95, 1.0, 1.03, 1.07, 1.12, 1.2, 1.4)):
    rows = []
    for c in cards:
        hard = float(rsb_point([c], [1.0]))
        for r in ratios:
            t = time.time()
            row = hs_point(c, r * hard)
            row['hard'] = hard
            row['ratio'] = r
            rows.append(row)
            print(f"  c={c} <k>={r * hard:.4f} ({r:.2f} x hard): rate {row['rate']:+.5f} "
                  f"({row['rate_a']:+.5f}, {row['rate_b']:+.5f})  shuffled {row['shuffled']:.2e}  "
                  f"entropy {row['entropy']:+.3f}({row['entropy_se']:.3f})  ({time.time() - t:.0f}s)",
                  flush=True)
    return rows


# ------------------------------------------------------ hitting set, regular

def scan_regular(cases=((2, 3), (3, 3), (4, 3), (2, 4), (3, 4), (2, 6), (3, 6), (4, 6))):
    rows = []
    for L, K in cases:
        t = time.time()
        m = BivariateHittingSet([K], [L], mu=60.0, size=SIZE, seed=5, regular=True)
        m.run(300)
        m.run_pair(600)
        d = m.history[-1]
        m2 = BivariateHittingSet([K], [L], mu=60.0, size=SIZE, seed=6, regular=True)
        m2.run(300)
        # the rate while the perturbation is still small: a growing one
        # saturates within a hundred sweeps and a late window reads zero
        m2.run_pair(120, perturb=0.06)
        r = rate(m2.history, 10, 70)
        rows.append({'L': L, 'K': K, 'entropy': float(regular_entropy(L, K)),
                     'shuffled': d, 'rate': r})
        print(f'  L={L} K={K}: entropy {regular_entropy(L, K):+.3f}  shuffled {d:.2e}  '
              f'rate {r:+.4f}  ({time.time() - t:.0f}s)', flush=True)
    return rows


def main(only=None):
    t0 = time.time()
    out = json.load(open(OUT)) if (only and OUT.exists()) else {}
    if only == 'regular':
        print('regular hitting sets:')
        out['regular'] = scan_regular()
        json.dump(out, open(OUT, 'w'), indent=1)
        print(f'rewrote the regular block of {OUT} in {time.time() - t0:.0f} s')
        return
    print('regular hitting sets:')
    out['regular'] = scan_regular()
    print('spin glass boundaries:')
    out['spin_glass'] = scan_spin_glass()
    print('Poisson hitting sets:')
    out['hitting_set'] = scan_hitting_set()
    OUT.parent.mkdir(exist_ok=True)
    json.dump(out, open(OUT, 'w'), indent=1)
    print(f'wrote {OUT} in {time.time() - t0:.0f} s')


if __name__ == '__main__':
    if len(sys.argv) > 1 and sys.argv[1] == 'scan':
        main(sys.argv[2] if len(sys.argv) > 2 else None)
    else:
        print(__doc__)
