"""Read probe/results/clique_moments.csv and decide TODO item 1.

Fits ``y ~ n^beta`` per (family, tau, kbar) for the quantities a chygraph needs
to exist:

    m2    second moment of the maximal-clique size distribution
    sbar  <c^2>/<c> - 1, the excess cardinality seen from a member
    c_max the clique number, whose exponent is known independently:
          Friedrich & Krohmer (2015) give Theta(n^{(3-tau)/2}) for 2 < tau < 3,
          so 0.45 / 0.25 / 0.05 at tau = 2.1 / 2.5 / 2.9.  That is the check
          that the measurement is reading the right thing.

Decision rule.  beta > 0 for the HRG and beta ~ 0 for the degree-matched
configuration model means the maximal-clique ensemble has no n-independent
limit, and the clique route to prediction 4 is dead.
"""

import csv
from collections import defaultdict
from pathlib import Path

import numpy as np

CSV = Path(__file__).parent / 'results' / 'clique_moments.csv'
FIT_FROM = 10_000          # fit the asymptotic range only


def load():
    rows = []
    with CSV.open() as fh:
        for r in csv.DictReader(fh):
            for k in ('tau', 'kbar_target', 'kbar', 'm1', 'm2', 'm3', 'sbar'):
                r[k] = float(r[k])
            for k in ('n', 'seed', 'c_max', 'degeneracy', 'n_cliques'):
                r[k] = int(r[k])
            rows.append(r)
    return rows


def fit(ns, ys):
    """Least-squares slope of log y against log n."""
    ns, ys = np.asarray(ns, float), np.asarray(ys, float)
    ok = (ns > 0) & (ys > 0)
    if ok.sum() < 3:
        return float('nan')
    return float(np.polyfit(np.log(ns[ok]), np.log(ys[ok]), 1)[0])


def main():
    rows = load()
    grouped = defaultdict(lambda: defaultdict(list))
    for r in rows:
        grouped[(r['family'], r['tau'], r['kbar_target'])][r['n']].append(r)

    print(f"Power-law exponents beta in y ~ n^beta, fitted for n >= {FIT_FROM}.")
    print("theory: c_max ~ n^((3-tau)/2) = 0.45 / 0.25 / 0.05 at tau = "
          "2.1 / 2.5 / 2.9  [Friedrich & Krohmer 2015]\n")
    hdr = (f"{'family':>7}{'tau':>6}{'kbar':>6}" +
           "".join(f"{k:>10}" for k in ('c_max', 'm2', 'm3', 'sbar')) +
           f"{'c_max@max n':>13}{'sbar@max n':>12}")
    print(hdr)
    print('-' * len(hdr))
    for key in sorted(grouped):
        family, tau, kbar = key
        byn = grouped[key]
        ns = sorted(n for n in byn if n >= FIT_FROM)
        allns = sorted(byn)
        if len(ns) < 3:
            continue
        mean = lambda n, k: float(np.mean([r[k] for r in byn[n]]))
        out = f"{family:>7}{tau:>6.1f}{kbar:>6.1f}"
        for k in ('c_max', 'm2', 'm3', 'sbar'):
            out += f"{fit(ns, [mean(n, k) for n in ns]):>10.3f}"
        out += (f"{mean(allns[-1], 'c_max'):>13.1f}"
                f"{mean(allns[-1], 'sbar'):>12.2f}")
        print(out)

    # -- the decisive statistic ------------------------------------------
    # kbar calibration drifts at heavy tails, so the n-trend *within* a family
    # is confounded by degree.  The HRG and its control share a degree sequence
    # at every (n, seed), so their ratio holds P(k) and assortativity fixed and
    # leaves clustering as the only difference.
    print("\nPaired ratio sbar_hrg / sbar_config, identical degree sequence.")
    ns_all = sorted({r['n'] for r in rows})
    print(f"{'tau':>5}{'kbar':>6}  " + "".join(f"{n:>9}" for n in ns_all)
          + f"{'beta':>8}")
    for tau in sorted({r['tau'] for r in rows}):
        for kbar in sorted({r['kbar_target'] for r in rows}):
            h, c = grouped[('hrg', tau, kbar)], grouped[('config', tau, kbar)]
            if not h or not c:
                continue
            ns, ratios, cells = sorted(h), [], ""
            for n in ns:
                pairs = [a['sbar'] / b[0]['sbar']
                         for a in h[n]
                         for b in [[x for x in c[n] if x['seed'] == a['seed']]]
                         if b and b[0]['sbar'] > 0]
                v = float(np.mean(pairs)) if pairs else float('nan')
                ratios.append(v)
                cells += f"{v:>9.2f}"
            big = [(n, v) for n, v in zip(ns, ratios) if n >= FIT_FROM]
            b = fit([n for n, _ in big], [v for _, v in big])
            print(f"{tau:>5.1f}{kbar:>6.0f}  {cells}{b:>8.3f}")

    print("\nMeasured kbar, HRG (calibration drifts at heavy tails; this is")
    print("why the paired ratio, not the within-family trend, is the statistic):")
    for tau in sorted({r['tau'] for r in rows}):
        for kbar in sorted({r['kbar_target'] for r in rows}):
            byn = grouped[('hrg', tau, kbar)]
            if byn:
                traj = " ".join(f"{np.mean([x['kbar'] for x in byn[n]]):>7.2f}"
                                for n in sorted(byn))
                print(f"  tau={tau} target={kbar:.0f}: {traj}")

    print("\nRaw trajectories of sbar = <c^2>/<c> - 1 (what the tensor needs):")
    for tau in sorted({r['tau'] for r in rows}):
        for kbar in sorted({r['kbar_target'] for r in rows}):
            for family in ('hrg', 'config'):
                byn = grouped[(family, tau, kbar)]
                if not byn:
                    continue
                traj = " ".join(
                    f"{np.mean([r['sbar'] for r in byn[n]]):>7.2f}"
                    for n in sorted(byn))
                print(f"  tau={tau} kbar={kbar:.0f} {family:>6}: {traj}")
    print(f"  (n = {', '.join(str(n) for n in sorted({r['n'] for r in rows}))})")


if __name__ == '__main__':
    main()
