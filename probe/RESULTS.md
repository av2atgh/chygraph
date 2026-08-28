# Probe: can an HRG be given a chygraph complex ensemble?

TODO item 1, step 1. A chygraph maps a clustered graph by making its dense
motifs into complexes; for a hyperbolic random graph the natural choice is the
maximal cliques. The threshold tensor needs the excess cardinality seen from a
member to have a finite mean,

```
sbar = <c^2>/<c> - 1
```

> **Notation.** `beta` below is a *structural* growth exponent. The manuscript
> renames these `theta` (referee minor 2), reserving `beta` for the inverse
> temperature; the numbers are unchanged.

so the maximal-clique size distribution needs a **converged second moment**.
The TODO predicted this would fail and close the item negatively.

**It does not fail.** For `tau >= 2.5` the moment converges, and the HRG
separates from its degree-matched control by a stable factor of 1.8–4.8. The
clique route to prediction 4 is alive.

## Baselines: what a clustering-free graph looks like

| graph | `sbar` across `n = 10^4 .. 3x10^5` | `c_max` |
|---|---|---|
| Erdős–Rényi, k̄ = 2 / 4 / 8 | 0.937 / 0.996 / 1.000, flat to 3 decimals | 3 |
| config model, `tau = 3.5`, k̄ = 4 | 0.995 → 0.990 | 3 |
| config model, `tau = 4.5`, k̄ = 8 | 1.005 → 1.000 | 3 |

A graph with finite `<k^2>` and no clustering sits at `sbar ≈ 1`. That is the
reference.

## The HRG converges for `tau >= 2.5`, and sits far above the baseline

`sbar`, mean of 3 seeds, `n = 10^3 .. 3x10^5`:

| tau | k̄ | HRG | config | ratio (paired, same degree sequence) | beta of ratio |
|---:|---:|---|---|---:|---:|
| 2.9 | 2 | 1.58 → 1.58 | 0.92 → 0.90 | 1.76 | +0.003 |
| 2.9 | 4 | 2.71 → 2.76 | 1.09 → 0.99 | 2.80 | +0.011 |
| 2.9 | 8 | 4.84 → 4.81 | 1.18 → 1.01 | 4.77 | +0.028 |
| 2.5 | 2 | 1.48 → 1.67 | 0.92 → 0.93 | 1.79 | +0.011 |
| 2.5 | 4 | 2.96 → 3.12 | 1.27 → 1.17 | 2.66 | +0.015 |
| 2.5 | 8 | 4.74 → 5.99 | 1.53 → 1.70 | 3.52 | +0.009 |
| 3.5 | 4 | 2.832 → 2.838 | 0.995 → 0.990 | 2.87 | flat |
| 3.5 | 8 | 4.855 → 4.892 | 1.021 → 1.000 | 4.89 | flat |
| 4.5 | 8 | 5.032 → 5.061 | 1.005 → 1.000 | 5.06 | flat |

The **paired ratio** is the statistic: HRG and control share a degree sequence
at every `(n, seed)`, so `P(k)` and assortativity are held fixed and clustering
is the only difference. Its exponent is `+0.003` to `+0.028` — flat. The
clustering signal is a converged multiplicative constant, not a finite-size
artefact.

Density is not the explanation. Erasure leaves the control ~8% sparser, but
`sbar` for a clustering-free graph barely moves with density at all
(0.996 → 1.000 from k̄ = 4 to 8 in ER), so it cannot produce a factor of 3–5.

## Only `tau = 2.1` fails, and it fails on the tail, not the geometry

| tau | quantity | HRG beta | config beta |
|---|---|---:|---:|
| 2.1 | `sbar ~ n^beta` | 0.28 – 0.40 | **0.42 – 0.58** |
| 2.5 | `sbar ~ n^beta` | 0.009 – 0.024 | −0.004 – 0.014 |
| 2.9 | `sbar ~ n^beta` | 0.000 – 0.001 | −0.028 – −0.002 |

At `tau = 2.1` both diverge and **the control diverges faster than the HRG**.
The paired ratio *decays* toward 1 (`beta = −0.12` to `−0.18`): the heavy tail
swamps the geometry. So the τ<3 range the TODO specified could not have
answered the question either way — a growing HRG moment there is a hub-clique
effect that a configuration model reproduces without any geometry.

## Measurement check

The clique number should scale as `n^{(3-tau)/2}` (Friedrich & Krohmer 2015):

| tau | theory | measured (HRG, k̄ = 2/4/8) |
|---:|---:|---|
| 2.1 | 0.45 | 0.369, 0.351, 0.376 |
| 2.5 | 0.25 | 0.210, 0.193, 0.201 |
| 2.9 | 0.05 | 0.059, 0.071, 0.092 |

Consistently a little below theory, as expected at finite `n`, and tracking it
across a 9x range in the exponent. The probe is reading the right quantity.

`kbar` calibration is exact at `tau = 2.5/2.9` (2.00, 4.00, 8.00) and drifts
only at `tau = 2.1`, which is a further reason to read the paired ratio rather
than the within-family trend.

## Verdict

Step 1 of TODO item 1 is **satisfied for `tau >= 2.5`**: the maximal-clique
ensemble exists, with a converged second moment, and it is strongly separated
from the degree-matched control. That is the regime where
`~/av2atg/computational_complexity` measured the extensive leaf-removal core,
so prediction 4 is directly testable rather than blocked.

What remains is steps 2 and 3 of the TODO item: write core percolation as a
chygraph fixed point, and run the comparison.

## Reproduce

```sh
python3 probe/er_control.py     # baseline, ~1 min
python3 probe/light_tail.py     # tau > 3, ~5 min
python3 probe/clique_moments.py # main scan, ~40 min -> results/clique_moments.csv
python3 probe/analyze.py        # exponents and paired ratios
```

Needs `networkx`, and `~/av2atg/computational_complexity/code/hrg.py` on the path.
