# Survey propagation on nested clauses, against the CNF flattening

**Status: unresolved.** Two of three parts work. The third does not produce a
number, and until it does this is not a result and is not in the book.

`nested_sp.py`. Self-contained (numpy only); nothing in `book/figs/` imports it
and nothing in the book depends on it.

## What it is for

Secs. 13.6–13.7 are the only place in the book where a complex contains another
complex — layer-2 clauses `D = NOT C OR l_1 ... OR l_k2` whose antecedent `C` is
itself a layer-1 clause. The book compares the nested formulation against its
CNF flattening at the level of the **warning** map, and finds the folds differ
by 7–22%, growing with `k1`, with an exact identity at `k2 = 1`.

Sec. 13.9 now does survey propagation and reaches `alpha_s` rather than
bracketing it. The obvious next experiment is the same nested-vs-flattened
comparison in that currency: a threshold on each side instead of a fold on each
side. That is what this file was written for, and what it does not yet deliver.

## What works

**Cross-level message passing.** At `k1 = 1` the sub-clause holds one literal,
`NOT C` is a single literal, and the three-level formula *is* flat
`(k2+1)`-SAT. The book establishes that identity at warning level exactly; here
both sides are population dynamics with independent draws, so the test is
whether the discrepancy is sampling noise. It is — clean `1/sqrt(N)`:

| population | mean \|nested − flat\| | ratio |
|---|---|---|
| 20 000 | 2.750e-03 | |
| 80 000 | 1.370e-03 | 2.01 |
| 320 000 | 6.988e-04 | 1.96 |

against 2.00 expected for a fourfold increase. No plateau, so no residual
asymmetry in how the two channels reaching a variable are merged.

**Cluster counting.** The hypothesis is that a layer-1 complex, being a
sub-expression rather than a constraint, contributes **no factor of its own** to
the count — its interior is already summed into the message its parent
receives. Then per variable

```
Sigma = <ln Z_i> + alpha <ln Z_D> - k2 alpha <ln Z_iD> - k1 alpha <ln Z_iC>
```

with two inclusion channels. At `k1 = 1` this must reduce to flat `(k2+1)`-SAT,
which requires both that the channels merge to exactly `(k2+1) alpha` and that
`Z_D = 1 - [prod_k2 pi] gamma` collapses to `1 - prod_{k2+1} pi`. Neither is
arranged. It reduces, again as noise:

| population | mean \|Sigma_nested − Sigma_flat\| | ratio |
|---|---|---|
| 20 000 | 6.036e-03 | |
| 80 000 | 3.507e-03 | 1.72 |
| 320 000 | 1.816e-03 | 1.93 |

So the counting hypothesis stands. **This is the part worth keeping**: it is the
expensive result in this file, and it licenses comparing complexities rather
than only branch onsets.

## What does not work

The nested system has no usable `Sigma = 0`. At `k1 = k2 = 2`:

| alpha | Sigma nested | ⟨eta⟩ | Sigma flat |
|---|---|---|---|
| 2.10 | +0.00000 | 0.0000 | +0.00472 |
| 2.15 | +0.04070 | 0.2049 | −0.00321 |
| 2.20 | +0.04327 | 0.2413 | −0.02448 |
| 2.35 | +0.04877 | 0.3005 | −0.04052 |
| 2.60 | +1.18 | 0.3406 | (saturated) |

The flattened side behaves and self-validates: it crosses zero between 2.10 and
2.15, which is 4.20–4.30 in its own clause density, bracketing the published
`alpha_s = 4.267`.

The nested side does not. Its branch appears **discontinuously** at a *finite*
`Sigma ≈ 0.041` and then stays flat rather than descending through zero, and by
`alpha ≈ 2.5` the `delta` channel saturates (`delta = prod pi -> 1` makes
`ln(1 - delta pi)` diverge) and the expression stops meaning anything. There is
no zero in the window between onset and breakdown.

Two candidate explanations, **not distinguished**, and separating them is the
open work:

1. **`m = 0` is the wrong parameter for this system.** A branch arriving
   discontinuously at finite `Sigma` is the signature of a transition that
   counting clusters at equal weight does not locate; that needs `Sigma(m)` with
   the Parisi parameter optimised. This is the more likely of the two.
2. **The saturation masks a crossing** that sits slightly higher and would be
   reachable with a better-conditioned `delta` update.

## Also incomplete

- **Only `(k1,k2) = (2,2)` was ever scanned.** `(3,2)` and `(2,3)` failed on
  badly chosen brackets and were not retried. This matters: `k1`-dependence is
  the interesting part, since Sec. 13.7's fold gap *grows* with `k1`.
- Two failure modes to avoid on any retry, both hit here. `Sigma` returns
  exactly `0.0` on the trivial branch, so a bisection bracket whose lower end
  sits below the branch onset is rejected rather than bracketing; the lower end
  must satisfy onset < lo < crossing. And `Sigma` is only meaningful near the
  crossing — evaluated well above it the logs saturate and return large positive
  values (`+361` per variable at one point, which is how the artefact announces
  itself). Check the magnitude before reading a difference.

## What is established, and how it would have to be stated

At the flattened threshold the nested complexity is still positive, `+0.04`: the
nested formulation has clusters where its own CNF flattening has none. That is
the right direction and matches Sec. 13.7's warning-level result. But it is a
one-sided comparison — a number on one side, a sign on the other — and not a
table.

One framing point for whenever this is written up. The flattened system violates
the treelike assumption that both SP and the Bethe counting rest on, *by
construction*: distributing `NOT C` produces `k1` clauses pairwise sharing `k2`
plain literals, which is exactly Sec. 14.1's condition. Its `Sigma = 0` is
therefore the tree approximation's answer to a non-treelike problem. That is the
content worth having, but it has to be stated as *what flattening costs a
treelike calculation*, not as two rival estimates of one truth.

## Running it

```sh
python3 probe/nested_sp.py        # the k1 = 1 relay identity
```

`nested_converge`, `nested_sigma` and their flat counterparts take
`size`/`sweeps`/`nsamp`; the defaults (20 000 / 200 / 300 000) put the noise
floor at a few parts in 10^3, which is comparable to `Sigma` itself near a
crossing. Raise them before trusting any threshold.
