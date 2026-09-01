# Set-valued survey propagation for proper colouring

**Status: the machinery is validated at cardinality two and produces nothing at
cardinality three.** The `c = 2` gate is the part worth keeping; the `c = 3`
application does not yield a threshold and is not in the book.

`setvalued_sp.py`. Self-contained (numpy only). Nothing in `book/figs/` imports
it and nothing in the book depends on it.

## What it is for

Sec. 12.7 says the proper colouring rule does not inherit the cardinality-two
one-step machinery: a complex of cardinality `c` can forbid up to `c-1` colours
**at once**, so the survey is a forbidden *set* rather than a scalar and the
prohibitions are correlated inside one complex as well as across complexes. That
is stated in the book as a reason not to attempt it. This file attempts it
anyway for `c = 3`, the smallest case, where the forbidden set has size at
most two.

The target was the colourability threshold of a **triangle-clustered graph**,
which nobody has: Mulet et al. do graphs, Gabrié et al. do hypergraphs, and
clique-structured graphs at 1RSB are neither. Sec. 12.5 gives the RS answer for
the same object; this would have given the 1RSB one.

## The interior, and why it is genuinely harder

At `c = 2`, colour `g` closes to member 0 exactly when its single neighbour is
*forced to* `g`. Scalar rule, one number.

At `c = 3` that is not enough. If `j` and `k` both hold available set `{g, d}`,
then `g` closes to member 0 — dropping `g` leaves both needing `d`, and they
must differ. **Neither j nor k is forced to anything.** The prohibition comes
from a correlation inside the complex that no scalar survey can carry. The
update is a Hall's-condition / SDR check on the other members' sets, and it is
precisely the failure Sec. 12.7 describes.

Colour symmetry still buys a lot: a survey is uniform within each size class, so
it is a `(q+1)`-vector `v[s] = P(|available set| = s)`, and the interior is
precomputed once as a tensor (`interior_tensor`, `_sat_tensor`) so the per-sample
cost is a small contraction.

## What works: the c = 2 gate

At `c = 2` the set-valued code and the scalar code in `book/figs/colouring.py`
compute the same object, so any disagreement is a bug here. Both the two-term
and the three-term complexity forms reproduce Mulet's `c_q`:

| q | two-term | three-term | `c_q` | scalar code (committed) |
|---|---|---|---|---|
| 3 | 4.685 (0.056) | 4.684 (0.026) | 4.69 | 4.683 |
| 4 | 8.930 (0.022) | 8.863 (0.027) | 8.90 | 8.901 |
| 5 | 13.734 (0.060) | 13.636 (0.026) | 13.69 | 13.660 |

All within 0.5%, and the two forms agree with each other — which is what
licenses taking the three-term form to `c = 3`, where the two-term collapse has
no justification.

## What does not work: c = 3

`Sigma` is negative from the moment the branch appears and never crosses zero:

| kappa | degree | Sigma | ⟨forced⟩ |
|---|---|---|---|
| 1.20 | 2.40 | +0.00000 | 0.0000 (branch off) |
| 1.60 | 3.20 | −0.14450 | 0.4979 |
| 2.40 | 4.80 | −0.41240 | 0.8638 |
| 3.00 | 6.00 | −0.64233 | 0.9411 |

Two signs that this is still the implementation rather than the physics:

- **Magnitudes.** `|Sigma| ~ 0.1-0.8` where the validated `c = 2` case near its
  threshold is `~0.01`. Every previous miscount in this line of work announced
  itself the same way.
- **Ordering.** For graphs the branch onset (4.42), the crossing (4.69) and the
  RS stability line (5) sit within 13% of one another. Here the branch appears
  at degree `~3` while Sec. 12.5's RS line for the same object is at degree `6`
  — a factor of two apart, which does not resemble the graph pattern.

The `c = 3` interior is the one component the `c = 2` gate cannot validate, by
construction, and it is where to look first.

**One explanation is ruled out.** All the scans above are at `q = 3`, where the
triangle is *saturated*: three members, three colours, so a proper colouring is
a permutation and the complex is as rigid as it can be. `q = c` could plausibly
be pathological rather than representative. It is not — `q = 4, c = 3` fails
identically, the branch appearing between degree 6 and 8 with `Sigma = -0.38`
already at degree 8 and no crossing below the RS line at 11:

| kappa | degree | Sigma | ⟨forced⟩ |
|---|---|---|---|
| 3.00 | 6.00 | +0.00000 | 0.0000 (branch off) |
| 4.00 | 8.00 | −0.37836 | 0.8022 |
| 5.50 | 11.00 | −0.68343 | 0.9383 |

So the fault is in the `c = 3` treatment itself, not in a degenerate choice of
`q`. Note also what the `c = 2` gate does and does not cover: `site`, `Z_a` and
`Z_ia` are each verified there, so the error is in what changes at `c = 3` --
the SDR satisfiability tensor `_sat_tensor`, the interior tensor, or the
`kappa/c` complex count.

## Five bugs, all caught by gates rather than by reading the code

Recorded because every one of them is invisible in `Sigma` and would recur:

1. **Missing `/q`** in the per-colour forbidding weight (in the hypergraph
   sibling, `book/figs/colouring.py:hyp_converge`). Saturates everything to
   `<e> = 1`.
2. **Collapsing the variable index.** Colour symmetry legitimately averages over
   colours; averaging over variables too is the scalar-closure error one level
   up. Diagnosed, fixed — and the fix changed the answer by 0.001, which
   **disproved the diagnosis** and pointed at (3).
3. **A mean substituted for a random product.** `<g^d> = exp(kappa(g-1))` uses
   only the mean of `g` and destroys the fluctuation the population exists to
   carry. The product must be taken over `d ~ Poisson` *actual* draws. This was
   the real cause of (2)'s symptom.
4. **Bracket floor.** `Sigma` returns exactly `0.0` on the trivial branch, so a
   bisection bracket below the branch onset is refused rather than bracketing.
   Already documented in `NESTED_SP.md`; hit again anyway. Scan for the onset
   before bracketing.
5. **A `c = 2`-only complexity form reused at `c = 3`.** `book/figs/colouring.py`
   uses the two-term `site - (kappa/2) ln Z_a`, which is a collapse valid when
   the edge *is* the complex. The general case needs three terms,
   `site + (kappa/c) ln Z_a - kappa ln Z_ia`, as the hypergraph sibling already
   did. Mirroring the collapsed form gave `Sigma ~ -1` and no crossing.

A sixth near-miss worth recording: the first `c = 2` gate compared a **regular**
chy-degree against Mulet's **Erdős–Rényi** thresholds and the brackets happened
to straddle them. That is not a test, and reading it as one would have validated
a wrong implementation.

## The gate that decides things here

Not the branch onset — that is a fold, and both this code and the scalar one
locate it to only a few per cent, disagreeing with each other by more than
either disagrees with `c_d`. Use `Sigma = 0`, which is a crossing: it is what
separates a real disagreement from a location problem, and it is what settled
the apparent 5% conflict between the two `c = 2` implementations.

## Running it

```sh
python3 probe/setvalued_sp.py      # the c = 2 gate
```

`converge`, `site_and_sigma` and `threshold` take `size`/`sweeps`/`seed`. The
defaults (4000 / 300) are enough for the `c = 2` crossing; the branch onset is
insensitive to population size, so raising it will not sharpen a fold.
