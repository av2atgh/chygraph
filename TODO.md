# TODO — open items

Work that is planned but not done. Items here are *not* covered by the tests;
what is done lives in the README's Results sections.

---

## 1. Test prediction 4 on hyperbolic random graphs

**The claim.** A chygraph whose complexes are the HRG's clustered motifs has a
non-empty leaf-removal core where the degree-matched configuration model does
not — i.e. the formalism sees the effect the `{p_d, e_dd'}` ensemble misses.

**Why it is the item that matters.** It is the whole motivating claim of this
repo (README, *Why now: the clustering gap*), and the one that would kill the
programme. Everything else here — WP1's reweighted tensor, WP2's fixed-point
stability, WP3's hitting set — stands or falls on its own and is about
higher-order constraint satisfaction. Only prediction 4 connects this repo to
the geometry question that `~/av2atg/computational_complexity` opened.

**What is already measured** (`computational_complexity`, no need to redo):

| family | core transition |
|---|---|
| Erdős–Rényi | `k̄_c = e = 2.71828` |
| configuration model, HRG's own `P(k)` | `k̄_c ≈ 4.50` (τ=2.9) … `≈10` (τ=2.3), none at 2.1 |
| HRG | none — core > 0 for all `k̄ > 0`, `core ∝ k̄^β`, β = 1.22–1.63 |

The HRG and its degree-matched configuration model have the same `P(k)` and
essentially the same assortativity (τ=2.3: r = −0.05 vs −0.045). WP3 recovers
`k̄_c = e` analytically, so the ER row is already reproduced here.

**Step 1 is done, and it came out the other way.**
`probe/` measures whether the maximal-clique ensemble converges; see
[`probe/RESULTS.md`](probe/RESULTS.md). It does, for `tau >= 2.5`:

| tau | HRG `sbar` beta | paired ratio HRG/control | ratio beta |
|---:|---:|---:|---:|
| 2.9 | 0.000 – 0.001 | 1.76 – 4.77 | +0.003 – +0.028 |
| 2.5 | 0.009 – 0.024 | 1.79 – 3.52 | +0.009 – +0.015 |
| 2.1 | 0.28 – 0.40 | decays to 1 | −0.12 – −0.18 |

So the blocker this item was written around does not exist in the regime that
matters. The HRG has a well-defined complex ensemble, and it is separated from
its degree-matched control by a converged factor of 1.8–4.8 with `P(k)` and
assortativity held fixed. That separation *is* clustering.

`tau = 2.1` is the exception and fails for the tail, not the geometry: there the
control's cliques diverge *faster* than the HRG's, so the paired ratio decays
toward 1. Anything measured at `tau = 2.1` cannot distinguish the two
mechanisms. Do not use it.

Validation: the measured clique-number exponents (0.37 / 0.20 / 0.07 at
`tau = 2.1 / 2.5 / 2.9`) track the `n^{(3-tau)/2}` of Friedrich & Krohmer across
a 9x range, and the clustering-free baseline (ER, and the control at `tau > 3`)
sits at `sbar = 1.00` as it should.

**What is left to build.**

1. ~~A complex ensemble for the HRG.~~ Done: maximal cliques, `tau >= 2.5`.
   The measured inputs a chygraph needs — clique-size distribution and the
   joint (clique-size, membership) distribution — come straight out of
   `probe/clique_moments.py`; only the joint distribution still needs
   recording, as a `JointChygraph`.
2. *Core percolation in the chygraph map.* Unchanged and now the critical path.
   The leaf-removal core is not the giant component; WP3's
   `sigma = Gbar_k(1 - Gbar_c(sigma))` is the hitting-set map, and the
   core-percolation order parameter is a different fixed point of a related
   anti-monotone map that needs writing down for chygraphs.
   `antimonotone.py` already supplies the solver.
3. *The comparison.* Same `P(k)`, same assortativity, complexes on / complexes
   off, at `tau = 2.5` and `2.9` where both the core effect and the ensemble
   exist. If the core fraction does not separate, prediction 4 is dead.

**How it fails now.** Not at the ensemble. Either at step 2 — core percolation
may not be a chygraph fixed point at all, in which case the formalism cannot
express the order parameter the question is about — or at step 3, where the
chygraph may reproduce the control's core rather than the HRG's, meaning
cliques are the wrong complexes even though they are a valid ensemble.

**Literature.**
- `~/av2atg/computational_complexity/plan.tex`, Observations 1–2, Conjectures.
- Soffer & Vázquez, Phys. Rev. E **71**, 057101 (2005) — clustering corrected
  for degree correlations; the `c̃` that separates core from hubs.
- arXiv:2607.09170 — the HRG family and its Thm-18 induced trees.
- Bauer & Golinelli, Eur. Phys. J. B **24**, 339 (2001) — core percolation at
  `c = e`.

---

## 2. WP4–WP6

Distributional messages, complexes as regions, and the Bethe free energy.
Stated in the README's Work packages; nothing to add here yet.
