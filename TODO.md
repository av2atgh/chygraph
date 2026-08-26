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

**What has to be built.**

1. *A complex ensemble for the HRG.* The blocker. Chygraphs need a distribution
   over complexes with finite moments; HRG clustering is geometric and
   scale-dependent, not a finite motif catalogue. Candidate reductions, in
   increasing order of honesty:
   - triangles only, with the measured triangle density → almost certainly too
     weak, since `GraphWithTriangles` is already in `chygraph` and a triangle
     ensemble is still locally treelike above the motif scale;
   - complexes = the maximal cliques of the HRG, with their measured size
     distribution and the joint (clique-size, membership) distribution as a
     `JointChygraph`;
   - complexes = the radial shells of the hyperbolic disk, which is where the
     Thm-18 induced trees of arXiv:2607.09170 live.
2. *Core percolation in the chygraph map.* The leaf-removal core is not the
   giant component. WP3's `sigma = Gbar_k(1 - Gbar_c(sigma))` is the hitting-set
   map; the core-percolation order parameter is a different fixed point of a
   related anti-monotone map and needs writing down for chygraphs.
3. *The comparison.* Same `P(k)`, same assortativity, complexes on / complexes
   off. If the core fraction does not separate, prediction 4 is dead.

**How it fails.** Most likely at step 1: no HRG complex ensemble with finite
moments. If so, say so in print — "chygraphs cannot represent hyperbolic
clustering" is a real result about the limits of the formalism, and it closes
Observation 1 of `computational_complexity/plan.tex` negatively rather than
leaving it open. The fallback framing for this repo is then WP1–3 alone:
higher-order constraint satisfaction, no bearing on geometry.

**Cheap first probe before any of the above.** Measure the maximal-clique size
distribution of the HRG at τ = 2.1/2.5/2.9 and `k̄ = 1..8` and check whether its
second moment converges with `n`. One afternoon with
`computational_complexity/code/hrg.py`. If it diverges, step 1 is already
answered and the item can be closed without building anything.

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
