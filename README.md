# chygraph_statmech

Statistical mechanics on chygraphs: extending the Bethe–Peierls treatment of
lattice-gas models on correlated random graphs
([Vázquez & Weigt, PRE **67**, 027101 (2003)](https://arxiv.org/abs/cond-mat/0207035))
from graphs to higher-order structures, by way of the chygraph formalism
([Vázquez, PRE **107**, 024316 (2023)](https://arxiv.org/abs/2308.00987);
`~/av2atg/chygraph`).

Status: **framing only.** No code yet. This README is the plan.

## The claim

The chygraph percolation map is already a Bethe–Peierls cavity recursion. It is
not an analogy: it is the same recursion, specialised to the one message type
that percolation needs — a scalar in `[0,1]` instead of a distribution over
effective fields.

VW03 Eq. (4), a message on a directed edge, updated from every incoming message
except the one it came from:

```
h^(i|j) = mu + sum_{k != j} u(h^(k|i))
```

`chygraph` Eqs. (Qm)–(Qp), a message on a directed *inclusion*, with the same
exclusion written as `[Fbar]_{k=m}`:

```
Q^{ml}_- = prod_k Phi^l_k(Q^{lk}_-) prod_k [Gbar^l_k]_{k=m}(Q^{lk}_+)
Q^{ml}_+ = prod_k [Phibar^l_k]_{k=m}(Q^{lk}_-) prod_k G^l_k(Q^{lk}_+)
```

Both are exact on a locally tree-like structure; both assume replica symmetry.

### Dictionary

| VW03 | chygraph |
|---|---|
| edge `(i,j)`, directed | inclusion, directed, plus the `i = ±` up/down role |
| excess degree class `d` | layer index `l` |
| `p(d\|d')` | the `<kappa>_{lk}`, `<s>_{lk}` coupling tables |
| `P_d(h)`, field distribution per class | `Q^{ml}_i`, scalar per class |
| `sum_{d_l} p(d_l\|d) int prod dh_l` — Eq. (10) | `prod_k Phi^l_k(.)` — the pgf *is* that convolution |
| RS instability | `Lambda = max eig(-A)` |
| Eq. (12), site + link free energy | **absent** |

The pgf product in `chygraph.giant.Chygraph.apply` and Eq. (10) of VW03 are the
same operation written two ways. A pgf `Phi(x) = sum_n p_n x^n` evaluated at a
scalar becomes, when the argument is a *measure on fields*, `sum_n p_n P^{*n}` —
the `n`-fold convolution. That single substitution — **argument scalar → argument
measure, product → convolution** — is the whole structural content of the
extension.

## Why now: the clustering gap

`~/av2atg/computational_complexity` measured the empirical side of this and
found the VW03 phase diagram fails where it should not. Hyperbolic random graphs
have an extensive leaf-removal core at **every** mean degree, while their
degree-matched configuration model — same `p_d`, same assortativity — is
core-free up to `kbar ~ 4.5–10`. The separating feature is clustering, which the
ensemble `{p_d, e_dd'}` does not encode and which the Bethe–Peierls derivation
explicitly assumes away.

Chygraphs are the construction that removes that assumption. A triangle, a
clique, any dense motif becomes a *complex*; the chygraph over complexes stays
locally tree-like even when the underlying graph is not. So the gap
Observation 1 of `computational_complexity/plan.tex` identifies is exactly the
gap this formalism is built to close.

That is the motivating claim of this repo, and it is the one most likely to be
wrong — see Risks.

## What transfers for free

**1. Degree–degree correlations are already in the formalism.** `e_dd'` needs no
new machinery, only layer refinement: layers `d = 0..D-1` for node degree
classes plus `D(D+1)/2` layers for edges typed by endpoint-degree pair. That is
literally the index count in `chygraph.percolation.InteractingHypergraphs`
(`L = g + g(g+1)/2`). `e_dd'` becomes `<s>` from an edge-layer down to its two
node-layers. The VW03 correlated ensemble is a chygraph constructible today.

**2. `PercolationMatrix` is the RS-stability tensor with the wrong weights.**
`A`'s entries are first moments because in percolation the message derivative
`du/dh` is identically the occupation probability, already folded into
`<kappa>_{0l} = p<k>` by thinning. For a general Hamiltonian the same `2L^2`
tensor works with

- `<kappa>_{lk} -> <kappa>_{lk} * <u'>` for the ferromagnetic / RS instability,
- `<kappa>_{lk} -> <kappa>_{lk} * <(u')^2>` for the AT / spin-glass line.

`theta()`, `eigenvals()`, `Lambda()` and the whole `chygraph.amplitude`
hierarchy (threshold ← first moments, amplitude ← second moments) then apply
verbatim. This is the cheapest real result available here: **closed-form RSB
thresholds for higher-order structures, from code that already exists.**

**3. The `mu -> inf` hard-core limit collapses back to scalars.** VW03 Eq. (16),

```
pi_d = sum_{d1} p(d1|d) (1 - pi_{d1})^{d1}
```

has no field distribution left in it — a scalar fixed-point map over layer
classes, the shape `chygraph.giant` already solves. It is *not* in the chygraph
class: the map is order-**reversing** where a pgf is order-preserving, so
`Chygraph.solve`'s monotone iteration from `Q = 0` will not converge. Everything
except that sign is shared.

The anti-monotonicity *is* the RSB detector VW03 uses ("an instability prevents
the program from convergence"), and it is the same phenomenon as the
period-doubling route in `chygraph/TODO.md` item 3. Both are loss of
monotonicity in the message map, and `Lambda` at the *non-trivial* fixed point
is the diagnostic for both. `chygraph.giant.Chygraph.jacobian` hardcodes
`subs({q: 1})`; parameterising that to evaluate at the solved fixed point is a
few lines and yields the local-stability criterion directly.

## What has to be built

**1. `G^l` does not survive.** The real obstruction, and the interesting one.
Percolation compresses a complex's entire internal structure into one
component-size pgf (`Gbar_{K_n}`, `Gbar_triangle`) because reachability
factorises. For a general Hamiltonian it does not: a `K_n` complex carries an
`n`-body interaction, and the model must be solved *inside* it and a full cavity
field emitted per member. The compression is percolation-specific.

The caveat under `manuscript_3` Eq. (excessderived) — that `Gbar` must be
supplied rather than derived when the entry vertex matters — is the first
visible symptom of exactly this.

The upside: complexes then become *regions* in the Kikuchi sense, and a chygraph
becomes an **ensemble-level region graph** — a cluster-variation method with a
controlled random ensemble and a layer index. That is the version of this worth
a paper rather than a section.

**2. Replica symmetry only.** VW03's punchline was RSB, and there is nothing
1RSB in the formalism. Survey propagation multiplies the message space again
(distributions of distributions, Parisi `m`). But the *AT condition* — where RS
breaks — is item 2 of the previous section and is cheap. Stop there.

**3. No interaction disorder.** VW03 calls this generalisation "evident" and it
is; it is orthogonal to the structural axis, and the two compose.

## Work packages

**WP1 — Reweighted stability matrix.** Generalise `PercolationMatrix` to carry
`<u'>`- and `<(u')^2>`-weighted moments. Sanity check: `u' = p` must return
`chygraph`'s percolation threshold identically, and `L = 1` must return the
known AT line for the Bethe lattice. *Cheapest path to a new result.*

**WP2 — Fixed-point Jacobian.** Parameterise `Chygraph.jacobian` off `Q = 1`.
Gives local stability at the non-trivial fixed point — the RSB detector — and
serves `chygraph/TODO.md` item 3 at the same time.

**WP3 — Hypergraph vertex cover / hitting set.** The named first target; see
below. Needs an anti-monotone solver (iterate `F∘F`, or damped iteration) beside
`Chygraph.solve`.

**WP4 — Distributional messages.** Population dynamics for `Q^{ml}_i(h)`,
finite temperature, general `w`. Only after WP1–3; the symbolic closed forms are
lost here and the payoff is lower.

**WP5 — Complexes as regions.** The Kikuchi reading. Exact solution inside a
complex, full field vector emitted. This is where clustering-that-VW03-cannot-see
actually gets computed, so it is the package that pays off Observation 1 of
`computational_complexity` — and the one with the most ways to fail.

**WP6 — Bethe free energy.** VW03 Eq. (12) has no chygraph counterpart. On a
chygraph the counting is complex-terms minus inclusion-overcounting, i.e. a
region-graph free energy. Needed for anything variational; not needed for
thresholds.

## First target

Minimal vertex cover / hitting set on a hypergraph with correlated hyperdegree
and cardinality.

- NP-hard, and the natural higher-order lift of the VW03 problem.
- `chygraph.applications.correlated_cardinality_hypergraph` already builds the
  ensemble.
- The `mu -> inf` limit stays scalar, so no population dynamics is needed.
- Generalised leaf removal lifts to hypergraphs cleanly, so the numerical
  control from `computational_complexity/code/leafremoval.py` is reusable.
- Hyperedge constraint: `e^w = 1 - prod_{i in e}(1 - x_i)`.

The physics question is sharp and, as far as I can tell, open: does *cardinality*
heterogeneity make hitting set easy the way degree heterogeneity makes vertex
cover easy — or does the AND-structure of a hyperedge constraint push the
ensemble into RSB earlier?

## Falsifiable predictions

1. WP1's reweighted matrix reduces to `chygraph.percolation.PercolationMatrix`
   identically at `u' = p`. *(If not, the correspondence claim is wrong.)*
2. Layer-refined `e_dd'` (free-transfer item 1) reproduces Fig. 1 of VW03 —
   `x_c(r)` at `gamma = 2.5` and `3.0` — from `chygraph` code with no new solver.
3. The hypergraph hitting-set RSB threshold moves *down* in cardinality
   heterogeneity, opposite to the degree-heterogeneity effect in VW03.
4. A chygraph whose complexes are the HRG's clustered motifs has a non-empty
   core where the degree-matched configuration model does not — i.e. the
   formalism sees the effect the `{p_d, e_dd'}` ensemble misses.

Prediction 4 is the whole motivating claim and the one to attack first, because
it is the one that would kill the programme.

## Risks

- **The HRG may not admit a useful chygraph mapping.** Clustering there is
  geometric and scale-dependent, not a finite catalogue of motifs. Chygraphs
  need a complex *ensemble*; hyperbolic geometry may not give one with finite
  moments. If so, prediction 4 fails and this repo reduces to WP1–3: a genuine
  but narrower contribution about higher-order constraint satisfaction, with no
  bearing on the geometry question.
- **WP5 may not close.** Region-graph free energies are not variational bounds
  in general and can fail to converge. The threshold results (WP1–3) do not
  depend on WP5 landing.
- **Overlap with existing higher-order cavity work.** Factor-graph BP for
  hypergraph constraint satisfaction is standard; the novelty claimed here is
  the *ensemble-level* index structure (layers, `2L^2` tensor, closed-form
  thresholds and amplitudes), not message passing on a given instance. Needs a
  literature pass before WP3 is written up.

## Related

| repo | what it holds |
|---|---|
| `~/av2atg/chygraph` | the percolation formalism, symbolic code, `manuscript_3` |
| `~/av2atg/computational_complexity` | VW03 on hyperbolic random graphs; the clustering gap |
| here | the statistical mechanics that would close it |

## References

- A. Vázquez, M. Weigt, *Computational complexity arising from degree
  correlations in networks*, Phys. Rev. E **67**, 027101 (2003),
  [cond-mat/0207035](https://arxiv.org/abs/cond-mat/0207035).
- A. Vázquez, *Percolation in higher order networks via mapping to chygraphs*,
  J. Complex Netw. (2024), [doi:10.1093/comnet/cnae047](https://doi.org/10.1093/comnet/cnae047),
  [arXiv:2308.00987](https://arxiv.org/abs/2308.00987).
- A. Vázquez, *The giant component of complex hypergraphs*, `chygraph/manuscript_3`.
- M. Weigt, A. K. Hartmann, Phys. Rev. Lett. **84**, 6118 (2000) — hard-sphere
  lattice-gas representation of vertex cover.
- M. Bauer, O. Golinelli, Eur. Phys. J. B **24**, 339 (2001) — leaf removal.
- M. Mézard, G. Parisi, Eur. Phys. J. B **20**, 217 (2001) — cavity method,
  Bethe lattice.
- J. S. Yedidia, W. T. Freeman, Y. Weiss, *Constructing free-energy
  approximations and generalized belief propagation algorithms*, IEEE Trans.
  Inf. Theory **51**, 2282 (2005) — region graphs, Kikuchi.
