# chygraph_statmech

Statistical mechanics on chygraphs: extending the Bethe–Peierls treatment of
lattice-gas models on correlated random graphs
([Vázquez & Weigt, PRE **67**, 027101 (2003)](https://arxiv.org/abs/cond-mat/0207035))
from graphs to higher-order structures, by way of the chygraph formalism
([Vázquez, PRE **107**, 024316 (2023)](https://arxiv.org/abs/2308.00987);
`~/av2atg/chygraph`).

Status: **the manuscript has been retired in favour of the book**, `book/` —
*Phase transitions on complex hypergraphs*, 202 pages, which carries everything
`main.tex` and `supplement.tex` did and a good deal more. Those two files, the
cover letter and the manuscript figures were deleted on 2026-08-28 and are
recoverable from git at commit `5ebd892`. WP1–WP6 complete; prediction 4 tested
and confirmed; generalised BP on the region graph implemented at the instance
level (`gbp.py`) and exact on the two-triangle example, with the ensemble lift
still open.

**Where to look now.** `book/README.md` is the live working document:
what is drafted, what is open, and which figure script generates which number.
The book's own *The software* chapter maps every computed equation to the
routine that evaluates it and the test that checks it — the successor to the
supplement's Sec. I.

Sections of this file that discuss `main.tex`, `supplement.tex` or the referee
response describe the retired submission. They are kept as a record of how the
results were arrived at; where a number in them disagrees with the book, **the
book is right** — it recomputed three of them and found the manuscript wrong.
See [Where this stands](#where-this-stands--handoff) for what is next.

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

**WP1 — Reweighted stability matrix.** ✅ Done. `stability.StabilityMatrix`
carries `<u'>`- and `<(u')^2>`-weighted moments; `cavity` solves a complex
exactly inside itself to get `u'`. See [Results](#results-wp1).

**WP2 — Fixed-point Jacobian.** ✅ Done. `fixedpoint.FixedPointStability`
evaluates the symbolic Jacobian at the fixed point actually reached, without
modifying `chygraph`. See [Results (WP2)](#results-wp2).

**WP3 — Hypergraph vertex cover / hitting set.** ✅ Done.
`antimonotone.py` solves order-reversing maps by `F∘F` bracketing;
`hittingset.py` is the hypergraph problem. See [Results (WP3)](#results-wp3).

**WP4 — Distributional messages.** ✅ Done. `population.py` carries
`P_l(h)` as a population of samples and solves the full non-linear problem at
finite temperature. See [Results (WP4)](#results-wp4).

**WP5 — Complexes as regions.** ✅ Done. `region.py` builds the region graph
and measures how far a chygraph is from treelike; `gbp.py` passes the
parent-to-child messages on it. The cost of ignoring overlap is measured on
HRGs, and on one instance it is now recovered rather than only measured. See
[Results (WP5)](#results-wp5).

**WP6 — Bethe free energy.** ✅ Done. `freeenergy.py`; it locates the
transition thermodynamically, a third route independent of WP1 and WP4. See
[Results (WP6)](#results-wp6).

## Where this stands — handoff

Manuscript under revision at PRE, **major revision** returned, revision
complete. `COVER_LETTER.md` answers the report point by point. `main.tex` is
**14 pages, the stated limit**, 0 overfull boxes, 0 undefined references, no
stuck floats; `supplement.tex` is 8 pages, same standard.

### What the last change was

Four things, in the order they matter.

**The Supplemental Material now exists** (`supplement.tex`). Sec. VIII of the
manuscript refers to it and previously it did not. It carries three things: the
map from all 37 numbered equations to the routine that evaluates each and the
test file that checks it (referee minor 9); nine derivations written out in full
(minor 14); and the region-graph calculation of `gbp.py` with its limits.

**One of those derivations is new physics.** The manuscript stated that the
simplicial spinodal approaches the Bragg–Williams limit as `1/k` except at
`q = 4`, where it is `1/k²`, and gave residuals with no coefficient. Expanding
the spinodal condition,

    T = T* − C_q/k + O(1/k²),   T* = q(q−1)/2^(q−1),   C_q = q(2q − 2^(q−1))/2^q,

and `C_q` vanishes exactly when `2q = 2^(q−1)`, i.e. at `q = 4` and nowhere
else. At `q = 3` it is `3/4`, reproducing the quoted 1.5e-2, 3.8e-3, 9.4e-4 at
`k = 50, 200, 800`. Pinned by `test_derivations.py`.

**Fig. 4 is new** (referee minor 7): hitting-set density, hard field against
soft field against Mézard–Tarzia, three panels. `probe/hittingset_density.py`
produces the data, `figures.py:hittingset_panel` the figure.

**GBP is implemented** (`gbp.py`), which was the standing item in the
conclusions. Exact on the manuscript's own two-triangle example; see
[Results (WP5)](#passing-the-messages--gbppy) for what it does and does not
settle.

Two defects surfaced while doing this and are fixed: Sec. VI wrote the
chy-degree generating function as `Ḡ` at cardinality two, where `Ḡ` is the
*intra-complex* function defined in Sec. II (it is `Φ̄`); and the rewired-control
paragraph carried a duplicated clause with a sentence restarting in lower case.
Table I of the previous version was removed — the paragraph after it stated
every row in prose, and the page was needed for Fig. 4.

Renumbering, since it changed: tables are now I `tab:tc`, II `tab:soft`,
III `tab:kikuchi`, IV `tab:overlap`, V `tab:pred4`, VI `tab:ensemble`; figures
1 schematic, 2 clustering, 3 simplicial, 4 hitting set, 5 prediction 4.
Equation numbers are unchanged.

### Next, in priority order

1. **Repo is private.** The acknowledgments cite
   `https://github.com/av2atgh/chygraph_statmech`. It must be made public before
   submission or the URL will not resolve for a referee. This is the only thing
   between the manuscript and resubmission.
2. **Referee minor 14 is now answered in the Supplemental Material rather than
   the body.** If the referee wanted the derivations in the main text, the
   page budget does not allow it and the cover letter says so.
3. **Open science, not review:** the *ensemble* lift of generalised BP — a
   message per region type over a region graph that is itself random. The
   instance-level calculation is done and exact where the clique structure is
   chordal; the ensemble average of it is what would close the 17–35% overlap
   deficit on hyperbolic random graphs. This is the standing item in the
   conclusions and it is stated that way there.

### Standing caution

Five claims in this paper were corrected during review, and every one had been
reported as a number without a derivation behind it: the clustering sign, the
`(q-1)` double count, the exponent agreement, `sigma` for the mixed-cardinality
example, and a convergence rate. The referee's closing request — that each
quantitative claim be independently checked — is the right instinct. Tests now
cover the derivations, not just the outputs; `tests/test_derivations.py` exists
for exactly that. Keep it that way.

## One class

`Chygraph` holds the structure — layers, cardinalities, chy-degrees — and every
calculation is a method on it, delegating to the modules below rather than
reimplementing them.

```python
from chygraph_statmech import Chygraph

g = Chygraph([2, 3], [4.0, 2.0])            # links and triangles, Poisson
g.critical_coupling()                        # Ising T_c
g.core().core_fraction()                     # leaf-removal core
g.hitting_set_bp(mu=60).run().density()      # hitting set, O(1) fields
Chygraph([2, 16], [4, 4], regular=True).simplicial([0.7, 0.3]).transition(scan=(0.96, 0.99))
```

| result | method |
|---|---|
| emitted field, any interior Hamiltonian | `emitted_field` |
| branching matrix `B` | `branching_matrix` |
| reweighted `2L²` tensor | `stability_tensor` |
| `u'` of a complex | `u_prime` |
| Ising `T_c` (continuous) | `critical_coupling` |
| de Almeida–Thouless line | `critical_coupling(squared=True)` |
| field distributions | `population().magnetisation()` |
| Bethe free energy | `free_energy().minus_beta_f()` |
| paramagnetic `−βf`, closed form | `paramagnetic_free_energy` |
| simplicial spinodal | `simplicial().spinodal()` |
| metastability limits `T*`, `T**` | `simplicial().coexistence()` |
| **first-order `T_c`** | `simplicial().transition()` |
| hard-field hitting set | `hitting_set()` |
| cover size at any `c` | `hitting_set_bp().run().density()` |
| RS validity | `hitting_set_bp().run().entropy()` |
| induced-graph vertex cover | `clique_cover().cover_size()` |
| core percolation | `core().core_fraction()` |
| measured ensemble | `from_samples(...).core_from_samples()` |
| Möbius counting numbers | `regions().counting` |
| distance from treelike | `overlap_profile` |

Three deliberate choices. `emitted_field` is **exposed**, so a new higher-order
interaction is a function passed in rather than a subclass — that's why the
simplicial model was a substitution. The two hitting-set methods are **separate**
because they answer different questions and agree only at `c = 2`. And
`critical_coupling` **refuses** the simplicial interaction, because there the
Perron condition gives the spinodal and returning it as `T_c` would be 10–50%
wrong.

```python
from chygraph import HypergraphPercolation, MultiplexHypergraph
## Prediction 4, tested

`probe/prediction4.py`, `n = 2×10⁵`, maximal cliques as complexes, ensemble
measured from each graph so nothing is fitted.

`θ` for structural exponents throughout, `β` being reserved for the inverse
temperature (referee minor 2). Uncertainties are bootstrap over seeds.

| τ | θ measured | θ chygraph | θ control |
|---:|---:|---:|---:|
| 2.5 | 1.590 ± 0.061 | **1.553 ± 0.035** | no core |
| 2.9 | 1.635 ± 0.008 | **1.573 ± 0.011** | no core |
| 2.1 | 1.493 ± 0.060 | 1.571 ± 0.013 | no core |
| spread | 0.143 | 0.020 | — |

The chygraph produces a power law `core ∝ k̄^θ` of about the right exponent
while the degree-matched control has **no core at all** across the fit range.
That is the qualitative separation the construction was built for — but two
qualifications, both the referee's and both right.

The predicted exponent **does not move with τ** (spread 0.020, consistent with
none) while the measured one spreads by 0.143 and rises monotonically; at
τ = 2.9 the two differ by roughly four standard deviations. So this establishes
that the chygraph gives a power law of about the right exponent, not that it
tracks the exponent's dependence on the tail. And a *positive* prediction
follows from the algebra the moment one triangle enters the ensemble — any layer
of cardinality ≥ 3 removes the core-free branch — so the qualitative half was
never really at risk. The sharp check is the rewired control below.

Magnitude at τ=2.9: ratio prediction/measurement is 0.97–0.99 at `k̄ ≤ 0.1`,
falls to 0.77 near `k̄ ~ 1.5–3`, recovers to 0.96 at `k̄ = 6` — tracking clique
overlap, negligible when sparse and saturating when nearly everything is core.

**The cleanest control:** at τ=2.9, `k̄=6` the configuration model *does* have a
core, and the chygraph gets it to **2%** (0.633 vs 0.621). Same map, same
pipeline, and its cliques barely overlap (`shared_2plus` 0.00–0.01 vs 0.08–0.32
for the HRG). So the 23% deficit on the HRG at that density is overlap, not the
map.

**τ=2.1 should not be used** — the chygraph *over*-predicts by up to 33% there,
consistent with the ensemble having no `n`-independent limit below τ=2.5.

**A claim of mine that was wrong.** I had written that a chygraph of independent
cliques predicts `1 − Φ(0)`, "which is not a power law in `k̄`". Both halves are
false. `1 − Φ(0)` holds only when *every* layer has cardinality ≥ 3, so no vertex
ever has degree 1 and leaf removal never fires. Real HRG ensembles have many
cardinality-2 cliques, leaf removal does fire, and the fixed point is non-trivial
— and it does give a power law, with nearly the right exponent.

## Absorbing the simplicial Ising model

[Son, Lee & Goh](https://doi.org/10.1038/s42005-026-02724-2) (arXiv:2411.19080)
study an Ising model where a hyperedge lowers the energy only when **all** its
members agree. It drops in without changing anything: `emitted_field` sums
whatever Hamiltonian sits inside a complex, and the chy-degree step never sees
it. Their Bethe–Peierls Eq. (9) *is* the intra-complex step with that energy.

The sum closes, so `simplicial.py` needs no enumeration at any cardinality:

```
Z(S0 = ±1) = (2 cosh h)^{q-1} + (e^a - 1) e^{± h(q-1)},   a = beta*J_q
u'         = (e^a - 1) / (2^{q-1} + e^a - 1)              per neighbour
```

`u'` is `tanh(beta J/2)` at `q = 2` — a simplicial pair is an ordinary bond of
half the coupling, since `delta_{S0 S1} = (1 + S0 S1)/2`.

**The convention, and a bug the referee found.** `u'` above is the derivative
with respect to **one** other member's field, matching `ising.clique_derivative`,
so the multiplicity `q − 1` is supplied once by the branching matrix and not
here. Differentiating instead with respect to a field *common* to all `q − 1`
others multiplies it by `q − 1`; that combination is what enters a symmetric
fixed point directly, and this file used to quote it while also feeding it to
Eq. (8), which counts the multiplicity twice. `test_simplicial.py` pins the
relation between the two conventions so the trap cannot reopen. No numerical
result changed — the code was consistent; the formula written here was not.

**The caveat that matters.** `det(I − B) = 0` locates a *linear* instability, so
where the transition is discontinuous it is the **spinodal**, not `T_c`. The
ordered branch has to be found by iterating from a magnetised start and the true
temperature by comparing free energies.

**Reproduced (their Fig. 4, `q=16`, `(J_q,J_2,k_q,k_2)=(0.3,0.7,4,4)`):**
continuous transition at `T = 1.0108`, equal to the pairwise `3 tanh(J_2/2T)=1`
to three digits because `u'` is `O(2^-q)`; `m_2` exceeds `m_q` by 2–3 orders of
magnitude between the transitions; and a genuinely first-order second
transition, located by free energy rather than hysteresis:

| | |
|---|---|
| continuous transition | `T = 1.01081` |
| coexistence window | `T* = 0.97701` … `T** = 0.98095` |
| **first-order transition** | **`T_c = 0.97898`** (free energies cross) |
| jump | `m: 0.548 → 0.834`, `Δm = 0.286` |

with `T* < T_c < T**` as required. Below it the two shares become comparable
(ratio ~3) and both jump. Every claim of theirs, recovered.

q-uniform at `k=4`: `T_c = 0.938, 0.714, 0.480` for `q = 6, 8, 12`, each inside
its metastability limits with `Δm > 0.9`; `q=3` has no coexistence window at all.

**Two traps.** The spinodal is the *lower* limit of metastability, not `T_c` —
quoting it would be wrong by 10–50% here. And just below a *continuous*
transition the near-zero branch converges slowly enough to mimic coexistence;
that's critical slowing down, and `branch_gap` documents it.

**Extended, and the normalisation has to be stated.** At *fixed* `J` on a
`q`-uniform Bethe hyperlattice with `k=4`, `3(q−1)u'=1` gives `e^a = 2` at both
`q=2` and `q=4` (so `T = 1/ln2 = 1.4427`) and `e^a = 9/5` at `q=3`
(`T = 1.7013`): the spinodal is **non-monotonic**, peaking at `q=3`. That is
their "ambivalent effect of group size" in closed form — but it is *not*
comparable with their Eq. (8), which fixes `rho_q J_q = 1` instead, i.e.
`J_q = q/k`. Imposing that normalisation and letting `k` grow,

```
T* -> q(q-1) / 2^(q-1)                        (their Eq. 8)
T   = T* - C_q/k + O(1/k^2),  C_q = q(2q - 2^(q-1)) / 2^q
```

so their maximum of `3/2` shared between `q = 3` and `q = 4` is recovered, and
at finite chy-degree the degeneracy is lifted. `C_q` vanishes exactly when
`2q = 2^(q-1)` — **at `q = 4` and nowhere else** — which is why the limit is
approached as `1/k^2` there and `1/k` everywhere else, and why `q = 4` is the
cardinality that is tricritical in the mean-field theory. At `q = 3`,
`C_3 = 3/4`, giving residuals `1.5e-2, 3.8e-3, 9.4e-4` at `k = 50, 200, 800`.
Derived in `supplement.tex` S4 B, checked in `test_derivations.py`.

**On the tricritical cardinality.** At `k = 3, 4` the transition is continuous
through `q=5`; from `k = 6` onward it is continuous only through `q=4`, the
Bragg–Williams value, and stays there to `k = 400`. The two treatments answer
different questions — Bragg–Williams is exact at infinite connectivity, Bethe on
a tree — and they agree in the limit where the first is exact. The shift is a
finite-connectivity effect confined to very sparse hyperlattices, not a
correction to their number.

## Can it solve the Ising model, and vertex cover?

**Ising, yes, generally.** `ising.py` gives `T_c` for any chygraph from
`det(I − B) = 0` with the `L × L` branching matrix

```
B_{lm} = [<kbar>_m if m == l else <kappa>_m] (c_m - 1) u'_m
```

The linear algebra is `L × L` however large the complexes are; `u'_m` is exact
by enumeration for any cardinality that can be enumerated. It reproduces every
closed form in Results (WP1) to `1e-11` and extends to cardinalities WP1 never
handled. Passing `squared=True` returns the AT line instead.

**Vertex cover: on graphs and hypergraphs yes; on chygraphs of larger cliques
the formalism computes a value and proves it untrustworthy.**

| problem | status |
|---|---|
| graph VC, correlated `e_dd'` | `vertexcover.py` — VW03 Fig. 1, Weigt–Hartmann to 10 digits |
| hypergraph hitting set, thresholds | `hittingset.py` — `k(c−1) = e`, 12 digits |
| hypergraph hitting set, density | `softfield.py` — MT Eq. (11), `rho = 1/K`, Weigt–Hartmann |
| VC of a chygraph's induced graph | `cover.py` — exact at `c = 2`, **never certified** at `c ≥ 3` |

`cover.py` is the other end of `hittingset.py`'s family: a hyperedge needs one
member taken, so `c−1` may be left out; a clique needs all but one, so **one**
may. The cavity differs in one place, and it is inside the complex — taking `i`
blocks the whole complex, but the complex could only ever have contributed one
member, so the cost is the **max** over it, not the sum:

```
sigma_m = Phibar^(m)(tau),   tau_l = (1 - sigma_l)^{c_l - 1}
```

against `tau_l = sigma_l^{c_l-1}` for hitting set. Both collapse at `c = 2`.

The `c ≥ 3` limitation is not an extra assumption — it is the core result: a
cardinality-≥3 layer has no core-free branch, so leaf removal leaves an
extensive core at every density and never proves a cover minimal. The symptom is
visible: on isolated triangles every vertex sits at `z = 0`, the graph degeneracy
rule returns `1/2`, and the truth is `2/3`. `certified()` reports this, and a
test asserts it equals `core.has_core_free_branch()`.

## First target — done, see [Results (WP3)](#results-wp3)

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

## Results (WP1)

`StabilityMatrix` is `PercolationMatrix` with the four moment tables multiplied
elementwise: `wkappa` on the inclusion channel, `ws` on the intra-complex one.
The `2L^2` index structure, the block layout, `theta = -det(A)` and
`Lambda = max eig(-A)` are untouched. Weight placement follows the two kinds of
step: going *up* a chy-degree a complex relays its own field (weight 1); going
*down* into a complex the interaction transforms it (weight `u'`).

`cavity.py` computes `u'` by exact enumeration inside the complex —
`u = ½ ln[Z(σ₀=+1)/Z(σ₀=−1)]`, differentiated at the trivial fixed point. That
is the step percolation does not need, and WP5 at the smallest scale.

### Recovered

| model | `theta = 0` | status |
|---|---|---|
| bond percolation, graph | `<kbar> p = 1` | Molloy–Reed; `A` identical to `chygraph` entry-for-entry |
| Ising, graph | `<kbar> tanh(beta J) = 1` | textbook |
| Ising AT line, graph | `<kbar> tanh²(beta J) = 1` | textbook |

### New

A triangle solved exactly inside itself transmits, per neighbour,

```
u_T = t / (1 - t + t²),      t = tanh(beta J)
```

which **exceeds `t` for all `0 < t < 1`** and reduces to `t` at leading order.
So a triangle is two independent edges at high temperature and strictly better
than two independent edges below that. For a graph with links and triangles the
threshold is the two-layer branching determinant with `u_T` in place of `t`:

```
K_L t + 2 K_T u_T + 2 t u_T (k_L k_T − K_L K_T) = 1
```

the third term vanishing for Poisson layers.

**The sign, and a claim of mine that was wrong.** This section used to say that
moving neighbours from links into triangles at fixed mean degree *raises* `T_c`,
by +20.6% at four neighbours, and that the two configurations are identical in
`{p_d, e_dd'}`. Both halves are false, and the referee caught it. Matching the
*mean* degree is not matching the degree distribution: `k_T` Poisson triangles
give `d = 2X` with `X ~ Poisson(k_T)`, so the degree is even, its variance is
twice a Poisson's, and the excess degree the branching matrix uses is `n + 1`
rather than `n`. The +20.6% was a degree effect wearing a clustering label.

Against a null that really is matched, clustering **lowers** `T_c`
(`examples/clustering_lowers_tc.py`, Table I of the manuscript):

| n | `t_c` links | `t_c` tri | `T_c` links | `T_c` tri | ΔT_c |
|---:|---:|---:|---:|---:|---:|
| 4 | 0.33333 | 0.38197 | 2.8854 | 2.4853 | **−13.9%** |
| 6 | 0.20000 | 0.20871 | 4.9326 | 4.7209 | −4.3% |
| 8 | 0.14286 | 0.14590 | 6.9521 | 6.8052 | −2.1% |
| 10 | 0.11111 | 0.11252 | 8.9628 | 8.8498 | −1.3% |
| 20 | 0.05263 | 0.05278 | 18.9824 | 18.9296 | −0.3% |

an `n`-regular graph against `n/2` triangles per vertex — identical
`p_d = δ_{d,n}`, identically neutral `e_dd'`, same edges per vertex. A Poisson
comparison against a configuration model carrying the triangle ensemble's own
excess degree `<kbar> = n + 1` gives the same sign.

**Why, mechanically.** A triangle transmits more *per traversal*, `u_T > t`. But
arriving at a degree-`n` vertex through one of its triangles leaves `n − 2`
branches rather than `n − 1`: one of the two neighbours in the triangle just
traversed is already accounted for, and its influence sits inside `u_T` rather
than being counted again. At `n = 4` the enhancement would have to carry `u_T`
from `1/3` to `1/2` to break even and reaches only `3/7`. **The lost branch
wins.** The effect is `O(t²)`, hence the `1/<k>` decay.

**Confirmed outside the formalism.** Wolff cluster updates on a 4-regular random
graph and on a network of two triangles per vertex, every vertex at degree
exactly 4 in both, put the Binder crossings at `2.894` and `2.482` against
`2.8854` and `2.4853` from the branching matrix — a shift of −14.2% against a
predicted −13.9%. `probe/ising_mc.py`; no cavity equation enters it.

### What this does and does not establish

It confirms the mechanical claim — the reweight is a pure elementwise multiply,
no structural change — and shows the machinery returns known answers where
known answers exist. Prediction 1 is true *by construction* once the reweight is
written that way; what carries content is that the **same** multiply, with
weights derived from a Hamiltonian rather than fitted, reproduces the Ising and
AT lines and then extends to a clustered structure that the VW03 ensemble cannot
represent at all.

It does **not** yet touch the HRG question. The triangle result shows the
formalism *can* see clustering the `{p_d, e_dd'}` ensemble cannot; whether it
can see the *specific* clustering that cores hyperbolic random graphs is
prediction 4, still open. The measurement in [`probe/`](probe/RESULTS.md) has
since removed the obstacle that looked most likely to stop it.

## Results (WP2)

`FixedPointStability` takes any `chygraph.Chygraph` and evaluates
`Matrix(F).jacobian(Q)` at the fixed point `solve` reaches, rather than at
`Q = 1`. `chygraph` is used as given — no upstream edit.

### Exchange of stability, made quantitative

| | Perron root / spectral radius |
|---|---|
| trivial fixed point `Q = 1` | `1 + Lambda` — the WP1 threshold diagnostic |
| physical fixed point `Q*` | `rho = 1 − Lambda + O(Lambda²)` |

Verified to `0.15%` at `Lambda = 1.1e-3`. The two roots move apart from 1 at
the *same rate*, in opposite directions. `rho < 1` everywhere in the percolating
phase and decreases monotonically with occupation.

One wart found and fixed: the sign of the leading eigenvalue does **not**
diagnose monotonicity. A chygraph's coupled core alternates between up and down
message roles, so `J` is a non-negative matrix of period 2 and Perron–Frobenius
puts eigenvalues at both `+rho` and `−rho` (verified exactly). `monotonicity()`
reads the entry signs instead, which is unambiguous: `J >= 0` for generating
functions, `J <= 0` for the vertex-cover map.

### The anti-monotone branch: VW03 Fig. 1

`vertexcover.py` implements VW03 Eqs. (16)–(18). Under the correlation model
`e_dd' = q_d[r delta + (1−r)q_d']` the map closes on a single scalar, so
bisecting on it reaches the fixed point **whether or not it is stable** — and
stability is then asked separately of the Jacobian there. That separation is
what WP2 buys. VW03 detects RSB *by* the iteration failing to converge; here the
fixed point is found and the diagnosis computed, two instruments instead of one.

RS breaking point, `p_d ~ d^-gamma`:

| `d_max` | `gamma = 2.5` | `gamma = 3.0` |
|---:|---:|---:|
| 200 | 0.7082 | 0.7805 |
| 400 | 0.7059 | 0.7803 |
| 800 | 0.7042 | 0.7802 |
| 1600 | 0.7031 | 0.7802 |

`gamma = 3.0` is converged; `gamma = 2.5` still drifts at `d_max = 1600`
because `<d>` itself converges slowly there — call it `0.70` and no more digits.
The same caveat applies to `x_c` at `gamma = 2.5`, which is still moving in the
third decimal at `d_max = 6400` (`~0.345`, drifting down by a factor `0.7` per
doubling); at `gamma = 3.0` it is converged to `0.40574`. Two decimals are real
at `gamma = 2.5`, four or five at `gamma = 3.0`. None of the *structural*
claims — curve ordering, monotonicity in `r`, RS at `r = 0` — depend on this.
The heavier tail breaks RS *earlier*, and both sit inside `(0,1)`, matching
"the RS solution breaks at a certain value of `r` that depends on `gamma`".

Cover size `x_c(r)` reproduces the structure of Fig. 1: rising with `r`,
`gamma = 2.5` the lower curve and `gamma = 3.0` the upper, uncorrelated
scale-free graphs replica-symmetric. `examples/vw03_figure1.py`.

### What validates it

The published figure has not been digitised, so the agreement above is
structural, not numerical. The numerical validation is independent and exact:
for uncorrelated Poisson graphs Eq. (17) must reduce to the Weigt–Hartmann
closed form `x_c = 1 − (2W + W²)/(2c)` with `W = LambertW(c)`. It does, to
**ten digits**, at `c = 1, 2, 3, 5, 10`. The rank-one secular criterion for the
instability is separately cross-checked against a dense eigensolve, which
reproduces every `r_RSB` digit above.

## Results (WP3)

`hittingset.py` derives and solves the `mu -> inf` cavity recursion for minimum
hitting set. With hyperedges split into layers by cardinality so a node's
participation can be correlated across them,

```
sigma_m = Phibar^(m)( 1 - sigma_1^{c_1-1}, ..., 1 - sigma_L^{c_L-1} )
x_c     = 1 - Phi(1-tau) - (1/2) sum_m <k_m> tau_m sigma_m,   tau_l = sigma_l^{c_l-1}
```

`Phi` being the joint hyperdegree generating function and `Phibar^(m)` its
inclusion-biased excess — the `JointChygraph` construction, with cardinality
layering from `correlated_cardinality_hypergraph`.

`antimonotone.py` supplies the solver the map needs. For order-reversing `F`,
`G = F∘F` is order-preserving; iterating `G` from `0` and from `1` gives its
least and greatest fixed points `a <= b`, with `F(a) = b`. If `a == b` there is
one stable fixed point and RS holds; if `a < b` that is a period-2 orbit and the
fixed point of `F` lies strictly between, unstable. **The bracket both finds the
solution and diagnoses it** — VW03's "the program fails to converge" becomes two
numbers saying by how much. It agrees with the Jacobian criterion at every point
tested.

### Recovered

| | | |
|---|---|---|
| `c = 2`, Poisson | `x_c = 1 − (2W + W²)/(2k)` | Weigt–Hartmann, 11 digits |
| `c = 2`, RSB point | `k = e = 2.718281828` | Bauer–Golinelli core percolation |

`c = e` is the Erdős–Rényi control in `~/av2atg/computational_complexity`, so
that repo's ER row is now reproduced analytically here.

### New

**Cardinality alone does nothing.** For fixed cardinality `c` and Poisson
hyperdegree the RSB point is exactly

```
k_RSB = e / (c − 1),      i.e.   k (c − 1) = e   for every c
```

verified to 12 digits at `c = 2..20`. Counted in *neighbours*, hyperedge size
does not move the transition at all.

**Spread in cardinality postpones it.** At fixed mean excess cardinality `2`:

| cardinality mix | `k <cbar>` at RSB |
|---|---:|
| all `c = 3` | 2.71828 |
| half `c = 2`, half `c = 4` | 4.84 |
| ¾ `c = 2`, ¼ `c = 6` | 6.83 |
| 9/10 `c = 2`, 1/10 `c = 12` | 6.04 |

So the README's open question — *does cardinality heterogeneity make hitting set
easy the way degree heterogeneity makes vertex cover easy, or does the
AND-structure of a hyperedge push into RSB earlier?* — answers **the first way.**
Heterogeneity keeps the problem easy longer, the same direction VW03 found for
degree, and the AND-structure does not offset it.

**Correlation is the axis that makes it harder.** Anti-correlating a node's
participation in small and large hyperedges, at identical marginals, brings RSB
forward monotonically: `6.83 -> 6.50 -> 5.20 -> 3.76 -> 2.64` as the spread
grows. That axis has no counterpart in VW03, whose ensemble has one edge type.


### Checked against Mézard–Tarzia (2007), and two things were wrong

Mézard & Tarzia, [PRE **76**, 041124](https://doi.org/10.1103/PhysRevE.76.041124),
solve hitting set on random regular hypergraphs with the full cavity method.
Comparing exposed two limits of the `mu -> inf` **hard-field** ansatz here
(warning propagation — the VW03 limit, which drops the `O(1)` part of the fields):

**1. The cover size is wrong at `c ≥ 3`.** The `1/2` weight VW03 gives a
degenerate `z = 0` vertex is right for a two-fold degeneracy only. Disjoint
3-hyperedges, one per vertex — *no interaction between complexes at all*, so RS
is trivially valid — need one vertex of every three taken, so the truth is `1/3`.
The rule returns `1/2`. Against MT Fig. 5, a regular hypergraph with 4 hyperedges
per vertex and 6 per hyperedge has `rho_cov ≃ 0.178`; this gives `0.252`.
`cover_bracket()` now brackets it honestly and `certified()` is False at `c ≥ 3`.

**2. The stability point is not the phase boundary at `c ≥ 3`.** `k(c−1) = e` is
a property of the hard-field map. At `c = 2` that coincides with the accepted RSB
point, which is why it returns `e`. On regular hypergraphs it reports broken
symmetry across almost the whole `(L, K)` grid where MT Fig. 4 finds a genuine RS
region. It *does* reproduce their result that `K = 2` breaks for every `L`.

Read `k(c−1) = e` as **where leaf-removal-style certification stops**, not where
the phase boundary lies. The σ recursion, the order-reversing structure, the
`F∘F` bracket, and every heterogeneity/correlation *threshold* result are
unaffected — they depend on the map, not on the degeneracy rule.

### Fixed: `softfield.py` keeps the O(1) fields

Both symptoms above come from scaling `h = mu*z` and keeping only the integer
`z`. Not scaling at all repairs them. Belief propagation on the chygraph,

```
v_{a->i} = -ln[ 1 - prod_{j in a\i} 1/(1 + e^{h_{j->a}}) ]
h_{i->a} = -mu + sum_{b ni i, b != a} v_{b->i}
```

with `rho = sigmoid(-mu + sum_b v)`. For a regular hypergraph symmetry closes it:

```
h_RS = -mu/L - ((L-1)/L) ln(K-1),    rho = 1/K
```

which is **MT Eq. (11)**, reproduced to 4 decimals. Note `h` carries a *fraction*
of `mu`, not an integer multiple — exactly what the hard ansatz cannot represent.

| cardinalities | `<k>` | hard | soft | exact |
|---|---|---:|---:|---|
| 3 (regular, L=1) | 1 | 0.500 | **0.333** | 1/3 |
| 2 | 1.0 | 0.2720 | 0.2721 | 0.27203 (Weigt–Hartmann) |
| 6 (regular, L=4) | 4 | 0.252 | 0.1667 | 0.178 (MT 1RSB) |
| 3 | 1.0 | 0.208 | 0.162 | — |
| 4 | 1.0 | 0.172 | 0.110 | — |
| {2,3} | {1.0,0.5} | 0.310 | 0.300 | — |

The correction grows with the weight on `c ≥ 3` — 56% for `c=4` alone, under 1%
once ordinary edges dominate. The last rows are outside MT's regular ansatz.

**Its own validity criterion**, replacing the hard-field instability above
`c = 2`. The Bethe free energy applies with `Z_i = 1 + e^H` and
`Z_a = prod_i (1 + e^{h_i}) - 1`, giving

```
s = sum_l n_l <ln Z_a> + <(1-k) ln Z_i> + mu*rho
```

which reduces to MT Eq. (13) on regular ensembles and reproduces it **to nine
digits**. At (4,6) `s = -0.157` and RS gives 0.1667 below MT's 0.178, as an
underestimate should; the hard-field 0.252 is above and further. At (6,12)
`s = -0.192` and the iteration stops converging — the same failure, visible.

Two honest qualifications. The estimator cancels two `O(mu)` terms to leave an
`O(1)` answer, so heterogeneous ensembles need `entropy_averaged()`; the regular
case is exact only because symmetry removes the sampling. And **the criterion is
sufficient, not necessary** — negative `s` proves RS wrong, positive `s` does not
prove it right. ER vertex cover breaks RS at mean degree `e`, yet `s = +0.093(8)`
at mean degree 1 and is still positive at 3. MT pair the entropy with a separate
stability criterion for this reason; only the entropy is implemented here.

**A bug worth recording.** Damping population dynamics by *averaging field
values* (`h ← (1-λ)h + λh_new`) looks right and is not: the two entries are
different messages, not two estimates of one, so averaging contracts the
distribution and *moves the fixed point*. It returned 0.172 for ER at `k=1`
where the truth is 0.272 (confirmed by exact leaf removal, empty core, `dx=0`).
Replacing a random *subset* of the population fixes it. `population.py` does
full parallel updates and was never affected.

### A confound worth naming

The obvious "positive correlation" construction — half the nodes at double the
hyperdegree, half at zero — gives an RSB point *exactly half* the independent
one. That is not correlation, it is dilution: the zero class is isolated
vertices, so the ensemble is the independent one at half density. Matched
constructions keep both classes populated, and `isolated_fraction()` reports
`Phi(0,...,0)` so the confound is visible rather than silent. It still grows at
extreme spread (0.47 at `s = 0.95` on the positive branch), which is why only
the negative branch is quoted as a clean measurement above.

## Results (WP4)

`population.py` replaces the scalar message with the distribution itself:
`P_l(h)` is the law of the field a node sends up into a layer-`l` complex,
carried as a population of samples. The up step is a Poisson-compound
convolution where the chygraph map has a generating-function product; the down
step is the exact enumeration inside the complex, evaluated numerically. That
is the README's substitution — *argument scalar → argument measure, product →
convolution* — actually implemented.

### It re-derives WP1 without sharing WP1's route

WP1 gets `T_c` by symbolic linearisation and never represents a field
distribution. WP4 solves the full non-linear stochastic problem and reads `T_c`
off the order parameter. At matched neighbour count (6 either way):

| system | WP1 closed form | WP4 population | rel err |
|---|---:|---:|---:|
| graph, `k_L = 6` | 0.168236 | 0.164628 | 2.14% |
| triangles, `k_T = 3` | 0.146947 | 0.144343 | 1.77% |

| | `T_c` gain from clustering |
|---|---:|
| WP1 closed form | 14.49% |
| WP4 population | 14.05% |

Both absolute values sit ~2% low for the same reason — with finite sweeps the
bisection calls a slowly decaying magnetisation "ordered" and stops just below
the true threshold. The bias is common to both systems, so the **ratio**, which
is the physical claim, is an order of magnitude more accurate than either
endpoint. `examples/wp4_validates_wp1.py`.

### The sharp check

The expensive comparison above is a few percent. The cheap one is exact: the
derivative of WP4's numeric `emitted` at zero field must equal WP1's symbolic
`cavity_derivative` for the same complex. It does to `1e-8` at `c = 2` and
`c = 3` across `beta J = 0.2 .. 2.0` — two independent enumerations of the same
object, pinned against each other. And `c = 2` reproduces
`atanh(tanh(beta J) tanh(h))` to `1e-12` at *every* field, not only at zero.

### What it cost

The symbolic closed forms are gone, as the work package predicted: `theta`,
`Lambda` and the amplitude hierarchy have no counterpart here, and a single
`T_c` costs ~60s of bisection against a `sympy` one-liner. The payoff is that
nothing in it is a linearisation, so it can be run away from the threshold and
at any temperature — which is what WP6 will need.

## Results (WP5)

Everything above treats a node's complexes as independent. That is the Bethe
approximation *on the chygraph*, exact when complexes meet in at most **one**
node. Two cliques sharing an edge are a loop the chygraph cannot see, just as a
triangle is a loop an ordinary graph cannot see. **Chygraphs move the problem up
one level; they do not remove it.**

`region.py` closes the complexes under intersection and counts by Möbius
inversion, `c_R = 1 − sum_{R' ⊃ R} c_{R'}`, so every node is counted once. When
complexes meet in ≤ 1 node the family is the complexes plus the single nodes and
`c_v = 1 − k_v` — exactly the Bethe counting WP1 and WP4 assume, so the Kikuchi
construction *contains* them. Two triangles on a shared edge separate: Kikuchi
puts `−1` on the **edge**, Bethe puts `−1` on each endpoint. Both count nodes
correctly; only Kikuchi subtracts the correlation.

### What ignoring overlap costs, on the HRG

Maximal cliques as complexes, ensemble measured from the graph itself,
`n = 3×10⁴`. "rewired" is that same clique ensemble arranged treelike:

| τ | k̄ | chygraph | rewired | real HRG | overlap gap |
|---:|---:|---:|---:|---:|---:|
| 2.5 | 2.00 | 0.0786 | 0.0780 | 0.1203 | 0.042 |
| 2.5 | 4.02 | 0.2040 | 0.1970 | 0.2993 | 0.102 |
| 2.9 | 2.01 | 0.1688 | 0.1742 | 0.2309 | 0.057 |
| 2.9 | 4.00 | 0.4298 | 0.4329 | 0.5167 | 0.084 |

**chygraph = rewired to 0.001–0.007**, so `core.py` is right on its own terms —
now checked on a measured, multi-layer, non-Poisson ensemble. The remaining
`real − rewired` is *purely* clique overlap, with the ensemble held fixed. The
chygraph accounts for **65–83%** of the HRG core; the rest is what region
counting would have to supply. The degree-matched control sits at 0.0000
throughout, chygraph and measured alike.

So the qualitative separation prediction 4 is about is reproduced. Its
**amount** is not, and the residual is now attributed to a named, measured
mechanism rather than left as a gap. `probe/overlap_cost.py`.

### A claim of mine that was too strong

I wrote earlier that a complex of cardinality ≥ 3 "is a core by itself". That is
wrong. Leaf removal deletes a leaf's *neighbour*, and that neighbour can be a
complex member — a triangle with one pendant edge on each vertex has an **empty**
core, every vertex removed (now a test). What the algebra actually gives is
narrower: a chygraph with any layer of cardinality ≥ 3 has no core-free *branch*,
so the core is strictly positive at every density — but it can be `1.8e-4`. The
strong reading holds only when *every* layer has cardinality ≥ 3, where no vertex
ever has degree 1, leaf removal never fires, and the core is exactly `1 − Φ(0)`.

### Passing the messages — `gbp.py`

`region.py` builds the object and measures the discrepancy. `gbp.py` passes the
messages: generalised belief propagation in the parent-to-child form of
[Yedidia–Freeman–Weiss (2005)](https://doi.org/10.1109/TIT.2005.850085).

Three things had to hold before any number it produces was worth quoting, and
`tests/test_gbp.py` checks all three.

**It contains BP.** On the two-layer region graph of a pairwise model — one
region per interaction plus one per variable — the denominator `D(P,R)` holds
only the edge being updated and the update *is* belief propagation. `gbp.GBP`
reproduces an independently written BP to `1e-10` on a chain, a star, a
four-cycle and a theta graph: the loopy error included, so the two agree on the
mistake as well as the answer.

**It is exact on a junction tree.** Two triangles sharing an edge — the
manuscript's own Table III — closes on `{0,1,2}, {1,2,3}, {1,2}` once the
zero-counting singletons are pruned, and that is a junction tree.

| βJ | exact ln Z | Bethe | Kikuchi (static) | GBP |
|---:|---:|---:|---:|---:|
| 0.2 | 2.8887 | +1.8e-2 | −1.4e-3 | <1e-13 |
| 0.5 | 3.5907 | +9.1e-2 | −2.9e-2 | <1e-13 |
| 1.0 | 5.7390 | +3.7e-1 | −6.6e-2 | <1e-13 |
| 2.0 | 10.6938 | +1.3 | −1.7e-2 | <1e-13 |

So the whole of what the Kikuchi column reports as left on the table is
recovered, and the residual overshoot of that column is an artefact of
evaluating a counting on *isolated* regions rather than of the counting itself.

**Its fixed point is a fixed point.** Where the region graph is not a junction
tree GBP is approximate, but the converged solution still satisfies
`sum_{x_P \ x_R} b_P = b_R`, so the residual error is the approximation's and
not the solver's. `GBP.consistency()` reports it.

**On real clique structures** (`probe/gbp_cliques.py`): 60 maximal-clique region
graphs of HRGs small enough to enumerate exactly, `n = 14, 18, 20`.

| | runs | GBP error | Kikuchi (static) | Bethe (static) |
|---|---:|---|---|---|
| chordal | 40 | ≤ 3e-12 | — | — |
| non-chordal, converged | 9 | 1.4e-9 – 5.4e-3 | 0.48 – 1.1 | 0.32 – 22 |
| non-chordal, did not converge | 11 | — | — | — |

Exact on every chordal instance, since there the maximal cliques and their
intersections *are* a junction tree. On the non-chordal ones it converges less
than half the time even at damping 0.999; where it does, it beats either static
counting by two to eight orders of magnitude and its marginals agree to `9e-8`;
where it does not, the residual says so and no number from it is worth quoting.
**That instability, not the accuracy, is what an ensemble treatment inherits** —
which is the useful thing this measurement says about the open problem.

**Two limits, stated rather than glossed.** It is not variationally bounded: on
a region graph built from complexes that *contain* one another — K4 covered by
its four triangles — it errs by `1e-2` where the static counting errs by
`4e-3`. That case does not arise for maximal cliques, which never nest. And it
is an *instance-level* calculation: every other module here works at the level
of an ensemble, where a message is a distribution over a chy-degree class rather
than a function on one region. Lifting the parent-to-child update to a message
per region **type**, over a region graph that is itself random, is what would
close the 17–35% gap on hyperbolic random graphs, and it is not done.

## Results (WP6)

VW03 Eq. (12) is a site-plus-link sum with weight `d − 1` per vertex and has no
chygraph counterpart. Parameterising the messages as `exp(h σ)` and `exp(u σ)`,
the overlap term of the general Bethe free energy collapses — because
`h_{i→a} + u_{a→i}` is the *same* full field `h_i` for every complex containing
`i` — leaving one term per complex and `1 − k` per node:

```
-beta f = sum_l n_l <ln Z_a^(l)> + <(1-k) ln Z_i>,     n_l = <k_l> / c_l
```

with `Z_a` the exact partition function *inside* a complex given its members'
cavity fields. **Those weights are WP5's Möbius counting numbers** in the
treelike case — `1` per complex, `1 − k_v` per node. WP5 computes the counting;
WP6 evaluates the free energy with it, and where complexes overlap both are
wrong in the same way. That identity is a test.

### Closed form and estimator

At zero cavity field, `-beta f = ln 2 + sum_l (<k_l>/c_l) ln(Z_c/2^{c_l})`,
which for a graph is the textbook `ln 2 + (c/2) ln cosh(beta J)`. The population
estimator reproduces it **to machine precision** — but only after a fix worth
recording: split `ln Z_i = ln 2 + ln cosh h` and take the `ln 2` piece exactly
rather than sampling it. Sampling costs a bias `~sqrt(<k>/n) ln 2 ≈ 5e-3` at
`<k> = 6`, `n = 8e4`, which is *larger than the free-energy differences the
module exists to resolve*. I hit that bias before spotting it: the error was
constant in `beta J`, which is what gave it away.

### The transition, located thermodynamically

Gap between the ordered and paramagnetic branches of `-beta f`:

| system | 0.8 βJc | 0.95 βJc | 1.05 βJc | 1.2 βJc | 1.5 βJc |
|---|---:|---:|---:|---:|---:|
| graph, `k = 6` | 1e-15 | -4e-14 | 5.1e-4 | 1.5e-2 | 7.8e-2 |
| triangles, `k = 3` | -1e-15 | 1e-16 | 2.2e-3 | 1.6e-2 | 6.4e-2 |

Zero to machine precision below the transition, positive above, lifting off at
WP1's `T_c` for both systems — using neither WP1's linearisation nor WP4's order
parameter. **Three independent routes to the same number.**

## Falsifiable predictions

1. ~~WP1's reweighted matrix reduces to `chygraph.percolation.PercolationMatrix`
   identically at `u' = p`.~~ **Confirmed**, and with it the Ising and AT lines
   on a configuration-model graph. See Results.
2. ~~Layer-refined `e_dd'` reproduces Fig. 1 of VW03 — `x_c(r)` at
   `gamma = 2.5` and `3.0` — from `chygraph` code with no new solver.~~
   **Half right, and the wrong half was mine.** Fig. 1 is *reproduced* (see
   Results, WP2), but "no new solver" contradicted free-transfer item 3 in this
   same README: the vertex-cover map is anti-monotone, so `Chygraph.solve`
   cannot reach its fixed point at all. A new solver was required and is in
   `vertexcover.solve`. The prediction should have read: *reproduces Fig. 1
   once the anti-monotone solver of free-transfer item 3 is supplied.*
3. ~~The hypergraph hitting-set RSB threshold moves *down* in cardinality
   heterogeneity, opposite to the degree-heterogeneity effect in VW03.~~
   **Wrong, and in the interesting direction.** It moves *up*: heterogeneity
   keeps hitting set easy longer, same as VW03's degree effect. What moves it
   down is *correlation* across cardinality layers, which the prediction did not
   consider. See Results (WP3).
4. ~~A chygraph whose complexes are the HRG's clustered motifs has a non-empty
   core where the degree-matched configuration model does not.~~ **Confirmed —
   but read the exponent claim carefully.** `core ~ kbar^theta` (θ for
   structural exponents; β is the inverse temperature) with `1.553 ± 0.035` and
   `1.573 ± 0.011` predicted against `1.590 ± 0.061` and `1.635 ± 0.008`
   measured at `tau = 2.5, 2.9`, where the degree-matched control has no core at
   all. The qualitative separation is real; the *exponent* agreement is weaker
   than the central values suggest. Bootstrapped over seeds the predicted
   exponent does not move with `tau` at all (spread 0.020) while the measured
   one spreads by 0.143 and rises monotonically, and at `tau = 2.9` the two
   differ by roughly four standard deviations. And a positive prediction follows
   from the algebra the moment one triangle enters the ensemble — Sec. VI shows
   any layer of cardinality ≥ 3 removes the core-free branch — so the
   qualitative half carries less weight than it looks. The **rewired control**,
   not the exponent, is the sharp check here. See below.

Prediction 4 is the whole motivating claim and the one to attack first, because
it is the one that would kill the programme. It is written up as a work item in
[`TODO.md`](TODO.md), including a cheap probe that could close it negatively
without building anything.

## Risks

- ~~**The HRG may not admit a useful chygraph mapping.**~~ **Measured, and it
  does** — for `tau >= 2.5`, where the maximal-clique second moment converges
  and the HRG separates from its degree-matched control by a converged factor
  of 1.8–4.8. See [`probe/RESULTS.md`](probe/RESULTS.md). The risk was real at
  `tau = 2.1`, where the ensemble does diverge — but there the *control*
  diverges faster, so that failure is the heavy tail, not the geometry. The
  remaining risk moved to step 2 of [`TODO.md`](TODO.md): core percolation may
  not be expressible as a chygraph fixed point at all.
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

## Layout

```
src/chygraph_statmech/
  stability.py   WP1  StabilityMatrix: PercolationMatrix, reweighted moments
  cavity.py      WP1  u' by exact enumeration inside a complex
  models.py      WP1  graph_percolation, graph_ising, graph_with_triangles_ising
  fixedpoint.py  WP2  Jacobian at the fixed point actually reached
  vertexcover.py WP2  VW03 Eqs. (16)-(18): the anti-monotone branch
  antimonotone.py WP3 order-reversing solver: F o F bracketing
  hittingset.py  WP3  minimum hitting set on a hypergraph
  core.py             leaf-removal core as a chygraph fixed point
  population.py  WP4  field distributions by population dynamics
  api.py              Chygraph: the single class; every method delegates below
  region.py      WP5  region graph, Mobius counting, overlap profile
  gbp.py         WP5  parent-to-child generalised belief propagation on it
  freeenergy.py  WP6  Bethe free energy; the transition thermodynamically
  ising.py            critical temperature for any chygraph
  cover.py            vertex cover of the induced graph
  softfield.py        hitting set with the O(1) fields kept (MT-validated)
  simplicial.py       the simplicial (unanimity) Ising model of Son-Lee-Goh
tests/           347 checks, each one a claim this README makes;
                 test_derivations.py covers the Supplemental Material's algebra
examples/        clustering_lowers_tc.py, vw03_figure1.py, hitting_set_rsb.py,
                 wp4_validates_wp1.py, softfield_vs_hardfield.py,
                 gbp_two_triangles.py
book/            the book that supersedes the manuscript; book/README.md is
                 the live working document, book/figs/ one script per chapter
TODO.md          open items; prediction 4 on hyperbolic random graphs
probe/           HRG clique moments, prediction 4, hitting-set density,
                 GBP on clique region graphs; see probe/RESULTS.md
```

Install with `pip install -e .` (needs `chygraph`, `numpy`, `scipy` on the
path); `pytest tests`.

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
