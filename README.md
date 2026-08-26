# chygraph_statmech

Statistical mechanics on chygraphs: extending the Bethe–Peierls treatment of
lattice-gas models on correlated random graphs
([Vázquez & Weigt, PRE **67**, 027101 (2003)](https://arxiv.org/abs/cond-mat/0207035))
from graphs to higher-order structures, by way of the chygraph formalism
([Vázquez, PRE **107**, 024316 (2023)](https://arxiv.org/abs/2308.00987);
`~/av2atg/chygraph`).

Status: **WP1–WP5 done** (WP5 counting only, no GBP), plus core percolation (`src/chygraph_statmech`, 56 tests). WP4–WP6 and the
HRG test of prediction 4 ([`TODO.md`](TODO.md)) are plan. Results:
[WP1](#results-wp1) · [WP2](#results-wp2) · [WP3](#results-wp3).

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

**WP5 — Complexes as regions.** ✅ Done for the counting; GBP on the region
graph is not implemented. `region.py` builds the region graph and measures how
far a chygraph is from treelike; the cost of ignoring overlap is measured on
HRGs. See [Results (WP5)](#results-wp5).

**WP6 — Bethe free energy.** VW03 Eq. (12) has no chygraph counterpart. On a
chygraph the counting is complex-terms minus inclusion-overcounting, i.e. a
region-graph free energy. Needed for anything variational; not needed for
thresholds.

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

the third term vanishing for Poisson layers. Holding the mean number of
neighbours fixed and moving them from links into triangles raises `T_c`
(`examples/clustering_raises_tc.py`):

| excess neighbours | as links | as triangles | `t_c` links | `t_c` tri | `T_c` gain |
|---:|---:|---:|---:|---:|---:|
| 4 | 4 | 2 | 0.25000 | 0.20871 | +20.6% |
| 6 | 6 | 3 | 0.16667 | 0.14590 | +14.5% |
| 8 | 8 | 4 | 0.12500 | 0.11252 | +11.2% |
| 10 | 10 | 5 | 0.10000 | 0.09167 | +9.1% |
| 20 | 20 | 10 | 0.05000 | 0.04773 | +4.8% |

These configurations are **identical** in `{p_d, e_dd'}`: same degree
distribution, same assortativity. The effect is `O(t²)`, hence the `1/<k>` decay.

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

### Not done

Generalised belief propagation on the region graph. `region.py` builds the
object and measures the discrepancy; it does not pass messages between regions.
That is what would close the 17–35% gap, and it is the natural next step.

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
4. A chygraph whose complexes are the HRG's clustered motifs has a non-empty
   core where the degree-matched configuration model does not — i.e. the
   formalism sees the effect the `{p_d, e_dd'}` ensemble misses.
   *Still open, but no longer blocked: the ensemble exists for `tau >= 2.5`
   ([`probe/RESULTS.md`](probe/RESULTS.md)), so the comparison can be run.*

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
  region.py      WP5  region graph, Mobius counting, overlap profile
tests/           116 checks, each one a claim this README makes
examples/        clustering_raises_tc.py, vw03_figure1.py, hitting_set_rsb.py,
                 wp4_validates_wp1.py
main.tex         manuscript in progress; results in place, prose to write
TODO.md          open items; prediction 4 on hyperbolic random graphs
probe/           HRG clique-moment measurement; see probe/RESULTS.md
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
