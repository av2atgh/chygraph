# Extensions: what the book should state that is new

Written 2026-09-24 from a reading of the outlook against the rest of the book.
The book's claim as it stands is methodological: one recursion, many models,
and a measured boundary where it fails. The three items below are, in rank
order, what would turn that into a stake in the ground: (i) an unsolved
problem solved, (ii) a new problem framed, (iii) two branches brought together.

## 1. An exactness theorem plus a new transition (thesis)

Chapter 15 finds the instance calculation exact "where the clique family is
chordal"; Chapter 16 finds merging exact "iff the merged incidence structure
is a forest". Both are one classical fact stated twice without its name:

- a graph is chordal iff its maximal-clique hypergraph is α-acyclic iff it has
  a join tree (Gavril 1974; Beeri, Fagin, Maier, Yannakakis 1983);
- belief propagation over a join tree is exact (Lauritzen and Spiegelhalter
  1988);
- GYO reduction, the certificate of α-acyclicity, is leaf removal on the
  incidence structure. The tool of Ch. 11 decides when Chs. 4–13 are exact.

The statement the book can make: **the chygraph cavity method is exact on an
ensemble iff its complex family is α-acyclic w.h.p., and the loss of
α-acyclicity is a percolation-type transition ("acyclicity percolation") whose
order parameter is the GYO residue**, sitting beside core percolation at
⟨k⟩ = e. That is (i), (ii) and (iii) in one move, and it answers the hierarchy
question of Sec. 17.5 for the ensembles where the answer is "one rung".

The ensemble to carry it cannot be Karrer–Newman: two placed triangles share
an edge O(1) times in the whole graph, so that ensemble is exact at every
branching ratio in the limit; the b = 1 threshold of Sec. 16.3 is the giant
component of the incidence structure, not the failure of the cavity method.
Extensive overlap needs geometry:

- random interval graphs (1D random geometric graphs) are chordal by
  construction, so the method is provably exact there for Ising, cover,
  colouring at any density, on a network with clustering of order one;
- this explains Ch. 15's "chordality follows provenance": hyperbolic random
  graphs are nearly one-dimensional in the angular coordinate;
- random geometric graphs in dimension d interpolate; the acyclicity
  transition against d and density is a well-posed unsolved problem with one
  end already answered.

Reframes Part IV from "where the method fails" to "the theory of when it is
exact".

## 2. Derive the hyperbolic random graph chygraph instead of measuring it (worked example)

Sec. 11.7 is the closest thing to a solved open problem, but its ensemble is
fitted from realisations and the exponent near 1.57 is admitted to be weak.
Bläsius, Friedrich and Krohmer (Algorithmica 2018, "Cliques in hyperbolic
random graphs") give the expected number of k-cliques and the clique number
of HRG as functions of τ and density: the cardinality and chy-degree layer
data of the clique chygraph in closed form. Plugging it into Eq. (11.4) makes
the decade-old failure of degree-based mean-field theory on HRG an analytic
prediction, exponent included. Bounded risk; feeds item 1 since the same
geometry controls the overlap.

## 3. State the replica conjecture as a conjecture, with its one solid result

The book "detects RSB and does not do it". `~/av2atg/replicas` holds a new
framing: a solution cluster's interior factorises into a frozen core plus
independent free blocks, so a cluster is itself a chygraph, and the sharpest
conjecture is that the Parisi m is a tree artefact whose chygraph analogue is
the structural ratio ⟨κ_l⟩/c_l. The withdrawn trend stays withdrawn; the
factorisation stands on exhaustive enumeration plus a control. A short section
in Ch. 8 stating the conjecture with its testable consequence. After 1 and 2,
not instead of them.

## Smaller bridges

- Dynamic message passing on chygraphs: the interior sum becomes a
  within-complex time course; closes the "no dynamics" gap of Sec. 17.4.
- Inference: the vertex-cover decomposition in `chygraph/vc` already treats
  SBM blocks as complexes.

## Recommendation

Item 1 as the book's thesis, item 2 as its worked example.

## Status 2026-09-24 (evening)

Work started on items 1 and 2. Draft chapter `book/exactness.tex` (not yet
included from `main.tex`); probes `statmech/probe/acyclicity.py` and
`real_acyclicity.py`; references appended to `references.bib`.

**Item 1, corrections found while measuring.**

- GYO residue (α-acyclicity) is an *instance* notion and keeps every cycle,
  long ones included. On real networks it is 23–100 % of the vertices and one
  connected piece (power grid 61 %), so it cannot be the ensemble order
  parameter. The ensemble notion is local: the fraction φ of vertices on a
  chordless cycle of length 4 or 5 (the loops between complexes that no
  complex contains), and "acyclicity percolation" is the percolation of that
  set. Both quantities are kept; the chapter says which is which.
- Radius-1 ball chordality ("wheel around v") is not Ch. 14's criterion
  (two complexes sharing ≥ 2 atoms): two triangles sharing an edge are
  chordal and the join tree handles them; a wheel with ≥ 4 spokes is not.
- 1D geometric graphs must be on a line: the ring closes one long chordless
  cycle that GYO keeps.
- On HRG the Ch. 16 merge closure produces a giant meta-complex (37–45 % of
  the vertices at τ = 2.5, k̄ = 4, n = 20000) while φ is 2 % and the GYO
  residue 3–4 %, both one connected piece: the centre. So the join-tree
  route, not merging, is the repair on HRG; Ch. 16's "always a forest" was
  measured on ego-networks, not whole graphs.
- 2D RGG: φ > 0 at every density (8.6 % at k̄ = 3, 43 % at 6, 84 % at 10),
  and the short-cycle set percolates well above the graph's own threshold
  4.51: giant fraction 0.1 %, 1 %, 84 % at k̄ = 3, 6, 10. There is a window
  above k̄_c where the graph has a giant component and the obstruction is
  finite islands. Finite-size sweep running to locate k̄_a.

**Item 2, what the literature actually gives.** Bläsius–Friedrich–Krohmer
give E[K_k] (all k-cliques, not maximal ones): n Θ(k)^{-k} for β ≥ 3, plus
n^{k(3−β)/2} Θ(k)^{-k} from the centre for β ∈ (2,3), and ω = Θ(n^{(3−β)/2})
resp. Θ(log n / log log n). Yamaji (arXiv:2303.06301) shows the number of
*maximal* cliques is exp(Θ(n^{(3−γ)/6})), all in the centre. So the layer data
of the clique chygraph (chy-degree by cardinality) is not in closed form
anywhere; it has to be derived from the geometry. At the periphery two nodes
connect iff |Δθ| < 2 e^{(R − r_i − r_j)/2}, a 1D graph with a *product*
kernel w_i w_j, which is not an interval graph (additive kernel), so
"nearly 1D ⇒ chordal" is a heuristic, not a theorem; the measured φ says the
periphery is exact and the centre is not. Item 2 is a derivation project, not
a lookup.

## Status 2026-09-24 (night): item 1 carried out as Ch. 17

`book/exactness.tex` is written, included from `main.tex` before the
Outlook, and built; `outlook.tex` and `metacomplex.tex` carry the
cross-references; `README.md` has the row and the figure-script paragraph.

Final numbers (n = 20000 unless stated; `statmech/probe/results/`):

- Karrer–Newman s = t = 1: vertices on short chordless cycles 157, 149, 129
  at n = 5000, 20000, 80000 — φ ∝ 1/n, exact at every b.
- RGG d = 1 (line): φ ≡ 0 at every density; merge closure largest 13 atoms
  at k̄ = 6.
- RGG d = 2: φ = 0.081, 0.17, 0.28, 0.42, 0.84 at k̄ = 3, 4, 5, 6, 10;
  largest island vanishes with n up to k̄ = 7, crosses at 8 (0.54, 0.54,
  0.29 for n = 5k, 20k, 80k), flat at 9 (0.75); k̄_a a little below 8 vs the
  graph's 4.512. Islands 11, 19, 48 at k̄ = 4, 5, 6.
- RGG d = 3: φ = 0.048, 0.14, 0.29, 0.62 at k̄ = 2, 3, 4, 6; island vanishes
  at 4.5, giant 0.60 at 6; k̄_a ∈ (4.5, 6) vs the graph's 2.736.
- HRG τ = 2.5: φ = 0.0002, 0.006, 0.028, 0.11 at k̄ = 1, 2, 4, 8, one piece
  (the centre); merge closure 0.03, 0.16, 0.46, 0.84. τ = 4: same φ, islands
  of 8 and 60 (no centre to gather them).
- Real networks (LCC, 5 rewirings): grouped ≪ control (co-authorship 0.077 vs
  0.77), interactomes ≈ control (yeast 0.18 vs 0.25), spatial ≫ control
  (power grid 0.29 vs 0.02, islands of 187). GYO residue 0.13–1 everywhere;
  wheel < 0.03 on every sparse network.

Open, in order: (a) the ensemble calculation past k̄_a (message indexed by
separator type over a random join tree) — now with the distribution of
separators and islands measured; (b) item 2, the HRG layer data from the
geometry; (c) item 3, the replica conjecture as a section of Ch. 8;
(d) a d = 2 finite-size scaling of k̄_a proper (three sizes here, no
collapse); (e) whether φ on the interactomes vanishes in a growing ensemble
with the same degree tail, which the rewired control suggests and does not
prove.

## Status 2026-09-25: item 3 carried out as Sec. 8.10

The original `~/av2atg/replicas` notes are no longer on disk (only the memory
note survives), so the section rests on a fresh enumeration:
`statmech/probe/cluster_blocks.py` (exhaustive 3-SAT, N = 16–24, α = 3.5,
3.8, 4.0, 30 instances each, 415 satisfiable of 450; clusters at Hamming
distance one; frozen set, pairwise-dependency blocks, exact product test,
size-matched random control) summarised by `book/figs/clusters.py`.

Numbers: 1219 clusters of ≥ 2 solutions, mean size 22, frozen fraction 0.77;
478 with ≥ 2 candidate blocks, of which 0.81 are exact products (0.80–0.82
across N, 0.79–0.85 across α, no trend); control: 0.66 of random sets split
pairwise, 0.25 of those are products; 2.3 blocks per multi-block cluster,
0.55 blocks per free variable.

Written into `statmech.tex` as Sec. 8.10 "A cluster is itself a chygraph":
the measurement, what it means for Eq. (rsbansatz), the conjecture (m is a
counting number ⟨κ⟩/c one level up) with its test (m from the m ≠ 0
one-step calculation against the measured block density), and the two things
the enumeration cannot say (inter-cluster geometry; the thermodynamic limit).
The withdrawn WalkSAT-vs-tempering trend is not mentioned: its data are gone.

Open: the m ≠ 0 population dynamics on 3-SAT that would test the conjecture;
a higher-order (beyond pairwise) block test for the remaining fifth.

## Status 2026-09-25 (later): the two open items of Sec. 8.10

**Beyond-pairwise block test.** `cluster_blocks.py` now also finds the
finest product partition (product decompositions are closed under common
refinement, so it exists; it is a merging of the pairwise candidate blocks,
found by testing the 2^k unions for separability). Of the non-product
fifth, 0.19 split at a coarser grain, against 0.17 of the random controls;
so 0.84 of multi-block clusters are products at the finest grain and the
rest are one coupled block. Written into Sec. 8.10.

**m ≠ 0 for 3-SAT.** `statmech/probe/onestep_sat.py`: entropic one-step
calculation as a population of populations (M = 1000 surveys × P = 200
messages, reweighting z^m at every update, potential by importance
weighting over 4e5 draws × 400 combinations). Anchors: BP = exact on a tree
instance to 1e-6; m = 1 potential = F_RS within population scatter (0.01);
below α_d the surveys collapse and F(m) = m F_RS. Σ(0) ≡ 0 with soft
messages (no anchor there). Results (Sec. 13.10, Fig. 13.2): Σ(m) of order
1e-3 in the condensed phase; Σ(1) first clearly negative at α = 4.0; m*
bracketed: [0.7, 0.9] at 4.0–4.05, [0.8, 0.9] at 4.1, [0.5, 0.6] at 4.15,
[0.6, 0.7] at 4.2, not monotone — noise (population scatter 8e-4, worst
4e-3, on a few 1e-3). Beyond m ≈ 0.6–0.9 at α ≥ 4.05 the soft population
freezes (spread < 0.08, Σ ~ 0.1–0.8 or divergent) and is discarded; at
4.25 everything above m = 0.3 freezes, so m* → 0 at α_s is not reached.

**The test.** Block density 0.55 (flat in α and N at N ≤ 24) vs m*: equal
only inside the 4.15 bracket, excluded at α ≤ 3.95 (m* = 1) and at 4.0–4.1.
Conjecture in the "m* is the block density" form rejected; what survives is
that the decomposition is structural, the exponent is not. Sec. 8.10 and
the Outlook say so.

Open: the hard fields at m → 0 (explicit delta weight in the survey, as
Montanari–Ricci-Tersenghi–Semerjian do) — that would make the soft
population stop freezing and reach α_s; k ≥ 4 where α_c ≠ α_d; a block
density at large N by equilibrium sampling, which is the only way the
conjecture's test could be reopened.

## 2026-09-25: Ch. 18, localization on a clustered graph (in progress)

Prompted by Tonetti, Cugliandolo, Tarzia, arXiv:2512.04037 (LLT vs Anderson
on the K+1 = 3 Bethe lattice; the lines differ, blamed on the absence of
loops). The chygraph version: the interior of a complex for the resolvent is
a Schur complement, Σ_{a→i} = t² 1ᵀ[(ε_J − z) − tA_J − Σ_{J→a}]⁻¹ 1,
polynomial in the cardinality; the landscape follows by linearity; LLT
percolation is a dependent layer with a threshold rule; the mobility edge is
the linearised imaginary map's growth rate min_β λ(β) = 1.

`statmech/probe/anderson.py` (checks, populations, `scan`); `book/figs/anderson.py`;
`book/localization.tex` drafted (theory sections; results pending the scan).

Verified: cavity exact on incidence trees (1e-15); random instances' error
falls with Im z (loops); landscape vs direct solve 1e-8 (trees), 3e-3 (densest
triangle ensemble, n = 3000); band bottoms −2√2, −2√3; band-centre λ(1/2) at
W = 18 rises 0.967 → 0.977 → 0.986 with P = 2e5, 1e6, 3e6 (W_c = 18.17 known);
(3,0) at W = 1.5: E_c^loc = −3.106 (paper −3.152), E'_c = 0.473 →
E_c^perc = 3.106 on the symmetric side (paper 3.165). LLT percolation also
checked on a 60000-site instance with an exact (CG) landscape: onset at
E' ≈ 0.45 (W = 1.5) and ≈ 2.1 (W = 6) vs population 0.47 and 2.30.

Lesson: the (Σ, η, p) triple must be regenerated jointly in the percolation
population; separate regeneration loses the smoothness of u and shifts the
threshold by several per cent.

First triangle result: (1,1) at W = 1.5: E_c^loc = −3.207, E_c^perc = −3.011,
gap 0.20, where (3,0) has gap 0.006. Triangles open the gap at a disorder
where the tree has none. Full scan (5 ensembles × 6 W, plus W_c at E = 0)
running.

## 2026-09-25 (night): Ch. 18 finished

`book/localization.tex` complete; `figs/localization.py` (renamed from
figs/anderson.py to avoid shadowing the probe module) draws Fig. 18.1 and
writes Tables 18.1–18.2. Scan: `statmech/probe/results/anderson_lines.csv`,
`anderson_wc.csv` (5 ensembles × 6 W; 1e6 elements for the edge, 3e5 for
percolation; ~3 h on 8 cores).

Results (bottom of the band; gap = E_c^perc − E_c^loc, positive = LLT misses
extended states):
- (3,0): gap 0.006, 0.30, 0.60, 0.88, 1.44, 1.37 at W = 1.5 … 12 — the paper's
  picture (lines together at weak disorder, apart at strong).
- (4,0): gap −0.42, −0.21, +0.03, +0.27, +0.75, +1.29 — the crossing near 4.3.
- (1,1): gap 0.20, 0.47, 0.74, 0.88; at W = 9, 12 the whole band is localised
  (W_c = 8.6 vs 17.4 for (3,0)) while the landscape still percolates.
- Degree 4: crossing at W ≈ 4.3 (4,0), 3.0 (2,1), 1.7 (0,2); gap at W = 6:
  0.27, 0.49, 0.75; at W = 12: 1.29, 1.49, 1.92. Percolation line moves
  toward the centre with triangles (0.4 at W = 6), mobility edge barely.
- W_c(E = 0): 17.4, 8.6, 32.3, 26.2, 18.1 — the √2 branching argument:
  (0,2) ≈ (3,0).

Open: exponents (not computed); explicit hard fields for the imaginary parts
would reduce the population bias; real-network clique chygraphs by layer.

## 2026-09-26: Ch. 18 reworked around figures

Table 18.1 (thirty rows of energies) dropped. The chapter now carries five
figures: the Schur-interior schematic (TikZ), the two lines measured from
the band edge (removes the −W/2 both share), the gap against W in one
panel, the band-centre growth rate λ(1/2) against W for the five ensembles
(the (0,2) curve lies on the (3,0) curve: the √2-per-triangle argument in
one picture), and the two criteria at work (λ(1/2) against E − E_min at
W = 3; percolation growth against E' with the triple drawn jointly and
with p drawn separately: threshold 0.47 vs 0.52). Side data from
`anderson.py side` and `anderson.py growth`; the latter uses a short
iteration from a uniform p because renormalised power iteration goes
extinct in the subcritical phase of a finite population. Error bars from
seven points computed twice: 0.002 on the threshold, 0.003–0.03 on the
mobility edge.

## 2026-09-25: the resolvent is a first-class message

`statmech/src/statmech/resolvent.py` now holds the message type of Ch. 18:
`transmit` (the Schur-complement down step, in Sherman–Morrison closed form
Σ = t²S/(1−tS), S = Σ_j 1/(m_j+t), linear in the cardinality; also returns
v = M⁻¹1 and the landscape row sum), `linearised_imaginary` (Eq. imlinear),
`cavity_instance`, `landscape_instance`, `hamiltonian`, `regular_instance`,
`incidence_tree`, `uniform_fixed_point`, `band_bottom`, `spectral_bottom`,
`clean_tree`, `kesten_mckay`, `incidence_branching`.  `Chygraph` gained
`instance`, `incidence_branching`, `band_bottom`, `spectral_bottom`,
`resolvent_instance`, `landscape_instance`.  `tests/test_resolvent.py` pins
all of it (26 tests).  `statmech/probe/anderson.py` keeps only the population
dynamics (`Ensemble`) and the scans, and routes every interior through
`transmit`; `anderson.py check` reproduces its numbers exactly and the
population regression at (3,0) W=18 agrees with the old code to 7 digits.
`figs/localization.py` imports the package, not the probe.  Sec. 18.2 gained
Eq. (sherman) and the cost sentence now says linear; the software chapter has
a Ch. 18 block and handle rows.  Note for the rerun: λ(1/2) at (3,0), W=18,
P=2e5, seed 0 is 0.9775 with both old and new code, not the 0.967 the text
quotes, so that number came from another run.

## 2026-09-25: Ch. 18 rerun from one code version (review item 2)

`anderson.py rerun` (108 min on 8 cores): scan at β = 1/2 only, repeats at
seed 1 (six points) and seeds 2–3 (cactus W = 6), the bias triple at three
seeds, a β scan 0.3–0.7 on tree and cactus, and the side data.  Percolation
thresholds identical to the previous scan at all 30 points (their scatter is
below the bisection grid of 0.003–0.007); mobility edges move by ≤ 0.09;
W_c unchanged to 0.1.  Seed scatter of the mobility edge: 0.005 at W = 1.5,
0.03 at 4.5, 0.06 on the cactus at W = 6 (−3.72, −3.69, −3.66, −3.66).
Bias triple: λ(1/2) = 0.98 ± 0.01 at 2e5, 1e6, 3e6 alike; no resolvable
trend with P.  β scan: flat to 0.002 on the cactus at its band-centre W_c;
on the tree the β-dependence (0.03) changes sign between seeds.  Secs.
18.4–18.7 rewritten from these; new CSVs anderson_repeats/bias/beta.

## 2026-09-25: Part V, the math connection

Six chapters written from `~/Downloads/chygraph_master_equation/draft.tex`
(the "master equation" note): Secs. 1–5 → Ch. 19 `chyequation.tex`; Secs.
6–10 → Chs. 20–24 `rde.tex`, `species.tex`, `operads.tex`, `marginal.tex`,
`mobius.tex`, one field each. Placed after Part IV, before the Outlook, which
still closes the book (now Ch. 25). Part V computes nothing and carries no
Checks sections; four TikZ figures (the factor graph of α+I, the two holes,
cyclic vs modular gluing, the Möbius numbers on two triangles).

Preamble gained `\Tint`, `\conv`, `\act`, `\push`, `\del`; `references.bib`
gained thirty entries (the draft's, checked; `foissy2011` renamed
`foissy2010`, its year; Peltre's thesis title and GSI pages least certain);
Notation's χ row lists the Euler characteristic and the kernel χ_a.

Corrections made to the draft while anchoring it on the book:
- Ch. 24: "exactness is collapsibility" weakened to an implication (the
  four-spoke wheel is collapsible and not α-acyclic); "exactness is b₁ = 0
  of the incidence structure" is the chygraph recursion's own condition
  (Ch. 16's forest), strictly stronger than Ch. 17's α-acyclicity (two
  triangles sharing an edge: α-acyclic, b₁ = 1); the ordering is
  b₁ = 0 ⇒ α-acyclic ⇒ collapsible with both converses false. The
  linearised recursion is a weighted non-backtracking operator, not a
  sheaf Laplacian; both come from one sheaf.
- Ch. 23: the sheaf class and the BP error are different diagnostics (a
  family can glue and still be wrong); whether Part IV's fixed points glue
  is not checked.
- Ch. 20: the AT-tensor = bivariate-uniqueness identification is stated for
  the linearised problem only; the two-population nonlinear test is the
  one computation Part V calls for.

Open, in order of cost: (a) the bivariate-uniqueness populations on the
models of Table 25.1; (b) r_C for the two-triangle loop from the loop
series, against Table 14.1 (the resummation claim of Sec. 24.2); (c) the
finite-component distribution by layer from Good's inversion, Eq. (21.7),
against simulation; (d) the level-wise Vorob'ev proviso (consistency on the
nested family vs on the flattened one); (e) whether Foissy's classification
covers the (layer, direction) decoration.
