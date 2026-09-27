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

## 2026-09-25: Part V, the math connection (superseded on the computing: see the three entries below, Part V now computes three things and Chs. 20, 21, 24 carry Checks)

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

## 2026-09-25: Part V reviewed (`~/Downloads/chygraph_book_review.md`)

Applied. A1: "exact when the factor graph is a tree" is Ch. 16's forest =
Berge-acyclicity (Fagin 1983), stated as such in Chs. 19, 22, 23, 24; Ch. 17's
three names are α-acyclicity and belong to the region-graph calculation. A2:
the Vorob'ev paragraph now runs on the region-graph fixed point (chygraph BP
agrees on single atoms only), gluing gives *a* law, Boltzmann needs the
junction-tree factorisation + Lauritzen–Spiegelhalter; "consistent with" the
240 runs. A3: the level-wise statement uses α-acyclicity of the flattened
family and is a corollary of Vorob'ev–Kellerer modulo the consistency-notion
check, which reduces to overlaps of size ≥ 2; no longer called a conjecture
(also in the Outlook). A4: smoothing transform — the α-stable fixed points
live on the localised side (min λ ≤ 1, characteristic exponent β* < 1/2),
the extended side is outside the class; AT line is m(1) = 1 at fixed exponent
on D² with weights u'², not a min over β. A5: Abramsky obstruction is
one-way (abramsky2012 for the class, abramsky2015 for paradoxes); support
bound, not dimension. A6: Aldous–Bandyopadhyay Thm 11(c) cited as the test.
B: arity c−1 + own leg; û' ≡ 1 / absent / (0,1) for the three cases; node vs
edge perspective = Φ vs Φ̄; Aldous–Lyons for unimodular; DMS conditions on
the fixed point; ν mass zero; "heuristically" on the power law; chygraph
readout for the component distribution; stochastic species: composite
exists (Pitman 2006, Diaconis–Pang–Ram 2014), fixed-point question new;
Giry 1982; character claim restricted to Com; coloured operad on the
positive cone; atoms as Com vertices explain the genus count; Getzler–
Kapranov 1995; Vallette 2007 and the properad section shortened
("disconnected operation"); hierarchical hypergraphs' strict nesting
flagged; two-triangle r_C is an identity check, the resummation claim is
about two series around different fixed points; Linial–Meshulam comparison
is with the giant of edge-adjacent triangles at c = 1/2; sheaf: the
F(j→a)→F(a)→F(a→i) factorisation named as the thing to check, Watanabe–
Fukumizu 2009 and Saade et al. 2014 cited; Peltre 2019 added. C: ∂ for
derivatives, `\del` only for members; χ no longer used for the Euler
characteristic (unlettered), kernel χ_a kept since Eq. (Xichy) already uses
χ_l; every "acyclic" qualified.

Open calculations, re-ranked by the review: (1) the two-population endogeny
test, Thm 11(c), hitting set first; (2) the finite-component distribution by
layer from Good's inversion against simulation; (3) the consistency-notion
check on overlaps of size ≥ 2 (turns the level-wise statement into a stated
corollary); (4) whether Part IV's region-graph fixed points are
Vorob'ev-consistent on pair overlaps (one-line check on the Ch. 15 runs);
(5) the term-by-term question of Sec. 24.2.

## 2026-09-25: the component distribution by layer, computed (Sec. 21.4)

`percolation/src/percolation/components.py`: the marked map of Eq. (21.6)
solved on a circle |x_l| = r < 1 (Cauchy's integral on a grid, forward
FFT / N / r^s; joint over two layers by a 2D grid), plus two exact routes
for polynomial pgfs (`symbolic_series`, rational iteration; `good_by_layer`,
Good's inversion with the rational determinant series-expanded). Tests
`tests/test_components.py` (5). `book/figs/components.py` (45 s): Borel to
1e-15; series = transform to 1e-17; Good 2/35, 6/245, 12/1225; households
{1:1/5, 3:1/2, 5:3/10}, p_H = 1/2, k = 2, T_c = 0.1569: theory = simulation
(12 × 1e5 households) at every s to 60, S 0.4895 vs 0.4893, households
touched 2.30 vs 2.31 at s = 8; nine real networks (Table 21.1): binary
interactomes S within 3% (clique closer on yeast), finite tail low by 2–4×
from s = 4; Collins (AP-MS) fails outright (S 0.99–1.00 vs 0.62) because its
finite components are isolated complexes — the cardinality–chy-degree
correlation of Sec. 5.5, on data. Marking does not support site thinning.

Open: measure the joint (cardinality, member chy-degree) distribution on the
interactomes and run the marked map on Sec. 5.5's joint construction.

## 2026-09-25: the loop series, computed (Sec. 24.3)

`statmech/src/statmech/loopseries.py`: `BinaryFactorGraph` (BP + exact Z
on any binary factor graph; constructors `pairwise` and `promoted` with
`assign='once'` — every bond in one clique, the true model on α+I — or
`'all'`, the recursion's double count), `generalised_loops` (backtracking
with last-factor pruning), `loop_term`, `loop_series`, `partial_sums`.
Tests `tests/test_loopseries.py` (6). Identity to 1e-15 everywhere.

Results (`figs/loopseries.py`, ~2 min): rings k = 4–6: promoted 1 term =
Table 16.1's error; pairwise 103/323/1019 terms, triangles 0.81–0.92 of the
sum at βJ = 0.3, 0.32–0.36 at 0.8. Forty clustered n = 10 instances: promoted
error 0.60× / 0.49× the pairwise, 285 vs 814 loops, smaller on 40/40.
Decomposition ln Z_BP^dc − ln Z = (ln Z^dc − ln Z) + (ln Z_BP^dc − ln Z^dc):
two triangles βJ = 0.5: +0.09 = +0.41 − 0.32; hyperbolic n = 14 (20 fresh
draws, own sampler; `computational_complexity/code/hrg.py` is gone from
disk): +1.29 = +2.02 − 0.73 (0.3), +5.54 = +6.27 − 0.73 (0.8); opposite
signs on all; **the double count dominates**. Assigning each bond to one
clique (BP on α+I for the true model) cuts the mean error to 0.32× / 0.16×
of the recursion's, smaller on 15/20 and 19/20 — a free repair Part IV did
not have. Caveat: on two triangles at weak coupling the assigned error is
larger than the recursion's because the two parts cancel there.

Bug found and fixed on the way (in the new module only; the probe's
`ChygraphBP._factor` was checked on asymmetric bond sets and is right): a
C-ordered 2^n table has position p on bit n−1−p.

Open: the ensemble loop series (terms indexed by the loops Ch. 17 counts);
the Bethe-Hessian check of Sec. 24.5; Part IV's tables rerun with bonds
assigned once (a one-flag change in cavity_clique.py) — the cheapest
follow-up and the one that changes Ch. 14's story.

## 2026-09-26: the endogeny test, run (Sec. 20.5)

`statmech/src/statmech/endogeny.py`: `BivariateSpinGlass` (±J Ising on a
Poisson/regular graph in a field, two copies sharing containers, neighbours
and coupling signs), `BivariateHittingSet` (the soft-field population of
Ch. 10 doubled, shared draws and shared damping masks), `endogeny_boundary`.
Tests `tests/test_endogeny.py` (6, ~90 s). Probe `probe/endogeny.py scan`
(~70 min, cache `results/endogeny.json`); `figs/endogeny.py` draws.

Results. Spin glass, degree 3: linearised boundary (perturbation carried
with the field) → exact T_c with a H^(2/3) law (a = 0.78, residuals 0.003);
shuffled-pair boundary 2–5% below it (fixed sweep budget); the annealed
line ⟨k̄⟩⟨u'²⟩ = 1 is 25–60% above and does not tend to T_c (1.52 at
H = 0.02): Sec. 9.7's substitution is exact only at the trivial fixed point.
Hitting set, Poisson c = 2, 3, 4: the stationary rate collapses in
⟨k⟩(c−1)/e and crosses zero at 1.034, 1.040, 1.029 × the hard-field line
(vertex cover: 2.81 vs e), entropy still positive there — the test sides
with the hard-field line; μ-independent (30–120), damping-independent,
plateau size-independent (5e4–4e5). Regular (Mézard–Tarzia, exact
entropy): K = 3 with s < 0 non-endogenous (distance 0.8–1.0, rate +0.08);
K = 4, 6 with s < 0 (L = 3,4) endogenous to the last digit, rate −0.06 to
−0.12: **endogenous and wrong** — the discontinuous case; only the entropy
sees it.

Method notes: the rate must be read from long stationary windows (short
windows read the transient and put the vertex-cover crossing at 2.95); the
"linear" perturbation in a near-hard-field population is flip-dominated and
noisy, so the two halves of the window are compared; the shuffled criterion
under-reads the boundary where decay is slow.

Open: the ensemble loop series; the Bethe-Hessian check of Sec. 24.5;
Part IV rerun with bonds assigned once.

## 2026-09-26: second Part V review (`~/Downloads/book_review.md`) applied

A1 Ch. 17's conclusions and opening now say the join-tree calculation is
exact iff chordal, the recursion on atoms iff Berge-acyclic. A2 Sec. 21.4's
interactome paragraph rewritten from Table 21.1 (giant within 4%, Vidal the
worst; small components scattered both ways; Figeys tenfold; road network a
10% miss, tied to Ch. 17's geometry). A3 hitting-set pair followed to 2000
sweeps at 0.95 and 1.00 of the hard-field line: converges at the line
(8e-6, 4e-5, 4e-5, rates −0.003 to −0.006/sweep), so the 600-sweep reading
was slow decay, not a plateau; Fig. 20.1b cut at 1.2 with a caption note on
saturation; Table 20.1 carries the entropy scatter, vertex cover's 0.11±0.08
called consistent with zero. A4 citations (Montanari–Ricci-Tersenghi 2003,
Rivoire et al. 2004, Zdeborová–Krzakala 2007, Bandyopadhyay 2011) and the
header comment; what is new stated narrowly. A5 "necessary for endogeny",
chain through 11(c) then 11(b). A6 β* ≤ 1 condition (Durrett–Liggett 1983),
inhomogeneous equation and Kesten–Goldie tail. A7 κ_a + 1 hyperedges. A8
Euler sum stated for two levels / the factor-graph count. A9 "intermediate
coupling". A10 truncation at s = 60 said. B19 pendant complex variables for
Part III; exactness sentence as a definition. B21 typo, weights on
structures, household check in the software table. B22 cyclic operad needs
an exchangeable kernel. B23 Kellerer → Math. Ann. 153 (1964). B24 Pakzad–
Anantharam for c_R = −μ; assignment spread measured (hyperbolic 0.10/0.19
against 0.46/0.76; clustered 0.04/0.60 against 0.15/1.02 — comparable at
strong coupling); the split run on Ch. 14's 120 real neighbourhoods
(`probe/loopseries_real.py`): the loop part is −ln 2 to three decimals on
five of six networks — the polarised sector, checked: assigned-once ln Z_BP
= exact ln Z with the centre fixed up to four digits — so beyond ln 2 the
real-network error was the double count and nothing else; Ch. 14's
conclusion and the Outlook's ladder now say so. C forward references from
Chs. 6, 10, 14, 17; figure widths 0.76; first Part V entry marked
superseded; marginal.tex comment fixed. Not done: the Bethe-Hessian check
(Sec. 24.5); an optimised bond-assignment rule.

## 2026-09-26: the Bethe Hessian, computed (Sec. 24.6; the notes above call it 24.5, before Collapsibility was inserted)

`statmech/src/statmech/bethehessian.py`: the linearised recursion in field
coordinates on the factor graph of α+I (`jacobian_block`: M^a_vu =
∂h_{a→v}/∂h_{u→a} = Cov_a(σ_v,σ_u)/(1−m_v²), by enumeration;
`linearised_operator`: non-backtracking on the incidence graph with a
*block* per complex), the product-form test (`factorise`, alternating least
squares off the diagonal), the Bethe Hessian on atoms + complexes from the
weighted Ihara–Bass identity (`bethe_hessian`), its Schur complement on the
atoms (`vertex_hessian`, `trivial_vertex_hessian`), `instance_threshold`
(bisection on ρ(T) = 1 or λ_min(H_V) = 0), `random_chygraph` (Poisson or
exactly regular layers). Tests `tests/test_bethehessian.py` (7, ~1 s).
`figs/bethehessian.py` (~100 s, cache `probe/results/bethehessian.json`).

Results. (1) At the trivial point of a homogeneous complex every block entry
is u'(c) (c = 2..5), so T is edge-weighted non-backtracking and
Watanabe–Fukumizu's identity gives det(I−T) = ∏(1−u')^{c−1}(1+(c−1)u') det H_V
with H_V = I + Σ_a u'/(1−u') [P_a − 1_a1_aᵀ/(1+(c−1)u')] — one rank-one term
per complex; on a graph this is t²/(1−t²) × Saade's (r²−1)I − rA + D at
r = 1/t (deviation 0). (2) Ferromagnet: H_V(0) = I and det > 0 until the
Perron root hits 1, so the trivial fixed point is stable iff H_V ≻ 0; the
instance threshold by ρ(T) and by λ_min(H_V) agree to 1e-9. (3) Instance
vs Eq. (branch): regular chygraphs exact at every n (all-ones Perron
vector; 1e-9 on 4-regular and two triangles per atom); Poisson graph ⟨k⟩=4,
Poisson triangles ⟨κ⟩=2, mixed links/triangles/4-cliques (1.5, 0.7, 0.3):
n = 3200 means within 0.02 / 0.45 / 0.07 % of the ensemble, seed spread
0.7–1.8 % shrinking ~n^{−1/2}. (4) The factorisation F(u→a)→F(a)→F(a→v)
(one-dimensional stalk) holds for c ≤ 3 at *any* fixed point (one cycle
condition, satisfied because M = D·C with C symmetric) and at the trivial
point of a homogeneous complex; from c = 4 it is the condition
C12C34 = C13C24 = C14C23 and fails: residual 1.2e-2 (c=4) / 3.7e-3 (c=5) at
field spread 0.4, 8e-2 / 1.2e-1 at 1.6; ±J at zero field 0.21 / 0.39. The
identity then fails by the same order (1.4% on a small instance).

Text: Sec. 24.6 rewritten around Eqs. (mb-block), (mb-iharabass),
(mb-hessian), Fig. 24.x and Table 24.x; "What it gives the book" has a
fourth item; the field paragraph names the matrix-weighted identity as the
missing piece; Checks; software table (three rows); README.

Open: the matrix-weighted Ihara–Bass identity for blocks that do not
factorise (what stalk dimension a 4-member complex needs); the ensemble
loop series; Part IV rerun with bonds assigned once.

## 2026-09-26: the stalk dimension and the ensemble loop series (Secs. 24.6, 24.3)

**Matrix-weighted identity** (`bethehessian.py`: `symmetric_block`,
`min_rank_completion`, `stalk_dimension`, `completion_factors`,
`matrix_bethe_hessian`; 1 more test). Up to a diagonal conjugation the
block is the belief's correlation matrix with the diagonal removed. Write
M = R − diag R with R = ABᵀ of rank d; the matrix determinant lemma on the
two-step form of I − T gives det(I−T) = ∏(1−r) det [[I + P diag(r) E Pᵀ,
−PEA], [−BᵀEPᵀ, I + BᵀEA]], E = diag(1−r)⁻¹: atoms + a d-dim stalk per
complex, symmetric after conjugating the stalk block by the completion's
eigenvalue signs (indefinite metric). Verified to 1e-12 with a full-rank
completion and with minimal ones on the 4-clique instance where rank one
failed by 1.4%. Stalk dimension = minimal rank of a real symmetric diagonal
completion; count: smallest d with (c−d)(c−d+1) ≤ 2c → 1,2,3,3,4 for
c = 3..7 (Table 24.x, `tab-stalks.tex`): found numerically (multistart
BFGS + Nelder–Mead over the diagonal) at random-field points: 1,2,3,3,4
(one c = 6 draw at 4); ±J zero field: at or below the count (gauge);
homogeneous zero field: 1 always. Metric: definite for c ≤ 3, indefinite
from 4 on.

**Ensemble loop series** (`ensembleloops.py`: `ensemble_series`,
`branching_power`, `cycle_sum`, `exact_minus_bethe`; 4 tests;
`figs/ensembleloops.py`, ~50 min, cache `probe/results/ensembleloops.json`).
Only unicyclic loops survive as n → ∞; r_C = ∏ u'(c_a) over the cycle
(μ_i = 1, μ_a = ⟨σ_vσ_u⟩_a = u'); Poisson cycle counts with weighted means
tr(Bˡ)/2ℓ, B = Eq. (branch). Closed form:
E[ln Z − ln Z_BP] → −½ Σ_j (−1)^{j+1}/j [ln det(I − B_j) + tr B_j], B_j with
u'^j. j = 1 is the Gaussian correction −½ ln det(I−B), diverging at
Eq. (det) = the Bethe Hessian on the ensemble; j = 2 carries the AT matrix;
the alternating sum turns −ln(1−r) into ln(1+r) (exact on a ring). Checks:
links+triangles (1.2, 0.5; β_c = 0.4208) and links+triangles+4-cliques
(1.0, 0.4, 0.15; β_c = 0.3739): n = 300 cycle sums (400 instances, ℓ ≤ 10)
match the closed form within errors up to ~0.6 β_c (0.136±0.004 vs 0.138
at βJ = 0.243); near β_c the truncation (0.654 vs 0.799 at 0.92 β_c) and
finite size (0.588 vs 0.654) both bite. n = 18 exact vs cycle sum on the
same instances: 2% to 0.6 β_c, 10% at 0.92 β_c (cyclomatic-two loops).
First attempt used a denser second ensemble (1.5, 0.7, 0.3): cycle
enumeration to ℓ = 10 with branching ≈ 3.8 is hopeless (>1 h, killed).

Text: Sec. 24.6 "The stalk dimension" paragraph + Table (stalks), closing
paragraph and field paragraph updated; Sec. 24.3 "The series on the
ensemble" + Fig. + Table; Checks; software (3 rows); README.

Open: O(1/n) corrections (cyclomatic-two loops) on the ensemble; the
ensemble series in the ordered phase (polarised fixed point, r_C with
nonzero m); Part IV rerun with bonds assigned once; an optimised
bond-assignment rule.

## 2026-09-26: Part IV rerun with bonds assigned once (Sec. 14.3) — and the fixed point

`statmech/probe/cavity_assigned.py` (results `cavity_assigned.json`, ~100 s):
`cavity_clique.ChygraphBP` gained `assign='all'|'once'` (default unchanged,
`validate()` still exact on trees of triangles; `'once'` reproduces
`loopseries.BinaryFactorGraph.promoted(..., 'once')` on two triangles). Same
Karrer–Newman and real instances as Ch. 14 (same caches, same seeds);
hyperbolic ones are fresh draws from `figs/loopseries.hyperbolic_graph`
(the generator in `computational_complexity/code/hrg.py` is gone with that
repo, and `figs/merge.py` imports it at module level — stubbed in the probe).

**Finding 1, the fixed point.** Ch. 14's solver starts symmetric and its
damping ladder starts at 0.0, so it converges in two sweeps to the
paramagnetic fixed point on every run (|m| = 0 exactly) whether stable or
not. Stability from `bethehessian.linearised_operator` at the trivial
blocks: unstable on 222/240 (double count), 191/240 (assigned). At the
stable fixed point (polarised start, damping ladder from 0.5) the
double-count medians are 2.45/8.90 (hyperbolic, 0.3/0.8), 2.57/9.70
(Karrer), 7.70/22.1 (real) vs the paramagnetic 0.32/4.04, 0.18/1.31,
1.60/13.9 quoted in the chapter. Story survives with larger numbers. This
also explains Table 24.3's "larger than 5.06" (caption corrected) and the
two non-unique hospital runs (rounding broke the symmetry).

**Finding 2, the assignment.** At the stable fixed point, each bond in its
largest clique: median over 240 runs 7.31 → 0.35, smaller on 214/240; max
at βJ = 0.8 is 0.70 = ln 2 (the sector) on every class; 227/240 below
ln 2 + 0.05. So Ch. 14 priced double count + fixed point together; the loop
part proper is a few tenths. Fig. 14.x (fig-cavity-assigned), new paragraph
in Sec. 14.3, Checks note, Ch. 16 conclusion qualified, software row, README.

**Consequence not yet applied.** `merge_lnz.py` (Ch. 16) uses the same
symmetric ChygraphBP iteration, so its 40 cyclic Karrer runs (median 0.18)
are paramagnetic-point errors too; the 200 acyclic ones are exact at any
fixed point (unique on a forest). Rerun Ch. 16 at the stable fixed point and
against the assigned recursion (merging vs assignment: what merging buys
beyond the free repair) — next item. `gbp_*` (Ch. 15) use static counting
and GBP, not this iteration.

Open: Ch. 16 rerun as above; an optimised bond-assignment rule.

## 2026-09-26: Ch. 16 merged-family recursion at the stable fixed point

`statmech/probe/merge_stable.py`: reruns the 60 Karrer–Newman runs of
`merge_lnz.py` on the merged family with the paramagnetic point's stability
and the polarised fixed point; patches the stable-point error into
`results/merge_lnz.json` (original kept as `merge_lnz_paramagnetic.json`;
the paramagnetic values reproduce the old file to 1e-8), and sets beside it
the unmerged double-count and assigned errors from `cavity_assigned.json`
(same instances). Hyperbolic and real merged families are acyclic → unique
fixed point → unchanged.

Paramagnetic point unstable on 9/60 merged runs (all cyclic): merging moves
the instance threshold. Cyclic 40: merged median 0.182 either way, max
2.25 → 0.71 (= ln 2, the sector). Merged vs unmerged at stable points: median
0.033 vs 5.28 (was "0.289 → 0.033" at the paramagnetic point); vs the free
assignment 0.406, merged smaller on 52/60. Text: Table 16.1, Fig. 16.x caption,
Sec. 16.2 paragraphs (fixed point; merge vs assignment), conclusions range.

Open: an optimised bond-assignment rule; the two Vorob'ev checks; the rest of
the list.

## 2026-09-26: the two Vorob'ev checks (Ch. 23)

**(4) Region-graph fixed points consistent on overlaps.** The Ch. 15 caches
record `consistency` (worst parent→child marginalisation violation; the
region graph is intersection-closed, so this is agreement on every overlap).
Split by residual < 1e-9: 186/240 settled — chordal ≤ 1.3e-12, non-chordal
≤ 8.5e-8 except one real run at 0.057; the 54 unsettled reach 0.5–1.0
(damping 0.999 makes the residual meaningless there). So the non-chordal
settled fixed points are consistent families that are not the true
marginals (ln Z off by up to 12). Side finding: `mean_abs_m` median 0 on
Karrer/real GBP runs → symmetric point; exact on junction trees, stability
unchecked elsewhere (open, Ch. 15).

**Chygraph fixed point not consistent on pairs.** `probe/pair_consistency.py`
(~2 min, `results/pair_consistency.json`) on the `cavity_assigned` instances
at the stable fixed point: single atoms agree to 3e-13; shared pairs differ
by up to 0.149 (double count; > 1e-6 on 101/240) and 0.328 (assigned;
> 1e-6 on 200/240, > 0.01 on 92); all at βJ = 0.3 — at 0.8 the beliefs are
concentrated on the ordered configuration.

**(3) Definition check.** Overlaps of {a}∪∂a and {b}∪∂b: {b}∪(∂a∩∂b) if b ∈
∂a, else ∂a∩∂b. Level-wise consistency = agreement on single vertices, so
the notions coincide iff every pairwise overlap of the flattened family is
a single vertex (no two containers share two members; no container holds a
member with one of its members); otherwise level-wise is strictly weaker
(two triangles sharing an edge, one level up). Stated corollary in Sec. 23.4;
for strictly layered chygraphs the condition is Part IV's pairwise
condition at every level. Ch. 23 now has a Checks section; software row;
README.

Open: Ch. 15's GBP fixed points' stability on non-chordal instances; an
optimised bond-assignment rule; Foissy's classification (Ch. 21); Ch. 17's
list; Sec. 8.10's list; Ch. 18's list.

## 2026-09-27: Ch. 15's GBP fixed points, stability (Sec. 15.x)

`statmech/probe/gbp_stable.py` (results `gbp_stable.json`, ~70 min with the
cap; a first attempt with full Jacobians on every run was killed after 5 h
— the hospital neighbourhoods have message vectors up to 1.3e5). Symmetric
ladder reproduced from the probes (130/139 cached settled runs to 1e-6; the
9 others are Football/Hospital runs where the symmetric iteration broke
symmetry in one run and not the other — the chapter's "run twice" runs);
stability = spectral radius of the undamped sweep by finite differences,
gauge eigenvalues at 1 excluded, formed when D ≤ 3000; polarised start
beside it; reported point = symmetric if stable (or, unassessed, if the
polarised start returns), else the polarised one. Arnoldi on FD matvecs and
a belief-growth test were tried and rejected (defective gauge cluster at 1;
nearby fixed points with different beliefs).

Findings: symmetric point stable on 66/78 assessed non-chordal runs (real
24/24, Karrer 34/43, hyperbolic 8/11) — the opposite of the chygraph
recursion (222/240 unstable); on the 12 unstable the polarised point's
error is ≤ ln 2 (median 0.27). 52 non-chordal runs have no fixed point
from either start. GBP vs BP at stable points on the 84 non-chordal pairs:
medians 0.030 vs 6.49 (was 0.31 vs 2.73), GBP worse on 6/84 (was 39/142),
corr −0.45; vs the assigned recursion 0.39, GBP worse on 23/84. Text: new
paragraph "The fixed point" in Sec. 15.6, verdict numbers, figure caption,
Checks; Outlook pitfall; software row (24 caches); README.

Open: an optimised bond-assignment rule; Foissy (Ch. 21); joint
(cardinality, chy-degree) on the interactomes (Ch. 21); Ch. 17's list;
Sec. 8.10's list; Ch. 18's list; O(1/n) corrections and the ordered phase
of the ensemble series (Ch. 24).

## 2026-09-27: the bond-assignment rule (Sec. 14.3)

`statmech/probe/assignment_rules.py` (results `assignment_rules.json`,
~30 min with the oracle capped at 40 evaluations; an uncapped first run was
killed after 2.5 h on the hospital neighbourhoods). `BinaryFactorGraph.promoted`
and `ChygraphBP` accept an explicit `{bond: clique}` assignment (test: the
exact model is unchanged; a first version of the hook dropped the unshared
bonds, which made the stability test pass wrongly — caught by the mismatch
with `cavity_assigned` on the runs with an unstable symmetric point).
Protocol = Sec. 14.3's stable point (the uniform-start BP of Sec. 24.3 sits
at the unstable symmetric point on these instances — its errors are the
paramagnetic ones). Medians over 240 runs: largest 0.343 (= cavity_assigned's
0.35), strongest indirect correlation 0.354 (picks the largest clique almost
always), smallest 0.438, weakest 0.454, balance 0.465, random median 0.440,
random best-of-10 0.285, oracle 0.190 (agrees with largest on 71% of shared
bonds; at βJ = 0.3 the oracle reaches 0.06, the sign of the error being
tunable; at 0.8 the ln 2 floor holds it at 0.31). Verdict: keep the
largest-clique rule; the assignment matters at the factor-of-two level and
nothing structural captures it. Text in Sec. 14.3 and Sec. 24.9; software
row (25 caches); README.

Open: Foissy (Ch. 21); joint (cardinality, chy-degree) on the interactomes
(Ch. 21); Ch. 17's list; Sec. 8.10's list; Ch. 18's list; O(1/n) corrections
and the ordered phase of the ensemble series (Ch. 24).
