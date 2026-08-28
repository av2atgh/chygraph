# Complex hypergraphs (book)

Book-length treatment of percolation and statistical mechanics on complex
hypergraphs. Production style follows *Local network growth* — 5×8 in trim,
grayscale, skippable "calculation" boxes, heavy on TikZ illustration, plain
language in the running text with the algebra boxed off.

## Contents

| | Source | Topic |
|---|--------|-------|
| — | `preface.tex` | **Drafted.** The tree assumption, the bargain, what the book claims |
| — | `software.tex` | **Drafted.** Repo links; equation-to-method-to-test table; how to reproduce a figure |
| **I** | | **Foundations** |
| 1 | `introduction.tex` | **Drafted.** Percolation; Ising, vertex cover, hitting set; mean field, cavity, replicas, BP; why real networks break treelikeness |
| 2 | `chygraphs.tex` | **Drafted.** The object: complexes whose vertices are complexes; every higher-order structure as one thing |
| 3 | `data.tex` | **Drafted.** Where chygraphs come from: papers, protein complexes, reactions, schedules; cliques of a clustered graph; how much overlap there is |
| **II** | | **Percolation** |
| 4 | `percolation.tex` | **Drafted.** The self-consistency map and the threshold tensor, structure by structure |
| 5 | `giant.tex` | **Drafted.** Order parameter, critical amplitude, the moment hierarchy, dependent layers, six constructions solved |
| 6 | `epidemics.tex` | **Drafted.** SIR as percolation, two levels of mixing, the household reproduction number, contagion inside a group |
| **III** | | **Statistical mechanics** |
| 7 | `potts.tex` | **Drafted.** Fortuin–Kasteleyn: percolation and Ising are one recursion at two values of $q$ |
| 8 | `statmech.tex` | **Drafted.** General theory: convolution up, exact interior sum down, the branching matrix |
| 9 | `ising.tex` | **Drafted.** Ising on chygraphs; clustering lowers $T_c$; the AT line; the unanimity interaction |
| 10 | `hittingset.tex` | **Drafted.** Hard fields, where they fail off the graph, soft fields, RSB |
| 11 | `cover.tex` | **Drafted.** Vertex cover, leaf removal, core percolation; hyperbolic random graphs |
| 12 | `colouring.tex` | **Drafted.** Proper against hypergraph colouring; `tau = -1/(q-1)` independent of cardinality; `(q-1)^2` |
| 13 | `satisfiability.tex` | **Drafted.** Clauses as complexes; `alpha = 1` exact at `k = 2`, no linear instability above it; **clauses whose members are clauses**, and what CNF flattening costs |
| **IV** | | **Limits** |
| 14 | `overlap.tex` | **Drafted.** The price of treelikeness; region graphs and GBP; what is open |
| 15 | `outlook.tex` | **Drafted.** One recursion, many models; the two running threads; what is not done |

**Chapters 12 and 13 have no manuscript behind them.** Every other chapter is
exposition of work that exists elsewhere; for these two the calculations are
done in the book's own figure scripts, `figs/colouring.py` and
`figs/satisfiability.py`, and checked against the published thresholds. PDFs of
the references are under `~/Downloads/chygraph_references/`.

## Status

Last updated 2026-08-28. `main.pdf` builds clean: **196 pages, 0 errors, 0
undefined references, 0 overfull boxes.**

**The book is drafted end to end.** Preface, the software chapter, and fifteen
numbered chapters, all with prose, figures and checks.

`software.tex` sits after the preface and is unnumbered, like the preface. It
carries the two repository URLs and **Table 1: every computed equation, the
routine that evaluates it, and the test or script that checks it** — the book's
version of the retired supplement's Sec. I. Every one of its 58 equation labels and
49 code names was verified to resolve before the table was written; if a routine
is renamed, that table is where it has to be fixed.

### What to do next

The remaining work is revision, not drafting:

| | what |
|---|---|
| 1 | **A read-through for continuity.** Nobody has yet read 1–15 in order. Watch for the two running threads being restated more often than they need to be, and for the "cardinality two is degenerate" observation, which now appears in Chs. 10, 11, 12 and 13 and should be stated once properly and referred back to. |
| 2 | ~~Fix the manuscripts~~ — **done by deletion.** The three errors the book found (the tricritical cardinality, the GBP convergence count, the Ch. 5 arithmetic slip) were in `main.tex`/`supplement.tex`, which no longer exist. The book carries the corrected values and says so; the entries below are kept as a record of what was wrong, in case either file is ever resurrected from git. |
| 3 | **Chs. 12 and 13 are the only chapters whose results are not backed by a manuscript.** If either is to be published separately, the calculations in `figs/colouring.py` and `figs/satisfiability.py` are the starting point. |
| 4 | Front matter: the dedication is still `\itshape Dedication to come.` in `main.tex`. |
| 5 | **Make the repositories public.** `software.tex` prints both URLs as if they resolve; `chygraph_statmech` is still private. This is `../TODO.md` item 3 and now blocks the book as well as the manuscript. |

Each chapter's `\include` in `main.tex` carries a comment naming its sources, so
that mapping does not have to be reconstructed.

### Decisions already taken — do not relitigate

- **Title is `Phase transitions on complex hypergraphs`.** No subtitle.
  (Changed 2026-08-28 from `Complex hypergraphs`.)
- **Contents lists parts and chapters only.** `\setcounter{tocdepth}{0}` sits
  immediately before `\tableofcontents` — it must go there, not in the preamble,
  where the packages override it. (The `.toc` *file* always contains every
  section; the filtering happens when the contents is typeset, so do not
  "fix" this by reading `main.toc`.)
- **Four `\part` divisions.** Not in *Local network growth*, added because
  thirteen chapters need them.
- **Chapter order.** Ch. 3 (data) sits before any machinery so a reader can
  build a chygraph before being asked to solve one. Ch. 6 (epidemics) follows
  percolation because SIR *is* percolation. Ch. 7 (Potts) is the hinge that
  makes Parts II and III one subject.

### Material in the book that is not in any manuscript

- **Sec. 4.5, "Do triangles help or hurt?"** resolves an apparent conflict
  between the published percolation result (triangles help, at fixed *link
  count*, `dtheta = 2q^2(1-q)<k>_tri > 0`) and Ch. 9's Ising result (clustering
  lowers `T_c`, at fixed *degree*). The link-matched comparison lets the degree
  distribution broaden; held at fixed degree the percolation sign reverses too,
  `q_c` going from `1/3` to `0.4030`. Verified three ways in
  `figs/triangles_vs_regular.py`.
- **Ch. 3's measurement on ten real networks** (Table 3.2, Fig. 3.2). Not in any
  manuscript. Its third finding — that a degree-matched control can carry *more*
  clique overlap than the real network, because heavy tails manufacture cliques
  at hubs — is a caution that applies to every degree-matched comparison in the
  book, Ch. 11 included.
- **Sec. 6.4, "Contagion inside a group."** The percolation machinery cannot do
  a threshold rule — a group that transmits only when several members are
  infected — and the chapter says so rather than pretending otherwise. It gives
  the reason (a product over onward routes assumes one number per inclusion),
  the consequence (with `r >= 2` the trivial fixed point stays stable, so
  `Lambda = 0` does not locate the transition — the same exclusion Sec. 5.4
  flags for interdependent networks), and then a self-consistent mean-field
  equation, `rho = 1 - exp(-b1 rho - b2 rho^2)`, with an exact tricritical point
  at `b2 = 1/2`. This does **not** contradict Sec. 5.4's proof that group
  structure alone gives continuous transitions: a threshold rule is outside the
  class that proof is about, and the chapter says so in as many words. Keep the
  two statements aligned if either is edited.
- **Chapter 7 in its entirety.** The manuscripts carry the Fortuin–Kasteleyn
  correspondence as one paragraph of the retired manuscript's Sec. I. The chapter turns it
  into a proved statement about chygraphs: `lim_{q->1}` of the Potts interior
  sum on a complex **is** Ch. 4's `Gbar`, and at `q = 2` it is Ch. 9's Ising
  cavity field. Both directions are verified in `figs/potts.py` by summing the
  interior at *symbolic* `q` over set partitions (a partition into `m` blocks
  admits a falling factorial `q(q-1)...(q-m+1)` of colourings, which is a
  polynomial in `q` and continues to real `q` on its own). The by-product is the
  single transmission factor `tau_c(q,v)`: `tau_2 = v/(q+v)`, which is the bond
  probability `p` at `q = 1` and `tanh(beta J)` at `q = 2`; `tau_3` is
  `sbar_tri/2` at `q = 1` and `t/(1-t+t^2)` at `q = 2`. So Ch. 4's triangle
  enumeration and Ch. 9's triangle transmission are one formula at two points —
  Ch. 9 is now written and Eqs. (7.7), (8.5) and (9.2) do agree; keep them
  consistent if any is edited. The
  chapter also records that above `q = 2` the transition is first order, so the
  linear condition is a spinodal (numbers from iterating the recursion, in the
  script).
- **Chapter 12 in its entirety.** No manuscript; the calculations are in
  `figs/colouring.py` and are done twice, by exact enumeration of a complex's
  interior in rational arithmetic and from a closed form. Two results worth
  guarding. **Proper colouring: `tau = -1/(q-1)`, independent of cardinality** —
  a clique of nine transmits what an edge transmits, because a traceless
  perturbation carries one fact ("not that colour") however many members repeat
  it. **Hypergraph colouring: `tau = -1/(q^(c-1) - 1)`**, which collapses like
  `q^-c`. So in neighbours per node the stability threshold is `(q-1)^2` at
  every cardinality for the proper rule and `(q^(c-1)-1)^2` for the hypergraph
  one — identical at `c = 2`, three orders apart by `c = 4`. The `c = 2` case
  reproduces **both** lines of Zdeborová & Krzakala Eq. (18), `(q-1)^2` for
  Poisson and `(q-1)^2 + 1` for regular chy-degree, which is a strong check
  because a formalism with the excess-degree bookkeeping wrong would get at most
  one of them. **Sec. 12.5 must not be softened:** what is computed is a local
  stability threshold, and for `q >= 4` it sits *above* the colourability
  threshold, not inside the colourable phase.
- **Chapter 13 in its entirety.** No manuscript; calculations in
  `figs/satisfiability.py`. Three results, of very different standing, and the
  chapter is careful about which is which. **Exact and convention-free:** the
  warning map's derivative at the trivial fixed point is `alpha` for `k = 2`,
  giving the rigorous 2-SAT threshold `alpha = 1`; and for `k >= 3` it is
  *exactly zero at every clause density*, so there is no linear instability
  anywhere — the extreme form of the book's "a linear instability is not always
  the transition" thread. **Convention-dependent:** the fold at which a
  non-trivial warning branch appears (`k = 3`: 5.386 strict, 4.667 net-field).
  Under *either* convention it lies above `alpha_s`, so the direction is robust
  and the number is not; the chapter quotes both. The structural point to keep:
  Eq. (13.2) is Ch. 10's hitting-set map with a factor 1/2 (negatable literals)
  and the opposite monotonicity — order-preserving, so no `F o F` bracket is
  needed. Also worth guarding: SAT has **no** uniform fixed point, because a
  clause forbids a *particular* assignment; that asymmetry is the whole reason
  Ch. 12 gets a closed form and Ch. 13 does not.
- **Secs. 13.6–13.7, clauses of clauses.** The only place in the book where a
  complex contains another complex, which is the freedom the chygraph
  definition exists for. Layer 2 clauses are implications `C => (l_1 v ... v
  l_k2)` whose antecedent `C` is a layer-1 clause; the layer-1 complex
  compresses its `2^k1` interior states into one number and the layer-2 complex
  consumes it exactly as it would a variable. Three results, all in
  `figs/satisfiability.py`:
  **(i)** at `k1 = 1` the middle layer is a relay and the system reduces to
  flat `(k2+1)`-SAT *identically, to the last digit* — the check that the
  degree bookkeeping is right, and the thing to re-run first if these equations
  are ever edited;
  **(ii)** distributing `NOT C` gives `k1` clauses that pairwise share all `k2`
  plain literals, so for `k2 >= 2` the flattening manufactures exactly Sec. 14.1's
  non-treelike condition — the overlap is created by the *encoding*, not found
  in the data;
  **(iii)** thresholds agree *exactly* at `k2 = 1` (both give `k1 alpha = 1`,
  as an identity, not numerically) and the folds differ by 7–22% at `k2 >= 2`,
  monotonically in `k1`. The agreement at `k2 = 1` is what makes (iii) evidence
  for (ii) rather than a coincidence: it is the rewired-control logic of
  Ch. 11 applied to an encoding instead of an ensemble.
- **Sec. 14.6, merging overlapping complexes.** Not in any manuscript. The
  obvious repair for Fig. 14.1 — merge two complexes that share 2+ atoms into
  one meta-complex and enumerate its interior — is *correct*, and legal because
  a complex may contain complexes. What decides it is cost, and **the merge
  closure is itself a percolation problem** on the graph whose nodes are
  complexes and whose links are shared pairs: affordable exactly while that
  graph has no giant component. So the feasibility of a Part IV repair is a
  Part II calculation. Measured on HRGs at `tau = 2.5`, the largest
  meta-complex is bounded (~5 atoms, flat in `n`) only below `<k> ~ 0.4`, and
  reaches 28% of all vertices at `<k> = 3`. Two details to keep if this is
  edited: the closure must be **iterated** (two meta-complexes from disjoint
  components can still share 2 atoms between their unions), and a merged
  complex inherits every inclusion of its parts, which is what makes the
  component grow.
- **Sec. 14.5's convergence count differed from the retired supplement.** It
  said 9 of the 20 non-chordal GBP runs converge and 11 do not, with
  errors `1.4e-9`–`5.4e-3` against Möbius `0.48`–`1.1` and Bethe `0.32`–`22`.
  Recomputed from the cached `../probe/results/gbp_cliques.json`, the residuals
  fall into two clusters separated by four orders of magnitude — everything at
  or below `3.2e-6`, then nothing until `2.0e-2` — so the threshold sits in the
  gap at `1e-4` and gives **14 converge, 6 do not**, with Möbius `0.48`–`1.8`
  and Bethe `0.26`–`30`. The GBP error range is unchanged. The book states the
  recomputed counts and says what criterion produced them.
  The supplement's numbers looked like an earlier run of the probe. It has since
  been deleted, so the book's values stand; this entry is the record.
- **Ch. 11 pairs two counterexamples that the retired manuscript kept apart.** Sec. 10.3
  (disjoint 3-hyperedges, hitting set, truth 1/3) and Sec. 11.2 (isolated
  triangles, induced-graph cover, truth 2/3) are the same hard-field defect
  producing the *same* wrong answer 1/2, because the half rule always returns
  1/2 when every vertex is degenerate. The book states the pairing; the
  manuscript has the two facts in different subsections without connecting
  them.
- **Sec. 11.5's `1 - e^{-<k>}`.** Above cardinality two the core of a pure
  clique network is exactly the non-isolated fraction, and it does not depend on
  cardinality at all — c = 3 and c = 4 give identical curves. That sharpening of
  "no core-free branch" is not in the manuscript and is what Figure 11.2 draws.
- **Sec. 10.3's counterexample is the book's only proof.** Disjoint
  3-hyperedges: no complex touches another, so the cavity equations are exact
  and replica symmetry cannot break, the density is 1/3 by counting, and the
  hard-field limit returns 1/2. It was one sentence of the manuscript; the book
  makes it the hinge of the chapter, because it settles the question by
  arithmetic rather than by numerics and its bracket `(0,1)` shows the
  hard-field scheme cannot self-diagnose. The mechanism to keep straight if this
  is edited: the integer message cannot separate the taken vertex from the two
  left out, the distinction sitting in the `O(1)` part that was scaled away —
  and Eq. (10.6)'s `h_RS = -mu/L - ((L-1)/L) ln(K-1)` carries a *fraction* of
  `mu`, which is precisely what an integer ansatz cannot hold.
- **Sec. 9.6 corrects the retired manuscript on the tricritical cardinality.**
  It said the continuous/discontinuous boundary for the unanimity
  interaction "sits between five and six only at k = 3 and k = 4 and has moved
  to between four and five by k = 6". Recomputed with a criterion that demands a
  coexistence window of non-negligible width **and** a free-energy crossing
  strictly inside it, the boundary is between q = 5 and 6 **only at k = 3**, and
  between q = 4 and 5 from k = 4 upward (checked at k = 4, 5, 6, 10, 40; at
  q = 5 the window runs 1.5% to 18% of `T*` with `dm` from 0.68 to 0.86, and at
  q = 4 the apparent window is ~5e-4 of `T*` with no crossing at any k). The
  manuscript's own warning is what explains the discrepancy: just above a
  continuous transition the ordered branch decays slowly enough that a finite
  iteration budget reports coexistence. The book states the corrected boundary
  and the reason; the manuscript has since been deleted. Probe:
  `figs/ising.py:check_tricritical`.
- **Ch. 8's cross-chapter identification.** `u'` of Eq. (8.5) and Ch. 7's
  `tau_c` are the same number; `figs/statmech.py` checks it for `c = 2, 3, 4` by
  taking Ch. 7's symbolic `q`-derivative to `q = 2` and comparing with
  `chygraph_statmech.cavity`'s Ising enumeration. This is what makes Eq. (7.9)
  the single-layer case of `det(I - B) = 0`, and it means Eqs. (7.7), (8.5) and
  whatever Ch. 9 quotes for `u'` must all agree: `t`, `t/(1-t+t^2)`,
  `(t+t^3)/(1-2t+3t^2)`.
- **Sec. 6.5's fixed-budget comparison.** Households help the epidemic when
  added on top of an unchanged global network (`T_c` 0.5 -> 0.1667) and *hurt*
  it at a fixed contact budget of four contacts per person (`T_c` 0.2500 ->
  0.2929, and 5.6% vs 31.4% infected at `T = 0.3`). Same running thread as
  Ch. 3, Sec. 4.5 and Sec. 5.5.

### The two running threads

Both are now gathered in **Sec. 15.2 of the Outlook**; this is the index of where
each instance lives, so that editing one does not leave the others behind.

Five chapters arrive at the same lesson from different directions: **name
what is held fixed before reading a comparison as a statement about
clustering.** Ch. 3 (a degree-matched control manufactures clustering), Ch. 4
(a link-matched control removes degree heterogeneity), Ch. 5 (identical
marginals, three different transitions), Ch. 6 (households help at a fixed
global network and hurt at a fixed contact budget), Ch. 9 (a Poisson link layer
matched only on mean degree reverses the sign of the clustering effect, because
the triangle construction carries `<kbar> = n + 1`, not `n`).

The second thread appeared four times: **a linear instability is not always the
transition.** Sec. 5.4
(interdependent networks appear by saddle-node while `Q = 1` is still stable),
Sec. 6.4 (threshold contagion, same exclusion, no invasion from one seed),
Sec. 7.4 (the ferromagnetic Potts transition is first order above `q = 2`, so
the branching condition is a spinodal — 5.7% to 25.5% off at `q = 3` to `6`),
and Sec. 9.6 (the unanimity interaction, where the free energy of Sec. 8.6 is
finally *used*: the q = 16 double transition has its first-order point at
`T_c = 0.9790` inside a window whose lower edge is the spinodal). Sec. 9.6 also
carries the practical corollary — critical slowing down mimics a coexistence
window, so the test needs a free-energy crossing, not a surviving branch. Keep
the four consistent; they are one statement.


## Underlying material

- `~/Dropbox/submissions/hyperabs.2022/hyperabs_v3.tex` — Phys. Rev. E **107**, 024316 (2023)
- `~/Dropbox/submissions/chygraph.2023/chygraph.tex` — J. Complex Netw. (2024), cnae047
- `~/av2atg/chygraph/manuscript_3/manuscript.tex` — the percolation extension
- the manuscript this book supersedes (`../main.tex`, `../supplement.tex`) —
  **deleted**, recoverable from git at commit `5ebd892`, the last commit before
  the deletion. Everything in it is in the book; the section numbers named in
  `main.tex`'s `\include` comments refer to that retired file.
- `../src/chygraph_statmech/`, `../examples/`, `../probe/` — every number in the book comes out of these

## Building

```sh
latexmk -pdf -interaction=nonstopmode main.tex     # -> main.pdf
```

Verify a build with

```sh
grep -c '^!' main.log                 # errors, must be 0
grep -i 'undefined' main.log          # citations and refs, must be empty
grep 'Overfull' main.log              # must be empty
grep -o 'Output written.*' main.log   # page count
```

Two environment notes. `latexmk` will hang rather than stop on a TikZ error even
under `-interaction=nonstopmode`; add `-halt-on-error`, or run `pdflatex`
directly, when a figure is new. GNU `timeout` is not on the path here.

To inspect a page at trim size, `pdftoppm -png -r 120 -f N -l M main.pdf out`
and read the PNG. **Always do this for a new figure** — compiling is not
evidence that a figure is legible at 5×8.

### Figure scripts

All live in `figs/` and write PDFs and `.tex` tables into that directory, as in
*Local network growth*. Rerun a script after editing it; the build picks up the
new PDF automatically.

`figs/real_chygraphs.py` generates
Table 3.2 and Figure 3.2 from the netzschleuder dumps cached under
`~/av2atg/LocalNetworkGrowth/figs/data/netzschleuder`;
`figs/triangles_vs_regular.py` generates Figure 4.2 and verifies
$\bar s_\triangle(q)$ by exact enumeration;
`figs/beyond_threshold.py` generates Figures 5.1–5.3 by importing the
`chygraph` package from `~/av2atg/chygraph/src`;
`figs/epidemics.py` generates Figures 6.2 and 6.3 and runs every check quoted in
Ch. 6 — that `Gbar^1_0` at `n = 3` is the bond-percolated triangle of Ch. 4,
that `theta + 1 = R*` for several household size distributions, the
fixed-budget comparison, and the tricritical condition and exponent of
Sec. 6.4. Figure 6.1 is TikZ in `epidemics.tex`.
`figs/potts.py` generates Figure 7.2 and runs Ch. 7's checks: the `q -> 1` limit
of the Potts interior sum against `Gbar` (cliques to `K_5`, paths, distinct `y_j`
per member, symbolic in `p`), the `q = 2` limit against the Ising cavity field of
`chygraph_statmech.cavity`, the transmission factors and their two limits, the
linear threshold, and the order of the transition. It takes a couple of minutes —
the `K_5` interior sum is the slow part. Figure 7.1 is TikZ in `potts.tex`.
`figs/statmech.py` generates Figure 8.2 and runs Ch. 8's checks: that `u'` is
Ch. 7's `tau_c` at `q = 2`; the size-biased average three ways (mixed layer,
split into one layer per cardinality — identical — and the mean-cardinality
substitution, wrong by 10.1%); the links-and-triangles determinant for Poisson
and regular layers, so the cross term is exercised both where it vanishes and
where it does not; the Bethe free energy against its closed form and the
textbook graph result; and the exchange of stability, taken along a sequence in
`Lambda` and stopped before rounding dominates. It imports `figs/potts.py`, so
the two live together. Figure 8.1 is TikZ in `statmech.tex`.
`figs/ising.py` generates Figures 9.1–9.3 and runs Ch. 9's checks: the closed
forms for `u'` at c = 2, 3, 4 against the enumeration; the fixed-degree
clustering family (monotone in f at every degree) and the regular-null table;
the AT line; the unanimity `u'`, the Bragg–Williams limit and its `C_q`
residual including the vanishing at q = 4; and the tricritical boundary. The
Monte Carlo points in Figure 9.2 are read from the cached
`../probe/results/ising_mc.log` rather than recomputed — that run is the one
test in the chapter from outside the formalism. `figure_unanimity` and
`check_tricritical` each solve coexistence windows and take a minute or two.
`figs/hittingset.py` generates Figure 10.2 and runs Ch. 10's checks: the
hard-field threshold `<k>(c-1) = e` to eleven digits for c = 2..20; the
disjoint-3-hyperedge counterexample (hard 0.5, soft 0.3333, truth 1/3);
Weigt–Hartmann at cardinality two; the regular closed forms `h_RS`, `rho = 1/K`
and the entropy; the mixed-cardinality and correlation scans with the
`isolated_fraction` confound reported alongside; and the entropy criterion.
Figure 10.1 is TikZ in `hittingset.tex`. Soft-field densities come from
population dynamics and carry a scatter of a few parts in a thousand — the
script averages over seeds and prints the spread, and the book quotes those
numbers to three decimals only.
`figs/cover.py` generates Figures 11.2 and 11.3 and runs Ch. 11's checks: the
`<k> = e` identification from core percolation and from Ch. 10's hard-field map
(agreeing to 2.5e-10); the induced-graph cover against Weigt–Hartmann at
cardinality two and its 1/2-vs-2/3 failure on isolated triangles; the
factorisation `f(delta) = -(1-delta)^2 sum (j+1) delta^j` symbolically for
c = 2..9; that the core is exactly `1 - e^{-<k>}` above cardinality two; and the
cached leaf-removal validation. Figure 11.1 is TikZ in `cover.tex`. The
hyperbolic-graph points are read from `../probe/results/prediction4.csv` and the
simulation check from `../probe/results/validate_core.txt`; neither is
recomputed by the book.
`figs/colouring.py` generates Figure 12.1 and runs Ch. 12's checks: both closed
forms against exact enumeration in rational arithmetic; that the linearised map
is a scalar on the traceless subspace (three directions, one answer); the graph
case against Zdeborová & Krzakala Eq. (18); and the comparison with Mulet's
published thresholds that shows the stability line is not the colourability
line. It is fast — everything is exact rational arithmetic over small
complexes.
`figs/satisfiability.py` generates Figure 13.1 and runs Ch. 13's checks: the
clause interior against enumeration in rational arithmetic (including the
non-uniform emitted message that rules out a symmetric fixed point); the `k = 2`
linearisation and the branch appearing continuously at `alpha = 1`; that the
`k >= 3` derivative is *exactly* zero, tested by the rate at which the finite
difference falls (`eps^(k-2)`, checked by halving `eps`) rather than by an
absolute bound, which would only show `eps` was small; and both conventions for
a contradicted variable.
`figs/overlap.py` generates Figure 14.2 and runs Ch. 14's checks: the Möbius and
Bethe counting numbers on two overlapping triangles, including the factor-coverage
test that is the whole point (shared bond covered once against twice); the
two-triangle table, recomputed here since it is a four-spin enumeration; the
60-instance summary from `../probe/results/gbp_cliques.json` with the
convergence threshold stated; and the clique-ensemble paired ratio read from
`../probe/results/analysis.txt`. Figure 14.1 is TikZ in `overlap.tex`.
`figs/merge.py` generates Figure 14.3 and measures Sec. 14.6: the merge closure
on hyperbolic random graphs. It generates the graphs and enumerates their
maximal cliques, so it is the slowest script in the book — a few minutes at
`n = 8000`.

## Conventions

- Drawing vocabulary is fixed once, in `main.tex`: `nd` (node), `ndf` (filled
  node), `hub`, `cx`/`cxb` (a complex drawn around a set of nodes), `lnk`,
  `msg` (a message along a directed inclusion), `lb`/`tb` (labels), and `\cxg{c}`
  for an inline glyph of a cardinality-$c$ complex. Every figure is built from
  these, so a reader who learns the first schematic can read the last one
  without its caption.
- Notation matches the papers (`\ave{}`, `\kbar`, `\sbar`), so a reader can move
  between book and paper without retranslating.
- Long derivations go in `calculation` boxes and must be genuinely skippable:
  the running text has to carry the argument on its own.
- Overlapping complexes are drawn as two translucent boxes of the same grey, so
  the shared atoms show as a darker patch. Established in Fig. 2.7 and used
  wherever overlap matters.
- Figures must be checked at the 5×8 trim, not just for compilation: more than
  three panels side by side is illegible on this page width. Stack into rows.
