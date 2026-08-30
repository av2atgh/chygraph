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
| 12 | `colouring.tex` | **Drafted.** Proper against hypergraph colouring; `tau = -1/(q-1)` independent of cardinality; `(q-1)^2`; a graph with triangles, where the graph calculation is wrong by one; survey propagation, and what the rule costs above cardinality two |
| 13 | `satisfiability.tex` | **Drafted.** Clauses as complexes; `alpha = 1` exact at `k = 2`, no linear instability above it; **clauses whose members are clauses**, and what CNF flattening costs |
| **IV** | | **Limits** |
| 14 | `overlap.tex` | **Drafted.** The price of treelikeness; region graphs and GBP; what is open |
| 15 | `outlook.tex` | **Drafted.** One recursion, many models; the two running threads; what is not done |

**Chapters 12 and 13 have no manuscript behind them.** Every other chapter is
exposition of work that exists elsewhere; for these two the calculations are
done in the book's own figure scripts, `figs/colouring.py` and
`figs/satisfiability.py`, and checked against the published thresholds --- which
since 2026-08-30 the last section of each chapter **computes** rather than
quotes, by survey propagation at `m = 0`. PDFs of
the references are under `~/Downloads/chygraph_references/`.

## Status

Last updated 2026-08-30. `main.pdf` builds clean: **228 pages, 0 errors, 0
undefined references, 0 multiply-defined labels, 0 overfull boxes, 0 underfull
vboxes.** 41 figures, 19 numbered tables, 23 calculation boxes, 104 numbered
equations, 70 references, an 85-term index.

**The book is drafted end to end.** Preface, the software chapter, and fifteen
numbered chapters, all with prose, figures and checks.

`software.tex` sits after the preface and is unnumbered, like the preface. It
carries the two repository URLs and **Table 1: every computed equation, the
routine that evaluates it, and the test or script that checks it** — the book's
version of the retired supplement's Sec. I. Every one of its equation labels and
all 92 of its `\code` names were re-verified to resolve on 2026-08-28 (the code
names across *both* repositories); if a routine is renamed, that table is where
it has to be fixed.

### What to do next

The remaining work is revision, not drafting:

| | what |
|---|---|
| 1 | ~~A read-through for continuity~~ — mechanical consistency passes were done 2026-08-28 and 2026-08-30, and a language pass 2026-08-28 (see below). What none of them could do is judge the argument: nobody has yet *read* 1–15 end to end for whether it persuades. |
| 2 | ~~Fix the manuscripts~~ — **done by deletion.** The three errors the book found (the tricritical cardinality, the GBP convergence count, the Ch. 5 arithmetic slip) were in `main.tex`/`supplement.tex`, which no longer exist. The book carries the corrected values and says so; the entries below are kept as a record of what was wrong, in case either file is ever resurrected from git. |
| 3 | **Chs. 12 and 13 are the only chapters whose results are not backed by a manuscript.** If either is to be published separately, the calculations in `figs/colouring.py` and `figs/satisfiability.py` are the starting point. |
| 4 | **Front matter: the dedication is still `\itshape Dedication to come.` in `main.tex`, and it prints.** The editorial review calls this a blocker. It needs the author's words; nothing else about the front matter is now outstanding. |
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
  linear condition is a spinodal. The chapter gives all three points of the
  coexistence window: the disordered spinodal `v_c = q`, the fold (ordered
  spinodal) at exactly `2 sqrt(q-1)` — **Eq. (7.11)** since 2026-08-30, proved by
  eliminating the message between the fixed-point equation and its derivative,
  which leaves `v^2 - 4q + 4 = 0` — and the transition itself, from a Bethe
  free-energy crossing between them (`figs/potts.py:check_transition_point`).
  `v_c` overstates the **transition** by 5.1 / 11.9 / 22.8% at `q = 3 / 4 / 6`;
  the window is wider, 5.7 / 13.4 / 25.5%. Quote the first set --- the second
  bounds the error rather than measuring it.
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

### Technical-correctness pass, 2026-08-28

A chapter-by-chapter check of the mathematics, not the prose: every closed form
re-derived independently where that was possible, every quoted number re-run
against the scripts. Chapters 1, 4, 5, 6 came through clean. What the rest
needed:

- **Ch. 2** contradicted itself five lines apart on `s` vs `sbar` — line 393
  says cardinality and component size coincide on ubergraphs (right), line 398
  told the hypergraph reader to read *s* as "cardinality minus the one I arrived
  by" (that is `sbar`). Fig. 2.5's caption miscounted reachable members the same
  way. Both fixed.
- **Ch. 3** attributed Friedrich & Krohmer's `n^((3-tau)/2)` to *the number of
  maximal cliques*. It is the **clique number**, which is what `probe/analyze.py`
  measures (`c_max`) and where the quoted 0.37/0.20/0.07 come from. The number of
  maximal cliques grows like `n^1.0` (881 -> 280848 over the range), so the
  sentence as written claimed ~12 maximal cliques at n = 3e5 — and it is the
  chapter's only independent check that the measurement reads the right thing.
  Fixed. Also "excess degree `<k^2>/<k>`" -> "degree heterogeneity" (Eq. 1.1
  defines the excess degree as `<k^2>/<k> - 1`).
- **Ch. 8, Eq. (8.6)** put the `u'` reweighting on the **wrong channel**:
  `<kappa> -> <kappa><u'>`, contradicting its own next sentence, Eq. (8.7), and
  the code (`stability.py`: "wkappa = 1 and ws carries the physics"; every
  constructor in `models.py` passes `ws=`). Now `<s> -> <s u'>`,
  `<sbar> -> <sbar u'>`. **The same slip was in `stability.py`'s docstring
  summary** — fixed at the source too. Also, "the whole difference between the
  two natural null models is one term that is easy to drop" was wrong: the
  Poisson/regular gap is mostly `<kbar>`, and deleting the cross term from the
  regular case moves T_c 6.8801 -> 5.2951 (the other way). Text now says so.
- **Ch. 9, Sec. 9.4** called the triangle construction's `n+1` a **chy-degree**
  excess. `n+1` is the *graph* excess degree `<kbar>`; the triangle layer's
  `<kappabar>` is `n/2`, and that is what Eq. (8.7) uses. Rewritten to give the
  real mechanism (a Poisson layer is its own excess, so
  `(c-1)<kappabar> = 2(n/2) = n` matches the link layer exactly and the lost
  branch is not lost), with the three T_c values added and a new
  `figs/ising.py:check_misleading_comparison` behind them. The Outlook repeated
  the same symbol error; fixed there too.
- **Ch. 10** claimed the negative-correlation branch is "clean everywhere" in a
  paragraph that had just said the isolated fraction must be reported. It runs
  0.033, 0.047, 0.198, **0.321**. Fixed. And "56% for c = 4 alone, under 1% once
  ordinary edges dominate" named no ensemble (the c = 4 gap runs 37–74% over the
  plotted range); traced to the retired manuscript's table — c = 4 Poisson at
  `<k> = 1` and the {2,6} mixture at `{1.5, 0.5}` — re-measured at 55.0% and
  0.77%, and the ensembles are now named.
- **Ch. 11** listed the measured HRG exponents as "rises monotonically, 1.59,
  1.64, 1.49" in tau order 2.5, 2.9, 2.1 — not monotone as printed. Both lists
  now run in ascending tau.
- **Ch. 12** had a botched edit in "What would be needed" (doubled em-dash,
  "survey propagation" twice). Fixed. **Both of its central results were
  re-derived from scratch** in exact rational arithmetic and are correct.
- **Ch. 13**: typo "and it little had to be added". The three-level warning map,
  the k1 = 1 reduction and the k2 = 1 identity were all re-derived by hand and
  are right.
- **Ch. 14**: "Couple every pair inside each triangle" would, under Eq. (9.1)'s
  convention, double the shared bond and give lnZ = 2.9607 at bJ = 0.2 rather
  than the table's 2.8887. The code deduplicates (`gbp.py:clique_edges`); the
  text now says the shared bond carries one coupling.
- **Ch. 15** omitted satisfiability from the "a linear instability is not the
  transition" thread (Ch. 13 enumerates five), and said the threads "turned up
  in three or four independent places" where the chapter itself says five twice.
  Both fixed, and the README's own count reconciled.

**Sec. 7.4 now quotes the transition, not just the spinodals.** Added
`figs/potts.py:check_transition_point`, a Bethe free-energy crossing on the
d-regular graph; at q = 2 all three points coincide, which is the check that the
free energy is right. `v_c` overstates the **transition** by 5.1/11.9/22.8% at
q = 3/4/6; the coexistence window is wider (5.7/13.4/25.5%). The ordered
spinodal is exactly `2 sqrt(q-1)`.

Verified and unchanged: all `\ref`/`\eqref` resolve; all 92 `\code` names in
`software.tex` resolve across both repositories; the numbers repeated across
chapters agree.

### Consistency pass, 2026-08-28

A mechanical check of the whole book. What it found and fixed:

- **Ch. 13 miscounted the "a linear instability is not the transition" thread** —
  said four, listing interdependent networks, threshold contagion, Potts above
  `q = 2` and satisfiability, but omitting Ch. 9's unanimity interaction. It is
  the fifth. The Outlook's version was stale for the same reason. **If a sixth
  instance is ever added, four places need updating: Secs. 7.4, 9.6, 13.4 and
  the Outlook's Sec. 15.2.**
- **The "cardinality two is degenerate" thread was unconnected across chapters.**
  Ch. 12 counts it to three (Secs. 10.3, 11.5, 12.4); Ch. 13's `k = 2` is a
  fourth in a different currency — a clause of two literals rather than a
  complex of two members — and now says so.
- **The software table had 13 computed equations missing**, all the Part II
  structure-by-structure thresholds and everything added after it was written
  (Chs. 12, 13's nested sections, 14's parent-to-child update). Its rows were
  also out of numeric order. Both fixed; coverage is now 70 of 98 equation
  labels, the remaining 28 being definitional or introductory.
- **Seven figures and tables were never referenced in the prose** (Chs. 1, 3 and
  12) while every other figure in the book is. All now have a pointer.
- The preface said "two of the chapters draw on published work", written when
  the book was shorter; it now names Part II as the published part and says
  Chs. 12–13 have no paper behind them.
- Five references added for Chs. 12–13 were never cited, and two pre-existing
  ones (Liu *et al.* on core percolation, Yoon *et al.* on Ising with motifs)
  belonged in Chs. 11 and 9. All now cited.

Checked and clean: terminology (no variant spellings of chy-degree, treelike,
hyperedge, higher-order), all `\ref`/`\eqref` resolve, no figure declared
without a file or present without a use, and the numbers repeated across
chapters agree — `u' = t/(1-t+t^2)` in Chs. 7, 8, 9; `-13.9%`/`-14.2%` in
Chs. 1, 4, 9, 15; the 77–100% HRG core in Chs. 11, 14, 15. (The last of these
was *not* in fact consistent — see the 2026-08-30 pass below, which found Chs. 2
and 3 quoting 65–83% for the same quantity.)

### Sec. 12.5, a graph with triangles, 2026-08-30

Author's suggestion, and it holds --- with a bound on it that the section states
rather than hides.

**The claim.** For *proper* colouring a triangle-clustered graph is an ordinary
graph: a triangle forbids exactly what its three edges already forbid. So the
graph calculation applies to it, and gets it wrong. On a regular network where
every vertex sits in `kappa` triangles, degree `n = 2 kappa`,

    chygraph:  n = (q-1)^2 + 2        graph reading of the same object:  n = (q-1)^2 + 1

off by one in degree at every `q`, because the graph calculation adds a vertex's
neighbours as independent when two of them are adjacent inside each triangle.
The chygraph puts that loop inside the interior sum, where it is exact.

**The internal check.** Both lines are one formula at its two ends. Links and
triangles together at fixed degree, with `f` the fraction of degree through
triangles, give `s(n-3) + 2 s^2 [n(1-f/2) - 1] = 1`, Eq. (12.6), which returns
`(q-1)^2 + 1` at `f = 0` and `(q-1)^2 + 2` at `f = 1` **exactly, in rational
arithmetic, with no layer dropped**. The cross term does all the interpolating.

**The bound, and it is the reason the section is worth having.** One in degree,
and *only for regular ensembles*. For Poisson layers `kbar = kappa`, the cross
term vanishes identically, and Eq. (12.6) returns `(q-1)^2` at every `f` --
no difference at all. So this is **not** "clustering moves the colouring
threshold". It is the regular ensemble's lost branch, the same one separating
`(q-1)^2` from `(q-1)^2 + 1` in Sec. 12.3, arriving where it could be mistaken
for physics. Set beside Ch. 9's 13.9% in `T_c` and Sec. 4.5's fifth of `q_c`,
this is the third and least dramatic sign the book has taken on that comparison,
and the section says so. **Do not restate it as a clustering effect.**

`check_triangles_at_fixed_degree` asserts all three: the endpoints exactly,
monotonicity in `f`, and the Poisson null. Any of the three breaking fails it.

### Colouring: what is ours and what is not, 2026-08-30

Prompted by the author asking whether anything in Ch. 12 is actually new. The
honest answer needed two corrections and produced one result.

**The hypergraph threshold is not ours.** Gabrié, Dani, Semerjian & Zdeborová,
*J. Phys. A* **50**, 505002 (2017), Eq. (43) gives
`l_stab = (q^(K-1) - 1)^2 / (K-1)` for Poisson degree, which is Ch. 12's
`(q^(c-1)-1)^2` neighbours per node. It reproduces **all six** of their tabulated
values (4.50, 16.33, 56.25, 32.00, 225.33, 112.50). Their Appendix B derives the
transmission at finite temperature and Eq. (12.2) is its `beta -> infinity`
limit. Cited now, and Sec. 12.4 says the agreement is a check rather than a
result. **The README previously listed this as one of two results "worth
guarding" as new. It was wrong to.**

**Sec. 12.6 recapitulates.** Survey propagation for colouring is Braunstein,
Mulet, Pagnani, Weigt & Zecchina (2003); the thresholds it returns are theirs.
The section now says so in its opening and again in its closing qualification,
because a section showing five-figure agreement is easily mistaken for one that
discovers something.

**What is left to Ch. 12** is the proper-colouring cardinality independence,
`tau = -1/(q-1)` at every `c`, and the contrast between the two rules — that one
substitution gives both and they part company only above cardinality two.

**Sec. 12.7, new.** The contrast has a consequence one level up, and this is the
chapter's own:

- **Hypergraph rule: a complex forbids at most ONE colour**, since its other
  members cannot all be forced to two at once. So the cardinality-two
  inclusion–exclusion carries over verbatim under `e/q -> h/q` with
  `h = q prod_{j<c}(e_j/q)`, Eq. (12.7). Cardinality enters the one-step
  calculation through one scalar and nothing else.
- **Proper rule: up to `c-1` at once**, so the message is a forbidden *set*,
  correlated within a complex as well as across complexes, and the node update
  needs that set's joint distribution. No substitution repairs it. Not
  attempted.

That is Sec. 12.4's "cardinality two is degenerate" thread arriving in the 1RSB
machinery: a two-member complex has no room for a partial violation and equally
none to forbid two things at once, and it is the second that decides whether the
machinery can be reused.

Checked against Gabrié's `l_col`, Table 12.3:

| q | c | Sigma = 0 at | published | error |
|---|---|---|---|---|
| 3 | 3 | 26.900 (0.000) | 26.92 | −0.07% |
| 4 | 3 | 63.276 (0.024) | 63.30 | −0.04% |
| 2 | 5 | 52.283 (0.043) | 52.32 | −0.07% |

All three errors agree in sign and size, so that is a converged finite-population
bias, not scatter — the caption says so rather than claiming agreement to a tenth
of a per cent.

**Two exclusions, both stated in the text so neither reads as an oversight.**
`q = 2` at `c = 3, 4` is omitted because Gabrié et al. mark those thresholds
invalid under an SP type-I instability; the code returns numbers there and they
would mean nothing. The proper rule is absent because of the argument above, not
because it was tried.

**Two bugs worth remembering**, since both are invisible in `Sigma` and show up
first in `<e>`: the per-colour weight is `h/q` and dropping the `/q` saturates
everything to `<e> = 1`; and the non-trivial branch is reached from **above**, so
initialising low gives `<e> = 0` and nothing else. `hyp_converge` starts at
0.97–1.0 and its docstring says why.

### Survey propagation, 2026-08-30

Secs. 12.6 and 13.9, new. The book previously said in six places that it did not
do the one-step calculation; it now does it at `m = 0` for the two chapters
whose messages are warnings, and the disclaimers are revised rather than
deleted, because for everything else they are still true.

**What is computed.** A survey is a distribution over warnings, taken over
clusters. Ch. 13's Eq. (13.4) already builds the numerator of the update — the
weight of "forced to falsify" — and SP divides it by the weight of the three
non-contradictory states. The threshold is not where the survey branch appears
but where the clusters run out, so both sections compute the complexity `Sigma`
as a Bethe count at `m = 0` and bisect `Sigma = 0`.

| | computed | published | error |
|---|---|---|---|
| `alpha_s`, k = 3 | 4.281 (0.008) | 4.267 | +0.3% |
| `alpha_s`, k = 4 | 9.906 (0.045) | 9.931 | −0.3% |
| `c_q`, q = 3 | 4.683 (0.012) | 4.69 | −0.1% |
| `c_q`, q = 4 | 8.901 (0.012) | 8.90 | +0.0% |
| `c_q`, q = 5 | 13.660 (0.008) | 13.69 | −0.2% |

Brackets are the spread over three seeds. No published value enters any of the
calculations. Ch. 12 also gets `c_d` as a by-product (4.46 vs 4.42, 8.49 vs
8.27) — a per cent and then three, against a tenth of a per cent for `c_q`,
because a vanishing complexity is a crossing and a branch appearing is a fold.

**Two things worth guarding.**

- **The scalar closure does not survive the step to surveys.** This is the whole
  reason the sections carry populations. Solving the SP equation with `eta` kept
  a single number gives `eta* = 0.1316` at k = 3, `alpha = 4.267`, against
  0.1828 for the distribution — 28% low, and 8% low at k = 4. `eta` stops being
  a probability and becomes a continuous variable whose spread enters the next
  update. Measured in `check_survey_needs_a_population`.
- **Colour symmetry collapses the survey to one scalar but does NOT decouple the
  colours.** A forced neighbour forbids exactly one colour, so "every other
  colour forbidden" is a coverage problem and needs the inclusion–exclusion of
  Eq. (12.5), not a product over colours. The independent-colour shortcut is
  wrong and a single threshold would hide it; the check asks for three.

**What is deliberately not claimed.** All three qualifications are stated in
both sections: it is `m = 0`, so condensation (k >= 4, and the rigidity
transitions) is outside it; Eq. (8.3)'s treelike condition is untouched, so
Ch. 14's obstruction is unaffected; and that the flat graph case reproduces
published constants is an implementation check, not a result. The Outlook's
1RSB entry is rewritten around `m != 0` rather than deleted, and Ch. 10's
disclaimer now says why a hitting-set field is the harder object — a
distribution over distributions rather than over a finite alphabet.

**Runtimes.** `figs/satisfiability.py`'s SP check is ~5 min, `figs/colouring.py`'s
~3.5 min. These are now the two slowest scripts after `merge.py`.

**Still open, and it is the honest gap:** nothing here is validated against
anything but flat `k`-SAT and flat `q`-COL. The nested clauses of Secs. 13.6–13.7
are the case with no flat counterpart, and running SP on the three-level chygraph
against SP on its CNF flattening is the experiment the machinery was built for.
It was attempted and is **unresolved** — see `../probe/NESTED_SP.md` and
`../probe/nested_sp.py`. Two of three parts work and are worth keeping: the
cross-level message passing and the cluster counting both validate against the
`k1 = 1` relay identity, each confirmed by clean `1/sqrt(N)` decay, the second
establishing that a layer-1 sub-expression contributes no factor of its own. The
third does not: the nested complexity appears discontinuously at a finite value
and has no usable zero before the `delta` channel saturates, so there is a
threshold on the flattened side and only a sign on the nested one. Nothing from
it is in the book.

### Editorial review, addressed 2026-08-30

`~/Downloads/chygraph_statmech_review.md` — an external editorial review of the
214-page build. Its verdict was accept-subject-to-revision, with the revisions
"almost entirely apparatus and production, not science"; it independently
re-derived about thirty of the book's analytical results and found no errors.
Every one of its line-specific claims that I checked was accurate, including a
third "next nine chapters" in `data.tex:6` that the consistency pass above had
missed. What was done:

**Front matter (its §3.2–3.5).**
- `\date{}` added. `\maketitle` had no date, so the title page was printing the
  LaTeX run date and changing it on every build.
- `\frontmatter` / `\mainmatter` / `\backmatter` added, and
  `\pagenumbering{gobble}` plus `introduction.tex`'s `\setcounter{page}{1}`
  removed — those two did the same job and would now fight. The preface, the
  software chapter and its Table 1 have roman page numbers and appear in the
  contents with them; Part I starts arabic 1.
- **An index**, `makeidx` plus 116 `\index` entries giving 82 terms, with the
  defining page in bold. Entries sit immediately after the `\label` of the
  section that owns them, so they move with the text. The `theindex`
  environment is redefined to one ragged-right column: the class default is two,
  which on a 5-inch page leaves ~1.7in of measure and produced five overfull
  lines on entries like "Fortuin--Kasteleyn correspondence".
- **A Notation chapter** (`notation.tex`), the review's §4.1: three tables —
  the structural symbols, the ten letters that carry more than one meaning with
  the chapters each meaning is in force in, and the chapter-local letters.

**Ch. 11's two symbol collisions (§4.1).** `\kappa` was the chy-degree
everywhere in the book *and* the "complex forces the entering node out"
probability in Eq. (11.4), which also contains `\bar\Phi`; `\lambda` was
Ch. 5's Perron root and Ch. 11's leaf-ready probability, including in a
displayed threshold condition. Renamed to `\xi` and `\psi` (neither was used
anywhere in the book), with a sentence saying that Ref. [liu2012] writes them
the other way.

**The nine orphan result tables (§4.2)** are now numbered floats with captions,
each with a pointer from the prose: Tables 6.1, 7.1, 9.1, 9.2, 12.1, 13.1,
13.2, 14.1, 14.2. The `-13.9%` clustering result is Table 9.2. Captions are
written to stand alone, since a float can drift from the paragraph that reads
its rows aloud.

**Bibliography (§4.4).** All twelve incomplete entries fixed, with nothing
invented — four resolved from the PDFs under `~/Downloads/chygraph_references/`,
three looked up against the published record, four were simply the wrong BibTeX
type. Three corrections are substantive and worth knowing about:
- **`cirigliano2025` had the wrong title.** The bib said "Percolation on
  clustered networks with static triadic closure"; the paper the book actually
  uses, and has on disk, is *How universal is the mean-field universality class
  for percolation in complex networks?*, arXiv:2506.17175.
- **`ha2025` had the wrong title.** Published as *Connected components in
  networks with higher-order interactions*, J. Phys. Complex. **6**, 045006.
- **`fujiki2024` is a 2018 paper** — Phys. Rev. E **97**, 062308. Key renamed
  `fujiki2018` and the one citation updated.
- `keating2026` had no identifier at all; it is arXiv:2511.15688, whose abstract
  says "Group effects alone, without long cycles, produce standard continuous
  phase transitions" — exactly what Sec. 5.4 attributes to it. Title corrected
  to *Loops, not groups*.
- `karp1972` was an `@article` carrying both `booktitle` and `journal`; it is an
  `@incollection`. `bianconi2021` was an `@article` whose journal was its
  publisher; it is a `@book`. `joslyn2017` and `leicht2009` had arXiv numbers in
  the `journal` field; both are `@misc` with the identifier in `note`.
- Left alone deliberately: `bianconi2017` and `montanari2008` have no volume
  because J. Stat. Mech. cites by article ID (034001, P04004), and `son2026` has
  no volume or pages because it is genuinely in press — its DOI is now recorded.

**The three substantive points (§4.5).**
- **Ch. 5's `5/3`.** The box said the published coefficient "should be 5/3" and
  then displayed `5q(1-q)^2`, so a reader could not check the correction. Both
  are right: 5/3 is the mean component size in the one-bond configuration (the
  three members sit in components of sizes 2, 2, 1) and 5/3 x 3q(1-q)^2 =
  5q(1-q)^2. The box now shows the uncollapsed form and says where the 5/3 comes
  from. **Still open for the author: whether an erratum has been filed against
  `vazquez2024comnet`.**
- **Ch. 7's percentages** are a fraction of `v_c`, not of the transition. The
  text now says so and gives the other normalisation (5.4 / 29.5% at q = 3 / 6)
  so the convention cannot be misread. Unchanged otherwise — see the decision
  above to quote the `v_c`-normalised set.
- **The merge table's `<k> = 1.5` row** decreases from n = 2000 to n = 4000 while
  labelled "grows". Text and caption now say that three seeds do not separate
  trend from scatter at that density and that the classification rests on the
  slope against n in Fig. 14.3, not on those three entries.

**Copy-editing (§7).** `software.tex`'s "is does"; the Outlook's "Thirteen
chapters" (it is Chapter 15, so fourteen precede it); `data.tex:6`'s third
"next nine chapters"; `chygraphs.tex`'s `Sec.~\ref{ch:giant}` (a chapter label
behind "Sec.", the book's only cross-reference type mismatch, now pointed at
`sec:joint` which is what the sentence means); `main.tex`'s "fourteen chapters"
comment; the orphan `figs/cl-1.png`. The preface now credits Ch. 2 as well as
Part II to the two papers, which the review noted. The ~50 unused `\label`s were
left as they are — they are defensive and harmless.

**Not done, and why.**
1. **`chygraph_statmech` is still private** — the review's one hard blocker, and
   `../TODO.md` item 3. `software.tex` prints the URL as if it resolves and the
   preface stakes the reproducibility claim on it. This needs the author: either
   make the repository public or rewrite the software chapter and preface.
2. **The dedication placeholder still prints.** It needs the author's words.
3. **No list of figures, list of tables, half-title, copyright/CIP page,
   acknowledgements or data-availability statement.** CIP and copyright are the
   publisher's. A list of figures and tables would contradict the standing
   decision that the contents lists parts and chapters only; that is the
   author's call to reopen.
4. **Chs. 12 and 13 have not been refereed** (§4.3). Informational; the review
   suggests routing them to a specialist referee.

### Consistency pass, 2026-08-30

A second mechanical pass over the whole book, run after the technical-correctness
pass. Everything below was found and fixed.

- **`sec:hrg` was defined twice** — Sec. 3.4 and Sec. 11.6 are both called
  "Hyperbolic random graphs" and both carried the label, so LaTeX resolved all
  seven references to the later one. The two references in `data.tex` that meant
  their own Sec. 3.4 printed **"Section 11.6"**. Ch. 3's copy is now
  `sec:hrg-data`; the five references that really do mean 11.6 are unchanged.
  Note that `main.log`'s `multiply defined` warning is **not** caught by the
  build recipe's `grep -i undefined` — the verify block above now checks for it.
- **The core-recovery headline was quoted at two values.** Chs. 11, 14 and 15
  said the chygraph accounts for **77–100%** of the leaf-removal core; Chs. 2 and
  3 said **65–83%** for the same quantity. Recomputed from
  `probe/results/prediction4.csv`, the per-seed ratios are 0.77–1.00 at
  `tau = 2.9`, 0.64–0.85 at `tau = 2.5` and 0.68–1.45 at `tau = 2.1` (excluded as
  uninformative by Sec. 3.4). So 77–100 is the `tau = 2.9` slice — which is what
  `figs/cover.py:check_hrg` prints, it hard-codes `tau - 2.9` — and 65–83 is
  roughly the `tau = 2.5` slice. **Standardised on 77–100 everywhere.** Live
  caveat: no site states the `tau`, and over both informative rows the range is
  64–100%. If that scoping is ever tightened, five places carry the number —
  `chygraphs.tex`, `data.tex`, `cover.tex`, `overlap.tex`, `outlook.tex`.
- **Ch. 8 decomposed the null-model gap the wrong way.** It said "Most of that
  gap is the excess chy-degrees ... The cross term is the rest ... a *further* 23
  per cent". The two effects oppose. Poisson 8.4256, regular 6.8801, regular
  without the cross term 5.2951: the gap is 1.55, the cross term alone is 1.58
  and it pushes *back up*, so the excess-degree change is 203% of the gap and
  the cross term −103% of it. The individual numbers and the 23% were right;
  only the framing was wrong. Rewritten to say the gap is a residue, not a sum.
- **Ch. 9 undercounted the "name what is held fixed" thread** — "three times
  already" where the Outlook counts five with Ch. 9 as the fifth. Sec. 5.5's
  identical-marginals demonstration was the one dropped. Now "four times".
- **`data.tex` said "eight of the ten networks"** at `shared_2+ <~ 0.02`.
  Table 3.2 has seven at `<= 0.020` and nine at `<= 0.022` (Vidal and Figeys both
  sit at 0.022); the paragraph three pages earlier puts nine in the group. Now
  nine.
- **Two chapter counts stopped one short of Ch. 13.** `chygraphs.tex`'s "the next
  ten chapters" (from Ch. 2) and `data.tex`'s "the next nine chapters" (from
  Ch. 3, twice) both reached only Ch. 12, though both sentences contrast the span
  with Ch. 14 as the obstruction. Now eleven and ten.
- **The preface had an ungrammatical sentence**, introduced in commit `e29852d`
  and never caught: "And once a complex is allowed to hold other complexes rather
  than only nodes, formulas that no flat encoding can hold without damage." No
  verb. Repaired, and satisfiability added to the list of models in the same
  sentence, which had omitted it.
- **Table 3.3's caption stated a notation convention the book breaks twice.** It
  said "`theta` is used for structural exponents throughout, `beta` being reserved
  for the inverse temperature." Ch. 5 uses `beta` for the order parameter exponent
  four times, and `theta` in Chs. 4–5 is the threshold determinant, not an
  exponent. The caption now says which letter means what where, rather than
  claiming a convention that does not hold.
- **Eq. (5.18) was missing from the software table**, whose caption is "every
  computed equation". It is a closed-form threshold checked against Cirigliano's
  exact result by `test_stc_threshold_is_cirigliano_eq12`, and it was the only
  `\label{eq:` added by the technical-correctness pass. Added, in numeric order.
- **`mezard2009` and `newman2001random` were in the bibliography and never
  cited**, though `main.tex`'s `\include` comment names both as the
  introduction's textbook background. Both now cited in Ch. 1 — Newman, Strogatz
  and Watts at the Molloy–Reed criterion, Mézard and Montanari at the cavity
  recursion.
- **Ch. 1 now cites Leone, Vázquez, Vespignani and Zecchina**, *Eur. Phys. J. B*
  **28**, 191 (2002), at the end of Sec. 1.3.1 — the exact Ising `T_c` and
  critical exponents on a random graph with arbitrary degree distribution. Added
  because Eq. (1.5) was stated with no source, and because the paper's divergent
  `<k^2>` case (ordered at every temperature) is the exact counterpart of the
  `p_c = 0` the chapter states for percolation two pages earlier, which is the
  parallel Sec. 1.3.1 is built on. PDF under `~/Downloads/chygraph_references/`.
- **Ch. 7 now states the fold in closed form.** The README already claimed the
  chapter gave the ordered spinodal "at exactly `2 sqrt(q-1)`"; it gave only the
  four numbers. The closed form is now **Eq. (7.11)**, proved by taking the
  resultant of the fixed-point equation and its derivative, which factors as
  `v^2 (q-1)^2 (q-v)^2 (q+v)^2 (v^2 - 4q + 4)` — the `(q-v)^2` factor being the
  disordered spinodal and `v^2 - 4q + 4` the fold. Verified numerically to
  1e-10 at `q = 2, 3, 4, 5, 6, 8, 10, 16`.

Software-table coverage is now **72 of 100** equation labels (was 70 of 98), the
remaining 28 being definitional or introductory.

Checked and clean this pass: all 223 `\ref`/`\eqref` resolve, no duplicate
labels of any kind, every one of the 49 figures and tables referenced in the
prose, all 20 `\includegraphics` files present with no orphan PDFs in `figs/`,
all 93 `\code` names in `software.tex` resolving across both repositories
(`joint.JointGiantComponent` is an alias at `chygraph/src/chygraph/joint.py:237`,
not a `def`, so a naive resolver will report it missing), no hardcoded internal
cross-reference numbers anywhere in the prose, terminology still single-valued,
and the remaining repeated numbers agreeing across chapters — `u'` in Chs. 7/8/9,
−13.9%/−14.2% in Chs. 1/4/9/15, `q_c` 1/3 → 0.4030 in Chs. 4/15, `T_c`
0.2500 → 0.2929 in Chs. 6/15, 5.1–22.8% in Chs. 7/15, six-of-twenty GBP
non-convergence in Chs. 14/15, the tricritical boundary, and the `(q-1)^2`
colouring thresholds. The three running counts in Secs. 7.4, 9.6 and 13.4 are
correct as running counts, and Secs. 5.4 and 6.4 remain aligned on continuity.

### Language pass, 2026-08-28

The register is: report what is done or achieved, no decoration around it,
without compromising understanding. Sixty-odd edits across all nineteen files,
removing three classes of construction:

- **"It is worth X-ing, because Y"** — 34 instances, now 8, and the survivors are
  literal uses (`two chapters' worth of argument`, `k1 alpha worth of clauses`).
  The fix is always the same: state the thing and drop the recommendation to
  find it interesting.
- **Staging and salesmanship in chapter openings** — "this is the chapter the
  whole book has been deferring to", "the reader is entitled to ask", "worth the
  chapter on their own", "saying which three is most of the work". Openings now
  state what the chapter does in its first sentence.
- **Evaluation standing in for statement** — "the sharpest statement in the
  chapter", "something uncomfortable", "the cleanest control and why it is the
  cleanest", "Stop and look at that", "Read that table carefully". Replaced by
  the statement itself.

**Keep this register when editing.** The distinction that matters: explanation
of *why* a step works is not decoration and should stay; instruction to the
reader on how to feel about it is. Connective tissue that makes an argument
followable was left alone.

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
the triangle construction carries graph excess degree `<k_bar> = n + 1`, not
`n`, while its Poisson chy-degree layer loses none of its branches).

The second thread appears five times: **a linear instability is not always the
transition.** Sec. 5.4
(interdependent networks appear by saddle-node while `Q = 1` is still stable),
Sec. 6.4 (threshold contagion, same exclusion, no invasion from one seed),
Sec. 7.4 (the ferromagnetic Potts transition is first order above `q = 2`, so
the branching condition is a spinodal — the transition is 5.1% to 22.8% below
it at `q = 3` to `6`),
Sec. 9.6 (the unanimity interaction, where the free energy of Sec. 8.6 is
finally *used*: the q = 16 double transition has its first-order point at
`T_c = 0.9790` inside a window whose lower edge is the spinodal), and Sec. 13.4
(satisfiability above `k = 2`, the extreme case — the warning map's derivative
at the trivial fixed point is *exactly zero at every clause density*, so there
is no linear instability anywhere to mistake for a transition). Sec. 9.6 also
carries the practical corollary — critical slowing down mimics a coexistence
window, so the test needs a free-energy crossing, not a surviving branch. Keep
the five consistent; they are one statement.


## Underlying material

- `~/Dropbox/submissions/hyperabs.2022/hyperabs_v3.tex` — Phys. Rev. E **107**, 024316 (2023)
- `~/Dropbox/submissions/chygraph.2023/chygraph.tex` — J. Complex Netw. (2024), cnae047
- `~/av2atg/chygraph/manuscript_3/manuscript.tex` — the percolation extension.
  **Will not be submitted** (author, 2026-08-30): the book is the first place
  Ch. 5's critical amplitude and moment hierarchy are reported, and Ch. 5's
  opening says so. Do not add a citation to it.
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
line. All of that is fast — exact rational arithmetic over small complexes.
The **survey-propagation sections added 2026-08-30** are not. Three bisections
of the complexity over three seeds each — Sec. 12.6's `c_q`, the clustering
by-product, and Sec. 12.7's hypergraph thresholds against Gabrié et al. — put
the whole script at about **fourteen minutes**, which makes it the slowest in
the book, ahead of `merge.py`.
`figs/satisfiability.py` generates Figure 13.1 and runs Ch. 13's checks: the
clause interior against enumeration in rational arithmetic (including the
non-uniform emitted message that rules out a symmetric fixed point); the `k = 2`
linearisation and the branch appearing continuously at `alpha = 1`; that the
`k >= 3` derivative is *exactly* zero, tested by the rate at which the finite
difference falls (`eps^(k-2)`, checked by halving `eps`) rather than by an
absolute bound, which would only show `eps` was small; and both conventions for
a contradicted variable. It also carries the **survey-propagation
section added 2026-08-30**: the SP update, the complexity `Sigma` as a Bethe
count of clusters at `m = 0`, and `Sigma = 0` bisected against the published
`alpha_s`. That check averages over three seeds and takes about five minutes,
which makes this the slowest script in the book after `merge.py`.
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
