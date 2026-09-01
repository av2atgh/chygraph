# Complex hypergraphs (book)

Book-length treatment of percolation and statistical mechanics on complex
hypergraphs. Production style follows *Local network growth* — 5×8 in trim,
grayscale, heavy on TikZ illustration. Derivations run in continuous text, each
ending on the routine that executes it; the `calculation` boxes that style
carried were removed on 2026-09-01, and the environment with them.

## Contents

| | Source | Topic |
|---|--------|-------|
| — | `preface.tex` | **Drafted.** The tree assumption, the bargain, what the book claims |
| — | `software.tex` | **Drafted.** Repo links; equation-to-method-to-test table; how to reproduce a figure |
| **I** | | **Foundations** |
| 1 | `introduction.tex` | **Drafted.** Percolation; Ising, vertex cover, hitting set; mean field, cavity, replicas, BP; why real networks break treelikeness |
| 2 | `chygraphs.tex` | **Drafted.** The object: complexes whose vertices are complexes; every higher-order structure as one thing |
| 3 | `data.tex` | **Drafted.** Chygraph representation of real systems: papers, protein complexes, reactions, schedules; cliques of a clustered graph; how much overlap there is |
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
| 12 | `colouring.tex` | **Drafted.** Proper against hypergraph colouring; `tau = -1/(q-1)` independent of cardinality; `(q-1)^2`; a graph with triangles, where the graph calculation is wrong by one; survey propagation; and **Sec. 12.9, extending Krzakala et al. to chygraphs to get `c_q` for a clustered graph** — the chapter's own result |
| 13 | `satisfiability.tex` | **Drafted.** Clauses as complexes; `alpha = 1` exact at `k = 2`, no linear instability above it; **clauses whose members are clauses**, and what CNF flattening costs |
| **IV** | | **Complexes with non-trivial overlap** |
| 14 | `overlap.tex` | **Drafted.** The price of treelikeness: what `BP_chi` loses to overlapping cliques, measured on three classes of network, with the rewired control that attributes it |
| 15 | `gbp.tex` | **Drafted.** The repair the literature already has: region graphs, Möbius counting, GBP over 240 runs — exact where the clique family is chordal, and **chordality follows provenance, not clustering** |
| 16 | `metacomplex.tex` | **Drafted.** The repair that stays inside the formalism: merge complexes sharing two or more atoms; exact iff the merged incidence structure is a forest; and **Sec. 16.3, the core-percolation transition at incidence branching one** — the chapter's own result |
| 17 | `outlook.tex` | **Drafted.** One recursion, many models; the two running threads; what is not done |

**Chapters 12 and 13 have no manuscript behind them.** Every other chapter is
exposition of work that exists elsewhere; for these two the calculations are
done in the book's own figure scripts, `figs/colouring.py` and
`figs/satisfiability.py`, and checked against the published thresholds --- which
since 2026-08-30 the last section of each chapter **computes** rather than
quotes, by survey propagation at `m = 0`. PDFs of
the references are under `~/Downloads/chygraph_references/`.

## Status

Last updated 2026-09-01. `main.pdf` builds with **0 errors, 0 undefined
references and 0 multiply-defined labels, across 284 pages.** Not box-clean:
**three overfull hboxes** — `cover.tex:239--249` (1.98pt, "Which
replica-symmetry-breaking point"), `overlap.tex:146--155` (6.37pt, the
`\paragraph{Real networks, where the triangles are whatever was recorded.}`
head) and `metacomplex.tex:347--353` (3.16pt) — and **33 underfull vboxes**,
every one of them `while \output is active` — page-breaking around floats, not
a line that runs into the margin.
55 figures, 25 numbered tables, 125 numbered equations of
which 116 carry labels, 73 references, a 96-term index. All 295 distinct
`\ref`/`\eqref` targets resolve, all 28 `\includegraphics` files are present,
and all 114 `\code` names in `software.tex` resolve across both repositories.
Software-table coverage is 85 of the 116 equation labels. Both checks under
*Two checks the build cannot make* print nothing.

These counts were re-measured on 2026-09-01. Six of them had drifted from what
this section claimed on 2026-08-31 without any content being added: numbered
equations 115 to 124, labelled 113 to 115, tables 23 to 24, references 72 to 73,
`\ref` targets 281 to 284, `\code` names 109 to 114. The 2026-08-31 pass
before that had found worse (228 pages, 0 overfull, 41 figures, 104 equations).
**Re-measure rather than trust them after any drafting pass** —
the recipe is under Building, below. The numbers inside the dated pass sections
further down are left as they were: they record what was true at the time.

**The book is drafted end to end.** Preface, the software chapter, and
seventeen numbered chapters, all with prose, figures and checks.

`software.tex` sits in the **back matter**, after Ch. 17 and before the index,
and is unnumbered — "The software", page 243, arabic. (This section used to say
it sits after the preface with roman page numbers; it has not since
`\backmatter` was added.) It carries the two repository URLs and **Table 1:
every computed equation, the routine that evaluates it, and the test or script
that checks it** — the book's version of the retired supplement's Sec. I. Every
one of its 82 equation labels and all 105 of its `\code` names were re-verified
to resolve on 2026-08-31, the code names across *both* repositories; if a
routine is renamed, that table is where it has to be fixed.

**It is Table 1 again, as of 2026-08-31, and the fix is in `main.tex`.** It had
been printing as **"Table 15.2"**: it is a `longtable` in the back matter, and
`\backmatter` stops numbering chapters but neither resets the float counters nor
clears `\thechapter`, so `\thetable` still expanded to `15.<n>`. The Outlook's
own table is 15.1, so the software map took 15.2 and its three
`Table~\ref{tab:map}` call-outs sent a reader into Chapter 15 to look for it.
(Those were the numbers when Part IV was one chapter; the Outlook is now
Ch. 17 and its table is 17.1.) Immediately after `\backmatter`, `main.tex` now
resets the `table` and `figure` counters and redefines `\thetable`/`\thefigure` to a plain `\arabic`, with a
comment saying why. **Add a float to the back matter and it will number plainly
— that is deliberate.** Ch. 12's caption cites `\citet{gabrie2017}, Table~1`,
which is somebody else's Table 1 and is attributed in place, so the two do not
collide.

### What to do next

The remaining work is revision, not drafting:

| | what |
|---|---|
| 1 | ~~A read-through for continuity~~ — mechanical consistency passes were done 2026-08-28 and 2026-08-30, and a language pass 2026-08-28 (see below). What none of them could do is judge the argument: nobody has yet *read* 1–17 end to end for whether it persuades. |
| 2 | ~~Fix the manuscripts~~ — **done by deletion.** The three errors the book found (the tricritical cardinality, the GBP convergence count, the Ch. 5 arithmetic slip) were in `main.tex`/`supplement.tex`, which no longer exist. The book carries the corrected values and says so; the entries below are kept as a record of what was wrong, in case either file is ever resurrected from git. |
| 3 | **Chs. 12 and 13 are the only chapters whose results are not backed by a manuscript.** If either is to be published separately, the calculations in `figs/colouring.py` and `figs/satisfiability.py` are the starting point. |
| 4 | **Front matter: the dedication is still `\itshape Dedication to come.` in `main.tex`, and it prints.** The editorial review calls this a blocker. It needs the author's words; nothing else about the front matter is now outstanding. |
| 5 | **Make the repositories public.** `software.tex` prints both URLs as if they resolve; `chygraph_statmech` is still private. This is `../statmech/TODO.md` item 3 and now blocks the book as well as the manuscript. |

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
  seventeen chapters need them. A Part IV called "Limitations",
  holding two chapters, was removed on 2026-08-31 (author's call). What stands
  now is a different object: **Part IV, "Complexes with non-trivial overlap"**,
  three chapters (14 failure, 15 GBP, 16 meta-complexes) added 2026-08-31, each
  carrying its own measurement, figures and checks. The Outlook (Ch. 17) closes
  the book rather than a part of it. Do not re-merge the three back into one.
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
  graph has no giant component. So the feasibility of Ch. 14's repair is a
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
  Recomputed from the cached `../statmech/probe/results/gbp_cliques.json`, the residuals
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

### Consistency and style pass, 2026-09-01

Read all twenty-one `.tex` files against each other and against this README.
Commit `8f6903a`, 19 files, +109/-99. The prose needed almost nothing: a grep
for `It is important to note`, `Of course`, `Indeed`, `Crucially`, `As we shall
see` and the rest returns **one** hit in 74k words, and British spelling is
uniform (`colour`/`neighbour`/`behaviour`, zero `-ize` forms). What it found was
text that had broken, terminology that had drifted, and cross-references that
Part IV's split into three chapters left behind.

**Four sentences were broken outright.** `giant.tex` had lost the end of a
clause mid-line — "group structure on its own produces only continuous" ran
straight into the next sentence, and "The qualifier is doing work" then pointed
at a qualifier the surviving text no longer isolated. `satisfiability.tex` had
`\paragraph{What does not survive the step.} The scalar closure does.`, which
reads as the opposite of what it means; `statmech.tex`'s parallel passage gets
it right with "The scalar closure goes." Also "Equation above is the derivative"
(`ising.tex`) and "\citet{krzakala2004} solve ... and locates"
(`colouring.tex`). **None of these is catchable mechanically.** LaTeX compiles a
truncated sentence without complaint, and there is no grep for a sentence that
means the opposite of its own heading; all four were found by reading.

**`chygraphs.tex` contradicted itself within a page.** The paragraph "The
pairwise test is necessary and not sufficient" establishes that sharing at most
one atom does *not* guarantee treelikeness — the ring of four triangles — and
eight lines later the chapter's bolded takeaway read "a chygraph is locally
treelike **when** its complexes meet in at most one atom". Now "only if", with
the case where it is also sufficient named.

**The condition itself was stated four ways** — atom, node, vertex, member —
across seventeen sites in nine files, sometimes alternating inside one file
(`metacomplex.tex` had node at :126 and atom at :186; `data.tex` had nodes in
the text and atoms in the caption of the same table). `notation.tex` defines
*atom*, so atom it is everywhere now.

**`software.tex` had not followed Part IV's split.** Table 1 filed the
region-graph and meta-complex equations under Ch. 14's heading; the region-graph
calculation was credited to Ch. 14 rather than 15; `figs/` was described as one
script per chapter (it is thirteen scripts covering Ch. 3 onward); three
chapters were said to have no manuscript behind them rather than two. Its
slow-script paragraph still named `ising.py` and `potts.py` as the only slow
ones — see below.

**`preface.tex`** described Part IV as one chapter that "shows what a
region-graph treatment recovers on a single instance"; it is three chapters and
240 runs. It also listed overlapping complexes inside "the whole of Part III".

Smaller: the preface's closing sentence about the field's predicament appeared
verbatim in `introduction.tex`; `overlap.tex` restated its own figure caption in
the following paragraph; `potts.tex` named one section twice in consecutive
sentences; `cover.tex` referred to "the fourth model" where Ch. 1 introduces
three; `satisfiability.tex` called itself the fourth hard-core family where
`colouring.tex` counts it as half of the third; `regime` was spelled three ways,
one of them the book's only non-ASCII byte; `gbp.tex` carried the only nine
unspaced em-dashes against 500-odd spaced ones elsewhere.

**On style**, four things were cut: `Our core idea is geometric before it is a
calculation` opened the cavity-method section, which is Bethe's idea and not
this book's; "The moral is..." had become a section-closing formula, four times;
the unanimity-interaction/threshold-contagion relationship was described three
times as being *in different clothes*; and four isolated flourishes went ("the
hinge of the book", "the kind of coincidence the whole book is arranged to
produce", "where that bill arrives", and the remark that the definition's
sentence is self-referential in the way its subject is).

**`figs/merge.py` was wrong in this README three times over.** It was called the
slowest script in the book, and ranked inconsistently against `colouring.py` and
`satisfiability.py` in three separate passages, none of them measured. Timed on
2026-09-01 it runs in **16.5 s**. It also **writes no figures** — it says so at
the top of the file — so "generates Figure 14.3" was wrong as well; what it
produces is the finite-size sweep behind Table 15.2. **Time a script before
calling it slow.**

### Sec. 14.7: the merge terminates on a placed ensemble, 2026-08-31

The author objected to Sec. 14.6, and was right. It measured the merge closure
on hyperbolic random graphs, found it affordable only when very sparse, and then
generalised: *"This is not special to hyperbolic graphs: clustering is exactly
what makes cliques share edges, so any ensemble with enough clustering to need
the repair tends to have enough to make it percolate."* That does not follow from
one ensemble, and **Ch. 2 already contradicted it** — `chygraphs.tex` says
triangles placed independently satisfy the treelike condition. The two
statements sat eleven chapters apart without meeting.

**The ensemble to test it on is Karrer & Newman**, *Phys. Rev. E* **82**, 066118
(2010), already in the bibliography and already cited in Ch. 2. Motifs are
placed by matching corners, so clustering is a free parameter and two motifs
share an edge only by coincidence.

**The measurement.** `shared_2+` is the wrong statistic here, because both
ensembles have `Theta(n)` intersecting pairs and dividing by that hides the
effect. Count the offending pairs instead:

| ensemble | C | pairs at n=1000 | at n=8000 | largest meta-complex |
|---|---|---|---|---|
| Karrer–Newman, `<k>`=3.0 | 0.25 | 14 | 11 | 7.0 → 6.0 |
| Karrer–Newman, `<k>`=8.9 | 0.16 | 340 | 351 | 23.7 → 9.0 |
| hyperbolic, `<k>`=2 | 0.28 | 511 | 5743 | 132 → 743 |
| hyperbolic, `<k>`=3 | 0.42 | 1563 | 15685 | 290 → 1556 |

**Flat against linear, at matched clustering (0.25 against 0.28), opposite
verdicts.** The mechanism is one line: two motifs meeting at a vertex share a
second one with probability `O(1/n)`, there are `O(n)` meeting pairs, so the
count is `O(1)` — hence `shared_2+ ~ 1/n`, measured slope −0.99.

**The best result is the diamond.** Karrer & Newman's own remedy for triangles
that share an edge is to place the *diamond* as a subgraph type — which is
Sec. 14.6's merge, done when the ensemble is written down. Run the closure blind
on such a graph and it returns the placed diamonds: 78.5% recovered exactly at
n = 1000, rising to 97.2% at n = 8000, the shortfall being the same O(1) residue
(~100 elements at every size). The merge inverts the modelling decision rather
than patching a broken mapping.

**The caveat, stated in the text and in Checks.** The O(1) constant grows with
density and at `<k>` ≈ 17 the closure percolates at n ≤ 4000 before collapsing to
19 atoms at 8000 and 13 at 16000. **A measurement at one n decides nothing; the
slope decides everything.** Also stated: the two ensembles are matched on
clustering and nothing else, and O(1) is an estimate supported by a flat
measurement, not a theorem.

**What the chapter now says**, replacing the claim above: the closure terminates
where the dense elements are *specified* (a modeller's motif set, or data that
records its groups — Ch. 3's route (a)) and fails where they are *emergent* (a
geometry producing overlaps at every scale). Clustering is compatible with
either. Sec. 14.6's verdict is scoped to its ensemble, and the Outlook's version
was stale for the same reason and is fixed.

**One bug caught before it reached the book.** The first `overlap_stats` counted
vertex *pairs covered by more than one clique* and divided by *clique* pairs —
two different denominators, and not Ch. 3's `shared_2+` at all. It now calls
`chygraph_statmech.region.overlap_profile`, the routine Ch. 3 measures real
networks with, so the two chapters count the same thing. The HRG numbers moved
by a factor of three when it was fixed.

New in `figs/merge.py`: `karrer_graph`, `diamond_graph`, `overlap_stats`,
`check_overlap_is_not_extensive`, `check_diamonds_are_recovered`,
`figure_motifs`. Runtime is up to about eight minutes; it was already the
slowest script in the book.

### Calculation-verification pass, 2026-09-01

Part IV (Chs. 14--16) had not been through one; it was drafted after the pass
below. Run: all thirteen scripts in `figs/`, `chygraph_statmech/tests` (347
passed), `chygraph/tests` (68 passed), and `latexmk` into a scratch `-outdir`.
Every quoted number in Chs. 3--13 still reproduces. **Twelve numbers in Part IV
did not, and are fixed**; all of them were transcription or rounding drift
against `probe/results/*.json`, none was a wrong calculation.

**Environment, because ten of the thirteen scripts fail without it.** They need
*both* `percolation/src` and `statmech/src` on `PYTHONPATH`. Since 2026-09-01
the scripts put both there themselves, resolved from the repository root rather
than from `$HOME`, so this is now a fallback rather than a requirement. Use

```sh
PYTHONPATH=~/av2atg/chygraph/percolation/src:~/av2atg/chygraph/statmech/src \
  ~/anaconda3/bin/python3 figs/<script>.py
```

**What was fixed.** In `gbp.tex`: "56 of 180" real runs chordal → `120`; the
residual-converged count `37` → `38` and the inconsistent count `7` → `8` in
three places; "margins as large as 2e-2" → `6e-2`; the grouped sound error
floor `7e-14` → `4e-14`; football ego 108's residual `5e-13` → `8e-13`;
"−ln 2 to eight decimal places" → `six` (the two runs sit `1.2e-7` and `2.1e-8`
from it). In `overlap.tex`: the real-class correlation `+0.83` → `+0.84`
(`0.8359`). In `metacomplex.tex`: sub-threshold core "to within 5e-4" → `1.5e-3`
(the `b=0.9` edges-only point at `n=1600` is `0.00151`); "`s≃1`, so `b` between
3 and 7" → `s = 1 or 1.5`, `b` between `3.5` and 7, which is what the
Conclusions already said and what `gbp_karrer.json` holds; "median four atoms"
→ a median of four complexes and four or five atoms (the atom median is `4.5`);
and "half of the twenty ... *are* the minimal ring" → half carry one cycle and
**eight** of those are the ring, which is what Fig. 16.5's caption says.

**The pair-correlation claim was dropped, then restored from the pipeline.**
`overlap.tex` reported that counting intersecting complex *pairs* rather than
doubled bonds gives `+0.53, +0.18, +0.68`. Nothing computed those:
`cavity_clique.solve` recorded `doubled_bonds` and no `shared_2plus`. It now
records both — `overlap_profile(cx)['shared_2plus']`, the same call and so the
same definition as `gbp_cliques`/`gbp_karrer`/`gbp_real`, which is the *fraction
of intersecting pairs meeting in two or more atoms* and not a raw count.
`probe/results/cavity_clique.json` is regenerated, and the correct triple is
**`+0.58, +0.52, +0.68`** against the bond measure's `+0.65, +0.63, +0.83`. The
bond measure wins on every class, which is the sentence's point; what does not
survive is the old claim that the pair measure gives "almost none" on the placed
ensemble, `+0.52` being a perfectly good signal. `figs/overlap.py:check_cavity`
now prints both correlations from the one file, so the comparison cannot drift
apart again.

**Sec. 16.3's transition is classical, and Ch. 16 now derives it.** The core in
Fig. 16.6 is the 2-core of the incidence structure, so its threshold is the
k-core transition at k = 2 (Pittel-Spencer-Wormald), equivalently Bauer and
Golinelli's alpha = 1. A calculation box gives the closed form —
`a = exp[-s(1-a) - t(1-a^2)]`, `mu = s(1-a) + t(1-a^2)`,
`P_C = 1 - e^-mu - mu e^-mu` — which yields `b = s + 2t` by differentiating at
`a = 1` and reproduces all three measured lines, threshold and amplitude, worst
discrepancy 0.007 against a standard error of 0.004.
`karrer_core_sweep.core_closed_form` evaluates it and `check_closed_form`
asserts the agreement, so Eqs. (16.2) and (16.3) are backed like every other
equation in Table 1. **Bauer and Golinelli's two removal rules are the book's
two cores** — leaf plus edge stops at mean degree one (Part IV), leaf plus
*neighbour* at `<k> = e` (Ch. 11) — which is why `P_C(Ising, .)` and
`P_C(VC, .)` differ on the same object. Sec. 14.2 says so now.

**Regenerating that file moved two of its 240 rows, and Secs. 14.3 and 14.5
now say so.** The running text carries the fact, where the error ranges are
quoted; the Checks section carries the detail. Do not push it back into the
checks alone.
Both are hospital-contact neighbourhoods at `beta J = 0.3`; their error in
`ln Z` goes from `116.6` to `58.1` and from `167.0` to `101.8`. The cavity
iteration has more than one fixed point there and the committed cache had found
the other one; two fresh regenerations agree with each other exactly, so this is
the environment rather than a seed. Nothing else in the file moves, no min,
median or max changes, and the only quoted number affected is the real-class
bond correlation — `0.8359` on the old file, `0.8340` on the new, which is why
it reads `+0.83` and not the `+0.84` this pass first corrected it to. **Do not
"fix" it back.** It is the same failure Sec. 15.4 documents for GBP on real
neighbourhoods, one method earlier.

**Two claims are cached-data-unbacked and were left alone.** `gbp.tex`'s
run-it-twice reproducibility check (86 runs to `2e-12`, largest disagreement
`14.19`, `−22.08` against `−7.88`): `probe/results/gbp_real_parts/` is identical
to `gbp_real.json`, so only one of the two runs survives. `−7.88` is in the data
(Hospital ego 64, βJ=0.3); `−22.08` is not. And `ising.tex`'s population-dynamics
`βJ_c = 0.164628` / `0.144343` — no script prints them.

**Verified and correct, do not "fix".** Table 14.1 and its `e^{1.3}≈3.7`; the
240-run split 98 chordal / 142 not; every Karrer--Newman number in Sec. 15.3
(43 sound of 58, gain `−0.07` to `1.70` orders, the `+1.018` against `−0.876`
instance); Table 15.1 to the digit; the ring table; Table 16.2 and the
`200` acyclic / `40` cyclic correspondence; `0.289 → 0.033`; the `b=1.6`
finite-size steps `+0.007, +0.012, +0.0007`; the treelike-validation `1.8e-15`;
`percolation.tex`'s `13.9` per cent, which is Ch. 9's degree-4 clustering shift
and not the percolation number; and all 109 `\code` names.

Also cleaned up: `covertriangles.aux`, left by the reverted Ch. 17 commit
(7826050), has no `.tex` behind it and can be deleted.

### Calculation-verification pass, 2026-08-31

Every number in the book that a script produces re-checked against that script,
and the closed forms behind them re-derived independently of both repositories.
Run: all thirteen scripts in `figs/`, `chygraph_statmech/tests` (347 passed),
`chygraph/tests` (68 passed), and `latexmk` into a scratch `-outdir`, which
before the edits below reproduced the committed `main.log` exactly.
**Every quoted number reproduced** — Table 3.2 exactly, the whole Potts
coexistence window, the AT line, Tables 12.2--12.4, the merge table, Ch. 13's
folds — except the five below, now fixed.

*Scope, so the next pass knows what this one did not do.* Values the book
**quotes from the literature** were checked for internal consistency and, where
a formula was given, re-derived — Gabrié's six, Zdeborová--Krzakala's two lines,
Mulet's `c_q`/`c_d`, Cirigliano's exact `p_c` from Eq. (5.18) — but a bare cited
constant with no formula behind it was not chased to its paper, and
Cirigliano's heterogeneous-mean-field `0.0635` is the one such number the book
sets its own result against. Nor were the bootstrap claims re-run: Ch. 11's
"roughly four standard deviations" at `tau = 2.9` rests on `probe` output that
was read, not regenerated.

**What was re-derived from scratch**, in sympy/mpmath, using neither repo:
`tau_3(q,v)` by brute-force Potts enumeration at `q = 2..6`; the Potts fold by
resultant, which factors as `v^2 (q-1)^2 (q-v)^2 (q+v)^2 (4q - v^2 - 4)` and so
gives Eq. (7.11); `u'(c)` for `c = 2..5` by Ising enumeration; both colouring
transmissions by symbolic perturbation of the interior up to `c = 6`; Gabrié's
six `l_stab` values; **all five of Ch. 8's critical temperatures including
5.2951, which no script prints**; `s_tri` with 5/3 (the published 3/2 fails);
`B_hyper = 4p/(1 + p<c>)` by solving the map near threshold; Eq. (12.6) at
`f = 0, 1` in exact arithmetic; `q_c = 0.403032`; Weigt--Hartmann; Eq. (5.18) at
`b = 3`; the `beta_1 = 0` saddle node (2.45541, 0.71533); the tricritical
coefficients; `T*` and `C_q` at `q = 2..8`. Also confirmed the unbacked "66% at
`<k> = 6`" in Sec. 14.6 (it is 66.0%) and Fig. 12.3's caption extrema (+0.064
and −0.400, so the caption's +0.06 is right and the script's docstring said
+0.07 — the docstring is what was wrong, and is fixed).

**1. The clustering effect decays as `1/n^2`, not `1/n`.** Three places:
`ising.tex`'s Table 9.2 caption and running text, and the Outlook's
"decays as `1/<k>`". The book's own table said so — a factor of 46 between
`n = 4` and `n = 20`, where `1/n` predicts 5 — and the stated reason was wrong
too. Solving the two conditions exactly, `(n-1)t = 1` for links and
`t^2 - (n-1)t + 1 = 0` for triangles (these reproduce the table's `T_c` column
to the digit), gives `dT_c/T_c = -1/(n-1)^2 + O(n^-3)`: −4.0, −2.0, −1.2, −0.28%
at `n = 6, 8, 10, 20` against the table's −4.3, −2.1, −1.3, −0.3. The mechanism
is that the lost branch and the transmission gain are **each** `O(1/n)` and
cancel at that order; the text now says so. It also now says that the
*absolute* shift does die as `1/n`, since `T_c ~ n-1` — which is probably where
the original claim came from, and separating the two is what stops it returning.

**2. Ch. 3's clique-number check overstated its own agreement.** "Consistently a
little below theory" for 0.37/0.20/0.07 against 0.45/0.25/0.05: the third is
*above*. Rewritten to report the two heavy-tailed rows as below and the third as
above, and to say why that row decides nothing — the predicted exponent (0.05)
is smaller than the measurement's spread across mean degrees (0.06--0.09). The
check still earns its place through the trend. **This is the chapter's only
independent check that the measurement reads the right quantity**, so it is the
one sentence there that has to be exactly right.

**3. Sec. 9.6 disagreed with itself by a digit.** "`m` jumping from 0.548 to
0.833" against the displayed `Delta m = 0.286`; `ising.py` prints 0.834.

**4. Table 3.3's `tau = 4.5` ratio was not a ratio.** 2.93 is the raw HRG
`sbar`; the column is `hrg/config` and the other three entries are ratios. Now
2.95.

**5. Ch. 1 called vertex cover "the third model"**, two lines after announcing
three and one subsection after the first. Now "the second".

**Checked and cleared, do not "fix" either side.** Ch. 12 cites Krzakala
*et al.* Eqs. (20) and (21) where the entry below says (16) and (20). Both are
right: the book cites the published PRE 70, 046705, whose numbering I confirmed
from the PDF, and the entry below refers to the arXiv source. Likewise
`giant.tex:597`'s "the map reports `S = 1` where simulation gives 1/2" does have
its counterexample in Sec. 5.6, and it is correct — the
`Phi = x_|^2/2 + x_tri^2/2` split leaves the link piece 2-regular, which is
exactly the `C = 0` case Sec. 5.5 flags.

**Two production defects surfaced by the re-measurement**, neither of them a
calculation and neither changing a number. The software chapter's map table was
printing as **Table 15.2**; that is fixed in `main.tex` and written up in the
paragraph on `software.tex` above. The overfull hbox in `cover.tex:239--249`
is **still open** and is recorded in Status.

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

### Set-valued SP for triangles: gate passed, application did not, 2026-08-30

Attempted after Sec. 12.5, to get the *colourability* threshold of a
triangle-clustered graph — which nobody has, Mulet being graphs and Gabrié
hypergraphs. See `../statmech/probe/SETVALUED_SP.md` and `../statmech/probe/setvalued_sp.py`.

**Validated at cardinality two.** Both a two-term and a three-term complexity
form reproduce Mulet's `c_q` — 4.684, 8.863, 13.636 against 4.69, 8.90, 13.69 —
and agree with the scalar code committed in `figs/colouring.py` (4.683, 8.901,
13.660, which still reproduce exactly). So the set-valued survey, the
intersection update and the complexity are all sound where they can be checked.
**The tolerance originally claimed here, "better than 0.5%", was one run's
luck**: re-run 2026-08-31 at population 4000 the same `threshold(q, 2, ...)`
returns 4.691, 8.825, 13.620, that is +0.02 / −0.84 / −0.51%. The gate holds to
about one per cent, which is all it needs to; nothing in the book depends on it.
Re-running it needs the bracket trap below — `lo` must sit **above** the branch
onset or the assertion refuses it, and `lo = 2.0` does.

**Produced nothing at cardinality three.** `Sigma` is negative from the moment
the branch appears and never crosses. Two signs it is the implementation and not
the physics: magnitudes of 0.1–0.8 where the validated `c = 2` case sits at
0.01, and an ordering unlike the graph's (branch onset at degree ~3 against
Sec. 12.5's RS line at 6, where for graphs onset/crossing/RS line sit within
13%). The `c = 3` interior is the one piece the `c = 2` gate cannot check.

**Sec. 12.7's statement is unaffected and remains the book's position:** the
proper rule does not carry over and no substitution repairs it. Nothing
committed depends on this work.

**Two things in that file are worth reading before any similar attempt.** Five
bugs are recorded, all caught by gates rather than by reading code, including a
diagnosis that the next run disproved within a minute and a `c = 2`-only
complexity form reused where it does not hold. And the gate note: **use
`Sigma = 0`, not the branch onset.** A crossing can be bisected; a fold cannot,
and both implementations locate `c_d` to only a few per cent while disagreeing
with each other by more — which briefly looked like a 5% conflict between them
and was not.

### Sec. 12.9: extending Krzakala et al. to chygraphs, 2026-08-30

The author read Ch. 12 and said, correctly, that it had no new result beyond
"for a graph with triangles the chygraph mapping is the way to go", then asked
whether Krzakala's approach could locate the transition. It can, and it does.

**The apparatus is theirs, verified against the arXiv source** (`cond-mat/0403725`,
`coloring.tex`), not the PDF text. With `eta = e/q`:

- `sp_update` **is** their Eq. (16), `eq_self_cons_q` — numerator, denominator,
  term for term.
- `sp_complexity` **is** their Eq. (20), `eq:sigma` — the site sum over the full
  degree `p_d`, and the link term `-(c/2) ln(1 - q eta_1 eta_2)`.

So Sec. 12.7 does not merely agree with the published calculation; it *is* the
published calculation with a chygraph in place of a graph. Worth stating in the
book, and now stated.

**The result.** `Sigma(y=inf)` counts ZERO-ENERGY clusters, so negative means
none exist and the colourable phase is bounded by where they run out. On a graph
that is a crossing; on triangles the branch arrives already negative, so the
phase ends at the lift-off.

| q | graph `c_q` | triangles | change |
|---|---|---|---|
| 3 | 4.69 | **3.078** | −34.4% |
| 4 | 8.90 | **6.763** | −24.0% |

(Re-run 2026-08-31; this entry previously recorded 3.083 / 6.771 and
−34.3% / −23.9%, from an earlier run of the same code. The difference is
population scatter at the third digit and the conclusion is untouched.
**Table 12.4 of the book carries the values above** — keep the two in step.)

**And it points the opposite way from Sec. 12.5.** The RS stability line goes
*up* with triangles (5 → 6, 10 → 11) while `c_q` goes *down*. On a graph the two
sit within 6% and one proxies the other; on triangles `c_q` is half the RS line.
**Anyone using the RS threshold as a proxy for colourability gets the SIGN of
the clustering effect wrong** — concludes clustering helps, when it costs a
third of the range.

**Three figures**, all added on the author's reading that the section asserted a
different technique without showing it. Fig. 12.2 draws what a survey *is*: one
cluster and one message, against many clusters each sending their own.
Fig. 12.3 draws `Sigma` against degree — the graph's arc above zero crossing at
`c_q`, the triangle branch starting below zero with no arc. Fig. 12.4 draws the
phase structure: a graph passes through one cluster, then exponentially many,
then none; triangles skip the shattered phase.

**The section was rewritten twice, both times on the author's reading, and the
lesson is worth keeping.** First it did not say this *extends* Krzakala et al.
to chygraphs — only that the apparatus was theirs, which understates it: their
Eqs. (16) and (20) are for graphs and have nowhere to put a complex. Second, it
still read as a recovery from confusion, because the original Sec. 12.9 was the
failed reading written down. The two sections are now one, opening with what was
done. **The discovery order is not the exposition order; do not restore the
narration.**

**What is still not done:** the energy above `c_q`. That is Krzakala's `e(y)` at
finite `y`, and four attempts at it failed (see the note below and
`../statmech/probe/FINITE_Y.md`). It answers *how badly* uncolourable, not *where* — the
threshold needed only the `y = inf` endpoint, which is gated.

### Sec. 12.9, the proper rule closes the window, 2026-08-30

The author asked for survey propagation on a graph with triangles and, when I
kept pointing at a failed attempt instead of doing it, told me to start over.
The restart found the result. **Three of my earlier conclusions were wrong**,
each to a check I should have run first:

- "`q = c` is degenerate, hence the failure" — no: `q = 4, c = 3` behaves the same.
- "the `c = 3` branch appears by a fold, unlike the graph" — no: **both** appear
  discontinuously, 0 → 0.74 and 0 → 0.68.
- "the `Sigma` magnitudes are too large, so it is a bug" — that compared `Sigma`
  near a crossing with `Sigma` far from one. Not a comparison.

**What is true.** A threshold needs an interval where a non-trivial survey
exists *and* `Sigma` is still positive. The graph has one; a triangle network
does not.

| | q=3 graph | q=3 triangles | q=4 graph | q=4 triangles |
|---|---|---|---|---|
| branch appears at (degree) | 4.45 | 3.20 | 8.20 | 6.80 |
| `Sigma` there | +0.026 | **−0.160** | +0.060 | **−0.216** |
| window | to 4.69 | none | to 8.90 | none |

**The gate that makes this physics and not a bug.** The set-valued apparatus is
checked twice against published numbers: at `c = 2` against Mulet's `c_q`
(three colour counts), and at `c = 3` against Gabrié's `l_col` with **only the
constraint rule swapped** and every `c = 3`-specific piece held fixed — the
interior tensor, the SDR tensor, the `kappa/c` count, the three-term `Sigma`.
The hypergraph rule at cardinality three *has* a window and returns 26.92 and
63.3. So the closure is the proper constraint, not the cardinality.

**State it as a negative result, not a threshold.** `m = 0` cannot locate the
colourability threshold of a triangle-clustered graph; it does not follow that
there is none. Survey propagation counts clusters without regard to size, which
is wrong once they have stopped being equally worth counting before any exist.
The repair is `Sigma(m)` with `m != 0` and is not done.

**Second occurrence of the same shape.** Sec. 13.9's nested clauses fail
identically — branch at finite `Sigma`, no crossing (`../statmech/probe/NESTED_SP.md`).
Two problems sharing nothing but a cardinality above two. **If a third turns
up, that is a thread and belongs in the Outlook.**

**Two traps, both of which cost real time.** The branch appears
*discontinuously* in every case, graph included, so a coarse scan steps over the
whole positive window — the graph's is about 8% wide in degree, and my `c = 3`
grid was 30% wide. And `Sigma` is exactly `0.0` on the trivial branch, so a
bracket below the onset is refused rather than bracketing. Both are in the
module header of `figs/colouring.py`.

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
It was attempted and is **unresolved** — see `../statmech/probe/NESTED_SP.md` and
`../statmech/probe/nested_sp.py`. Two of three parts work and are worth keeping: the
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
   `../statmech/TODO.md` item 3. `software.tex` prints the URL as if it resolves and the
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
- `~/av2atg/chygraph/percolation/manuscript_3/manuscript.tex` — the percolation extension.
  **Will not be submitted** (author, 2026-08-30): the book is the first place
  Ch. 5's critical amplitude and moment hierarchy are reported, and Ch. 5's
  opening says so. Do not add a citation to it.
- the manuscript this book supersedes (`main.tex`, `supplement.tex`, at the
  root of the then-separate `chygraph_statmech` repository) —
  **deleted**, recoverable from git at commit `5ebd892`, the last commit before
  the deletion. Everything in it is in the book; the section numbers named in
  `main.tex`'s `\include` comments refer to that retired file.
- `../statmech/src/statmech/`, `../statmech/examples/`, `../statmech/probe/` — every number in the book comes out of these

## Building

```sh
latexmk -pdf -interaction=nonstopmode main.tex     # -> main.pdf
```

Verify a build with

```sh
grep -c '^!' main.log                    # errors, must be 0
grep -i 'undefined' main.log             # citations and refs, must be empty
grep -i 'multiply defined' main.log      # must be empty; NOT caught by the line above
grep -c 'Overfull \\hbox' main.log       # 3 known, listed in Status; a 4th is new
grep -o 'Output written.*' main.log      # page count
```

Underfull vboxes are not a defect here: all 33 are `while \output is active`,
which is page-breaking around floats. An **overfull hbox** is, since it puts ink
in the margin; the three that survive are recorded in Status.

### Two checks the build cannot make

A clean build is not evidence of either of these, because both fail silently.
Run them after any pass that moves a figure or renumbers a chapter.

**Every included figure is in the repository.** `.gitignore` carries `*.pdf`, so
a figure that was never force-added is invisible to `git status` and present
only on the machine that drew it — the book builds there and nowhere else. This
found two, in Sept 2026, after an earlier pass had found six by eye and stopped:

```sh
cd book
for f in $(grep -rho 'fig-[a-z-]*\.pdf' *.tex | sort -u); do
    git ls-files --error-unmatch "figs/$f" >/dev/null 2>&1 ||
        echo "UNTRACKED (but included): figs/$f"
done
```

Reverse it to find orphans — figures still generated but included nowhere — by
looping over `figs/fig-*.pdf` and grepping the `.tex` for each. Both directions
should print nothing.

**Every section reference in the code still resolves.** The scripts cite the
book by number (`Sec. 15.4`, `Chapter 16`) rather than by label, so nothing
breaks when a chapter is renumbered; the comments just start lying. Part IV
splitting into three chapters left eleven stale references behind. This
compares them against the numbers LaTeX actually assigned:

```sh
cd book && python3 - <<'EOF'
import re, glob, pathlib
sec, cha = {}, {}
for f in glob.glob('*.aux'):
    txt = open(f).read()
    for m in re.finditer(r'\\newlabel\{[^}]*\}\{\{(\d+(?:\.\d+)*)\}\{\d+\}\{(.*?)\}\{(section|chapter|subsection)\.', txt):
        (cha if m.group(3) == 'chapter' else sec)[m.group(1)] = m.group(2)[:52]
    for m in re.finditer(r'\\contentsline \{(subsection|section)\}\{\\numberline \{(\d+(?:\.\d+)*)\}(.*?)\}\{', txt):
        sec.setdefault(m.group(2), m.group(3)[:52])
for r in ('figs', '../statmech/probe', '../statmech/src', '../statmech/tests',
          '../statmech/examples', '../percolation/src', '../percolation/tests'):
    for p in sorted(pathlib.Path(r).rglob('*.py')):
        if '__pycache__' in str(p):
            continue
        for i, line in enumerate(p.read_text(errors='ignore').split('\n'), 1):
            for m in re.finditer(r'\b(?:Sec\.?|Section)s?\s+(\d+\.\d+(?:\.\d+)?)', line):
                if m.group(1) not in sec:
                    print(f'{p}:{i}  Sec. {m.group(1)}  {line.strip()[:70]}')
            for m in re.finditer(r'\b(?:Ch\.?|Chapter)s?\s+(\d+)\b', line):
                if m.group(1) not in cha:
                    print(f'{p}:{i}  Ch. {m.group(1)}  {line.strip()[:70]}')
EOF
```

It should print nothing. Two caveats learned the hard way: index **subsections**
as well as sections, or `Sec. 5.6.2` reads as broken, and match the title
non-greedily, or `Sec. 9.4` does — its title contains `\texorpdfstring`. And
resolving is weaker than being right: a reference can point at a section that
exists and is the wrong one, which is how `Sec. 14.2's claim` survived in
`overlap.py` after 14.2 became a different section. The audit narrows the
reading; it does not do it.

Build into a scratch directory --- `latexmk -pdf -interaction=nonstopmode
-outdir=/tmp/bookbuild main.tex` --- when the point is to verify rather than to
refresh `main.pdf`; the output is identical and the working tree is left alone.

The counts in Status are re-measured, not maintained by hand. Figures, tables
are `grep -c '\\begin{figure}'` and friends over `*.tex`;
numbered equations are `\begin{equation}` plus the numbered rows of `align`;
`\ref`/`\eqref` resolution and `\code`-name resolution are worth scripting
against both repositories, remembering that `joint.JointGiantComponent` is an
alias rather than a `def` and that `\allowbreak` has to be stripped out of a
`\code` argument before it is looked up.

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
`percolation` package from `../percolation/src`;
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
`../statmech/probe/results/ising_mc.log` rather than recomputed — that run is the one
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
hyperbolic-graph points are read from `../statmech/probe/results/prediction4.csv` and the
simulation check from `../statmech/probe/results/validate_core.txt`; neither is
recomputed by the book.
`figs/cover.py` also writes **Table 11.1** (`table_real_core`) from
`../statmech/probe/real_core.json`: leaf removal on sixteen real networks --
Ch. 3's ten and Ch. 14's six, from the same netzschleuder cache -- against
Eq. (11.4) on each graph's own clique ensemble and against a degree-matched
control averaged over twenty rewirings. `statmech/probe/real_core.py` produces
the JSON and takes about two minutes; it needs `numba`, through
`computational_complexity/code/leafremoval.py`.
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
the book.
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
which makes this the second slowest script in the book.
`figs/overlap.py` plots **all of Part IV** — it is the one script serving three
chapters. Figures 14.2, 14.3, 15.2–15.4 and 16.3–16.6, from the cached probe
outputs (`gbp_cliques.json`, `gbp_real.py`, `merge_lnz.py`, `core_fraction.py`,
`karrer_core_sweep.py`); Figures 14.1, 15.1, 16.1 and 16.2 are TikZ in the
chapters themselves. It also runs the checks: the Möbius and Bethe counting
numbers on two overlapping triangles, including the factor-coverage test that is
the whole point (shared bond covered once against twice); the two-triangle
table, recomputed here since it is a four-spin enumeration; the 60-instance
summary from `../statmech/probe/results/gbp_cliques.json` with the convergence threshold
stated; and the clique-ensemble paired ratio read from
`../statmech/probe/results/analysis.txt`.
`figs/merge.py` **writes no figures** — it says so at the top of the file, and
it used to write `fig-merge` and `fig-motifs`. What it produces now are numbers:
the finite-size sweep behind **Table 15.2** (`check_placed_finite_size`), the
merge closure on the six real networks, and the two-triangle and diamond checks.
It is **not** slow: timed on 2026-09-01 it runs in **16.5 s**. Three separate
passages of this README used to call it the slowest script in the book, ranked
against each other inconsistently; none of the three had been measured. It needs
`PYTHONPATH` to carry both `src` trees, or it fails at import.

## Conventions

- Drawing vocabulary is fixed once, in `main.tex`: `nd` (node), `ndf` (filled
  node), `hub`, `cx`/`cxb` (a complex drawn around a set of nodes), `lnk`,
  `msg` (a message along a directed inclusion), `lb`/`tb` (labels), and `\cxg{c}`
  for an inline glyph of a cardinality-$c$ complex. Every figure is built from
  these, so a reader who learns the first schematic can read the last one
  without its caption.
- Notation matches the papers (`\ave{}`, `\kbar`, `\sbar`), so a reader can move
  between book and paper without retranslating.
- **No `calculation` boxes.** Reversed 2026-09-01: the boxes were an obstruction
  to reading rather than an aid, because a derivation the reader is told to skip
  is a derivation they cannot follow when they want it. Derivations run in
  continuous text, where they are needed, and **each ends by naming the routine
  that executes it** --- `cover.CliqueCover.solve` and the like --- so the path
  from the algebra to the code is in the sentence rather than only in
  `software.tex`'s Table 1. All 25 boxes are gone, and the `calculation`
  environment and the `tcolorbox` dependency went with them, so a new box will
  not compile.
- **Long `\code` names carry `\allowbreak` at their separators.** `\code` is a
  bare `\texttt` with no break points, and a name like
  `statmech/probe/karrer_core_sweep.core_closed_form` overflows the 4.05in
  measure by 232pt without them. The rule is applied everywhere, including
  `software.tex`'s table, and is idempotent: insert after `.`, `/` and `\_`
  where one is not already there.
- Overlapping complexes are drawn as two translucent boxes of the same grey, so
  the shared atoms show as a darker patch. Established in Fig. 2.7 and used
  wherever overlap matters.
- **Section and subsection headings are ragged right**, patched in `main.tex`
  from `book.cls`'s own definitions with nothing changed but the addition of
  `\raggedright`. The class sets them justified, which on this 4.05in measure
  makes a long title run into the margin instead of breaking — renaming Sec. 1.6
  to "Replica symmetry and symmetry breaking" overflowed by 15.5pt. Chapter
  heads were already ragged. Do not revert this to get "tidier" headings; a
  single-line heading looks identical either way.
- Figures must be checked at the 5×8 trim, not just for compilation: more than
  three panels side by side is illegible on this page width. Stack into rows.
