# TODO — open items

Work that is planned but not done. Items here are *not* covered by the tests;
what is done lives in the README's Results sections.

---

## 1. Test prediction 4 on hyperbolic random graphs — **RESOLVED, confirmed**

**The claim.** A chygraph whose complexes are the HRG's clustered motifs has a
non-empty core where the degree-matched configuration model does not.

**Answer: yes — with the exponent claim stated more carefully than it was.**
`probe/prediction4.py`, `n = 2x10^5`, maximal cliques as complexes, ensemble
measured from each graph so nothing is fitted. `theta` throughout, `beta` being
reserved for the inverse temperature (referee minor 2).

| tau | theta measured | theta chygraph | theta control |
|---:|---:|---:|---:|
| 2.5 | 1.590 ± 0.061 | 1.553 ± 0.035 | no core |
| 2.9 | 1.635 ± 0.008 | 1.573 ± 0.011 | no core |
| 2.1 | 1.493 ± 0.060 | 1.571 ± 0.013 | no core |

The chygraph produces a power law `core ~ kbar^theta` of about the right
exponent while the degree-matched control has **no core at all** over the fit
range. What it does *not* do is track the exponent's dependence on the tail:
bootstrapped over seeds the predicted exponent has spread 0.020, consistent with
no `tau` dependence, while the measured one spreads by 0.143 and rises
monotonically, and at `tau = 2.9` the two differ by roughly four standard
deviations. The referee was right that the earlier "including the exponent"
reading overstated it.

The sharp check in this item is the **rewired control**, not the exponent: the
same measured clique ensemble arranged treelike, which the prediction matches to
`1e-3`. That separates whether the map is right from whether treating complexes
as independent is, and it says the map is right on its own terms.

Magnitude, tau = 2.9: the ratio prediction/measurement is 0.97-0.99 at
`kbar <= 0.1`, falls to 0.77 near `kbar ~ 1.5-3`, and returns to 0.96 at
`kbar = 6`. That tracks clique overlap, which is negligible when the graph is
sparse and saturates when almost everything is core.

**The cleanest control.** At tau = 2.9, `kbar = 6` the configuration model does
have a core, and the chygraph predicts it to 2% (0.633 against 0.621). Same map,
same pipeline, and the control's cliques barely overlap (`shared_2plus` 0.00-0.01
against 0.08-0.32 for the HRG). So the 23% deficit on the HRG at the same
density is overlap and not the map.

**tau = 2.1 behaves differently and should not be used**: the chygraph
*over*-predicts there, by up to 33%. That is consistent with step 1, which found
the clique ensemble has no `n`-independent limit below tau = 2.5, so the measured
ensemble is not converged and the prediction built on it is not either.

**A claim of mine that was wrong.** This item previously said "a chygraph of
independent cliques predicts `1 - Phi(0)`, which is not a power law in `kbar`".
Both halves are false. `1 - Phi(0)` is the prediction only when *every* layer has
cardinality >= 3, so that no vertex ever has degree 1 and leaf removal never
fires. Real HRG clique ensembles contain many cardinality-2 cliques, leaf removal
does fire, and the fixed point is non-trivial — and it does give a power law,
with very nearly the right exponent.

**What is left.** Only the magnitude, and it is attributed rather than open:
generalised belief propagation on the region graph of Sec. VIII, which is the
standing item in the paper's conclusions.

---

## 2. Generalised belief propagation at the level of the *ensemble*

**Done, at the level of one instance.** `gbp.py` passes parent-to-child messages
on the region graph of an explicit complex list. It contains belief propagation
exactly (checked against an independent implementation, loopy error included),
it is exact on the two-triangle example — recovering the whole of
what the book's Sec. 14.4 reports as left on the table — and where it is approximate its
fixed point still satisfies `sum_{x_P \ x_R} b_P = b_R`.

**Open: the ensemble.** Every other module here carries a message per chy-degree
*class*, not per object. The lift needs a message indexed by region **type**,
over a region graph that is itself a random object: what is the distribution of
intersection patterns among the maximal cliques of an HRG, and what does the
parent-to-child update become when averaged over it? That is what would close
the 17–35% overlap deficit of Sec. VII, and it is the standing item in the
manuscript's conclusions.

Two facts to build on, both measured in `probe/gbp_cliques.py` over 60
maximal-clique region graphs of HRGs at `n = 14, 18, 20`. The instance
calculation is exact exactly when the clique structure is chordal — all 40
chordal runs to `3e-12` — so the ensemble question is really about the
distribution of chordality-violating motifs. And of the 20 non-chordal runs only
**9 converge**, even at damping 0.999; where they do the error is `1.4e-9` to
`5.4e-3` against `0.48–1.1` for the static Möbius counting, and where they do
not the residual says so. An ensemble treatment inherits that instability, and
damping will not be enough — it needs a schedule or a different message
parameterisation.

## 3. Repo visibility

The book's *The software* chapter prints
`https://github.com/av2atgh/chygraph` and
`https://github.com/av2atgh/chygraph_statmech` as if both resolve, and tells the
reader that every number in the book comes out of them. The second is private.
It must be public before the book goes anywhere.

(This was the manuscript's acknowledgments requirement; the manuscript has been
retired, and the book inherits the obligation in a stronger form, since it gives
the reader an equation-to-routine table to follow.)

## 4. Boolean networks — assessed and **not pursued**, 2026-08-28

**Decision: skipped.** Kept as a record of what was checked, so the question
does not have to be reopened from scratch. The findings below stand; what is
missing is the directed extension of item 7, which everything else waits on.

**Question.** Would the chygraph construction help with random Boolean networks
(Kauffman): the order–chaos transition, canalyzing rules, attractor scaling?

**Verdict: yes for three things, one of which is new; no for the thing the field
most wants.** Details below, with what was actually checked.

### What was verified

`influence = u'`. The per-input influence of a Boolean function — the
probability that flipping that input flips the output — plays exactly the role
of `u'` in Eq. (8.5). Measured by enumeration: a random function of bias `b` has
per-input influence `2b(1-b)` (0.5081 vs 0.5000 at `b = 1/2`, 0.4217 vs 0.4200
at `b = 0.3`, K = 3..5). The branching condition `K <influence> = 1` is then
Bastolla & Parisi's critical line `p_c = 1/K` and the Derrida–Pomeau annealed
result `2b(1-b)K = 1`.

Total sensitivity, measured:

| K | random `f` | nested canalyzing |
|---:|---:|---:|
| 2 | 1.015 | 0.797 |
| 4 | 1.985 | 0.951 |
| 6 | 2.992 | 0.984 |
| 8 | 3.979 | 1.030 |

Random functions give `K/2`, so chaotic for `K > 2`. Nested canalyzing functions
**saturate near 1** — marginally critical whatever `K`. (Measured over one
random-NCF ensemble; consistent with Kauffman *et al.* 2004.)

### Where the construction helps

1. **The order–chaos transition is already a chygraph calculation.** Substitute
   the influence for `u'` and `det(I - B) = 0` is the multi-layer generalisation
   of `K<s> = 1`. **But the one-layer case is Derrida–Pomeau (1986) and
   multi-type versions exist (Pomerance *et al.* 2009).** This is
   systematisation, not novelty, and on its own would not justify a chapter.

2. **It is the successor to Markert, Baas, Levine & Vázquez (2010).** That paper
   already introduced higher-order Boolean networks with an explicit distinction
   between cell components and got stability from an eigenvalue problem, using
   Baas's hyperstructures. Chygraphs are the finished formalism for exactly
   that: arbitrary depth, the `L x L` determinant of Eq. (8.9), correlated
   layers (Sec. 5.5), mixed cardinalities (Eq. 8.8). **This is the most
   defensible claim**, because the prior work is the same author reaching for
   the same thing without the machinery.

3. **The new part: motif-rich regulation breaks the annealed approximation.**
   Derrida–Pomeau assumes the damage-spreading graph is locally treelike. Real
   gene regulatory networks are motif-rich — feed-forward loops, bi-fans — which
   is exactly what treelikeness forbids. Promoting a motif to a complex and
   summing its interior exactly is the book's central move, and it has no
   counterpart in the Boolean-network literature that I can find. **This is
   where a chapter would earn its place.**

4. **Canalyzing rules lift the `2^c` ceiling.** Sec. 8.7 names the interior cost
   `2^c` as the formalism's one hard limit. A nested canalyzing function is a
   *chain of conditions* — a nested object, which is what a chygraph is — and
   its interior can be evaluated in `O(K)` rather than `2^K`. So the one hard
   limit does not bite for the function class biology actually uses. Worth
   saying out loud.

5. **The frozen core is a lead worth following.** Bastolla & Parisi 1998b:
   "the phase transition in random boolean networks can also be described as a
   percolation transition", with relevant elements grouped into functionally
   independent modules. Relevant elements are found by iteratively decimating
   constant-output nodes — structurally leaf removal — and Ch. 11's three-state
   message (`L`/`D`/`C`) is built for that class of algorithm. Ch. 11's result
   that the core-free branch dies above cardinality two would be the thing to
   test for group regulation.

### Where it does not help

6. **Attractor number and length scaling.** This is what the field most wants
   and the book has no dynamics at all — the Outlook lists that as not done.
   Attractor counts are not fixed points of a message equation. The one indirect
   route is Bastolla's (relevant elements → modules → attractors): if the
   chygraph can predict the module-size distribution, which is a percolation
   calculation, it might reach the `sqrt(N)` scaling of relevant elements at
   criticality. The superpolynomial attractor count of Samuelsson & Troein needs
   more than a branching matrix. **Speculative; do not promise it.**

7. **The formalism would need a directed extension, and does not have one.**
   Boolean networks are directed: damage flows along regulation, and the
   branching ratio is *out*-degree times influence. The book's inclusions carry
   a direction of arrival (up/down) which is *not* causal direction. In-degree
   and out-degree distributions would both be needed, and the four moment
   matrices of Sec. 2.8 would have to be split. This is real work, not a
   substitution, and it is the main technical risk.

8. Synchronous update, and the annealed approximation itself, are assumptions
   the chygraph does not improve.

### If it goes ahead

Order of work: (i) directed chy-degrees, since everything else waits on it;
(ii) reproduce `2b(1-b)K = 1` and Kauffman 2004's canalyzing stability from
`det(I - B) = 0`, as the benchmark; (iii) the motif calculation of item 3, which
is the actual contribution; (iv) the frozen core of item 5. References are under
`~/Downloads/chygraph_references/` (`boolean_networks` in the name).
