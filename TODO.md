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
