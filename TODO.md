# TODO — open items

Work that is planned but not done. Items here are *not* covered by the tests;
what is done lives in the README's Results sections.

---

## 1. Test prediction 4 on hyperbolic random graphs — **RESOLVED, confirmed**

**The claim.** A chygraph whose complexes are the HRG's clustered motifs has a
non-empty core where the degree-matched configuration model does not.

**Answer: yes, and the shape as well as the sign.** `probe/prediction4.py`,
`n = 2x10^5`, maximal cliques as complexes, ensemble measured from each graph so
nothing is fitted.

| tau | beta measured | beta chygraph | beta control |
|---:|---:|---:|---:|
| 2.5 | 1.584 | 1.551 | no core |
| 2.9 | 1.635 | 1.573 | no core |
| 2.1 | 1.491 | 1.571 | no core |

The chygraph reproduces the power law `core ~ kbar^beta` that
`computational_complexity` measured (1.50 and 1.63 at tau = 2.5, 2.9), to within
0.03 and 0.06 in the exponent, while the degree-matched control has **no core at
all** over the whole fit range. That was the open half of the item: a nonzero
prediction was never in doubt after WP5, the exponent was.

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

## 2. WP4–WP6

Distributional messages, complexes as regions, and the Bethe free energy.
Stated in the README's Work packages; nothing to add here yet.
