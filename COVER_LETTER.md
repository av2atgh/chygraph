# Response to the referee

I am grateful for a report that engaged with the calculations rather than the
prose. One of its objections is correct and overturns a headline result; I set
that out first.

## The essential item: §4.1

**The referee is right, and the claim had the wrong sign.** I reproduced the
regular comparison by hand and independently in code, and my numbers agree with
the referee's to every digit given: for an `n`-regular graph against `n/2`
triangles per vertex — identical `p_d = δ_{d,n}`, identical `e_dd'`, identical
edge count — the transition is at `t_c = 1/(n-1)` for links and at
`(n-2) u'_T = 1` for triangles, giving `T_c` **lower** by 13.9% at `n = 4`.

The diagnosis is the referee's. Table II compared a Poisson link layer against a
Poisson triangle layer. In the latter the node degree is `d = 2X` with
`X ~ Poisson(n/2)`, so `⟨d²⟩ = n² + 2n` and the excess degree that Eq. (6)
actually uses is `κ̄ = n+1` rather than `n`. The +14.5% was a degree effect
mislabelled as clustering, and "identical degree distribution and assortativity"
was false: the construction forces even degrees with twice the variance.

Changes: Table II is recomputed against both nulls, the regular one (which is
unambiguous) and a Poisson null carrying the triangle ensemble's own excess
degree. Sec. V B now derives the mechanism the referee identifies — the
traversed complex removes a branch, so `u'_T` would have to reach 1/2 and
reaches 3/7 — and states that the Poisson comparison is retained *because* it
misleads. Abstract and conclusions corrected.

**Monte Carlo (§4.8).** Wolff cluster updates on a 4-regular random graph and on
a network of two triangles per vertex, both verified to have every vertex at
degree exactly 4. Binder crossings at n = 3000 and 12000:

| | Monte Carlo | Eq. (6) |
|---|---|---|
| 4-regular graph | 2.894 | 2.8854 |
| two triangles per vertex | 2.482 | 2.4853 |
| change | **−14.2%** | **−13.9%** |

No cavity equation enters the simulation. Sec. V C is retitled and now separates
implementation checks from this one physical test.

## Other major points

**§4.3 — a real bug, fixed.** `u'` for the simplicial interaction was the
derivative with respect to a field common to all `q-1` members and so already
contained `(q-1)`; I confirmed the ratio is exactly `q-1` by independent-field
enumeration. It is now per-neighbour, matching Eq. (12), with the multiplicity
supplied once by Eq. (6). A regression test pins the relation between the two
conventions so the trap cannot reopen. No numerical result changes.

**§4.4.** Stated: a layer is single-cardinality by construction throughout, so
`c_m - 1` is correct as used; the size-biased `s̄_m = ⟨c²⟩/⟨c⟩ - 1` is given for
the general case, and the clique ensembles of Sec. VIII are split into one layer
per clique size.

**§4.5.** The ensemble-level map is now displayed, Eqs. (3)–(4): the action of a
generating function on a measure, and the convolution form of the chy-degree
step, with the reductions to the percolation map and to Eq. (5) stated.

**§4.7 — the referee understated this, and the revision says so.** Bootstrapping
over seeds gives a predicted exponent of 1.571 ± 0.013, 1.553 ± 0.035,
1.573 ± 0.011 at τ = 2.1, 2.5, 2.9: spread 0.020, consistent with no τ
dependence. The measured exponent spreads by 0.143 and rises monotonically, and
at τ = 2.9 the two differ by roughly four standard deviations. Table IV now
carries uncertainties, and the text states that the comparison shows the
chygraph produces a power law of about the right exponent and *not* that it
tracks the exponent's dependence on the tail — together with the referee's point
that a positive prediction follows from the algebra once one triangle enters, so
the qualitative separation carries less weight than the previous draft implied.

**§4.9 — demonstrated rather than asserted.** `ln Z` enumerated exactly for the
two-triangle example: the Möbius counting is closer than the Bethe counting by a
factor of 3 to 80, and the two err in opposite directions, Bethe overestimating
because it subtracts the shared vertices independently. New table and test.

**§4.2.** Newman (2009), Karrer & Newman (2010), Yoon *et al.* (2011),
Cantwell & Newman (2019), Liu *et al.* (2012) and Coutinho *et al.* (2020) are
added, with a paragraph in the introduction positioning the work against them:
the structural idea is theirs, and what is added here is that Eq. (3) leaves the
interior Hamiltonian free, so one determinant serves four models. Against
Liu *et al.* the novelty is the extension to cardinality ≥ 3, where the
core-free branch disappears, not the three-state message itself. WH00, MP01 and
YFW05 are now cited in the text.

**§4.6.** The abstract now says what `⟨k⟩(c-1) = e` is a threshold *of* — the
loss of stability of warning propagation, i.e. of leaf-removal certification —
in the same sentence, and the caveat is stated where the mixed-cardinality and
correlation results are presented rather than only afterwards.

## Minor points

1 (H overloaded), 2 (structural exponents renamed θ, β reserved for inverse
temperature), 3 (Eq. 5 as one substitution with a note), 4 (the held-fixed and
reported quantities defined explicitly), 5 (a paragraph distinguishing
implementation checks from physical validation), 6 (the residual divided by Λ²
now shown converging to −0.61 over a range, rather than a bound at one point),
8 (new figure: both branches, the coexistence window, and the free-energy
crossing), 9 (Sec. IX compressed to a paragraph, method table to Supplemental),
11 (continuity and monotonicity conditions stated), 12 (ensemble specified),
13 (the rewired control given its own paragraph). 7 is addressed by three new figures: the
clustering fraction at fixed degree (which makes §4.1 visible at a glance, with
the Monte Carlo on the same axes), the simplicial double transition with the
free-energy difference that locates it, and the hyperbolic-random-graph
comparison. 14: the derivations the report named as compressed past
reproducibility are restored — the collapse of Eq. (18) to one scalar and the
elimination giving `⟨k⟩ = e/(c−1)`, the multi-layer instability condition with
the three-to-one mixture worked through to `σ = 0.3777`, `⟨k⟩ = 3.415`, and the
branch construction behind the `q = 16` free-energy comparison. Each is now
checked by a test as well as displayed.

10: confirmed. The accepted manuscript carries "Article in Press", received
21 December 2025 and accepted 5 June 2026, with the publisher's own recommended
citation giving no volume or article number. The reference is now "Commun. Phys.
(2026), in press (accepted 5 June 2026)" with the DOI; the volume should be
filled at proof.

Consulting the accepted version also improved Sec. IV C. The two treatments fix
different normalisations --- Ref. [SLG26] sets rho_q J_q = 1, this paper had
been quoting fixed J --- so the earlier statement that the Bethe spinodal peaks
at q = 3, while correct at fixed J, was not comparable with their Eq. (8).
Imposing their normalisation, J_q = q/k, the Bethe spinodal converges to
T* = q(q-1)/2^(q-1) with a residual falling as 1/k, and their maximum of 3/2
shared between q = 3 and q = 4 is recovered. The exception is q = 4 itself,
where the leading finite-connectivity correction vanishes and the approach is
1/k^2. That is a stronger check than the manuscript previously contained, and it
resolves what would otherwise have read as a disagreement.

## §4.15

Each quantitative claim has been checked. The corrections above were reached by
re-deriving the disputed result by hand before touching the code, and the
Monte Carlo of §4.8 was written to be independent of the formalism it tests.
