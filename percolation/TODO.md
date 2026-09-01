# TODO — extensions deferred to separate work

Constructions that need a change to the chygraph formalism itself, not just a
new mapping. Everything addressable by a mapping alone lives in
`src/chygraph/applications.py`; this file records what does not.

---

## 1. Directed chygraphs

**Why the current formalism cannot do it.** Inclusion is a symmetric relation in
`chygraph.giant`: a complex is entered from below (`i = -`) or from above
(`i = +`), and both directions carry the same information. Directed higher-order
interactions have an intrinsic asymmetry — damage to a reactant suppresses a
reaction, but not conversely — so each inclusion needs an orientation and the
component structure splits into in-, out- and strongly connected components.

**What the extension looks like.** Split each layer index into in/out roles,
doubling the index set from `2L²` to `4L²`, with two coupled families of
generating functions in place of `Phi^l` and `G^l`. The undirected case must
come back when the two coincide. Directed percolation on ordinary graphs is the
`L = 1` sanity check.

**Caveat that must be handled.** Ha, Neri & Annibale show that for hypergraphs
under AND-logic the strongly connected component is *not* the intersection of
the in- and out-components, which it is for ordinary directed graphs. The naive
product ansatz is therefore wrong and needs care.

**Literature.**
- Sun, Liu & Bianconi, *Directionality and node heterogeneity reshape criticality
  in hypergraph percolation*, arXiv:2601.20726 (2026) — message passing,
  anomalous critical exponents.
- Ha, Neri & Annibale, *Connected components in networks with higher-order
  interactions*, arXiv:2504.03060 — AND/OR logic, directed case.
- Traversa et al., *Robustness and complexity of directed and weighted metabolic
  hypergraphs* (2023) — the driving application.

---

## 2. Discontinuous and double transitions

**Why the current formalism cannot do it.** The curvature
`C = l . M[r, r]` of `chygraph.amplitude` is a sum of non-negative terms, so
within the chygraph class the transition is continuous with `beta = 1` whenever
`C > 0`, and `Lambda = 0` locates it. Clique percolation with overlap `l > 1`
has a *discontinuous* node fraction alongside a continuous clique-cluster
fraction, and dense simplicial complexes show double transitions. Those cannot
be branching-process extinction transitions of the present map.

**What to do.** Two separate questions, worth keeping apart:
1. Is the discontinuity an artefact of which layer's order parameter is being
   watched — `S_0` versus `S_l`? The formalism gives both, so this is checkable
   now and may not need an extension at all.
2. For genuine hybrid transitions (interdependent networks), the giant component
   appears by saddle-node bifurcation while `Q = 1` is still stable. Then
   `Lambda = 0` does not locate the transition and the physical fixed point is
   reached by a cascade from above, not by monotone iteration from `Q = 0`.
   This needs a different criterion than `theta` or `Lambda`.

**Literature.**
- Zhao, Li, Peng, Zhong & Wang, Appl. Math. Comput. 431, 127330 (2022) — clique
  percolation on random graphs (already cited in cnae047).
- Higher-order percolation in simplicial complexes, Chaos Solitons Fractals
  (2021) — double transitions.
- Buldyrev et al., Nature 464, 1025 (2010); Radicchi, Nat. Phys. 11, 597 (2015)
  — interdependent networks.
- *Higher-order interdependent percolation on hypergraphs*, Chaos Solitons
  Fractals 177 (2023).

---

## 3. Triadic percolation

**Why the current formalism cannot do it.** The chygraph map is static. Triadic
percolation iterates percolation with signed regulation: a third complex
up- or down-regulates whether an inclusion is active, so the generating
functions change between steps and the giant component becomes a dynamical
variable with a logistic-map route to chaos.

**What the extension looks like.** The chygraph map supplies the static kernel
`S(theta)` of that iteration exactly, so the dynamics could be built on top of
it rather than rederived. One concrete lead worth checking cheaply: the route to
chaos is governed by the derivative of the iterated map at its fixed point, and
near the percolation threshold that derivative involves `dS/dLambda = B`, the
critical amplitude computed in `chygraph.amplitude`. If that connection holds,
the amplitude predicts where the period doubling starts.

**Literature.**
- Millán, Sun, Torres & Bianconi, *Triadic percolation induces dynamical
  topological patterns in higher-order networks*, PNAS Nexus 3, pgae270 (2024).
- *Higher-order triadic percolation on random hypergraphs*, Phys. Rev. E 110,
  064315 (2024), arXiv:2407.14213.
- *Triadic percolation on multilayer networks*, arXiv:2510.09341.

---

## Out of class entirely

Not extensions — these are not branching-process extinction problems and should
not be forced into the formalism:

- Threshold / complex contagion with group activation, where a hyperedge fires
  only above a fraction of infected members.
- Bootstrap and k-core percolation; hypergraph k-core (Phys. Rev. E 109, 014307);
  hyper-cores (arXiv:2301.04235).

Keating & Hébert-Dufresne (arXiv:2511.15688) argue that group structure alone
gives only continuous transitions and that discontinuity requires intersecting
transmission chains. The non-negativity of `C` is a short proof of that claim
for anything that maps to a chygraph, and delimits what the formalism can
represent.
