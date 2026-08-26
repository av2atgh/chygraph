# chygraph
Symbolic calculation of percolation in chygraphs, including random networks, hypergraphs and higher order networks.

Sources:

Alexei Vazquez, Percolation in higher order networks via mapping to chygraphs
(percolation threshold), https://arxiv.org/abs/2308.00987

Alexei Vazquez, The giant component of complex hypergraphs: automated
generating function calculations (order parameter, critical amplitude,
dependent layers), `manuscript_3/`

## Installation

Install in editable mode from the repository:

```bash
pip install -e .
```

To also install the dependencies used by the plotting examples and the tests:

```bash
pip install -e ".[examples,test]"
```

## Usage

### Percolation threshold

```python
from chygraph import HypergraphPercolation, MultiplexHypergraph, InteractingHypergraphs, GraphWithTriangles

H = HypergraphPercolation()
print(H.A.theta())
print(H.A.eigenvals())
```

### One object, three orders

The threshold tensor `A` is the Jacobian, at its trivial fixed point, of a
non-linear self-consistency map.  A chygraph is specified once, by its
generating functions, and then answers all three questions:

```python
from chygraph import hypergraph_giant

M = hypergraph_giant()

M.theta()                 # c*k*p*q - 1              threshold, symbolic
M.Lambda()                # sqrt(c*k*p*q) - 1        order parameter of the threshold
M.amplitude()             # 4*p/(c*p + sqrt(c*k*p*q))  S = B*Lambda + O(Lambda^2)
M.node_fraction({'k': 3, 'c': 3, 'p': 0.6, 'q': 0.8})   # 0.5082...  S itself
```

`theta`, `Lambda` and `amplitude` are closed-form symbolic; `node_fraction`
solves the map by monotone iteration from `Q = 0`, since the fixed point has no
closed form in general.  `amplitude_at_threshold` reduces `B` on the critical
manifold, `curvature` returns the `C` that must be finite and positive for `B`
to mean anything, and `amplitude_numeric` re-derives `B` by numerical linear
algebra on the full index set as an independent check.

The input is four tables of generating functions -- the same data the threshold
needs, as functions rather than as their first derivatives at 1.  For Poisson
distributions they are fixed by the means already in use.

### Dependent layers

When a complex's participation in different layers is correlated, specify the
chygraph by joint generating functions instead.  `JointChygraph` derives the
excess functions from them and carries inclusion-biased second moments where the
published tensor carries unconditional first moments; the rest of the surface
above is unchanged.

```python
from chygraph import JointChygraph
```

### Constructions from the literature

`chygraph.applications` holds published percolation and spreading problems,
each with its chygraph mapping written out in the docstring: which layer holds
what, what `Phi^l` is, and what `G^l` is.

```python
from chygraph import (and_or_hypergraph, correlated_cardinality_hypergraph,
                      household_epidemic, clique_network)

and_or_hypergraph('and').theta()      # hypergraph vs factor graph percolation
household_epidemic({3: 1}).theta()    # = R* - 1, the household reproduction number
```

| construction | layers | reference |
|---|---|---|
| OR / factor graph percolation | 0 nodes, 1 hyperedges | Bianconi & Dorogovtsev, PRE 109, 014306 (2024) |
| AND / hypergraph percolation | 0 nodes, 1 hyperedges | same; Ha, Neri & Annibale (2025) |
| hyperdegree-cardinality correlation | 0 nodes, 1..L hyperedges by cardinality | Ha et al.; Valdez & La Rocca (2026) |
| two levels of mixing | 0 people, 1 households, 2 global contacts | Ball, Mollison & Scalia-Tomba (1997) |
| network of cliques / bipartite projection | 0 nodes, 1 cliques | Fujiki & Mizutaka (2024) |
| static triadic closure / extended range R=2 | 0 occupied, 1 vacant, 2 backbone edges | Cirigliano, EPL 152, 31002 (2025) |

See the `examples/` folder for more, `manuscript_3/` for the derivations, and
`TODO.md` for constructions that need an extension to the formalism rather than
just a mapping (directed chygraphs, discontinuous transitions, triadic
percolation).
