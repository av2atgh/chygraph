# chygraph
Symbolic calculation of percolation in chygraphs, including random networks, hypergraphs and higher order networks.

Sources:

Alexei Vazquez, Percolation in higher order networks via mapping to chygraphs
(percolation threshold), https://arxiv.org/abs/2308.00987

Alexei Vazquez, The giant component of complex hypergraphs (order parameter,
critical amplitude, dependent layers), `manuscript_3/`

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

### Giant component fraction

The threshold tensor `A` is the Jacobian, at its trivial fixed point, of a
non-linear self-consistency map.  Solving that map instead of linearising it
gives the order parameter.  The inputs are generating functions rather than
their first derivatives; for Poisson distributions they are fixed by the same
means used above.

```python
from chygraph import hypergraph_giant, graph_with_triangles_giant

G = hypergraph_giant()
print(G.theta())                                        # same as above
print(G.node_fraction({'k': 3, 'c': 3, 'p': 0.6, 'q': 0.8}))
```

### Critical amplitude

`S = B * Lambda + O(Lambda^2)` near the threshold, with `B` in closed form from
second moments alone.

```python
from chygraph import CriticalAmplitude
from sympy import symbols

C = CriticalAmplitude(hypergraph_giant())
print(C.amplitude_at_threshold(0, symbols('q')))        # 4*p/(c*p + 1)
```

### Dependent layers

When a complex's participation in different layers is correlated, the excess
generating functions are derivatives of the joint one and the threshold tensor
picks up inclusion-biased second moments.

```python
from chygraph import JointGiantComponent
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
