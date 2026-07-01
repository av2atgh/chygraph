# chygraph
Symbolic calculation of percolation in chygraphs, including random networks, hypergraphs and higher order networks.

Source:

Alexei Vazquez, Percolation in higher order networks via mapping to chygraphs

https://arxiv.org/abs/2308.00987

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

```python
from chygraph import HypergraphPercolation, MultiplexHypergraph, InteractingHypergraphs, GraphWithTriangles

H = HypergraphPercolation()
print(H.A.theta())
print(H.A.eigenvals())
```

See the `examples/` folder for more.
