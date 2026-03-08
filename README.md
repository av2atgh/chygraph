# chygraph
Symbolic calculation of percolation in chygraphs, including random networks, hypergraphs and higher order networks.

Source:

Alexei Vazquez, Percolation in higher order networks via mapping to chygraphs

https://arxiv.org/abs/2308.00987

## Installation

Install from the source distribution:

```bash
pip install chygraph-0.1.0.tar.gz
```

Or install in editable mode from the repository:

```bash
pip install -e .
```

## Usage

```python
from chygraph import HypergraphPercolation, MultiplexHypergraph, InteractingHypergraphs, GraphWithTriangles

H = HypergraphPercolation()
print(H.A.theta())
print(H.A.eigenvals())
```

See the `examples/` folder for more.
