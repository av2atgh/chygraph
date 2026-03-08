"""chygraph — symbolic percolation in Chygraphs, hypergraphs, and higher-order networks.

References:
    Alexei Vazquez, "Percolation in higher order networks via mapping to chygraphs"
    https://doi.org/10.1093/comnet/cnae047
    https://arxiv.org/abs/2308.00987
"""

from chygraph.percolation import (
    PercolationMatrix,
    vec2A,
    HypergraphPercolation,
    MultiplexHypergraph,
    InteractingHypergraphs,
    GraphWithTriangles,
)

__all__ = [
    "PercolationMatrix",
    "vec2A",
    "HypergraphPercolation",
    "MultiplexHypergraph",
    "InteractingHypergraphs",
    "GraphWithTriangles",
]
