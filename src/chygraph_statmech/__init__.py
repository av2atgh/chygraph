"""chygraph_statmech — statistical mechanics on chygraphs.

Extends the Bethe-Peierls treatment of Vazquez & Weigt, Phys. Rev. E 67, 027101
(2003) from correlated random graphs to chygraphs.  See the README for the
programme; this package is WP1, the reweighted stability tensor.
"""

from chygraph_statmech.cavity import (cavity_derivative, emitted_field,
                                      in_tanh, ising_clique,
                                      ising_edge_derivative,
                                      ising_triangle_derivative, tanh_of)
from chygraph_statmech.models import (graph_ising, graph_percolation,
                                      graph_with_triangles_ising)
from chygraph_statmech.stability import (StabilityMatrix, reweight,
                                         uniform_weights)

__all__ = [
    "StabilityMatrix", "reweight", "uniform_weights",
    "emitted_field", "cavity_derivative", "ising_clique",
    "ising_edge_derivative", "ising_triangle_derivative", "in_tanh", "tanh_of",
    "graph_percolation", "graph_ising", "graph_with_triangles_ising",
]
