"""statmech — statistical mechanics on chygraphs.

Extends the Bethe-Peierls treatment of Vazquez & Weigt, Phys. Rev. E 67, 027101
(2003) from correlated random graphs to chygraphs.  See the README for the
programme; this package is WP1, the reweighted stability tensor.
"""

from statmech.api import Chygraph
from statmech.cavity import (cavity_derivative, emitted_field,
                                      in_tanh, ising_clique,
                                      ising_edge_derivative,
                                      ising_triangle_derivative, tanh_of)
from statmech.fixedpoint import FixedPointStability
from statmech.population import CavityPopulation, critical_coupling
from statmech.models import (graph_ising, graph_percolation,
                                      graph_with_triangles_ising)
from statmech import (antimonotone, core, cover, freeenergy,
                               gbp, hittingset, ising, region, simplicial,
                               softfield, vertexcover)
from statmech.gbp import GBP
from statmech.stability import (StabilityMatrix, reweight,
                                         uniform_weights)

__all__ = [
    "Chygraph",
    "StabilityMatrix", "reweight", "uniform_weights",
    "emitted_field", "cavity_derivative", "ising_clique",
    "ising_edge_derivative", "ising_triangle_derivative", "in_tanh", "tanh_of",
    "graph_percolation", "graph_ising", "graph_with_triangles_ising",
    "FixedPointStability", "vertexcover", "hittingset", "antimonotone",
    "CavityPopulation", "critical_coupling", "core", "cover", "ising",
    "region", "freeenergy", "softfield", "simplicial", "gbp", "GBP",
]
