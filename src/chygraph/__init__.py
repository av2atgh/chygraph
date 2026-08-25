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
from chygraph.giant import (
    GiantComponent,
    poisson_pgf,
    thin,
    finite_pgf,
    constant_pgf,
    moment_pgf,
    hypergraph_giant,
    multiplex_hypergraph_giant,
    graph_with_triangles_giant,
)
from chygraph.amplitude import CriticalAmplitude
from chygraph.joint import JointGiantComponent
from chygraph.applications import (
    clique_cluster_distribution,
    clique_excess_pgf,
    clique_pgf,
    size_biased,
    and_or_hypergraph,
    correlated_cardinality_hypergraph,
    two_class_joint_degree,
    household_epidemic,
    clique_network,
    stc_percolation,
)

__all__ = [
    "PercolationMatrix",
    "vec2A",
    "HypergraphPercolation",
    "MultiplexHypergraph",
    "InteractingHypergraphs",
    "GraphWithTriangles",
    "GiantComponent",
    "poisson_pgf",
    "thin",
    "finite_pgf",
    "constant_pgf",
    "moment_pgf",
    "hypergraph_giant",
    "multiplex_hypergraph_giant",
    "graph_with_triangles_giant",
    "CriticalAmplitude",
    "JointGiantComponent",
    "clique_cluster_distribution",
    "clique_excess_pgf",
    "clique_pgf",
    "size_biased",
    "and_or_hypergraph",
    "correlated_cardinality_hypergraph",
    "two_class_joint_degree",
    "household_epidemic",
    "clique_network",
    "stc_percolation",
]
