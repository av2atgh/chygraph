"""Examples: percolation in Hypergraphs and Graphs.

Demonstrates HypergraphPercolation for both hypergraph and graph cases.
"""

from chygraph import HypergraphPercolation


def example_hypergraph():
    H = HypergraphPercolation()
    print("theta:", H.A.theta())
    print("eigenvalues:", H.A.eigenvals())
    print("Lambda:", H.A.eigenvals()[1])


def example_graph():
    G = HypergraphPercolation(graph=True)
    print("theta:", G.A.theta())
    print("eigenvalues:", G.A.eigenvals())
    print("Lambda:", G.A.eigenvals()[1])


if __name__ == "__main__":
    example_hypergraph()
    example_graph()
