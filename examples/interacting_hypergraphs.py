"""Examples: percolation in Interacting Hypergraphs and Graphs.

Demonstrates InteractingHypergraphs with Poisson and general degree distributions.

WARNING: eigenvalue calculations for the general-degree interacting graph
can consume significant time and memory.
"""

from chygraph import InteractingHypergraphs


def example_interacting_hypergraphs_poisson():
    IHP2 = InteractingHypergraphs(g=2, poisson=True)
    print("theta:", IHP2.A.theta())
    print("eigenvalues:", IHP2.A.eigenvals())
    print("Lambda:", IHP2.A.eigenvals()[4])


def example_interacting_graphs_poisson():
    IGP2 = InteractingHypergraphs(g=2, graph=True, poisson=True)
    print("theta:", IGP2.A.theta())
    print("eigenvalues:", IGP2.A.eigenvals())
    print("Lambda:", IGP2.A.eigenvals()[4])


def example_interacting_graphs_general():
    IG2 = InteractingHypergraphs(g=2, graph=True)
    print("theta:", IG2.A.theta())
    # WARNING: the following eigenvalue calculations can be slow and memory-intensive
    print("eigenvalues:", IG2.A.eigenvals())
    print("Lambda:", IG2.A.eigenvals()[4])


if __name__ == "__main__":
    example_interacting_hypergraphs_poisson()
    example_interacting_graphs_poisson()
    example_interacting_graphs_general()
