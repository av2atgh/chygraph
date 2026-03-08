"""Examples: percolation in Multiplex Hypergraphs.

Demonstrates MultiplexHypergraph with 2 types, 3 types, and Poisson distributions.
"""

from chygraph import MultiplexHypergraph


def example_multiplex_2types():
    MH = MultiplexHypergraph(number_of_types=2)
    print("theta:", MH.A.theta())
    print("eigenvalues:", MH.A.eigenvals())
    print("Lambda:", MH.A.eigenvals()[3])


def example_multiplex_3types():
    MH3 = MultiplexHypergraph(number_of_types=3)
    print("theta:", MH3.A.theta())


def example_multiplex_poisson():
    MHP5 = MultiplexHypergraph(number_of_types=4, poisson=True)
    ev = MHP5.A.eigenvals()
    print("theta:", MHP5.A.theta())
    print("eigenvalues:", ev)
    print("Lambda:", ev[2])


if __name__ == "__main__":
    example_multiplex_2types()
    example_multiplex_3types()
    example_multiplex_poisson()
