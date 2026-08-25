"""Example: from the percolation threshold to the giant component fraction.

``theta`` and ``Lambda`` are properties of the linearisation of a non-linear
self-consistency map; this example solves the map itself.
"""

import math

from sympy import symbols, nsolve, expand, simplify

from chygraph import hypergraph_giant, graph_with_triangles_giant


def example_hypergraph():
    G = hypergraph_giant()
    print("theta:", simplify(G.theta()))
    print("k=3, c=3, p=0.6, q=0.8 ->  S =",
          G.node_fraction({'k': 3, 'c': 3, 'p': 0.6, 'q': 0.8}))


def example_er_graph():
    G = hypergraph_giant(graph=True)
    print("\nER graph, S vs the textbook S = 1 - exp(-kS):")
    for k in (1.5, 2.0, 3.0):
        S = G.node_fraction({'k': k, 'p': 1, 'q': 1})
        print(f"  k={k}: S={S:.8f}   1-exp(-kS)={1 - math.exp(-k * S):.8f}")
    print("  critical amplitude B =",
          G.amplitude_numeric({'k': 1, 'p': 1, 'q': 1}), "(exact: 4)")


def example_triangles():
    q, kL, kT = symbols('q k_L k_T')
    T = graph_with_triangles_giant()
    print("\nGraph with links and triangles:")
    print("  theta:", expand(simplify(T.theta())))
    base = {kL: 1, kT: 0.5}
    qc = float(nsolve(T.theta().subs(base), q, 0.4))
    print(f"  bond percolation threshold q_c = {qc:.6f}")
    print(f"  critical amplitude B = {T.amplitude_numeric({'k_L': 1, 'k_T': 0.5, 'q': qc}):.6f}")
    for qq in (0.5, 0.7, 1.0):
        print(f"  q={qq}: S = {T.node_fraction({'k_L': 1, 'k_T': 0.5, 'q': qq}):.6f}")


if __name__ == "__main__":
    example_hypergraph()
    example_er_graph()
    example_triangles()
