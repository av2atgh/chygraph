"""Example: percolation problems from the literature, solved by mapping.

Each construction is a chygraph; the threshold, the order parameter and the
critical amplitude then come from the generic machinery.  The mappings are
documented in chygraph/applications.py and in manuscript_3.
"""

from sympy import Rational, exp, expand, nsolve, simplify, symbols

from chygraph import (and_or_hypergraph, correlated_cardinality_hypergraph,
                      two_class_joint_degree, household_epidemic,
                      clique_network, clique_cluster_distribution,
                      stc_percolation)

p, q, pH, T, k, c = symbols('p q p_H T k c')


def and_versus_or():
    """Bianconi & Dorogovtsev: one chygraph, two intra-complex generating
    functions.  AND-logic (the hyperedge dies if any member is removed) always
    has the higher threshold."""
    print('=== AND vs OR hypergraph site percolation ===')
    for logic, guess in (('or', 0.2), ('and', 0.6)):
        M = and_or_hypergraph(logic)
        pc = float(nsolve(M.theta().subs({k: 2, c: 3}), p, guess))
        print(f'  {logic.upper():3s}  theta = {simplify(M.theta())}')
        print(f'       p_c(<k>=2,<c>=3) = {pc:.5f}   '
              f'S(p=1) = {M.node_fraction({"k": 2, "c": 3, "p": 1}):.5f}')


def hyperdegree_cardinality_correlation():
    """Three hypergraphs with identical marginal hyperdegree distributions and
    different hyperdegree-cardinality correlation."""
    print('\n=== hyperdegree-cardinality correlation (identical marginals) ===')
    card, ma, mb, sp = [2, 5], Rational(1), Rational(3, 5), Rational(4, 5)
    print(f'{"w":>6}{"<kbar^(1)>_02":>16}{"p_c":>10}{"S(0.5)":>10}{"S(1)":>10}')
    for w in (0, Rational(1, 2), 1):
        M = correlated_cardinality_hypergraph(card, two_class_joint_degree(ma, mb, sp, w))
        pc = float(nsolve(M.theta(), p, 0.3))
        print(f'{float(w):>6}{str(M.kappa_bar(1, 0, 2)):>16}{pc:>10.5f}'
              f'{M.node_fraction({"p": 0.5}):>10.5f}'
              f'{M.node_fraction({"p": 1.0}):>10.5f}')
    print('  <kappa>_02 is the same for all three, so the marginals alone '
          'give the middle row.')


def two_levels_of_mixing():
    """Ball, Mollison & Scalia-Tomba: households plus a global network.  The
    threshold is their household reproduction number R*."""
    print('\n=== SIR with households and a global network ===')
    H = household_epidemic({3: 1})
    print('  theta =', expand(H.theta()))
    print('  (= R* - 1 with R* = T[<kbar> + mu_H <k>])')
    print(f'{"p_H":>6}{"T":>6}{"k":>5}{"final size S":>15}')
    for a, b, d in ((0.0, 0.6, 2.0), (0.5, 0.3, 2.0), (0.8, 0.2, 2.0)):
        print(f'{a:>6}{b:>6}{d:>5}'
              f'{H.node_fraction({"p_H": a, "T": b, "k": d}):>15.5f}')
    print('  within-household final size for n = 3 (= Fig. 3 of cnae047):')
    for j, pr in sorted(clique_cluster_distribution(3, pH).items()):
        print(f'    {j} infected: {expand(pr)}')


def cliques():
    """A network of cliques: bond percolation, SIR with clique-dependent
    transmission, and bipartite projections at p_bond = 1."""
    print('\n=== network of cliques ===')
    for n in (2, 3, 4):
        M = clique_network({n: 1}, p_bond=q)
        print(f'  clique size {n}: theta = {expand(M.theta())}')


def triadic_closure():
    """Cirigliano: site percolation on a static triadic closure graph.

    The chygraph lives on the treelike backbone, not on the clustered graph.
    Mapping the clustered graph itself gives an incidence structure with
    four-cycles and the prediction fails.
    """
    print('\n=== static triadic closure (Poisson backbone) ===')
    M = stc_percolation()
    print('  theta =', expand(simplify(M.theta())))
    print('  phi_c = [(1+b) - sqrt((1+b)^2-4)] / 2b,  b = <k(k-1)>/<k>')
    b = symbols('k')
    pc = float(nsolve(M.theta().subs(b, 3), symbols('phi'), 0.1))
    print(f'  <k>=3:  phi_c = {pc:.6f}   (naive mapping of G1 gives 0.0711,'
          f' hMF gives 0.0635)')
    for v in (0.15, 0.25, 0.40):
        print(f'  phi={v}: S = {M.node_fraction({"phi": v, "k": 3.0}):.6f}')


if __name__ == '__main__':
    and_versus_or()
    hyperdegree_cardinality_correlation()
    two_levels_of_mixing()
    cliques()
    triadic_closure()
