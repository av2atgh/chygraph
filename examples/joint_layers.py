"""Example: chygraphs whose layers are statistically dependent.

Three joint distributions of the number of links and triangles a node
participates in, all with the same marginals, give three different percolation
thresholds and three different order parameters.  The published threshold tensor
sees only the marginals and predicts the independent case for all three.
"""

from sympy import Rational, symbols, nsolve

from chygraph import JointChygraph

q = symbols('q')
HALF = Rational(1, 2)


def triangle_parts():
    """Intra-complex generating functions of the link and triangle layers.

    The triangle excess must be supplied: a triangle's excess component size
    under bond percolation is not the size-biased cardinality, so deriving it by
    differentiation would be the wrong relation.
    """
    tri_bar = lambda y: ((3 * q**2 - 2 * q**3) * y[0]**2
                         + 2 * q * (1 - q)**2 * y[0]
                         + (1 - q)**3 + q * (1 - q)**2)
    Gbar = [[None] * 3 for _ in range(3)]
    Gbar[1][0] = lambda y: (1 - q) + q * y[0]
    Gbar[2][0] = tri_bar
    G = [None,
         lambda y: (1 - q) * y[0] + q * y[0]**2,
         lambda y: y[0] * tri_bar(y)]
    return G, Gbar


# marginals shared by all three: k_| in {1, 3}, k_T in {0, 2}, each w.p. 1/2
PHI = {
    'correlated':      lambda x: HALF * x[1] + HALF * x[1]**3 * x[2]**2,
    'anti-correlated': lambda x: HALF * x[1] * x[2]**2 + HALF * x[1]**3,
    'independent':     lambda x: (HALF * x[1] + HALF * x[1]**3) * (HALF + HALF * x[2]**2),
}


def main():
    G, Gbar = triangle_parts()
    print(f"{'model':<17}{'<kbar^(1)>_0T':>15}{'q_c':>10}{'S(q=1)':>10}{'B':>10}")
    for name, Phi in PHI.items():
        M = JointChygraph(Phi=[Phi, None, None], G=G, Gbar=Gbar)
        qc = float(nsolve(M.theta(), q, 0.3))
        print(f'{name:<17}{str(M.kappa_bar(1, 0, 2)):>15}{qc:>10.5f}'
              f'{M.node_fraction({"q": 1.0}):>10.5f}'
              f'{float(M.amplitude().subs(q, qc)):>10.5f}')
    print('\n<kappa>_0T = 1 for all three: the published tensor uses that value '
          'in\nthe off-diagonal slot and cannot distinguish them.')


if __name__ == '__main__':
    main()
