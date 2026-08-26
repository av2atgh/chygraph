"""Clustering raises the Ising critical temperature at fixed mean degree.

A node in <k_T> triangles has 2<k_T> neighbours, the same count as 2<k_T>
links.  The two are indistinguishable to the {p_d, e_dd'} ensemble of
Vazquez & Weigt (2003).  They are not indistinguishable here, because a
triangle is solved exactly inside itself and transmits u_T = t/(1-t+t^2)
per neighbour rather than t.
"""

import sympy as sp

import chygraph_statmech as cs

t = sp.Symbol('t', positive=True)
K_L, K_T, k_L, k_T = sp.symbols('K_L K_T k_L k_T')

theta = sp.factor(sp.simplify(
    cs.in_tanh(sp.expand(cs.graph_with_triangles_ising().theta()))))
print('theta =', theta, '\n')


def t_c(kL, kT):
    """Poisson (K = k), so the layer coupling term drops out."""
    roots = sp.solve(theta.subs({K_L: kL, k_L: kL, K_T: kT, k_T: kT}), t)
    return min(float(r) for r in roots if r.is_real and 0 < float(r) < 1)


print(f'{"excess nbrs":>12} {"links":>8} {"triangles":>10} '
      f'{"t_c links":>10} {"t_c tri":>9} {"T_c gain":>9}')
for n in (4, 6, 8, 10, 20):
    a, b = t_c(n, 0), t_c(0, n // 2)
    # T_c / J = 1 / atanh(t_c)
    Ta, Tb = 1 / sp.atanh(a), 1 / sp.atanh(b)
    print(f'{n:>12} {n:>8} {n//2:>10} {a:>10.5f} {b:>9.5f} '
          f'{float(Tb / Ta - 1) * 100:>8.2f}%')
