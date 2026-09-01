"""Clustering *lowers* the Ising critical temperature at fixed degree.

This is Table I of the manuscript, and the sign matters: an earlier version of
this file claimed the opposite, and the claim was wrong for the reason the
referee gave.  A Poisson layer of ``n`` links against a Poisson layer of ``n/2``
triangles matches the *mean* degree and nothing else -- the triangle
construction gives ``d = 2X`` with ``X ~ Poisson(n/2)``, so the degree is even,
its variance is ``2n`` rather than ``n``, and the excess degree the branching
matrix actually uses is ``n + 1`` rather than ``n``.  Comparing against a link
layer at ``<kbar> = n`` therefore compares two different degree distributions
and reverses the apparent sign.

Two nulls are run here.  The unambiguous one is regular: an ``n``-regular graph
against a network of ``n/2`` triangles per vertex.  Those have identical
``p_d = delta_{d,n}``, identically neutral ``e_dd'`` and the same number of
edges per vertex, and differ only in whether those edges close into triangles.
The second is Poisson against a configuration model carrying the triangle
ensemble's own excess degree.  Both give the same sign.

The mechanism: a triangle transmits ``t/(1-t+t^2) > t`` per traversal, so per
neighbour clustering helps.  But arriving at a degree-``n`` vertex through one
of its triangles leaves only ``n-2`` branches rather than ``n-1``, because one of
the two neighbours in the triangle just traversed is already accounted for and
its influence sits inside ``u'``.  At ``n = 4`` the enhancement would have to
carry ``u'_T`` from ``1/3`` to ``1/2`` to break even and it reaches only ``3/7``.
The lost branch wins.

    python examples/clustering_lowers_tc.py
"""

import numpy as np

from chygraph_statmech import Chygraph, ising


def T_c(cardinalities, degrees, excess):
    return 1.0 / Chygraph(cardinalities, degrees, excess=excess).critical_coupling()


if __name__ == '__main__':
    print(__doc__.split('\n\n')[0])
    print("\nRegular nulls: an n-regular graph against n/2 triangles per vertex."
          "\nIdentical p_d, identical e_dd', same edges per vertex.\n")
    print(f"{'n':>4}{'t_c links':>12}{'t_c tri':>10}{'T_c links':>12}"
          f"{'T_c tri':>10}{'Delta T_c':>11}")
    for n in (4, 6, 8, 10, 20):
        # links: <kbar> = n - 1, condition <kbar> t = 1
        tL = 1.0 / (n - 1)
        # triangles: <kbar>_T = n/2 - 1, condition 2 <kbar>_T u'_T = 1
        TL = T_c([2], [float(n)], [n - 1.0])
        TT = T_c([3], [n / 2.0], [n / 2.0 - 1.0])
        tT = np.tanh(1.0 / TT)
        print(f"{n:>4}{tL:>12.5f}{tT:>10.5f}{TL:>12.4f}{TT:>10.4f}"
              f"{(TT / TL - 1) * 100:>10.1f}%")

    print("\nPoisson layers, against a null carrying the triangle ensemble's own"
          "\nexcess degree <kbar> = n + 1.  Same sign.\n")
    print(f"{'n':>4}{'matched null':>14}{'t_c tri':>10}{'T_c null':>11}"
          f"{'T_c tri':>10}{'Delta T_c':>11}")
    for n in (4, 6, 8, 10, 20):
        Tnull = T_c([2], [n + 1.0], [n + 1.0])          # <kbar> = n + 1
        TT = T_c([3], [n / 2.0], [n / 2.0])             # Poisson is its own excess
        print(f"{n:>4}{np.tanh(1.0 / Tnull):>14.5f}{np.tanh(1.0 / TT):>10.5f}"
              f"{Tnull:>11.4f}{TT:>10.4f}{(TT / Tnull - 1) * 100:>10.1f}%")

    print("\nAnd the comparison that misleads, kept because it does: matching"
          "\nonly the mean degree, a link layer at <kbar> = n.\n")
    for n in (4, 6):
        Tbad = T_c([2], [float(n)], [float(n)])
        TT = T_c([3], [n / 2.0], [n / 2.0])
        print(f"  n = {n}: link layer at <kbar> = {n} gives T_c = {Tbad:.4f},"
              f" triangles {TT:.4f}  ->  {(TT / Tbad - 1) * 100:+.1f}%,"
              f" the wrong sign")

    print("\nThe transmission, per traversal, is the half that does favour order:")
    for bJ in (0.1, 0.3, 0.6):
        t = np.tanh(bJ)
        print(f"  t = {t:.4f}:  edge {ising.clique_derivative(2, bJ):.5f}"
              f"   triangle {ising.clique_derivative(3, bJ):.5f}"
              f"   = t/(1-t+t^2) = {t/(1-t+t*t):.5f}")
    print("\n  Both effects are real; Eq. (8) carries the multiplicity as well"
          "\n  as the transmission, and the multiplicity wins.")
