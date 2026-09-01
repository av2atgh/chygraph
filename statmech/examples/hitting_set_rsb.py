"""Where minimum hitting set stops being easy.

Three questions, in order of how much they move the answer:
does hyperedge cardinality matter, does spread in cardinality matter, and does
correlation between a node's participation across cardinalities matter.

Everything is measured in *neighbours* per node, k <cbar>, so the comparisons
are at matched mean degree in the underlying graph sense.
"""

import numpy as np

from chygraph_statmech import hittingset as hs

E = np.e
CS, M = [2, 6], [0.75, 0.25]
CBAR = sum(m * (c - 1) for m, c in zip(M, CS)) / sum(M)

print("1. Fixed cardinality: k_RSB = e/(c-1), so k(c-1) = e for every c.")
print(f"{'c':>5}{'k_RSB':>12}{'e/(c-1)':>12}{'k(c-1)':>16}")
for c in (2, 3, 4, 5, 10, 20):
    k = hs.rsb_point([c], [1.0])
    print(f"{c:>5}{k:>12.6f}{E/(c-1):>12.6f}{k*(c-1):>16.12f}")
print("  c=2 is Bauer-Golinelli core percolation at c = e.")

print("\n2. Spread in cardinality, at fixed <cbar> = 2.  Heterogeneity")
print("   postpones RSB, the same direction degree heterogeneity does in VW03.")
print(f"{'mix':>24}{'k_RSB':>10}{'k<cbar>':>10}")
for label, w, cs in (("all c=3", [1.0], [3]),
                     ("half c=2, half c=4", [.5, .5], [2, 4]),
                     ("3/4 c=2, 1/4 c=6", [.75, .25], [2, 6]),
                     ("9/10 c=2, 1/10 c=12", [.9, .1], [2, 12])):
    cb = sum(wi * (c - 1) for wi, c in zip(w, cs))
    print(f"{label:>24}{hs.rsb_point(cs, w):>10.5f}"
          f"{hs.rsb_point(cs, w)*cb:>10.5f}")


def matched(sign, s):
    """Two equal-weight node classes, identical marginals, correlated across
    layers.  s = 0 is the uncorrelated Poisson ensemble."""
    def f(t):
        if sign > 0:
            a, b = [t*m*(1+s) for m in M], [t*m*(1-s) for m in M]
        else:
            a = [t*M[0]*(1+s), t*M[1]*(1-s)]
            b = [t*M[0]*(1-s), t*M[1]*(1+s)]
        return hs.HittingSet(CS, hs.two_class_phi(a, b, 0.5))
    return f


print(f"\n3. Correlation across cardinality layers, cardinalities {CS},")
print(f"   marginal split {M}, <cbar> = {CBAR}.  Negative correlation brings")
print("   RSB forward -- an axis with no counterpart in VW03.")
print(f"{'spread s':>9}{'positive':>11}{'isolated':>10}"
      f"{'negative':>11}{'isolated':>10}")
for s in (0.0, 0.2, 0.4, 0.6, 0.8, 0.95):
    out = f"{s:>9.2f}"
    for sign in (+1, -1):
        f = matched(sign, s)
        t = hs.rsb_scale(f)
        out += f"{t*CBAR:>11.5f}{f(t).isolated_fraction():>10.3f}"
    print(out)
print("\n  'isolated' is Phi(0..0), the fraction of vertices in no hyperedge.")
print("  Where it grows, the shift is dilution rather than correlation: the")
print("  positive family's low class empties as s -> 1, so only its small-s")
print("  entries are a clean correlation measurement.  The negative family")
print("  keeps both classes populated and moves monotonically throughout.")
