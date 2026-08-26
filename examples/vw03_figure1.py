"""Reproduce Fig. 1 of Vazquez & Weigt, PRE 67, 027101 (2003).

Minimal vertex cover on scale-free graphs with the correlations of Eq. (18),
e_dd' = q_d [ r delta_dd' + (1-r) q_d' ], which interpolates from uncorrelated
(r=0) to fully assorted (r=1).

The paper's curves stop where RS breaks, because the iteration stops
converging there.  Here the fixed point is found by bisection on the scalar the
map closes on, so it exists on both sides, and stability is asked separately.
Rows past the break are marked: the number is a fixed point, not the answer.
"""

from chygraph_statmech import vertexcover as vc

DMAX = 800

print("RS breaking point (VW03: 'a certain value of r that depends on gamma')")
for gamma in (2.5, 3.0):
    print(f"  gamma = {gamma}:  r_RSB = {vc.rsb_point(gamma, DMAX):.4f}")

print("\nMinimal vertex cover size x_c")
print(f"{'r':>6}{'gamma=2.5':>14}{'gamma=3.0':>14}")
for i in range(11):
    r = i / 10
    row = f"{r:>6.1f}"
    for gamma in (2.5, 3.0):
        p = vc.scale_free(gamma, DMAX)
        e, q, _ = vc.excess(p)
        pi = vc.solve(r, e, q)
        row += f"{vc.cover_size(p, pi):>13.6f}"
        row += "*" if vc.is_unstable(r, e, q, pi) else " "
    print(row)
print("\n* RS unstable; the RS value is no longer the minimal cover size.")
print("gamma=2.5 is the lower curve and gamma=3.0 the upper, as in the paper.")

print("\nTruncation dependence of r_RSB")
print(f"{'d_max':>8}{'gamma=2.5':>12}{'gamma=3.0':>12}")
for dmax in (200, 400, 800, 1600):
    print(f"{dmax:>8}{vc.rsb_point(2.5, dmax):>12.4f}"
          f"{vc.rsb_point(3.0, dmax):>12.4f}")
