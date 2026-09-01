"""Reproduce Fig. 1 of Vazquez & Weigt, PRE 67, 027101 (2003).

Minimal vertex cover on scale-free graphs with the correlations of Eq. (18),
e_dd' = q_d [ r delta_dd' + (1-r) q_d' ], which interpolates from uncorrelated
(r=0) to fully assorted (r=1).

The paper's curves stop where RS breaks, because the iteration stops
converging there.  Here the fixed point is found by bisection on the scalar the
map closes on, so it exists on both sides, and stability is asked separately.
Rows past the break are marked: the number is a fixed point, not the answer.
"""

from statmech import vertexcover as vc

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
        row += f"{vc.cover_size(p, pi):>13.4f}"
        row += "*" if vc.is_unstable(r, e, q, pi) else " "
    print(row)
print("\n* RS unstable; the RS value is no longer the minimal cover size.")
print("gamma=2.5 is the lower curve and gamma=3.0 the upper, as in the paper.")

print("\nTruncation dependence.  gamma=3.0 is converged; gamma=2.5 is not,")
print("because <d> itself converges slowly for a 2.5 tail.  Two decimals of")
print("x_c and of r_RSB are real at gamma=2.5, four at gamma=3.0.")
print(f"{'d_max':>8}{'r_RSB 2.5':>11}{'r_RSB 3.0':>11}"
      f"{'x_c 2.5':>10}{'x_c 3.0':>10}   (x_c at r=0)")
for dmax in (200, 400, 800, 1600, 3200):
    xs = []
    for gamma in (2.5, 3.0):
        pd = vc.scale_free(gamma, dmax)
        e, q, _ = vc.excess(pd)
        xs.append(vc.cover_size(pd, vc.solve(0.0, e, q)))
    print(f"{dmax:>8}{vc.rsb_point(2.5, dmax):>11.4f}"
          f"{vc.rsb_point(3.0, dmax):>11.4f}{xs[0]:>10.4f}{xs[1]:>10.4f}")
