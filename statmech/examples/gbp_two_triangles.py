"""Generalised belief propagation closes what the Bethe counting leaves open.

Section VII of the manuscript closes the complex family under intersection and
assigns Mobius counting numbers, then evaluates that counting on *isolated*
regions -- an estimate, not a fixed point.  This runs the messages instead.

Two triangles sharing an edge is the manuscript's own example (Table IV).  Once
the zero-counting singletons are pruned the region family is
{0,1,2}, {1,2,3}, {1,2}, which is a junction tree, so GBP is exact there.

    python examples/gbp_two_triangles.py
"""

import numpy as np

from statmech import Chygraph
from statmech.gbp import (GBP, clique_edges, exact_log_Z,
                                   ising_factors, static_log_Z)
from statmech.region import RegionGraph

CX = [[0, 1, 2], [1, 2, 3]]

if __name__ == '__main__':
    rg = RegionGraph(CX)
    print("regions and Mobius counting numbers")
    for r, c in sorted(rg.counting.items(), key=lambda kv: (-len(kv[0]), sorted(kv[0]))):
        print(f"    {str(sorted(r)):>12}  c_R = {c:+d}"
              + ("   (pruned: contributes nothing)" if c == 0 else ""))
    print(f"\n  Bethe counting puts -1 on each shared *vertex*; Mobius puts -1"
          f"\n  on the shared *edge*.  Both count nodes once; only one counts"
          f"\n  the bond once, and GBP refuses to run on the other.\n")

    edges = clique_edges(CX)
    print(f"{'beta J':>8}{'exact':>10}{'Bethe err':>12}{'Kikuchi err':>13}"
          f"{'GBP err':>12}{'sweeps':>8}")
    for bJ in (0.2, 0.5, 1.0, 2.0):
        f = ising_factors(edges, bJ)
        ex = exact_log_Z(f, range(4))
        bet = static_log_Z(rg.bethe_counting(), f) - ex
        kik = static_log_Z(rg.counting, f) - ex
        g = Chygraph.region_gbp(CX, bJ, damping=0.0).run()
        print(f"{bJ:>8.1f}{ex:>10.4f}{bet:>+12.2e}{kik:>+13.2e}"
              f"{g.log_Z()-ex:>+12.2e}{g.sweeps:>8}")

    print("\n  The Bethe counting overestimates: it subtracts the two shared"
          "\n  vertices independently and so never removes the correlation"
          "\n  between them.  The Mobius counting removes it and slightly"
          "\n  overshoots -- but only because it is evaluated on isolated"
          "\n  regions.  Iterating the messages closes the gap entirely.")

    g = Chygraph.region_gbp(CX, 1.0, damping=0.0).run()
    print(f"\n  marginal consistency at the fixed point: {g.consistency():.1e}")
    print(f"  magnetisation (zero field, so zero): "
          f"{max(abs(v) for v in g.magnetisation().values()):.1e}")

    print("\n  Where the region graph is NOT a junction tree, GBP is"
          "\n  approximate and needs damping.  K4 covered by its four"
          "\n  triangles -- a cover no maximal-clique family ever produces,"
          "\n  since maximal cliques do not nest:")
    K4 = [[0, 1, 2], [0, 1, 3], [0, 2, 3], [1, 2, 3]]
    rg4 = RegionGraph(K4)
    for bJ in (0.2, 0.5):
        f = ising_factors(clique_edges(K4), bJ)
        ex = exact_log_Z(f, range(4))
        g = Chygraph.region_gbp(K4, bJ, damping=0.99).run(20_000)
        print(f"    beta J = {bJ}:  GBP {g.log_Z()-ex:+.2e} (residual "
              f"{g.residual:.0e}), static Kikuchi "
              f"{static_log_Z(rg4.counting, f)-ex:+.2e}")
