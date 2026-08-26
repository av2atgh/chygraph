"""Separating geometry from the heavy tail.

The main probe runs tau = 2.1/2.5/2.9, all below 3, where <k^2> diverges.  A
degree-matched configuration model then grows maximal cliques among its hubs
whether or not there is any geometry, so a growing clique moment in the HRG
there does not tell clustering from the tail.

Above tau = 3 the second moment of the degree is finite and the control's
cliques should converge, exactly as they do for Erdos-Renyi.  The HRG is still
geometric and still clustered.  So this is the run that isolates the question:

    tau > 3, HRG clique moments grow, control converges  -> geometry, and the
        clique route to prediction 4 is worth pursuing at light tails
    tau > 3, both converge                               -> the clique ensemble
        is fine above 3, and the failure below 3 is the tail, not geometry
"""

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path.home() / 'av2atg' / 'computational_complexity' / 'code'))
from hrg import erased_configuration_model, hrg_calibrated  # noqa: E402

sys.path.insert(0, str(Path(__file__).parent))
from er_control import moments  # noqa: E402

TAUS = (3.5, 4.5)
KBARS = (4.0, 8.0)
NS = (10_000, 30_000, 100_000, 300_000)
SEEDS = (1, 2, 3)

if __name__ == '__main__':
    print("tau > 3: <k^2> finite, so the control should converge like ER.")
    print(f"{'tau':>5}{'kbar':>6}{'family':>8}" +
          "".join(f"{n:>10}" for n in NS) + "   (sbar; c_max in brackets)")
    for tau in TAUS:
        for kbar in KBARS:
            rows = {'hrg': [], 'config': []}
            for n in NS:
                acc = {'hrg': [], 'config': []}
                for seed in SEEDS:
                    src, dst, r, th, R, C = hrg_calibrated(
                        n, tau=tau, kbar=kbar, rng=seed, tol=0.01, max_iter=25)
                    deg = np.bincount(np.concatenate((src, dst)), minlength=n)
                    acc['hrg'].append(moments(n, src, dst))
                    acc['config'].append(
                        moments(n, *erased_configuration_model(deg, rng=seed)))
                for fam in rows:
                    cm = np.mean([a[0] for a in acc[fam]])
                    sb = np.mean([a[3] for a in acc[fam]])
                    rows[fam].append((sb, cm))
            for fam in ('hrg', 'config'):
                cells = "".join(f"{sb:>7.3f}[{cm:>2.0f}]" for sb, cm in rows[fam])
                print(f"{tau:>5.1f}{kbar:>6.0f}{fam:>8}{cells}", flush=True)
