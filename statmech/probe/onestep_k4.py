"""Sec. 8.10 / 13.10 at k = 4, where alpha_d and alpha_c differ (9.38 and
9.55 against alpha_s = 9.93): the one-step scan of onestep_sat.py with k = 4."""
import sys
from pathlib import Path
sys.path.append(str(Path(__file__).resolve().parent))
from onestep_sat import scan, OUT  # noqa: E402
if __name__ == '__main__':
    scan(alphas=(9.0, 9.2, 9.35, 9.45, 9.55, 9.65, 9.8, 9.9), ms=(0.2, 0.4, 0.6, 0.8, 1.0),
         seeds=(0, 1), M=1000, P=200, sweeps=300, nsamp=400000, procs=6,
         out=OUT.with_name('onestep_sat_k4.csv'), k=4)
