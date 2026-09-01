"""Finite-y 1RSB on a chygraph: the reweighting, and the kappa x y plane.

The extension is a substitution, not a derivation. At reweighting y each cluster
carries e^{-yE} with E the number of violated complexes, and in the cavity
recursion the energy shift of adding one complex is 1 if that complex is
violated. So the interior sum of Eq. (8.4) runs over ALL configurations with the
violating ones weighted e^{-y} instead of excluded -- which is Ch. 8's -beta H
taking a different argument. Nothing upstream or downstream changes.

  y -> infinity   the hard constraint, and every result in Secs. 12.6-12.9
  finite y        metastable states at energy density e(y)

This lives in probe/ and not in book/figs/ because the finite-y numbers have no
external benchmark here. The y -> infinity reduction IS checked (exactly, to
every digit, at c = 2 and c = 3) and that check is in figs/colouring.py.

Reference: Krzakala, Pagnani & Weigt, Phys. Rev. E 70, 046705 (2004), which sets
out the apparatus for graphs -- reweighting y, complexity Sigma(y), energy e(y),
and y = infinity as the ground-state limit. Its stability criterion is what
Sec. 12.3 already borrows.
"""
import sys
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'book' / 'figs'))
import colouring as C            # noqa: E402  (the machinery lives there)

YS = (1.5, 2.0, 3.0, 4.0, 6.0, 10.0, np.inf)


def plane(q, kappas, c=3, ys=YS):
    """Sigma on the kappa x y plane. Degrees are 2*kappa for triangles."""
    print(f'=== q={q}, c={c}: Sigma on the kappa x y plane ===')
    head = '  '.join(f'y={"inf" if y == np.inf else y:>4}' for y in ys)
    print(f'  kappa  deg  {head}')
    out = {}
    for k in kappas:
        row = [C._sv_sigma(q, c, k, yreweight=y)[0] for y in ys]
        out[k] = row
        print(f'  {k:5.2f} {(2 if c == 3 else 1) * k:5.1f}  '
              + '  '.join(f'{s:+7.3f}' for s in row))
    return out


def check_reduction():
    """y = infinity must reproduce the hard-constraint results exactly."""
    for q, c, k in ((3, 2, 4.5), (4, 2, 8.3), (3, 3, 1.6), (4, 3, 3.6)):
        a = C._sv_sigma(q, c, k)[0]
        b = C._sv_sigma(q, c, k, yreweight=np.inf)[0]
        assert a == b, (q, c, k, a, b)
    print('  y = inf reproduces the hard constraint exactly at c = 2 and c = 3')


if __name__ == '__main__':
    check_reduction()
    print()
    plane(4, (3.5, 3.8, 4.2, 4.6, 5.0, 5.5))
    print()
    plane(3, (1.6, 1.8, 2.0, 2.3, 2.6, 3.0))
