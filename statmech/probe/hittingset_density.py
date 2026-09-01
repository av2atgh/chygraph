"""Data for Fig. 5: hitting-set density, hard field against soft field.

Table III of the manuscript states the comparison at eight ensembles; the
referee asked for it made visible.  This produces the three series the figure
needs and writes them to ``probe/results/hittingset.json``.

    python probe/hittingset_density.py

Panel A: Poisson layers, ``rho`` against chy-degree at cardinality 2, 3 and 4,
both treatments, with the Weigt-Hartmann closed form as the exact curve at
``c = 2`` where hard and soft must agree.

Panel B: the regular hypergraphs of Mezard & Tarzia [Phys. Rev. E 76, 041124
(2007)], where the soft field must return their ``rho = 1/K`` and the hard field
does not.  The replica-symmetric entropy is carried alongside, since where it is
negative the RS answer is an underestimate and the point should be read as such.

Panel C: the relative correction against the fraction of a node's complexes that
have cardinality three or more, at fixed total chy-degree.  Layers of
cardinality 2 and ``c`` mixed at ``<k> = 1``, so the right-hand end of each
curve is the pure-``c`` row of Table III.

Roughly ten minutes; the soft-field points are population dynamics.
"""

import json
from pathlib import Path

import numpy as np
from scipy.special import lambertw

import statmech.hittingset as hs
from statmech.hittingset import layer_symbols
from statmech.softfield import HittingSetBP, regular_entropy

OUT = Path(__file__).parent / 'results' / 'hittingset.json'

MU, SIZE, SWEEPS = 60.0, 200_000, 600
KS = [0.5, 0.75, 1.0, 1.5, 2.0, 2.5, 3.0]
REGULAR = ((1, 3), (2, 3), (3, 4), (2, 6), (4, 6), (6, 12))
MIX = [0.0, 0.15, 0.3, 0.5, 0.7, 0.85, 1.0]


def soft(cardinalities, degrees, regular=False, size=SIZE, sweeps=SWEEPS):
    return _soft_run(cardinalities, degrees, regular, size, sweeps)[0]


def _soft_run(cardinalities, degrees, regular=False, size=SIZE, sweeps=SWEEPS,
              tail=40):
    """``(density, spread)``, the second being the drift over ``tail`` more sweeps.

    Population dynamics carries sampling noise of order ``1/sqrt(size)`` whatever
    happens, so a converged run still moves a little; a run that has not settled
    moves by orders of magnitude more.  At ``L = 6``, ``K = 12`` the entropy is
    negative and the iteration does not settle at all, and the manuscript says
    so --- the spread is what says it here, and the figure draws that point
    hollow.
    """
    m = HittingSetBP(cardinalities, degrees, regular=regular, mu=MU,
                     size=size, seed=1, damping=0.5).run(sweeps)
    vals = []
    for _ in range(tail):
        m.sweep()
        vals.append(m.density())
    v = np.asarray(vals)
    return float(v.mean()), float(v.max() - v.min())


def main():
    out = {}

    A = {}
    for c in (2, 3, 4):
        A[c] = {'k': KS,
                'hard': [hs.poisson([c], [k]).cover_size() for k in KS],
                'soft': [soft([c], [k]) for k in KS]}
        print(c, ['%.4f' % v for v in A[c]['soft']], flush=True)
    A['WH'] = [float(1 - (2 * lambertw(k).real + lambertw(k).real ** 2) / (2 * k))
               for k in KS]
    out['A'] = A

    x = layer_symbols(1)[0]
    out['B'] = []
    for L, K in REGULAR:
        rho, spread = _soft_run([K], [L], regular=True, size=100_000, sweeps=500)
        out['B'].append({
            'L': L, 'K': K,
            'hard': float(hs.HittingSet([K], x ** L).cover_size()),
            'soft': rho, 'spread': spread,
            'mt': 1.0 / K,
            's': float(regular_entropy(L, K))})
        print(out['B'][-1], flush=True)

    C = {}
    for c in (3, 4, 6):
        rel = []
        for f in MIX:
            cs, ks = [], []
            if 1 - f > 1e-12:
                cs.append(2); ks.append(1.0 - f)
            if f > 1e-12:
                cs.append(c); ks.append(f)
            h = hs.poisson(cs, ks).cover_size()
            s = soft(cs, ks)
            rel.append((h - s) / s)
            print(c, f, '%.4f %.4f %.3f' % (h, s, rel[-1]), flush=True)
        C[c] = {'x': MIX, 'rel': rel}
    out['C'] = C

    OUT.parent.mkdir(parents=True, exist_ok=True)
    json.dump(out, OUT.open('w'), indent=1)
    print('wrote', OUT)


if __name__ == '__main__':
    main()
