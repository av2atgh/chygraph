"""Sec. 8.10: is a solution cluster itself a chygraph?  Writes no figures.

Summarises ../statmech/probe/results/cluster_blocks.csv, the cached output of
../statmech/probe/cluster_blocks.py (exhaustive 3-SAT at N = 16..24, clusters
at Hamming distance one, block decomposition of each cluster against a
size-matched random control), and prints the numbers the section quotes:
over the clusters with two or more blocks, the fraction that are exactly a
product of their block projections, the same for the control, blocks per free
variable, the frozen fraction, and the additivity of the entropy over blocks.
"""
import csv
from collections import defaultdict
from pathlib import Path

import numpy as np

PROBE = Path(__file__).resolve().parents[2] / 'statmech' / 'probe' / 'results'


def summarise(rows, key):
    groups = defaultdict(list)
    for r in rows:
        groups[key(r)].append(r)
    out = []
    for k in sorted(groups):
        g = groups[k]
        multi = [r for r in g if int(r['nblocks']) >= 2]
        ctrl_multi = [r for r in g if int(r['ctrl_nblocks']) >= 2]
        prod = np.mean([int(r['product']) for r in multi]) if multi else np.nan
        # control: a random set of the same size; the product test is meaningful
        # only when its dependency graph has >= 2 components, else it is trivially
        # a product with one block.  Quote both the fraction of controls that split
        # at all and the fraction that are products among those.
        ctrl_split = np.mean([int(r['ctrl_nblocks']) >= 2 for r in g])
        ctrl_prod = np.mean([int(r['ctrl_product']) for r in ctrl_multi]) if ctrl_multi else np.nan
        bpf = np.mean([int(r['nblocks']) / int(r['free']) for r in g if int(r['free']) > 0])
        frozen = np.mean([int(r['frozen']) / int(r['n']) for r in g])
        add = np.mean([abs(float(r['entropy']) - float(r['block_entropy'])) < 1e-9 for r in multi]) if multi else np.nan
        # beyond pairwise: the finest product partition (unions of candidate blocks)
        nonprod = [r for r in multi if not int(r['product'])]
        fine_multi = np.mean([int(r['finest_nblocks']) >= 2 for r in multi]) if multi else np.nan
        fine_rescued = np.mean([int(r['finest_nblocks']) >= 2 for r in nonprod]) if nonprod else np.nan
        ctrl_fine = np.mean([int(r['ctrl_finest_nblocks']) >= 2 for r in g if int(r['free']) > 0])
        out.append(dict(key=k, clusters=len(g), multi=len(multi), product=prod,
                        finest_multi=fine_multi, rescued=fine_rescued, ctrl_finest_multi=ctrl_fine,
                        ctrl_split=ctrl_split, ctrl_product=ctrl_prod,
                        blocks_per_free=bpf, frozen=frozen, additive=add,
                        mean_size=np.mean([int(r['size']) for r in g]),
                        mean_blocks=np.mean([int(r['nblocks']) for r in multi]) if multi else np.nan))
    return out


def check_decomposition():
    """A hand-built product is recovered with exactly its blocks."""
    import sys
    sys.path.insert(0, str(PROBE.parent))
    import cluster_blocks as cb
    A = np.array([[0, 0], [1, 1]], bool)
    B = np.array([[0, 0, 0], [0, 1, 1], [1, 0, 1]], bool)
    S = np.array([np.concatenate([a, b]) for a in A for b in B])
    frozen, blocks, prod, sizes = cb.block_decomposition(S)
    assert prod and blocks == [[0, 1], [2, 3, 4]] and sizes == [2, 3], (blocks, sizes)
    print(f'hand-built product recovered with blocks {blocks} of sizes {sizes}')


def main():
    check_decomposition()
    rows = list(csv.DictReader(open(PROBE / 'cluster_blocks.csv')))
    print(f'{len(rows)} clusters of >= 2 solutions')
    for name, key in (('by n', lambda r: int(r['n'])), ('by alpha', lambda r: float(r['alpha'])),
                      ('all', lambda r: 'all')):
        print(f'--- {name}')
        print('key clusters multi product finest_multi rescued ctrl_finest_multi ctrl_split ctrl_product blocks/free frozen additive mean_size mean_blocks')
        for o in summarise(rows, key):
            print(' '.join(f'{v:.3g}' if isinstance(v, float) else str(v) for v in o.values()))


if __name__ == '__main__':
    main()
