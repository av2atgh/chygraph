# chygraph

Percolation and statistical mechanics on **chygraphs** — complex hypergraphs,
whose complexes may contain other complexes. One repository, three parts.

| | what | |
|---|---|---|
| [`percolation/`](percolation/) | the `percolation` package: self-consistency map, threshold tensor, critical amplitude, correlated layers. Symbolic, on `sympy`. | Part II of the book |
| [`statmech/`](statmech/) | the `statmech` package: cavity recursion with a general interior, branching matrix, Ising / hitting-set / core-percolation solvers, Bethe free energy, region graphs and GBP. Plus `probe/`, the expensive measurements and their cached results. | Parts III and IV |
| [`book/`](book/) | *Phase transitions on complex hypergraphs* — 278 pages, seventeen chapters. Every number in it comes out of the two packages. | |

`statmech` builds on `percolation`: `statmech.stability` imports
`PercolationMatrix`. That dependency is why the two live together.

A third package, **`chygraph`** (`src/chygraph/`), is one namespace over both:

```python
from chygraph import PercolationMatrix, CriticalAmplitude, GBP
```

It re-exports all 58 public names of the two packages, adds nothing and takes
nothing away — `chygraph.percolation` and `chygraph.statmech` remain the real
packages, and importing either directly is unaffected. One name is deliberately
withheld: both define a class called `Chygraph`, and they are different objects
(`percolation.giant.Chygraph` is the self-consistency map, `statmech.api.Chygraph`
the calculation handle). Re-exporting both would let one silently shadow the
other, so `chygraph.Chygraph` raises with a message naming the alternatives, and
the two are available as **`ChygraphPercolation`** and **`ChygraphStatmech`**.
Any *undeclared* collision between the two packages raises on import rather than
being resolved by import order.

Each package installs on its own and carries its own test suite and examples:

```sh
pip install -e percolation
pip install -e statmech
pip install -e .                     # the chygraph namespace over both
pytest tests percolation/tests statmech/tests
```

Or without installing anything, which is what the figure scripts do for
themselves:

```sh
PYTHONPATH=$PWD/src:$PWD/percolation/src:$PWD/statmech/src \
  pytest tests percolation/tests statmech/tests
```

**422 tests, all passing** (2026-09-01, ~6 min). The suite needs `sympy`,
`numpy`, `scipy` and `pytest`, and `statmech/tests/test_core.py` additionally
needs `numba`, because it checks the core against an independent leaf-removal
implementation in a separate repository (`~/av2atg/computational_complexity`).
Without `numba` that one test errors on import and the rest still pass.

The book builds with `latexmk -pdf -interaction=nonstopmode main.tex` from
`book/`; `book/README.md` is the working log and carries the build checks, the
figure-script map and a dated record of every drafting and revision pass.

Both packages were separate repositories until 2026-09-01 —
`av2atgh/chygraph` and `av2atgh/chygraph_statmech`. The latter is archived; its
history is merged into this one, so every commit either README cites still
resolves here.
