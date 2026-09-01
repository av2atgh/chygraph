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

Each package installs on its own and carries its own test suite and examples:

```sh
pip install -e percolation
pip install -e statmech
pytest percolation/tests statmech/tests
```

The book builds with `latexmk -pdf -interaction=nonstopmode main.tex` from
`book/`; `book/README.md` is the working log and carries the build checks, the
figure-script map and a dated record of every drafting and revision pass.

Both packages were separate repositories until 2026-09-01 —
`av2atgh/chygraph` and `av2atgh/chygraph_statmech`. The latter is archived; its
history is merged into this one, so every commit either README cites still
resolves here.
