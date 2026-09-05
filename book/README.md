# Phase transitions on complex hypergraphs (book)

Book-length treatment of percolation and statistical mechanics on complex
hypergraphs. Production style follows *Local network growth* — 5×8 in trim,
grayscale, heavy on TikZ illustration. Derivations run in continuous text, each
ending on the routine that executes it.

## Contents

| | Source | Topic |
|---|--------|-------|
| — | `preface.tex` | The tree assumption, the bargain, what the book claims |
| — | `notation.tex` | Structural symbols, the letters that carry more than one meaning, chapter-local letters |
| **I** | | **Foundations** |
| 1 | `introduction.tex` | Percolation; Ising, vertex cover, hitting set; leaf removal and the core fraction $P_C$; mean field, cavity, replicas, BP; why real networks break treelikeness |
| 2 | `chygraphs.tex` | The object: complexes whose vertices are complexes; every higher-order structure as one thing |
| 3 | `data.tex` | Chygraph representation of real systems: papers, protein complexes, reactions, schedules; cliques of a clustered graph; how much overlap there is |
| **II** | | **Percolation** |
| 4 | `percolation.tex` | The self-consistency map and the threshold tensor, structure by structure |
| 5 | `giant.tex` | Order parameter, critical amplitude, the moment hierarchy, dependent layers, six constructions solved |
| 6 | `epidemics.tex` | SIR as percolation, two levels of mixing, the household reproduction number, contagion inside a group |
| **III** | | **Statistical mechanics** |
| 7 | `potts.tex` | Fortuin–Kasteleyn: percolation and Ising are one recursion at two values of $q$ |
| 8 | `statmech.tex` | General theory: convolution up, exact interior sum down, the branching matrix |
| 9 | `ising.tex` | Ising on chygraphs; clustering lowers $T_c$; the AT line; the unanimity interaction |
| 10 | `hittingset.tex` | Hard fields, where they fail off the graph, soft fields, RSB |
| 11 | `cover.tex` | Vertex cover, leaf removal, core percolation; hyperbolic random graphs |
| 12 | `colouring.tex` | Proper against hypergraph colouring; `tau = -1/(q-1)` independent of cardinality; `(q-1)^2`; a graph with triangles, where the graph calculation is wrong by one; survey propagation; and Sec. 12.9, extending Krzakala et al. to chygraphs to get `c_q` for a clustered graph |
| 13 | `satisfiability.tex` | Clauses as complexes; `alpha = 1` exact at `k = 2`, no linear instability above it; clauses whose members are clauses, and what CNF flattening costs |
| **IV** | | **Complexes with non-trivial overlap** |
| 14 | `overlap.tex` | The price of treelikeness: what `BP_chi` loses to overlapping cliques, measured on three classes of network, with the rewired control that attributes it |
| 15 | `gbp.tex` | The repair the literature already has: region graphs, Möbius counting, GBP over 240 runs — exact where the clique family is chordal, and chordality follows provenance, not clustering |
| 16 | `metacomplex.tex` | The repair that stays inside the formalism: merge complexes sharing two or more atoms; exact iff the merged incidence structure is a forest; and Sec. 16.3, the core-percolation transition at incidence branching one |
| 17 | `outlook.tex` | One recursion, many models; the two running threads; what is not done |
| — | `software.tex` | Back matter. Repo links; equation-to-method-to-test table; how to reproduce a figure |

**Chapters 12 and 13 have no manuscript behind them.** Every other chapter is
exposition of work that exists elsewhere; for these two the calculations are
done in the book's own figure scripts, `figs/colouring.py` and
`figs/satisfiability.py`, and checked against the published thresholds — which
the last section of each chapter **computes** rather than quotes, by survey
propagation at `m = 0`. PDFs of the references are under
`~/Downloads/chygraph_references/`.

## Status

`main.pdf` builds with **0 errors, 0 undefined references and 0
multiply-defined labels, across 364 pages.** Not box-clean: **two overfull
hboxes** — `cover.tex:483--493` (1.99pt, "Which replica-symmetry-breaking
point") and `metacomplex.tex:383--389` (3.16pt) — and 52 underfull vboxes,
every one of them `while \output is active`, which is page-breaking around
floats and not a line running into the margin.

60 figures, 31 numbered tables, 163 numbered equations, 74 references and a
109-term index. Both checks under *Two checks the build cannot make* print
nothing.

**These counts are re-measured, not maintained by hand** — the recipes are
under Building, below. Re-measure rather than trust them after any edit.

`software.tex` sits in the **back matter**, after Ch. 17 and before the index,
and is unnumbered — "The software", page 325, arabic. It carries the two
repository URLs and **Table 1: every computed equation, the routine that
evaluates it, and the test or script that checks it**. If a routine is renamed,
that table is where it has to be fixed.

**Table 1 numbers plainly, and that is deliberate.** It is a `longtable` in the
back matter, and `\backmatter` stops numbering chapters but neither resets the
float counters nor clears `\thechapter`, so `\thetable` used to expand to
`15.<n>` and the software map printed as "Table 15.2". Immediately after
`\backmatter`, `main.tex` now resets the `table` and `figure` counters and
redefines `\thetable`/`\thefigure` to a plain `\arabic`, with a comment saying
why. Add a float to the back matter and it will number plainly. Ch. 12's
caption cites `\citet{gabrie2017}, Table~1`, which is somebody else's Table 1
and is attributed in place, so the two do not collide.

## Underlying material

- `~/Dropbox/submissions/hyperabs.2022/hyperabs_v3.tex` — Phys. Rev. E **107**, 024316 (2023)
- `~/Dropbox/submissions/chygraph.2023/chygraph.tex` — J. Complex Netw. (2024), cnae047
- `~/av2atg/chygraph/percolation/manuscript_3/manuscript.tex` — the percolation extension.
  **Will not be submitted**: the book is the first place Ch. 5's critical
  amplitude and moment hierarchy are reported, and Ch. 5's opening says so.
  Do not add a citation to it.
- the manuscript this book supersedes (`main.tex`, `supplement.tex`, at the
  root of the then-separate `chygraph_statmech` repository) —
  **deleted**, recoverable from git at commit `5ebd892`, the last commit before
  the deletion. Everything in it is in the book; the section numbers named in
  `main.tex`'s `\include` comments refer to that retired file.
- `../statmech/src/statmech/`, `../statmech/examples/`, `../statmech/probe/` — every number in the book comes out of these

## Building

```sh
latexmk -pdf -interaction=nonstopmode main.tex     # -> main.pdf
```

Verify a build with

```sh
grep -c '^!' main.log                    # errors, must be 0
grep -i 'undefined' main.log             # citations and refs, must be empty
grep -i 'multiply defined' main.log      # must be empty; NOT caught by the line above
grep -c 'Overfull \\hbox' main.log       # 2 known, listed in Status; a 3rd is new
grep -o 'Output written.*' main.log      # page count
```

Underfull vboxes are not a defect here: every one is `while \output is
active`, which is page-breaking around floats. An **overfull hbox** is, since it
puts ink in the margin; the two that survive are recorded in Status.

### Two checks the build cannot make

A clean build is not evidence of either of these, because both fail silently.
Run them after any pass that moves a figure or renumbers a chapter.

**Every included figure is in the repository.** `.gitignore` carries `*.pdf`, so
a figure that was never force-added is invisible to `git status` and present
only on the machine that drew it — the book builds there and nowhere else:

```sh
cd book
for f in $(grep -rho 'fig-[a-z-]*\.pdf' *.tex | sort -u); do
    git ls-files --error-unmatch "figs/$f" >/dev/null 2>&1 ||
        echo "UNTRACKED (but included): figs/$f"
done
```

Reverse it to find orphans — figures still generated but included nowhere — by
looping over `figs/fig-*.pdf` and grepping the `.tex` for each. Both directions
should print nothing.

**Every section reference in the code still resolves.** The scripts cite the
book by number (`Sec. 15.4`, `Chapter 16`) rather than by label, so nothing
breaks when a chapter is renumbered; the comments just start lying. This
compares them against the numbers LaTeX actually assigned:

```sh
cd book && python3 - <<'EOF'
import re, glob, pathlib
sec, cha = {}, {}
for f in glob.glob('*.aux'):
    txt = open(f).read()
    for m in re.finditer(r'\\newlabel\{[^}]*\}\{\{(\d+(?:\.\d+)*)\}\{\d+\}\{(.*?)\}\{(section|chapter|subsection)\.', txt):
        (cha if m.group(3) == 'chapter' else sec)[m.group(1)] = m.group(2)[:52]
    for m in re.finditer(r'\\contentsline \{(subsection|section)\}\{\\numberline \{(\d+(?:\.\d+)*)\}(.*?)\}\{', txt):
        sec.setdefault(m.group(2), m.group(3)[:52])
for r in ('figs', '../statmech/probe', '../statmech/src', '../statmech/tests',
          '../statmech/examples', '../percolation/src', '../percolation/tests'):
    for p in sorted(pathlib.Path(r).rglob('*.py')):
        if '__pycache__' in str(p):
            continue
        for i, line in enumerate(p.read_text(errors='ignore').split('\n'), 1):
            for m in re.finditer(r'\b(?:Sec\.?|Section)s?\s+(\d+\.\d+(?:\.\d+)?)', line):
                if m.group(1) not in sec:
                    print(f'{p}:{i}  Sec. {m.group(1)}  {line.strip()[:70]}')
            for m in re.finditer(r'\b(?:Ch\.?|Chapter)s?\s+(\d+)\b', line):
                if m.group(1) not in cha:
                    print(f'{p}:{i}  Ch. {m.group(1)}  {line.strip()[:70]}')
EOF
```

It should print nothing. Two caveats learned the hard way: index **subsections**
as well as sections, or `Sec. 5.6.2` reads as broken, and match the title
non-greedily, or `Sec. 9.4` does — its title contains `\texorpdfstring`. And
resolving is weaker than being right: a reference can point at a section that
exists and is the wrong one, which is how `Sec. 14.2's claim` survived in
`overlap.py` after 14.2 became a different section. The audit narrows the
reading; it does not do it.

Build into a scratch directory --- `latexmk -pdf -interaction=nonstopmode
-outdir=/tmp/bookbuild main.tex` --- when the point is to verify rather than to
refresh `main.pdf`; the output is identical and the working tree is left alone.

The counts in Status are re-measured, not maintained by hand. The reliable
route is the `.aux` files, which carry what LaTeX actually numbered:
`grep -h newlabel *.aux | grep -c 'figure\.'` and likewise `table\.` and
`equation\.`; index terms are `grep -c '^  \\item' main.ind` and references
`grep -c '\\bibitem' main.bbl`. Counting `\begin{figure}` and friends over
`*.tex` overcounts, since not every float is labelled.
`\ref`/`\eqref` resolution and `\code`-name resolution are worth scripting
against both repositories, remembering that `joint.JointGiantComponent` is an
alias rather than a `def` and that `\allowbreak` has to be stripped out of a
`\code` argument before it is looked up.

Two environment notes. `latexmk` will hang rather than stop on a TikZ error even
under `-interaction=nonstopmode`; add `-halt-on-error`, or run `pdflatex`
directly, when a figure is new. GNU `timeout` is not on the path here.

To inspect a page at trim size, `pdftoppm -png -r 120 -f N -l M main.pdf out`
and read the PNG. **Always do this for a new figure** — compiling is not
evidence that a figure is legible at 5×8.

### Figure scripts

All live in `figs/` and write PDFs and `.tex` tables into that directory, as in
*Local network growth*. Rerun a script after editing it; the build picks up the
new PDF automatically.

`figs/real_chygraphs.py` generates
Table 3.2 and Figure 3.2 from the netzschleuder dumps cached under
`~/av2atg/LocalNetworkGrowth/figs/data/netzschleuder`;
`figs/triangles_vs_regular.py` generates Figure 4.2 and verifies
$\bar s_\triangle(q)$ by exact enumeration;
`figs/beyond_threshold.py` generates Figures 5.1–5.3 by importing the
`percolation` package from `../percolation/src`;
`figs/epidemics.py` generates Figures 6.2 and 6.3 and runs every check quoted in
Ch. 6 — that `Gbar^1_0` at `n = 3` is the bond-percolated triangle of Ch. 4,
that `theta + 1 = R*` for several household size distributions, the
fixed-budget comparison, and the tricritical condition and exponent of
Sec. 6.4. Figure 6.1 is TikZ in `epidemics.tex`.
`figs/potts.py` generates Figure 7.2 and runs Ch. 7's checks: the `q -> 1` limit
of the Potts interior sum against `Gbar` (cliques to `K_5`, paths, distinct `y_j`
per member, symbolic in `p`), the `q = 2` limit against the Ising cavity field of
`chygraph_statmech.cavity`, the transmission factors and their two limits, the
linear threshold, and the order of the transition. It takes a couple of minutes —
the `K_5` interior sum is the slow part. Figure 7.1 is TikZ in `potts.tex`.
`figs/statmech.py` generates Figure 8.2 and runs Ch. 8's checks: that `u'` is
Ch. 7's `tau_c` at `q = 2`; the size-biased average three ways (mixed layer,
split into one layer per cardinality — identical — and the mean-cardinality
substitution, wrong by 10.1%); the links-and-triangles determinant for Poisson
and regular layers, so the cross term is exercised both where it vanishes and
where it does not; the Bethe free energy against its closed form and the
textbook graph result; and the exchange of stability, taken along a sequence in
`Lambda` and stopped before rounding dominates. It imports `figs/potts.py`, so
the two live together. Figure 8.1 is TikZ in `statmech.tex`.
`figs/ising.py` generates Figures 9.1–9.3 and runs Ch. 9's checks: the closed
forms for `u'` at c = 2, 3, 4 against the enumeration; the fixed-degree
clustering family (monotone in f at every degree) and the regular-null table;
the AT line; the unanimity `u'`, the Bragg–Williams limit and its `C_q`
residual including the vanishing at q = 4; and the tricritical boundary. The
Monte Carlo points in Figure 9.2 are read from the cached
`../statmech/probe/results/ising_mc.log` rather than recomputed — that run is the one
test in the chapter from outside the formalism. `figure_unanimity` and
`check_tricritical` each solve coexistence windows and take a minute or two.
`figs/hittingset.py` generates Figure 10.2 and runs Ch. 10's checks: the
hard-field threshold `<k>(c-1) = e` to eleven digits for c = 2..20; the
disjoint-3-hyperedge counterexample (hard 0.5, soft 0.3333, truth 1/3);
Weigt–Hartmann at cardinality two; the regular closed forms `h_RS`, `rho = 1/K`
and the entropy; the mixed-cardinality and correlation scans with the
`isolated_fraction` confound reported alongside; and the entropy criterion.
Figure 10.1 is TikZ in `hittingset.tex`. Soft-field densities come from
population dynamics and carry a scatter of a few parts in a thousand — the
script averages over seeds and prints the spread, and the book quotes those
numbers to three decimals only.
`figs/cover.py` generates Figures 11.2 and 11.3 and runs Ch. 11's checks: the
`<k> = e` identification from core percolation and from Ch. 10's hard-field map
(agreeing to 2.5e-10); the induced-graph cover against Weigt–Hartmann at
cardinality two and its 1/2-vs-2/3 failure on isolated triangles; the
factorisation `f(delta) = -(1-delta)^2 sum (j+1) delta^j` symbolically for
c = 2..9; that the core is exactly `1 - e^{-<k>}` above cardinality two; and the
cached leaf-removal validation. Figure 11.1 is TikZ in `cover.tex`. The
hyperbolic-graph points are read from `../statmech/probe/results/prediction4.csv` and the
simulation check from `../statmech/probe/results/validate_core.txt`; neither is
recomputed by the book.
`figs/cover.py` also writes **Table 11.1** (`table_real_core`) from
`../statmech/probe/real_core.json`: leaf removal on sixteen real networks --
Ch. 3's ten and Ch. 14's six, from the same netzschleuder cache -- against
Eq. (11.4) on each graph's own clique ensemble and against a degree-matched
control averaged over twenty rewirings. `statmech/probe/real_core.py` produces
the JSON and takes about two minutes; it needs `numba`, through
`computational_complexity/code/leafremoval.py`.
`figs/colouring.py` generates Figure 12.1 and runs Ch. 12's checks: both closed
forms against exact enumeration in rational arithmetic; that the linearised map
is a scalar on the traceless subspace (three directions, one answer); the graph
case against Zdeborová & Krzakala Eq. (18); and the comparison with Mulet's
published thresholds that shows the stability line is not the colourability
line. All of that is fast — exact rational arithmetic over small complexes.
The **survey-propagation sections** are not. Three bisections
of the complexity over three seeds each — Sec. 12.6's `c_q`, the clustering
by-product, and Sec. 12.7's hypergraph thresholds against Gabrié et al. — put
the whole script at about **fourteen minutes**, which makes it the slowest in
the book.
`figs/satisfiability.py` generates Figure 13.1 and runs Ch. 13's checks: the
clause interior against enumeration in rational arithmetic (including the
non-uniform emitted message that rules out a symmetric fixed point); the `k = 2`
linearisation and the branch appearing continuously at `alpha = 1`; that the
`k >= 3` derivative is *exactly* zero, tested by the rate at which the finite
difference falls (`eps^(k-2)`, checked by halving `eps`) rather than by an
absolute bound, which would only show `eps` was small; and both conventions for
a contradicted variable. It also carries the **survey-propagation section**:
the SP update, the complexity `Sigma` as a Bethe
count of clusters at `m = 0`, and `Sigma = 0` bisected against the published
`alpha_s`. That check averages over three seeds and takes about five minutes,
which makes this the second slowest script in the book.
`figs/overlap.py` plots **all of Part IV** — it is the one script serving three
chapters. Figures 14.2, 14.3, 15.2–15.4 and 16.3–16.6, from the cached probe
outputs (`gbp_cliques.json`, `gbp_real.py`, `merge_lnz.py`, `core_fraction.py`,
`karrer_core_sweep.py`); Figures 14.1, 15.1, 16.1 and 16.2 are TikZ in the
chapters themselves. It also runs the checks: the Möbius and Bethe counting
numbers on two overlapping triangles, including the factor-coverage test that is
the whole point (shared bond covered once against twice); the two-triangle
table, recomputed here since it is a four-spin enumeration; the 60-instance
summary from `../statmech/probe/results/gbp_cliques.json` with the convergence threshold
stated; and the clique-ensemble paired ratio read from
`../statmech/probe/results/analysis.txt`.
`figs/merge.py` **writes no figures** — it says so at the top of the file. What
it produces are numbers: the finite-size sweep behind Table 15.2
(`check_placed_finite_size`), the merge closure on the six real networks, and
the two-triangle and diamond checks. It runs in about 17 s. It needs
`PYTHONPATH` to carry both `src` trees, or it fails at import.

## Conventions

- Drawing vocabulary is fixed once, in `main.tex`: `nd` (node), `ndf` (filled
  node), `hub`, `cx`/`cxb` (a complex drawn around a set of nodes), `lnk`,
  `msg` (a message along a directed inclusion), `lb`/`tb` (labels), and `\cxg{c}`
  for an inline glyph of a cardinality-$c$ complex. Every figure is built from
  these, so a reader who learns the first schematic can read the last one
  without its caption.
- Notation matches the papers (`\ave{}`, `\kbar`, `\sbar`), so a reader can move
  between book and paper without retranslating.
- **No `calculation` boxes.** A derivation the reader is told to skip is a
  derivation they cannot follow when they want it. Derivations run in continuous
  text, where they are needed, and **each ends by naming the routine that
  executes it** --- `cover.CliqueCover.solve` and the like --- so the path from
  the algebra to the code is in the sentence rather than only in
  `software.tex`'s Table 1. The `calculation` environment and the `tcolorbox`
  dependency are gone, so a new box will not compile.
- **Long `\code` names carry `\allowbreak` at their separators.** `\code` is a
  bare `\texttt` with no break points, and a name like
  `statmech/probe/karrer_core_sweep.core_closed_form` overflows the 4.05in
  measure by 232pt without them. The rule is applied everywhere, including
  `software.tex`'s table, and is idempotent: insert after `.`, `/` and `\_`
  where one is not already there.
- Overlapping complexes are drawn as two translucent boxes of the same grey, so
  the shared atoms show as a darker patch. Established in Fig. 2.7 and used
  wherever overlap matters.
- **Section and subsection headings are ragged right**, patched in `main.tex`
  from `book.cls`'s own definitions with nothing changed but the addition of
  `\raggedright`. The class sets them justified, which on this 4.05in measure
  makes a long title run into the margin instead of breaking — renaming Sec. 1.6
  to "Replica symmetry and symmetry breaking" overflowed by 15.5pt. Chapter
  heads were already ragged. Do not revert this to get "tidier" headings; a
  single-line heading looks identical either way.
- **Contents lists parts and chapters only.** `\setcounter{tocdepth}{0}` sits
  immediately before `\tableofcontents` — it must go there, not in the preamble,
  where the packages override it. (The `.toc` *file* always contains every
  section; the filtering happens when the contents is typeset, so do not "fix"
  this by reading `main.toc`.)
- **The index is one ragged-right column.** `theindex` is redefined in
  `main.tex`; the class default is two columns, which on a 5-inch page leaves
  ~1.7in of measure and overflows on entries like "Fortuin--Kasteleyn
  correspondence".
- Figures must be checked at the 5×8 trim, not just for compilation: more than
  three panels side by side is illegible on this page width. Stack into rows.
