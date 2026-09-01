# Finite-y 1RSB on a chygraph

**Status: the extension is verified, the numbers it produces are not gated.**
The reduction to the hard constraint is exact and lives as a check in
`book/figs/colouring.py`; the finite-`y` plane below has no external benchmark
and is therefore here rather than in the book.

`finite_y_chygraph.py`. Imports the machinery from `book/figs/colouring.py`
rather than duplicating it.

## The extension, which is a substitution

At reweighting `y` each cluster carries `e^{-yE}`, `E` the number of violated
complexes. In the cavity recursion the energy shift of adding one complex is 1
if that complex is violated, so the interior sum of Eq. (8.4) runs over **all**
configurations with the violating ones weighted `e^{-y}` instead of excluded.

That is Ch. 8's `-beta H` taking a different argument. The chy-degree
convolution, the survey structure and the three-term complexity are unchanged:

```
y -> infinity   the hard constraint, and every result in Secs. 12.6-12.9
finite y        metastable states at energy density e(y)
```

Reference: Krzakala, Pagnani & Weigt, *Phys. Rev. E* **70**, 046705 (2004), which
sets out the apparatus for graphs — reweighting `y`, complexity `Sigma(y)`,
energy `e(y)`, `y = infinity` as the ground-state limit. Sec. 12.3 already
borrows its stability criterion.

## What is checked

`y -> infinity` reproduces the hard-constraint results **exactly, to every
digit**, at `c = 2` and `c = 3`. That is the check that the substitution is the
right one, and it is in `figs/colouring.py` where the rest of the gated work
lives.

## What is not checked

Every finite-`y` number below. Krzakala et al. give `Sigma(y)` for graphs and
Gabrié et al. give `m = 1` condensation and rigidity thresholds; either would
gate this at `c = 2` before it is turned on triangles. **That has not been
done**, and until it is, nothing in the plane should be quoted.

## The plane, for triangles (c = 3)

`Sigma` against chy-degree and reweighting. Degree is `2 kappa`.

**q = 4**

| kappa | deg | y=1.5 | y=2 | y=3 | y=4 | y=6 | y=10 | y=inf |
|---|---|---|---|---|---|---|---|---|
| 3.50 | 7.0 | −0.752 | −0.326 | **+0.091** | **+0.136** | −0.086 | −0.271 | −0.281 |
| 3.80 | 7.6 | −2.057 | −1.167 | **+0.011** | **+0.488** | **+0.340** | −0.221 | −0.293 |
| 4.20 | 8.4 | −4.299 | −2.914 | −0.748 | **+0.505** | **+0.938** | −0.075 | −0.374 |
| 4.60 | 9.2 | −6.743 | −4.970 | −1.981 | **+0.039** | **+1.360** | **+0.200** | −0.445 |
| 5.00 | 10.0 | −9.577 | −7.487 | −3.769 | −0.952 | **+1.528** | **+0.598** | −0.545 |
| 5.50 | 11.0 | −13.513 | −11.053 | −6.479 | −2.677 | **+1.526** | **+1.345** | −0.644 |

**q = 3**

| kappa | deg | y=1.5 | y=2 | y=3 | y=4 | y=6 | y=10 | y=inf |
|---|---|---|---|---|---|---|---|---|
| 1.60 | 3.2 | −0.160 | −0.108 | −0.079 | −0.097 | −0.141 | −0.159 | −0.160 |
| 1.80 | 3.6 | −0.396 | −0.192 | **+0.015** | **+0.035** | −0.108 | −0.257 | −0.270 |
| 2.00 | 4.0 | −0.649 | −0.276 | **+0.188** | **+0.336** | **+0.171** | −0.187 | −0.248 |
| 2.30 | 4.6 | −1.814 | −1.142 | −0.134 | **+0.404** | **+0.507** | −0.139 | −0.363 |
| 2.60 | 5.2 | −3.262 | −2.332 | −0.769 | **+0.303** | **+1.034** | **+0.192** | −0.425 |
| 3.00 | 6.0 | −6.170 | −4.888 | −2.568 | −0.732 | **+1.155** | **+0.750** | −0.461 |

## Two readings, and only the first is safe

**Sec. 12.9 survives the test that could have killed it.** The `y = inf` column
is negative at *every* `kappa`, for both `q`, and grows *more* negative with
`kappa` rather than approaching zero from below. So the claim there — that
promoting triangles to complexes leaves no window at the ground-state endpoint
— is not an artefact of the `kappa` values sampled. This is the strongest thing
the plane says and it is a negative result, which is why it can be trusted
without a finite-`y` benchmark: it concerns the `y = inf` column, which *is*
gated.

**The positive ridge is not a colourability threshold**, and should not be
written up as one. It sits at intermediate `y`, widening with `kappa`, and
describes metastable states at finite energy density. Krzakala et al. locate
`c_q` at `y = infinity`, which is exactly where the ridge is absent. Reading a
threshold off the ridge would be quoting the wrong endpoint of the family.

## What would finish it

Gate the finite-`y` machinery at `c = 2` against Krzakala's `Sigma(y)` or
Gabrié's `m = 1` thresholds, then extract `e(y)` alongside `Sigma(y)` so the
ridge can be placed on the energy axis rather than the reweighting axis. Only
then is there a statement to make about what the metastable states of a
triangle-clustered graph are.
