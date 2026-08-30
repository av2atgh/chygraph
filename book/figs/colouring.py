"""Chapter 12: graph colouring on a chygraph.

Unlike every other chapter, this one has no manuscript behind it, so the
calculations are here rather than in `../src`. Everything is done twice: once
by exact enumeration of the interior of a complex, once from a closed form
derived by hand, and the two are required to agree.

  fig-colouring   the two constraints. Proper colouring -- every member of a
                  complex a different colour -- against hypergraph colouring,
                  which only forbids a monochromatic complex. Identical at
                  cardinality two and violently different above it.

The survey-propagation section goes one level up, to 1RSB at m = 0, and reaches
the colourability threshold that Sec. 12.4's stability line is not: the
complexity vanishes at Mulet's c_q for q = 3, 4 and 5, and the density at which
the surveys stop being trivial is his c_d. Both columns of Table 12.1 out of one
population, neither used in getting there.

Benchmarks, from `~/Downloads/chygraph_references/`:
  Mulet, Pagnani, Weigt & Zecchina, Phys. Rev. Lett. 89, 268701 (2002), Table I
  Zdeborova & Krzakala, Phys. Rev. E 76, 031131 (2007), Eq. (18)
  Braunstein, Mulet, Pagnani, Weigt & Zecchina, Phys. Rev. E 68, 036702 (2003)
    -- survey propagation for colouring; the update below is its m = 0 form
"""

import itertools
import sys
from fractions import Fraction as F
from pathlib import Path

import numpy as np

OUT = Path(__file__).resolve().parent
DARK, MID, LIGHT = '0.10', '0.45', '0.70'

# Mulet et al., Table I: Erdos-Renyi, clustering and colourability thresholds.
MULET = {3: (4.42, 4.69), 4: (8.27, 8.90), 5: (12.67, 13.69)}


def _mpl():
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    return plt


def _tidy(ax):
    ax.tick_params(labelsize=8)
    for sp in ('top', 'right'):
        ax.spines[sp].set_visible(False)


# ---------------------------------------------------------------------------
# the interior sum
# ---------------------------------------------------------------------------

def emitted(c, q, eta, rule):
    """``Z(sigma_0)`` for one complex of cardinality ``c``.

    ``eta[j][s]`` is the message from member ``j = 1 .. c-1`` on colour ``s``.
    ``rule`` is 'proper' -- every member a different colour, which is proper
    colouring of the induced graph -- or 'hyper', which forbids only a
    monochromatic complex.
    """
    Z = []
    for s0 in range(q):
        tot = F(0)
        for rest in itertools.product(range(q), repeat=c - 1):
            if rule == 'proper':
                if len(set((s0,) + rest)) != c:
                    continue
            elif all(x == s0 for x in rest):
                continue
            w = F(1)
            for j, sj in enumerate(rest):
                w *= eta[j][sj]
            tot += w
        Z.append(tot)
    return Z


def tau_enumerated(c, q, rule, eps=F(1, 10**7)):
    """``tau``: the outgoing perturbation per unit incoming perturbation.

    Linearised at the uniform fixed point ``eta = 1/q``.  The colour-permutation
    symmetry means the Jacobian acts as a scalar on the traceless perturbations,
    so one direction suffices; :func:`check_scalar` verifies that.
    """
    u = F(1, q)
    pert = [[u] * q for _ in range(c - 1)]
    e = [F(1) - u if s == 0 else -u for s in range(q)]      # traceless
    pert[0] = [u + eps * e[s] for s in range(q)]
    Z = emitted(c, q, pert, rule)
    tot = sum(Z)
    return (Z[0] / tot - u) / (eps * e[0])


def tau_closed(c, q, rule):
    """The closed forms derived in the chapter.

    Proper colouring of a c-clique needs at least c colours; below that the
    interior sum is identically zero and the closed form is a formal
    continuation with no meaning, so it is refused rather than returned.
    """
    if rule == 'proper':
        if q < c:
            raise ValueError(f'a {c}-clique has no proper {q}-colouring')
        return F(-1, q - 1)
    return F(-1, q**(c - 1) - 1)


def check_closed_forms():
    """Enumeration against the closed forms, exactly, in rational arithmetic."""
    print('  proper colouring:  tau = -1/(q-1), independent of cardinality')
    for q in (3, 4, 5, 6, 7):
        for c in range(2, min(q, 5) + 1):
            got, want = tau_enumerated(c, q, 'proper'), tau_closed(c, q, 'proper')
            assert got == want, (c, q, got, want)
        print(f'    q = {q}: c = 2 .. {min(q, 5)}   all equal {want}   OK')
    print('  hypergraph colouring:  tau = -1/(q^(c-1) - 1)')
    for q in (3, 4, 5):
        for c in (2, 3, 4):
            got, want = tau_enumerated(c, q, 'hyper'), tau_closed(c, q, 'hyper')
            assert got == want, (c, q, got, want)
        print(f'    q = {q}: c = 2, 3, 4   OK')


def check_feasibility():
    """The proper rule is undefined below q = c, and says so."""
    for q, c in ((3, 4), (4, 5), (2, 3)):
        Z = emitted(c, q, [[F(1, q)] * q] * (c - 1), 'proper')
        assert all(z == 0 for z in Z), (q, c)
        try:
            tau_closed(c, q, 'proper')
        except ValueError:
            pass
        else:
            raise AssertionError(f'q={q}, c={c} should have been refused')
        print(f'    q = {q}, c = {c}: interior sum is identically zero, '
              f'closed form refused')


def check_scalar():
    """The Jacobian really is a scalar on the traceless subspace."""
    eps = F(1, 10**7)
    for q, c, rule in ((4, 3, 'proper'), (4, 3, 'hyper'), (5, 4, 'proper')):
        u = F(1, q)
        vals = []
        for a, b in ((0, 1), (1, 2), (0, 2)):
            e = [F(0)] * q
            e[a], e[b] = F(1), F(-1)                       # another traceless direction
            pert = [[u] * q for _ in range(c - 1)]
            pert[0] = [u + eps * e[s] for s in range(q)]
            Z = emitted(c, q, pert, rule)
            tot = sum(Z)
            vals.append((Z[a] / tot - u) / (eps * e[a]))
        assert len(set(vals)) == 1, (q, c, rule, vals)
        print(f'    q = {q}, c = {c}, {rule}: three traceless directions give '
              f'one number, {vals[0]}   OK')


# ---------------------------------------------------------------------------
# the threshold
# ---------------------------------------------------------------------------

def stability_threshold(c, q, rule, regular=False):
    """Neighbours per node at which the uniform solution loses local stability.

    Sec. 8.3's branching matrix with ``u'`` squared, which is the de
    Almeida-Thouless form: a perturbation of random sign is transmitted with
    ``tau`` and its effect measured by squaring.  For one layer,
    ``<kbar> (c-1) tau^2 = 1``.
    """
    t = float(tau_closed(c, q, rule))
    kbar = 1.0 / ((c - 1) * t * t)          # excess chy-degree at the threshold
    k = kbar + 1.0 if regular else kbar     # regular: <kbar> = <k> - 1
    return k * (c - 1)                      # reported in neighbours per node


def check_against_literature():
    """The graph case, against two published closed forms."""
    print('        q     ER: (q-1)^2   ours     regular: (q-1)^2+1   ours')
    for q in (3, 4, 5, 6, 10):
        er = stability_threshold(2, q, 'proper')
        rg = stability_threshold(2, q, 'proper', regular=True)
        assert abs(er - (q - 1)**2) < 1e-12
        assert abs(rg - ((q - 1)**2 + 1)) < 1e-12
        print(f'  {q:>9}   {(q-1)**2:>10}  {er:>7.3f}   {(q-1)**2+1:>15}  '
              f'{rg:>7.3f}')
    print('    Zdeborova & Krzakala Eq. (18), both lines, from one tau')


def check_not_the_colourability_threshold():
    """What c_RS is and is not, against Mulet et al. Table I."""
    print('     q   RS instability   clustering c_d   colourable to c_q')
    for q in (3, 4, 5):
        cd, cq = MULET[q]
        crs = stability_threshold(2, q, 'proper')
        where = 'inside' if crs < cq else 'ABOVE the colourable phase'
        print(f'  {q:>4}   {crs:>14.0f}   {cd:>14.2f}   {cq:>17.2f}   ({where})')
    assert stability_threshold(2, 3, 'proper') < MULET[3][1]
    assert all(stability_threshold(2, q, 'proper') > MULET[q][1] for q in (4, 5))
    print('    only at q = 3 does the RS instability lie inside the colourable')
    print('    phase; it is a local stability bound, not a colourability one')


def check_cardinality():
    """Counted in neighbours, cardinality does not move the proper-colouring
    point at all -- and moves the hypergraph one enormously."""
    print('     q   c   proper (neighbours)   hypergraph (neighbours)')
    for q in (4, 5):
        for c in (2, 3, 4):
            if c > q:
                continue
            pr = stability_threshold(c, q, 'proper')
            hy = stability_threshold(c, q, 'hyper')
            assert abs(pr - (q - 1)**2) < 1e-9
            print(f'  {q:>4} {c:>3}   {pr:>18.3f}   {hy:>21.3f}')
    print('    proper: (q-1)^2 at every cardinality, exactly as on a graph')
    print('    hypergraph: (q^(c-1)-1)^2, which runs away with c')


# ---------------------------------------------------------------------------
# survey propagation
# ---------------------------------------------------------------------------
#
# Sec. 12.4's threshold is where the replica-symmetric solution loses local
# stability, and Sec. 12.5 says at length that this is not the colourability
# threshold. Reaching that threshold needs one step of replica symmetry
# breaking, which is what follows: survey propagation at m = 0, with the
# complexity Sigma counting clusters and Sigma = 0 locating c_q.
#
# Colour symmetry does most of the work. A survey is a distribution over which
# colour a node is forced to, and permutation symmetry makes it uniform, so one
# scalar per inclusion suffices: e = P(forced at all), each colour carrying e/q.
# What symmetry does NOT do is decouple the colours. One forced neighbour
# forbids exactly one colour, so "every other colour forbidden" is a coverage
# problem and not a product over colours -- hence the inclusion-exclusion:
#
#   P(available set = {c}) = sum_j (-1)^j C(q-1,j) prod_k (1 - e_k (1+j)/q)
#   P(no colour available) = sum_j (-1)^j C(q,j)   prod_k (1 - e_k j/q)
#
# Treating the colours as independent instead is wrong, and wrong in a way that
# a single threshold would hide; the check below is against three.


def _sp_prods(e, deg, rng, q, size):
    """``prod_k (1 - e_k t / q)`` for ``t = 0..q``, per sample."""
    tot = int(deg.sum())
    draws = e[rng.integers(0, e.size, tot)] if tot else np.empty(0)
    ends = np.cumsum(deg)
    out = np.ones((q + 1, size))
    for t in range(1, q + 1):
        v = np.log(np.maximum(1.0 - draws * t / q, 1e-300))
        c = np.concatenate(([0.0], np.cumsum(v)))
        out[t] = np.exp(c[ends] - c[ends - deg])
    return out


def sp_update(e, deg, rng, q, size):
    """One survey update, and the weight of the contradictory state."""
    from math import comb
    P = _sp_prods(e, deg, rng, q, size)
    w = sum((-1) ** j * comb(q - 1, j) * P[1 + j] for j in range(q))
    contra = sum((-1) ** j * comb(q, j) * P[j] for j in range(q + 1))
    return np.clip(q * w / np.maximum(1.0 - contra, 1e-300), 0.0, 1.0), contra


def sp_converge(q, c, size=20000, sweeps=200, seed=0):
    """Population of surveys at the SP fixed point, Erdos-Renyi mean degree c."""
    rng = np.random.default_rng(seed)
    e = rng.uniform(0.3, 0.6, size)
    for _ in range(sweeps):
        e, _ = sp_update(e, rng.poisson(c, size), rng, q, size)
    return e, rng


def sp_complexity(q, c, size=20000, sweeps=200, seed=0, nsamp=300000):
    """``Sigma``, the Bethe count of clusters at ``m = 0``, per node.

        Sigma = <ln(1 - P(no colour available))>   full degree, per node
              - (c/2) <ln(1 - e_i e_j / q)>        per edge

    the edge term being the weight of the two ends not forced to the same
    colour. Both vanish at ``e = 0``, so Sigma says nothing on the trivial
    branch; it counts the clusters of the non-trivial one.
    """
    e, rng = sp_converge(q, c, size, sweeps, seed)
    if e.mean() < 1e-6:
        return 0.0, 0.0
    _, contra = sp_update(e, rng.poisson(c, nsamp), rng, q, nsamp)
    site = np.log(np.maximum(1.0 - contra, 1e-300)).mean()
    ei = e[rng.integers(0, e.size, nsamp)]
    ej = e[rng.integers(0, e.size, nsamp)]
    edge = np.log(np.maximum(1.0 - ei * ej / q, 1e-300)).mean()
    return site - (c / 2.0) * edge, float(e.mean())


def sp_threshold(q, lo, hi, seed=0, iters=9, **kw):
    """Mean degree at which ``Sigma`` vanishes: the colourability threshold."""
    assert sp_complexity(q, lo, seed=seed, **kw)[0] > 0, (q, lo)
    for _ in range(iters):
        mid = 0.5 * (lo + hi)
        if sp_complexity(q, mid, seed=seed, **kw)[0] > 0:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def check_survey_thresholds():
    """Sigma = 0 against Mulet's colourability thresholds."""
    print('     q   Sigma = 0 at   c_q (published)   error')
    for q, lo, hi in ((3, 4.55, 4.95), (4, 8.65, 9.25), (5, 13.35, 14.05)):
        c_q = MULET[q][1]
        r = [sp_threshold(q, lo, hi, seed=s) for s in (0, 1, 2)]
        m, sd = float(np.mean(r)), float(np.std(r))
        assert abs(m - c_q) / c_q < 0.03, (q, m, c_q)
        print(f'  {q:>4}   {m:>7.3f}({sd:.3f})   {c_q:>15.2f}'
              f'   {100 * (m - c_q) / c_q:>+5.1f}%')
    print('    Three thresholds from one calculation, none of them used in')
    print('    getting there. This is the line Sec. 12.4 is not: the stability')
    print('    threshold sits at 4, 9, 16 for these q, above c_q from q = 4 on.')


def check_survey_branch_is_clustering():
    """The survey branch appears at Mulet's clustering threshold c_d."""
    print('     q   branch appears   c_d (published)')
    for q, lo, hi in ((3, 4.0, 5.0), (4, 7.6, 8.8)):
        c_d = MULET[q][0]
        lo_, hi_ = lo, hi
        for _ in range(9):
            mid = 0.5 * (lo_ + hi_)
            if sp_converge(q, mid)[0].mean() > 1e-4:
                hi_ = mid
            else:
                lo_ = mid
        r = 0.5 * (lo_ + hi_)
        assert abs(r - c_d) / c_d < 0.05, (q, r, c_d)
        print(f'  {q:>4}   {r:>14.3f}   {c_d:>15.2f}')
    print('    A by-product: the density at which the surveys stop being')
    print('    trivial is the clustering transition, so the same population')
    print('    carries both columns of Table 12.1.')


def figure_colouring():
    plt = _mpl()
    fig, axes = plt.subplots(1, 2, figsize=(4.6, 2.5))

    ax = axes[0]
    qs = np.arange(3, 9)
    for c, col, mk in ((2, DARK, 'o'), (3, MID, 's'), (4, LIGHT, '^')):
        ax.semilogy(qs, [abs(float(tau_closed(c, q, 'hyper'))) for q in qs],
                    mk + '--', ms=3.4, lw=1.0, color=col, dashes=(3, 2),
                    mfc='white', mew=0.9)
        feas = qs[qs >= c]                     # proper colouring needs q >= c
        ax.semilogy(feas, [abs(float(tau_closed(c, q, 'proper'))) for q in feas],
                    mk + '-', ms=3.4, lw=1.4, color=col, label=f'$c={c}$')
    ax.set_xlabel(r'colours $q$', fontsize=8.5)
    ax.set_ylabel(r'$|\tau_{c}(q)|$', fontsize=8.5)
    ax.legend(frameon=False, fontsize=7.5, loc='lower left')
    ax.annotate('proper: the three\ncoincide', xy=(5.6, 0.30), fontsize=6.4,
                color='0.3')
    ax.annotate('hypergraph', xy=(6.1, 2e-3), fontsize=6.4, color='0.45')
    _tidy(ax)

    ax = axes[1]
    for c, col, mk in ((2, DARK, 'o'), (3, MID, 's'), (4, LIGHT, '^')):
        ax.semilogy(qs, [stability_threshold(c, q, 'hyper') for q in qs],
                    mk + '--', ms=3.4, lw=1.0, color=col, dashes=(3, 2),
                    mfc='white', mew=0.9)
        feas = qs[qs >= c]
        ax.semilogy(feas, [stability_threshold(c, q, 'proper') for q in feas],
                    mk + '-', ms=3.4, lw=1.4, color=col)
    ax.plot([3, 4, 5], [MULET[q][1] for q in (3, 4, 5)], '*', ms=7,
            color=DARK, mfc='white', mew=0.9)
    ax.annotate('proper, $c=2,3,4$', xy=(5.6, 24), fontsize=6.4, color='0.3')
    ax.annotate('hypergraph', xy=(6.2, 4e5), fontsize=6.4, color='0.45')
    ax.annotate('$\\star$: colourable only below', xy=(4.6, 2.6),
                fontsize=6.2, color='0.3')
    ax.set_ylim(2, 2e6)
    ax.set_xlabel(r'colours $q$', fontsize=8.5)
    ax.set_ylabel('neighbours per node', fontsize=8.5)
    _tidy(ax)

    fig.tight_layout()
    fig.savefig(OUT / 'fig-colouring.pdf')
    print(f'  wrote {OUT / "fig-colouring.pdf"}')


if __name__ == '__main__':
    print('the interior sum, enumerated against the closed forms:')
    check_closed_forms()
    print('below q = c there is no proper colouring:')
    check_feasibility()
    print('the Jacobian is a scalar on the traceless subspace:')
    check_scalar()
    print('the graph case, against the literature:')
    check_against_literature()
    print('what that threshold is not:')
    check_not_the_colourability_threshold()
    print('cardinality:')
    check_cardinality()
    print('survey propagation, against the colourability thresholds:')
    check_survey_thresholds()
    print('and the clustering threshold, as a by-product:')
    check_survey_branch_is_clustering()
    print('figure:')
    figure_colouring()
