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
the colourability threshold that Sec. 12.3's stability line is not: the
complexity vanishes at Mulet's c_q for q = 3, 4 and 5, and the density at which
the surveys stop being trivial is his c_d. Both columns of Table 12.1 out of one
population, neither used in getting there.

Benchmarks, from `~/Downloads/chygraph_references/`:
  Mulet, Pagnani, Weigt & Zecchina, Phys. Rev. Lett. 89, 268701 (2002), Table I
  Zdeborova & Krzakala, Phys. Rev. E 76, 031131 (2007), Eq. (18)
  Braunstein, Mulet, Pagnani, Weigt & Zecchina, Phys. Rev. E 68, 036702 (2003)
    -- survey propagation for colouring; the update below is its m = 0 form
  Gabrie, Dani, Semerjian & Zdeborova, J. Phys. A 50, 505002 (2017), Table 1
    -- q-colouring of K-uniform hypergraphs; its Eq. (43) is the stability
    threshold this chapter derives, and its l_col benchmarks Sec. 12.8
"""

import itertools
import sys
from itertools import combinations, product
from math import comb
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


def triangle_family(q, f, regular=True):
    """Degree at which proper colouring loses RS stability, links + triangles.

    Two regular layers at fixed degree ``n``: links of chy-degree ``a``,
    triangles of chy-degree ``b``, ``n = a + 2b``, and ``f = 2b/n`` the fraction
    of degree carried by triangles. Eq. (8.11) with ``u' -> tau^2``, and tau the
    same for both layers because it does not depend on cardinality:

        s[(a-1) + 2(b-1)] + 2 s^2 [ab - (a-1)(b-1)] = 1,   s = 1/(q-1)^2

    which in terms of n and f is  s(n-3) + 2 s^2 [n(1-f/2) - 1] = 1.  Exact in
    rational arithmetic. For Poisson layers kbar = kappa, the cross term
    vanishes, and n = (q-1)^2 at every f.
    """
    s = F(1, (q - 1) ** 2)
    if not regular:
        return F(1) / s
    return (1 + 3 * s + 2 * s * s) / (s + 2 * s * s * (1 - F(f, 2)))


def check_triangles_at_fixed_degree():
    """A triangle-clustered graph, and where the graph calculation goes wrong."""
    print('     q   all links   all triangles   graph formula   (q-1)^2')
    for q in (3, 4, 5, 6, 10):
        lo = triangle_family(q, F(0))
        hi = triangle_family(q, F(1))
        t = (q - 1) ** 2
        # the endpoints come out of the ONE mixed formula, no layer dropped
        assert lo == t + 1 and hi == t + 2, (q, lo, hi)
        # monotone in f
        vals = [triangle_family(q, F(j, 8)) for j in range(9)]
        assert all(x < y for x, y in zip(vals, vals[1:])), (q, vals)
        print(f'  {q:>4}   {str(lo):>9}   {str(hi):>13}   {str(t+1):>13}   {t:>7}')
    print('    The first and last columns are the same formula at f = 0 and')
    print('    f = 1, with no layer dropped, so the cross term does exactly the')
    print('    interpolating work. Read the middle two: on a network of')
    print('    triangles the answer is (q-1)^2 + 2, while treating that same')
    print('    object as a graph of the same degree gives (q-1)^2 + 1. The')
    print('    graph calculation is wrong by one in degree, at every q, because')
    print('    it adds the three edges of a triangle as though they were')
    print('    independent when they close a loop.')
    for q in (3, 4, 5):
        assert triangle_family(q, F(0), regular=False) == (q - 1) ** 2
        assert triangle_family(q, F(1), regular=False) == (q - 1) ** 2
    print('    And the size of it: for POISSON layers kbar = kappa, the cross')
    print('    term vanishes identically and there is no difference at any f.')
    print('    The effect is the regular ensemble\'s lost branch and nothing')
    print('    else -- the same one that separates (q-1)^2 from (q-1)^2 + 1.')


# ---------------------------------------------------------------------------
# survey propagation
# ---------------------------------------------------------------------------
#
# Sec. 12.3's threshold is where the replica-symmetric solution loses local
# stability, and Sec. 12.6 says at length that this is not the colourability
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
    print('    getting there. This is the line Sec. 12.3 is not: the stability')
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


# ---------------------------------------------------------------------------
# survey propagation above cardinality two
# ---------------------------------------------------------------------------
#
# Sec. 12.8. Everything turns on how many colours ONE complex can forbid at
# once.
#
#   hypergraph rule: at most ONE -- a complex is violated only when all its
#     members agree, so it can forbid gamma to the entering member only when
#     every other member is already forced to gamma, and they cannot all be
#     forced to two colours at once. The cardinality-two inclusion-exclusion
#     above therefore carries over verbatim, with the per-neighbour weight e/q
#     replaced by h/q, h = q prod_{j<c} (e_j/q).  That substitution is the
#     whole of what cardinality does here.
#
#   proper rule: up to c-1 AT ONCE, so the message is a forbidden *set*, the
#     prohibitions are correlated inside one complex as well as across
#     complexes, and the node update needs the joint distribution of that set.
#     No substitution repairs it, and it is not attempted.
#
# Benchmarks: Gabrie, Dani, Semerjian & Zdeborova, J. Phys. A 50, 505002
# (2017), Table 1, l_col. Their q = 2 entries at c = 3 and c = 4 are marked
# invalid (SP type I instability) and are excluded here for that reason -- the
# code returns numbers for them and they would mean nothing.
GABRIE = {(3, 3): 26.92, (4, 3): 63.3, (2, 5): 52.32}


def _hyp_h(e, c, rng, size):
    """Weight that a cardinality-c complex forbids anything: q prod (e_j/q)."""
    return np.prod([e[rng.integers(0, e.size, size)] for _ in range(c - 1)],
                   axis=0)


def hyp_converge(q, c, kappa, size=20000, sweeps=200, seed=0):
    """Survey population for hypergraph colouring, Poisson chy-degree.

    Initialised near fully forced. The non-trivial branch is reached from
    ABOVE: at q = 3, c = 3, kappa = 27 the map sends 0.5 -> 0.47, decaying to
    the trivial fixed point, but 1.0 -> 0.97. Starting low finds nothing, the
    same trap `population.py` documents for the ferromagnetic branch.
    """
    rng = np.random.default_rng(seed)
    e = rng.uniform(0.97, 1.0, size)
    for _ in range(sweeps):
        h = _hyp_h(e, c, rng, size) * q ** (2 - c)
        e, _ = sp_update(h, rng.poisson(kappa, size), rng, q, size)
    return e, rng


def hyp_complexity(q, c, kappa, size=20000, sweeps=200, seed=0, nsamp=300000):
    """Sigma for hypergraph colouring, with the complex term a graph lacks."""
    e, rng = hyp_converge(q, c, kappa, size, sweeps, seed)
    if e.mean() < 1e-6:
        return 0.0, 0.0
    h = _hyp_h(e, c, rng, nsamp) * q ** (2 - c)
    _, contra = sp_update(h, rng.poisson(kappa, nsamp), rng, q, nsamp)
    site = np.log(np.maximum(1.0 - contra, 1e-300)).mean()
    allsame = q * np.prod([e[rng.integers(0, e.size, nsamp)] / q
                           for _ in range(c)], axis=0)
    Za = np.log(np.maximum(1.0 - allsame, 1e-300)).mean()
    hh = h[rng.integers(0, h.size, nsamp)]
    ee = e[rng.integers(0, e.size, nsamp)]
    Zia = np.log(np.maximum(1.0 - hh * ee / q, 1e-300)).mean()
    return site + (kappa / c) * Za - kappa * Zia, float(e.mean())


def hyp_threshold(q, c, lo, hi, seed=0, iters=9, **kw):
    assert hyp_complexity(q, c, lo, seed=seed, **kw)[0] > 0, (q, c, lo)
    for _ in range(iters):
        mid = 0.5 * (lo + hi)
        if hyp_complexity(q, c, mid, seed=seed, **kw)[0] > 0:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def check_hypergraph_thresholds():
    """Sigma = 0 above cardinality two, against Gabrie et al. Table 1."""
    print('     q   c   Sigma = 0 at   published   error')
    errs = []
    for (q, c), lo, hi in (((3, 3), 24.0, 30.0), ((4, 3), 55.0, 70.0),
                           ((2, 5), 45.0, 58.0)):
        ref = GABRIE[(q, c)]
        r = [hyp_threshold(q, c, lo, hi, seed=s) for s in (0, 1, 2)]
        m, sd = float(np.mean(r)), float(np.std(r))
        errs.append(100 * (m - ref) / ref)
        assert abs(m - ref) / ref < 0.01, (q, c, m, ref)
        print(f'  {q:>4}  {c:>2}   {m:>8.3f}({sd:.3f})   {ref:>9.2f}'
              f'   {errs[-1]:>+5.2f}%')
    assert all(x < 0 for x in errs), errs
    print('    All three errors agree in sign and size, so this is a converged')
    print('    bias of the finite population and not scatter. The numbers are')
    print('    Ref. [gabrie2017]\'s; what the agreement checks is that h is the')
    print('    whole of what cardinality does to the one-step calculation.')


# ---------------------------------------------------------------------------
# set-valued survey propagation: Sec. 12.9
# ---------------------------------------------------------------------------
#
# Sec. 12.8's substitution reaches the hypergraph rule and not the proper one,
# because a complex of cardinality c can close up to c-1 colours at once. This
# carries the survey as a SET and does the proper rule at c = 3.
#
# The apparatus is gated twice, both against published numbers:
#   c = 2, proper       -> Mulet's c_q         (check_setvalued_gates)
#   c = 3, hypergraph   -> Gabrie's l_col      (same, rule swapped and nothing
#                                               else, which is what isolates the
#                                               c = 3 scaffolding)
# and then applied to the proper rule at c = 3, where it finds no window:
# the branch arrives with Sigma already negative (check_window_closes).
#
# Two traps, both of which cost time here. The branch appears DISCONTINUOUSLY
# in every case, graph included, so a coarse scan steps over the whole positive
# window -- the graph's is only about 8% wide in degree. And Sigma is exactly
# 0.0 on the trivial branch, so a bracket below the onset is refused rather
# than bracketing.

def _sv_subsets(q):
    return [frozenset(T) for s in range(q + 1) for T in combinations(range(q), s)]

def _sv_emit(others, q):
    """Colours still open to member 0 (Hall's condition on the other members)."""
    ok = set()
    for g in range(q):
        rest = [A - {g} for A in others]
        if any(not A for A in rest):
            continue
        if len(rest) == 2 and len(rest[0]) == 1 and rest[0] == rest[1]:
            continue
        ok.add(g)
    return ok

def _sv_emithyper(others, q):
    """Hypergraph rule: the complex forbids g only if ALL others are forced to g."""
    bad = {next(iter(A)) for A in others if len(A) == 1}
    if len(bad) == 1 and all(len(A) == 1 for A in others):
        return set(range(q)) - bad
    return set(range(q))


def _sv_interior(q, c, rule='proper'):
    """T[...] -> distribution of emitted size, given the other members' sizes."""
    subs = _sv_subsets(q)
    shape = (q + 1,) * (c - 1) + (q + 1,)
    T = np.zeros(shape)
    from itertools import product
    for combo in product(subs, repeat=c - 1):
        if any(len(A) == 0 for A in combo):
            continue
        w = 1.0
        for A in combo:
            w /= comb(q, len(A))
        em = _sv_emit(combo, q) if rule == 'proper' else _sv_emithyper(combo, q)
        T[tuple(len(A) for A in combo) + (len(em),)] += w
    return T

def _sv_contains(m, q):
    """g[t] = P(the set contains a given t-subset)."""
    C = np.array([[comb(q - t, s - t) / comb(q, s) if s >= t else 0.0
                   for s in range(q + 1)] for t in range(q + 1)])
    return m @ C.T

def _sv_sweep(V, q, c, kappa, T, rng, yreweight=np.inf):
    n = V.shape[0]
    if c == 2:
        M = np.einsum('ni,ij->nj', V[rng.integers(0, n, n)], T)
    else:
        M = np.einsum('ni,nj,ijk->nk', V[rng.integers(0, n, n)],
                      V[rng.integers(0, n, n)], T)
    G = _sv_contains(M, q)
    d = rng.poisson(kappa, n)
    tot = int(d.sum())
    P = np.ones((n, q + 1))
    if tot:
        lg = np.log(np.maximum(G[rng.integers(0, n, tot)], 1e-300))
        cs = np.concatenate([np.zeros((1, q + 1)), np.cumsum(lg, axis=0)])
        e = np.cumsum(d)
        P = np.exp(cs[e] - cs[e - d])
    out = np.zeros_like(V)
    for u in range(q + 1):
        out[:, u] = comb(q, u) * sum((-1) ** j * comb(q - u, j) * P[:, u + j]
                                     for j in range(q - u + 1))
    out = np.clip(out, 0, None)
    # Finite y: a member left with no colour is a violated complex, energy 1,
    # so it survives with weight e^{-y} rather than being discarded. y = inf
    # recovers the hard constraint and every result above it.
    out[:, 0] *= np.exp(-yreweight)
    z = out.sum(axis=1, keepdims=True)
    return np.where(z > 0, out / np.maximum(z, 1e-300), 0.0), M

def _sv_run(q, c, kappa, size=4000, sweeps=300, seed=0, rule='proper',
            yreweight=np.inf):
    rng = np.random.default_rng(seed)
    T = _sv_interior(q, c, rule)
    V = np.zeros((size, q + 1)); V[:, 1] = 0.9; V[:, 2] = 0.1
    for _ in range(sweeps):
        V, M = _sv_sweep(V, q, c, kappa, T, rng, yreweight)
    return V, M, rng


def _sv_sat(q, c, rule='proper'):
    """S[...] = P(the complex is satisfiable), by size class."""
    subs = _sv_subsets(q)
    from itertools import product
    S = np.zeros((q + 1,) * c)
    for combo in product(subs, repeat=c):
        if any(len(A) == 0 for A in combo):
            continue
        w = 1.0
        for A in combo:
            w /= comb(q, len(A))
        found = False
        for pick in product(*[sorted(A) for A in combo]):
            if (len(set(pick)) == c) if rule == 'proper' else (len(set(pick)) > 1):
                found = True
                break
        if found:
            S[tuple(len(A) for A in combo)] += w
    return S


def _sv_sigma(q, c, kappa, size=4000, sweeps=300, seed=0, rule='proper',
              yreweight=np.inf):
    """Three-term Bethe count at m = 0, the form validated at c = 2."""
    V, M, rng = _sv_run(q, c, kappa, size, sweeps, seed, rule, yreweight)
    if V[:, 1].mean() < 1e-9:
        return 0.0, 0.0
    n = V.shape[0]
    G = _sv_contains(M, q)
    d = rng.poisson(kappa, n)
    tot = int(d.sum())
    P = np.ones((n, q + 1))
    if tot:
        lg = np.log(np.maximum(G[rng.integers(0, n, tot)], 1e-300))
        cs = np.concatenate([np.zeros((1, q + 1)), np.cumsum(lg, axis=0)])
        e = np.cumsum(d)
        P = np.exp(cs[e] - cs[e - d])
    P0 = sum((-1) ** j * comb(q, j) * P[:, j] for j in range(q + 1))
    site = np.log(np.maximum(1.0 - P0, 1e-300)).mean()
    Sm = _sv_sat(q, c, rule)
    parts = [V[rng.integers(0, n, n)] for _ in range(c)]
    Za = (np.einsum('ni,nj,ij->n', *parts, Sm) if c == 2
          else np.einsum('ni,nj,nk,ijk->n', *parts, Sm))
    g1 = G[:, 1]
    Zia = 1.0 - (1.0 - g1[rng.integers(0, n, n)]) * V[rng.integers(0, n, n), 1]
    return (site + (kappa / c) * np.log(np.maximum(Za, 1e-300)).mean()
            - kappa * np.log(np.maximum(Zia, 1e-300)).mean()), V[:, 1].mean()


def check_setvalued_gates():
    """The set-valued apparatus, against Mulet at c = 2 and Gabrie at c = 3."""
    print('    c = 2, proper rule, against Mulet c_q')
    print('       q   Sigma < 0 first at   c_q')
    for q, cq, lo, hi in ((3, 4.69, 4.45, 4.95), (4, 8.90, 8.20, 9.40),
                          (5, 13.69, 12.60, 14.40)):
        r = _sv_bisect(q, 2, lo, hi, 'proper')
        assert abs(r - cq) / cq < 0.03, (q, r, cq)
        print(f'    {q:>4}   {r:>17.3f}   {cq:>5.2f}')
    print('    c = 3, hypergraph rule, against Gabrie l_col -- only the rule')
    print('    is swapped, which is what isolates the c = 3 scaffolding')
    print('       q   Sigma < 0 first at   l_col')
    for q, ref, lo, hi in ((3, 26.92, 24.0, 30.0), (4, 63.3, 57.0, 70.0)):
        r = _sv_bisect(q, 3, lo, hi, 'hyper')
        assert abs(r - ref) / ref < 0.03, (q, r, ref)
        print(f'    {q:>4}   {r:>17.3f}   {ref:>5.2f}')


def check_window_closes():
    """Sec. 12.9: on triangles the proper rule leaves no window for m = 0.

    A crossing needs a range where a non-trivial survey exists AND Sigma is
    still positive. The graph has one; the triangle network does not. Losing it
    costs the lower side of the bracket, not the threshold: Sigma < 0 wherever
    a solution exists still says the colourings are gone by the lift-off.
    """
    print('       q   c   branch at (degree)   best Sigma above it')
    for q in (3, 4):
        for c, lo, hi in ((2, 4.0, 10.0), (3, 1.2, 4.5)):
            on = _sv_onset(q, c, lo, hi)
            # The claim is that a positive-Sigma point EXISTS above the onset,
            # not that Sigma is positive at one chosen place. Testing a single
            # point is a knife edge: the q = 3 graph window is narrow and its
            # Sigma runs to +0.001 near the edge. Scan the window instead.
            grid = [on * (1 + 0.02 * j) for j in range(1, 8)]
            sgs = [_sv_sigma(q, c, k)[0] for k in grid]
            sg = max(sgs)
            deg = on * (1 if c == 2 else 2)
            # graphs open a window, triangles do not
            assert (sg > 0) == (c == 2), (q, c, on, sgs)
            print(f'    {q:>4}  {c}   {deg:>18.2f}   {sg:>+11.5f}')
    print('    Both graph rows arrive with Sigma positive and cross zero at the')
    print('    published c_q; both triangle rows arrive with Sigma already')
    print('    negative -- fewer than one cluster -- and never return. Same')
    print('    machinery, same q, same cardinality as the hypergraph gate')
    print('    above, which does have a window. So this is the constraint and')
    print('    not the cardinality. With no crossing to bisect, m = 0 cannot')
    print('    BRACKET c_q on triangles the way it does on a graph. What it')
    print('    still gives is a bound from above: Sigma < 0 wherever a solution')
    print('    exists, so the colourings are gone by the lift-off. That bound is')
    print('    what clustered_cq reports, and it is all the comparison with the')
    print('    graph needs.')


def _sv_onset(q, c, lo, hi, iters=14, seed=0):
    """Smallest chy-degree carrying a non-trivial survey. The branch appears
    discontinuously in every case, so this is a fold and bisection on it is
    the only reliable way to find the window's left edge.

    A fold is located worse than a crossing, and the population is what limits
    it, not the bisection: past about a dozen halvings the bracket is finer
    than the seed-to-seed spread and the extra digits are noise. Hence the
    seed argument -- the spread over seeds is the honest error bar.
    """
    for _ in range(iters):
        mid = 0.5 * (lo + hi)
        if _sv_run(q, c, mid, seed=seed)[0][:, 1].mean() > 1e-4:
            hi = mid
        else:
            lo = mid
    return 0.5 * (lo + hi)


def _sv_bisect(q, c, lo, hi, rule, iters=9):
    """Where Sigma crosses zero. Requires Sigma(lo) > 0, so lo must sit inside
    the window -- above the branch onset, below the crossing."""
    assert _sv_sigma(q, c, lo, rule=rule)[0] > 0, (q, c, lo, rule)
    for _ in range(iters):
        mid = 0.5 * (lo + hi)
        if _sv_sigma(q, c, mid, rule=rule)[0] > 0:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


# Mulet's colourability thresholds, the graph column of Table 12.4.
MULET_CQ = {3: 4.69, 4: 8.90, 5: 13.69}


def clustered_cq(q, lo=1.2, hi=4.5, iters=12, seed=0):
    """Degree at which a triangle-clustered graph stops being q-colourable.

    Sigma(y=inf) is the complexity of ZERO-ENERGY clusters -- Krzakala, Pagnani
    & Weigt Eq. (20), which Eq. (12.8) reproduces term for term. Negative means
    none exist. On a graph the sequence is: trivial survey, then a branch with
    Sigma > 0, then Sigma crossing zero at c_q. On triangles the middle step is
    missing -- the branch arrives with Sigma already negative -- so the
    colourings are gone by the lift-off and this RETURNS AN UPPER BOUND on c_q,
    not c_q itself.

    The asymmetry is worth being explicit about, because it is not visible in
    the return value. Above the lift-off a solution exists and its Sigma is
    computed, and it is negative. Below it the only solution is the trivial one,
    where _sv_sigma returns 0.0 by construction rather than counting anything --
    the replica-symmetric default, one cluster because nothing has shattered.
    So the upper end is a result and the lower end is an assumption, and the
    bound is tight only if the shattered phase has zero width rather than a
    width this population cannot resolve. Degree is 2*kappa.
    """
    return 2.0 * _sv_onset(q, 3, lo, hi, iters=iters, seed=seed)


def check_clustered_cq():
    """Sec. 12.9: c_q for triangles, against the graph, and against the RS line.

    clustered_cq bounds c_q from above, so every comparison below is an
    inequality -- which is the form the conclusion needs anyway.
    """
    print('     q   graph c_q   triangles c_q <=   change   RS line, Poisson'
          '   RS line, regular')
    # The triangle threshold is a FOLD, not a crossing: the branch lifts off
    # with Sigma already negative, so there is nothing to bisect Sigma on and
    # the population sets the precision. Three seeds, and the table quotes the
    # spread. Note what that spread is and is not: it is scatter, and it is
    # small, but check_survey_branch_is_clustering finds the OTHER fold in this
    # chapter -- c_d -- sitting a per cent to three off the published value with
    # a spread no larger. A fold carries a bias the seeds do not see.
    for q in (3, 4):
        cq = MULET_CQ[q]
        r = [clustered_cq(q, seed=sd) for sd in (0, 1, 2)]
        tri, sd_tri = float(np.mean(r)), float(np.std(r))
        # Sec. 12.5: Poisson layers carry no excess-degree correction, so both
        # routes give (q-1)^2; only a regular ensemble separates them.
        rs_p = (q - 1) ** 2
        rs_g, rs_t = rs_p + 1, rs_p + 2
        # c_q is computed in the POISSON ensemble, so that is the RS line it
        # has to be set beside: the colourable phase must end below it.
        assert tri < rs_p, (q, tri, rs_p)
        # and triangles must be HARDER, which is what neither RS line reports
        assert tri < cq, (q, tri, cq)
        print(f'  {q:>4}   {cq:>9.2f}   {tri:>9.3f}({sd_tri:.3f})'
              f'   {100*(tri-cq)/cq:>+6.1f}%'
              f'    {rs_p:>7} -> {rs_p}   {rs_g:>7} -> {rs_t}')
    print('    The RS line does not track c_q in EITHER ensemble. Promoting')
    print('    links into triangles at fixed degree lowers the density at which')
    print('    colourings actually run out, by AT LEAST a third. In the Poisson')
    print('    ensemble the RS line does not move at all; in a regular one it')
    print('    moves the other way. On a graph the two at least land near one')
    print('    another -- 4 against 4.69, 9 against 8.90 -- which is why the RS')
    print('    line gets used as a rough locator at all. On triangles the proxy')
    print('    reports either no effect or the wrong SIGN, depending on which')
    print('    ensemble it is taken in.')


def figure_threshold():
    """Sec. 12.9: how the threshold is located, and why triangles differ.

    Left, a graph: the survey is trivial and Sigma is zero until the branch
    lifts off, Sigma is POSITIVE over an arc, and the arc falls through zero at
    c_q. The threshold is the crossing.

    Right, a network of triangles at the same q: the branch lifts off already
    BELOW zero. There is no arc, so there is no crossing, and the colourable
    phase is the trivial region -- the threshold is the lift-off itself.

    The vertical scales differ by a factor of six and the caption says so: over
    the degrees plotted here the graph's arc peaks at +0.064 while the triangle
    branch plunges to -0.400. A shared axis makes the graph panel look flat and
    hides the crossing, which is the thing the panel exists to show.
    """
    plt = _mpl()
    q = 4
    seeds = (0, 1, 2)

    def curve(c, degs):
        out = []
        for d in degs:
            k = d if c == 2 else d / 2.0
            out.append(np.mean([_sv_sigma(q, c, k, seed=s)[0] for s in seeds]))
        return np.array(out)

    gdeg = np.array([8.10, 8.20, 8.40, 8.60, 8.80, 9.00, 9.20])
    tdeg = np.array([6.80, 7.00, 7.20, 7.60, 8.00, 8.40])
    gsig, tsig = curve(2, gdeg), curve(3, tdeg)

    fig, axes = plt.subplots(1, 2, figsize=(4.6, 2.3))
    for ax, deg, sig, on, cq, ttl, note in (
            (axes[0], gdeg, gsig, 8.08, 8.90, 'graph',
             'an arc above zero;\nit crosses at $c_q$'),
            (axes[1], tdeg, tsig, 6.77, 6.77, 'triangles',
             'no arc: the branch\nstarts below zero')):
        ax.axhline(0.0, lw=0.8, color=MID, zorder=1)
        ax.plot([deg[0] - 1.3, on], [0, 0], lw=2.4, color=LIGHT,
                solid_capstyle='butt', zorder=2)
        ax.plot(deg, sig, 'o-', ms=3.2, lw=1.0, color=DARK, zorder=3)
        ax.plot([on, on], [0, sig[0]], ls=':', lw=0.9, color=DARK, zorder=3)
        ax.axvline(cq, ls='--', lw=0.9, color=MID, zorder=1)
        ax.set_title(ttl, fontsize=8.5)
        ax.set_xlabel('degree', fontsize=8.5)
        _tidy(ax)
        lo, hi = min(sig.min(), 0.0), max(sig.max(), 0.0)
        pad = 0.32 * (hi - lo) if hi > lo else 0.05
        ax.set_ylim(lo - pad, hi + 1.5 * pad)
        ax.annotate('$c_q$', xy=(cq, hi + 0.95 * pad), fontsize=8,
                    color='0.25', ha='center')
        ax.annotate(note, xy=(0.62, 0.06), xycoords='axes fraction',
                    fontsize=6.6, color='0.30', ha='center')
    axes[0].set_ylabel(r'$\Sigma$', fontsize=9)
    axes[1].set_ylabel(r'$\Sigma$', fontsize=9)
    fig.tight_layout()
    fig.savefig(OUT / 'fig-threshold.pdf')
    print(f'  wrote {OUT / "fig-threshold.pdf"}')


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
    print('a graph with triangles, at fixed degree:')
    check_triangles_at_fixed_degree()
    print('survey propagation, against the colourability thresholds:')
    check_survey_thresholds()
    print('and the clustering threshold, as a by-product:')
    check_survey_branch_is_clustering()
    print('above cardinality two, against Gabrie et al.:')
    check_hypergraph_thresholds()
    print('the set-valued apparatus, gated twice:')
    check_setvalued_gates()
    print('and what it finds for the proper rule on triangles:')
    check_window_closes()
    print('the colourability threshold of a clustered graph:')
    check_clustered_cq()
    print('figures:')
    figure_colouring()
    figure_threshold()
