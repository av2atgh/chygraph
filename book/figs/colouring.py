"""Chapter 12: graph colouring on a chygraph.

Unlike every other chapter, this one has no manuscript behind it, so the
calculations are here rather than in `../src`. Everything is done twice: once
by exact enumeration of the interior of a complex, once from a closed form
derived by hand, and the two are required to agree.

  fig-colouring   the two constraints. Proper colouring -- every member of a
                  complex a different colour -- against hypergraph colouring,
                  which only forbids a monochromatic complex. Identical at
                  cardinality two and violently different above it.

Benchmarks, from `~/Downloads/chygraph_references/`:
  Mulet, Pagnani, Weigt & Zecchina, Phys. Rev. Lett. 89, 268701 (2002), Table I
  Zdeborova & Krzakala, Phys. Rev. E 76, 031131 (2007), Eq. (18)
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
    print('figure:')
    figure_colouring()
