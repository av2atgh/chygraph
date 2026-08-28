"""Chapter 8: the general theory, and the checks the chapter quotes.

  fig-exchange   the two Perron roots leaving 1 in opposite directions at the
                 same rate: the trivial fixed point loses stability exactly
                 where the non-trivial one gains it.

Everything else here is a check.  Chapter 8's other figure, the two steps, is
TikZ in `statmech.tex`.
"""

import sys
from pathlib import Path

import numpy as np
from scipy.optimize import brentq

sys.path.insert(0, str(Path.home() / 'av2atg' / 'chygraph_statmech' / 'src'))
sys.path.insert(0, str(Path.home() / 'av2atg' / 'chygraph' / 'src'))
from chygraph import hypergraph_giant  # noqa: E402
from chygraph_statmech import Chygraph  # noqa: E402
from chygraph_statmech.freeenergy import paramagnetic  # noqa: E402
from chygraph_statmech.fixedpoint import FixedPointStability  # noqa: E402

OUT = Path(__file__).resolve().parent
DARK, MID, LIGHT = '0.10', '0.45', '0.70'


def _mpl():
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    return plt


def _tidy(ax):
    ax.tick_params(labelsize=8)
    for sp in ('top', 'right'):
        ax.spines[sp].set_visible(False)


# ------------------------------------------------- (1) u' is Chapter 7's tau
def check_uprime_is_tau():
    """Chapter 7's transmission at q = 2 is this chapter's u'.

    Chapter 7 defined tau_c(q, v) as the derivative of the Potts emitted ratio
    at the symmetric fixed point; Sec. 8.3 defines u' as the derivative of the
    emitted cavity field.  For the Ising interaction they are the same number,
    which is what makes Eq. (7.9) the single-layer case of det(I - B) = 0.
    """
    from sympy import Rational, cancel, diff, simplify, symbols
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    from potts import clique, emitted_ratio
    from chygraph_statmech.ising import clique_derivative

    q, v, t = symbols('q v t')
    for c in (2, 3, 4):
        members, edges = clique(c)
        ys = {j: symbols(f'y{j}') for j in members}
        rho = emitted_ratio(members, edges, ys, q, v)
        at1 = {ys[j]: 1 for j in members if j != 0}
        tau = cancel(diff(rho, ys[1]).subs(at1)).subs(q, 2)
        # Ising coupling J_I = J_P / 2, so v = exp(2 beta J_I) - 1
        tau_t = simplify(tau.subs(v, 2 * t / (1 - t)))
        for bJ in (0.15, 0.4, 0.9):
            want = clique_derivative(c, bJ)
            got = float(tau_t.subs(t, np.tanh(bJ)))
            assert abs(got - want) < 1e-12, (c, bJ, got, want)
        print(f'  c = {c}: tau_c(2, v) = u\'_c   OK   [{tau_t}]')


# ------------------------------------------------ (2) the size-biased average
def check_size_biased(kappa=4.0, sizes=(3, 5)):
    """A mixed-cardinality layer, three ways.

    The branching matrix needs <sbar u'>, the size-biased average of
    (c-1)u'(c).  Splitting the layer into one layer per cardinality with
    chy-degrees in the ratio c p_c must give the same transition; replacing the
    layer by its mean cardinality must not, and does not.
    """
    a, b = sizes
    mixed = Chygraph([{a: 0.5, b: 0.5}], [kappa])
    cbar = 0.5 * (a + b)
    split = Chygraph([a, b], [kappa * a * 0.5 / cbar, kappa * b * 0.5 / cbar])
    naive = Chygraph([int(round(cbar))], [kappa])

    sbar = mixed.excess_cardinality_layer(0)
    print(f'  layer mixing c = {a} and {b} equally, <kappa> = {kappa:g}')
    print(f'    <sbar> = <c^2>/<c> - 1 = {sbar:g},  not <c> - 1 = {cbar - 1:g}')
    Tm, Ts, Tn = (m.critical_temperature()
                  for m in (mixed, split, naive))
    assert abs(Tm - Ts) < 1e-9, (Tm, Ts)
    print(f'    T_c mixed layer      {Tm:.6f}')
    print(f'    T_c split into two   {Ts:.6f}   (identical, as it must be)')
    print(f'    T_c at mean c = {int(round(cbar))}      {Tn:.6f}   '
          f'(wrong by {100 * (Tm - Tn) / Tm:.1f}%)')
    # and <sbar u'> is not <sbar><u'>
    bJ = mixed.critical_coupling()
    prod = mixed.transmission(0, bJ)
    fact = mixed.excess_cardinality_layer(0) * mixed.u_prime(0, bJ)
    print(f'    <sbar u\'> = {prod:.6f}, <sbar><u\'> = {fact:.6f}   '
          f'(differ by {100 * abs(prod - fact) / prod:.2f}%)')
    assert abs(prod - fact) > 1e-6


# ------------------------------------------- (3) links and triangles, Eq. (15)
def check_links_triangles(kL=4.0, kT=2.0):
    """The determinant, and the cross term that vanishes only for Poisson."""
    for regular in (False, True):
        g = Chygraph([2, 3], [kL, kT], regular=regular)
        bJ = g.critical_coupling()
        B = g.branching_matrix(bJ)
        uL, uT = g.u_prime(0, bJ), g.u_prime(1, bJ)
        cross = (g.k[0] * g.k[1] - g.kbar[0] * g.kbar[1])
        det = np.linalg.det(np.eye(2) - B)
        # Eq. (8.9) written out
        closed = (1 - g.kbar[0] * uL - 2 * g.kbar[1] * uT
                  - 2 * uL * uT * cross)
        assert abs(det) < 1e-10 and abs(closed) < 1e-9, (det, closed)
        name = 'regular' if regular else 'Poisson'
        print(f'  {name:<8} kbar = {list(g.kbar)}:  T_c = '
              f'{1 / bJ:.6f},  cross term coefficient = {cross:g}')


# ------------------------------------------------------ (4) the free energy
def check_free_energy(kappa=4.0):
    """Population dynamics against the paramagnetic closed form."""
    g = Chygraph([2], [kappa])
    for bJ in (0.10, 0.20):
        num = g.free_energy(bJ, sweeps=200).minus_beta_f()
        closed = paramagnetic([2], [kappa], bJ)
        book = np.log(2) + (kappa / 2) * np.log(np.cosh(bJ))
        assert abs(num - closed) < 1e-9 and abs(closed - book) < 1e-12
        print(f'  beta J = {bJ:.2f}:  -beta f = {num:.10f} (population) '
              f'= {closed:.10f} (closed form) = {book:.10f} (textbook)')


# --------------------------------------------- (5) the exchange of stability
MODELS = (('graph', 'Poisson graph <k> = 3',
           lambda: hypergraph_giant(graph=True), {'q': 1.0, 'k': 3.0}, DARK),
          ('hypergraph', 'Poisson hypergraph <k> = 2, <c> = 3',
           hypergraph_giant, {'q': 1.0, 'k': 2.0, 'c': 3.0}, MID))


def _roots(F, sub0, p):
    lam = F.trivial_perron_root({**sub0, 'p': p})
    rho = F.spectral_radius({**sub0, 'p': p})
    return lam - 1.0, rho


def check_exchange():
    """rho(J(Q*)) = 1 - Lambda + O(Lambda^2), and the order is quadratic."""
    for _, plain, make, sub0, _ in MODELS:
        F = FixedPointStability(make())
        pc = brentq(lambda x: F.trivial_perron_root({**sub0, 'p': x}) - 1,
                    1e-3, 1 - 1e-9)
        print(f'  {plain}:  p_c = {pc:.10f}')
        rows = []
        for d in (3e-2, 1e-2, 3e-3, 1e-3):
            L, rho = _roots(F, sub0, pc * (1 + d))
            rows.append((L, rho, 1 - rho - L))
            print(f'    Lambda = {L:.3e}   rho = {rho:.12f}   '
                  f'residual = {1 - rho - L:+.3e}   /Lambda^2 = '
                  f'{(1 - rho - L) / L**2:+.4f}')
        # the residual must fall faster than Lambda: quadratic, not linear
        for (L0, _, r0), (L1, _, r1) in zip(rows, rows[1:]):
            assert abs(r1 / r0) < 0.6 * (L1 / L0), (r0, r1)
        print('    residual falls faster than Lambda at every step   OK')
        print('    (closer in, the residual reaches 1e-8 and rounding takes '
              'over; the sequence is stopped before it does)')


def figure_exchange():
    """Left: the two roots crossing.  Right: where the models actually differ."""
    plt = _mpl()
    fig, axes = plt.subplots(1, 2, figsize=(4.6, 2.4))

    # (a) the exchange itself, on one model
    name, _, make, sub0, _ = MODELS[0]
    F = FixedPointStability(make())
    pc = brentq(lambda x: F.trivial_perron_root({**sub0, 'p': x}) - 1,
                1e-3, 1 - 1e-9)
    ds = np.linspace(-0.05, 0.05, 61)
    lam = [F.trivial_perron_root({**sub0, 'p': pc * (1 + d)}) for d in ds]
    up = ds > 0
    rho = [F.spectral_radius({**sub0, 'p': pc * (1 + d)}) for d in ds[up]]
    ax = axes[0]
    ax.plot(ds, lam, '-', lw=1.4, color=DARK)
    ax.plot(ds[up], rho, '--', lw=1.4, color=DARK, dashes=(4, 2))
    ax.axhline(1.0, ls=':', lw=0.8, color='0.6')
    ax.axvline(0.0, ls=':', lw=0.8, color='0.6')
    ax.annotate(r'$1+\Lambda$', xy=(0.043, lam[-3]), fontsize=7.5,
                color='0.25', ha='right', va='bottom')
    ax.annotate(r'$\rho(J(Q^{*}))$', xy=(0.043, rho[-3]), fontsize=7.5,
                color='0.25', ha='right', va='top')
    ax.set_xlabel(r'$(p-p_{c})/p_{c}$', fontsize=8.5)
    ax.set_ylabel('Perron root', fontsize=8.5)
    ax.set_ylim(0.972, 1.028)
    _tidy(ax)

    # (b) the residual, scaled by Lambda^2: this is where the models differ
    ax = axes[1]
    for label, _, make, sub0, col in MODELS:
        F = FixedPointStability(make())
        pc = brentq(lambda x: F.trivial_perron_root({**sub0, 'p': x}) - 1,
                    1e-3, 1 - 1e-9)
        L, R = [], []
        for d in np.logspace(np.log10(3e-3), np.log10(1e-1), 22):
            lam_, rho_ = _roots(F, sub0, pc * (1 + d))
            L.append(lam_)
            R.append((1 - rho_ - lam_) / lam_**2)
        ax.semilogx(L, R, 'o-', ms=2.6, lw=1.0, color=col, label=label)
    ax.set_xlabel(r'$\Lambda$', fontsize=8.5)
    ax.set_ylabel(r'$(1-\rho-\Lambda)/\Lambda^{2}$', fontsize=8.5)
    ax.legend(frameon=False, fontsize=7.5, loc='center')
    _tidy(ax)

    fig.tight_layout()
    fig.savefig(OUT / 'fig-exchange.pdf')
    print(f'  wrote {OUT / "fig-exchange.pdf"}')


if __name__ == '__main__':
    print("Chapter 7's tau at q = 2 against Chapter 8's u':")
    check_uprime_is_tau()
    print('the size-biased average:')
    check_size_biased()
    print('links and triangles:')
    check_links_triangles()
    print('the free energy:')
    check_free_energy()
    print('the exchange of stability:')
    check_exchange()
    print('figure:')
    figure_exchange()
