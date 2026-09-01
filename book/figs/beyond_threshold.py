"""Chapter 5: what the threshold tensor does not know.

Three figures, all computed with the `chygraph` package so that the book quotes
the same numbers the software returns.

  fig-degeneracy   two degree distributions sharing every entry of the
                   threshold tensor, and therefore the threshold, but differing
                   in the order parameter everywhere above it and in the
                   critical amplitude;
  fig-andor        AND-logic against OR-logic hyperedge percolation -- one
                   chygraph, one generating function apart;
  fig-correlated   three joint (link, triangle) degree distributions with
                   identical marginals and three different transitions,
                   with Monte Carlo.
"""

import sys
from pathlib import Path

import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
from sympy import Rational, nsolve, symbols

sys.path.insert(0, str(Path.home() / 'av2atg' / 'chygraph' / 'src'))
from chygraph import (CriticalAmplitude, JointGiantComponent,  # noqa: E402
                      and_or_hypergraph, correlated_cardinality_hypergraph,
                      finite_pgf, hypergraph_giant, two_class_joint_degree)

OUT = Path(__file__).resolve().parent
RNG = np.random.default_rng(5)
DARK, MID, LIGHT = '0.10', '0.45', '0.70'


def _mpl():
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    return plt


def excess_pgf(d):
    """PGF of the excess degree of a finite degree distribution."""
    m = sum(k * pr for k, pr in d.items())
    return finite_pgf({k - 1: k * pr / m for k, pr in d.items() if k > 0})


# ------------------------------------------------------------ (1) degeneracy
def figure_degeneracy():
    plt = _mpl()
    p, q = symbols('p q')
    dists = ((r'$P_A$', {0: 0.5, 4: 0.5}, DARK, '-'),
             (r'$P_B$', {0: 0.2, 1: 0.5, 5: 0.3}, MID, '-'))
    fig, ax = plt.subplots(figsize=(4.3, 2.9))
    for name, d, col, ls in dists:
        G = hypergraph_giant(degree=finite_pgf(d), excess_degree=excess_pgf(d),
                             graph=True, poisson=False)
        ps = np.linspace(0.05, 1.0, 96)
        ax.plot(ps, [G.node_fraction({'p': v, 'q': 1}) for v in ps], ls,
                lw=1.4, color=col, label=name)
        C = CriticalAmplitude(G)
        Bv = float(C.amplitude_at_threshold(0, p).subs(q, 1))
        lam = C.perron_root().subs(q, 1)
        pl = np.linspace(1 / 3, 0.62, 40)
        ax.plot(pl, [Bv * (float(lam.subs(p, v)) - 1) for v in pl], '--',
                lw=1.0, color=col, dashes=(4, 2))
        print(f'  {name}: theta = {G.theta().subs(q, 1)},  B = {Bv:.4f},  '
              f'S(p=1) = {G.node_fraction({"p": 1, "q": 1}):.4f}')
    ax.axvline(1 / 3, ls=':', lw=0.8, color='0.6')
    ax.text(0.355, 0.56, r'$p_c=1/3$, both', fontsize=7.5, color='0.35')
    ax.text(0.99, 0.30, r'dashed: $B\Lambda$', fontsize=7.5, color='0.35',
            ha='right')
    ax.set_xlabel(r'site occupation probability $p$')
    ax.set_ylabel(r'giant component $S$')
    ax.legend(frameon=False, fontsize=8.5, loc='upper left')
    _tidy(ax); fig.tight_layout(); fig.savefig(OUT / 'fig-degeneracy.pdf')


# ---------------------------------------------------------------- (2) AND/OR
def figure_andor():
    plt = _mpl()
    p = symbols('p')
    fig, ax = plt.subplots(figsize=(4.3, 2.9))
    for logic, col, mark in (('or', MID, '-'), ('and', DARK, '-')):
        G = and_or_hypergraph(logic=logic, p=p)
        Gs = G.__class__(G.phi, G.phibar, G.g, G.gbar,
                         root_occupation=G.root_occupation)
        sub = {'k': 2, 'c': 3}
        ps = np.linspace(0.02, 1.0, 99)
        S = [Gs.node_fraction({**sub, 'p': v}) for v in ps]
        pc = float(nsolve(Gs.theta().subs(sub), p, 0.3 if logic == 'or' else 0.6))
        ax.plot(ps, S, mark, lw=1.4, color=col,
                label=f'{logic.upper()}-logic  ($p_c={pc:.4f}$)')
        ax.axvline(pc, ls=':', lw=0.8, color=col)
        print(f'  {logic.upper()}: p_c = {pc:.4f},  S(p=1) = {S[-1]:.4f}')
    ax.set_xlabel(r'site occupation probability $p$')
    ax.set_ylabel(r'giant component $S$')
    ax.legend(frameon=False, fontsize=8, loc='upper left')
    _tidy(ax); fig.tight_layout(); fig.savefig(OUT / 'fig-andor.pdf')


# ----------------------------------------------------------- (3) correlation
HALF = Rational(1, 2)
PAIRS = {'correlated': [(1, 0), (3, 2)], 'anti-correlated': [(1, 2), (3, 0)]}
PHI = {'correlated': lambda x: HALF * x[1] + HALF * x[1]**3 * x[2]**2,
       'anti-correlated': lambda x: HALF * x[1] * x[2]**2 + HALF * x[1]**3,
       'independent': lambda x: (HALF * x[1] + HALF * x[1]**3)
                                * (HALF + HALF * x[2]**2)}


def _triangle_parts():
    q = symbols('q')
    tb = lambda y: ((3 * q**2 - 2 * q**3) * y[0]**2 + 2 * q * (1 - q)**2 * y[0]
                    + (1 - q)**3 + q * (1 - q)**2)
    Gbar = [[None] * 3 for _ in range(3)]
    Gbar[1][0] = lambda y: (1 - q) + q * y[0]
    Gbar[2][0] = tb
    return [None, lambda y: (1 - q) * y[0] + q * y[0]**2,
            lambda y: y[0] * tb(y)], Gbar


def mc_joint(kind, n, qv, reps=2):
    """Configuration model with prescribed joint (link, triangle) degrees."""
    out = []
    for _ in range(reps):
        on = np.zeros(n, bool); on[RNG.permutation(n)[:n // 2]] = True
        if kind == 'independent':
            on2 = np.zeros(n, bool); on2[RNG.permutation(n)[:n // 2]] = True
            st, co = np.where(on, 3, 1), np.where(on2, 2, 0)
        else:
            (a1, b1), (a2, b2) = PAIRS[kind]
            st, co = np.where(on, a2, a1), np.where(on, b2, b1)
        st, co = st.astype(int), co.astype(int)
        if st.sum() % 2:
            st[0] += 1
        L = RNG.permutation(np.repeat(np.arange(n), st)).reshape(-1, 2)
        T = RNG.permutation(np.repeat(np.arange(n), co)).reshape(-1, 3)
        e = np.vstack([L, T[:, [0, 1]], T[:, [1, 2]], T[:, [0, 2]]])
        e = e[RNG.random(len(e)) < qv]
        g = coo_matrix((np.ones(len(e)), (e[:, 0], e[:, 1])), shape=(n, n))
        _, lab = connected_components(g, directed=False)
        out.append(np.bincount(lab).max() / n)
    return float(np.mean(out)), float(np.std(out))


def figure_correlated(mc=True):
    plt = _mpl()
    q = symbols('q')
    G, Gbar = _triangle_parts()
    style = {'correlated': (DARK, 'o'), 'anti-correlated': (LIGHT, '^'),
             'independent': (MID, 's')}
    fig, ax = plt.subplots(figsize=(4.3, 3.1))
    qs = np.linspace(0.10, 1.0, 91)
    mq = [0.25, 0.4, 0.6, 1.0]
    for name, Phi in PHI.items():
        col, mk = style[name]
        M = JointGiantComponent(Phi=[Phi, None, None], G=G, Gbar=Gbar)
        qc = float(nsolve(M.theta(), q, 0.3))
        ax.plot(qs, [M.node_fraction({'q': v}) for v in qs], '-', lw=1.4,
                color=col, label=f'{name} ($q_c={qc:.4f}$)')
        ax.axvline(qc, ls=':', lw=0.7, color=col)
        if mc:
            mm = [mc_joint(name, 150000, v, reps=2) for v in mq]
            ax.errorbar(mq, [m[0] for m in mm], yerr=[m[1] for m in mm],
                        fmt=mk, ms=4, color=col, mfc='white', mew=0.9, lw=0.8)
        print(f'  {name}: q_c = {qc:.4f},  S(q=1) = '
              f'{M.node_fraction({"q": 1.0}):.4f}')
    ax.set_xlabel(r'bond occupation probability $q$')
    ax.set_ylabel(r'giant component $S$')
    ax.legend(frameon=False, fontsize=7.5, loc='lower right')
    _tidy(ax); fig.tight_layout(); fig.savefig(OUT / 'fig-correlated.pdf')


def check_cardinality_correlation():
    """Sec. 5.6.2: hyperdegree-cardinality correlation at fixed marginals.

    Two cardinality classes, c = 2 and c = 5, and a one-parameter family in
    which a node is one of two equally likely types; ``w`` pairs the high mean
    in one class with the high or the low mean in the other.  Every marginal,
    and so every entry of the published threshold tensor, is the same for every
    ``w``.
    """
    from sympy import Rational, nsolve
    p_ = symbols('p')
    card, ma, mb, sp = [2, 5], Rational(1), Rational(3, 5), Rational(4, 5)
    print('    w      p_c      S(1/2)    S(1)')
    for w in (0, Rational(1, 2), 1):
        M = correlated_cardinality_hypergraph(card,
                                              two_class_joint_degree(ma, mb, sp, w))
        pc = float(nsolve(M.theta(), p_, 0.2))
        print(f'  {float(w):<5}  {pc:.4f}   {M.node_fraction({"p": 0.5}):.4f}'
              f'    {M.node_fraction({"p": 1.0}):.4f}')
    print('    positive correlation lowers the threshold; the order parameter')
    print('    is not ordered the same way at every p, so the curves cross')


def _tidy(ax):
    ax.tick_params(labelsize=8)
    for sp in ('top', 'right'):
        ax.spines[sp].set_visible(False)


if __name__ == '__main__':
    print('cardinality correlation:'); check_cardinality_correlation()
    print('degeneracy:');   figure_degeneracy()
    print('AND/OR:');       figure_andor()
    print('correlation:');  figure_correlated()
    print('wrote three PDFs')
