"""Figures for manuscript_3: the giant component in chygraphs."""

import math
import numpy as np
import matplotlib.pyplot as plt
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
from sympy import symbols, nsolve, Rational

from percolation.giant import hypergraph_giant, graph_with_triangles_giant, finite_pgf
from percolation.joint import JointGiantComponent
from percolation.amplitude import CriticalAmplitude
from percolation.applications import (and_or_hypergraph, household_epidemic,
                                   correlated_cardinality_hypergraph,
                                   two_class_joint_degree)

plt.rcParams.update({'font.size': 9, 'axes.linewidth': 0.8,
                     'xtick.direction': 'in', 'ytick.direction': 'in',
                     'figure.dpi': 160})
rng = np.random.default_rng(11)


def mc_triangles(n, kL, kT, q, reps=3):
    """Configuration model with Poisson link stubs and triangle corners."""
    out = []
    for _ in range(reps):
        st = rng.poisson(kL, n); st[0] += st.sum() % 2
        co = rng.poisson(kT, n); co[0] += (-co.sum()) % 3
        L = rng.permutation(np.repeat(np.arange(n), st)).reshape(-1, 2)
        T = rng.permutation(np.repeat(np.arange(n), co)).reshape(-1, 3)
        e = np.vstack([L, T[:, [0, 1]], T[:, [1, 2]], T[:, [0, 2]]])
        e = e[rng.random(len(e)) < q]
        g = coo_matrix((np.ones(len(e)), (e[:, 0], e[:, 1])), shape=(n, n))
        _, lab = connected_components(g, directed=False)
        out.append(np.bincount(lab).max() / n)
    return np.mean(out), np.std(out)


def figure_triangles():
    """Order parameter for a graph with over-represented triangles (Sec. III D)."""
    q, kL, kT = symbols('q k_L k_T')
    T = graph_with_triangles_giant()
    C = CriticalAmplitude(T)
    base = {kL: 1, kT: 0.5}
    qc = float(nsolve(C.Lambda().subs(base), q, 0.4))

    fig, ax = plt.subplots(figsize=(3.4, 2.6))
    qs = np.linspace(0.30, 1.0, 71)
    ax.plot(qs, [T.node_fraction({'k_L': 1, 'k_T': 0.5, 'q': v}) for v in qs],
            '-', lw=1.4, color='C0', label='chygraph map')
    mq = [0.35, 0.45, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
    mm = [mc_triangles(300000, 1.0, 0.5, v) for v in mq]
    ax.errorbar(mq, [m[0] for m in mm], yerr=[m[1] for m in mm], fmt='o',
                ms=3.5, color='k', mfc='none', lw=0.8, label='simulation')
    ax.axvline(qc, ls=':', lw=0.8, color='0.5')
    ax.text(qc + 0.014, 0.30, r'$q_c$', fontsize=8, color='0.4')
    ax.set_xlabel(r'bond occupation $q$')
    ax.set_ylabel(r'$S$')
    ax.legend(frameon=False, fontsize=8, loc='upper left')
    fig.tight_layout()
    fig.savefig('fig_triangles.pdf'); fig.savefig('fig_triangles.png')
    print(f'fig_triangles  q_c={qc:.6f}')


def figure_amplitude():
    """S/Lambda approaching the closed-form amplitude (Sec. V D)."""
    q, kL, kT = symbols('q k_L k_T')
    T = graph_with_triangles_giant()
    C = CriticalAmplitude(T)
    base = {kL: 1, kT: 0.5}
    qc = float(nsolve(C.Lambda().subs(base), q, 0.4))
    B = float(C.amplitude(0).subs({**base, q: qc}))

    lam, ratio = [], []
    for dq in np.geomspace(2e-4, 8e-2, 26):
        v = qc + dq
        s = T.node_fraction({'k_L': 1, 'k_T': 0.5, 'q': v},
                            tol=1e-16, maxiter=4_000_000)
        Lm = float(C.perron_root().subs({**base, q: v})) - 1
        lam.append(Lm); ratio.append(s / Lm)

    fig, ax = plt.subplots(figsize=(3.4, 2.6))
    ax.semilogx(lam, ratio, 'o-', ms=3, lw=1.0, color='C0')
    ax.axhline(B, ls='--', lw=1.0, color='C3')
    ax.text(3e-4, B + 0.015, rf'$B={B:.4f}$', fontsize=8, color='C3')
    ax.set_xlabel(r'$\Lambda$')
    ax.set_ylabel(r'$S/\Lambda$')
    fig.tight_layout()
    fig.savefig('fig_amplitude.pdf'); fig.savefig('fig_amplitude.png')
    print(f'fig_amplitude  q_c={qc:.6f}  B={B:.6f}')


def figure_degeneracy():
    """Two degree distributions with identical first moments."""
    A = {0: 0.5, 4: 0.5}
    Bd = {0: 0.2, 1: 0.5, 5: 0.3}

    def excess(d):
        m = sum(k * p for k, p in d.items())
        return finite_pgf({k - 1: k * p / m for k, p in d.items() if k > 0})

    fig, ax = plt.subplots(figsize=(3.4, 2.6))
    p = symbols('p')
    for (name, d, col) in ((r'$P_A$', A, 'C0'), (r'$P_B$', Bd, 'C1')):
        G = hypergraph_giant(degree=finite_pgf(d), excess_degree=excess(d),
                             graph=True, poisson=False)
        ps = np.linspace(0.05, 1.0, 96)
        S = [G.node_fraction({'p': v, 'q': 1}) for v in ps]
        ax.plot(ps, S, '-', lw=1.4, color=col, label=name)
        # leading order S = B * Lambda from the closed-form amplitude
        Cc = CriticalAmplitude(G)
        Bv = float(Cc.amplitude_at_threshold(0, p).subs(symbols('q'), 1))
        lam = Cc.perron_root().subs(symbols('q'), 1)
        pl = np.linspace(1 / 3, 0.62, 40)
        ax.plot(pl, [Bv * (float(lam.subs(p, v)) - 1) for v in pl],
                '--', lw=1.0, color=col, dashes=(4, 2))
    ax.axvline(1 / 3, ls=':', lw=0.8, color='0.5')
    ax.text(1.0, 0.055, r'$p_c=1/3$ (both);  dashed: $B\Lambda$',
            fontsize=7.5, color='0.4', ha='right')
    ax.set_xlabel(r'site occupation $p$')
    ax.set_ylabel(r'$S$')
    ax.legend(frameon=False, fontsize=8, loc='upper left')
    fig.tight_layout()
    fig.savefig('fig_degeneracy.pdf'); fig.savefig('fig_degeneracy.png')

    for name, d in (('A', A), ('B', Bd)):
        G = hypergraph_giant(degree=finite_pgf(d), excess_degree=excess(d),
                             graph=True, poisson=False)
        Cc = CriticalAmplitude(G)
        print(f'fig_degeneracy  {name}: theta={G.theta()}  '
              f'B={float(Cc.amplitude_at_threshold(0, symbols("p")).subs(symbols("q"), 1)):.4f}  '
              f'S(p=1)={G.node_fraction({"p":1,"q":1}):.6f}')


HALF = Rational(1, 2)
PAIRS = {'correlated': [(1, 0), (3, 2)], 'anti-correlated': [(1, 2), (3, 0)]}
PHI = {'correlated':     lambda x: HALF*x[1] + HALF*x[1]**3*x[2]**2,
       'anti-correlated': lambda x: HALF*x[1]*x[2]**2 + HALF*x[1]**3,
       'independent':    lambda x: (HALF*x[1] + HALF*x[1]**3)*(HALF + HALF*x[2]**2)}


def _triangle_parts():
    from sympy import symbols
    q = symbols('q')
    tb = lambda y: ((3*q**2 - 2*q**3)*y[0]**2 + 2*q*(1-q)**2*y[0]
                    + (1-q)**3 + q*(1-q)**2)
    Gbar = [[None]*3 for _ in range(3)]
    Gbar[1][0] = lambda y: (1-q) + q*y[0]
    Gbar[2][0] = tb
    return [None, lambda y: (1-q)*y[0] + q*y[0]**2, lambda y: y[0]*tb(y)], Gbar


def mc_joint(kind, n, qv, reps=2):
    """Configuration model with prescribed joint (link, triangle) degrees."""
    out = []
    for _ in range(reps):
        on = np.zeros(n, bool); on[rng.permutation(n)[:n//2]] = True
        if kind == 'independent':
            on2 = np.zeros(n, bool); on2[rng.permutation(n)[:n//2]] = True
            st, co = np.where(on, 3, 1), np.where(on2, 2, 0)
        else:
            (a1, b1), (a2, b2) = PAIRS[kind]
            st, co = np.where(on, a2, a1), np.where(on, b2, b1)
        st, co = st.astype(int), co.astype(int)
        if st.sum() % 2: st[0] += 1
        L = rng.permutation(np.repeat(np.arange(n), st)).reshape(-1, 2)
        T = rng.permutation(np.repeat(np.arange(n), co)).reshape(-1, 3)
        e = np.vstack([L, T[:, [0,1]], T[:, [1,2]], T[:, [0,2]]])
        e = e[rng.random(len(e)) < qv]
        g = coo_matrix((np.ones(len(e)), (e[:,0], e[:,1])), shape=(n,n))
        _, lab = connected_components(g, directed=False)
        out.append(np.bincount(lab).max()/n)
    return np.mean(out), np.std(out)


def figure_correlated():
    """Correlated layers: identical marginals, three different transitions."""
    from sympy import symbols
    q = symbols('q')
    G, Gbar = _triangle_parts()
    fig, ax = plt.subplots(figsize=(3.4, 2.6))
    cols = {'correlated': 'C0', 'anti-correlated': 'C3', 'independent': 'C2'}
    qs = np.linspace(0.10, 1.0, 91)
    mq = [0.2, 0.3, 0.4, 0.5, 0.7, 1.0]
    for name, Phi in PHI.items():
        M = JointGiantComponent(Phi=[Phi, None, None], G=G, Gbar=Gbar)
        qc = float(nsolve(M.theta(), q, 0.3))
        ax.plot(qs, [M.node_fraction({'q': v}) for v in qs], '-', lw=1.4,
                color=cols[name], label=f'{name}  ($q_c={qc:.3f}$)')
        mm = [mc_joint(name, 300000, v) for v in mq]
        ax.errorbar(mq, [m[0] for m in mm], yerr=[m[1] for m in mm], fmt='o',
                    ms=3.5, color=cols[name], mfc='none', lw=0.8)
        ax.axvline(qc, ls=':', lw=0.7, color=cols[name])
        print(f'fig_correlated  {name}: q_c={qc:.6f}  '
              f'S(1)={M.node_fraction({"q":1.0}):.5f}  '
              f'<kbar^(1)>_02={M.kappa_bar(1,0,2)}')
    ax.set_xlabel(r'bond occupation $q$')
    ax.set_ylabel(r'$S$')
    ax.legend(frameon=False, fontsize=7, loc='lower right')
    fig.tight_layout()
    fig.savefig('fig_correlated.pdf'); fig.savefig('fig_correlated.png')


def _giant(n, e):
    if len(e) == 0:
        return 0.0
    e = np.asarray(e)
    g = coo_matrix((np.ones(len(e)), (e[:, 0], e[:, 1])), shape=(n, n))
    _, lab = connected_components(g, directed=False)
    return np.bincount(lab).max() / n


def _pairs(m):
    a, b = np.triu_indices(len(m), 1)
    return np.stack([m[a], m[b]], axis=1)


def mc_andor(n, kmean, cmean, p, logic, reps=2):
    out = []
    for _ in range(reps):
        deg = rng.poisson(kmean, n)
        m = max(1, int(round(deg.sum() / cmean)))
        card = rng.poisson(cmean, m)
        d, c = deg.sum(), card.sum()
        if d > c:
            card[rng.integers(0, m, d - c)] += 1
        elif c > d:
            idx = rng.permutation(np.repeat(np.arange(m), card))[:d]
            card = np.bincount(idx, minlength=m)
        stubs = rng.permutation(np.repeat(np.arange(n), deg))
        present = rng.random(n) < p
        ed, at = [], 0
        for j in range(m):
            mem = stubs[at:at + card[j]]; at += card[j]
            if len(mem) < 2:
                continue
            if logic == 'and':
                if present[mem].all():
                    ed.append(_pairs(mem))
            else:
                g = mem[present[mem]]
                if len(g) >= 2:
                    ed.append(_pairs(g))
        out.append(_giant(n, np.vstack(ed) if ed else np.empty((0, 2), int)))
    return np.mean(out), np.std(out)


def mc_corr(n, w, p, card, ma, mb, sp, reps=2):
    ap, am = float(ma * (1 + sp)), float(ma * (1 - sp))
    bp, bm = float(mb * (1 + sp)), float(mb * (1 - sp))
    out = []
    for _ in range(reps):
        hi = rng.random(n) < 0.5
        same = rng.random(n) < w
        d1 = rng.poisson(np.where(hi, ap, am))
        d2 = rng.poisson(np.where(hi == same, bp, bm))
        present = rng.random(n) < p
        ed = []
        for deg, c in ((d1, card[0]), (d2, card[1])):
            s = deg.sum(); s -= s % c
            grp = rng.permutation(np.repeat(np.arange(n), deg))[:s].reshape(-1, c)
            for gg in grp:
                a = gg[present[gg]]
                if len(a) >= 2:
                    ed.append(_pairs(a))
        out.append(_giant(n, np.vstack(ed) if ed else np.empty((0, 2), int)))
    return np.mean(out), np.std(out)


def mc_house(n, hsize, kmean, pH, T, reps=2):
    out = []
    for _ in range(reps):
        perm = rng.permutation((n // hsize) * hsize).reshape(-1, hsize)
        ed = []
        for h in perm:
            pr = _pairs(h)
            ed.append(pr[rng.random(len(pr)) < pH])
        deg = rng.poisson(kmean, n)
        if deg.sum() % 2:
            deg[0] += 1
        st = rng.permutation(np.repeat(np.arange(n), deg)).reshape(-1, 2)
        ed.append(st[rng.random(len(st)) < T])
        out.append(_giant(n, np.vstack(ed)))
    return np.mean(out), np.std(out)


def figure_applications():
    """Three literature constructions mapped to chygraphs."""
    fig, ax = plt.subplots(1, 3, figsize=(6.9, 2.5))

    # (a) AND vs OR hypergraph site percolation
    ps = np.linspace(0.05, 1.0, 49)
    for logic, col, lab in (('or', 'C0', 'OR (factor graph)'),
                            ('and', 'C3', 'AND (hypergraph)')):
        M = and_or_hypergraph(logic)
        ax[0].plot(ps, [M.node_fraction({'k': 2, 'c': 3, 'p': v}) for v in ps],
                   '-', lw=1.4, color=col, label=lab)
        mp = [0.3, 0.5, 0.7, 0.9, 1.0]
        mm = [mc_andor(150000, 2.0, 3.0, v, logic) for v in mp]
        ax[0].errorbar(mp, [m[0] for m in mm], yerr=[m[1] for m in mm], fmt='o',
                       ms=3, color=col, mfc='none', lw=0.8)
    ax[0].set_xlabel(r'node occupation $p$'); ax[0].set_ylabel(r'$S$')
    ax[0].legend(frameon=False, fontsize=6.5, loc='upper left')
    ax[0].set_title('(a)', loc='left', fontsize=9)

    # (b) hyperdegree-cardinality correlation, identical marginals
    card, ma, mb, sp = [2, 5], Rational(1), Rational(3, 5), Rational(4, 5)
    ps = np.linspace(0.10, 1.0, 46)
    for w, col, lab in ((0, 'C3', r'$w=0$ (anti)'),
                        (Rational(1, 2), 'C2', r'$w=1/2$ (indep.)'),
                        (1, 'C0', r'$w=1$ (corr.)')):
        M = correlated_cardinality_hypergraph(card, two_class_joint_degree(ma, mb, sp, w))
        ax[1].plot(ps, [M.node_fraction({'p': v}) for v in ps], '-', lw=1.4,
                   color=col, label=lab)
        mp = [0.3, 0.5, 0.7, 1.0]
        mm = [mc_corr(150000, float(w), v, card, ma, mb, sp) for v in mp]
        ax[1].errorbar(mp, [m[0] for m in mm], yerr=[m[1] for m in mm], fmt='o',
                       ms=3, color=col, mfc='none', lw=0.8)
        print(f'fig_applications  corr w={float(w)}: S(0.5)={M.node_fraction({"p":0.5}):.5f} '
              f'S(1)={M.node_fraction({"p":1.0}):.5f}')
    ax[1].set_xlabel(r'node occupation $p$'); ax[1].set_ylabel(r'$S$')
    ax[1].legend(frameon=False, fontsize=6.5, loc='upper left')
    ax[1].set_title('(b)', loc='left', fontsize=9)

    # (c) households plus a global network
    H = household_epidemic({3: 1})
    Ts = np.linspace(0.02, 0.8, 40)
    for pH, col in ((0.0, 'C2'), (0.5, 'C0'), (0.9, 'C3')):
        ax[2].plot(Ts, [H.node_fraction({'p_H': pH, 'T': v, 'k': 2.0}) for v in Ts],
                   '-', lw=1.4, color=col, label=rf'$p_H={pH}$')
        mp = [0.1, 0.2, 0.3, 0.5, 0.7]
        mm = [mc_house(150000, 3, 2.0, pH, v) for v in mp]
        ax[2].errorbar(mp, [m[0] for m in mm], yerr=[m[1] for m in mm], fmt='o',
                       ms=3, color=col, mfc='none', lw=0.8)
    ax[2].set_xlabel(r'global transmission $T$'); ax[2].set_ylabel(r'$S$')
    ax[2].legend(frameon=False, fontsize=6.5, loc='lower right')
    ax[2].set_title('(c)', loc='left', fontsize=9)

    fig.tight_layout()
    fig.savefig('fig_applications.pdf'); fig.savefig('fig_applications.png')


if __name__ == '__main__':
    import os
    os.chdir(os.path.dirname(os.path.abspath(__file__)))
    figure_triangles()
    figure_amplitude()
    figure_degeneracy()
    figure_correlated()
    figure_applications()
