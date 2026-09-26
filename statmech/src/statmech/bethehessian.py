"""The Bethe Hessian of a chygraph: the linearised recursion, the weighted
Ihara--Bass identity on the incidence graph, and what it needs.

Sec. 24.6 of the book says the linearised chygraph recursion is a weighted
non-backtracking operator on the incidence graph, that the Bethe Hessian
is the symmetric operator on the same graph, and that the two are tied by a
graph zeta function -- Watanabe and Fukumizu's weighted Ihara--Bass identity
(NIPS 22, 2009), which Saade, Krzakala and Zdeborova (NIPS 27, 2014) use to
read the stability of a fixed point off the Bethe Hessian on a graph.  This
module carries that out on the factor graph of alpha + I.

The recursion in field coordinates.  A message is a field,
``h = (ln m(+) - ln m(-)) / 2``.  A variable relays, ``h_{v->a} = sum_{b != a}
h_{b->v}``, derivative one on every leg.  A complex ``a`` sends
``h_{a->v}`` as a function of the fields on its other members, and

    M^a_{vu} = d h_{a->v} / d h_{u->a}
             = [<s_u>_{s_v=+} - <s_u>_{s_v=-}] / 2
             = Cov_a(s_v, s_u) / (1 - m_v^2),                     (u != v)

the covariance under the factor belief ``b_a``.  The linearised map ``T`` on
directed inclusions has ``T[(a->v), (u->a)] = M^a_{vu}`` and
``T[(v->a), (b->v)] = 1`` for ``b != a``: non-backtracking on the incidence
graph, with a *block* on each complex.  At the trivial fixed point of a
homogeneous complex every entry of the block is the one number ``u'(c)`` of
Eq. (8.x), :func:`statmech.ising.clique_derivative`.

When the block factorises.  ``T`` is similar to a non-backtracking operator
with one weight per *directed edge*, ``w_{a->v} = alpha^a_v`` and
``w_{v->a} = beta^a_v``, exactly when every block has the product form
``M^a_{vu} = alpha^a_v beta^a_u`` off its diagonal.  That is the
factorisation ``F(u->a) -> F(a) -> F(a->v)`` through a one-dimensional
stalk at the complex that the book names as the thing to check.  It holds
for any block of size two or three (six off-diagonal entries, the one
cycle condition ``M_12 M_23 M_31 = M_13 M_32 M_21`` satisfied because
``M = D C`` with ``C`` symmetric), and at the trivial fixed point of a
homogeneous complex of any size.  From cardinality four on it is a genuine
condition, ``C_12 C_34 = C_13 C_24 = C_14 C_23`` on the belief's covariance,
and :func:`factorise` measures how far a block is from it.

The identity.  With factorised weights, on the incidence graph
``I = V + C`` (a vertex per variable and per complex, an edge per inclusion),

    det(I - T) = prod_{(v,a)} (1 - alpha^a_v beta^a_v) det(I + D - A),
    A_{va} = beta^a_v / (1 - alpha^a_v beta^a_v),
    A_{av} = alpha^a_v / (1 - alpha^a_v beta^a_v),
    D_{vv} = sum_{a ∋ v} alpha^a_v beta^a_v / (1 - alpha^a_v beta^a_v),
    D_{aa} = sum_{v ∈ a} the same,

and ``H = I + D - A`` is the Bethe Hessian of the chygraph.  It is similar
to a symmetric matrix whenever every ``alpha beta`` is positive (conjugate
by the square roots), and eliminating the complex vertices -- ``H_CC`` is
diagonal, since no complex of Part III is a member of another -- gives an
operator on the variables alone.  At the trivial fixed point, with
``u = u'(c_a)``,

    H_V = I + sum_a u/(1-u) [ P_a - 1_a 1_a^T / (1 + (c_a - 1) u) ],
    det(I - T) = prod_a (1-u)^{c_a-1} (1 + (c_a-1) u) det(H_V),

``P_a`` the projector on the members of ``a``: each complex couples its
members all-to-all through one rank-one term.  For a graph, ``c_a = 2`` and
``u = t = tanh(beta J)``, this is ``t^2/(1-t^2)`` times Saade's
``(r^2 - 1) I - r A + D`` at ``r = 1/t``.

What follows from it.  ``H(beta) = I`` at ``beta = 0``; on a ferromagnet
``T`` is non-negative and its Perron root is real, so ``det(I - T) > 0``
until the root reaches one, and by continuity no eigenvalue of ``H``
crosses zero before that: the trivial fixed point of the recursion is
stable iff the chygraph Bethe Hessian is positive definite, and the
instance's threshold is where its gap closes.  The ensemble version is
Eq. (8.branch): the Perron root of ``T`` on a random chygraph converges to
the Perron root of the branching matrix, and so does the point where the
gap closes.
"""

from itertools import combinations

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
from scipy.special import logsumexp

from statmech.loopseries import SPIN, BinaryFactorGraph


# ---------------------------------------------------------------------------
# instances
# ---------------------------------------------------------------------------

def random_chygraph(n, cardinalities, means, rng, regular=False):
    """A random chygraph on ``n`` atoms with one Poisson layer per
    cardinality: chy-degrees drawn per atom, stubs matched at random,
    complexes with a repeated atom dropped.  ``regular=True`` gives every
    atom exactly ``means[l]`` complexes of each layer, the matching redrawn
    until no complex repeats an atom, so that the chygraph is exactly
    regular.  Returns a list of tuples."""
    complexes = []
    for c, kappa in zip(cardinalities, means):
        if regular:
            deg = np.full(n, int(round(kappa)))
        else:
            deg = rng.poisson(kappa, n)
        stubs = np.repeat(np.arange(n), deg)
        for attempt in range(1000):
            rng.shuffle(stubs)
            cut = stubs[: (len(stubs) // c) * c].reshape(-1, c)
            good = [tuple(sorted(int(v) for v in row)) for row in cut
                    if len(set(row.tolist())) == c]
            if not regular or len(good) == len(cut):
                break
        else:
            raise RuntimeError('no regular matching without a repeated atom')
        complexes.extend(good)
    return complexes


def clique_factor_graph(complexes, beta_J, field=0.0, fields=None):
    """The factor graph with one full-clique factor per complex (the
    recursion's own model; a bond shared by two complexes is counted in
    both).  ``field`` is a uniform external field, ``fields`` a per-atom
    one, each spread evenly over the factors at the atom."""
    cx = [tuple(sorted(int(v) for v in c)) for c in complexes]
    deg = {}
    for c in cx:
        for v in c:
            deg[v] = deg.get(v, 0) + 1
    factors = []
    for c in cx:
        n = len(c)
        s = np.array([[1.0 if (i >> (n - 1 - b)) & 1 == 0 else -1.0
                       for b in range(n)] for i in range(2 ** n)])
        e = np.zeros(2 ** n)
        for p, q in combinations(range(n), 2):
            e += s[:, p] * s[:, q]
        table = (beta_J * e).reshape((2,) * n)
        for p in range(n):
            hv = field + (0.0 if fields is None else float(fields[c[p]]))
            if hv:
                sh = [1] * n
                sh[p] = 2
                table = table + (hv * SPIN).reshape(sh) / deg[c[p]]
        factors.append((c, table))
    return BinaryFactorGraph(factors)


# ---------------------------------------------------------------------------
# the linearised recursion
# ---------------------------------------------------------------------------

def _field(m):
    return 0.5 * float(m[0] - m[1])


def jacobian_block(fg, a):
    """``M^a_{vu} = d h_{a->v} / d h_{u->a}`` at the current messages, as a
    ``c x c`` array with zero diagonal, by enumeration of the factor."""
    sc, t = fg.factors[a]
    n = len(sc)
    M = np.zeros((n, n))
    for i, v in enumerate(sc):
        acc = t.copy()
        for j, u in enumerate(sc):
            if u == v:
                continue
            sh = [1] * n
            sh[j] = 2
            acc = acc + fg.m_va[(u, a)].reshape(sh)
        # conditional on s_v = +1 (index 0) and -1 (index 1)
        for j, u in enumerate(sc):
            if u == v:
                continue
            means = []
            for state in (0, 1):
                sub = np.take(acc, state, axis=i)
                w = np.exp(sub - logsumexp(sub))
                jj = j if j < i else j - 1
                ax = tuple(k for k in range(n - 1) if k != jj)
                p = w.sum(axis=ax) if ax else w
                means.append(float(p @ SPIN))
            M[i, j] = 0.5 * (means[0] - means[1])
    return M


def belief_covariance(fg, a):
    """The covariance matrix of the factor belief ``b_a``; ``M^a`` is its
    off-diagonal part divided by ``1 - m_v^2`` row by row."""
    sc, _ = fg.factors[a]
    b = fg.factor_belief(a)
    n = len(sc)
    m = np.zeros(n)
    C = np.zeros((n, n))
    for i in range(n):
        ax = tuple(k for k in range(n) if k != i)
        m[i] = float(b.sum(axis=ax) @ SPIN)
    for i, j in combinations(range(n), 2):
        ax = tuple(k for k in range(n) if k not in (i, j))
        pij = b.sum(axis=ax) if ax else b
        C[i, j] = C[j, i] = float(SPIN @ pij @ SPIN) - m[i] * m[j]
    for i in range(n):
        C[i, i] = 1.0 - m[i] ** 2
    return C, m


def directed_inclusions(fg):
    """Index the directed inclusions: ``(a, v, +1)`` is ``a -> v`` and
    ``(a, v, -1)`` is ``v -> a``."""
    idx = {}
    for a, (sc, _) in enumerate(fg.factors):
        for v in sc:
            idx[(a, v, +1)] = len(idx)
            idx[(a, v, -1)] = len(idx)
    return idx


def linearised_operator(fg, blocks=None):
    """The Jacobian ``T`` of the message map on directed inclusions, sparse.

    ``blocks`` overrides the complex blocks (a list of ``c x c`` arrays);
    by default they are computed at the current messages."""
    idx = directed_inclusions(fg)
    rows, cols, vals = [], [], []
    for a, (sc, _) in enumerate(fg.factors):
        M = jacobian_block(fg, a) if blocks is None else blocks[a]
        for i, v in enumerate(sc):
            for j, u in enumerate(sc):
                if u != v:
                    rows.append(idx[(a, v, +1)])
                    cols.append(idx[(a, u, -1)])
                    vals.append(M[i, j])
    for v in fg.nodes:
        for a in fg.of[v]:
            for b in fg.of[v]:
                if b != a:
                    rows.append(idx[(a, v, -1)])
                    cols.append(idx[(b, v, +1)])
                    vals.append(1.0)
    m = len(idx)
    return sp.csr_matrix((vals, (rows, cols)), shape=(m, m))


def spectral_radius(T, k=1):
    """The largest real part among the eigenvalues of ``T`` (its Perron root
    when ``T`` is non-negative); dense below 400, Arnoldi above."""
    if T.shape[0] <= 400:
        ev = np.linalg.eigvals(T.toarray())
        return float(ev.real.max())
    ev = spla.eigs(T, k=k, which='LR', return_eigenvectors=False, tol=1e-12)
    return float(ev.real.max())


def trivial_blocks(fg):
    """The blocks at the trivial fixed point (all messages zero)."""
    saved = dict(fg.m_va)
    for key in fg.m_va:
        fg.m_va[key] = np.zeros(2)
    blocks = [jacobian_block(fg, a) for a in range(len(fg.factors))]
    fg.m_va.update(saved)
    return blocks


# ---------------------------------------------------------------------------
# the factorisation through the complex
# ---------------------------------------------------------------------------

def factorise(M, iters=500, tol=1e-15):
    """Best product form ``M_vu ~ alpha_v beta_u`` on the off-diagonal
    entries, by alternating least squares from the leading singular pair.

    Returns ``(alpha, beta, residual)`` with the residual the Frobenius
    norm of the off-diagonal misfit over that of ``M``.  Exact (residual
    at rounding) for two and three members and at the trivial fixed point
    of a homogeneous complex."""
    M = np.asarray(M, float)
    n = M.shape[0]
    off = ~np.eye(n, dtype=bool)
    norm = np.linalg.norm(M[off])
    if n == 2:
        return np.array([M[0, 1], M[1, 0]]), np.ones(2), 0.0
    if norm == 0:
        return np.zeros(n), np.zeros(n), 0.0
    U, S, Vt = np.linalg.svd(M)
    alpha = U[:, 0] * np.sqrt(S[0])
    beta = Vt[0] * np.sqrt(S[0])
    last = np.inf
    for _ in range(iters):
        for v in range(n):
            b = np.delete(beta, v)
            alpha[v] = (np.delete(M[v], v) @ b) / (b @ b) if b @ b else 0.0
        for u in range(n):
            a = np.delete(alpha, u)
            beta[u] = (np.delete(M[:, u], u) @ a) / (a @ a) if a @ a else 0.0
        R = M - np.outer(alpha, beta)
        res = np.linalg.norm(R[off]) / norm
        if abs(last - res) < tol:
            break
        last = res
    return alpha, beta, res


def factorisation_residuals(fg, blocks=None):
    """The residual of :func:`factorise` on every complex."""
    out = []
    for a in range(len(fg.factors)):
        M = jacobian_block(fg, a) if blocks is None else blocks[a]
        out.append(factorise(M)[2])
    return np.array(out)


# ---------------------------------------------------------------------------
# the Bethe Hessian
# ---------------------------------------------------------------------------

def edge_weights(fg, blocks=None):
    """``alpha^a_v`` (on ``a -> v``) and ``beta^a_v`` (on ``v -> a``) from the
    factorised blocks; returns two dicts keyed by ``(a, v)`` and the largest
    factorisation residual met."""
    alpha, beta, worst = {}, {}, 0.0
    for a, (sc, _) in enumerate(fg.factors):
        M = jacobian_block(fg, a) if blocks is None else blocks[a]
        al, be, res = factorise(M)
        worst = max(worst, res)
        for i, v in enumerate(sc):
            alpha[(a, v)] = float(al[i])
            beta[(a, v)] = float(be[i])
    return alpha, beta, worst


def bethe_hessian(fg, alpha, beta, symmetric=False):
    """``H = I + D - A`` on the incidence graph (variables first, then the
    complexes), sparse, and the prefactor ``prod (1 - alpha beta)`` of the
    identity.  ``symmetric=True`` returns the conjugate symmetric matrix,
    which needs every ``alpha beta`` non-negative."""
    nv = len(fg.nodes)
    pos = {v: i for i, v in enumerate(fg.nodes)}
    n = nv + len(fg.factors)
    rows, cols, vals = [], [], []
    diag = np.ones(n)
    pref = 1.0
    for a, (sc, _) in enumerate(fg.factors):
        for v in sc:
            p = alpha[(a, v)] * beta[(a, v)]
            pref *= 1.0 - p
            diag[pos[v]] += p / (1.0 - p)
            diag[nv + a] += p / (1.0 - p)
            if symmetric:
                if p < 0:
                    raise ValueError('alpha beta < 0: no symmetric form')
                w = np.sqrt(p) / (1.0 - p)
                rows += [pos[v], nv + a]
                cols += [nv + a, pos[v]]
                vals += [-w, -w]
            else:
                rows += [pos[v], nv + a]
                cols += [nv + a, pos[v]]
                vals += [-beta[(a, v)] / (1.0 - p), -alpha[(a, v)] / (1.0 - p)]
    H = sp.csr_matrix((vals, (rows, cols)), shape=(n, n)) + sp.diags(diag)
    return H.tocsr(), pref


def vertex_hessian(fg, alpha, beta):
    """The Schur complement of the Bethe Hessian on the variables, and the
    determinant of the eliminated diagonal block ``H_CC``."""
    nv = len(fg.nodes)
    pos = {v: i for i, v in enumerate(fg.nodes)}
    diag = np.ones(nv)
    rows, cols, vals = [], [], []
    det_cc = 1.0
    for a, (sc, _) in enumerate(fg.factors):
        ps = np.array([alpha[(a, v)] * beta[(a, v)] for v in sc])
        hcc = 1.0 + np.sum(ps / (1.0 - ps))
        det_cc *= hcc
        for i, v in enumerate(sc):
            diag[pos[v]] += ps[i] / (1.0 - ps[i])
            for j, u in enumerate(sc):
                rows.append(pos[v])
                cols.append(pos[u])
                vals.append(-(beta[(a, v)] / (1.0 - ps[i]))
                            * (alpha[(a, u)] / (1.0 - ps[j])) / hcc)
    S = sp.csr_matrix((vals, (rows, cols)), shape=(nv, nv)) + sp.diags(diag)
    return S.tocsr(), det_cc


def trivial_vertex_hessian(fg, uprime):
    """The closed form at the trivial fixed point, ``uprime`` a list of
    ``u'(c_a)`` per complex, together with the prefactor
    ``prod (1-u)^{c-1} (1 + (c-1) u)``."""
    nv = len(fg.nodes)
    pos = {v: i for i, v in enumerate(fg.nodes)}
    diag = np.ones(nv)
    rows, cols, vals = [], [], []
    pref = 1.0
    for a, (sc, _) in enumerate(fg.factors):
        c, u = len(sc), uprime[a]
        pref *= (1.0 - u) ** (c - 1) * (1.0 + (c - 1) * u)
        w = u / (1.0 - u)
        for v in sc:
            diag[pos[v]] += w
            for x in sc:
                rows.append(pos[v])
                cols.append(pos[x])
                vals.append(-w / (1.0 + (c - 1) * u))
    S = sp.csr_matrix((vals, (rows, cols)), shape=(nv, nv)) + sp.diags(diag)
    return S.tocsr(), pref


def smallest_eigenvalue(H, k=1):
    """The smallest eigenvalue of a symmetric sparse matrix."""
    if H.shape[0] <= 400:
        return float(np.linalg.eigvalsh(H.toarray()).min())
    ev = spla.eigsh(H, k=k, which='SA', return_eigenvectors=False, tol=1e-12)
    return float(ev.min())


def log_det(A):
    """``ln |det A|`` and the sign, dense."""
    sign, ld = np.linalg.slogdet(A.toarray() if sp.issparse(A) else A)
    return float(sign), float(ld)


# ---------------------------------------------------------------------------
# thresholds on an instance
# ---------------------------------------------------------------------------

def instance_threshold(complexes, lo=0.02, hi=3.0, tol=1e-9, via='T'):
    """``beta J`` at which the trivial fixed point loses stability on one
    instance: ``via='T'`` bisects the Perron root of the linearised
    recursion through one, ``via='H'`` the smallest eigenvalue of the
    vertex Bethe Hessian through zero."""
    from statmech.ising import clique_derivative
    fg = clique_factor_graph(complexes, 1.0)
    cards = [len(sc) for sc, _ in fg.factors]

    def gap(bj):
        u = [clique_derivative(c, bj) for c in cards]
        if via == 'T':
            blocks = []
            for a, (sc, _) in enumerate(fg.factors):
                c = len(sc)
                blocks.append(u[a] * (np.ones((c, c)) - np.eye(c)))
            return spectral_radius(linearised_operator(fg, blocks)) - 1.0
        S, _ = trivial_vertex_hessian(fg, u)
        return smallest_eigenvalue(S)

    from scipy.optimize import brentq
    return brentq(gap, lo, hi, xtol=tol)


# ---------------------------------------------------------------------------
# the matrix-weighted identity: stalks of dimension d at the complexes
# ---------------------------------------------------------------------------
#
# Write a block as a completion minus its diagonal, ``M = R - diag(R)`` with
# ``R = A B^T`` of rank ``d``.  On the two spaces X (atom-to-complex fields)
# and Y (complex-to-atom fields), both indexed by the inclusions, the map is
# ``T = [[0, G], [F, 0]]`` with ``F = R - D`` (``D = diag R``) and
# ``G = P^T P - I``, ``P`` the atom incidence.  With ``E = (I - D)^{-1}`` the
# matrix determinant lemma gives
#
#     det(I - T) = det(I - D) det [[ I + P D E P^T,  -P E A     ],
#                                  [ -B^T E P^T,     I + B^T E A ]],
#
# a matrix on the atoms plus a ``d_a``-dimensional stalk at every complex.
# For a symmetric completion ``R = Q L Q^T`` take ``A = Q |L|^{1/2}`` and
# ``B = Q sign(L) |L|^{1/2}``; conjugating the stalk block by ``sign(L)``
# makes the matrix symmetric, with the signature of the completion as the
# metric on the stalk.  Rank one is the edge-weighted identity above.

from scipy.optimize import minimize


def symmetric_block(fg, a):
    """The block conjugated to symmetric form: the correlation matrix of the
    factor belief with its diagonal removed.  ``T`` built from these blocks
    is conjugate to the one built from :func:`jacobian_block`."""
    C, m = belief_covariance(fg, a)
    s = 1.0 / np.sqrt(1.0 - m ** 2)
    R = s[:, None] * C * s[None, :]
    np.fill_diagonal(R, 0.0)
    return R


def _tail_singular(Ms, d, x):
    sv = np.linalg.svd(Ms + np.diag(x), compute_uv=False)
    return float(np.sum(sv[d:] ** 2))


def min_rank_completion(Ms, d, rng=None, starts=12, tol=1e-18):
    """The diagonal ``x`` that brings ``Ms + diag(x)`` closest to rank ``d``
    (sum of squares of the trailing singular values), by multistart BFGS
    polished with Nelder--Mead; returns ``(x, misfit)`` with the misfit
    relative to ``|Ms|_F^2``."""
    rng = np.random.default_rng(0) if rng is None else rng
    n = Ms.shape[0]
    norm = float(np.sum(Ms ** 2))
    best, bx = np.inf, None
    f = lambda x: _tail_singular(Ms, d, x)  # noqa: E731
    for k in range(starts):
        x0 = np.zeros(n) if k == 0 else rng.normal(0, 0.5, n)
        r = minimize(f, x0, method='BFGS', options=dict(gtol=1e-14))
        r = minimize(f, r.x, method='Nelder-Mead',
                     options=dict(xatol=1e-13, fatol=1e-24, maxiter=4000))
        if r.fun < best:
            best, bx = r.fun, r.x
        if best / norm < tol:
            break
    return bx, best / norm


def stalk_dimension(Ms, tol=1e-12, rng=None):
    """The smallest ``d`` admitting a real diagonal completion of rank
    ``d`` (misfit below ``tol``), with that completion."""
    n = Ms.shape[0]
    for d in range(1, n + 1):
        x, mis = min_rank_completion(Ms, d, rng=rng)
        if mis < tol:
            return d, x
    raise RuntimeError('unreachable')


def completion_factors(R, d):
    """``A, B, S`` with ``R = A B^T`` of rank ``d`` from the symmetric
    eigendecomposition, ``S`` the signs of the kept eigenvalues."""
    w, Q = np.linalg.eigh(R)
    keep = np.argsort(-np.abs(w))[:d]
    w, Q = w[keep], Q[:, keep]
    A = Q * np.sqrt(np.abs(w))
    S = np.sign(w)
    return A, A * S, S


def matrix_bethe_hessian(fg, blocks, completions):
    """``H`` of the matrix-weighted identity and the prefactor
    ``prod (1 - r)``, from symmetric ``blocks`` and ``completions``: per
    complex a diagonal ``x`` and a rank ``d`` (``R = M + diag(x)``).  Also
    returns the block dimensions."""
    nv = len(fg.nodes)
    pos = {v: i for i, v in enumerate(fg.nodes)}
    dims = [d for _, d in completions]
    off = nv + np.concatenate([[0], np.cumsum(dims)[:-1]]).astype(int)
    N = nv + sum(dims)
    H = np.eye(N)
    pref = 1.0
    for a, (sc, _) in enumerate(fg.factors):
        x, d = completions[a]
        R = blocks[a] + np.diag(x)
        A, B, _ = completion_factors(R, d)
        r = np.diag(R)
        e = 1.0 / (1.0 - r)
        pref *= float(np.prod(1.0 - r))
        idx = [pos[v] for v in sc]
        sl = slice(off[a], off[a] + d)
        for i, v in enumerate(sc):
            H[idx[i], idx[i]] += r[i] * e[i]
            H[idx[i], sl] -= e[i] * A[i]
            H[sl, idx[i]] -= e[i] * B[i]
        H[sl, sl] += B.T @ (e[:, None] * A)
    return H, pref, dims
