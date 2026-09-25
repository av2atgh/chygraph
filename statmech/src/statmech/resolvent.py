"""The resolvent as a message (Chapter 18).

Every other module of this package passes a probability, a field, a survey or
a count.  Here the message is a Green's function.  For a quantum particle on
the atoms of a chygraph, hopping with amplitude ``hop`` between any two members
of a complex and with site energies ``eps``,

    H = sum_i eps_i |i><i| - hop sum_{a} sum_{i != j in a} |i><j|,

the object is the diagonal of the resolvent ``G(z) = (H - z)^{-1}`` at
``z = E + i0+``: its imaginary part over pi is the local density of states,
zero at almost every energy in the localised phase and finite in the extended
one.  On an incidence tree ``G_ii`` obeys the two steps of Sec. 8.2 with a
self-energy for a message:

    up:    G_ii^{-1} = eps_i - z - sum_{a ni i} Sigma_{a->i}            (18.9)
    down:  Sigma_{a->i} = hop^2 1^T M_a^{-1} 1,
           M_a = diag(eps_J - z - Sigma_{J->a}) - hop A_J                (18.5)

where ``J = a \\ i`` are the other members of the complex, ``A_J`` the adjacency
of the clique among them and ``Sigma_{J->a}`` the diagonal of their cavity
self-energies from their other complexes.  The down step is a Schur
complement: every path from ``i`` into the complex and back, the members' own
branches folded into ``Sigma_{J->a}``.  At cardinality two it is the tree
recursion of Abou-Chacra, Anderson and Thouless; at three the closed ``2 x 2``
form of Eq. (18.6).

**The clique interior is linear in the cardinality.**  ``A_J = 1 1^T - I``, so
``M_a = diag(m_J + hop) - hop 1 1^T`` is a rank-one update of a diagonal and
Sherman--Morrison gives, with ``d_j = 1 / (m_j + hop)`` and ``S = sum_j d_j``,

    v = M_a^{-1} 1 = d / (1 - hop S),
    Sigma_{a->i} = hop^2 S / (1 - hop S),
    eta_{a->i}   = hop 1^T M_a^{-1} eta_J = hop (d . eta_J) / (1 - hop S).

No matrix is inverted, at any cardinality.  The pole ``1 = hop S`` is the
resonance of the complex.  :func:`transmit` is that formula and is the only
place the interior is computed; :func:`cavity_instance` runs the recursion on
an explicit complex list, :func:`landscape_instance` carries the Localization
Landscape ``u = (H - E_min)^{-1} 1`` by the same linearity (Eq. 18.7), and
:func:`linearised_imaginary` is the linear map of Eq. (18.10) whose growth
decides the mobility edge.  The population dynamics that runs these on a
random ensemble lives with the scans in ``statmech/probe/anderson.py``.

Reference: Tonetti, Cugliandolo and Tarzia, Phys. Rev. Lett. 137, 136302
(2026), for the tree; Abou-Chacra, Thouless and Anderson, J. Phys. C 6, 1734
(1973), for the criterion.
"""

from collections import namedtuple

import numpy as np

Transmitted = namedtuple('Transmitted', 'sigma v eta')
Transmitted.__doc__ = """What a complex transmits to one member.

``sigma`` is the self-energy ``Sigma_{a->i}``; ``v = M_a^{-1} 1`` carries the
weights ``hop^2 v_j^2`` of the linearised imaginary map; ``eta`` is the landscape
row sum ``eta_{a->i}`` when ``eta_J`` was given, else ``None``."""


# --------------------------------------------------------------- the down step

def transmit(m_J, hop=1.0, eta_J=None):
    """The intra-complex step, Eq. (18.5), for a clique of any cardinality.

    Args:
        m_J: ``(..., m)`` array of ``m_j = eps_j - z - Sigma_{j->a}`` over the
            other members ``J`` of the complex; real or complex.
        hop: the hopping amplitude.
        eta_J: optional ``(..., m)`` array of the members' cavity row sums, for
            the landscape.

    Returns:
        :class:`Transmitted` with ``sigma`` of shape ``(...)``, ``v`` of shape
        ``(..., m)`` and ``eta`` of shape ``(...)`` or ``None``.

    The formula is the Sherman--Morrison closed form of the module docstring,
    which equals ``hop^2 1^T M_a^{-1} 1`` with ``M_a = diag(m_J) - hop A_J``
    (checked against the explicit inverse in ``test_resolvent``).  With one
    other member it is ``hop^2 / m_1``, the tree recursion; with two it is
    ``hop^2 (m_1 + m_2 + 2 hop) / (m_1 m_2 - hop^2)``, Eq. (18.6).
    """
    m_J = np.asarray(m_J)
    d = 1.0 / (m_J + hop)
    S = d.sum(axis=-1)
    denom = 1.0 - hop * S
    v = d / denom[..., None]
    sigma = hop * hop * S / denom
    eta = None
    if eta_J is not None:
        eta = hop * (d * np.asarray(eta_J)).sum(axis=-1) / denom
    return Transmitted(sigma, v, eta)


def linearised_imaginary(v, im_J, hop=1.0):
    """Eq. (18.10): ``Im Sigma_{a->i} = hop^2 sum_j v_j^2 Im Sigma_{j->a}``.

    The imaginary parts to first order, with the real parts frozen: a linear
    map with positive random weights ``hop^2 v_j^2``, whose growth rate per
    generation of a fractional moment decides the mobility edge.
    """
    return hop * hop * (np.asarray(v) ** 2 * np.asarray(im_J)).sum(axis=-1)


# --------------------------------------------------------------- instances

def hamiltonian(n, complexes, eps, hop=1.0):
    """The Anderson Hamiltonian on an explicit complex list, as a dense matrix."""
    H = np.diag(np.asarray(eps, float)).astype(complex)
    for c in complexes:
        for a in c:
            for b in c:
                if a != b:
                    H[a, b] -= hop
    return H


def _incidence(n, complexes):
    inc = [[] for _ in range(n)]
    for ci, c in enumerate(complexes):
        for p, a in enumerate(c):
            inc[a].append((ci, p))
    return inc


def _sum_messages(n, complexes, msg, dtype, base=0.0):
    tot = np.full(n, base, dtype)
    for (ci, p), val in msg.items():
        tot[complexes[ci][p]] += val
    return tot


def cavity_instance(n, complexes, eps, z, hop=1.0, sweeps=200, damping=0.0):
    """The two steps on one instance: ``G_ii`` for every atom.

    Messages ``Sigma_{a->i}`` are indexed by ``(complex index, position of i in
    the complex)`` and updated in parallel.  Exact on an incidence tree once
    ``sweeps`` exceeds its depth; on a graph with long loops the error decays
    with ``Im z``.  Returns ``(G, messages, totals)`` with ``totals[i]`` the
    full self-energy ``sum_a Sigma_{a->i}``.
    """
    eps = np.asarray(eps, float)
    sig = {(ci, p): 0j for ci, c in enumerate(complexes) for p in range(len(c))}
    for _ in range(sweeps):
        tot = _sum_messages(n, complexes, sig, complex)
        new = {}
        for ci, c in enumerate(complexes):
            for p, i in enumerate(c):
                J = [(q, j) for q, j in enumerate(c) if q != p]
                m_J = np.array([eps[j] - z - (tot[j] - sig[(ci, q)]) for q, j in J])
                new[(ci, p)] = transmit(m_J, hop).sigma
        if damping:
            sig = {k: damping * sig[k] + (1 - damping) * new[k] for k in sig}
        else:
            sig = new
    tot = _sum_messages(n, complexes, sig, complex)
    G = 1.0 / (eps - z - tot)
    return G, sig, tot


def landscape_instance(n, complexes, eps, emin, hop=1.0, sweeps=200):
    """The Localization Landscape ``u = (H - E_min)^{-1} 1`` by cavity, Eq. (18.7).

    ``eta_{a->i} = hop 1^T M_a^{-1} eta_{J->a}``, ``eta_i = 1 + sum_a eta_{a->i}``,
    ``u_i = G_ii(E_min) eta_i``, in real arithmetic: at ``E_min`` below the
    spectrum nothing needs regularising.
    """
    eps = np.asarray(eps, float)
    G, sig, tot = cavity_instance(n, complexes, eps, emin + 0j, hop, sweeps)
    eta = {k: 0.0 for k in sig}
    for _ in range(sweeps):
        etot = _sum_messages(n, complexes, eta, float, base=1.0)
        new = {}
        for ci, c in enumerate(complexes):
            for p, i in enumerate(c):
                J = [(q, j) for q, j in enumerate(c) if q != p]
                m_J = np.array([(eps[j] - emin - (tot[j] - sig[(ci, q)])).real for q, j in J])
                eta_J = np.array([etot[j] - eta[(ci, q)] for q, j in J])
                new[(ci, p)] = transmit(m_J, hop, eta_J).eta
        eta = new
    etot = _sum_messages(n, complexes, eta, float, base=1.0)
    return G.real * etot


def regular_instance(n, cardinalities, degrees, rng):
    """A random regular chygraph: ``n`` atoms, each in ``degrees[l]`` complexes
    of cardinality ``cardinalities[l]``.

    Configuration model on stubs, one layer at a time.  Degenerate complexes
    (a repeated atom) are redrawn a few times and then dropped, which costs a
    vanishing fraction at large ``n``.  Returns the list of complexes as
    tuples of atoms.
    """
    complexes = []
    for c, k in zip(cardinalities, degrees):
        c, k = int(c), int(k)
        if k == 0:
            continue
        stubs = np.repeat(np.arange(n), k)
        for _ in range(20):
            rng.shuffle(stubs)
            groups = stubs[:(len(stubs) // c) * c].reshape(-1, c)
            bad = np.zeros(len(groups), bool)
            for i in range(c):
                for j in range(i + 1, c):
                    bad |= groups[:, i] == groups[:, j]
            if not bad.any():
                break
        complexes += [tuple(int(x) for x in g) for g in groups if len(set(g)) == c]
    return complexes


def incidence_tree(cardinalities, degrees, depth):
    """The incidence tree of a regular chygraph to a given depth: the root in
    ``degrees[l]`` complexes of each layer, every other atom in ``degrees[l]-1``
    beyond the one it was reached through.  Returns ``(n, complexes)``.  This is
    the structure on which the recursion is exact.
    """
    complexes, n, frontier = [], 1, [(0, 0)]
    while frontier:
        a, g = frontier.pop()
        if g >= depth:
            continue
        for c, k in zip(cardinalities, degrees):
            for _ in range(int(k) if g == 0 else int(k) - 1):
                members = tuple(range(n, n + int(c) - 1))
                complexes.append((a,) + members)
                frontier += [(b, g + 1) for b in members]
                n += int(c) - 1
    return n, complexes


# --------------------------------------------------------------- the clean ensemble

def uniform_fixed_point(cardinalities, degrees, E, hop=1.0, iters=5000, tol=1e-13):
    """The real fixed point of the recursion at ``W = 0`` and energy ``E`` on a
    regular chygraph, one message per layer, or ``None`` if there is none.

    Without disorder every message on a layer is one number: ``sigma_l`` with
    ``m_l = -E - sum_m (k_m - delta_lm) sigma_m`` and the down step of
    :func:`transmit` over ``c_l - 1`` copies of ``m_l``.  A real solution
    exists below the band and not inside it, which is what
    :func:`band_bottom` bisects on; the imaginary part it acquires inside the
    band is Eq. (18.3) on the tree.
    """
    c = np.asarray(cardinalities, int)
    k = np.asarray(degrees, float)
    L = len(c)
    sig = np.zeros(L)
    for _ in range(iters):
        new = np.empty(L)
        for l in range(L):
            cav = float(((k - np.eye(L)[l]) * sig).sum())
            m = -E - cav
            d = 1.0 / (m + hop)
            S = (c[l] - 1) * d
            if not (m + hop > 0 and 1.0 - hop * S > 0):
                return None
            new[l] = hop * hop * S / (1.0 - hop * S)
        if not np.all(np.isfinite(new)):
            return None
        if np.abs(new - sig).sum() < tol:
            return new
        sig = 0.5 * (sig + new)
    return None


def band_bottom(cardinalities, degrees, hop=1.0, lo=-12.0, hi=0.0, iters=60):
    """Bottom ``e_min`` of the bulk spectrum of the pure hopping problem: the
    largest ``E`` below which :func:`uniform_fixed_point` exists.

    ``-2 sqrt(K) hop`` on the tree of degree ``K + 1``; lower with triangles,
    since a triangle's lowest eigenvalue is ``-2 hop`` against an edge's
    ``-hop``.  The isolated eigenvalue of the uniform vector, at ``-hop`` times
    the degree, is not part of the bulk and is not found here.
    """
    for _ in range(iters):
        mid = 0.5 * (lo + hi)
        if uniform_fixed_point(cardinalities, degrees, mid, hop) is not None:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def spectral_bottom(cardinalities, degrees, W, hop=1.0):
    """``E_min(W)``: the band bottom minus ``W/2``, or the isolated eigenvalue
    ``-hop * degree`` when that lies below it (``W < W_min``)."""
    degree = float(sum(k * (c - 1) for c, k in zip(cardinalities, degrees)))
    return min(band_bottom(cardinalities, degrees, hop) - W / 2.0, -hop * degree)


def clean_tree(K, z, hop=1.0):
    """The clean tree of degree ``K + 1`` in closed form, Eqs. (18.3)-(18.4).

    Returns ``(G_cavity, G_ii)`` on the branch with ``Im G > 0`` for
    ``Im z > 0``: ``G_cavity`` is the root of ``G^{-1} = -z - K hop^2 G`` and
    ``G_ii^{-1} = -z - (K + 1) hop^2 G_cavity``.  ``Im G_ii / pi`` is the
    Kesten--McKay density of states on ``|E| <= 2 sqrt(K) hop``.
    """
    z = complex(z)
    r = np.sqrt(z * z - 4 * K * hop * hop)
    roots = ((-z + r) / (2 * K * hop * hop), (-z - r) / (2 * K * hop * hop))
    if abs(z.imag) > 1e-14:
        g = roots[0] if roots[0].imag > 0 else roots[1]
    else:
        # real z outside the band: both roots real, the physical one is the
        # branch continuous with G ~ -1/z at infinity, the smaller in modulus
        g = min(roots, key=abs)
    return g, 1.0 / (-z - (K + 1) * hop * hop * g)


def kesten_mckay(K, E, hop=1.0):
    """The density of states of the clean tree of degree ``K + 1``, Eq. (18.4)."""
    E = np.asarray(E, float)
    inside = np.abs(E) < 2 * np.sqrt(K) * hop
    rho = np.zeros_like(E)
    Ei = E[inside]
    rho[inside] = ((K + 1) / (2 * np.pi) * np.sqrt(4 * K * hop * hop - Ei * Ei)
                   / ((K + 1) ** 2 * hop * hop - Ei * Ei))
    return rho


# --------------------------------------------------------------- the incidence tree's branching

def incidence_branching(cardinalities, degrees):
    """Eq. (18.13): how the numbers of messages per layer grow per generation
    of the incidence tree of a regular chygraph.

    An atom reached through a layer-``l`` complex passes a message on through
    ``k_m - delta_lm`` complexes of layer ``m``, each of which reaches
    ``c_m - 1`` atoms.  The leading eigenvalue is the branching the imaginary
    parts propagate on, and it orders the band-centre critical disorders: a
    triangle spends two of a site's neighbours on each other, so at degree
    three it cuts the branching from 2 to ``sqrt 2``.
    """
    c = np.asarray(cardinalities, float)
    k = np.asarray(degrees, float)
    L = len(c)
    B = np.empty((L, L))
    for l in range(L):
        for m in range(L):
            B[l, m] = (k[m] - (1.0 if m == l else 0.0)) * (c[m] - 1.0)
    return B
