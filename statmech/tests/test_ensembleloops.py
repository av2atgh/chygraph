"""Tests for the loop series on the ensemble (statmech.ensembleloops)."""

import math

import numpy as np

from statmech.bethehessian import clique_factor_graph, random_chygraph
from statmech.ensembleloops import (
    branching_power, cycle_sum, ensemble_series, exact_minus_bethe,
)
from statmech.ising import branching_matrix, clique_derivative
from statmech.loopseries import loop_series


def test_cycle_term_is_product_of_uprimes_and_exact_on_a_ring():
    # a ring of three triangles, each pair sharing one atom: one cycle
    cx = [(0, 1, 2), (2, 3, 4), (4, 5, 0)]
    bj = 0.4
    fg = clique_factor_graph(cx, bj)
    ex = exact_minus_bethe(fg)
    u = clique_derivative(3, bj)
    assert abs(ex - math.log1p(u ** 3)) < 1e-12
    assert abs(cycle_sum(cx, bj)[0] - ex) < 1e-12


def test_full_series_equals_exact_and_cycles_when_disjoint():
    rng = np.random.default_rng(1)
    seen = 0
    for _ in range(12):
        cx = random_chygraph(14, [2, 3], [1.2, 0.5], rng)
        if not 3 <= len(cx) <= 12:
            continue
        fg = clique_factor_graph(cx, 0.25)
        ex = exact_minus_bethe(fg)
        fg.bp()
        full = math.log1p(sum(r for _, _, r in loop_series(fg)))
        assert abs(full - ex) < 1e-10
        seen += 1
    assert seen >= 3


def test_closed_form_is_the_limit_of_the_truncated_sum():
    cards, means = [2, 3], [1.2, 0.5]
    clos, first = ensemble_series(cards, means, 0.3)
    prev = None
    for lmax in (8, 16, 32):
        tr, _ = ensemble_series(cards, means, 0.3, lmax=lmax)
        if prev is not None:
            assert abs(tr - clos) < abs(prev - clos)
        prev = tr
    assert abs(prev - clos) < 1e-3
    # single layer: (lambda^l / 2l) ln(1 + u'^l) summed
    lam = 2 * 0.8
    u = clique_derivative(3, 0.3)
    hand = sum(lam ** l / (2 * l) * math.log1p(u ** l) for l in range(2, 400))
    assert abs(ensemble_series([3], [0.8], 0.3)[0] - hand) < 1e-10
    B2 = branching_power([3], [0.8], 0.3, 2)
    assert abs(B2[0, 0] - lam * u ** 2) < 1e-14
    assert np.allclose(branching_power(cards, means, 0.3, 1),
                       branching_matrix(cards, means, 0.3))


def test_series_diverges_at_the_branching_threshold():
    from statmech.ising import critical_coupling
    bc = critical_coupling([2, 3], [1.2, 0.5])
    below = ensemble_series([2, 3], [1.2, 0.5], 0.999 * bc)[0]
    assert math.isfinite(below) and below > 1.0
    assert ensemble_series([2, 3], [1.2, 0.5], 1.001 * bc)[0] == math.inf
