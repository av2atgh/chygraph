"""Examples: percolation in graphs with links and triangles.

The graph is mapped to a multiplex hypergraph with two layers:
  - l=1: links
  - l=2: triangles

Covers both general and Poisson degree distributions, plus numerical evaluation.
"""

import matplotlib.pyplot as plt
from sympy import symbols
from chygraph import GraphWithTriangles

plt.rcParams.update({'font.size': 14})


def make_theta_evaluator(theta_expr):
    """Return a callable that substitutes (kL, kT, q) into a theta expression."""
    def ftheta(kL, kT, q_val):
        return theta_expr.subs([
            (symbols('k_L'), kL), (symbols('k_T'), kT),
            (symbols('K_L'), kL), (symbols('K_T'), kT),
            (symbols('q'),   q_val),
        ])
    return ftheta


def make_lambda_evaluator(lambda_expr):
    """Return a callable that substitutes (kL, kT, q) into a Lambda expression."""
    def fLambda(kL, kT, q_val):
        return lambda_expr.subs([
            (symbols('k_L'), kL), (symbols('k_T'), kT),
            (symbols('K_L'), kL), (symbols('K_T'), kT),
            (symbols('q'),   q_val),
        ])
    return fLambda


def plot_theta_lambda(ftheta, fLambda, kL=1, kT=0.5, n_points=100):
    """Plot theta and Lambda as functions of q.

    Args:
        ftheta:   Callable(kL, kT, q) returning theta.
        fLambda:  Callable(kL, kT, q) returning Lambda.
        kL:       Mean link degree.
        kT:       Mean triangle degree.
        n_points: Number of q samples in [0, 1].

    Returns:
        (fig, ax) matplotlib objects.
    """
    x = [i / n_points for i in range(0, n_points + 1)]

    fig, ax = plt.subplots(2, 1)

    label1 = r" $\langle k\rangle_|=\langle \bar{k}\rangle_|=1,\ \langle k\rangle_\Delta=\langle \bar{k}\rangle_\Delta=0.5$"
    label2 = r" $\langle k\rangle_|=\langle \bar{k}\rangle_|=2,\ \langle k\rangle_\Delta=\langle \bar{k}\rangle_\Delta=0$"

    y1 = [ftheta(kL, kT, xi) for xi in x]
    y2 = [ftheta(kL + 2*kT, 0, xi) for xi in x]
    ax[0].plot(x, y1, c="orange", ls="-",  label=label1)
    ax[0].plot(x, y2, c="green",  ls="-.", label=label2)
    ax[0].plot([min(x), max(x)], [0, 0], "k--")
    ax[0].set(xlabel=r"$q$", ylabel=r"$\theta$")
    ax[0].legend(loc='upper left', frameon=0, handletextpad=0, fontsize=14)
    ax[0].set_title("a)", x=-0.18, y=1)

    z1 = [fLambda(kL, kT, xi) for xi in x]
    z2 = [fLambda(kL + 2*kT, 0, xi) for xi in x]
    ax[1].plot(x, z1, c="orange", ls="-",  label=label1)
    ax[1].plot(x, z2, c="green",  ls="-.", label=label2)
    ax[1].plot([min(x), max(x)], [0, 0], "k--")
    ax[1].set(xlabel=r"$q$", ylabel=r"$\Lambda$")
    ax[1].legend(loc='lower right', frameon=0, handletextpad=0, fontsize=14)
    ax[1].set_title("b)", x=-0.18, y=1)

    plt.subplots_adjust(bottom=0, top=2, hspace=0.25)
    return fig, ax


def example_general_degree():
    """Any degree distribution."""
    NM = GraphWithTriangles(poisson=False)
    print("theta:", NM.A.theta())
    print("eigenvalues:", NM.A.eigenvals())
    print("Lambda:", NM.A.eigenvals()[3])
    return NM


def example_poisson_degree():
    """Poisson degree distribution."""
    NMP = GraphWithTriangles(poisson=True)
    print("theta:", NMP.A.theta())
    print("eigenvalues:", NMP.A.eigenvals())
    print("Lambda:", NMP.A.eigenvals()[2])


def example_numerical_evaluation(NM):
    """Numerical evaluation and plotting for the general-degree model."""
    theta_expr  = NM.A.theta()
    lambda_expr = NM.A.eigenvals()[3]
    ftheta  = make_theta_evaluator(theta_expr)
    fLambda = make_lambda_evaluator(lambda_expr)
    fig, ax = plot_theta_lambda(ftheta, fLambda, kL=1, kT=0.5)
    plt.show()


if __name__ == "__main__":
    NM = example_general_degree()
    example_poisson_degree()
    example_numerical_evaluation(NM)
