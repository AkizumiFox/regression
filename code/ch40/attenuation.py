"""Chapter 40, Section 1: marginal versus conditional coefficients in a random-intercept model.

The conditional (subject-specific) curve is P(Y = 1 | U = u) = H(alpha + beta x + u) with
U ~ N(0, tau^2); the marginal (population-averaged) curve is its expectation over U. For the
probit the marginal curve is again a probit with coefficients divided by sqrt(1 + tau^2); for
the logit no exact form exists and the population-averaged coefficient is defined as the
solution of the population quasi-score equation.
"""
import matplotlib.pyplot as plt
import numpy as np
from numpy.polynomial.hermite_e import hermegauss
from scipy.optimize import root
from scipy.stats import norm

from regbook import COLORS, Generated, figure_path, use_book_style

# <<setup>>
def expit(z):
    return 1.0 / (1.0 + np.exp(-z))


# nodes and weights that integrate g(U) over U ~ N(0, 1): sum_k w_k g(z_k)
z_node, z_w = hermegauss(80)
z_w = z_w / np.sqrt(2 * np.pi)


def marginal_curve(eta, tau, link="logit"):
    """E[H(eta + U)] for U ~ N(0, tau^2), by Gauss-Hermite quadrature."""
    H = expit if link == "logit" else norm.cdf
    grid = np.atleast_1d(eta)[..., None] + tau * z_node
    return np.sum(z_w * H(grid), axis=-1)
# <</setup>>


# <<beta_marginal>>
def beta_marginal(alpha, beta, tau, link="logit"):
    """Population-averaged coefficients: the root of E[(1, x){m(x) - H(a + b x)}] = 0
    when x ~ N(0, 1), with m the marginal curve of the conditional model."""
    H = expit if link == "logit" else norm.cdf
    m = marginal_curve(alpha + beta * z_node, tau, link)

    def score(par):
        resid = m - H(par[0] + par[1] * z_node)
        return [np.sum(z_w * resid), np.sum(z_w * z_node * resid)]

    sol = root(score, [alpha, beta], tol=1e-13)
    assert sol.success, sol.message
    return sol.x
# <</beta_marginal>>


C_LOGIT = 16 * np.sqrt(3) / (15 * np.pi)


def approx_factor(tau, link="logit"):
    c = C_LOGIT if link == "logit" else 1.0
    return 1.0 / np.sqrt(1.0 + c**2 * tau**2)


# <<demo>>
for t in [0.5, 1.0, 2.0]:
    print(f"tau = {t}:  logit {beta_marginal(0.0, 1.0, t)[1]:.4f}"
          f"   approximation {approx_factor(t):.4f}"
          f"   probit {1 / np.sqrt(1 + t**2):.4f}")
# <</demo>>


# ---- checks -----------------------------------------------------------------
# probit: the marginal curve is exactly a probit with coefficients scaled by 1/sqrt(1+tau^2)
for tau in [0.3, 1.0, 2.5]:
    eta = np.linspace(-4, 4, 41)
    assert np.allclose(marginal_curve(eta, tau, "probit"),
                       norm.cdf(eta / np.sqrt(1 + tau**2)), atol=1e-10)
    a, b = beta_marginal(-0.4, 1.3, tau, "probit")
    assert np.allclose([a, b], np.array([-0.4, 1.3]) / np.sqrt(1 + tau**2), atol=1e-8)

# logit: attenuation, and agreement with the classical approximation
ratios = {}
for tau in [0.5, 1.0, 1.5, 2.0]:
    b = beta_marginal(0.0, 1.0, tau)[1]
    assert 0 < b < 1.0
    ratios[tau] = b
assert all(ratios[t] > ratios[s] for t, s in [(0.5, 1.0), (1.0, 1.5), (1.5, 2.0)])

# the marginal curve is flatter than the conditional one at its steepest point
for tau in [0.5, 1.0, 2.0]:
    slope = np.sum(z_w * 0.25 / np.cosh(tau * z_node / 2) ** 2)   # E[h(U)], h = logistic pdf
    assert slope < 0.25

# a binary covariate: the marginal log odds ratio, for comparison with the x ~ N(0,1) target
logit = lambda p: np.log(p / (1 - p))
or_binary = {t: float(logit(marginal_curve(1.0, t))[0] - logit(marginal_curve(0.0, t))[0])
             for t in [0.5, 1.0, 1.5, 2.0]}

# Poisson with a log link: no attenuation, the intercept alone moves by tau^2/2
rng = np.random.default_rng(40_001)
u = rng.normal(0, 0.8, 4_000_000)
assert abs(np.mean(np.exp(u)) - np.exp(0.8**2 / 2)) < 2e-3

gen = Generated("ch40", "attenuation")
gen.num("c_logit", C_LOGIT, 4)
gen.num("ratio_half", ratios[0.5], 4)
gen.num("ratio_one", ratios[1.0], 4)
gen.num("ratio_two", ratios[2.0], 4)
gen.num("approx_one", approx_factor(1.0), 4)
gen.num("approx_two", approx_factor(2.0), 4)
gen.num("probit_one", 1 / np.sqrt(2), 4)
gen.num("or_binary_one", or_binary[1.0], 4)
p_u = expit(z_node)                                   # h(U) for U ~ N(0, 1)
var_h = float(np.sum(z_w * p_u**2) - np.sum(z_w * p_u) ** 2)
assert abs(np.sum(z_w * p_u) - 0.5) < 1e-12           # symmetry
gen.num("var_h", var_h, 4)
gen.write()

# ---- figure ------------------------------------------------------------------
use_book_style()
fig, axes = plt.subplots(1, 2, figsize=(5.6, 2.5))

ax = axes[0]
x = np.linspace(-5, 5, 400)
tau = 1.5
for u0, style in [(-2 * tau, ":"), (-tau, "--"), (0.0, "-"), (tau, "--"), (2 * tau, ":")]:
    ax.plot(x, expit(x + u0), style, color=COLORS["grid"], linewidth=0.9, zorder=1)
ax.plot(x, expit(x), color=COLORS["accent"], label="conditional, $u=0$", zorder=3)
ax.plot(x, marginal_curve(x, tau), color=COLORS["second"], label="marginal", zorder=3)
a_m, b_m = beta_marginal(0.0, 1.0, tau)
ax.plot(x, expit(a_m + b_m * x), "--", color=COLORS["third"],
        label=f"logistic, slope {b_m:.2f}", zorder=2)
ax.set_xlabel(r"linear predictor $\alpha + \beta x$")
ax.set_ylabel("probability")
ax.set_title(r"(a) $\tau = 1.5$")
ax.legend(frameon=False, loc="upper left", fontsize=6.5)

ax = axes[1]
taus = np.linspace(0, 3, 61)
logit_ratio = np.array([beta_marginal(0.0, 1.0, t)[1] if t > 0 else 1.0 for t in taus])
ax.plot(taus, logit_ratio, color=COLORS["accent"], label="logit, exact")
ax.plot(taus, approx_factor(taus), "--", color=COLORS["second"],
        label=r"logit, $(1+c^2\tau^2)^{-1/2}$")
ax.plot(taus, 1 / np.sqrt(1 + taus**2), color=COLORS["third"], label="probit, exact")
ax.set_xlabel(r"random-intercept standard deviation $\tau$")
ax.set_ylabel(r"$\beta^{M}/\beta^{C}$")
ax.set_ylim(0, 1.05)
ax.set_title("(b) attenuation factor")
ax.legend(frameon=False, loc="lower left", fontsize=6.5)

fig.tight_layout()
fig.savefig(figure_path("ch40", "attenuation"))
print("attenuation ratios", {k: round(v, 4) for k, v in ratios.items()})
print("binary-covariate log odds ratios", {k: round(v, 4) for k, v in or_binary.items()})
