"""Chapter 45, Section 2: linear quantile regression on Engel's budget survey.

Solves the check-loss problem as a linear program from first principles (scipy's
HiGHS solver), checks it against statsmodels' implementation, verifies the counting
identity for the signs of the residuals and the equivariance properties, and draws
the fitted quantile lines.
"""
import matplotlib.pyplot as plt
import numpy as np
import statsmodels.api as sm
import statsmodels.formula.api as smf
from scipy.optimize import linprog

from regbook import COLORS, Generated, figure_path, use_book_style

data = sm.datasets.engel.load_pandas().data
TAUS = (0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95)

# <<lp>>
income = data["income"].to_numpy()
food = data["foodexp"].to_numpy()
n = len(food)
X = np.column_stack([np.ones(n), income])


def check_loss(u, tau):
    """The check function rho_tau applied elementwise."""
    return u * (tau - (u < 0))


def qreg(X, y, tau):
    """Quantile regression by linear programming.

    Variables are (b+, b-, u, v), all nonnegative, with beta = b+ - b- and
    u - v = y - X beta the positive and negative parts of the residual.
    """
    n, p = X.shape
    cost = np.concatenate([np.zeros(2 * p), tau * np.ones(n), (1 - tau) * np.ones(n)])
    A = np.hstack([X, -X, np.eye(n), -np.eye(n)])
    out = linprog(cost, A_eq=A, b_eq=y, bounds=(0, None), method="highs")
    return out.x[:p] - out.x[p:2 * p]


beta = {tau: qreg(X, food, tau) for tau in TAUS}
for tau in (0.10, 0.50, 0.90):
    print(f"tau = {tau:.2f}:  intercept {beta[tau][0]:8.2f}   slope {beta[tau][1]:.4f}")
# <</lp>>

# the same fits from statsmodels, to three decimals
for tau in TAUS:
    ref = smf.quantreg("foodexp ~ income", data).fit(q=tau).params.to_numpy()
    assert np.allclose(ref, beta[tau], atol=1e-3), (tau, ref, beta[tau])

# <<counts>>
for tau in TAUS:
    r = food - X @ beta[tau]
    n_zero = np.sum(np.abs(r) < 1e-8)
    n_neg = np.sum(r < -1e-8)
    assert n_zero == X.shape[1]                      # p residuals vanish exactly
    assert n_neg <= tau * n <= n_neg + n_zero        # the counting identity
print("at every tau: 2 residuals are zero and #negative <= tau*n <= #negative + 2")
# <</counts>>

# equivariance: scale, shift, and reparameterization
a, gamma = 2.5, np.array([30.0, -0.2])
b25, b75 = qreg(X, food, 0.25), qreg(X, food, 0.75)
assert np.allclose(qreg(X, a * food, 0.25), a * b25, atol=1e-6)
assert np.allclose(qreg(X, -food, 0.25), -b75, atol=1e-6)
assert np.allclose(qreg(X, food + X @ gamma, 0.25), b25 + gamma, atol=1e-6)
A = np.array([[1.0, 500.0], [0.0, 1000.0]])          # income recentred and rescaled
assert np.allclose(A @ qreg(X @ A, food, 0.25), b25, atol=1e-6)

# Koenker and Machado's R1(tau): the check loss relative to the intercept-only fit
r1 = {}
for tau in TAUS:
    full = check_loss(food - X @ beta[tau], tau).sum()
    null = check_loss(food - np.quantile(food, tau), tau).sum()
    r1[tau] = 1 - full / null

# the least squares slope, for comparison
ols = np.linalg.lstsq(X, food, rcond=None)[0]
median_share = beta[0.50] @ [1, 1000.0] / 1000.0     # median food share at 1000 francs
low_share = beta[0.10] @ [1, 1000.0] / 1000.0

gen = Generated("ch45", "quantile_engel")
gen.int("n", n)
for tau in TAUS:
    key = f"{round(100 * tau):02d}"
    gen.num(f"slope{key}", beta[tau][1], 4)
    gen.num(f"int{key}", beta[tau][0], 1)
    gen.num(f"r1_{key}", r1[tau], 3)
gen.num("ols_slope", ols[1], 4)
gen.num("slope_ratio", beta[0.95][1] / beta[0.05][1], 2)
gen.num("median_share", median_share, 3)
gen.num("low_share", low_share, 3)
gen.num("fit_at_500_10", beta[0.10] @ [1, 500.0], 1)
gen.num("fit_at_500_90", beta[0.90] @ [1, 500.0], 1)
gen.num("fit_at_2000_10", beta[0.10] @ [1, 2000.0], 1)
gen.num("fit_at_2000_90", beta[0.90] @ [1, 2000.0], 1)
gen.write()

# ---- figure ----------------------------------------------------------------
use_book_style()
fig, axes = plt.subplots(1, 2, figsize=(5.6, 2.6))

ax = axes[0]
ax.scatter(income, food, s=8, color=COLORS["muted"], alpha=0.6, linewidths=0)
grid = np.linspace(income.min(), 3000, 2)
shades = plt.cm.viridis(np.linspace(0.12, 0.88, len(TAUS)))
for tau, colour in zip(TAUS, shades):
    ax.plot(grid, beta[tau][0] + beta[tau][1] * grid, color=colour, linewidth=1.0)
    ax.annotate(f"{tau:g}", (3000, beta[tau][0] + beta[tau][1] * 3000),
                color=colour, fontsize=6.5, xytext=(2, -2), textcoords="offset points")
ax.plot(grid, ols[0] + ols[1] * grid, color=COLORS["second"], linewidth=1.2, linestyle="--")
ax.set_xlim(0, 3400)
ax.set_ylim(0, 2200)
ax.set_xlabel("household income")
ax.set_ylabel("food expenditure")
ax.set_title("(a) seven fitted quantile lines")

ax = axes[1]
taus_fine = np.linspace(0.05, 0.95, 19)
slopes = [qreg(X, food, t)[1] for t in taus_fine]
ax.plot(taus_fine, slopes, color=COLORS["accent"], marker="o", markersize=2.5)
ax.axhline(ols[1], color=COLORS["second"], linestyle="--", linewidth=1.0)
ax.annotate("least squares", (0.08, ols[1]), color=COLORS["second"], fontsize=7,
            xytext=(0, 4), textcoords="offset points")
ax.set_xlabel(r"$\tau$")
ax.set_ylabel("slope on income")
ax.set_title("(b) the slope depends on the quantile")
fig.tight_layout()
fig.savefig(figure_path("ch45", "engel_quantiles"))
