"""Chapter 45, Section 3: standard errors, intervals and crossing for quantile regression.

Three routes to a standard error are compared on Engel's data: the iid-error formula
with an estimated sparsity function, Powell's sandwich for independent but not
identically distributed errors, and the pairs bootstrap of chapter 23. The crossing
of fitted quantile curves is located exactly, and rearrangement repairs it.
"""
import matplotlib.pyplot as plt
import numpy as np
import statsmodels.api as sm
from scipy import stats
from scipy.optimize import linprog

from regbook import COLORS, Generated, figure_path, use_book_style

data = sm.datasets.engel.load_pandas().data
TAUS = (0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95)


def qreg(X, y, tau):
    """Quantile regression by linear programming (as in Section 45.2)."""
    n, p = X.shape
    cost = np.concatenate([np.zeros(2 * p), tau * np.ones(n), (1 - tau) * np.ones(n)])
    A = np.hstack([X, -X, np.eye(n), -np.eye(n)])
    out = linprog(cost, A_eq=A, b_eq=y, bounds=(0, None), method="highs")
    return out.x[:p] - out.x[p:2 * p]


# <<sparsity>>
income = data["income"].to_numpy()
food = data["foodexp"].to_numpy()
n = len(food)
X = np.column_stack([np.ones(n), income])


def hall_sheather(n, tau, alpha=0.05):
    """Bandwidth in the tau direction for a difference quotient of the quantile function."""
    z, za = stats.norm.ppf(tau), stats.norm.ppf(1 - alpha / 2)
    return n ** (-1 / 3) * za ** (2 / 3) * (1.5 * stats.norm.pdf(z) ** 2 / (2 * z ** 2 + 1)) ** (1 / 3)


def sparsity(resid, tau, n):
    """Siddiqui's estimate of s(tau) = 1 / f(F^{-1}(tau)): a difference quotient."""
    h = min(hall_sheather(n, tau), tau, 1 - tau)
    hi, lo = np.quantile(resid, [min(tau + h, 1.0), max(tau - h, 0.0)])
    return (hi - lo) / (2 * h)


def se_iid(X, y, tau):
    """Standard errors under independent and identically distributed errors."""
    beta = qreg(X, y, tau)
    s = sparsity(y - X @ beta, tau, len(y))
    return beta, np.sqrt(tau * (1 - tau) * s ** 2 * np.diag(np.linalg.inv(X.T @ X)))


def se_sandwich(X, y, tau):
    """Powell's sandwich: a kernel estimate of the matrix sum_i f_i(0) x_i x_i'."""
    beta = qreg(X, y, tau)
    n = len(y)
    r = y - X @ beta
    h = sparsity(r, tau, n) * hall_sheather(n, tau)    # a window on the residual scale
    dens = np.exp(-0.5 * (r / h) ** 2) / (h * np.sqrt(2 * np.pi))
    D1 = np.linalg.inv(X.T @ (dens[:, None] * X))
    return beta, np.sqrt(tau * (1 - tau) * np.diag(D1 @ (X.T @ X) @ D1))


def se_bootstrap(X, y, tau, B, rng):
    """The pairs bootstrap of chapter 23: resample (y_i, x_i) with replacement and refit."""
    draws = np.array([qreg(X[i], y[i], tau)
                      for i in rng.integers(0, len(y), (B, len(y)))])
    return draws.std(axis=0, ddof=1)


rng = np.random.default_rng(4545)
beta50, se50_iid = se_iid(X, food, 0.50)
_, se50_sand = se_sandwich(X, food, 0.50)
se50_small = se_bootstrap(X, food, 0.50, 199, np.random.default_rng(4545))
print(f"slope at tau = 0.5: {beta50[1]:.4f}")
print(f"  se: iid errors {se50_iid[1]:.4f}, sandwich {se50_sand[1]:.4f},"
      f" bootstrap with B = 199 {se50_small[1]:.4f}")
# <</sparsity>>

# the table in the text uses B = 999; the cell above uses 199 so that it runs quickly
se50_boot = se_bootstrap(X, food, 0.50, 999, rng)
print(f"  bootstrap se at B = 999: {se50_boot[1]:.4f}")

assert se50_sand[1] > 2 * se50_iid[1]           # the iid formula is far too optimistic here
assert abs(se50_boot[1] - se50_sand[1]) < 0.35 * se50_sand[1]

beta90, se90_iid = se_iid(X, food, 0.90)
_, se90_sand = se_sandwich(X, food, 0.90)
se90_boot = se_bootstrap(X, food, 0.90, 999, rng)

# coefficient path with bootstrap intervals, for the figure
taus_fine = np.round(np.linspace(0.05, 0.95, 19), 3)
path = np.array([qreg(X, food, t) for t in taus_fine])
path_se = np.array([se_bootstrap(X, food, t, 199, rng) for t in taus_fine])
ols = np.linalg.lstsq(X, food, rcond=None)[0]
ols_se = np.sqrt(np.sum((food - X @ ols) ** 2) / (n - 2) * np.diag(np.linalg.inv(X.T @ X)))

# ---- crossing --------------------------------------------------------------
# <<crossing>>
beta = {tau: qreg(X, food, tau) for tau in TAUS}
crossings = {}
for lo, hi in zip(TAUS[:-1], TAUS[1:]):
    d0, d1 = beta[hi] - beta[lo]
    crossings[(lo, hi)] = -d0 / d1                 # where two fitted lines meet
print(f"incomes in the data run from {income.min():.0f} to {income.max():.0f}")
print("neighbouring fitted lines meet at income",
      " ".join(f"{v:.0f}" for v in crossings.values()))
# <</crossing>>

assert max(crossings.values()) < income.min()      # every crossing is outside the data

# a quadratic fit crosses inside the range of the data, where the data thin out
Xq = np.column_stack([np.ones(n), income, income ** 2])
grid = np.linspace(income.min(), np.quantile(income, 0.99), 400)
Gq = np.column_stack([np.ones_like(grid), grid, grid ** 2])
curves = np.array([Gq @ qreg(Xq, food, t) for t in TAUS])
bad = np.any(np.diff(curves, axis=0) < 0, axis=0)
cross_frac = bad.mean()
first_cross = grid[np.argmax(bad)]
assert 0 < cross_frac < 0.5

# rearrangement repairs the ordering, and moves the curves by at most:
rearranged = np.sort(curves, axis=0)
assert np.all(np.diff(rearranged, axis=0) >= 0)
max_move = np.abs(rearranged - curves).max()

gen = Generated("ch45", "qr_inference")
gen.num("slope50", beta50[1], 4)
gen.num("se50_iid", se50_iid[1], 4)
gen.num("se50_sand", se50_sand[1], 4)
gen.num("se50_boot", se50_boot[1], 4)
gen.num("se_ratio50", se50_sand[1] / se50_iid[1], 1)
gen.num("slope90", beta90[1], 4)
gen.num("se90_iid", se90_iid[1], 4)
gen.num("se90_sand", se90_sand[1], 4)
gen.num("se90_boot", se90_boot[1], 4)
gen.num("ols_se", ols_se[1], 4)
gen.num("sparsity50", sparsity(food - X @ beta50, 0.50, n), 1)
gen.num("sparsity90", sparsity(food - X @ beta90, 0.90, n), 1)
gen.num("cross_max", max(crossings.values()), 0)
gen.num("income_min", income.min(), 0)
gen.num("cross_frac", cross_frac, 3)
gen.num("first_cross", first_cross, 0)
gen.num("max_move", max_move, 1)
gen.write()

# ---- figure ----------------------------------------------------------------
use_book_style()
fig, axes = plt.subplots(1, 2, figsize=(5.6, 2.5))

ax = axes[0]
ax.plot(taus_fine, path[:, 1], color=COLORS["accent"], marker="o", markersize=2.5)
ax.fill_between(taus_fine, path[:, 1] - 1.96 * path_se[:, 1],
                path[:, 1] + 1.96 * path_se[:, 1],
                color=COLORS["accent"], alpha=0.15, linewidth=0)
ax.axhline(ols[1], color=COLORS["second"], linestyle="--", linewidth=1.0)
ax.fill_between([0.02, 0.98], ols[1] - 1.96 * ols_se[1], ols[1] + 1.96 * ols_se[1],
                color=COLORS["second"], alpha=0.12, linewidth=0)
ax.annotate("least squares", (0.30, ols[1]), color=COLORS["second"], fontsize=7,
            xytext=(0, -12), textcoords="offset points")
ax.set_xlim(0.02, 0.98)
ax.set_xlabel(r"$\tau$")
ax.set_ylabel("slope on income")
ax.set_title("(a) slope path, bootstrap bands")

ax = axes[1]
shades = plt.cm.viridis(np.linspace(0.12, 0.88, len(TAUS)))
low = np.linspace(0, 800, 2)
ax.scatter(income, food, s=8, color=COLORS["muted"], alpha=0.5, linewidths=0, zorder=0)
for tau, colour in zip(TAUS, shades):
    ax.plot(low, beta[tau][0] + beta[tau][1] * low, color=colour, linewidth=0.9)
ax.axvspan(0, income.min(), color=COLORS["second"], alpha=0.12, linewidth=0, zorder=0)
ax.annotate("no data here", (10, 560), color=COLORS["second"], fontsize=7)
ax.set_xlim(0, 800)
ax.set_ylim(80, 620)
ax.set_xlabel("household income")
ax.set_ylabel("food expenditure")
ax.set_title("(b) the fan closes and crosses")
fig.tight_layout()
fig.savefig(figure_path("ch45", "qr_inference"))
