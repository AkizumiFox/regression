"""Chapter 45, Section 4: expectile regression on Engel's budget survey.

Expectiles are fitted by iterated weighted least squares, checked against the
defining equation and against least squares at tau = 1/2, and the level of the
expectile that matches each quantile is computed for the fitted residuals: the
correspondence depends on the distribution, which is the price of the smooth loss.
"""
import matplotlib.pyplot as plt
import numpy as np
import statsmodels.api as sm
from scipy import stats

from regbook import COLORS, Generated, figure_path, use_book_style

data = sm.datasets.engel.load_pandas().data
TAUS = (0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95)

# <<iwls>>
income = data["income"].to_numpy()
food = data["foodexp"].to_numpy()
n = len(food)
X = np.column_stack([np.ones(n), income])


def expectile_reg(X, y, tau, tol=1e-12, max_iter=100):
    """Expectile regression by iterated weighted least squares.

    The weight of observation i is tau if it sits above the current fit and
    1 - tau if it sits below; each step is a weighted least squares fit.
    """
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    for step in range(1, max_iter + 1):
        w = np.where(y >= X @ beta, tau, 1 - tau)
        new, *_ = np.linalg.lstsq(X * np.sqrt(w)[:, None], y * np.sqrt(w), rcond=None)
        if np.max(np.abs(new - beta)) <= tol * (1 + np.max(np.abs(beta))):
            return new, step
        beta = new
    raise RuntimeError("iterated weighted least squares did not converge")


fits = {tau: expectile_reg(X, food, tau) for tau in TAUS}
for tau in (0.10, 0.50, 0.90):
    beta, steps = fits[tau]
    print(f"tau = {tau:.2f}: intercept {beta[0]:7.2f}  slope {beta[1]:.4f}  ({steps} steps)")
# <</iwls>>

ols, _ = np.linalg.lstsq(X, food, rcond=None)[:2]
assert np.allclose(fits[0.50][0], ols)            # tau = 1/2 is least squares
for tau in TAUS:                                  # the defining equation holds
    r = food - X @ fits[tau][0]
    w = np.where(r >= 0, tau, 1 - tau)
    assert np.max(np.abs(X.T @ (w * r))) < 1e-8 * np.max(np.abs(X).T @ np.abs(w * r))
for lo, hi in zip(TAUS[:-1], TAUS[1:]):           # expectile curves are ordered here
    assert np.all(X @ fits[hi][0] > X @ fits[lo][0])
max_steps = max(steps for _, steps in fits.values())

# <<level>>
def expectile(y, tau):
    """The sample tau-expectile of a vector: the scalar case of the same iteration."""
    m = y.mean()
    for _ in range(200):
        w = np.where(y >= m, tau, 1 - tau)
        new = (w * y).sum() / w.sum()
        if abs(new - m) < 1e-12 * (1 + abs(m)):
            return new
        m = new
    return m


def level(y, mu):
    """tau = E(mu - Y)_+ / E|Y - mu|: the expectile level of the point mu."""
    return np.mean(np.maximum(mu - y, 0)) / np.mean(np.abs(y - mu))


resid = food - X @ ols
for tau in (0.1, 0.5, 0.9):
    assert abs(level(resid, expectile(resid, tau)) - tau) < 1e-9
print("the identity tau = E(mu - Y)_+ / E|Y - mu| holds at every level tried")
# <</level>>

# which quantile level does each expectile level correspond to?
def matching_quantile(y, tau):
    e = expectile(y, tau)
    return np.mean(y <= e)


alphas_engel = [matching_quantile(resid, t) for t in np.linspace(0.02, 0.98, 49)]
normal = stats.norm.rvs(size=200000, random_state=np.random.default_rng(1))
alphas_normal = [matching_quantile(normal, t) for t in np.linspace(0.02, 0.98, 49)]

# a heavier right tail moves the correspondence the other way
skewed = stats.lognorm.rvs(0.8, size=200000, random_state=np.random.default_rng(2))
alphas_skew = [matching_quantile(skewed, t) for t in np.linspace(0.02, 0.98, 49)]

# the standard normal 0.9-expectile, from the defining equation
from scipy.optimize import brentq
normal_eq = lambda m, t: (t * (stats.norm.pdf(m) - m * (1 - stats.norm.cdf(m)))
                          - (1 - t) * (stats.norm.pdf(m) + m * stats.norm.cdf(m)))
e90_normal = brentq(lambda m: normal_eq(m, 0.9), 0.0, 5.0)
assert abs(stats.norm.cdf(e90_normal) - matching_quantile(normal, 0.90)) < 2e-3

gen = Generated("ch45", "expectiles")
for tau in TAUS:
    key = f"{round(100 * tau):02d}"
    gen.num(f"slope{key}", fits[tau][0][1], 4)
    gen.num(f"int{key}", fits[tau][0][0], 1)
gen.int("max_steps", max_steps)
gen.num("q90_of_e90", matching_quantile(resid, 0.90), 3)
gen.num("q10_of_e10", matching_quantile(resid, 0.10), 3)
gen.num("normal_q90", matching_quantile(normal, 0.90), 3)
gen.num("skew_q90", matching_quantile(skewed, 0.90), 3)
gen.num("normal_q50", matching_quantile(normal, 0.50), 3)
gen.num("normal_e90", e90_normal, 4)
gen.num("normal_e90_level", stats.norm.cdf(e90_normal), 3)
gen.write()

# ---- figure ----------------------------------------------------------------
use_book_style()
fig, axes = plt.subplots(1, 2, figsize=(5.6, 2.5))

ax = axes[0]
grid = np.linspace(income.min(), 3000, 2)
shades = plt.cm.viridis(np.linspace(0.12, 0.88, len(TAUS)))
ax.scatter(income, food, s=8, color=COLORS["muted"], alpha=0.55, linewidths=0, zorder=0)
for tau, colour in zip(TAUS, shades):
    beta = fits[tau][0]
    ax.plot(grid, beta[0] + beta[1] * grid, color=colour, linewidth=1.0)
    ax.annotate(f"{tau:g}", (3000, beta[0] + beta[1] * 3000), color=colour,
                fontsize=6.5, xytext=(2, -2), textcoords="offset points")
ax.set_xlim(0, 3400)
ax.set_ylim(0, 2200)
ax.set_xlabel("household income")
ax.set_ylabel("food expenditure")
ax.set_title("(a) seven fitted expectile lines")

ax = axes[1]
ts = np.linspace(0.02, 0.98, 49)
ax.plot(ts, ts, color=COLORS["grid"], linewidth=0.8)
ax.plot(ts, alphas_normal, color=COLORS["accent"], label="normal")
ax.plot(ts, alphas_engel, color=COLORS["second"], label="Engel residuals")
ax.plot(ts, alphas_skew, color=COLORS["third"], label="lognormal")
ax.legend(frameon=False, loc="upper left")
ax.set_xlabel(r"expectile level $\tau$")
ax.set_ylabel("matching quantile level")
ax.set_title("(b) the match depends on the law")
fig.tight_layout()
fig.savefig(figure_path("ch45", "engel_expectiles"))
