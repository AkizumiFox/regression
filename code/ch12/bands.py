"""Chapter 12, Section 4: the Working-Hotelling band for a straight line, compared with pointwise
intervals and prediction intervals.

Stack loss against air flow for Brownlee's 21 days of plant operation
(statsmodels.datasets.stackloss, public domain).
"""
import matplotlib.pyplot as plt
import numpy as np
import statsmodels.api as sm
from scipy import stats

from regbook import COLORS, Generated, figure_path, use_book_style

data = sm.datasets.stackloss.load_pandas().data

# <<band>>
x = data["AIRFLOW"].to_numpy()
y = data["STACKLOSS"].to_numpy()
n = len(y)
X = np.column_stack([np.ones(n), x])
C = np.linalg.inv(X.T @ X)
beta_hat = C @ X.T @ y
s = np.sqrt(np.sum((y - X @ beta_hat) ** 2) / (n - 2))

grid = np.linspace(45, 85, 401)
Z = np.column_stack([np.ones_like(grid), grid])
fit = Z @ beta_hat
se_mean = s * np.sqrt(np.einsum("ij,jk,ik->i", Z, C, Z))   # s * sqrt(1/n + (x - xbar)^2 / Sxx)
t_pt = stats.t.ppf(0.975, n - 2)                            # pointwise multiplier
w_wh = np.sqrt(2 * stats.f.ppf(0.95, 2, n - 2))             # Working-Hotelling multiplier
pointwise = (fit - t_pt * se_mean, fit + t_pt * se_mean)
band = (fit - w_wh * se_mean, fit + w_wh * se_mean)
se_new = np.sqrt(s**2 + se_mean**2)
prediction = (fit - t_pt * se_new, fit + t_pt * se_new)
print(f"slope {beta_hat[1]:.4f}, s = {s:.3f}; multipliers: t {t_pt:.3f}, Working-Hotelling {w_wh:.3f}")
# <</band>>

xbar, Sxx = x.mean(), np.sum((x - x.mean()) ** 2)
assert np.allclose(se_mean, s * np.sqrt(1 / n + (grid - xbar) ** 2 / Sxx))
i_min = np.argmin(band[1] - band[0])
assert abs(grid[i_min] - xbar) < 0.1                        # narrowest at the mean of x
half = lambda x0: w_wh * s * np.sqrt(1 / n + (x0 - xbar) ** 2 / Sxx)

# ---- the supremum identity behind the band (one simulated data set with known beta) -----------
rng = np.random.default_rng(1204)
beta_true = beta_hat.copy()
sigma = s
y_sim = X @ beta_true + sigma * rng.normal(size=n)
b_sim = C @ X.T @ y_sim
d = b_sim - beta_true
quad = d @ X.T @ X @ d                                      # ||X(b - beta)||^2 = ||M y - X beta||^2
ratio = (Z @ d) ** 2 / np.einsum("ij,jk,ik->i", Z, C, Z)
assert ratio.max() <= quad + 1e-9
u = X.T @ X @ d                                             # the maximizing direction
x_star = u[1] / u[0]
z_star = np.array([1.0, x_star])
assert np.isclose((z_star @ d) ** 2 / (z_star @ C @ z_star), quad)

# ---- simultaneous coverage: pointwise intervals versus the band -------------------------------
reps = 40_000
rng = np.random.default_rng(1205)
Y = X @ beta_true + sigma * rng.normal(size=(reps, n))
B = Y @ (C @ X.T).T
S = np.sqrt(np.sum((Y - B @ X.T) ** 2, axis=1) / (n - 2))
D = B - beta_true
rng_grid = np.linspace(x.min(), x.max(), 301)
Zr = np.column_stack([np.ones_like(rng_grid), rng_grid])
Tabs = np.abs(D @ Zr.T) / (S[:, None] * np.sqrt(np.einsum("ij,jk,ik->i", Zr, C, Zr))[None, :])
cov_pointwise_range = np.mean(Tabs.max(axis=1) <= t_pt)
cov_wh_range = np.mean(Tabs.max(axis=1) <= w_wh)
sup_all = np.einsum("ij,jk,ik->i", D, X.T @ X, D) / S**2   # sup over the whole line, attained a.s.
cov_wh_line = np.mean(sup_all <= w_wh**2)
cov_single = np.mean(Tabs[:, 150] <= t_pt)
print(f"simultaneous coverage over the observed range: pointwise {cov_pointwise_range:.4f}, "
      f"band {cov_wh_range:.4f}; band over the whole line {cov_wh_line:.4f}")
tol = 4 * np.sqrt(0.95 * 0.05 / reps)
assert abs(cov_wh_line - 0.95) < tol and abs(cov_single - 0.95) < tol
assert cov_wh_range >= cov_wh_line - 1e-12 and cov_pointwise_range < 0.93

# ---- multipliers sqrt(r F(r, nu)) relative to t -----------------------------------------------
dfs = [5, 10, 20, 60, np.inf]
dims = [1, 2, 3, 5, 10]


def mult(r, nu):
    if np.isinf(nu):
        return np.sqrt(stats.chi2.ppf(0.95, r))
    return np.sqrt(r * stats.f.ppf(0.95, r, nu))


table = {(r, nu): mult(r, nu) for r in dims for nu in dfs}
for nu in dfs[:-1]:
    assert np.isclose(table[(1, nu)], stats.t.ppf(0.975, nu))
for r in dims[1:]:
    for nu in dfs:
        assert table[(r, nu)] > table[(r - 1 if r - 1 in dims else 1, nu)]

# large-nu ratios of the band multiplier to the pointwise one, at two levels
ratio_inf = {a: np.sqrt(stats.chi2.ppf(1 - a, 2)) / stats.norm.ppf(1 - a / 2) for a in (0.05, 0.01)}
assert np.isclose(stats.chi2.ppf(0.99, 2), -2 * np.log(0.01))
assert ratio_inf[0.01] < ratio_inf[0.05]

gen = Generated("ch12", "bands", prefix="band")
gen.num("rinf05", ratio_inf[0.05], 3)
gen.num("rinf01", ratio_inf[0.01], 3)
gen.num("w01", np.sqrt(stats.chi2.ppf(0.99, 2)), 3)
gen.num("z01", stats.norm.ppf(0.995), 3)
gen.int("n", n)
gen.num("b0", beta_hat[0], 3)
gen.num("b1", beta_hat[1], 4)
gen.num("s", s, 3)
gen.num("xbar", xbar, 2)
gen.num("Sxx", Sxx, 1)
gen.num("t", t_pt, 3)
gen.num("wh", w_wh, 3)
gen.num("ratio", w_wh / t_pt, 3)
for x0 in (60, 80):
    gen.num(f"half{x0}", half(x0), 2)
    gen.num(f"pt{x0}", t_pt * s * np.sqrt(1 / n + (x0 - xbar) ** 2 / Sxx), 2)
    gen.num(f"pred{x0}", t_pt * s * np.sqrt(1 + 1 / n + (x0 - xbar) ** 2 / Sxx), 2)
    gen.num(f"fit{x0}", beta_hat[0] + beta_hat[1] * x0, 2)
gen.num("cov:ptrange", cov_pointwise_range, 4)
gen.num("cov:whrange", cov_wh_range, 4)
gen.num("cov:whline", cov_wh_line, 4)
gen.num("cov:single", cov_single, 4)
gen.int("reps", reps)
for (r, nu), v in table.items():
    gen.num(f"m:{r}:{'inf' if np.isinf(nu) else int(nu)}", v, 3)
gen.write()

# ---- figure ----------------------------------------------------------------------------------
use_book_style()
fig, ax = plt.subplots(figsize=(4.6, 3.2))
ax.fill_between(grid, *prediction, color=COLORS["grid"], alpha=0.6, linewidth=0,
                label="95% prediction intervals")
ax.plot(grid, band[0], color=COLORS["accent"], label="95% Working–Hotelling band")
ax.plot(grid, band[1], color=COLORS["accent"])
ax.plot(grid, pointwise[0], color=COLORS["second"], linestyle="--", linewidth=0.9,
        label="95% pointwise intervals")
ax.plot(grid, pointwise[1], color=COLORS["second"], linestyle="--", linewidth=0.9)
ax.plot(grid, fit, color=COLORS["ink"], linewidth=0.8)
ax.scatter(x, y, s=12, color=COLORS["ink"], zorder=3, linewidths=0)
ax.set_xlabel("air flow")
ax.set_ylabel("stack loss")
ax.set_xlim(grid[0], grid[-1])
ax.legend(loc="upper left", frameon=False)
fig.tight_layout()
fig.savefig(figure_path("ch12", "working_hotelling"))
