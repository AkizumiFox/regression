"""Chapter 31, Section 2: equicorrelated clusters of unequal size.

With V = (1 - rho) I + rho Z Z' and a model [1, w] whose second column is centred inside
every cluster, the within-cluster slope is exactly best under ordinary least squares while
the overall level is not. Everything here is checked against closed forms.
"""
import matplotlib.pyplot as plt
import numpy as np

from regbook import COLORS, Generated, figure_path, use_book_style

# <<setup>>
import numpy as np

sizes = [2, 3, 6]
w = np.array([-1.0, 1.0, -1.0, 0.0, 1.0, -3.0, -2.0, -1.0, 1.0, 2.0, 3.0])
n = len(w)
Z = np.zeros((n, len(sizes)))                       # cluster indicators
Z[np.arange(n), np.repeat(np.arange(len(sizes)), sizes)] = 1.0
X = np.column_stack([np.ones(n), w])
rho = 0.5
V = (1 - rho) * np.eye(n) + rho * Z @ Z.T

Ginv = np.linalg.inv(X.T @ X)
Cov_ols = Ginv @ X.T @ V @ X @ Ginv                 # sigma^2 = 1 throughout
Cov_gls = np.linalg.inv(X.T @ np.linalg.solve(V, X))
print("efficiency of OLS:  level %.4f   slope %.4f"
      % (Cov_gls[0, 0] / Cov_ols[0, 0], Cov_gls[1, 1] / Cov_ols[1, 1]))
# <</setup>>

assert np.allclose(Z.T @ w, 0)                      # w is centred inside every cluster
assert np.isclose(Cov_gls[1, 1], Cov_ols[1, 1])     # the slope loses nothing
assert Cov_gls[0, 0] < Cov_ols[0, 0] - 1e-8         # the level does

M = X @ Ginv @ X.T
for j, best in ((0, False), (1, True)):
    c = X @ Ginv[:, j]                              # lambda = e_j gives this vector c
    assert np.allclose((np.eye(n) - M) @ V @ c, 0, atol=1e-12) == best

# closed forms: the GLS level is a weighted mean of the cluster means
weights = np.array([g / (1 - rho + g * rho) for g in sizes])
ybar_weights = weights / weights.sum()
b_gls_row = np.linalg.solve(X.T @ np.linalg.solve(V, X), X.T @ np.linalg.inv(V))[0]
assert np.allclose(b_gls_row @ Z, ybar_weights)   # the cluster totals of the weights
assert np.isclose(Cov_gls[0, 0], 1 / weights.sum())
assert np.isclose(Cov_ols[0, 0], ((1 - rho) * n + rho * sum(g * g for g in sizes)) / n**2)

eff_level = Cov_gls[0, 0] / Cov_ols[0, 0]
lam = np.linalg.eigvalsh(V)
kantorovich = 4 * lam[0] * lam[-1] / (lam[0] + lam[-1]) ** 2
Es2 = np.trace((np.eye(n) - M) @ V) / (n - 2)
se_ratio_level = np.sqrt(Es2 * Ginv[0, 0] / Cov_ols[0, 0])
se_ratio_slope = np.sqrt(Es2 * Ginv[1, 1] / Cov_ols[1, 1])
print("Kantorovich bound %.4f, attained efficiency %.4f" % (kantorovich, eff_level))
print("reported/true s.e.: level %.3f  slope %.3f" % (se_ratio_level, se_ratio_slope))
assert kantorovich < eff_level

gen = Generated("ch31", "unequal_clusters", prefix="clus")
gen.int("n", n)
gen.text("sizes", ", ".join(str(g) for g in sizes))
gen.num("rho", rho, 1)
gen.num("eff_level", eff_level, 3)
gen.num("var_gls_level", Cov_gls[0, 0], 4)
gen.num("var_ols_level", Cov_ols[0, 0], 4)
gen.num("var_slope", Cov_gls[1, 1], 4)
gen.num("sum_weights", weights.sum(), 4)
gen.num("w1", weights[0], 3)
gen.num("w2", weights[1], 3)
gen.num("w3", weights[2], 3)
gen.num("kantorovich", kantorovich, 3)
gen.num("lam1", lam[-1], 3)
gen.num("lamn", lam[0], 3)
gen.num("se_ratio_level", se_ratio_level, 3)
gen.num("se_ratio_slope", se_ratio_slope, 3)
gen.num("Es2", Es2, 3)

# ---- how bad can the imbalance get? -----------------------------------------
grid = np.linspace(0.0, 0.98, 99)
layouts = {"2, 3, 6": [2, 3, 6], "1, 1, 20": [1, 1, 20], "4, 4, 4": [4, 4, 4]}
curves, bounds = {}, {}
for name, gs in layouts.items():
    nn = sum(gs)
    ww = np.array([g / (1 - grid + g * grid) for g in gs])
    var_gls = 1 / ww.sum(axis=0)
    var_ols = ((1 - grid) * nn + grid * sum(g * g for g in gs)) / nn**2
    curves[name] = var_gls / var_ols
    Zg = np.zeros((nn, len(gs)))
    Zg[np.arange(nn), np.repeat(np.arange(len(gs)), gs)] = 1.0
    bs = []
    for r in grid:
        ev = np.linalg.eigvalsh((1 - r) * np.eye(nn) + r * Zg @ Zg.T)
        bs.append(4 * ev[0] * ev[-1] / (ev[0] + ev[-1]) ** 2)
    bounds[name] = np.array(bs)

assert np.allclose(curves["4, 4, 4"], 1.0)                 # a balanced layout loses nothing
assert curves["1, 1, 20"].min() < 0.42
for name in layouts:
    assert np.all(curves[name] >= bounds[name] - 1e-12)     # the Kantorovich bound holds

gen.num("eff_worst", curves["1, 1, 20"].min(), 3)
gen.num("eff_worst_rho", grid[int(np.argmin(curves["1, 1, 20"]))], 2)
gen.write()

use_book_style()
fig, ax = plt.subplots(figsize=(3.6, 2.4))
for name, colour in zip(layouts, (COLORS["accent"], COLORS["second"], COLORS["third"])):
    ax.plot(grid, curves[name], color=colour, label=f"sizes {name}")
    ax.plot(grid, bounds[name], color=colour, linewidth=0.7, linestyle=":")
ax.set_xlabel(r"$\rho$")
ax.set_ylabel("efficiency of the OLS level")
ax.set_ylim(0, 1.05)
ax.legend(loc="lower left", frameon=False)
fig.tight_layout()
fig.savefig(figure_path("ch31", "unequal_clusters"))
