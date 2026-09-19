"""Chapter 7, Section 5: sampling distributions in the normal linear model, for a rank-deficient X.

A one-way layout with three groups, an intercept and a covariate (p = 5, rank 4).
Under normal errors: an estimable contrast is normal, SSE / sigma^2 ~ chi^2(n - r),
the two are independent, and the standardized contrast is t(n - r).
Under Laplace errors: the moments of the contrast are unchanged, Var(SSE) follows the
kurtosis formula of Chapter 2, and the contrast and SSE are uncorrelated but dependent.
"""
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats

from regbook import COLORS, Generated, figure_path, use_book_style

# <<design>>
import numpy as np
from scipy import stats

groups = np.repeat([0, 1, 2], [4, 5, 6])
z = np.array([2.1, 3.4, 1.8, 4.0, 2.9, 3.3, 5.1, 2.2, 3.8, 1.5, 4.4, 2.7, 3.9, 5.0, 3.1])
n = len(z)
X = np.column_stack([np.ones(n), np.eye(3)[groups], z])   # p = 5 columns, rank 4
r = np.linalg.matrix_rank(X)
G = np.linalg.pinv(X.T @ X)                                # one generalized inverse
M = X @ G @ X.T
lam = np.array([0.0, 1.0, -1.0, 0.0, 0.0])                 # group 1 minus group 2: estimable
beta, sigma = np.array([1.0, 0.5, -0.5, 0.0, 0.8]), 2.0
var_lam = sigma ** 2 * lam @ G @ lam                       # Var(lambda^T beta_hat)
print(f"n = {n}, p = {X.shape[1]}, rank = {r}, Var(lambda^T beta_hat) = {var_lam:.4f}")
# <</design>>

assert r == 4
assert np.allclose(lam @ G @ X.T @ X, lam)                 # lambda in C(X^T): estimable

# <<simulate>>
rng = np.random.default_rng(75)
reps = 100_000
Y = X @ beta + sigma * rng.normal(size=(reps, n))          # one data set per row
est = Y @ (X @ G @ lam)                                    # lambda^T beta_hat = lambda^T G X^T y
sse = np.sum((Y - Y @ M) ** 2, axis=1)                     # y^T (I - M) y
s2 = sse / (n - r)
t = (est - lam @ beta) / np.sqrt(s2 * (lam @ G @ lam))
print(f"mean {est.mean():.4f}, variance {est.var():.4f}")
print(f"E(SSE)/sigma^2 = {np.mean(sse) / sigma**2:.3f}   (n - r = {n - r})")
print(f"corr(estimate, s^2) = {np.corrcoef(est, s2)[0, 1]:.4f}")
print("KS p-value, SSE/sigma^2 vs chi^2(n-r):", round(stats.kstest(sse / sigma**2, "chi2", args=(n - r,)).pvalue, 3))
print("KS p-value, t vs t(n-r):              ", round(stats.kstest(t, "t", args=(n - r,)).pvalue, 3))
# <</simulate>>

k = n - r
assert abs(est.mean() - lam @ beta) < 5 * np.sqrt(var_lam / reps)
assert abs(est.var() / var_lam - 1) < 0.02
assert abs(np.mean(sse) / sigma ** 2 / k - 1) < 0.01
assert stats.kstest(sse / sigma ** 2, "chi2", args=(k,)).pvalue > 0.001
assert stats.kstest(t, "t", args=(k,)).pvalue > 0.001
assert stats.kstest(t, "norm").pvalue < 1e-6               # t(11) is visibly not N(0,1)
assert abs(np.corrcoef(est, s2)[0, 1]) < 0.01
assert abs(np.corrcoef((est - lam @ beta) ** 2, s2)[0, 1]) < 0.01
cover = np.mean(np.abs(t) <= stats.t.ppf(0.975, k))
assert abs(cover - 0.95) < 0.003
out196 = np.mean(np.abs(t) > 1.96)                         # tail mass beyond the normal 97.5% point
t_out196 = 2 * stats.t.sf(1.96, k)
assert abs(out196 - t_out196) < 0.003
assert out196 > 0.05 + 0.02                                # clearly more than N(0,1) allows

# joint statement: two contrasts at once, quadratic form ~ chi^2(2)
L = np.array([[0.0, 1.0, -1.0, 0.0, 0.0], [0.0, 1.0, 0.0, -1.0, 0.0]])
D = Y @ (X @ G @ L.T) - L @ beta
W = L @ G @ L.T
q = np.einsum("ij,jk,ik->i", D, np.linalg.inv(W), D) / sigma ** 2
assert stats.kstest(q, "chi2", args=(2,)).pvalue > 0.001

# <<laplace>>
E_lap = rng.laplace(scale=sigma / np.sqrt(2), size=(reps, n))   # variance sigma^2, mu_4 = 6 sigma^4
Y_lap = X @ beta + E_lap
est_lap = Y_lap @ (X @ G @ lam)
sse_lap = np.sum((Y_lap - Y_lap @ M) ** 2, axis=1)
h = np.diag(np.eye(n) - M)
var_sse_lap = (6 - 3) * sigma**4 * h @ h + 2 * sigma**4 * (n - r)
print(f"Var(SSE): normal theory {2 * sigma**4 * (n - r):.1f}, Laplace formula {var_sse_lap:.1f}, "
      f"simulated {sse_lap.var():.1f}")
print(f"corr(squared error of estimate, SSE) under Laplace errors: "
      f"{np.corrcoef((est_lap - lam @ beta) ** 2, sse_lap)[0, 1]:.3f}")
# <</laplace>>

assert abs(est_lap.var() / var_lam - 1) < 0.02
assert abs(sse_lap.var() / var_sse_lap - 1) < 0.03
corr_lap = np.corrcoef((est_lap - lam @ beta) ** 2, sse_lap)[0, 1]
assert corr_lap > 0.1                                     # dependent, though uncorrelated
assert abs(np.corrcoef(est_lap, sse_lap)[0, 1]) < 0.01

gen = Generated("ch07", "sampling_distributions", prefix="sd")
gen.int("n", n)
gen.int("r", r)
gen.int("k", k)
gen.num("var_lam", var_lam, 4)
gen.num("sim_var", est.var(), 4)
gen.num("mean_sse", np.mean(sse) / sigma ** 2, 3)
gen.num("corr", np.corrcoef(est, s2)[0, 1], 4)
gen.num("cover", cover, 4)
gen.num("tq", stats.t.ppf(0.975, k), 3)
gen.num("var_sse_normal", 2 * sigma ** 4 * k, 1)
gen.num("var_sse_lap", var_sse_lap, 1)
gen.num("sim_var_sse_lap", sse_lap.var(), 1)
gen.num("corr_lap", corr_lap, 3)
gen.num("sumh2", h @ h, 3)
gen.num("out196", out196, 4)
gen.num("t_out196", t_out196, 4)
gen.write()

use_book_style()
fig, axes = plt.subplots(1, 3, figsize=(6.2, 2.2))
ax = axes[0]
u = np.linspace(0, 35, 300)
ax.hist(sse / sigma ** 2, bins=80, range=(0, 35), density=True, color=COLORS["grid"])
ax.plot(u, stats.chi2.pdf(u, k), color=COLORS["accent"])
ax.set_xlabel(r"$\mathrm{SSE}/\sigma^2$")
ax.set_title(r"(a) $\chi^2(%d)$" % k)
ax = axes[1]
v = np.linspace(-6, 6, 300)
counts, edges = np.histogram(t, bins=60, range=(-6, 6))
dens = counts / (reps * np.diff(edges))
mid = 0.5 * (edges[:-1] + edges[1:])
keep = dens > 0
ax.plot(mid[keep], dens[keep], ls="none", marker="o", ms=2.2, color=COLORS["grid"], mec="0.45", mew=0.3,
        label="simulated")
ax.plot(v, stats.t.pdf(v, k), color=COLORS["accent"], label=r"$t(%d)$" % k)
ax.plot(v, stats.norm.pdf(v), color=COLORS["second"], ls="--", lw=0.9, label=r"$\mathrm{N}(0,1)$")
ax.set_yscale("log")
ax.set_ylim(1e-5, 1)
ax.set_xlabel(r"standardized contrast $T$")
ax.set_ylabel("density (log scale)")
ax.set_title(r"(b) $t(%d)$" % k)
ax.legend(frameon=False, fontsize=6, loc="lower center", handlelength=1.5)
ax = axes[2]
ax.scatter(est[:1500], s2[:1500], s=2, color=COLORS["accent"], alpha=0.5, linewidths=0)
ax.set_xlabel(r"$\boldsymbol{\lambda}^{\top}\hat{\boldsymbol{\beta}}$")
ax.set_ylabel(r"$s^2$")
ax.set_title("(c) independent")
fig.tight_layout()
fig.savefig(figure_path("ch07", "sampling_distributions"))
