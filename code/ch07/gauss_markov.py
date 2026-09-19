"""Chapter 7, Sections 1-2: linear unbiased estimators of a slope, and the Gauss-Markov theorem.

Three linear unbiased estimators of the slope in a straight-line model are compared,
exactly (variance = sigma^2 ||a||^2) and by simulation with skewed, non-normal errors.
Then a nonlinear unbiased estimator (the midrange) beats the sample mean when the
errors are uniform: the Gauss-Markov theorem says nothing about nonlinear estimators.
"""
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats

from regbook import COLORS, Generated, figure_path, use_book_style

# <<design>>
import numpy as np

x = np.array([1.0, 2.0, 2.5, 3.0, 4.0, 5.0, 5.5, 7.0, 8.0, 9.0, 10.5, 12.0])
n = len(x)
X = np.column_stack([np.ones(n), x])
lam = np.array([0.0, 1.0])                      # target: the slope beta_1
M = X @ np.linalg.solve(X.T @ X, X.T)           # projection onto C(X)

# coefficient vectors a of three linear estimators a^T y of the slope
a_ls = X @ np.linalg.solve(X.T @ X, lam)        # least squares: a = X (X^T X)^{-1} lambda
a_end = np.zeros(n)
a_end[[0, -1]] = [-1.0, 1.0]
a_end /= x[-1] - x[0]                           # slope through the two end points
low, high = np.arange(n) < n // 2, np.arange(n) >= n // 2
a_grp = (high / high.sum() - low / low.sum()) / (x[high].mean() - x[low].mean())  # group means

for name, a in [("least squares", a_ls), ("end points", a_end), ("group means", a_grp)]:
    print(f"{name:14s} X^T a = {X.T @ a}   Var/sigma^2 = {a @ a:.5f}   Ma = a_ls: {np.allclose(M @ a, a_ls)}")
# <</design>>

for a in [a_ls, a_end, a_grp]:
    assert np.allclose(X.T @ a, lam)            # unbiased for the slope, whatever beta is
    assert np.allclose(M @ a, a_ls)             # projection proof: every LUE projects to a_ls
    assert a @ a >= a_ls @ a_ls - 1e-15         # Gauss-Markov
Sxx = np.sum((x - x.mean()) ** 2)
assert np.isclose(a_ls @ a_ls, 1 / Sxx)

# <<simulate>>
estimators = {"ls": a_ls, "end": a_end, "grp": a_grp}
rng = np.random.default_rng(7)
beta, sigma, reps = np.array([2.0, 0.5]), 1.0, 200_000
E = sigma * (rng.exponential(size=(reps, n)) - 1.0)   # skewed errors, mean 0, variance sigma^2
Y = X @ beta + E                                      # one simulated data set per row
sims = {name: Y @ a for name, a in estimators.items()}
for name, est in sims.items():
    print(f"{name:4s} mean {est.mean():.4f}   variance {est.var():.4f}")
# <</simulate>>

for name, a in estimators.items():
    est = sims[name]
    se_mean = np.sqrt(sigma ** 2 * (a @ a) / reps)
    assert abs(est.mean() - beta[1]) < 5 * se_mean
    assert abs(est.var() / (sigma ** 2 * (a @ a)) - 1) < 0.03

# Matrix form: Cov(A y) - Cov(beta_hat) is nonnegative definite for any LUE matrix A of beta.
B = rng.normal(size=(2, n))
A_other = np.linalg.solve(X.T @ X, X.T) + B @ (np.eye(n) - M)
assert np.allclose(A_other @ X, np.eye(2))
diff = A_other @ A_other.T - np.linalg.inv(X.T @ X)
assert np.all(np.linalg.eigvalsh(diff) > -1e-12)

# Exercise exr-opt-weighing (b): two +-1 weighing designs for k = 3 objects in n = 4 weighings
X1 = np.array([[1, 1, 1], [1, -1, -1], [-1, 1, -1], [-1, -1, 1]], dtype=float)
X2 = np.array([[1, 1, 1], [1, 1, -1], [1, -1, 1], [1, 1, 1]], dtype=float)
assert np.allclose(X1.T @ X1, 4 * np.eye(3))
C2 = np.linalg.inv(X2.T @ X2)
assert np.allclose(32 * C2, [[16, -8, -8], [-8, 12, 4], [-8, 4, 12]])
assert np.isclose(np.trace(C2), 5 / 4) and np.isclose(np.linalg.det(C2), 1 / 32)
assert np.isclose(np.trace(np.linalg.inv(X1.T @ X1)), 3 / 4)

# ---- nonlinear unbiased estimators can win: the location model with uniform errors
# <<midrange>>
ns = np.array([5, 10, 20, 50, 100, 200])
var_mean = 1 / (3 * ns)                          # errors uniform on (-1, 1): variance 1/3
var_mid = 2 / ((ns + 1) * (ns + 2))              # exact variance of the midrange
for m, v1, v2 in zip(ns, var_mean, var_mid):
    print(f"n = {m:3d}   Var(mean) = {v1:.5f}   Var(midrange) = {v2:.6f}   ratio = {v1 / v2:.1f}")
# <</midrange>>

sim_mid = {}
for m in [5, 20]:
    U = rng.uniform(-1, 1, size=(100_000, m)) + 3.0
    mid = (U.min(axis=1) + U.max(axis=1)) / 2
    assert abs(mid.mean() - 3.0) < 0.005
    assert abs(mid.var() / (2 / ((m + 1) * (m + 2))) - 1) < 0.03
    assert abs(U.mean(axis=1).var() / (1 / (3 * m)) - 1) < 0.03
    sim_mid[m] = mid.var()

gen = Generated("ch07", "gauss_markov", prefix="gm")
gen.int("n", n)
gen.num("sxx", Sxx, 3)
for name, a in estimators.items():
    gen.num(f"var:{name}", a @ a, 5)
    gen.num(f"eff:{name}", (a_ls @ a_ls) / (a @ a), 3)
    gen.num(f"simvar:{name}", sims[name].var(), 5)
    gen.num(f"simmean:{name}", sims[name].mean(), 4)
for m in [20, 100]:
    k = int(np.where(ns == m)[0][0])
    gen.num(f"varmean{m}", var_mean[k], 4)
    gen.num(f"varmid{m}", var_mid[k], 5)
    gen.num(f"ratio{m}", var_mean[k] / var_mid[k], 1)
gen.write()

# ---- figure ------------------------------------------------------------------
use_book_style()
fig, axes = plt.subplots(1, 2, figsize=(5.8, 2.5))
ax = axes[0]
grid = np.linspace(0.1, 0.9, 400)
styles = [("end", "end points", COLORS["second"], "--", 1.0), ("grp", "group means", COLORS["third"], "-.", 1.0),
          ("ls", "least squares", COLORS["accent"], "-", 1.8)]
for key, label, color, ls, lw in styles:              # least squares drawn last, heaviest
    kde = stats.gaussian_kde(sims[key][:40_000])
    ax.plot(grid, kde(grid), color=color, ls=ls, lw=lw,
            label=f"{label} (var {estimators[key] @ estimators[key]:.5f})")
for (key, label, color, ls, lw), ybar in zip(styles, [6.85, 6.4, 5.95]):   # mean +- 2 SD markers
    sd = np.sqrt(estimators[key] @ estimators[key])
    ax.errorbar(beta[1], ybar, xerr=2 * sd, color=color, lw=lw, capsize=2, capthick=lw, ls="none")
ax.text(0.9, 7.15, r"$\pm2$ SD", fontsize=6.5, ha="right", color=COLORS["muted"])
ax.plot([beta[1], beta[1]], [0, 5.9], color=COLORS["muted"], lw=0.6, ls=":")
ax.set_xlabel(r"estimate of the slope $\beta_1$")
ax.set_ylabel("density")
ax.set_title("(a) three linear unbiased estimators")
ax.legend(frameon=False, fontsize=6.5, loc="upper left")
ax.set_xlim(0.1, 0.9)
ax.set_ylim(0, 10.2)
ax = axes[1]
ax.loglog(ns, var_mean, "o-", ms=3, color=COLORS["accent"], label="sample mean (BLUE)")
ax.loglog(ns, var_mid, "s-", ms=3, color=COLORS["second"], label="midrange (nonlinear)")
ax.set_xlabel("sample size $n$")
ax.set_ylabel("variance")
ax.set_title("(b) uniform errors")
ax.legend(frameon=False, fontsize=7)
fig.tight_layout()
fig.savefig(figure_path("ch07", "gauss_markov"))
