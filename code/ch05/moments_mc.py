"""Chapter 5, Section 5: the mean and covariance of the least squares estimator.

(i) A quadratic regression on twelve design points with skewed (centred exponential)
    errors: the Monte Carlo mean and covariance of beta_hat match beta and
    sigma^2 (X^T X)^{-1}, which need only second moments.
(ii) Simple regression with n = 20 points in [0, 10] under three designs. The slope
    variance is sigma^2 / S_xx, so spreading the design out pays.
"""
import matplotlib.pyplot as plt
import numpy as np

from regbook import COLORS, Generated, figure_path, use_book_style

gen = Generated("ch05", "moments_mc", prefix="mc")

# <<quadratic>>
rng = np.random.default_rng(5)
x = np.arange(1.0, 13.0)                        # twelve design points
X = np.column_stack([np.ones_like(x), x, x**2])
beta = np.array([2.0, 0.5, -0.03])
sigma = 1.5
reps = 100_000

A = np.linalg.solve(X.T @ X, X.T)               # beta_hat = A y
E = sigma * (rng.exponential(size=(reps, len(x))) - 1.0)   # skewed, mean 0, sd sigma
B = (X @ beta + E) @ A.T                        # one estimate per row

print("mean of estimates:", np.round(B.mean(axis=0), 4), " true:", beta)
print("Monte Carlo covariance:\n", np.round(np.cov(B.T), 5))
print("sigma^2 (X^T X)^{-1}:\n", np.round(sigma**2 * np.linalg.inv(X.T @ X), 5))
# <</quadratic>>

cov_theory = sigma**2 * np.linalg.inv(X.T @ X)
cov_mc = np.cov(B.T)
se_mean = np.sqrt(np.diag(cov_theory) / reps)
assert np.all(np.abs(B.mean(axis=0) - beta) < 4 * se_mean)
rel = np.abs(cov_mc - cov_theory) / np.sqrt(np.outer(np.diag(cov_theory), np.diag(cov_theory)))
assert rel.max() < 0.02
# correlation of the estimates: nearly collinear columns x and x^2
corr_theory = cov_theory / np.sqrt(np.outer(np.diag(cov_theory), np.diag(cov_theory)))
assert corr_theory[1, 2] < -0.9

gen.text("reps", f"{reps:,}")
gen.num("relmax", rel.max(), 3)
gen.num("corr12", corr_theory[1, 2], 3)
for j in range(3):
    gen.num(f"sd{j}", np.sqrt(cov_theory[j, j]), 4)
    gen.num(f"sdmc{j}", np.sqrt(cov_mc[j, j]), 4)
    gen.num(f"mean{j}", B.mean(axis=0)[j], 4)

# <<designs>>
n = 20
designs = {
    "equally spaced": np.linspace(0, 10, n),
    "two ends": np.repeat([0.0, 10.0], n // 2),
    "middle only": np.linspace(4, 6, n),
}
for name, xd in designs.items():
    Sxx = np.sum((xd - xd.mean()) ** 2)
    print(f"{name:15s} S_xx = {Sxx:7.2f}   sd(slope) = sigma * {1 / np.sqrt(Sxx):.4f}")
# <</designs>>

sims = {}
rng2 = np.random.default_rng(55)
for name, xd in designs.items():
    Sxx = np.sum((xd - xd.mean()) ** 2)
    xc = xd - xd.mean()
    Y = 1.0 + 0.5 * xd + rng2.standard_normal((20_000, n))
    slopes = (Y - Y.mean(axis=1, keepdims=True)) @ xc / Sxx
    assert abs(slopes.std() * np.sqrt(Sxx) - 1) < 0.03
    sims[name] = slopes
    key = {"equally spaced": "eq", "two ends": "ends", "middle only": "mid"}[name]
    gen.num("Sxx" + key, Sxx, 2)
    gen.num("sd" + key, 1 / np.sqrt(Sxx), 4)
S = {k: np.sum((v - v.mean()) ** 2) for k, v in designs.items()}
assert S["two ends"] == max(S.values()) and np.isclose(S["two ends"], n * 25)
gen.num("ratio", np.sqrt(S["middle only"] / S["two ends"]), 3)
gen.write()

# ---- figure -------------------------------------------------------------------
use_book_style()
fig, axes = plt.subplots(1, 2, figsize=(5.6, 2.4), gridspec_kw={"width_ratios": [1, 1.25]})
ax = axes[0]
for k, (name, xd) in enumerate(designs.items()):
    ax.scatter(xd, np.full(n, 2 - k) + np.random.default_rng(k).uniform(-0.12, 0.12, n), s=6,
               color=[COLORS["accent"], COLORS["second"], COLORS["third"]][k], linewidths=0)
ax.set_yticks([2, 1, 0])
ax.set_yticklabels(list(designs))
ax.set_ylim(-0.6, 2.6)
ax.set_xlabel("x")
ax.set_title("(a) three designs, n = 20")
ax = axes[1]
bins = np.linspace(-0.5, 1.5, 81)
for k, (name, sl) in enumerate(sims.items()):
    ax.hist(sl, bins=bins, histtype="step", density=True,
            color=[COLORS["accent"], COLORS["second"], COLORS["third"]][k], label=name)
ax.axvline(0.5, color=COLORS["muted"], linewidth=0.6)
ax.set_xlabel("least squares slope")
ax.set_yticks([])
ax.set_title(r"(b) sampling distributions, $\beta_1=0.5$")
ax.legend(frameon=False, fontsize=7, loc="upper left")
fig.tight_layout()
fig.savefig(figure_path("ch05", "moments_designs"))
