"""Chapter 29, Section 1: optimism of the training error and degrees of freedom.

(a) Least squares with a WRONG mean (a quadratic trend fitted by a straight line plus two
    irrelevant regressors): the expected optimism is still 2 sigma^2 p / n.
(b) Ridge regression, a linear smoother: optimism 2 sigma^2 tr(S) / n.
(c) Thresholding in the orthonormal (sequence) model: soft thresholding has
    df = E #{|Y_i| > lambda} (Stein); hard thresholding has the extra jump term.
(d) Best subset of size k (keep the k largest |Y_j|): its degrees of freedom exceed k.
"""
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import norm

from regbook import COLORS, Generated, figure_path, use_book_style

gen = Generated("ch29", "optimism", prefix="opt")

# <<penalty>>
import numpy as np
rng = np.random.default_rng(2901)
n, sigma, reps = 30, 1.0, 20_000

t = np.linspace(-1, 1, n)
X = np.column_stack([np.ones(n), t, rng.normal(size=(n, 2))])   # line + 2 irrelevant columns
mu = 1 + t + 1.5 * t**2                                         # the true mean is curved
p = X.shape[1]
M = X @ np.linalg.solve(X.T @ X, X.T)                           # projection onto C(X)

Y = mu + sigma * rng.normal(size=(reps, n))                     # one data set per row
fit = Y @ M                                                     # M is symmetric
train = np.mean((Y - fit) ** 2, axis=1)                         # training error
err_in = sigma**2 + np.mean((fit - mu) ** 2, axis=1)            # in-sample prediction error
print(f"average optimism {np.mean(err_in - train):.4f},  2 sigma^2 p / n = {2 * sigma**2 * p / n:.4f}")
# <</penalty>>

omega_ls = np.mean(err_in - train)
assert abs(omega_ls - 2 * p / n) < 0.01
bias2 = np.sum(((np.eye(n) - M) @ mu) ** 2)
assert bias2 > 1.0                                             # the model is genuinely wrong
# exact expectations
assert abs(np.mean(train) - (sigma**2 * (n - p) + bias2) / n) < 0.01
assert abs(np.mean(err_in) - (sigma**2 * (n + p) + bias2) / n) < 0.01
gen.int("n", n)
gen.int("p", p)
gen.num("theory", 2 * p / n, 4)
gen.num("omega_ls", omega_ls, 4)
gen.num("bias2", bias2 / n, 4)

# <<ridge>>
lam = 5.0
S = X @ np.linalg.solve(X.T @ X + lam * np.eye(p), X.T)         # ridge smoother
fit_r = Y @ S.T
train_r = np.mean((Y - fit_r) ** 2, axis=1)
err_in_r = sigma**2 + np.mean((fit_r - mu) ** 2, axis=1)
omega_r = np.mean(err_in_r - train_r)
print(f"ridge: optimism {omega_r:.4f},  2 sigma^2 tr(S) / n = {2 * np.trace(S) / n:.4f}")
# <</ridge>>
assert abs(omega_r - 2 * np.trace(S) / n) < 0.01
d = np.linalg.svd(X, compute_uv=False)
assert np.isclose(np.trace(S), np.sum(d**2 / (d**2 + lam)))
gen.num("ridge_trS", np.trace(S), 3)
gen.num("ridge_omega", omega_r, 4)

# ---- (c) thresholding in the sequence model --------------------------------------
lam_t = np.sqrt(2.0)


def hard_df(theta, lam):
    """Per-coordinate degrees of freedom of hard thresholding at lam (sigma = 1)."""
    return norm.sf(lam - theta) + norm.cdf(-lam - theta) + lam * (norm.pdf(lam - theta) + norm.pdf(lam + theta))


Z = rng.normal(size=2_000_000)
for theta in (0.0, 1.0, 3.0):
    Yt = theta + Z
    hard = Yt * (np.abs(Yt) > lam_t)
    soft = np.sign(Yt) * np.maximum(np.abs(Yt) - lam_t, 0)
    cov_h = np.mean(hard * Z)                                   # Cov(mu_hat, Y) = E[mu_hat (Y - theta)]
    cov_s = np.mean(soft * Z)
    assert abs(cov_h - hard_df(theta, lam_t)) < 0.005
    assert abs(cov_s - np.mean(np.abs(Yt) > lam_t)) < 0.005   # Stein: E[derivative]
h0 = hard_df(0.0, lam_t)
gen.num("hard_df0", h0, 3)
gen.num("hard_sel0", 2 * norm.sf(lam_t), 3)
gen.num("hard_jump0", 2 * lam_t * norm.pdf(lam_t), 3)
gen.num("hard_df1", hard_df(1.0, lam_t), 3)
gen.num("hard_df3", hard_df(3.0, lam_t), 3)

# ---- (d) best subset of size k in the orthonormal model -----------------------------
# <<subset>>
m, reps_s = 20, 20_000
Zs = rng.normal(size=(reps_s, m))


def subset_df(theta, k):
    """Monte Carlo df of 'keep the k largest |Y_j|' when Y ~ N(theta, I)."""
    Yb = theta + Zs
    keep = np.argsort(-np.abs(Yb), axis=1)[:, :k]
    fit = np.zeros_like(Yb)
    np.put_along_axis(fit, keep, np.take_along_axis(Yb, keep, axis=1), axis=1)
    return np.mean(np.sum(fit * Zs, axis=1))                   # sum_j Cov(fit_j, Y_j)


null = np.zeros(m)
sparse = np.r_[np.full(5, 6.0), np.zeros(m - 5)]               # five large means
for k in (1, 5, 10):
    print(f"k = {k:2d}: df null {subset_df(null, k):5.2f}, df sparse {subset_df(sparse, k):5.2f}")
# <</subset>>

ks = np.arange(0, m + 1)
df_null = np.array([subset_df(null, k) for k in ks])
df_sparse = np.array([subset_df(sparse, k) for k in ks])
assert abs(df_null[0]) < 1e-12 and abs(df_null[-1] - m) < 0.2
assert np.all(df_null[1:m] > ks[1:m])                            # search costs extra df
assert abs(df_sparse[5] - 5) < 0.3                               # obvious choice: little extra
assert df_sparse[8] > 8 + 3
gen.int("m", m)
gen.num("df_null1", df_null[1], 2)
gen.num("df_null5", df_null[5], 2)
gen.num("df_null10", df_null[10], 2)
gen.num("df_sparse5", df_sparse[5], 2)
gen.num("df_sparse6", df_sparse[6], 2)
gen.num("df_sparse10", df_sparse[10], 2)
gen.write()

use_book_style()
fig, ax = plt.subplots(figsize=(4.2, 2.9))
ax.plot(ks, ks, color=COLORS["muted"], linewidth=0.8, linestyle="--", label="df = k (fixed subset)")
ax.plot(ks, df_null, "o-", color=COLORS["accent"], markersize=3, label="all means zero")
ax.plot(ks, df_sparse, "s-", color=COLORS["second"], markersize=3, label="five means equal to 6")
ax.set_xlabel("subset size k")
ax.set_ylabel("degrees of freedom")
ax.set_xticks([0, 5, 10, 15, 20])
ax.legend(frameon=False, loc="lower right")
fig.tight_layout()
fig.savefig(figure_path("ch29", "subset_df"))
