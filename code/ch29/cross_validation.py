"""Chapter 29, Section 3: cross-validation.

(a) The leave-one-out shortcut for ridge regression (intercept unpenalized) on the state
    data, against brute-force refitting; generalized cross-validation.
(b) Bias and variance of K-fold cross-validation for least squares with k = 10 normal
    regressors and n = 40: E(CV_K) equals the expected error of a fit to n - n/K cases, given
    by the normal-theory formula of Chapter 14; and CV_K hardly tracks the error of the fit at hand.
"""
import matplotlib.pyplot as plt
import numpy as np

from regbook import COLORS, Generated, figure_path, use_book_style

gen = Generated("ch29", "cross_validation", prefix="cv")

# <<shortcut>>
import numpy as np
import statsmodels.api as sm
data = sm.datasets.statecrime.load_pandas().data.drop(index="District of Columbia")
names = ["hs_grad", "poverty", "single", "white", "urban"]
y = np.log(data["violent"].to_numpy())
Zs = data[names].to_numpy()
Zs = (Zs - Zs.mean(axis=0)) / Zs.std(axis=0)               # standardized regressors
X = np.column_stack([np.ones(len(y)), Zs])
n, p = X.shape
Omega = np.diag([0.0] + [1.0] * (p - 1))                    # ridge penalty, intercept left free


def ridge_fit(X, y, lam):
    return np.linalg.solve(X.T @ X + lam * Omega, X.T @ y)


lam = 5.0
S = X @ np.linalg.solve(X.T @ X + lam * Omega, X.T)          # the smoother matrix
resid = y - S @ y
loo_short = resid / (1 - np.diag(S))                         # the shortcut
loo_brute = np.array([y[i] - X[i] @ ridge_fit(np.delete(X, i, 0), np.delete(y, i), lam)
                      for i in range(n)])                    # n refits
print("largest discrepancy:", np.max(np.abs(loo_short - loo_brute)))
cv_loo = np.mean(loo_short**2)
gcv = np.mean(resid**2) / (1 - np.trace(S) / n) ** 2
print(f"CV_n = {cv_loo:.5f}, GCV = {gcv:.5f}, tr(S) = {np.trace(S):.3f}")
# <</shortcut>>
assert np.max(np.abs(loo_short - loo_brute)) < 1e-10
assert np.all(np.diag(S) < 1)
gen.num("lam", lam, 1)
gen.num("cv_loo", cv_loo, 4)
gen.num("gcv", gcv, 4)
gen.num("trS", np.trace(S), 3)

# least squares (lam = 0): CV_n = PRESS / n, GCV = n s^2 / (n - p)
H = X @ np.linalg.solve(X.T @ X, X.T)
e = y - H @ y
press = np.sum((e / (1 - np.diag(H))) ** 2)
gcv_ls = np.mean(e**2) / (1 - p / n) ** 2
assert np.isclose(gcv_ls, n * (e @ e / (n - p)) / (n - p))
gen.num("press_ls", press / n, 4)
gen.num("gcv_ls", gcv_ls, 4)

# ---- (b) bias and variance of K-fold CV ----------------------------------------------------------
# <<kfold>>
rng = np.random.default_rng(2905)
n_s, k_s, sigma = 40, 10, 1.0
beta = np.r_[1.0, np.full(k_s, 0.5)]                         # intercept and ten slopes


def err_theory(m):
    """Expected error of a least squares fit to m cases (Chapter 14, normal regressors)."""
    return sigma**2 * (1 + 1 / m) * (m - 2) / (m - k_s - 2)


def cv_kfold(X, y, K):
    folds = np.arange(len(y)) % K                            # equal folds (cases are exchangeable)
    sq = np.empty(len(y))
    for f in range(K):
        tr, te = folds != f, folds == f
        b, *_ = np.linalg.lstsq(X[tr], y[tr], rcond=None)
        sq[te] = (y[te] - X[te] @ b) ** 2
    return sq.mean()


Ks = [2, 5, 10, n_s]


def simulate(reps):
    cv = np.empty((reps, len(Ks)))
    err_cond = np.empty(reps)                                 # error of the fit at hand
    for r in range(reps):
        Xr = np.column_stack([np.ones(n_s), rng.normal(size=(n_s, k_s))])
        yr = Xr @ beta + sigma * rng.normal(size=n_s)
        b, *_ = np.linalg.lstsq(Xr, yr, rcond=None)
        err_cond[r] = sigma**2 + np.sum((b - beta) ** 2)     # new case X0 ~ N(0, I)
        cv[r] = [cv_kfold(Xr, yr, K) for K in Ks]
    return cv, err_cond


cv, err_cond = simulate(300)                                 # the book's figures use 4000
for j, K in enumerate(Ks):
    print(f"K = {K:2d}: mean {cv[:, j].mean():.3f} (theory {err_theory(n_s - n_s // K):.3f}), sd {cv[:, j].std():.3f}")
print(f"error of the fit at hand: mean {err_cond.mean():.3f} (theory {err_theory(n_s):.3f})")
print("correlation of CV_10 with it:", np.corrcoef(cv[:, 2], err_cond)[0, 1])
# <</kfold>>
reps = 4000
cv, err_cond = simulate(reps)
for j, K in enumerate(Ks):
    th = err_theory(n_s - n_s // K)
    assert abs(cv[:, j].mean() / th - 1) < 0.03
    gen.num(f"mean_{K}", cv[:, j].mean(), 3)
    gen.num(f"theory_{K}", th, 3)
    gen.num(f"sd_{K}", cv[:, j].std(), 3)
assert abs(err_cond.mean() / err_theory(n_s) - 1) < 0.02
corr = np.corrcoef(cv[:, 2], err_cond)[0, 1]
corr_loo = np.corrcoef(cv[:, 3], err_cond)[0, 1]
assert abs(corr) < 0.3
gen.num("err_n", err_theory(n_s), 3)
gen.num("err_mean", err_cond.mean(), 3)
gen.num("err_sd", err_cond.std(), 3)
gen.num("corr10", corr, 3)
gen.num("corr_loo", corr_loo, 3)
gen.int("n_s", n_s)
gen.int("k_s", k_s)
gen.int("reps", reps)
gen.write()

use_book_style()
fig, axes = plt.subplots(1, 2, figsize=(5.8, 2.7))
ax = axes[0]
xs = np.arange(len(Ks))
means, sds = cv.mean(axis=0), cv.std(axis=0)
ax.errorbar(xs, means, yerr=sds, fmt="o", color=COLORS["accent"], capsize=3, markersize=4, label="CV (mean, sd)")
ax.plot(xs, [err_theory(n_s - n_s // K) for K in Ks], "x", color=COLORS["second"], markersize=6,
        label=r"Err$(n-n/K)$")
ax.axhline(err_theory(n_s), color=COLORS["muted"], linestyle="--", linewidth=0.8)
ax.text(2.45, err_theory(n_s) - 0.25, r"Err$(n)$", fontsize=7, color=COLORS["muted"])
ax.set_xticks(xs, ["2", "5", "10", "n = 40"])
ax.set_xlabel("K")
ax.set_ylabel("estimated prediction error")
ax.set_title("(a) bias and spread")
ax.legend(frameon=False, loc="upper right")
ax = axes[1]
ax.scatter(err_cond, cv[:, 2], s=3, color=COLORS["accent"], alpha=0.35, linewidths=0)
lim = [0.9, 4.0]
ax.plot(lim, lim, color=COLORS["muted"], linewidth=0.8, linestyle="--")
ax.set_xlim(lim)
ax.set_ylim(0.3, 4.5)
ax.set_xlabel("error of the fitted model")
ax.set_ylabel(r"CV$_{10}$")
ax.set_title("(b) what CV tracks")
fig.tight_layout()
fig.savefig(figure_path("ch29", "cv_bias_variance"))
