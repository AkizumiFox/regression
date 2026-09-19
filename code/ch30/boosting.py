"""Chapter 30, Section 4: componentwise L2 boosting.

Simulated data: n = 60 training cases, p = 30 correlated regressors (correlation 0.5^|j-k|), five
nonzero coefficients, and the exact prediction risk computed from the design law. The script checks the one-step
decrease formula and the convergence to least squares (full column rank) at the guaranteed
geometric rate, traces the test error along the iterations for two step sizes, chooses the
stopping iteration by 5-fold cross-validation, and compares two measures of degrees of freedom:
the trace of the boosting operator (selection held fixed) and the covariance degrees of freedom
estimated by simulation with the design and mean held fixed.
"""
import matplotlib.pyplot as plt
import numpy as np

from regbook import COLORS, Generated, figure_path, use_book_style


# <<boost>>
def l2_boost(X, y, nu, steps):
    """Componentwise L2 boosting from b = 0; returns the coefficient path (steps + 1) x p."""
    n, p = X.shape
    norms2 = np.sum(X ** 2, axis=0)
    b = np.zeros(p)
    r = y.astype(float).copy()
    path = [b.copy()]
    for _ in range(steps):
        corr = X.T @ r
        j = np.argmax(corr ** 2 / norms2)               # largest drop in the residual sum of squares
        step = nu * corr[j] / norms2[j]
        b[j] += step
        r -= step * X[:, j]
        path.append(b.copy())
    return np.array(path)
# <</boost>>


# <<data>>
rng = np.random.default_rng(3004)
n, p, sigma = 60, 30, 1.5
C = 0.5 ** np.abs(np.subtract.outer(np.arange(p), np.arange(p)))
L = np.linalg.cholesky(C)
beta = np.zeros(p)
beta[[0, 4, 9, 14, 19]] = [2.0, -1.5, 1.0, 1.0, -0.5]

def draw(m):
    Z = rng.normal(size=(m, p)) @ L.T
    return Z, Z @ beta + sigma * rng.normal(size=m)

X, y = draw(n)
mx, sx, my = X.mean(axis=0), X.std(axis=0), y.mean()
Xs, ys = (X - mx) / sx, y - my                          # centre and scale by the training data

def risk(path):
    """Expected squared error on a new case, for each row b of path (exact, from the design law)."""
    d = beta[None, :] - path / sx                       # coefficient error on the original scale
    offset = path @ (mx / sx) - my                      # error in the fitted intercept
    return np.sum((d @ C) * d, axis=1) + offset ** 2 + sigma ** 2
# <</data>>

b_ls = np.linalg.lstsq(Xs, ys, rcond=None)[0]
Q_star = 0.5 * np.sum((ys - Xs @ b_ls) ** 2)

# ---- convergence to least squares, with the step-decrease identity and the geometric bound ------
nu = 0.1
long_path = l2_boost(Xs, ys, nu, 50000)
Qs = 0.5 * np.sum((ys[None, :] - long_path @ Xs.T) ** 2, axis=1)
norms2 = np.sum(Xs ** 2, axis=0)
# decrease = nu (1 - nu/2) max_j (x_j'r)^2/||x_j||^2, checked at every step of the path
R_path = ys[None, :] - long_path[:-1] @ Xs.T
drops = nu * (1 - nu / 2) * np.max((R_path @ Xs) ** 2 / norms2, axis=1)
assert np.allclose(Qs[:-1] - Qs[1:], drops, rtol=1e-6, atol=1e-12)
s_min = np.linalg.svd(Xs, compute_uv=False).min()
kappa = 2 * nu * (1 - nu / 2) * s_min ** 2 / (p * norms2.max())
ms = np.arange(len(Qs))
assert np.all(Qs - Q_star <= (1 - kappa) ** ms * (Qs[0] - Q_star) * (1 + 1e-9) + 1e-12)
assert np.max(np.abs(long_path[-1] - b_ls)) < 1e-5
dist_ls = {m: np.linalg.norm(long_path[m] - b_ls) / np.linalg.norm(b_ls) for m in [100, 1000, 50000]}

# ---- test error along the path ----------------------------------------------------------------
# <<test>>
M = 1500
test_err = {}
for step in [0.1, 1.0]:
    test_err[step] = risk(l2_boost(Xs, ys, step, M))
err_ls = risk(b_ls[None, :])[0]
m_best = int(np.argmin(test_err[0.1]))
print(f"least squares: risk {err_ls:.3f}")
print(f"nu = 0.1: best m = {m_best}, risk {test_err[0.1][m_best]:.3f}")
print(f"nu = 1.0: best m = {np.argmin(test_err[1.0])}, risk {test_err[1.0].min():.3f}")
# <</test>>
assert test_err[0.1][m_best] < err_ls
near = {s: np.flatnonzero(test_err[s] <= 1.05 * test_err[s].min()) for s in test_err}
near_lo = {s: int(near[s].min()) for s in near}
near_hi = {s: int(near[s].max()) for s in near}
assert near_hi[0.1] - near_lo[0.1] > 10 * (near_hi[1.0] - near_lo[1.0])
assert test_err[0.1][-1] > test_err[0.1][m_best] and np.isclose(test_err[1.0][-1], err_ls, rtol=0.01)
active_best = int(np.sum(np.abs(l2_boost(Xs, ys, 0.1, m_best)[-1]) > 0))

# ---- 5-fold cross-validation for the stopping iteration ------------------------------------------
folds = np.arange(n) % 5
rng.shuffle(folds)
cv = np.zeros(M + 1)
for k in range(5):
    tr, te = folds != k, folds == k
    mxk, sxk, myk = X[tr].mean(axis=0), X[tr].std(axis=0), y[tr].mean()
    pk = l2_boost((X[tr] - mxk) / sxk, y[tr] - myk, 0.1, M)
    cv += np.sum((y[te][None, :] - myk - pk @ ((X[te] - mxk) / sxk).T) ** 2, axis=1)
cv /= n
m_cv = int(np.argmin(cv))
err_cv = test_err[0.1][m_cv]
assert err_cv < err_ls

# ---- degrees of freedom: trace of the boosting operator vs covariance df -------------------------
steps_df = 120
checkpoints = [5, 10, 20, 40, 80, 120]
mu = Xs @ (beta * sx)                                    # the true mean, in the centred scale
mu -= mu.mean()
R = 2000
fits = np.zeros((R, len(checkpoints), n))
ys_all = np.zeros((R, n))
traces = np.zeros((R, len(checkpoints)))
for rep in range(R):
    yr = mu + sigma * rng.normal(size=n)
    ys_all[rep] = yr
    path = l2_boost(Xs, yr, nu, steps_df)
    # boosting operator B_m = I - prod_k (I - nu H_{j_k}) for the selected sequence
    js = [int(np.flatnonzero(path[k + 1] != path[k])[0]) for k in range(steps_df)]
    Rm = np.eye(n)                                       # R_m = I - B_m
    for k, j in enumerate(js, start=1):
        x = Xs[:, j]
        Rm -= nu * np.outer(x, x @ Rm) / (x @ x)
        if k in checkpoints:
            c = checkpoints.index(k)
            traces[rep, c] = n - np.trace(Rm)
            fits[rep, c] = (np.eye(n) - Rm) @ yr
            assert np.allclose(fits[rep, c], Xs @ path[k])
yc = ys_all - ys_all.mean(axis=0)
df_cov = np.array([np.sum(np.mean(yc * (fits[:, c] - fits[:, c].mean(axis=0)), axis=0)) / sigma ** 2
                   for c in range(len(checkpoints))])
df_trace = traces.mean(axis=0)
# Monte Carlo standard error of the covariance df at each checkpoint
se_cov = np.array([np.std(np.sum(yc * (fits[:, c] - fits[:, c].mean(axis=0)), axis=1)) / sigma ** 2
                   / np.sqrt(R) for c in range(len(checkpoints))])
print("m       ", checkpoints)
print("trace   ", np.round(df_trace, 1))
print("cov df  ", np.round(df_cov, 1))
assert np.all(df_cov[:4] - df_trace[:4] > 3 * se_cov[:4])  # selection adds degrees of freedom

# a small version for the runnable cell
# <<df_small>>
def operator_trace(X, path, nu):
    """Trace of B_m = I - prod_k (I - nu H_{j_k}), with the selected columns held fixed."""
    n = X.shape[0]
    Rm = np.eye(n)
    for k in range(len(path) - 1):
        j = int(np.flatnonzero(path[k + 1] != path[k])[0])
        x = X[:, j]
        Rm -= nu * np.outer(x, x @ Rm) / (x @ x)
    return n - np.trace(Rm)

reps, m_df = 300, 20
rng_df = np.random.default_rng(30041)
Y = mu + sigma * rng_df.normal(size=(reps, n))            # repeated responses, same design and mean
paths = [l2_boost(Xs, Y[r], 0.1, m_df) for r in range(reps)]
F = np.array([Xs @ P[-1] for P in paths])                 # fitted values after m_df steps
cov_df = np.sum(np.mean((Y - Y.mean(axis=0)) * (F - F.mean(axis=0)), axis=0)) / sigma ** 2
trace_df = np.mean([operator_trace(Xs, P, 0.1) for P in paths])
print(f"m = {m_df}: mean trace {trace_df:.1f}, covariance df {cov_df:.1f}")
# <</df_small>>
assert cov_df > trace_df

# the full least squares base learner: B_m = (1 - (1 - nu)^m) M
M_hat = Xs @ np.linalg.pinv(Xs)
fit = np.zeros(n)
for _ in range(7):
    fit += nu * M_hat @ (ys - fit)
assert np.allclose(fit, (1 - (1 - nu) ** 7) * M_hat @ ys)

gen = Generated("ch30", "boosting")
gen.int("n", n)
gen.int("p", p)
gen.num("sigma", sigma, 1)
gen.num("kappa", kappa, 2, sci=True)
gen.num("dist_100", dist_ls[100], 3)
gen.num("dist_1000", dist_ls[1000], 3)
gen.num("dist_50000", dist_ls[50000], 1, sci=True)
gen.num("err_ls", err_ls, 3)
gen.int("m_best", m_best)
gen.num("err_best", test_err[0.1][m_best], 3)
gen.int("m_best_nu1", int(np.argmin(test_err[1.0])))
gen.num("err_best_nu1", test_err[1.0].min(), 3)
gen.int("active_best", active_best)
gen.int("m_cv", m_cv)
gen.num("err_cv", err_cv, 3)
gen.num("err_bayes", sigma ** 2, 2)
gen.num("small_trace", trace_df, 1)
gen.num("small_cov", cov_df, 1)
for st, tag in [(0.1, "01"), (1.0, "1")]:
    gen.int(f"near_lo_{tag}", near_lo[st])
    gen.int(f"near_hi_{tag}", near_hi[st])
for c, m in enumerate(checkpoints):
    gen.num(f"df_trace_{m}", df_trace[c], 1)
    gen.num(f"df_cov_{m}", df_cov[c], 1)
gen.write()

use_book_style()
fig, axes = plt.subplots(1, 3, figsize=(6.2, 2.2))
ax = axes[0]
path = l2_boost(Xs, ys, 0.1, M)
for j in range(p):
    ax.plot(np.arange(M + 1), path[:, j], color=COLORS["accent"] if beta[j] != 0 else COLORS["muted"],
            linewidth=1.0 if beta[j] != 0 else 0.6)
ax.axvline(m_cv, color=COLORS["second"], linewidth=0.8, linestyle="--")
ax.set_xscale("symlog", linthresh=10)
ax.set_xlim(0, M)
ax.set_xlabel("iteration $m$")
ax.set_ylabel("coefficient")
ax.set_title(r"(a) path, $\nu=0.1$")
ax = axes[1]
ax.plot(np.arange(M + 1), test_err[0.1], color=COLORS["accent"], label=r"$\nu=0.1$")
ax.plot(np.arange(M + 1), test_err[1.0], color=COLORS["second"], label=r"$\nu=1$")
ax.axhline(err_ls, color=COLORS["muted"], linestyle="--", linewidth=0.8, label="least squares")
ax.set_xscale("symlog", linthresh=10)
ax.set_xlim(0, M)
ax.set_ylim(sigma ** 2 * 0.9, 8)
ax.set_xlabel("iteration $m$")
ax.set_title("(b) test error")
ax.legend(frameon=False, fontsize=7, loc="upper right")
ax = axes[2]
ax.plot(checkpoints, df_trace, "o-", color=COLORS["accent"], markersize=3, label="trace")
ax.plot(checkpoints, df_cov, "s-", color=COLORS["second"], markersize=3, label="covariance")
ax.set_xlabel("iteration $m$")
ax.set_title(r"(c) degrees of freedom")
ax.legend(frameon=False, fontsize=7, loc="lower right")
fig.tight_layout()
fig.savefig(figure_path("ch30", "boosting"))
