"""Chapter 27, Section 4: the lasso by coordinate descent, checked three ways, on Longley's data.

The criterion is (1/2)||y - X b||^2 + lam ||b||_1 with the regressors centred and scaled to unit
standard deviation and the response centred (intercept unpenalized). Coordinate descent is checked
against the KKT conditions, against a brute-force search over active sets and sign patterns, and
against a bound-constrained quasi-Newton solver on the split variables b = u - v.
"""
import itertools

import matplotlib.pyplot as plt
import numpy as np
import statsmodels.api as sm
from scipy.optimize import minimize

from regbook import COLORS, Generated, figure_path, use_book_style

# <<cd>>
def soft(z, t):
    """Soft thresholding: sign(z) * max(|z| - t, 0)."""
    return np.sign(z) * np.maximum(np.abs(z) - t, 0.0)


def lasso_cd(X, y, lam, b=None, tol=1e-10, max_sweeps=100_000):
    """Cyclic coordinate descent for (1/2)||y - X b||^2 + lam ||b||_1."""
    p = X.shape[1]
    b = np.zeros(p) if b is None else b.copy()
    r = y - X @ b                                   # current residual
    sq = np.sum(X ** 2, axis=0)                     # ||x_j||^2
    for sweep in range(1, max_sweeps + 1):
        biggest = 0.0
        for j in range(p):
            z = X[:, j] @ r + sq[j] * b[j]          # x_j^T (partial residual without j)
            new = soft(z, lam) / sq[j]
            if new != b[j]:
                r -= X[:, j] * (new - b[j])
                biggest = max(biggest, abs(new - b[j]) * np.sqrt(sq[j]))
                b[j] = new
        if biggest < tol:
            break
    return b, sweep


def kkt_violation(X, y, b, lam):
    """Largest violation of the KKT conditions X^T (y - X b) = lam * s, s in the subdifferential."""
    g = X.T @ (y - X @ b)
    active = b != 0
    return max(np.max(np.abs(g[active] - lam * np.sign(b[active])), initial=0.0),
               np.max(np.abs(g[~active]) - lam, initial=0.0))
# <</cd>>


def objective(X, y, b, lam):
    return 0.5 * np.sum((y - X @ b) ** 2) + lam * np.sum(np.abs(b))


def brute_force(X, y, lam):
    """Try every active set A and sign pattern s: solve the KKT equations on A, keep consistent ones."""
    p = X.shape[1]
    best, best_obj = np.zeros(p), objective(X, y, np.zeros(p), lam)
    for k in range(1, p + 1):
        for A in itertools.combinations(range(p), k):
            XA = X[:, A]
            G = XA.T @ XA
            for s in itertools.product([-1.0, 1.0], repeat=k):
                bA = np.linalg.solve(G, XA.T @ y - lam * np.array(s))
                if np.all(np.sign(bA) == s):
                    b = np.zeros(p)
                    b[list(A)] = bA
                    obj = objective(X, y, b, lam)
                    if obj < best_obj:
                        best, best_obj = b, obj
    return best


def split_solver(X, y, lam):
    """L-BFGS-B on b = u - v with u, v >= 0, a smooth bound-constrained reformulation."""
    p = X.shape[1]

    def f(w):
        u, v = w[:p], w[p:]
        r = y - X @ (u - v)
        g = -X.T @ r
        return 0.5 * r @ r + lam * np.sum(w), np.r_[g + lam, -g + lam]

    res = minimize(f, np.zeros(2 * p), jac=True, method="L-BFGS-B", bounds=[(0, None)] * (2 * p),
                   options={"ftol": 1e-15, "gtol": 1e-12, "maxiter": 100_000, "maxfun": 100_000})
    return res.x[:p] - res.x[p:]


# <<longley>>
data = sm.datasets.longley.load_pandas().data
names = ["GNPDEFL", "GNP", "UNEMP", "ARMED", "POP", "YEAR"]
Z = data[names].to_numpy()
y = data["TOTEMP"].to_numpy()
n, p = Z.shape
X = (Z - Z.mean(axis=0)) / Z.std(axis=0, ddof=1)
yc = y - y.mean()

lam_max = np.max(np.abs(X.T @ yc))                  # smallest lam with b = 0
lams = lam_max * np.logspace(0, -3, 61)
path, b = [], np.zeros(p)
for lam in lams:                                     # warm starts along a decreasing grid
    b, _ = lasso_cd(X, yc, lam, b, tol=1e-6)
    path.append(b.copy())
path = np.array(path)
first = {names[j]: np.argmax(path[:, j] != 0) for j in range(p) if np.any(path[:, j] != 0)}
entry = sorted(first, key=first.get)                  # variables in order of entry
print(f"lambda_max = {lam_max:.0f}; order of entry:", entry)
print("GNP coefficient along the path:", np.round(path[::6, 1], 0))
# <</longley>>

gnp_peak = np.argmax(path[:, 1])
assert path[gnp_peak, 1] > 3000 and path[-1, 1] == 0      # GNP enters first, later leaves

# ---- three checks at several values of lambda -------------------------------------------------
checks = []
for frac in [0.5, 0.1, 0.02, 0.005, 0.001]:
    lam = frac * lam_max
    b_cd, sweeps = lasso_cd(X, yc, lam)
    b_bf = brute_force(X, yc, lam)
    b_sp = split_solver(X, yc, lam)
    kkt = kkt_violation(X, yc, b_cd, lam)
    assert kkt < 1e-6 * lam
    assert np.allclose(b_cd, b_bf, atol=1e-6 * np.abs(b_bf).max())
    rel = (objective(X, yc, b_sp, lam) - objective(X, yc, b_cd, lam)) / objective(X, yc, b_cd, lam)
    assert rel > -1e-12 and rel < 1e-7              # quasi-Newton is never better, and close
    checks.append((frac, sweeps, np.count_nonzero(b_cd)))
    print(f"lam = {frac:5.3f} lam_max: {sweeps:6d} sweeps, {np.count_nonzero(b_cd)} nonzero, "
          f"KKT violation {kkt:.1e}, relative objective gap of L-BFGS-B {rel:.1e}")

# lam >= lam_max gives zero; the path is piecewise linear between changes of the active set
assert np.all(lasso_cd(X, yc, lam_max * 1.0001)[0] == 0)
lam_a, lam_b = 0.30 * lam_max, 0.25 * lam_max
ba, bb = lasso_cd(X, yc, lam_a)[0], lasso_cd(X, yc, lam_b)[0]
bm = lasso_cd(X, yc, 0.5 * (lam_a + lam_b))[0]
if np.array_equal(np.sign(ba), np.sign(bb)):
    assert np.allclose(bm, 0.5 * (ba + bb), atol=1e-6)

# leave-one-out choice of lambda (refit with the case removed, recentring the training data)
# <<loo>>
grid = lam_max * np.logspace(-1, -4, 31)
loo = np.zeros(len(grid))
for i in range(n):
    keep = np.arange(n) != i
    Xi, yi = X[keep], y[keep]
    xm, ym = Xi.mean(axis=0), yi.mean()
    b = np.zeros(p)
    for g, lam in enumerate(grid):
        b, _ = lasso_cd(Xi - xm, yi - ym, lam, b, tol=1e-8)
        loo[g] += (y[i] - ym - (X[i] - xm) @ b) ** 2 / n
lam_loo = grid[np.argmin(loo)]
b_loo, _ = lasso_cd(X, yc, lam_loo)
print(f"LOO chooses lambda = {lam_loo:.0f} ({lam_loo / lam_max:.4f} lam_max), "
      f"root mean squared error {np.sqrt(loo.min()):.0f}; coefficients", np.round(b_loo, 0))
# <</loo>>



# honest assessment of the tuned lasso: nested leave-one-out. Within each outer fold the fifteen
# training years are re-standardized and lambda is chosen by their own leave-one-out error.
# A looser tolerance keeps the 16 x 15 path fits affordable (a few minutes).
def lasso_tuned_predict(Zt, yt, z0, tol=1e-3):
    nt = len(yt)
    mt, st = Zt.mean(axis=0), Zt.std(axis=0, ddof=1)
    Xt = (Zt - mt) / st
    ytc = yt - yt.mean()
    grid_t = np.max(np.abs(Xt.T @ ytc)) * np.logspace(-1, -4, 31)
    err = np.zeros(len(grid_t))
    for i in range(nt):
        keep = np.arange(nt) != i
        Xi, yi = Xt[keep], yt[keep]
        xm, ym = Xi.mean(axis=0), yi.mean()
        b = np.zeros(p)
        for g, lam in enumerate(grid_t):
            b, _ = lasso_cd(Xi - xm, yi - ym, lam, b, tol=tol)
            err[g] += (yt[i] - ym - (Xt[i] - xm) @ b) ** 2
    bt, _ = lasso_cd(Xt, ytc, grid_t[np.argmin(err)], tol=tol)
    return yt.mean() + ((z0 - mt) / st) @ bt


nested = [y[i] - lasso_tuned_predict(Z[np.arange(n) != i], y[np.arange(n) != i], Z[i]) for i in range(n)]
rmse_nested = np.sqrt(np.mean(np.square(nested)))
print(f"lasso: minimized LOO root mean squared error {np.sqrt(loo.min()):.0f}, nested {rmse_nested:.0f}")
assert rmse_nested > np.sqrt(loo.min())

use_book_style()
fig, axes = plt.subplots(1, 2, figsize=(6.2, 2.6))
ax = axes[0]
z = np.linspace(-3, 3, 601)
ax.plot(z, z, color=COLORS["grid"], linewidth=0.8)
ax.plot(z, soft(z, 1.0), color=COLORS["accent"], label="lasso: soft")
ax.plot(z, np.where(np.abs(z) > 1.0, z, 0.0), color=COLORS["second"], label="subset: hard")
ax.plot(z, z / 2, color=COLORS["third"], label="ridge: proportional")
ax.set_xlabel(r"least squares coefficient $\hat\beta_j$")
ax.set_ylabel("estimate")
ax.set_title("(a) orthonormal design")
ax.legend(frameon=False, fontsize=6.5, loc="upper left")
ax = axes[1]
ls = np.linalg.lstsq(X, yc, rcond=None)[0]
shrink = np.sum(np.abs(path), axis=1) / np.sum(np.abs(ls))
cols = [COLORS["accent"], COLORS["second"], COLORS["third"], COLORS["thread"], COLORS["muted"], COLORS["ink"]]
for j, name in enumerate(names):
    ax.plot(shrink, path[:, j] / 1000, color=cols[j], label=name)
ax.axvline(np.sum(np.abs(b_loo)) / np.sum(np.abs(ls)), color=COLORS["grid"], linewidth=0.8, zorder=0)
ax.set_xlabel(r"$\|\hat{b}_\lambda\|_1 / \|\hat\beta\|_1$")
ax.set_ylabel("coefficient (millions per SD)")
ax.set_title("(b) lasso path, Longley")
ax.legend(frameon=False, fontsize=6.5, ncol=2, loc="upper left")
fig.tight_layout()
fig.savefig(figure_path("ch27", "lasso_path"))

gen = Generated("ch27", "lasso_cd")
gen.num("lam_max", lam_max, 0)
gen.text("entry", ", ".join(entry))
gen.num("lam_loo", lam_loo, 0)
gen.num("lam_loo_frac", lam_loo / lam_max, 4)
gen.num("rmse_loo", np.sqrt(loo.min()), 0)
gen.num("rmse_nested", rmse_nested, 0)
gen.int("nonzero_loo", np.count_nonzero(b_loo))
for j, name in enumerate(names):
    gen.num(f"loo_{name}", b_loo[j], 0)
for frac, sweeps, nz in checks:
    gen.int(f"sweeps_{frac}", sweeps)
gen.write()
