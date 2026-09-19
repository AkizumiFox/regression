"""Chapter 27, Section 2: ridge regression on Longley's employment data.

Regressors centred and scaled to unit standard deviation (divisor n - 1); the response, total
employment in thousands, centred. The intercept is not penalized. The script computes ridge
through the SVD, checks it against a direct solve, the augmented-data least squares form and the
conjugate posterior mean, computes the trace, effective degrees of freedom, leave-one-out and
generalized cross-validation, and the exact mean squared error of ridge in the "plug-in world"
where the true coefficients are the least squares estimates and sigma^2 = s^2.
Public-domain data shipped with statsmodels (statsmodels.datasets.longley).
"""
import matplotlib.pyplot as plt
import numpy as np
import statsmodels.api as sm
from scipy import optimize

from regbook import COLORS, Generated, figure_path, use_book_style

# <<setup>>
data = sm.datasets.longley.load_pandas().data
names = ["GNPDEFL", "GNP", "UNEMP", "ARMED", "POP", "YEAR"]
Z = data[names].to_numpy()
y = data["TOTEMP"].to_numpy()
n, p = Z.shape

X = (Z - Z.mean(axis=0)) / Z.std(axis=0, ddof=1)     # centred, unit standard deviation
yc = y - y.mean()                                     # centred response (the intercept is ybar)

U, d, Vt = np.linalg.svd(X, full_matrices=False)      # X = U diag(d) V^T
V = Vt.T
c = U.T @ yc                                          # coordinates of y along u_1, ..., u_p


def ridge(lam):
    """Ridge coefficients for the centred problem, through the SVD."""
    return V @ (d / (d ** 2 + lam) * c)


def edf(lam):
    """Effective degrees of freedom of the slopes: sum of the shrinkage factors."""
    return np.sum(d ** 2 / (d ** 2 + lam))


beta_ls = V @ (c / d)
print("singular values squared:", np.round(d ** 2, 4))
print("least squares:", np.round(beta_ls, 1))
for lam in [0.01, 0.1, 1.0]:
    print(f"lambda = {lam:5.2f}  edf = {edf(lam):.2f}  ridge =", np.round(ridge(lam), 1))
# <</setup>>

# ---- checks of the closed form ----------------------------------------------------------------
lam0 = 0.05
direct = np.linalg.solve(X.T @ X + lam0 * np.eye(p), X.T @ yc)
assert np.allclose(ridge(lam0), direct)
aug = np.linalg.lstsq(np.vstack([X, np.sqrt(lam0) * np.eye(p)]),
                      np.concatenate([yc, np.zeros(p)]), rcond=None)[0]
assert np.allclose(aug, direct)                         # ridge is least squares on augmented data
assert np.allclose(ridge(1e-12), beta_ls, rtol=1e-6)    # lambda -> 0 gives least squares
norms = [np.linalg.norm(ridge(l)) for l in np.logspace(-4, 3, 60)]
assert np.all(np.diff(norms) < 0)                       # the length of the estimate decreases

# intercept unpenalized: fitting [1, X] with penalty diag(0, lam I) gives ybar and the centred ridge
W = np.column_stack([np.ones(n), X])
Omega = np.diag(np.r_[0.0, lam0 * np.ones(p)])
full = np.linalg.solve(W.T @ W + Omega, W.T @ y)
assert np.isclose(full[0], y.mean()) and np.allclose(full[1:], direct)

# Bayesian reading: posterior mean under beta ~ N(0, (sigma^2/lam) I), NIG form of Chapter 7
V0 = np.eye(p) / lam0
Vn = np.linalg.inv(np.linalg.inv(V0) + X.T @ X)
assert np.allclose(Vn @ (X.T @ yc), direct)

# ---- selection of lambda: leave-one-out and GCV ---------------------------------------------
# <<cv>>
def hat_diag(lam):
    """Leverages of the ridge smoother with an unpenalized intercept."""
    f = d ** 2 / (d ** 2 + lam)
    return 1 / n + np.sum(U ** 2 * f, axis=1)


def loo_and_gcv(lam):
    resid = yc - X @ ridge(lam)
    h = hat_diag(lam)
    loo = np.mean((resid / (1 - h)) ** 2)                  # leave-one-out shortcut
    gcv = np.mean(resid ** 2) / (1 - (1 + edf(lam)) / n) ** 2
    return loo, gcv


grid = np.logspace(-4, 2, 601)
scores = np.array([loo_and_gcv(l) for l in grid])
lam_loo = grid[np.argmin(scores[:, 0])]
lam_gcv = grid[np.argmin(scores[:, 1])]
print(f"LOO chooses lambda = {lam_loo:.4f} (edf {edf(lam_loo):.2f})")
print(f"GCV chooses lambda = {lam_gcv:.4f} (edf {edf(lam_gcv):.2f})")
print("ridge at the LOO choice:", np.round(ridge(lam_loo), 1))
# <</cv>>

# brute-force leave-one-out agrees with the shortcut
lam1 = 0.02
errs = []
for i in range(n):
    keep = np.arange(n) != i
    Wi = W[keep]
    bi = np.linalg.solve(Wi.T @ Wi + np.diag(np.r_[0.0, lam1 * np.ones(p)]), Wi.T @ y[keep])
    errs.append(y[i] - W[i] @ bi)
assert np.isclose(np.mean(np.square(errs)), loo_and_gcv(lam1)[0])

# the GNP coefficient changes sign along the ridge trace
lam_sign = optimize.brentq(lambda l: ridge(l)[1], 1e-5, 1.0)
assert beta_ls[1] < 0 and ridge(lam_sign / 2)[1] < 0 and ridge(2 * lam_sign)[1] > 0
print(f"GNP coefficient changes sign at lambda = {lam_sign:.4f} (edf {edf(lam_sign):.2f})")


# honest assessment of the tuned procedure: nested leave-one-out, the choice of lambda by the
# leave-one-out shortcut redone inside each fold (training data re-standardized in each fold)
def ridge_tuned_predict(Zt, yt, z0):
    nt = len(yt)
    mt, st = Zt.mean(axis=0), Zt.std(axis=0, ddof=1)
    Xt = (Zt - mt) / st
    Ut, dt, Vtt = np.linalg.svd(Xt, full_matrices=False)
    ct = Ut.T @ (yt - yt.mean())

    def loo_t(lam):
        f = dt ** 2 / (dt ** 2 + lam)
        r = yt - yt.mean() - Ut @ (f * ct)
        return np.mean((r / (1 - 1 / nt - np.sum(Ut ** 2 * f, axis=1))) ** 2)

    lam = grid[np.argmin([loo_t(l) for l in grid])]
    bt = Vtt.T @ (dt / (dt ** 2 + lam) * ct)
    return yt.mean() + ((z0 - mt) / st) @ bt


nested = [y[i] - ridge_tuned_predict(Z[np.arange(n) != i], y[np.arange(n) != i], Z[i]) for i in range(n)]
rmse_nested = np.sqrt(np.mean(np.square(nested)))
rmse_min = np.sqrt(scores[:, 0].min())
print(f"ridge: minimized LOO root mean squared error {rmse_min:.0f}, nested {rmse_nested:.0f}")
assert rmse_nested > rmse_min

# ---- mean squared error in the plug-in world --------------------------------------------------
# <<mse>>
resid_ls = yc - X @ beta_ls
s2 = resid_ls @ resid_ls / (n - p - 1)                    # residual df n - p - 1 (intercept)
alpha = Vt @ beta_ls                                      # canonical coefficients


def ridge_mse(lam, alpha=alpha, sigma2=s2):
    """Exact total MSE of ridge: variance and squared bias, summed over components."""
    var = sigma2 * np.sum(d ** 2 / (d ** 2 + lam) ** 2)
    bias2 = np.sum((lam * alpha / (d ** 2 + lam)) ** 2)
    return var, bias2


mse_grid = np.logspace(-5, 1, 1201)
total = np.array([sum(ridge_mse(l)) for l in mse_grid])
lam_star = mse_grid[np.argmin(total)]
mse_ls = s2 * np.sum(1 / d ** 2)
hk_bound = 2 * s2 / np.max(alpha ** 2)                    # Hoerl-Kennard: better for lam below this
print(f"MSE of least squares {mse_ls:.0f}, best ridge {total.min():.0f} at lambda {lam_star:.4f}")
print(f"guaranteed improvement for 0 < lambda < {hk_bound:.5f}")
# <</mse>>

assert np.isclose(sum(ridge_mse(0.0)), mse_ls)
inside = np.linspace(hk_bound / 1000, hk_bound * 0.999, 200)
assert all(sum(ridge_mse(l)) < mse_ls for l in inside)
# the sufficient bound is conservative: the improvement goes on well beyond it
lam_break = mse_grid[np.nonzero(total < mse_ls)[0].max()]
print(f"ridge beats least squares for lambda up to {lam_break:.4f}")
assert lam_break > hk_bound

# Theobald's matrix comparison: MSE(LS) - MSE(ridge) nnd for lam <= 2 sigma^2/||beta||^2
S = X.T @ X
theo = 2 * s2 / (beta_ls @ beta_ls)
for lam in [theo * 0.5, theo]:
    A = np.linalg.solve(S + lam * np.eye(p), S)
    bias = A @ beta_ls - beta_ls
    diff = s2 * np.linalg.inv(S) - s2 * A @ np.linalg.inv(S) @ A.T - np.outer(bias, bias)
    assert np.linalg.eigvalsh(diff).min() > -1e-8 * np.abs(diff).max()

# ---- figure: ridge trace and the MSE decomposition -------------------------------------------
use_book_style()
fig, axes = plt.subplots(1, 2, figsize=(6.2, 2.7))
ax = axes[0]
lams = np.logspace(-5, 3, 400)
path = np.array([ridge(l) for l in lams])
dfs = np.array([edf(l) for l in lams])
cols = [COLORS["accent"], COLORS["second"], COLORS["third"], COLORS["thread"], COLORS["muted"], COLORS["ink"]]
for j, name in enumerate(names):
    ax.plot(dfs, path[:, j] / 1000, color=cols[j], label=name)
ax.axvline(edf(lam_loo), color=COLORS["grid"], linewidth=0.8, zorder=0)
ax.axhline(0, color=COLORS["grid"], linewidth=0.6, zorder=0)
ax.set_xlabel("effective degrees of freedom")
ax.set_ylabel("coefficient (millions per SD)")
ax.set_title("(a) ridge trace")
ax.legend(frameon=False, fontsize=6.5, ncol=2, loc="lower left")
ax = axes[1]
parts = np.array([ridge_mse(l) for l in mse_grid])
ax.loglog(mse_grid, parts[:, 0] / 1e6, color=COLORS["accent"], label="variance")
ax.loglog(mse_grid, parts[:, 1] / 1e6, color=COLORS["second"], label="squared bias")
ax.loglog(mse_grid, total / 1e6, color=COLORS["ink"], label="total")
ax.axhline(mse_ls / 1e6, color=COLORS["muted"], linestyle="--", linewidth=0.8, label="least squares")
ax.set_ylim(1e-3, 30)
ax.set_xlabel(r"$\lambda$")
ax.set_ylabel(r"MSE (millions)")
ax.set_title("(b) MSE if the truth were the LS fit")
ax.legend(frameon=False, fontsize=6.5, loc="lower left")
fig.tight_layout()
fig.savefig(figure_path("ch27", "ridge_longley"))

gen = Generated("ch27", "ridge_longley")
gen.int("n", n)
gen.num("d2_max", d[0] ** 2, 2)
gen.num("d2_min", d[-1] ** 2, 4)
for j, name in enumerate(names):
    gen.num(f"ls_{name}", beta_ls[j], 0)
    gen.num(f"loo_{name}", ridge(lam_loo)[j], 0)
gen.num("lam_loo", lam_loo, 4)
gen.num("lam_gcv", lam_gcv, 4)
gen.num("edf_loo", edf(lam_loo), 2)
gen.num("edf_gcv", edf(lam_gcv), 2)
gen.num("s", np.sqrt(s2), 1)
gen.num("mse_ls", mse_ls / 1e6, 2)
gen.num("mse_best", total.min() / 1e6, 2)
gen.num("mse_ratio", total.min() / mse_ls, 3)
gen.num("lam_star", lam_star, 4)
gen.num("edf_star", edf(lam_star), 2)
gen.num("hk_bound", hk_bound, 5)
gen.num("theobald", theo, 5)
gen.num("lam_break", lam_break, 4)
gen.num("lam_sign", lam_sign, 4)
gen.num("edf_sign", edf(lam_sign), 1)
gen.num("rmse_min", rmse_min, 0)
gen.num("rmse_nested", rmse_nested, 0)
gen.write()
