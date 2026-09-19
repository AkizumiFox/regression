"""Chapter 19, Section 2: least squares when Cov(e) = sigma^2 V.

(a) AR(1) errors, v_st = rho^|s-t|, n = 20, a straight line in a trend regressor and in an
    alternating regressor. The ratio E(estimated variance)/(true variance) of the slope, with
    the eigenvalue bounds of Swindel; coverage of the nominal 95% t interval.
(b) Efficiency of ordinary relative to generalized least squares, against the Kantorovich
    bound 4 lam_1 lam_n / (lam_1 + lam_n)^2, and a design that attains the bound.
(c) Equicorrelated errors: slopes keep correct standard errors, the intercept does not.
"""
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats

from regbook import COLORS, Generated, figure_path, use_book_style

gen = Generated("ch19", "wrong_covariance", prefix="wc")

# <<ratio>>
import numpy as np
n = 20
t = np.arange(1, n + 1)
designs = {"trend": t.astype(float), "alternating": (-1.0) ** t}

def ar1(rho):
    return rho ** np.abs(t[:, None] - t[None, :])

def slope_summary(x, V):
    """True variance, mean of the usual estimate, GLS variance of the slope (sigma = 1)."""
    X = np.column_stack([np.ones(n), x])
    G = np.linalg.inv(X.T @ X)
    A = G @ X.T                                     # beta_hat = A y
    true_var = (A @ V @ A.T)[1, 1]                  # the sandwich
    M = X @ A
    Es2 = np.trace((np.eye(n) - M) @ V) / (n - 2)   # E(s^2)
    usual = Es2 * G[1, 1]                           # E of s^2 [(X'X)^{-1}]_{11}
    gls = np.linalg.inv(X.T @ np.linalg.solve(V, X))[1, 1]
    return true_var, usual, gls

rho = 0.6
for name, x in designs.items():
    true_var, usual, gls = slope_summary(x, ar1(rho))
    print(f"{name:12s} usual/true = {usual / true_var:.3f}   GLS/OLS variance = {gls / true_var:.3f}")
# <</ratio>>

def swindel_bounds(V, k):
    lam = np.sort(np.linalg.eigvalsh(V))
    return lam[:k].mean() / lam[-1], lam[-k:].mean() / lam[0]

res = {}
for name, x in designs.items():
    tv, us, gl = slope_summary(x, ar1(rho))
    res[name] = (tv, us, gl)
    lo, hi = swindel_bounds(ar1(rho), n - 2)
    assert lo <= us / tv <= hi
    assert gl <= tv + 1e-12                          # Aitken: GLS never worse
lam = np.linalg.eigvalsh(ar1(rho))
kant = 4 * lam[-1] * lam[0] / (lam[-1] + lam[0]) ** 2
for name in designs:
    assert res[name][2] / res[name][0] >= kant - 1e-12
lo, hi = swindel_bounds(ar1(rho), n - 2)
gen.num("rho", rho, 1)
gen.int("n", n)
gen.num("ratio_trend", res["trend"][1] / res["trend"][0], 3)
gen.num("ratio_alt", res["alternating"][1] / res["alternating"][0], 3)
gen.num("se_ratio_trend", np.sqrt(res["trend"][1] / res["trend"][0]), 3)
gen.num("se_ratio_alt", np.sqrt(res["alternating"][1] / res["alternating"][0]), 3)
gen.num("eff_trend", res["trend"][2] / res["trend"][0], 3)
gen.num("eff_alt", res["alternating"][2] / res["alternating"][0], 3)
gen.num("swindel_lo", lo, 3)
gen.num("swindel_hi", hi, 2)
gen.num("kant", kant, 4)
gen.num("lam_max", lam[-1], 3)
gen.num("lam_min", lam[0], 4)

# <<coverage>>
rng = np.random.default_rng(1902)
reps = 100_000
L = np.linalg.cholesky(ar1(rho))
E = rng.normal(size=(reps, n)) @ L.T              # rows are AR(1) error vectors
tq = stats.t.ppf(0.975, n - 2)
for name, x in designs.items():
    X = np.column_stack([np.ones(n), x])
    G = np.linalg.inv(X.T @ X)
    B = E @ (G @ X.T).T                             # beta_hat - beta, one row per data set
    s2 = np.sum((E - B @ X.T) ** 2, axis=1) / (n - 2)
    cover = np.mean(np.abs(B[:, 1]) <= tq * np.sqrt(s2 * G[1, 1]))
    print(f"{name:12s} coverage of the nominal 95% interval: {cover:.4f}")
# <</coverage>>
cov = {}
for name, x in designs.items():
    X = np.column_stack([np.ones(n), x])
    G = np.linalg.inv(X.T @ X)
    B = E @ (G @ X.T).T
    s2 = np.sum((E - B @ X.T) ** 2, axis=1) / (n - 2)
    cov[name] = np.mean(np.abs(B[:, 1]) <= tq * np.sqrt(s2 * G[1, 1]))
    # the simulated variance of the slope matches the sandwich
    assert abs(B[:, 1].var() / res[name][0] - 1) < 0.02
assert cov["trend"] < 0.85 and cov["alternating"] > 0.99
gen.num("cover_trend", cov["trend"], 3)
gen.num("cover_alt", cov["alternating"], 4)
gen.int("reps", reps)

# ---- (b) the Kantorovich bound is attained ---------------------------------------------
lam_all, Q = np.linalg.eigh(ar1(rho))
x_worst = (Q[:, 0] + Q[:, -1]) / np.sqrt(2)          # half on each extreme eigenvector
V = ar1(rho)
eff_worst = (x_worst @ x_worst) ** 2 / ((x_worst @ V @ x_worst) * (x_worst @ np.linalg.solve(V, x_worst)))
assert np.isclose(eff_worst, kant)

# ---- (c) equicorrelation: slope standard error right, intercept wrong --------------------
rho_e = 0.4
Ve = (1 - rho_e) * np.eye(n) + rho_e * np.ones((n, n))
X = np.column_stack([np.ones(n), t.astype(float), np.log(t)])
G = np.linalg.inv(X.T @ X)
A = G @ X.T
true_cov = A @ Ve @ A.T
Es2 = np.trace((np.eye(n) - X @ A) @ Ve) / (n - 3)
assert np.isclose(Es2, 1 - rho_e)
e1 = np.zeros(3); e1[0] = 1
assert np.allclose(true_cov, (1 - rho_e) * G + rho_e * np.outer(e1, e1))
M = X @ A
assert np.allclose(M @ Ve @ X, Ve @ X)               # Kruskal's condition C(VX) in C(X)
gen.write()

# ---- figure -----------------------------------------------------------------------------
use_book_style()
rhos = np.linspace(-0.9, 0.9, 73)
out = {k: [] for k in designs}
bl, bh, kb = [], [], []
for rh in rhos:
    V = ar1(rh)
    for name, x in designs.items():
        tv, us, gl = slope_summary(x, V)
        out[name].append((np.sqrt(us / tv), gl / tv))
    lo_, hi_ = swindel_bounds(V, n - 2)
    bl.append(np.sqrt(lo_)); bh.append(np.sqrt(hi_))
    lm = np.linalg.eigvalsh(V)
    kb.append(4 * lm[-1] * lm[0] / (lm[-1] + lm[0]) ** 2)
fig, axes = plt.subplots(1, 2, figsize=(5.8, 2.5))
ax = axes[0]
ax.fill_between(rhos, bl, bh, color=COLORS["grid"], alpha=0.6, lw=0, label="eigenvalue bounds")
ax.plot(rhos, [v[0] for v in out["trend"]], color=COLORS["accent"], label="trend")
ax.plot(rhos, [v[0] for v in out["alternating"]], color=COLORS["second"], label="alternating")
ax.axhline(1, color=COLORS["muted"], lw=0.6, ls=":")
ax.set_yscale("log")
ax.set_ylim(0.08, 12)
ax.set_xlabel(r"autocorrelation $\rho$")
ax.set_ylabel(r"$\sqrt{\mathrm{E}(v)\,/\,\mathrm{Var}(\hat\beta_1)}$")
ax.set_title("(a) the usual standard error")
ax.legend(frameon=False, loc="upper center")
ax = axes[1]
ax.plot(rhos, [v[1] for v in out["trend"]], color=COLORS["accent"], label="trend")
ax.plot(rhos, [v[1] for v in out["alternating"]], color=COLORS["second"], label="alternating")
ax.plot(rhos, kb, color=COLORS["muted"], ls="--", label="Kantorovich bound")
ax.set_yscale("log")
ax.set_xlabel(r"autocorrelation $\rho$")
ax.set_ylabel("Var(GLS) / Var(OLS)")
ax.set_title("(b) efficiency of least squares")
ax.legend(frameon=False, loc="lower center")
fig.tight_layout()
fig.savefig(figure_path("ch19", "wrong_covariance"))
