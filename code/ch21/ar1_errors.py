"""Chapter 21, Section 6: regression with first-order autoregressive errors.

(a) The quarterly Okun's-law regression of Section 21.5 (statsmodels.datasets.macrodata, public
domain): Prais-Winsten and iterated Cochrane-Orcutt estimates, and Newey-West standard errors.
(b) Exact efficiency of OLS relative to GLS, and the ratio of the true to the nominal OLS
variance, for n = 100 and AR(1) errors, with three regressors: a linear trend, a fixed
draw of white noise, and a fixed draw of an AR(1) series with coefficient 0.8.
(c) Coverage of nominal 95% intervals by simulation.
"""
import matplotlib.pyplot as plt
import numpy as np
import statsmodels.api as sm
from scipy import stats
from statsmodels.stats.sandwich_covariance import cov_hac

from regbook import COLORS, Generated, figure_path, use_book_style

macro = sm.datasets.macrodata.load_pandas().data
du = np.diff(macro["unemp"].to_numpy())
growth = 400 * np.diff(np.log(macro["realgdp"].to_numpy()))
n = len(du)
X = np.column_stack([np.ones(n), growth])


# <<pw>>
def prais_winsten(X, y, rho):
    """Whiten AR(1) errors: row 1 times sqrt(1 - rho^2), row t - rho row t-1."""
    Xs, ys = X.copy(), y.copy()
    Xs[1:], ys[1:] = X[1:] - rho * X[:-1], y[1:] - rho * y[:-1]
    Xs[0], ys[0] = np.sqrt(1 - rho**2) * X[0], np.sqrt(1 - rho**2) * y[0]
    return Xs, ys


def ols(X, y):
    b, *_ = np.linalg.lstsq(X, y, rcond=None)
    e = y - X @ b
    cov = (e @ e / (len(y) - X.shape[1])) * np.linalg.inv(X.T @ X)
    return b, e, cov


b_ols, e, V_ols = ols(X, du)
# regress e_t on e_(t-1)
rho = np.sum(e[1:] * e[:-1]) / np.sum(e[:-1] ** 2)
# two-step Prais-Winsten
b_pw, e_pw, V_pw = ols(*prais_winsten(X, du, rho))
print(f"rho-hat {rho:.3f}")
print(f"OLS slope {b_ols[1]:.4f} (se {np.sqrt(V_ols[1, 1]):.4f}); "
      f"Prais-Winsten {b_pw[1]:.4f} (se {np.sqrt(V_pw[1, 1]):.4f})")
# <</pw>>


# <<co>>
def cond_ss(b, r):
    """Conditional sum of squares S(beta, rho) over t = 2, ..., n."""
    u = (du[1:] - r * du[:-1]) - (X[1:] - r * X[:-1]) @ b
    return u @ u


b_co, r_co, path = b_ols, 0.0, []
for it in range(100):
    res = du - X @ b_co
    r_new = np.sum(res[1:] * res[:-1]) / np.sum(res[:-1] ** 2)   # rho step
    path.append(cond_ss(b_co, r_new))
    # beta step
    b_co, *_ = np.linalg.lstsq(X[1:] - r_new * X[:-1], du[1:] - r_new * du[:-1],
                               rcond=None)
    path.append(cond_ss(b_co, r_new))
    if abs(r_new - r_co) < 1e-10:
        break
    r_co = r_new
print(f"Cochrane-Orcutt after {it + 1} iterations: rho {r_co:.4f}, "
      f"slope {b_co[1]:.4f}")
# <</co>>
# each half-step lowers S
assert np.all(np.diff(path) <= 1e-12)


# <<nw>>
def newey_west(X, e, L):
    """Bartlett-weighted HAC covariance of the OLS estimate, truncation lag L."""
    v = X * e[:, None]
    S = v.T @ v
    for k in range(1, L + 1):
        G = v[k:].T @ v[:-k]
        S += (1 - k / (L + 1)) * (G + G.T)
    B = np.linalg.inv(X.T @ X)
    return B @ S @ B


V_nw = newey_west(X, e, 4)
print(f"Newey-West (L = 4) se of slope {np.sqrt(V_nw[1, 1]):.4f}")
# <</nw>>
fit = sm.OLS(du, X).fit()
assert np.allclose(V_nw, cov_hac(fit, nlags=4, use_correction=False))
assert np.min(np.linalg.eigvalsh(V_nw)) > 0

# <<dynamic>>
# the AR(1) error model is the dynamic regression
#   du_t = a + rho du_(t-1) + b g_t + c g_(t-1) + u_t
#   with the restriction c = -rho b
Xd = np.column_stack([np.ones(n - 1), du[:-1], growth[1:], growth[:-1]])
bd, ed, Vd = ols(Xd, du[1:])
# zero under the restriction
restr = bd[3] + bd[1] * bd[2]
grad = np.array([0.0, bd[2], bd[1], 1.0])                      # delta method
z_cf = restr / np.sqrt(grad @ Vd @ grad)
print("dynamic regression: const, du(t-1), g(t), g(t-1) =", np.round(bd, 4))
print(f"common-factor restriction c + rho b = {restr:.4f}, z = {z_cf:.1f}")
# <</dynamic>>
lagsd = np.column_stack([np.concatenate([np.zeros(k), ed[:-k]]) for k in range(1, 5)])
Zd = np.column_stack([Xd, lagsd])
cd, *_ = np.linalg.lstsq(Zd, ed, rcond=None)
bg_dyn = (n - 1) * (1 - np.sum((ed - Zd @ cd) ** 2) / np.sum((ed - ed.mean()) ** 2))

# residual autocorrelation after whitening
e_star = e_pw[1:]
r1_after = np.sum(e_star[1:] * e_star[:-1]) / np.sum(e_star**2)

gen = Generated("ch21", "ar1_errors", prefix="ar")
gen.num("rho", rho, 3)
gen.num("ols_slope", b_ols[1], 4)
gen.num("ols_se", np.sqrt(V_ols[1, 1]), 4)
gen.num("pw_slope", b_pw[1], 4)
gen.num("pw_se", np.sqrt(V_pw[1, 1]), 4)
gen.num("pw_int", b_pw[0], 3)
gen.num("co_rho", r_co, 3)
gen.num("co_slope", b_co[1], 4)
gen.int("co_iter", it + 1)
gen.num("nw_se", np.sqrt(V_nw[1, 1]), 4)
gen.num("r1_after", r1_after, 3)
for k, name in enumerate(("const", "lagdu", "g", "lagg")):
    gen.num(f"dyn_{name}", bd[k], 4)
    gen.num(f"dyn_se_{name}", np.sqrt(Vd[k, k]), 4)
gen.num("dyn_implied", -bd[1] * bd[2], 4)
gen.num("dyn_z", z_cf, 1)
gen.num("dyn_bg", bg_dyn, 2)
gen.num("dyn_bg_p", stats.chi2.sf(bg_dyn, 4), 2)
assert abs(z_cf) > 3 and stats.chi2.sf(bg_dyn, 4) > 0.1

# ---- (b) exact efficiency and variance ratios ----------------------------------------------
# <<efficiency>>
def ar1_cov(n, rho):
    t = np.arange(n)
    return rho ** np.abs(t[:, None] - t[None, :]) / (1 - rho**2)


def slope_variances(X, rho):
    """True OLS variance, GLS variance, and the expected nominal OLS variance, of the slope."""
    n, p = X.shape
    V = ar1_cov(n, rho)
    B = np.linalg.inv(X.T @ X)
    true_ols = (B @ X.T @ V @ X @ B)[1, 1]
    gls = np.linalg.inv(X.T @ np.linalg.solve(V, X))[1, 1]
    M = X @ B @ X.T
    es2 = np.trace((np.eye(n) - M) @ V) / (n - p)                 # E(s^2)
    return true_ols, gls, es2 * B[1, 1]


rng = np.random.default_rng(2109)
m = 100
x_ar = np.zeros(m)
for t in range(1, m):
    x_ar[t] = 0.8 * x_ar[t - 1] + rng.normal()
regressors = {"trend": np.arange(m, dtype=float), "white noise": rng.normal(size=m), "AR(0.8)": x_ar}
for name, x in regressors.items():
    true_ols, gls, nominal = slope_variances(np.column_stack([np.ones(m), x]), 0.6)
    print(f"{name:11s} rho = 0.6: efficiency of OLS {gls / true_ols:.3f}, "
          f"true / nominal OLS variance {true_ols / nominal:.2f}")
# <</efficiency>>

rhos = np.linspace(0, 0.95, 39)
curves = {name: np.array([slope_variances(np.column_stack([np.ones(m), x]), r) for r in rhos])
          for name, x in regressors.items()}
tags = {"trend": "trend", "white noise": "wn", "AR(0.8)": "arx"}
for name, x in regressors.items():
    for r in (0.6, 0.9):
        true_ols, gls, nominal = slope_variances(np.column_stack([np.ones(m), x]), r)
        gen.num(f"eff_{tags[name]}_{int(r * 10)}", gls / true_ols, 3)
        gen.num(f"infl_{tags[name]}_{int(r * 10)}", true_ols / nominal, 2)
        assert gls <= true_ols * (1 + 1e-12)
phi_hat = np.sum(x_ar[1:] * x_ar[:-1]) / np.sum(x_ar**2)
gen.num("phi_hat", phi_hat, 2)
gen.num("approx_infl_6", (1 + 0.6 * phi_hat) / (1 - 0.6 * phi_hat), 2)
gen.num("wn_limit_6", (1 - 0.36) / (1 + 0.36), 3)
gen.num("wn_limit_9", (1 - 0.81) / (1 + 0.81), 3)

# ---- (c) coverage by simulation --------------------------------------------------------------
# <<coverage>>
def ar1_draws(n, rho, reps, rng):
    u = rng.normal(size=(reps, n))
    x = np.empty((reps, n))
    x[:, 0] = u[:, 0] / np.sqrt(1 - rho**2)
    for t in range(1, n):
        x[:, t] = rho * x[:, t - 1] + u[:, t]
    return x


def coverage(n, reps, rng):
    """Coverage of 95% slope intervals: OLS, Newey-West, Prais-Winsten, GLS."""
    x = ar1_draws(n, 0.8, 1, rng)[0]
    Xc = np.column_stack([np.ones(n), x])
    L = int(4 * (n / 100) ** (2 / 9))
    q = stats.t.ppf(0.975, n - 2)
    hits = np.zeros(4)
    for y in ar1_draws(n, 0.6, reps, rng):          # true slope 0
        b, res, V = ols(Xc, y)
        r = np.sum(res[1:] * res[:-1]) / np.sum(res[:-1] ** 2)
        b_f, _, V_f = ols(*prais_winsten(Xc, y, r))
        b_g, _, V_g = ols(*prais_winsten(Xc, y, 0.6))
        ses = [V[1, 1], newey_west(Xc, res, L)[1, 1], V_f[1, 1], V_g[1, 1]]
        hits += [abs(bb[1]) <= q * np.sqrt(v)
                 for bb, v in zip((b, b, b_f, b_g), ses)]
    return hits / reps


for n_sim in (50, 200):
    print(n_sim, "OLS, Newey-West, Prais-Winsten, GLS:",
          np.round(coverage(n_sim, 500, rng), 3))
# <</coverage>>

cov = {n_sim: coverage(n_sim, 10_000, rng) for n_sim in (50, 200)}
for n_sim, c in cov.items():
    for k, name in enumerate(("ols", "nw", "pw", "gls")):
        gen.num(f"cov_{name}_{n_sim}", c[k], 3)
    assert c[0] < 0.85 and c[1] < c[2] and abs(c[3] - 0.95) < 0.01
gen.write()

# ---- figure -------------------------------------------------------------------------------------
use_book_style()
fig, axes = plt.subplots(1, 2, figsize=(5.8, 2.4))
colors = {"trend": COLORS["accent"], "white noise": COLORS["second"], "AR(0.8)": COLORS["third"]}
for name, c in curves.items():
    axes[0].plot(rhos, c[:, 1] / c[:, 0], color=colors[name], label=name)
    axes[1].plot(rhos, c[:, 0] / c[:, 2], color=colors[name], label=name)
axes[0].set_ylim(0, 1.03)
axes[0].set_xlabel(r"error autocorrelation $\rho$")
axes[0].set_ylabel("GLS variance / OLS variance")
axes[0].set_title("(a) efficiency of OLS")
axes[0].legend(frameon=False, loc="lower left")
axes[1].set_yscale("log")
axes[1].axhline(1, color=COLORS["grid"], linewidth=0.8, zorder=0)
axes[1].set_xlabel(r"error autocorrelation $\rho$")
axes[1].set_ylabel("true / nominal variance")
axes[1].set_title("(b) what the OLS standard error misses")
fig.tight_layout()
fig.savefig(figure_path("ch21", "ar1_efficiency"))
