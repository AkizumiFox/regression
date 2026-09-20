"""Chapter 33, Section 4: covariance models for longitudinal data, fitted by REML.

Grunfeld's panel of 11 US firms observed in each of the 20 years 1935-1954
(public domain, shipped with statsmodels). The response is log gross investment,
the mean model an intercept, a linear trend and log market value.
"""
import matplotlib.pyplot as plt
import numpy as np
import statsmodels.api as sm
from scipy.optimize import minimize

from regbook import COLORS, Generated, figure_path, use_book_style

# <<data>>
df = sm.datasets.grunfeld.load_pandas().data
firms = sorted(df["firm"].unique())
years = np.sort(df["year"].unique())
N, m = len(firms), len(years)                        # 11 firms, 20 years
t = (years - years.mean()) / 10                      # centred and scaled time

Ylist, Xlist = [], []
for f in firms:
    d = df[df["firm"] == f].sort_values("year")
    Ylist.append(np.log(d["invest"].to_numpy()))
    Xlist.append(np.column_stack([np.ones(m), t, np.log(d["value"].to_numpy())]))
p = Xlist[0].shape[1]
n = N * m
# <</data>>
assert all(len(y) == m for y in Ylist) and n == 220


# <<reml>>
def reml(Sigma, Ylist, Xlist):
    """-2 x the restricted log likelihood, and the GLS estimate, for a common Sigma."""
    L = np.linalg.cholesky(Sigma)
    logdet = 2 * np.sum(np.log(np.diag(L)))
    solve = lambda B: np.linalg.solve(Sigma, B)
    XtWX = sum(X.T @ solve(X) for X in Xlist)
    XtWy = sum(X.T @ solve(y) for X, y in zip(Xlist, Ylist))
    beta = np.linalg.solve(XtWX, XtWy)
    quad = sum((y - X @ beta) @ solve(y - X @ beta) for X, y in zip(Xlist, Ylist))
    s, ld_info = np.linalg.slogdet(XtWX)
    n, p = len(Ylist) * Sigma.shape[0], len(beta)
    return len(Ylist) * logdet + quad + ld_info + (n - p) * np.log(2 * np.pi), beta


lag = np.abs(np.subtract.outer(np.arange(m), np.arange(m)))
Zrc = np.column_stack([np.ones(m), t])               # random intercept and slope basis


def build(name, th):
    """Sigma for each covariance model, from an unconstrained parameter vector."""
    if name == "independence":
        return np.exp(th[0]) * np.eye(m)
    if name == "compound symmetry":
        return np.exp(th[0]) * np.eye(m) + np.exp(th[1]) * np.ones((m, m))
    if name == "AR(1)":
        return np.exp(th[0]) * np.tanh(th[1])**lag
    if name == "AR(1) + random intercept":
        return (np.exp(th[0]) * np.tanh(th[1])**lag + np.exp(th[2]) * np.ones((m, m)))
    if name == "random intercept and slope":
        D = np.array([[np.exp(th[1]), th[3]], [th[3], np.exp(th[2])]])
        return np.exp(th[0]) * np.eye(m) + Zrc @ D @ Zrc.T
    if name == "Toeplitz(4)":
        a = np.concatenate([[1.0], th[1:4]])         # a moving average of order 3 ...
        g = np.array([a[:4 - u] @ a[u:] for u in range(4)])
        c = np.concatenate([g / g[0], np.zeros(m - 4)])   # ... so the band is nonneg. definite
        return np.exp(th[0]) * c[lag]
    raise ValueError(name)


MODELS = {"independence": 1, "compound symmetry": 2, "AR(1)": 2,
          "AR(1) + random intercept": 3, "random intercept and slope": 4,
          "Toeplitz(4)": 4}


def fit(name, start):
    def obj(th):
        S = build(name, th)
        if np.min(np.linalg.eigvalsh(S)) <= 1e-9:
            return 1e6
        return reml(S, Ylist, Xlist)[0]
    best = min((minimize(obj, s, method="Nelder-Mead",
                         options={"maxiter": 20000, "xatol": 1e-9, "fatol": 1e-9})
                for s in start), key=lambda r: r.fun)
    S = build(name, best.x)
    d2l, beta = reml(S, Ylist, Xlist)
    q = MODELS[name]
    return {"m2ll": d2l, "aic": d2l + 2 * q, "bic": d2l + q * np.log(n - p),
            "q": q, "Sigma": S, "beta": beta}


starts = {"independence": [[-1.0]],
          "compound symmetry": [[-1.0, -1.0]],
          "AR(1)": [[-1.0, 1.0]],
          "AR(1) + random intercept": [[-1.5, 1.0, -1.5]],
          "random intercept and slope": [[-2.0, -1.0, -2.0, 0.0], [-3.0, -2.0, -3.0, 0.0]],
          "Toeplitz(4)": [[-1.0, 1.5, 1.0, 0.5], [-1.0, 2.5, 2.0, 1.0]]}
res = {name: fit(name, starts[name]) for name in MODELS}
for name in MODELS:
    r = res[name]
    print(f"{name:28s} q={r['q']}  -2logL_R={r['m2ll']:9.3f}  "
          f"AIC={r['aic']:9.3f}  BIC={r['bic']:9.3f}")
# <</reml>>

best_name = min(MODELS, key=lambda k: res[k]["aic"])
assert best_name == min(MODELS, key=lambda k: res[k]["bic"])
assert res["independence"]["aic"] == max(r["aic"] for r in res.values())
# nesting: adding a parameter cannot decrease the restricted likelihood
assert res["AR(1) + random intercept"]["m2ll"] <= res["AR(1)"]["m2ll"] + 1e-6
assert res["AR(1)"]["m2ll"] <= res["independence"]["m2ll"] + 1e-6
assert res["compound symmetry"]["m2ll"] <= res["independence"]["m2ll"] + 1e-6
assert res["Toeplitz(4)"]["m2ll"] <= res["independence"]["m2ll"] + 1e-6
# an unstructured covariance would need this many parameters, from 11 firms
n_unstructured = m * (m + 1) // 2

# <<standard-errors>>
Xall = np.vstack(Xlist)
yall = np.concatenate(Ylist)
Sbest = res[best_name]["Sigma"]
A = np.linalg.inv(sum(X.T @ np.linalg.solve(Sbest, X) for X in Xlist))
beta = res[best_name]["beta"]
resid = [y - X @ beta for X, y in zip(Xlist, Ylist)]
meat = sum(np.outer(X.T @ np.linalg.solve(Sbest, r), X.T @ np.linalg.solve(Sbest, r))
           for X, r in zip(Xlist, resid))
cov_robust = A @ meat @ A * N / (N - 1)              # empirical ("sandwich") covariance
XtXinv = np.linalg.inv(Xall.T @ Xall)
beta_ols = XtXinv @ Xall.T @ yall
s2_ols = np.sum((yall - Xall @ beta_ols)**2) / (n - p)
for j, lab in enumerate(["intercept", "time/10", "log value"]):
    print(f"{lab:10s} GLS {beta[j]: .4f}  model se {np.sqrt(A[j, j]):.4f}  "
          f"empirical se {np.sqrt(cov_robust[j, j]):.4f}  "
          f"| OLS {beta_ols[j]: .4f} se {np.sqrt(s2_ols * XtXinv[j, j]):.4f}")
# <</standard-errors>>

# ---- checking the fitted covariance -----------------------------------------
R = np.array([y - X @ beta_ols for X, y in zip(Xlist, Ylist)])   # N x m residual array
emp_corr = np.array([np.mean([np.corrcoef(R[:, j], R[:, j + u])[0, 1]
                              for j in range(m - u)]) for u in range(m)])
variogram = np.array([np.mean([(R[:, j] - R[:, j + u])**2 for j in range(m - u)]) / 2
                      for u in range(1, m)])
total_var = np.mean([np.var(R[i]) for i in range(N)]) + np.var(R.mean(axis=1))


def model_corr(S):
    D = np.sqrt(np.diag(S))
    C = S / np.outer(D, D)
    return np.array([np.mean([C[j, j + u] for j in range(m - u)]) for u in range(m)])


assert emp_corr[0] == 1.0
assert emp_corr[1] > 0.4                              # strong correlation between adjacent years

gen = Generated("ch33", "covariance_models")
gen.int("N", N)
gen.int("m", m)
gen.int("n", n)
gen.int("p", p)
gen.int("n_unstructured", n_unstructured)
short = {"independence": "ind", "compound symmetry": "cs", "AR(1)": "ar1",
         "AR(1) + random intercept": "ar1ri", "random intercept and slope": "rcs",
         "Toeplitz(4)": "toep"}
for name, key in short.items():
    gen.num("m2ll_" + key, res[name]["m2ll"], 2)
    gen.num("aic_" + key, res[name]["aic"], 2)
    gen.num("bic_" + key, res[name]["bic"], 2)
    gen.int("q_" + key, res[name]["q"])
gen.text("best", best_name)
gen.num("delta_aic", res["compound symmetry"]["aic"] - res[best_name]["aic"], 2)
for j, key in enumerate(["b0", "b_time", "b_value"]):
    gen.num(key, beta[j], 4)
    gen.num("se_model_" + key, np.sqrt(A[j, j]), 4)
    gen.num("se_robust_" + key, np.sqrt(cov_robust[j, j]), 4)
    gen.num("ols_" + key, beta_ols[j], 4)
    gen.num("se_ols_" + key, np.sqrt(s2_ols * XtXinv[j, j]), 4)
for u in [1, 2, 5, 10]:
    gen.num(f"corr_lag{u}", emp_corr[u], 3)
gen.num("total_var", total_var, 4)
gen.num("variogram_plateau", float(np.mean(variogram[4:])), 3)
gen.write()

# ---- figure ------------------------------------------------------------------
use_book_style()
fig, axes = plt.subplots(1, 2, figsize=(5.6, 2.4))
ax = axes[0]
ax.plot(np.arange(1, m), variogram, "o", color=COLORS["ink"], markersize=3,
        label="empirical")
ax.axhline(total_var, color=COLORS["grid"], linewidth=0.8, zorder=0)
for name, col in [("compound symmetry", COLORS["second"]),
                  (best_name, COLORS["accent"])]:
    S = res[name]["Sigma"]
    g = np.array([np.mean([(S[j, j] + S[j + u, j + u]) / 2 - S[j, j + u]
                           for j in range(m - u)]) for u in range(1, m)])
    ax.plot(np.arange(1, m), g, color=col, label=name)
ax.set_xlabel("lag in years")
ax.set_ylabel("semi-variance")
ax.set_title("(a) variogram of residuals")
ax.set_ylim(0, 0.62)
ax.legend(fontsize=6.5, frameon=False, loc="upper left")
ax.text(18.5, total_var - 0.035, "total variance", ha="right", va="top",
        fontsize=6.5, color=COLORS["muted"])
ax = axes[1]
ax.plot(np.arange(m), emp_corr, "o", color=COLORS["ink"], markersize=3, label="empirical")
for name, col in [("compound symmetry", COLORS["second"]),
                  ("AR(1)", COLORS["third"]), (best_name, COLORS["accent"])]:
    ax.plot(np.arange(m), model_corr(res[name]["Sigma"]), color=col, label=name)
ax.axhline(0, color=COLORS["grid"], linewidth=0.6, zorder=0)
ax.set_xlabel("lag in years")
ax.set_ylabel("correlation")
ax.set_title("(b) average correlation by lag")
ax.set_ylim(0, 1.15)
ax.legend(fontsize=6.5, frameon=False, loc="lower left")
fig.tight_layout()
fig.savefig(figure_path("ch33", "covariance_models"))
