"""Chapter 40, Section 4: generalized estimating equations.

A compact GEE for balanced clusters, applied to the travel-mode choice data (public domain,
statsmodels.datasets.modechoice: 210 travellers, each with one record per transport mode), and
a simulation of the coverage of the model-based, sandwich and bias-corrected covariances when
the working correlation is wrong.
"""
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import statsmodels.api as sm
from numpy.polynomial.hermite_e import hermegauss
from scipy.optimize import root
from scipy.stats import norm, t

from regbook import COLORS, Generated, figure_path, use_book_style


# <<gee>>
def expit(z):
    return 1.0 / (1.0 + np.exp(-z))


def working_correlation(kind, r, p, scale_fixed=False):
    """Moment estimate of the working correlation from the Pearson residuals r (m by n)."""
    m, n = r.shape
    phi = 1.0 if scale_fixed else np.sum(r**2) / (m * n - p)
    if kind == "independence":
        return np.eye(n), 0.0, phi
    if kind == "exchangeable":
        off = (np.sum(r.sum(1) ** 2) - np.sum(r**2)) / 2
        a = off / (phi * (m * n * (n - 1) / 2 - p))
        return np.eye(n) + a * (1 - np.eye(n)), a, phi
    if kind == "ar1":
        a = np.sum(r[:, :-1] * r[:, 1:]) / (phi * (m * (n - 1) - p))
        lag = np.abs(np.arange(n)[:, None] - np.arange(n))
        return a**lag, a, phi
    if kind == "unstructured":
        R = (r.T @ r) / (phi * m)
        R = R / np.sqrt(np.outer(np.diag(R), np.diag(R)))
        return R, np.nan, phi
    raise ValueError(kind)


def gee(y, X, kind="independence", family="binomial", scale_fixed=False, maxit=200):
    """Solve sum_i D_i' V_i^{-1} (y_i - mu_i) = 0 for balanced clusters y (m, n), X (m, n, p).

    With S_i = diag(mu'_ij / sqrt(V(mu_ij))) and Pearson residuals r_i, the Fisher-scoring
    step is beta += (sum_i X_i'S_i R^{-1} S_i X_i)^{-1} sum_i X_i'S_i R^{-1} r_i."""
    m, n, p = X.shape
    beta, R, alpha, phi = np.zeros(p), np.eye(n), 0.0, 1.0
    for _ in range(maxit):
        eta = X @ beta
        if family == "binomial":
            mu = expit(eta)
            dmu, V = mu * (1 - mu), mu * (1 - mu)
        else:                                            # log link, V(mu) = mu^2 (gamma)
            mu = np.exp(eta)
            dmu, V = mu, mu**2
        r = (y - mu) / np.sqrt(V)
        R, alpha, phi = working_correlation(kind, r, p, scale_fixed)
        Rinv = np.linalg.inv(R)
        S = dmu / np.sqrt(V)
        SX = S[..., None] * X                            # (m, n, p)
        A = np.einsum("mjp,jk,mkq->pq", SX, Rinv, SX)
        g = np.einsum("mjp,jk,mk->p", SX, Rinv, r)
        step = np.linalg.solve(A, g)
        beta = beta + step
        if np.max(np.abs(step)) < 1e-11:
            break
    u = np.einsum("mjp,jk,mk->mp", SX, Rinv, r)          # cluster contributions
    Ainv = np.linalg.inv(A)
    robust = Ainv @ (u.T @ u) @ Ainv                     # the sandwich covariance
    naive = phi * Ainv                                   # the model-based covariance
    return dict(beta=beta, naive=naive, robust=robust, alpha=alpha, phi=phi, R=R,
                A=A, Ainv=Ainv, S=S, r=r, Rinv=Rinv, X=X)
# <</gee>>


# <<mancl>>
def mancl_derouen(fit):
    """Bias-corrected sandwich: inflate each cluster's residual by (I - H_i)^{-1}."""
    SX, Rinv, Ainv, r = fit["S"][..., None] * fit["X"], fit["Rinv"], fit["Ainv"], fit["r"]
    n = r.shape[1]
    H = np.einsum("mjp,pq,mkq,kl->mjl", SX, Ainv, SX, Rinv)   # (m, n, n)
    adj = np.linalg.solve(np.eye(n) - H, r[..., None])[..., 0]
    u = np.einsum("mjp,jk,mk->mp", SX, Rinv, adj)
    return Ainv @ (u.T @ u) @ Ainv
# <</mancl>>


# ---- the travel-mode data ----------------------------------------------------
# <<data>>
d = sm.datasets.modechoice.load_pandas().data.sort_values(["individual", "mode"])
m, n = d["individual"].nunique(), 4                       # 210 travellers, 4 modes each
y = d["choice"].to_numpy().reshape(m, n)
D = pd.get_dummies(d["mode"].astype(int), drop_first=True).to_numpy(dtype=float)
cols = [np.ones(m * n)] + [D[:, j] for j in range(3)] + [d["ttme"] / 10, d["gc"] / 10]
X = np.stack([np.asarray(c, dtype=float).reshape(m, n) for c in cols], axis=2)

# a Bernoulli response carries no information about a dispersion parameter, so phi = 1
fit = gee(y, X, "independence", scale_fixed=True)
se_naive, se_robust = np.sqrt(np.diag(fit["naive"])), np.sqrt(np.diag(fit["robust"]))
names = ["intercept", "train", "bus", "car", "ttme/10", "gc/10"]
for k, name in enumerate(names):
    print(f"{name:10s} {fit['beta'][k]:8.4f}  naive {se_naive[k]:.4f}  robust {se_robust[k]:.4f}")
print("exchangeable alpha:", round(gee(y, X, 'exchangeable', scale_fixed=True)['alpha'], 4))
# <</data>>

logistic = sm.GLM(y.ravel(), X.reshape(m * n, -1), family=sm.families.Binomial()).fit()
assert np.allclose(fit["beta"], logistic.params, atol=1e-7)     # independence GEE = ordinary GLM
assert np.allclose(se_naive, logistic.bse, atol=1e-6)
phi_pearson = logistic.pearson_chi2 / logistic.df_resid         # not usable: Bernoulli data
ratio = se_robust / se_naive
assert ratio.max() > 1.5 and ratio.min() < 1.0     # the correction goes both ways
assert np.argmax(ratio) == 4 and np.argmin(ratio) == 5

alpha_ex = gee(y, X, "exchangeable", scale_fixed=True)["alpha"]
assert -1 / (n - 1) < alpha_ex < -0.25             # near the structural lower bound -1/(n-1)
bound = np.eye(n) - (1 - np.eye(n)) / (n - 1)
assert abs(np.linalg.eigvalsh(bound).min()) < 1e-12
assert np.linalg.eigvalsh(gee(y, X, "exchangeable", scale_fixed=True)["R"]).min() < 0.1

# cross-check the point estimates and the sandwich against statsmodels
smgee = sm.GEE(y.ravel(), X.reshape(m * n, -1), groups=d["individual"],
               family=sm.families.Binomial()).fit()
assert np.allclose(smgee.params, fit["beta"], atol=1e-6)
assert np.allclose(smgee.bse, se_robust, rtol=0.02)

# ---- coverage when the working correlation is wrong --------------------------
# One covariate varies within the cluster, one is fixed for the whole cluster (a treatment
# assigned to clusters). The independence working correlation is wrong for both.
# <<simulation>>
TAU, BETA = 1.5, np.array([-0.5, 1.0, 0.8])
node, w = hermegauss(60)
w = w / np.sqrt(2 * np.pi)


def marginal_target():
    """The population-averaged coefficients: the root of E[x {m(x) - expit(x'b)}] = 0,
    with x1 ~ N(0,1) within clusters and x2 ~ Bernoulli(1/2) between them."""
    x1 = np.repeat(node, 2)
    x2 = np.tile([0.0, 1.0], len(node))
    weight = np.repeat(w, 2) / 2
    eta = BETA[0] + BETA[1] * x1 + BETA[2] * x2
    mm = np.sum(w * expit(eta[:, None] + TAU * node), axis=1)      # the marginal curve
    Z = np.column_stack([np.ones_like(x1), x1, x2])

    def score(par):
        return Z.T @ (weight * (mm - expit(Z @ par)))

    sol = root(score, BETA, tol=1e-13)
    assert sol.success
    return sol.x


target = marginal_target()
assert np.all(np.abs(target[1:]) < np.abs(BETA[1:]))
# <</simulation>>

rng = np.random.default_rng(40_005)
ms, reps, nclus = [15, 30, 60, 120, 240], 600, 8
cover = {k: [] for k in ["naive", "robust", "corrected"]}
slopes = None
z975 = norm.ppf(0.975)
for m_sim in ms:
    hit = {k: 0 for k in cover}
    draws = []
    tq = t.ppf(0.975, m_sim - 3)
    for _ in range(reps):
        x1 = rng.normal(size=(m_sim, nclus))
        x2 = np.repeat(rng.integers(0, 2, m_sim)[:, None], nclus, axis=1).astype(float)
        Xs = np.stack([np.ones((m_sim, nclus)), x1, x2], axis=2)
        u = rng.normal(0, TAU, m_sim)
        ys = (rng.uniform(size=(m_sim, nclus)) < expit(Xs @ BETA + u[:, None])).astype(float)
        f = gee(ys, Xs, "independence", scale_fixed=True)
        draws.append(f["beta"][2])
        err = abs(f["beta"][2] - target[2])
        hit["naive"] += err < z975 * np.sqrt(f["naive"][2, 2])
        hit["robust"] += err < z975 * np.sqrt(f["robust"][2, 2])
        hit["corrected"] += err < tq * np.sqrt(mancl_derouen(f)[2, 2])
    for k in cover:
        cover[k].append(hit[k] / reps)
    if m_sim == ms[-1]:
        slopes = np.array(draws)
    print(f"m = {m_sim:4d}: " + "  ".join(f"{k} {cover[k][-1]:.3f}" for k in cover))

assert cover["naive"][-1] < 0.90                    # the model-based interval never recovers
assert cover["robust"][-1] > 0.93
assert cover["corrected"][0] > cover["robust"][0]   # the correction helps most when m is small
assert abs(slopes.mean() - target[2]) < 0.05        # consistent for the MARGINAL coefficient
assert abs(slopes.mean() - BETA[2]) > 0.05          # and not for the conditional one

gen = Generated("ch40", "gee_modechoice", prefix="gee")
gen.int("m", m)
gen.int("n", n)
gen.num("beta_ttme", fit["beta"][4], 4)
gen.num("beta_gc", fit["beta"][5], 4)
gen.num("se_naive_ttme", se_naive[4], 4)
gen.num("se_robust_ttme", se_robust[4], 4)
gen.num("se_naive_gc", se_naive[5], 4)
gen.num("se_robust_gc", se_robust[5], 4)
gen.num("ratio_min", ratio.min(), 3)
gen.num("ratio_max", ratio.max(), 3)
gen.num("alpha_ex", alpha_ex, 4)
gen.num("ratio_ttme", ratio[4], 3)
gen.num("ratio_gc", ratio[5], 3)
gen.num("eig_min", float(np.linalg.eigvalsh(gee(y, X, "exchangeable", scale_fixed=True)["R"]).min()), 4)
gen.num("phi_pearson", phi_pearson, 3)
gen.num("target_slope", target[2], 4)
gen.num("target_within", target[1], 4)
gen.num("mean_slope", float(slopes.mean()), 4)
gen.num("cover_naive", cover["naive"][-1], 3)
gen.num("cover_robust", cover["robust"][-1], 3)
gen.num("cover_robust_small", cover["robust"][0], 3)
gen.num("cover_corrected_small", cover["corrected"][0], 3)
gen.int("reps", reps)
gen.write()

# ---- figure ------------------------------------------------------------------
use_book_style()
fig, axes = plt.subplots(1, 2, figsize=(5.6, 2.5))

ax = axes[0]
for k, colour, marker in [("naive", COLORS["second"], "s"), ("robust", COLORS["accent"], "o"),
                          ("corrected", COLORS["third"], "^")]:
    ax.plot(ms, cover[k], marker + "-", color=colour, markersize=3.5, label=k)
ax.axhline(0.95, color=COLORS["grid"], linewidth=0.8, zorder=0)
ax.set_xscale("log")
ax.minorticks_off()
ax.set_xticks(ms)
ax.set_xticklabels([str(k) for k in ms])
ax.set_xlabel("clusters $m$")
ax.set_ylabel("coverage")
ax.set_ylim(0.6, 1.0)
ax.set_title("(a) 95% intervals, cluster-level effect")
ax.legend(frameon=False, fontsize=7, loc="lower right")

ax = axes[1]
ax.hist(slopes, bins=30, color=COLORS["grid"], edgecolor="none")
ax.axvline(target[2], color=COLORS["accent"], label=rf"marginal $\beta^M={target[2]:.2f}$")
ax.axvline(BETA[2], color=COLORS["second"], linestyle="--",
           label=rf"conditional $\beta^C={BETA[2]:.2f}$")
ax.set_xlabel(r"$\hat\beta_2$ from the independence GEE")
ax.set_ylabel("count")
ax.set_title(f"(b) {reps} samples, $m={ms[-1]}$")
ax.legend(frameon=False, fontsize=7, loc="upper left")

fig.tight_layout()
fig.savefig(figure_path("ch40", "gee_coverage"))
