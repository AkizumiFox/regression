"""Chapter 33, Section 2: equicorrelation within clusters, the design effect, and
what a cluster-robust standard error repairs.

Simulated clusters of equal size with one cluster-level regressor and one regressor
that varies only within clusters. The first has its variance inflated by the factor
1 + (m-1)rho, the second deflated by 1 - rho.
"""
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats

from regbook import COLORS, Generated, figure_path, use_book_style

# <<structure>>
m, rho, sigma2 = 5, 0.4, 1.0                         # cluster size, intraclass correlation
S = sigma2 * ((1 - rho) * np.eye(m) + rho * np.ones((m, m)))
w = np.linalg.eigvalsh(S)                            # one large, m-1 equal small eigenvalues
print("eigenvalues of the cluster covariance:", w.round(4))
print("large:", sigma2 * (1 + (m - 1) * rho), "  small:", sigma2 * (1 - rho))
# <</structure>>
assert np.allclose(sorted(w), sorted([sigma2 * (1 + (m - 1) * rho)]
                                     + [sigma2 * (1 - rho)] * (m - 1)))
Sinv = np.linalg.inv(S)
J = np.ones((m, m)) / m
assert np.allclose(Sinv, (np.eye(m) - (m * rho / (1 + (m - 1) * rho)) * J)
                   / (sigma2 * (1 - rho)))


def design(G, m, seed):
    """Cluster indicator, a cluster-level regressor and a within-cluster regressor."""
    rng = np.random.default_rng(seed)
    g = np.repeat(np.arange(G), m)
    xb = np.repeat(rng.normal(0, 1, G), m)           # constant inside each cluster
    xw = rng.normal(0, 1, G * m)
    xw = xw - np.repeat(xw.reshape(G, m).mean(axis=1), m)   # sums to zero inside each cluster
    X = np.column_stack([np.ones(G * m), xb, xw])
    return g, X


# <<deff>>
G = 40
g, X = design(G, m, 2024)
V = np.kron(np.eye(G), S)                            # block diagonal, one block per cluster
XtXinv = np.linalg.inv(X.T @ X)
cov_true = XtXinv @ X.T @ V @ X @ XtXinv             # the real covariance of the OLS estimate
cov_iid = sigma2 * XtXinv                            # what the usual formula reports
deff = np.diag(cov_true) / np.diag(cov_iid)
print("design effects (intercept, between, within):", deff.round(4))
print("predicted:", [1 + (m - 1) * rho, 1 + (m - 1) * rho, 1 - rho])
# <</deff>>
assert np.allclose(deff, [1 + (m - 1) * rho, 1 + (m - 1) * rho, 1 - rho])
assert np.allclose(cov_true[1, 2], 0, atol=1e-12)

# the two regressors are eigenvectors of V ------------------------------------
assert np.allclose(V @ X[:, 1], sigma2 * (1 + (m - 1) * rho) * X[:, 1])
assert np.allclose(V @ X[:, 2], sigma2 * (1 - rho) * X[:, 2])
# so ordinary least squares is already the generalized least squares fit
Vinv = np.linalg.inv(V)
beta_gls_map = np.linalg.inv(X.T @ Vinv @ X) @ X.T @ Vinv
assert np.allclose(beta_gls_map, XtXinv @ X.T)

# how badly the usual standard error is wrong ----------------------------------
n, p = X.shape
Mres = np.eye(n) - X @ XtXinv @ X.T
Es2 = np.trace(Mres @ V) / (n - p)                   # expectation of the usual s^2
ratio_between = Es2 * XtXinv[1, 1] / cov_true[1, 1]
ratio_within = Es2 * XtXinv[2, 2] / cov_true[2, 2]
print(f"E(v)/Var: between {ratio_between:.4f}, within {ratio_within:.4f}")
# the closed form: p_b = 2 columns in the whole-cluster stratum, p_w = 1 in the other
lam_w, lam_s = sigma2 * (1 + (m - 1) * rho), sigma2 * (1 - rho)
p_b, p_w = 2, 1
Es2_formula = (lam_w * (G - p_b) + lam_s * (n - G - p_w)) / (n - p)
assert np.isclose(Es2, Es2_formula)
assert np.isclose(ratio_between, Es2_formula / lam_w)
assert np.isclose(ratio_within, Es2_formula / lam_s)


# <<robust>>
def cluster_sandwich(X, resid, g):
    """The cluster-robust covariance estimate: one outer product per cluster."""
    XtXinv = np.linalg.inv(X.T @ X)
    meat = np.zeros((X.shape[1], X.shape[1]))
    for k in np.unique(g):
        s = X[g == k].T @ resid[g == k]               # the cluster's score contribution
        meat += np.outer(s, s)
    return XtXinv @ meat @ XtXinv


rng_one = np.random.default_rng(55)
e_one = (np.linalg.cholesky(S) @ rng_one.standard_normal((m, G))).T.reshape(-1)
y_one = X[:, 1] + e_one                              # the between coefficient is 1
b_one = XtXinv @ X.T @ y_one
res_one = y_one - X @ b_one
s2_one = res_one @ res_one / (n - p)
se_naive_one = np.sqrt(s2_one * XtXinv[1, 1])
se_cr_one = np.sqrt(cluster_sandwich(X, res_one, g)[1, 1] * G / (G - 1))
se_true_one = np.sqrt(cov_true[1, 1])
print(f"between coefficient {b_one[1]:.4f}: naive se {se_naive_one:.4f}, "
      f"cluster-robust {se_cr_one:.4f}, true sd {se_true_one:.4f}")
# <</robust>>
assert se_cr_one > 1.4 * se_naive_one

# ---- simulated coverage ------------------------------------------------------
def coverage(G, B=3000, seed=7, m=m, rho=rho):
    S = (1 - rho) * np.eye(m) + rho * np.ones((m, m))
    L = np.linalg.cholesky(S)
    g, X = design(G, m, 1000 + G)
    n, p = X.shape
    XtXinv = np.linalg.inv(X.T @ X)
    Vinv = np.kron(np.eye(G), np.linalg.inv(S))
    A = np.linalg.inv(X.T @ Vinv @ X)
    gls_map = A @ X.T @ Vinv
    corr = G / (G - 1) * (n - 1) / (n - p)           # the usual small-sample correction
    idx = [np.flatnonzero(g == k) for k in range(G)]
    rng = np.random.default_rng(seed)
    hit = np.zeros(3)
    for _ in range(B):
        e = (L @ rng.standard_normal((m, G))).T.reshape(-1)
        y = X[:, 1] + e                              # true coefficient of the between regressor is 1
        bhat = XtXinv @ X.T @ y
        res = y - X @ bhat
        s2 = res @ res / (n - p)
        se_naive = np.sqrt(s2 * XtXinv[1, 1])
        meat = sum(np.outer(X[i].T @ res[i], X[i].T @ res[i]) for i in idx)
        se_cr = np.sqrt(corr * (XtXinv @ meat @ XtXinv)[1, 1])
        bg = gls_map @ y
        rg = y - X @ bg
        s2g = rg @ Vinv @ rg / (n - p)
        se_gls = np.sqrt(s2g * A[1, 1])
        for i, (b, se, dfree) in enumerate([(bhat[1], se_naive, n - p),
                                            (bhat[1], se_cr, G - 1),
                                            (bg[1], se_gls, n - p)]):
            hit[i] += abs(b - 1) <= stats.t.ppf(0.975, dfree) * se
    return hit / B


Gs = [8, 12, 20, 40, 80, 160]
cov_rows = np.array([coverage(Gv) for Gv in Gs])
for Gv, row in zip(Gs, cov_rows):
    print(f"G = {Gv:3d}   naive {row[0]:.3f}   cluster-robust {row[1]:.3f}   GLS {row[2]:.3f}")
assert cov_rows[:, 0].max() < 0.90                   # the naive interval never comes close
assert cov_rows[0, 1] < cov_rows[-1, 1]              # cluster-robust improves with more clusters
assert cov_rows[-1, 1] > 0.93
assert cov_rows[:, 2].min() > 0.93                   # GLS is right throughout

gen = Generated("ch33", "design_effect")
gen.int("m", m)
gen.num("rho", rho, 1)
gen.num("deff_between", deff[1], 2)
gen.num("deff_within", deff[2], 2)
gen.num("Es2", Es2, 4)
gen.num("ratio_between", ratio_between, 3)
gen.num("ratio_within", ratio_within, 3)
gen.num("se_factor_between", np.sqrt(1 / ratio_between), 2)
gen.int("G", G)
gen.num("b_one", b_one[1], 4)
gen.num("se_naive_one", se_naive_one, 4)
gen.num("se_cr_one", se_cr_one, 4)
gen.num("se_true_one", se_true_one, 4)
for Gv, row in zip(Gs, cov_rows):
    gen.num(f"cov_naive_{Gv}", row[0], 3)
    gen.num(f"cov_cr_{Gv}", row[1], 3)
    gen.num(f"cov_gls_{Gv}", row[2], 3)
gen.write()

# ---- figure ------------------------------------------------------------------
use_book_style()
fig, axes = plt.subplots(1, 2, figsize=(5.6, 2.4))
ax = axes[0]
rr = np.linspace(0, 0.9, 200)
for mm, style in [(3, ":"), (5, "-"), (10, "--")]:
    ax.plot(rr, 1 + (mm - 1) * rr, style, color=COLORS["accent"], label=f"between, m={mm}")
ax.plot(rr, 1 - rr, color=COLORS["second"], label="within, any m")
ax.axhline(1, color=COLORS["grid"], linewidth=0.6, zorder=0)
ax.set_xlabel(r"intraclass correlation $\rho$")
ax.set_ylabel("design effect")
ax.set_title("(a) variance inflation")
ax.legend(fontsize=6.5, frameon=False)
ax = axes[1]
ax.plot(Gs, cov_rows[:, 0], "o-", color=COLORS["second"], markersize=3, label="naive")
ax.plot(Gs, cov_rows[:, 1], "s-", color=COLORS["accent"], markersize=3, label="cluster-robust")
ax.plot(Gs, cov_rows[:, 2], "^-", color=COLORS["third"], markersize=3, label="GLS")
ax.axhline(0.95, color=COLORS["grid"], linewidth=0.6, zorder=0)
ax.set_xscale("log")
ax.set_xticks(Gs)
ax.set_xticklabels([str(v) for v in Gs])
ax.minorticks_off()
ax.set_xlabel("number of clusters")
ax.set_ylabel("coverage")
ax.set_ylim(0.55, 1.0)
ax.set_title("(b) coverage of a 95% interval")
ax.legend(fontsize=6.5, frameon=False, loc="lower right")
fig.tight_layout()
fig.savefig(figure_path("ch33", "design_effect"))
