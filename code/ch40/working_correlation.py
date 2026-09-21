"""Chapter 40, Section 5: choosing the working correlation.

(i) A simulation of the efficiency of four working correlations under two truths, separately
for a covariate that varies within the cluster and one that is fixed for the whole cluster.
(ii) QIC on Grunfeld's investment panel (public domain, statsmodels.datasets.grunfeld), for the
mean model and for the working correlation.
"""
import matplotlib.pyplot as plt
import numpy as np
import statsmodels.api as sm

from regbook import COLORS, Generated, figure_path, use_book_style

KINDS = ["independence", "exchangeable", "ar1", "unstructured"]


# <<working>>
def correlation(kind, r, p, phi):
    """Moment estimates of the working correlation from Pearson residuals r (m by n)."""
    m, n = r.shape
    if kind == "independence":
        return np.eye(n)
    if kind == "exchangeable":
        off = (np.sum(r.sum(1) ** 2) - np.sum(r**2)) / 2
        a = off / (phi * (m * n * (n - 1) / 2 - p))
        return np.eye(n) + a * (1 - np.eye(n))
    if kind == "ar1":
        a = np.sum(r[:, :-1] * r[:, 1:]) / (phi * (m * (n - 1) - p))
        return a ** np.abs(np.arange(n)[:, None] - np.arange(n))
    if kind == "unstructured":
        R = (r.T @ r) / (phi * m)
        return R / np.sqrt(np.outer(np.diag(R), np.diag(R)))
    raise ValueError(kind)
# <</working>>


def gee(y, X, kind, family="poisson", maxit=100):
    m, n, p = X.shape
    beta = np.zeros(p)
    beta[0] = np.log(max(y.mean(), 1e-3)) if family != "binomial" else 0.0
    for _ in range(maxit):
        eta = X @ beta
        mu = np.exp(eta)
        dmu, V = (mu, mu) if family == "poisson" else (mu, mu**2)
        r = (y - mu) / np.sqrt(V)
        phi = np.sum(r**2) / (m * n - p)
        R = correlation(kind, r, p, phi)
        Rinv = np.linalg.pinv(R)
        SX = (dmu / np.sqrt(V))[..., None] * X
        A = np.einsum("mjp,jk,mkq->pq", SX, Rinv, SX)
        step = np.linalg.solve(A, np.einsum("mjp,jk,mk->p", SX, Rinv, r))
        beta = beta + step
        if np.max(np.abs(step)) < 1e-10:
            break
    u = np.einsum("mjp,jk,mk->mp", SX, Rinv, r)
    Ainv = np.linalg.inv(A)
    return beta, Ainv @ (u.T @ u) @ Ainv, phi, R


# ---- efficiency simulation ---------------------------------------------------
M_SIM, N_SIM, REPS = 60, 6, 400
BETA = np.array([0.4, 0.5, 0.3])       # intercept, within-cluster slope, cluster-level slope
TAU, RHO = 0.6, 0.8


def draw(rng, truth):
    x = rng.normal(size=(M_SIM, N_SIM))                      # varies within the cluster
    z = np.repeat(rng.normal(size=(M_SIM, 1)), N_SIM, axis=1)  # fixed within the cluster
    X = np.stack([np.ones((M_SIM, N_SIM)), x, z], axis=2)
    if truth == "exchangeable":
        u = np.repeat(rng.normal(0, TAU, (M_SIM, 1)), N_SIM, axis=1)
    else:                                                    # a Gaussian AR(1) latent process
        e = rng.normal(0, TAU, (M_SIM, N_SIM))
        u = np.empty_like(e)
        u[:, 0] = e[:, 0]
        for j in range(1, N_SIM):
            u[:, j] = RHO * u[:, j - 1] + np.sqrt(1 - RHO**2) * e[:, j]
    y = rng.poisson(np.exp(X @ BETA + u)).astype(float)
    return y, X


rng = np.random.default_rng(40_006)
sd = {}
for truth in ["exchangeable", "ar1"]:
    draws = {k: [] for k in KINDS}
    for _ in range(REPS):
        y, X = draw(rng, truth)
        for k in KINDS:
            draws[k].append(gee(y, X, k)[0][1:])
    sd[truth] = {k: np.std(np.array(v), axis=0, ddof=1) for k, v in draws.items()}
    for k in KINDS:
        print(f"{truth:13s} {k:14s} sd(within) {sd[truth][k][0]:.4f}  "
              f"sd(cluster) {sd[truth][k][1]:.4f}")

# the working correlation that matches the truth is the most efficient for the within-cluster
# covariate; for the cluster-level covariate all four are within a few per cent of each other
for truth, best in [("exchangeable", "exchangeable"), ("ar1", "ar1")]:
    within = {k: sd[truth][k][0] for k in KINDS}
    assert min(within, key=within.get) in (best, "unstructured")
    assert within["independence"] > 1.10 * within[best]
    cluster = np.array([sd[truth][k][1] for k in KINDS])
    assert cluster.max() / cluster.min() < 1.10
    assert cluster.argmin() == 0            # independence is not beaten at the cluster level

eff = {t: {k: (sd[t]["independence"][0] / sd[t][k][0]) ** 2 for k in KINDS} for t in sd}
eff_c = {t: {k: (sd[t]["independence"][1] / sd[t][k][1]) ** 2 for k in KINDS} for t in sd}

# ---- QIC on the investment panel ---------------------------------------------
# <<qic>>
d = sm.datasets.grunfeld.load_pandas().data.sort_values(["firm", "year"])
m, n = d["firm"].nunique(), 20
y = d["invest"].to_numpy().reshape(m, n)
lv = np.log(d["value"].to_numpy()).reshape(m, n)
lc = np.log(d["capital"].to_numpy()).reshape(m, n)


def qic(y, X, kind, phi):
    """Pan's QIC: the quasi-likelihood of the independence model at the fitted mean, plus a
    penalty built from the sandwich. For V(mu) = mu^2 the quasi-likelihood is -y/mu - log mu.
    The dispersion phi is held fixed across the models compared, as in Mallows's Cp."""
    beta, robust, _, _ = gee(y, X, kind, family="gamma")
    mu = np.exp(X @ beta)
    Q = np.sum(-y / mu - np.log(mu)) / phi
    omega = np.einsum("mjp,mjq->pq", X, X) / phi     # independence information, log link
    return -2 * Q + 2 * np.trace(omega @ robust)


X_full = np.stack([np.ones((m, n)), lv, lc], axis=2)
X_small = np.stack([np.ones((m, n)), lv], axis=2)
phi_ref = gee(y, X_full, "independence", family="gamma")[2]     # from the largest model
for name, XX in [("value + capital", X_full), ("value only", X_small)]:
    for kind in KINDS:
        print(f"{name:16s} {kind:14s} QIC {qic(y, XX, kind, phi_ref):10.2f}")
# <</qic>>

qic_full = {k: qic(y, X_full, k, phi_ref) for k in KINDS}
qic_small = {k: qic(y, X_small, k, phi_ref) for k in KINDS}
assert min(qic_full.values()) < min(qic_small.values())         # capital earns its place
best_kind = min(qic_full, key=qic_full.get)
# QIC ranks the mean models sensibly; as a chooser of the working correlation it merely
# rewards the structure whose sandwich is smallest, and here that is independence
assert best_kind == "independence"
assert qic_full["unstructured"] == max(qic_full.values())
_, _, _, R_ar = gee(y, X_full, "ar1", family="gamma")
assert R_ar[0, 1] > 0.7                                          # strong serial correlation

# the same beta under two working correlations: by thm-gmm-gee(a) they estimate the same
# thing when the mean model is right, so a large gap accuses the mean model
b_ind = gee(y, X_full, "independence", family="gamma")[0]
b_exc, rob_exc, _, _ = gee(y, X_full, "exchangeable", family="gamma")
se_exc = np.sqrt(np.diag(rob_exc))
print("independence beta", b_ind.round(4), " exchangeable beta", b_exc.round(4))
assert abs(b_ind[1] - b_exc[1]) > 1.0 * se_exc[1]

# the Mundlak remedy of exr-gmm-within-between: add the firm means, so that lv and lc carry
# only within-firm variation. The gap closes exactly -- both working correlations then give
# the same estimate to machine precision, because the two sources of variation are separated.
mv = np.repeat(lv.mean(1)[:, None], n, axis=1)
mc = np.repeat(lc.mean(1)[:, None], n, axis=1)
X_mund = np.stack([np.ones((m, n)), lv, lc, mv, mc], axis=2)
b_mund_ind = gee(y, X_mund, "independence", family="gamma")[0]
b_mund_exc = gee(y, X_mund, "exchangeable", family="gamma")[0]
assert np.max(np.abs(b_mund_ind - b_mund_exc)) < 1e-8
print("Mundlak beta (independence)", b_mund_ind.round(4))

gen = Generated("ch40", "working_correlation", prefix="wc")
gen.int("m_sim", M_SIM)
gen.int("n_sim", N_SIM)
gen.int("reps", REPS)
for t in sd:
    for k in KINDS:
        gen.num(f"eff_{t}_{k}", eff[t][k], 3)
        gen.num(f"effc_{t}_{k}", eff_c[t][k], 3)
gen.num("rho_ar", R_ar[0, 1], 3)
gen.num("b_ind_value", b_ind[1], 4)
gen.num("b_ind_capital", b_ind[2], 4)
gen.num("b_exc_value", b_exc[1], 4)
gen.num("b_exc_capital", b_exc[2], 4)
gen.num("se_exc_value", se_exc[1], 4)
gen.num("b_mund_value", b_mund_ind[1], 4)
gen.num("b_mund_capital", b_mund_ind[2], 4)
for k in KINDS:
    gen.num(f"qic_full_{k}", qic_full[k], 1)
    gen.num(f"qic_small_{k}", qic_small[k], 1)
gen.text("best_kind", best_kind)
gen.write()

# ---- figure ------------------------------------------------------------------
use_book_style()
fig, axes = plt.subplots(1, 2, figsize=(5.6, 2.5), sharey=True)
width, pos = 0.38, np.arange(len(KINDS))
short = ["indep.", "exch.", "AR(1)", "unstr."]
for ax, truth, title in [(axes[0], "exchangeable", "(a) exchangeable truth"),
                         (axes[1], "ar1", r"(b) AR(1) truth, $\rho=0.8$")]:
    ax.bar(pos - width / 2, [eff[truth][k] for k in KINDS], width,
           color=COLORS["accent"], label="within-cluster covariate")
    ax.bar(pos + width / 2, [eff_c[truth][k] for k in KINDS], width,
           color=COLORS["second"], label="cluster-level covariate")
    ax.axhline(1.0, color=COLORS["grid"], linewidth=0.8, zorder=0)
    ax.set_xticks(pos)
    ax.set_xticklabels(short, fontsize=7.5)
    ax.set_title(title)
axes[0].set_ylim(0, 2.35)
axes[0].set_ylabel("efficiency relative to independence")
axes[0].legend(frameon=False, fontsize=7, loc="upper left")
fig.tight_layout()
fig.savefig(figure_path("ch40", "efficiency"))
