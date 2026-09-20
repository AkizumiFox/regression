"""Chapter 31, Section 4: feasible generalized least squares with equicorrelated clusters.

Clusters of unequal size share a random shock. The covariance is known up to the
intraclass correlation rho, which is estimated from the ordinary least squares residuals
and then treated as if it were known. The simulation measures what that costs.
"""
import matplotlib.pyplot as plt
import numpy as np
import statsmodels.api as sm
from scipy import stats

from regbook import COLORS, Generated, figure_path, use_book_style

PATTERN = (2, 4, 8)           # the cluster sizes repeat in this order


def design(blocks):
    """Cluster sizes, indicator matrix and model matrix [1, w, z] for `blocks` repeats."""
    sizes = np.array(PATTERN * blocks)
    n, G = sizes.sum(), len(sizes)
    Z = np.zeros((n, G))
    Z[np.arange(n), np.repeat(np.arange(G), sizes)] = 1.0
    z = Z @ np.tile([-1.0, 0.0, 1.0], blocks)                   # a cluster-level regressor
    w = np.concatenate([np.arange(s) - (s - 1) / 2 for s in sizes])   # centred in the cluster
    return sizes, Z, np.column_stack([np.ones(n), w, z])


# <<helpers>>
def quasi_demean(A, Z, sizes, rho):
    """Whitening for equicorrelation rho: rows of A are vectors in R^n."""
    c = 1.0 - np.sqrt((1 - rho) / (1 - rho + sizes * rho))
    return (A - (A @ Z / sizes * c) @ Z.T) / np.sqrt(1 - rho)


def fit(X, Y, Z, sizes, rho):
    """GLS coefficients and the model-based covariance factor, for many data sets at once."""
    Xw = quasi_demean(X.T, Z, sizes, rho).T
    Yw = quasi_demean(Y, Z, sizes, rho)
    A = np.linalg.inv(Xw.T @ Xw)
    B = Yw @ Xw @ A
    R = Yw - B @ Xw.T
    s2 = np.sum(R * R, axis=1) / (X.shape[0] - X.shape[1])
    return B, A, s2


def moment_rho(R, Z, sizes):
    """ANOVA estimate of the intraclass correlation from residual vectors R."""
    G, n = len(sizes), R.shape[1]
    means = R @ Z / sizes
    between = np.sum(means**2 * sizes, axis=1) / G
    within = (np.sum(R * R, axis=1) - np.sum(means**2 * sizes, axis=1)) / (n - G)
    mbar = (n - np.sum(sizes**2) / n) / (G - 1)
    est = (between - within) / (between + (mbar - 1) * within)
    return np.clip(est, 0.0, 0.95)
# <</helpers>>


def experiment(blocks, rho, reps, seed):
    """One row of the table: variances and coverages for OLS, GLS and feasible GLS."""
    sizes, Z, X = design(blocks)
    n, p, G = X.shape[0], X.shape[1], len(sizes)
    beta = np.array([1.0, 0.5, 0.25])
    rng = np.random.default_rng(seed + 1)
    out = {k: [] for k in ("ols", "gls", "fgls", "cov_gls", "cov_naive", "cov_robust")}
    tq = stats.t.ppf(0.975, n - p)
    for start in range(0, reps, 2000):
        r = min(2000, reps - start)
        E = (np.sqrt(rho) * rng.normal(size=(r, G)) @ Z.T
             + np.sqrt(1 - rho) * rng.normal(size=(r, n)))
        Y = X @ beta + E
        B_ols, _, _ = fit(X, Y, Z, sizes, 0.0)
        B_gls, A_gls, s2_gls = fit(X, Y, Z, sizes, rho)
        rho_hat = moment_rho(Y - B_ols @ X.T, Z, sizes)
        B_f = np.empty_like(B_ols)
        se_f = np.empty(r)
        se_rob = np.empty(r)
        for i in range(r):
            b, A, s2 = fit(X, Y[i:i + 1], Z, sizes, rho_hat[i])
            B_f[i] = b[0]
            se_f[i] = np.sqrt(s2[0] * A[2, 2])
            # the cluster-robust covariance OF THE FEASIBLE GLS ESTIMATE: the bread is
            # A = (X' V(rho-hat)^{-1} X)^{-1} and the meat is built from whitened quantities
            Xw = quasi_demean(X.T, Z, sizes, rho_hat[i]).T
            ew = quasi_demean(Y[i:i + 1], Z, sizes, rho_hat[i])[0] - Xw @ b[0]
            tot = Z.T @ (ew[:, None] * Xw)                       # one row per cluster
            meat = tot.T @ tot * G / (G - 1)
            Vr = A @ meat @ A
            se_rob[i] = np.sqrt(Vr[2, 2])
        out["ols"].append(B_ols[:, 2])
        out["gls"].append(B_gls[:, 2])
        out["fgls"].append(B_f[:, 2])
        out["cov_gls"].append(np.abs(B_gls[:, 2] - beta[2])
                              <= tq * np.sqrt(s2_gls * A_gls[2, 2]))
        out["cov_naive"].append(np.abs(B_f[:, 2] - beta[2]) <= tq * se_f)
        out["cov_robust"].append(np.abs(B_f[:, 2] - beta[2])
                                 <= stats.t.ppf(0.975, G - 1) * se_rob)
    res = {k: np.concatenate(v) for k, v in out.items()}
    var_gls = res["gls"].var()
    return {
        "G": G, "n": n,
        "eff_ols": var_gls / res["ols"].var(),
        "eff_fgls": var_gls / res["fgls"].var(),
        "cov_gls": res["cov_gls"].mean(),
        "cov_naive": res["cov_naive"].mean(),
        "cov_robust": res["cov_robust"].mean(),
    }


rho_true = 0.6
rows = [experiment(b, rho_true, 20_000, seed=100 * b) for b in (1, 2, 4, 8, 16, 32)]
for r in rows:
    print("G = %3d  n = %3d   eff(OLS) %.3f  eff(FGLS) %.3f   coverage: GLS %.3f  "
          "naive %.3f  robust %.3f"
          % (r["G"], r["n"], r["eff_ols"], r["eff_fgls"], r["cov_gls"], r["cov_naive"],
             r["cov_robust"]))

assert all(0.93 < r["cov_gls"] < 0.97 for r in rows)          # known rho: exact
assert rows[0]["cov_naive"] < 0.87                            # three clusters: far too liberal
assert rows[-1]["cov_naive"] > 0.93                           # ninety-six clusters: about right
assert rows[0]["eff_fgls"] < rows[-1]["eff_fgls"]
assert rows[-1]["eff_fgls"] > 0.97
# the robust interval is never worse than the model-based one, and neither reaches 0.95
assert all(r["cov_robust"] >= r["cov_naive"] - 0.001 for r in rows)
assert all(r["cov_robust"] < 0.95 for r in rows)
assert rows[-1]["cov_robust"] - rows[-1]["cov_naive"] < 0.005   # indistinguishable at G = 96

gen = Generated("ch31", "feasible_gls", prefix="fgls")
gen.num("rho", rho_true, 1)
gen.text("pattern", ", ".join(str(s) for s in PATTERN))
for tag, r in (("small", rows[0]), ("mid", rows[2]), ("large", rows[-1])):
    gen.int(f"G_{tag}", r["G"])
    gen.int(f"n_{tag}", r["n"])
    gen.num(f"eff_ols_{tag}", r["eff_ols"], 3)
    gen.num(f"var_ratio_ols_{tag}", 1.0 / r["eff_ols"], 2)
    gen.num(f"eff_fgls_{tag}", r["eff_fgls"], 3)
    gen.num(f"cov_naive_{tag}", r["cov_naive"], 3)
    gen.num(f"cov_robust_{tag}", r["cov_robust"], 3)
    gen.num(f"cov_gls_{tag}", r["cov_gls"], 3)

# ---- the El Nino series with the intraclass correlation estimated ------------
# <<elnino>>
import numpy as np
import statsmodels.api as sm

frame = sm.datasets.elnino.load_pandas().data
temps = frame.iloc[:, 1:].to_numpy()
G, m = temps.shape
n = G * m
y = temps.ravel()
month = np.tile(np.arange(1, m + 1), G)
time = (np.repeat(frame["YEAR"].to_numpy(), m) - 1980.0) + (month - 6.5) / 12.0
cols = [np.ones(n), time / 10.0]
for k in (1, 2, 3):
    cols += [np.cos(2 * np.pi * k * month / m), np.sin(2 * np.pi * k * month / m)]
X = np.column_stack(cols)
Zn = np.zeros((n, G))
Zn[np.arange(n), np.repeat(np.arange(G), m)] = 1.0
sizes_n = np.full(G, m)

b_ols = np.linalg.solve(X.T @ X, X.T @ y)
rho_hat = float(moment_rho((y - X @ b_ols)[None, :], Zn, sizes_n)[0])
b_f, A_f, s2_f = fit(X, y[None, :], Zn, sizes_n, rho_hat)
se_f = np.sqrt(s2_f[0] * A_f[1, 1])
print(f"El Nino: rho_hat = {rho_hat:.3f}, trend {b_f[0, 1]:.4f} (s.e. {se_f:.4f})")
# <</elnino>>

res = (y - X @ b_ols)[:, None] * X
tot = Zn.T @ res
A_ols = np.linalg.inv(X.T @ X)
V_rob = A_ols @ (tot.T @ tot * G / (G - 1)) @ A_ols
se_rob = np.sqrt(V_rob[1, 1])
se_naive = np.sqrt(np.sum((y - X @ b_ols) ** 2) / (n - X.shape[1]) * A_ols[1, 1])
print(f"trend s.e.: usual {se_naive:.4f}, cluster-robust {se_rob:.4f}, feasible GLS {se_f:.4f}")
assert se_rob > 2 * se_naive and abs(se_rob / se_f - 1) < 0.2

gen.num("nino_rho", rho_hat, 3)
gen.num("nino_trend", b_f[0, 1], 4)
gen.num("nino_se", se_f, 4)
gen.num("nino_se_naive", se_naive, 4)
gen.num("nino_se_robust", se_rob, 4)
gen.write()

# ---- figure ------------------------------------------------------------------
use_book_style()
Gs = np.array([r["G"] for r in rows], dtype=float)
fig, axes = plt.subplots(1, 2, figsize=(5.8, 2.4))
ax = axes[0]
ax.plot(Gs, [r["eff_fgls"] for r in rows], "o-", color=COLORS["accent"], ms=3,
        label="feasible GLS")
ax.plot(Gs, [r["eff_ols"] for r in rows], "s-", color=COLORS["second"], ms=3, label="OLS")
ax.axhline(1, color=COLORS["grid"], linewidth=0.6, zorder=0)
ax.set_xscale("log")
ax.set_xticks(Gs)
ax.set_xticklabels([f"{int(g)}" for g in Gs])
ax.set_ylim(0.85, 1.02)
ax.set_xlabel("number of clusters")
ax.set_ylabel("efficiency")
ax.set_title("(a) efficiency relative to GLS")
ax.legend(loc="lower right", frameon=False)

ax = axes[1]
ax.plot(Gs, [r["cov_gls"] for r in rows], "^-", color=COLORS["third"], ms=3,
        label=r"GLS, $\rho$ known")
ax.plot(Gs, [r["cov_naive"] for r in rows], "o-", color=COLORS["accent"], ms=3,
        label="feasible GLS, usual s.e.")
ax.plot(Gs, [r["cov_robust"] for r in rows], "s-", color=COLORS["second"], ms=3,
        label="feasible GLS, robust s.e.")
ax.axhline(0.95, color=COLORS["grid"], linewidth=0.6, zorder=0)
ax.set_xscale("log")
ax.set_xticks(Gs)
ax.set_xticklabels([f"{int(g)}" for g in Gs])
ax.set_xlabel("number of clusters")
ax.set_ylabel("coverage")
ax.set_title("(b) nominal 95% intervals")
ax.legend(loc="lower right", frameon=False)
fig.tight_layout()
fig.savefig(figure_path("ch31", "feasible_gls"))
