"""Chapter 23, Section 5: permutation tests of one coefficient in the presence of nuisance regressors.

(a) Brownlee's stack loss data (statsmodels.datasets.stackloss, public domain): test the acid-concentration
    coefficient with air flow and water temperature in the model, by the F test and five permutation schemes.
(b) A fixed-seed simulation of the size of the six tests: n = 12, normal errors, a nuisance regressor
    correlated with the regressor of interest, which has one high-leverage value; the nuisance
    coefficient gamma ranges from 0 to 30.
"""
import matplotlib.pyplot as plt
import numpy as np
import statsmodels.api as sm
from scipy import stats

from regbook import COLORS, Generated, figure_path, use_book_style

data = sm.datasets.stackloss.load_pandas().data

# <<schemes>>
def proj(Z):
    """Orthogonal projection onto C(Z), from a QR factorization."""
    Q, _ = np.linalg.qr(Z)
    return Q @ Q.T


def perm_rows(v, B, rng):
    """B random permutations of the vector v, one per row."""
    return v[np.argsort(rng.random((B, len(v))), axis=1)]


def perm_pvalues(y, x, Z, B, rng):
    """p-values for H0: beta = 0 in E(y) = Z gamma + x beta, by six methods."""
    n, p = len(y), Z.shape[1] + 1
    MZ = proj(Z)
    xt = x - MZ @ x                               # x residualized on Z
    yt = y - MZ @ y                               # reduced-model residuals
    e_full = yt - xt * (xt @ yt) / (xt @ xt)      # full-model residuals (FWL)

    def F(Y):
        """F statistic for beta = 0, for each row of Y (model matrix [Z, x] fixed)."""
        num = (Y @ xt) ** 2 / (xt @ xt)
        rss_Z = np.sum(Y ** 2, axis=-1) - np.sum((Y @ MZ) * Y, axis=-1)
        return num / ((rss_Z - num) / (n - p))

    F_obs = F(y)
    out = {"F test": stats.f.sf(F_obs, 1, n - p)}

    def pval(F_star, F0=F_obs):
        return (1 + np.sum(F_star >= F0)) / (len(F_star) + 1)

    out["Manly"] = pval(F(perm_rows(y, B, rng)))                    # permute y
    Xs = perm_rows(x, B, rng)                                       # permute x (Draper-Stoneman)
    xts = Xs - Xs @ MZ
    num = (xts @ y) ** 2 / np.sum(xts ** 2, axis=1)
    out["Draper-Stoneman"] = pval(num / ((yt @ yt - num) / (n - p)))
    Yk = perm_rows(yt, B, rng)                                      # Kennedy: no Z in the refit
    num = (Yk @ xt) ** 2 / (xt @ xt)
    out["Kennedy"] = pval(num / ((np.sum(Yk ** 2, axis=1) - num) / (n - p)))
    out["Freedman-Lane"] = pval(F(perm_rows(yt, B, rng)))           # permute reduced residuals
    out["ter Braak"] = pval(F(perm_rows(e_full, B, rng)))           # permute full residuals
    return out, F_obs
# <</schemes>>


# <<stackloss>>
rng = np.random.default_rng(2307)
y = data["STACKLOSS"].to_numpy()
x = data["ACIDCONC"].to_numpy()
Z = np.column_stack([np.ones(len(y)), data[["AIRFLOW", "WATERTEMP"]]])
pvals, F_obs = perm_pvalues(y, x, Z, B=19_999, rng=rng)
print(f"F = {F_obs:.3f}")
for m, pv in pvals.items():
    print(f"{m:16s} p = {pv:.4f}")
# <</stackloss>>

# the Freedman-Lane statistic equals the F statistic of the permuted reduced residuals, i.e. refitting
# y* = M_Z y + P e_0 gives the same F as regressing P e_0 (Proposition bs-freedman-lane (a))
MZ = proj(Z)
e0 = y - MZ @ y
perm = rng.permutation(len(y))
Xf = np.column_stack([Z, x])
def F_refit(yy):
    full = sm.OLS(yy, Xf).fit()
    red = sm.OLS(yy, Z).fit()
    return (red.ssr - full.ssr) / (full.ssr / (len(yy) - Xf.shape[1]))
assert np.isclose(F_refit(MZ @ y + e0[perm]), F_refit(e0[perm]))
assert np.isclose(F_refit(y), F_obs)
assert np.isclose(sm.OLS(y, Xf).fit().pvalues[-1], pvals["F test"])
for m in ["Manly", "Draper-Stoneman", "Kennedy", "Freedman-Lane", "ter Braak"]:
    assert abs(pvals[m] - pvals["F test"]) < 0.05
# with the same permutations, Kennedy's statistic never exceeds Freedman-Lane's (Proposition bs-kennedy)
xt = x - MZ @ x
Mfull = proj(Xf)
for _ in range(200):
    u = e0[rng.permutation(len(y))]
    num = (u @ xt) ** 2 / (xt @ xt)
    F_fl = num / ((u @ u - u @ Mfull @ u) / (len(y) - 4))
    F_k = num / ((u @ u - num) / (len(y) - 4))
    assert F_k <= F_fl + 1e-12 and np.isclose(F_fl, F_refit(u))

gen = Generated("ch23", "freedman_lane", prefix="fl")
gen.num("F", F_obs, 3)
for m, pv in pvals.items():
    gen.num(f"stack:{m.replace(' ', '').replace('-', '')}", pv, 3)

# ---- (b) size simulation --------------------------------------------------------------------
# <<size>>
def size_study(gamma, reps, B, seed, n=12, alpha=0.05):
    """Rejection rates at level alpha when beta = 0 and the nuisance coefficient is gamma."""
    rng = np.random.default_rng(seed)
    z = np.linspace(-1, 1, n)
    x = 0.8 * z + 0.6 * np.random.default_rng(1).normal(size=n)   # fixed, correlated with z
    x[0] = 6.0                                    # one high-leverage value of x
    Z = np.column_stack([np.ones(n), z])
    rejections = {}
    for _ in range(reps):
        y = 1.0 + gamma * z + rng.normal(size=n)  # the null hypothesis beta = 0 holds
        for m, pv in perm_pvalues(y, x, Z, B, rng)[0].items():
            rejections[m] = rejections.get(m, 0) + (pv <= alpha)
    return {m: r / reps for m, r in rejections.items()}


for gamma in (0.0, 30.0):
    print(gamma, size_study(gamma, reps=300, B=99, seed=2308))
# <</size>>

reps_s, B_s = 5000, 199
gammas = [0.0, 1.0, 3.0, 10.0, 30.0]
sizes = {g: size_study(g, reps_s, B_s, seed=2309) for g in gammas}
for g, v in sizes.items():
    print(f"gamma {g:5.1f}: " + "  ".join(f"{m} {r:.3f}" for m, r in v.items()))
mc = np.sqrt(0.05 * 0.95 / reps_s)
methods = list(sizes[0.0])
# the F test and four permutation schemes depend on y only through (I - M_Z) y, which does not
# involve gamma: with the same random numbers their rejections are identical for every gamma
for m in ["F test", "Draper-Stoneman", "Kennedy", "Freedman-Lane", "ter Braak"]:
    assert len({sizes[g][m] for g in gammas}) == 1
# the raw permutation (Manly) is not: its size grows with the nuisance effect
manly = [sizes[g]["Manly"] for g in gammas]
assert manly[-1] > 0.05 + 4 * mc and manly[-1] > manly[0] + 3 * mc
assert abs(sizes[0.0]["F test"] - 0.05) < 4 * mc
assert abs(sizes[0.0]["Freedman-Lane"] - 0.05) < 4 * mc + 0.005
assert abs(sizes[0.0]["ter Braak"] - 0.05) < 4 * mc
assert sizes[0.0]["Kennedy"] > sizes[0.0]["Freedman-Lane"]
gen.int("reps", reps_s)
gen.int("B", B_s)
gen.num("mc", mc, 4)
for g, v in sizes.items():
    for m, r in v.items():
        gen.num(f"size:{int(g)}:{m.replace(' ', '').replace('-', '')}", r, 3)
gen.write()

# ---- figure --------------------------------------------------------------------------------------
use_book_style()
fig, ax = plt.subplots(figsize=(4.8, 2.7))
styles = {"F test": (COLORS["ink"], "o"), "Manly": (COLORS["second"], "s"),
          "Draper-Stoneman": (COLORS["muted"], "v"), "Kennedy": (COLORS["thread"], "D"),
          "Freedman-Lane": (COLORS["accent"], "^"), "ter Braak": (COLORS["third"], "x")}
xpos = np.arange(len(gammas))
for m in methods:
    col, mk = styles[m]
    ax.plot(xpos, [sizes[g][m] for g in gammas], color=col, marker=mk, markersize=4, linewidth=0.9, label=m)
ax.axhspan(0.05 - 2 * mc, 0.05 + 2 * mc, color=COLORS["grid"], alpha=0.7, linewidth=0)
ax.axhline(0.05, color=COLORS["ink"], linewidth=0.5, linestyle=":")
ax.set_xticks(xpos)
ax.set_xticklabels([f"{g:g}" for g in gammas])
ax.set_xlabel(r"nuisance coefficient $\gamma$")
ax.set_ylabel("rejection rate")
ax.legend(frameon=False, fontsize=7, ncol=2, loc="upper left")
ax.set_ylim(0.03, 0.095)
fig.tight_layout()
fig.savefig(figure_path("ch23", "freedman_lane_size"))
