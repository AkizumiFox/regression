"""Chapter 29, Section 5: selection bias and honest inference.

(a) Winner's curse: the largest of ten unbiased estimates of equal means.
(b) Freedman's screening paradox: 50 pure-noise regressors, keep those with p < 0.25, refit.
(c) Coverage after selection: the conditional coverage of a reported z interval, and the
    unconditional coverage of the usual interval for beta_1 after a pretest on beta_2.
(d) Honest intervals for the state data: Scheffe, Bonferroni and the PoSI constant.
(e) The truncated-normal (selective) interval in one dimension.
"""
import itertools

import matplotlib.pyplot as plt
import numpy as np
from scipy import stats
from scipy.integrate import quad
from scipy.optimize import brentq

from regbook import COLORS, Generated, figure_path, use_book_style

gen = Generated("ch29", "selection_bias", prefix="sb")
rng = np.random.default_rng(2912)

# ---- (a) winner's curse ----------------------------------------------------------------------------
Zw = rng.normal(size=(400_000, 10))
emax = Zw.max(axis=1).mean()
exact = quad(lambda z: z * 10 * stats.norm.pdf(z) * stats.norm.cdf(z) ** 9, -np.inf, np.inf)[0]
assert abs(emax - exact) < 0.01
gen.num("emax10", exact, 3)

# ---- (b) Freedman's paradox ----------------------------------------------------------------------------
# <<freedman>>
import numpy as np
from scipy import stats
rng = np.random.default_rng(2913)
n, m = 100, 50


def ols(X, y):
    """Coefficients, t statistics (slopes), overall F p-value; X without intercept."""
    A = np.column_stack([np.ones(len(y)), X])
    b, *_ = np.linalg.lstsq(A, y, rcond=None)
    e = y - A @ b
    df = len(y) - A.shape[1]
    s2 = e @ e / df
    t = b / np.sqrt(s2 * np.diag(np.linalg.inv(A.T @ A)))
    ssr = np.sum((A @ b - y.mean()) ** 2)
    F = (ssr / X.shape[1]) / s2
    return t[1:], df, stats.f.sf(F, X.shape[1], df)


def screen_and_refit():
    X = rng.normal(size=(n, m))
    y = rng.normal(size=n)                                  # no regressor matters
    t, df, _ = ols(X, y)
    keep = np.abs(t) > stats.t.ppf(1 - 0.25 / 2, df)        # first pass: p < 0.25
    t2, df2, p_F = ols(X[:, keep], y)                       # second pass on the survivors
    return keep.sum(), np.sum(np.abs(t2) > stats.t.ppf(0.975, df2)), p_F


kept, signif, p_F = screen_and_refit()
print(f"one data set: {kept} kept, {signif} significant at 5%, overall F p-value {p_F:.4f}")
runs = np.array([screen_and_refit() for _ in range(500)])
print("over 500 data sets: mean kept", runs[:, 0].mean(), " mean significant", runs[:, 1].mean(),
      " Pr(F significant)", np.mean(runs[:, 2] < 0.05))
# <</freedman>>
runs = np.vstack([runs, [screen_and_refit() for _ in range(1500)]])
assert np.mean(runs[:, 2] < 0.05) > 0.5
gen.int("kept1", kept)
gen.int("signif1", signif)
gen.num("pF1", p_F, 4)
gen.num("mean_kept", runs[:, 0].mean(), 1)
gen.num("mean_signif", runs[:, 1].mean(), 1)
gen.num("prob_F", np.mean(runs[:, 2] < 0.05), 2)

# ---- (c) coverage after selection ------------------------------------------------------------------------
# <<coverage>>
z = stats.norm.ppf(0.975)


def conditional_coverage(theta, c=z):
    """Pr(|Z - theta| <= z | |Z| > c) for Z ~ N(theta, 1): the reported interval."""
    lo, hi = theta - z, theta + z
    inside_sel = (stats.norm.cdf(max(hi, c) - theta) - stats.norm.cdf(max(lo, c) - theta)) \
        + (stats.norm.cdf(min(hi, -c) - theta) - stats.norm.cdf(min(lo, -c) - theta))
    return inside_sel / (stats.norm.sf(c - theta) + stats.norm.cdf(-c - theta))


def pretest_coverage(delta, r, reps=200_000):
    """Coverage of the usual 95% interval for beta_1 after a 5% pretest of beta_2 = 0.
    Two standardized regressors with correlation r, sigma known; delta = beta_2 / sd(beta_2 hat)."""
    sd2 = 1 / np.sqrt(1 - r**2)                             # sd of either long-model slope
    C = np.array([[1, -r], [-r, 1]]) / (1 - r**2)           # covariance of (b1, b2) in the long model
    b = rng.normal(size=(reps, 2)) @ np.linalg.cholesky(C).T + [0.0, delta * sd2]
    keep = np.abs(b[:, 1] / sd2) > z                        # the pretest keeps x2
    b1_short = b[:, 0] + r * b[:, 1]                        # short slope (sd 1)
    cover_long = np.abs(b[:, 0]) <= z * sd2                 # true beta_1 = 0
    cover_short = np.abs(b1_short) <= z
    return np.mean(np.where(keep, cover_long, cover_short))


for r, d_min in ((0.5, 2.2), (0.8, 1.8), (0.95, 1.2)):              # the minimizing delta of the figure
    print(f"r = {r}: coverage at delta = 0, {d_min}, 5:",
          ", ".join(f"{pretest_coverage(d, r, 50_000):.3f}" for d in (0.0, d_min, 5.0)))
# <</coverage>>

assert conditional_coverage(0.0) == 0.0
assert abs(conditional_coverage(10.0) - 0.95) < 1e-6
thetas = np.linspace(0, 5, 101)
cc = np.array([conditional_coverage(t) for t in thetas])
# check the conditional coverage by simulation at theta = 1.5
Zc = 1.5 + rng.normal(size=1_000_000)
sel = np.abs(Zc) > z
assert abs(np.mean(np.abs(Zc[sel] - 1.5) <= z) - conditional_coverage(1.5)) < 0.003
gen.num("cc1", conditional_coverage(1.0), 3)
gen.num("cc2", conditional_coverage(2.0), 3)
gen.num("cc3", conditional_coverage(3.0), 3)

deltas = np.linspace(0, 6, 31)
rs = (0.5, 0.8, 0.95)
cov = {r: np.array([pretest_coverage(d, r) for d in deltas]) for r in rs}
for r in rs:
    assert abs(cov[r][0] - 0.95) < 0.01 or cov[r][0] > 0.9
    assert cov[r].min() < 0.95
    gen.num(f"min_{int(r * 100)}", cov[r].min(), 3)
    gen.num(f"argmin_{int(r * 100)}", deltas[cov[r].argmin()], 1)
assert cov[0.95].min() < cov[0.5].min()

# ---- (d) honest multipliers for the state data -----------------------------------------------------
# <<posi>>
import itertools
import statsmodels.api as sm
rng = np.random.default_rng(2914)
data = sm.datasets.statecrime.load_pandas().data.drop(index="District of Columbia")
names = ["hs_grad", "poverty", "single", "white", "urban"]
ys = np.log(data["violent"].to_numpy())
Xs = np.column_stack([np.ones(len(ys)), data[names].to_numpy()])
ns, P = Xs.shape
nu = ns - P

# unit vectors u with  coefficient of x_j in submodel M  =  (u' y) / ||v||,  for every M and j in M
U = []
for k in range(1, 6):
    for M in itertools.combinations(range(1, 6), k):
        XM = Xs[:, [0, *M]]
        V = XM @ np.linalg.inv(XM.T @ XM)                   # column i is v for coefficient i
        for i in range(1, len(M) + 1):
            U.append(V[:, i] / np.linalg.norm(V[:, i]))
U = np.array(U)                                              # 80 unit vectors, all in C(X) minus 1
Qb, _ = np.linalg.qr(Xs)
Qb = Qb[:, 1:]                                               # orthonormal basis of C(X) orthogonal to 1
W = U @ Qb                                                   # coordinates of the u's in that basis
sims = 100_000
Zp = rng.normal(size=(sims, P - 1))
s = np.sqrt(rng.chisquare(nu, size=sims) / nu)
K_posi = np.quantile(np.max(np.abs(Zp @ W.T), axis=1) / s, 0.95)
K_t = stats.t.ppf(0.975, nu)
K_bonf = stats.t.ppf(1 - 0.025 / len(U), nu)
K_sch = np.sqrt((P - 1) * stats.f.ppf(0.95, P - 1, nu))
print(f"{len(U)} coefficients; multipliers: t {K_t:.3f}, PoSI {K_posi:.3f}, "
      f"Bonferroni {K_bonf:.3f}, Scheffe {K_sch:.3f}")
# <</posi>>
assert len(U) == 80 and np.allclose(U @ np.ones(ns), 0)
assert K_t < K_posi < min(K_bonf, K_sch)
gen.int("ncoef", len(U))
gen.int("nu", nu)
gen.num("K_t", K_t, 3)
gen.num("K_posi", K_posi, 2)
gen.num("K_bonf", K_bonf, 3)
gen.num("K_sch", K_sch, 3)

# intervals for the AIC choice (single + urban)
XM = Xs[:, [0, 3, 5]]
bM = np.linalg.lstsq(XM, ys, rcond=None)[0]
s_full = np.sqrt(np.sum((ys - Xs @ np.linalg.lstsq(Xs, ys, rcond=None)[0]) ** 2) / nu)
s_M = np.sqrt(np.sum((ys - XM @ bM) ** 2) / (ns - 3))
se_M = s_M * np.sqrt(np.diag(np.linalg.inv(XM.T @ XM)))
se_full_s = s_full * np.sqrt(np.diag(np.linalg.inv(XM.T @ XM)))
for i, nm in ((1, "single"), (2, "urban")):
    gen.num(f"t_{nm}", bM[i] / se_M[i], 2)
    gen.num(f"tfull_{nm}", bM[i] / se_full_s[i], 2)
t_sin, t_urb = bM[1] / se_full_s[1], bM[2] / se_full_s[2]
assert t_sin > K_bonf and abs(t_urb) < K_t                 # single survives every correction
gen.num("s_full", s_full, 4)
gen.num("s_M", s_M, 4)

# ---- (e) the truncated-normal interval ---------------------------------------------------------------
# <<selective>>
from scipy.optimize import brentq


def trunc_cdf(x, theta, c=z):
    """Distribution function at x of N(theta, 1) conditioned on |Z| > c."""
    den = stats.norm.sf(c - theta) + stats.norm.cdf(-c - theta)
    if x <= -c:
        num = stats.norm.cdf(x - theta)
    elif x < c:
        num = stats.norm.cdf(-c - theta)
    else:
        num = stats.norm.cdf(-c - theta) + stats.norm.cdf(x - theta) - stats.norm.cdf(c - theta)
    return num / den


z_obs = 2.5
lower = brentq(lambda th: trunc_cdf(z_obs, th) - 0.975, -10, z_obs)
upper = brentq(lambda th: trunc_cdf(z_obs, th) - 0.025, -10, 10)
print(f"naive interval [{z_obs - z:.2f}, {z_obs + z:.2f}], selective interval [{lower:.2f}, {upper:.2f}]")
# <</selective>>
assert lower < z_obs - z and upper < z_obs + z + 0.1
# conditional coverage of the selective interval, by simulation at theta = 0.5
th0 = 0.5
Zs_ = th0 + rng.normal(size=400_000)
Zs_ = Zs_[np.abs(Zs_) > z][:4000]
covered = [trunc_cdf(zz, th0) > 0.025 and trunc_cdf(zz, th0) < 0.975 for zz in Zs_]
assert abs(np.mean(covered) - 0.95) < 0.012
gen.num("z_obs", z_obs, 1)
gen.num("sel_lo", lower, 2)
gen.num("sel_hi", upper, 2)
gen.num("naive_lo", z_obs - z, 2)
gen.num("naive_hi", z_obs + z, 2)
gen.write()

# ---- figure ------------------------------------------------------------------------------------------
use_book_style()
fig, axes = plt.subplots(1, 2, figsize=(5.8, 2.7))
ax = axes[0]
ax.plot(thetas, cc, color=COLORS["accent"])
ax.axhline(0.95, color=COLORS["muted"], linestyle="--", linewidth=0.8)
ax.set_xlabel(r"$\theta$")
ax.set_ylabel("coverage given selection")
ax.set_title(r"(a) report only if $|Z|>1.96$")
ax.set_ylim(0, 1)
ax = axes[1]
for r, col in zip(rs, (COLORS["accent"], COLORS["third"], COLORS["second"])):
    ax.plot(deltas, cov[r], color=col, label=f"r = {r}")
ax.axhline(0.95, color=COLORS["muted"], linestyle="--", linewidth=0.8)
ax.set_xlabel(r"$\beta_2$ / sd$(\hat\beta_2)$")
ax.set_ylabel(r"coverage for $\beta_1$")
ax.set_title("(b) after a pretest on $\\beta_2$")
ax.set_ylim(0, 1)
ax.legend(frameon=False, loc="lower right")
fig.tight_layout()
fig.savefig(figure_path("ch29", "post_selection"))
