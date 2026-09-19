"""Chapter 23, Section 4: permutation tests for a regression slope.

(a) Brownlee's stack loss data (statsmodels.datasets.stackloss, public domain): the simple regression of
    stack loss on acid concentration. The permutation p-value from 19,999 random permutations, and a
    check that the slope, the correlation and the t statistic give identical permutation p-values.
(b) A fixed-seed simulation of the size of the permutation test of a zero slope when the errors are
    uncorrelated with x but their variance depends on x: the raw slope gives an invalid test, the
    slope studentized by a heteroscedasticity-consistent standard error an approximately valid one.
"""
import matplotlib.pyplot as plt
import numpy as np
import statsmodels.api as sm
from scipy import stats

from regbook import COLORS, Generated, figure_path, use_book_style

data = sm.datasets.stackloss.load_pandas().data

# <<stackloss>>
rng = np.random.default_rng(2304)
y = data["STACKLOSS"].to_numpy()
x = data["ACIDCONC"].to_numpy()
n = len(y)
xc = x - x.mean()


def slope_stats(Y):
    """Slope, correlation and t statistic for each row of Y, regressed on the fixed x."""
    Yc = Y - Y.mean(axis=-1, keepdims=True)
    b = Yc @ xc / (xc @ xc)
    r = Yc @ xc / np.sqrt((xc @ xc) * np.sum(Yc ** 2, axis=-1))
    t = r * np.sqrt(n - 2) / np.sqrt(1 - r ** 2)
    return b, r, t


b_obs, r_obs, t_obs = slope_stats(y)
p_t = 2 * stats.t.sf(abs(t_obs), n - 2)

B = 19_999
Yperm = y[np.argsort(rng.random((B, n)), axis=1)]  # B random permutations of y
b_perm, r_perm, t_perm = slope_stats(Yperm)
p_perm = (1 + np.sum(np.abs(b_perm) >= abs(b_obs))) / (B + 1)
print(f"slope {b_obs:.4f}, t = {t_obs:.3f}, t-test p = {p_t:.4f}, permutation p = {p_perm:.4f}")
# <</stackloss>>

# the three statistics order the permutations identically (Proposition bs-simple-equivalence)
tol = 1e-9
same_b = np.abs(b_perm) >= abs(b_obs) - tol
same_r = np.abs(r_perm) >= abs(r_obs) - tol
same_t = np.abs(t_perm) >= abs(t_obs) - tol
assert np.array_equal(same_b, same_r) and np.array_equal(same_r, same_t)
mc_sd = np.sqrt(p_perm * (1 - p_perm) / B)
assert abs(p_perm - p_t) < 0.03                   # close to the normal-theory answer
# number of distinct permutations
n_perms = float(np.prod(np.arange(1, n + 1, dtype=float)))

gen = Generated("ch23", "permutation", prefix="perm")
gen.int("n", n)
gen.int("B", B)
gen.num("slope", b_obs, 4)
gen.num("t", t_obs, 3)
gen.num("r", r_obs, 3)
gen.num("p_t", p_t, 4)
gen.num("p_perm", p_perm, 4)
gen.num("mc_sd", mc_sd, 4)
gen.num("nperms", n_perms, 2, sci=True)

# ---- (b) heteroscedastic null -------------------------------------------------------------
# <<hetero>>
def perm_sizes(n, reps, B, rng, alpha=0.05):
    """Rejection rates of permutation tests of zero slope when y = x * eps (x, eps iid normal)."""
    reject = {"raw": 0, "studentized": 0}
    for _ in range(reps):
        x = rng.normal(size=n)
        y = x * rng.normal(size=n)                # E(y | x) = 0, Var(y | x) = x^2
        Y = np.vstack([y, y[np.argsort(rng.random((B, n)), axis=1)]])   # row 0: the data
        xc = x - x.mean()
        Yc = Y - Y.mean(axis=1, keepdims=True)
        b = Yc @ xc / (xc @ xc)
        e = Yc - b[:, None] * xc
        se_hc0 = np.sqrt((e ** 2) @ xc ** 2) / (xc @ xc)
        for name, T in [("raw", np.abs(b)), ("studentized", np.abs(b / se_hc0))]:
            p = np.mean(T >= T[0])                # (1 + #{T_b >= T_0}) / (B + 1)
            reject[name] += p <= alpha
    return {k: v / reps for k, v in reject.items()}


rng = np.random.default_rng(2305)
for n_sim in (20, 200):
    print(n_sim, perm_sizes(n_sim, reps=400, B=199, rng=rng))
# <</hetero>>

rng = np.random.default_rng(2306)
reps_h, B_h = 4000, 399
sizes = {m: perm_sizes(m, reps_h, B_h, rng) for m in (20, 200)}
for m, v in sizes.items():
    print(f"n = {m}: raw {v['raw']:.3f}, studentized {v['studentized']:.3f}")
# limit of the raw test: sqrt(n) r -> N(0, 3) while its permutation distribution -> N(0, 1)
limit_raw = 2 * stats.norm.sf(stats.norm.ppf(0.975) / np.sqrt(3))
mc = np.sqrt(0.05 * 0.95 / reps_h)
assert abs(sizes[200]["raw"] - limit_raw) < 0.03
assert sizes[20]["raw"] > 0.15
assert abs(sizes[200]["studentized"] - 0.05) < 4 * mc
assert sizes[20]["studentized"] < sizes[20]["raw"] / 2
gen.int("reps_h", reps_h)
gen.int("B_h", B_h)
gen.num("limit_raw", limit_raw, 3)
gen.num("z_over_sqrt3", stats.norm.ppf(0.975) / np.sqrt(3), 4)
for m, v in sizes.items():
    gen.num(f"size:{m}:raw", v["raw"], 3)
    gen.num(f"size:{m}:stud", v["studentized"], 3)
gen.write()

# ---- figure: the permutation distribution of t for the stack loss data --------------------------
use_book_style()
fig, ax = plt.subplots(figsize=(4.4, 2.5))
ax.hist(t_perm, bins=80, density=True, histtype="stepfilled", color=COLORS["grid"],
        edgecolor=COLORS["muted"], linewidth=0.5, label="permutation")
tt = np.linspace(-4.5, 4.5, 300)
ax.plot(tt, stats.t.pdf(tt, n - 2), color=COLORS["accent"], label=f"t({n - 2})")
ax.axvline(t_obs, color=COLORS["second"], linewidth=1.0)
ax.axvline(-t_obs, color=COLORS["second"], linewidth=1.0, linestyle="--")
ax.set_xlabel("t statistic for acid concentration")
ax.set_yticks([])
ax.legend(frameon=False, loc="upper left")
fig.tight_layout()
fig.savefig(figure_path("ch23", "stackloss_permutation"))
