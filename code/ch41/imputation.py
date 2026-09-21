"""Chapter 41, Section 5: multiple imputation and Rubin's rules.

A covariate is missing at random, with a probability that depends on the response
and on the other covariate. Five analyses are compared over 2000 replicates:
complete cases, mean imputation, one regression imputation with noise, and proper
multiple imputation with M = 5 and M = 20. What is recorded is bias, the average
reported standard error against the actual spread, and the coverage of the
nominal 95% interval. A second part demonstrates chained equations when two
variables have gaps, and a third draws the efficiency of a finite M.
"""
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats

from regbook import COLORS, Generated, figure_path, use_book_style

n = 200
reps = 2000
rho = 0.5
beta_true = np.array([1.0, 0.5, 0.6])                 # intercept, x1, x2


def make_data(rng):
    x1 = rng.normal(size=n)
    x2 = rho * x1 + np.sqrt(1 - rho**2) * rng.normal(size=n)
    y = beta_true[0] + beta_true[1] * x1 + beta_true[2] * x2 + rng.normal(size=n)
    eta = -0.3 + 0.8 * x1 + 2.0 * (y - y.mean())      # depends on observed values only: MAR
    seen = rng.uniform(size=n) > 1 / (1 + np.exp(-eta))
    return x1, x2, y, seen


def fit(X, y):
    """Least squares with the coefficient of x2 and its squared standard error."""
    b, *_ = np.linalg.lstsq(X, y, rcond=None)
    dof = len(y) - X.shape[1]
    s2 = np.sum((y - X @ b) ** 2) / dof
    v = s2 * np.linalg.inv(X.T @ X)[2, 2]
    return b[2], v, dof


# <<impute>>
def impute_once(x1, x2, y, seen, rng):
    """One draw from the posterior predictive distribution of the missing x2."""
    Z = np.column_stack([np.ones(n), x1, y])          # imputation model: x2 on x1 and y
    Zo, k = Z[seen], int(seen.sum())
    g, *_ = np.linalg.lstsq(Zo, x2[seen], rcond=None)
    sse = float(np.sum((x2[seen] - Zo @ g) ** 2))
    s2_star = sse / rng.chisquare(k - Z.shape[1])     # draw the variance ...
    cov = s2_star * np.linalg.inv(Zo.T @ Zo)
    g_star = rng.multivariate_normal(g, cov)          # ... then the coefficients
    draw = Z @ g_star + np.sqrt(s2_star) * rng.normal(size=n)
    return np.where(seen, x2, draw)                   # observed values are kept


def rubin(estimates, variances):
    """Combine M analyses: the average, the total variance and the degrees of freedom."""
    M = len(estimates)
    q_bar, u_bar = np.mean(estimates), np.mean(variances)
    b = np.var(estimates, ddof=1)
    total = u_bar + (1 + 1 / M) * b
    gamma = (1 + 1 / M) * b / total                   # fraction of missing information
    nu = (M - 1) / gamma**2
    return q_bar, total, nu


demo_rng = np.random.default_rng(4105)
x1d, x2d, yd, seend = make_data(demo_rng)             # one data set, x2 missing at random
qs, us = [], []
for _ in range(5):                                    # M = 5 completions of it
    xi = impute_once(x1d, x2d, yd, seend, demo_rng)
    qi, ui, _ = fit(np.column_stack([np.ones(n), x1d, xi]), yd)
    qs.append(qi)
    us.append(ui)
q_bar, total, nu = rubin(qs, us)
print(f"M = 5: estimate {q_bar:.3f}, total variance {total:.4f}, "
      f"degrees of freedom {nu:.1f}")
# <</impute>>

assert abs(q_bar - beta_true[2]) < 0.25 and total > np.mean(us)


def one_replicate(rng, M_list=(5, 20)):
    x1, x2, y, seen = make_data(rng)
    out = {}
    Xc = np.column_stack([np.ones(int(seen.sum())), x1[seen], x2[seen]])
    q, u, d = fit(Xc, y[seen])
    out["complete cases"] = (q, u, d)
    x_mean = np.where(seen, x2, x2[seen].mean())
    out["mean imputation"] = fit(np.column_stack([np.ones(n), x1, x_mean]), y)
    rng_single = rng
    x_one = impute_once(x1, x2, y, seen, rng_single)
    out["single imputation"] = fit(np.column_stack([np.ones(n), x1, x_one]), y)
    for M in M_list:
        qs, us = [], []
        for _ in range(M):
            xi = impute_once(x1, x2, y, seen, rng)
            qi, ui, _ = fit(np.column_stack([np.ones(n), x1, xi]), y)
            qs.append(qi)
            us.append(ui)
        out[f"MI, M = {M}"] = rubin(qs, us)
    return out, int(seen.sum())


rng = np.random.default_rng(4106)
methods = ["complete cases", "mean imputation", "single imputation", "MI, M = 5", "MI, M = 20"]
res = {k: [] for k in methods}
kept = []
for _ in range(reps):
    out, k = one_replicate(rng)
    kept.append(k)
    for name in methods:
        res[name].append(out[name])

summary = {}
for name in methods:
    q, u, d = np.array(res[name]).T
    half = stats.t.ppf(0.975, d) * np.sqrt(u)
    cover = np.mean(np.abs(q - beta_true[2]) < half)
    summary[name] = (q.mean() - beta_true[2], np.sqrt(u).mean(), q.std(ddof=1), 100 * cover)
    print(f"{name:18s} bias {summary[name][0]:+.4f}  mean se {summary[name][1]:.4f}"
          f"  sd {summary[name][2]:.4f}  coverage {summary[name][3]:.1f}%")
print(f"complete cases: {np.mean(kept):.0f} of {n} on average")

assert summary["complete cases"][3] < 85           # the interval misses far too often
assert summary["mean imputation"][3] < 60
assert summary["single imputation"][1] < 0.85 * summary["single imputation"][2]
for M in (5, 20):
    assert abs(summary[f"MI, M = {M}"][0]) < 0.02
    assert 92.5 < summary[f"MI, M = {M}"][3] < 97.0

# ---- chained equations with gaps in two variables ----------------------------
# <<chained>>
rng2 = np.random.default_rng(4107)
x1, x2, y, seen2 = make_data(rng2)
seen1 = rng2.uniform(size=n) > 0.25                   # x1 now has gaps of its own
x1_work = np.where(seen1, x1, x1[seen1].mean())       # a crude start
x2_work = np.where(seen2, x2, x2[seen2].mean())
for sweep in range(10):                               # one sweep = one pass over the variables
    Z = np.column_stack([np.ones(n), x2_work, y])     # x1 given the others
    g, *_ = np.linalg.lstsq(Z[seen1], x1[seen1], rcond=None)
    s = np.std(x1[seen1] - Z[seen1] @ g, ddof=3)
    x1_work = np.where(seen1, x1, Z @ g + s * rng2.normal(size=n))
    Z = np.column_stack([np.ones(n), x1_work, y])     # x2 given the others
    g, *_ = np.linalg.lstsq(Z[seen2], x2[seen2], rcond=None)
    s = np.std(x2[seen2] - Z[seen2] @ g, ddof=3)
    x2_work = np.where(seen2, x2, Z @ g + s * rng2.normal(size=n))
b_chained, *_ = np.linalg.lstsq(np.column_stack([np.ones(n), x1_work, x2_work]), y, rcond=None)
print("after chained equations:", np.round(b_chained, 3))
# <</chained>>

assert abs(b_chained[2] - beta_true[2]) < 0.2

gen = Generated("ch41", "imputation")
gen.int("n", n)
gen.int("reps", reps)
gen.int("kept", int(round(np.mean(kept))))
keys = {"complete cases": "cc", "mean imputation": "mean", "single imputation": "single",
        "MI, M = 5": "mi5", "MI, M = 20": "mi20"}
for name, key in keys.items():
    gen.num(f"bias_{key}", summary[name][0], 3)
    gen.num(f"se_{key}", summary[name][1], 3)
    gen.num(f"sd_{key}", summary[name][2], 3)
    gen.num(f"cover_{key}", summary[name][3], 1)
nu5 = np.array([r[2] for r in res["MI, M = 5"]])
gamma5 = np.sqrt(4 / nu5)                             # gamma recovered from nu = (M-1)/gamma^2
gen.num("gamma5", float(gamma5.mean()), 3)
gen.num("nu5", float(np.median(nu5)), 1)
gen.num("infl_g3_m5", np.sqrt(1 + 0.3 / 5), 4)
gen.num("infl_g5_m5", np.sqrt(1 + 0.5 / 5), 4)
gen.num("infl_g5_m20", np.sqrt(1 + 0.5 / 20), 4)
gen.write()

# ---- figure ------------------------------------------------------------------
use_book_style()
fig, axes = plt.subplots(1, 2, figsize=(5.6, 2.3))
ax = axes[0]
ax.bar(np.arange(len(methods)), [summary[m][3] for m in methods],
       color=[COLORS["second"]] * 3 + [COLORS["third"]] * 2, width=0.62, linewidth=0)
ax.axhline(95, color=COLORS["ink"], linewidth=0.7, linestyle="--")
ax.set_xticks(np.arange(len(methods)))
ax.set_xticklabels(["complete\ncases", "mean\nimput.", "single\nimput.",
                    "MI\n$M=5$", "MI\n$M=20$"], fontsize=6.5)
ax.set_ylabel("coverage of the 95% interval (%)")
ax.set_ylim(0, 100)
ax.set_title("(a) coverage")

ax = axes[1]
Ms = np.arange(1, 31)
for g, c in zip([0.1, 0.3, 0.5, 0.8], [COLORS["accent"], COLORS["third"],
                                       COLORS["second"], COLORS["thread"]]):
    ax.plot(Ms, np.sqrt(1 + g / Ms), color=c, label=f"$\\gamma={g}$")
ax.axhline(1.0, color=COLORS["grid"], linewidth=0.6, zorder=0)
ax.set_xlabel("number of imputations $M$")
ax.set_ylabel("standard error inflation")
ax.set_title("(b) the cost of a finite $M$")
ax.legend(frameon=False, fontsize=7)
fig.tight_layout()
fig.savefig(figure_path("ch41", "mi_coverage"))
