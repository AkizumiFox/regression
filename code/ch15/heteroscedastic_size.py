"""Chapter 15, Section 5: the size of three tests of equal means when the variances differ.

Four normal groups with equal means and standard deviations growing geometrically from 1 to R
(R = largest / smallest). Three designs with 40 observations: balanced (10 each); positive
pairing (the largest group has the largest variance); negative pairing (the smallest group has
the largest variance). Tests: the classical F test, Welch's test and the Brown-Forsythe test
for means, all at nominal level 0.05. Simulation with a fixed seed.
"""
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats
from statsmodels.stats.oneway import anova_oneway

from regbook import COLORS, Generated, figure_path, use_book_style

# <<tests>>
def three_tests(ybar, s2, n_k):
    """p-values of the classical F, Welch and Brown-Forsythe tests; arrays of shape (reps, g)."""
    g, n = n_k.size, n_k.sum()
    grand = ybar @ n_k / n
    ss_b = ((ybar - grand[:, None]) ** 2) @ n_k
    ss_w = s2 @ (n_k - 1)
    p_f = stats.f.sf(ss_b / (g - 1) / (ss_w / (n - g)), g - 1, n - g)
    # Welch (1951): weights n_k / s_k^2 and a Satterthwaite-type denominator df
    w = n_k / s2
    mu_w = np.sum(w * ybar, axis=1) / w.sum(axis=1)
    h = (1 - w / w.sum(axis=1)[:, None]) ** 2 / (n_k - 1)
    a = np.sum(w * (ybar - mu_w[:, None]) ** 2, axis=1) / (g - 1)
    b = 1 + 2 * (g - 2) / (g ** 2 - 1) * h.sum(axis=1)
    p_w = stats.f.sf(a / b, g - 1, (g ** 2 - 1) / (3 * h.sum(axis=1)))
    # Brown-Forsythe (1974): SS_between over its null expectation under unequal variances
    d = (1 - n_k / n) * s2
    c = d / d.sum(axis=1)[:, None]
    p_bf = stats.f.sf(ss_b / d.sum(axis=1), g - 1, 1 / np.sum(c ** 2 / (n_k - 1), axis=1))
    return p_f, p_w, p_bf
# <</tests>>

# the formulas agree with statsmodels on one data set
rng = np.random.default_rng(2151)
n_k = np.array([4, 6, 10, 20])
samples = [rng.normal(0, sd, k) for sd, k in zip([3, 2, 1.5, 1], n_k)]
ybar = np.array([[x.mean() for x in samples]])
s2 = np.array([[x.var(ddof=1) for x in samples]])
p_f, p_w, p_bf = three_tests(ybar, s2, n_k.astype(float))
assert np.isclose(p_f[0], stats.f_oneway(*samples).pvalue)
assert np.isclose(p_w[0], anova_oneway(samples, use_var="unequal").pvalue)
assert np.isclose(p_bf[0], anova_oneway(samples, use_var="bf").pvalue2)   # the 1974 df

# <<simulate>>
def size(n_k, sds, reps=40000, seed=0, alpha=0.05):
    """Rejection rates under equal means, from group means and variances drawn exactly."""
    rng = np.random.default_rng(seed)
    n_k, sds = np.asarray(n_k, float), np.asarray(sds, float)
    ybar = rng.normal(0, sds / np.sqrt(n_k), (reps, n_k.size))
    s2 = sds ** 2 * rng.chisquare(n_k - 1, (reps, n_k.size)) / (n_k - 1)
    return [np.mean(p < alpha) for p in three_tests(ybar, s2, n_k)]

designs = {"balanced": [10, 10, 10, 10],
           "positive pairing": [4, 6, 10, 20],      # big groups get the big variances
           "negative pairing": [20, 10, 6, 4]}      # small groups get the big variances
ratios = [1, 1.5, 2, 3, 4]
res = {}
for name, sizes in designs.items():                 # sds grow geometrically from 1 to R
    res[name] = np.array([size(sizes, R ** (np.arange(4) / 3), seed=int(10 * R)) for R in ratios])
    for R, row in zip(ratios, res[name]):
        print(f"{name:17s} R = {R:3.1f}: F, Welch, BF sizes = {np.round(row, 3)}")
# <</simulate>>

# the expected mean squares of the theorem, under H0, for R = 3
def ems_ratio(n_k, var):
    n_k, var = np.asarray(n_k, float), np.asarray(var, float)
    n, g = n_k.sum(), n_k.size
    return (np.sum((1 - n_k / n) * var) / (g - 1)) / (np.sum((n_k - 1) * var) / (n - g))

v3 = 3.0 ** (2 * np.arange(4) / 3)
# the difference of the expected mean squares under H0 (the theorem of Section 15.5, part (c))
for sizes in designs.values():
    nk = np.array(sizes, float)
    n_, g_ = nk.sum(), nk.size
    diff = (np.sum((1 - nk / n_) * v3) / (g_ - 1)) - np.sum((nk - 1) * v3) / (n_ - g_)
    assert np.isclose(diff, g_ * (n_ - 1) / (n_ * (g_ - 1) * (n_ - g_)) * np.sum((n_ / g_ - nk) * v3))
ratio = {name: ems_ratio(sizes, v3) for name, sizes in designs.items()}
assert np.isclose(ratio["balanced"], 1.0)
assert ratio["positive pairing"] < 1 < ratio["negative pairing"]
# with equal variances the F test is exact and Brown-Forsythe close; Welch is a little liberal
for name in designs:
    assert np.all(np.abs(res[name][0, [0, 2]] - 0.05) < 0.006)
# the classical F test: conservative with positive pairing, liberal with negative pairing
assert res["positive pairing"][3, 0] < 0.03 and res["negative pairing"][3, 0] > 0.12
assert abs(res["balanced"][3, 0] - 0.05) < 0.03          # balance: a modest excess only
# Welch stays between 0.045 and 0.075 throughout
welch = np.array([res[k][:, 1] for k in designs])
assert np.all((welch > 0.045) & (welch < 0.075))

gen = Generated("ch15", "heteroscedastic_size")
for key, name in (("bal", "balanced"), ("pos", "positive pairing"), ("neg", "negative pairing")):
    gen.num(f"ems_{key}", ratio[name], 3)
    for j, R in enumerate(ratios):
        for t, test in enumerate(["F", "W", "BF"]):
            gen.num(f"{test}_{key}_{int(10 * R)}", res[name][j, t], 3)
gen.write()

# ---- figure ------------------------------------------------------------------------------
use_book_style()
fig, axes = plt.subplots(1, 3, figsize=(6.0, 2.3), sharey=True)
styles = [("classical F", COLORS["second"], "-", "o"), ("Welch", COLORS["accent"], "-", "s"),
          ("Brown–Forsythe", COLORS["third"], "--", "^")]
for ax, (name, sizes) in zip(axes, designs.items()):
    for t, (label, col, ls, mk) in enumerate(styles):
        ax.plot(ratios, res[name][:, t], color=col, ls=ls, marker=mk, ms=3, label=label)
    ax.axhline(0.05, color=COLORS["muted"], lw=0.7, ls=":")
    ax.set_title(f"{name}\n$n_k$ = {', '.join(map(str, sizes))}", fontsize=8)
    ax.set_xlabel("largest / smallest sd")
    ax.set_xticks(ratios)
axes[0].set_ylabel("rejection rate under $H_0$")
axes[0].legend(frameon=False, fontsize=7, loc="upper left")
fig.tight_layout()
fig.savefig(figure_path("ch15", "heteroscedastic_size"))
