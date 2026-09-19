"""Chapter 15, Section 5: an unbalanced one-way layout with unequal variances.

Synthetic data (fixed seed): setting times (minutes) of a repair mortar made with four admixture
formulations A-D, tested 12, 9, 6 and 4 times. The newer formulations were tested least and are
also the most variable (true standard deviations 0.8, 1.2, 2.0 and 3.0 minutes). Values are
rounded to 0.1 minute, as recorded. The classical F test, Welch's test and the Brown-Forsythe
test disagree; residual plots show why.
"""
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats
from statsmodels.stats.oneway import anova_oneway

from regbook import COLORS, Generated, figure_path, use_book_style

# <<data>>
labels = ["A", "B", "C", "D"]
n_k = np.array([12, 9, 6, 4])
rng = np.random.default_rng(1520)
samples = [np.round(rng.normal(mu, sd, k), 1)
           for mu, sd, k in zip([20.0, 20.6, 21.0, 22.2], [0.8, 1.2, 2.0, 3.0], n_k)]
# <</data>>

# <<summary>>
g, n = len(samples), n_k.sum()
means = np.array([x.mean() for x in samples])
s2_k = np.array([x.var(ddof=1) for x in samples])        # group variances
s2 = np.sum((n_k - 1) * s2_k) / (n - g)                  # pooled within-group variance
grand = means @ n_k / n
ss_b = np.sum(n_k * (means - grand) ** 2)
F = ss_b / (g - 1) / s2
print("means:", means.round(3), " sds:", np.sqrt(s2_k).round(2))
print(f"classical F = {F:.2f}, p = {stats.f.sf(F, g - 1, n - g):.4f}")
# null expectation of MS_between if the variances really are s2_k (the theorem of this section)
print(f"MS_within = {s2:.3f};  estimated null E(MS_between) = "
      f"{np.sum((1 - n_k / n) * s2_k) / (g - 1):.3f}")
# <</summary>>

# <<robust>>
welch = anova_oneway(samples, use_var="unequal")
bf = anova_oneway(samples, use_var="bf")
print(f"Welch: F = {welch.statistic:.2f} on ({g - 1}, {welch.df_denom:.1f}) df, p = {welch.pvalue:.3f}")
print(f"Brown-Forsythe: F = {bf.statistic:.2f} on ({g - 1}, {bf.df2[1]:.1f}) df, p = {bf.pvalue2:.3f}")
lev = stats.levene(*samples, center="median")           # dispersion test on |y - median|
print(f"variance test on absolute deviations from medians: p = {lev.pvalue:.4f}")
# <</robust>>

p_f = stats.f.sf(F, g - 1, n - g)
assert p_f < 0.01 and welch.pvalue > 0.1 and bf.pvalue2 > 0.1
assert np.isclose(F, stats.f_oneway(*samples).statistic)
null_ems = np.sum((1 - n_k / n) * s2_k) / (g - 1)
assert null_ems > s2
# Brown-Forsythe statistic = SS_between / estimated null expectation of SS_between
assert np.isclose(bf.statistic, ss_b / np.sum((1 - n_k / n) * s2_k))
# the unweighted and weighted averages of the group means differ in an unbalanced layout
unweighted = means.mean()
# Levene's test by hand: one-way F on absolute deviations from the group medians
z = [np.abs(x - np.median(x)) for x in samples]
assert np.isclose(lev.statistic, stats.f_oneway(*z).statistic)

# Welch t intervals for D - A, against the pooled interval
d = means[3] - means[0]
se_pool = np.sqrt(s2 * (1 / n_k[3] + 1 / n_k[0]))
se_w = np.sqrt(s2_k[3] / n_k[3] + s2_k[0] / n_k[0])
df_w = se_w ** 4 / ((s2_k[3] / n_k[3]) ** 2 / (n_k[3] - 1) + (s2_k[0] / n_k[0]) ** 2 / (n_k[0] - 1))
ci_pool = d + np.array([-1, 1]) * stats.t.ppf(0.975, n - g) * se_pool
ci_w = d + np.array([-1, 1]) * stats.t.ppf(0.975, df_w) * se_w
assert se_w > se_pool

gen = Generated("ch15", "mortar_trial")
for k in range(g):
    gen.text(f"data{k}", ", ".join(f"{v:.1f}" for v in samples[k]))
    gen.num(f"mean{k}", means[k], 3)
    gen.num(f"sd{k}", np.sqrt(s2_k[k]), 2)
gen.num("grand", grand, 3)
gen.num("unweighted", unweighted, 3)
gen.num("s2", s2, 3)
gen.num("s", np.sqrt(s2), 3)
gen.num("ss_b", ss_b, 3)
gen.num("F", F, 2)
gen.num("p_F", p_f, 4)
gen.num("null_ems", null_ems, 3)
gen.num("F_welch", welch.statistic, 2)
gen.num("df_welch", welch.df_denom, 1)
gen.num("p_welch", welch.pvalue, 3)
gen.num("F_bf", bf.statistic, 2)
gen.num("df_bf", bf.df2[1], 1)
gen.num("p_bf", bf.pvalue2, 3)
gen.num("p_levene", lev.pvalue, 4)
gen.num("d_DA", d, 2)
gen.num("se_pool", se_pool, 3)
gen.num("se_welch", se_w, 3)
gen.num("df_welch_t", df_w, 1)
gen.num("ci_pool_lo", ci_pool[0], 2)
gen.num("ci_pool_hi", ci_pool[1], 2)
gen.num("ci_w_lo", ci_w[0], 2)
gen.num("ci_w_hi", ci_w[1], 2)
gen.write()

# ---- figure: residuals against fitted values, and a normal quantile plot -----------------
use_book_style()
fitted = np.concatenate([np.full(k, mu) for k, mu in zip(n_k, means)])
resid = np.concatenate([x - x.mean() for x in samples])
groups = np.repeat(np.arange(g), n_k)
fig, axes = plt.subplots(1, 2, figsize=(5.8, 2.5))
ax = axes[0]
cols = [COLORS["accent"], COLORS["third"], COLORS["thread"], COLORS["second"]]
for k in range(g):
    sel = groups == k
    ax.scatter(fitted[sel] + np.linspace(-0.03, 0.03, sel.sum()), resid[sel], s=11,
               color=cols[k], linewidths=0, label=f"{labels[k]} ($n$={n_k[k]})")
ax.axhline(0, color=COLORS["muted"], lw=0.7)
ax.set_xlabel("fitted value (group mean, min)")
ax.set_ylabel("residual (min)")
ax.set_title("(a) residuals against fitted values")
ax.legend(frameon=False, fontsize=7, loc="center", bbox_to_anchor=(0.62, 0.5))
ax = axes[1]
std_resid = resid / np.sqrt(s2 * (1 - 1 / n_k[groups]))    # internally studentized
osm, osr = stats.probplot(std_resid, dist="norm", fit=False)
ax.scatter(osm, osr, s=10, color=COLORS["accent"], linewidths=0)
ax.plot([-2.5, 2.5], [-2.5, 2.5], color=COLORS["muted"], lw=0.7)
ax.set_xlabel("normal quantile")
ax.set_ylabel("studentized residual")
ax.set_title("(b) normal quantile plot")
fig.tight_layout()
fig.savefig(figure_path("ch15", "mortar_residuals"))
