"""Chapter 13: the running one-way example.

Where respondents to the 1996 American National Election Study place Bill Clinton on a
seven-point left-right scale (1 = extremely liberal, 7 = extremely conservative), by
the respondent's education, coded 1 (least) to 7 (most). Public-domain data shipped with
statsmodels (statsmodels.datasets.anes96). Used in Sections 13.2, 13.3, 13.4, 13.5 and 13.6.
"""
import itertools

import matplotlib.pyplot as plt
import numpy as np
import statsmodels.api as sm
from scipy import stats
from statsmodels.stats.multitest import multipletests

from regbook import COLORS, Generated, figure_path, use_book_style

anes = sm.datasets.anes96.load_pandas().data

# <<data>>
y = anes["ClinLR"].to_numpy()
level = anes["educ"].to_numpy().astype(int) - 1       # 0, ..., 6
g = 7
n_k = np.bincount(level).astype(float)
means = np.bincount(level, weights=y) / n_k
n = len(y)
nu = n - g                                            # error degrees of freedom
s2 = np.sum((y - means[level]) ** 2) / nu             # pooled variance estimate
s = np.sqrt(s2)
between = np.sum(n_k * (means - y.mean()) ** 2)
F = between / (g - 1) / s2
print("group sizes", n_k.astype(int))
print("means      ", means.round(3))
print(f"s = {s:.4f} on {nu} df;  F = {F:.3f},  p = {stats.f.sf(F, g - 1, nu):.2g}")
# <</data>>

alpha = 0.05

# <<scheffe>>
c_sch = np.sqrt((g - 1) * stats.f.ppf(1 - alpha, g - 1, nu))   # Scheffe multiplier, q = g - 1

def contrast(c):
    """Estimate, standard error and sum of squares of the contrast sum_k c_k mu_k."""
    est = c @ means
    v = np.sum(c ** 2 / n_k)
    return est, s * np.sqrt(v), est ** 2 / v

c_max = n_k * (means - y.mean())                      # the most significant contrast
c_lohi = np.array([0.5, 0.5, 0, 0, 0, -0.5, -0.5])    # levels 1-2 against levels 6-7
for name, c in (("maximal", c_max), ("1-2 vs 6-7", c_lohi)):
    est, se, ss = contrast(c)
    print(f"{name:11s} estimate {est:8.3f}  t = {est / se:6.2f}  SS = {ss:8.2f}"
          f"  Scheffe interval [{est - c_sch * se:.3f}, {est + c_sch * se:.3f}]")
print(f"between-groups SS = {between:.2f};  Scheffe multiplier = {c_sch:.3f}")
# <</scheffe>>

est_max, se_max, ss_max = contrast(c_max)
assert np.isclose(ss_max, between)                    # the maximal contrast carries it all
assert np.isclose((est_max / se_max) ** 2, (g - 1) * F)
# every contrast's SS is at most the between SS (random contrasts)
rng = np.random.default_rng(20260919)
C = rng.normal(size=(5000, g))
C -= C.mean(axis=1, keepdims=True)
ss_rand = (C @ means) ** 2 / (C ** 2 / n_k).sum(axis=1)
assert ss_rand.max() <= between * (1 + 1e-12)
est_lh, se_lh, ss_lh = contrast(c_lohi)

# <<pairwise>>
pairs = list(itertools.combinations(range(g), 2))
m = len(pairs)                                        # 21 comparisons
diff = np.array([means[k] - means[l] for k, l in pairs])
se_d = np.array([s * np.sqrt(1 / n_k[k] + 1 / n_k[l]) for k, l in pairs])
t = diff / se_d
pval = 2 * stats.t.sf(np.abs(t), nu)

mult = {                                              # half-width = multiplier x se
    "unadjusted": stats.t.ppf(1 - alpha / 2, nu),
    "Bonferroni": stats.t.ppf(1 - alpha / (2 * m), nu),
    "Sidak": stats.t.ppf((1 + (1 - alpha) ** (1 / m)) / 2, nu),
    "Tukey-Kramer": stats.studentized_range.ppf(1 - alpha, g, nu) / np.sqrt(2),
    "Scheffe": c_sch,
}
for name, c in mult.items():
    print(f"{name:13s} multiplier {c:.3f}   pairs declared different: {np.sum(np.abs(t) > c)}")
# <</pairwise>>

# <<stepwise>>
def holm(p, alpha):
    """Holm's step-down procedure: reject the i smallest p-values, i maximal with
    p_(j) <= alpha / (m - j + 1) for all j <= i."""
    m = len(p)
    order = np.argsort(p)
    reject = np.zeros(m, dtype=bool)
    for j, idx in enumerate(order):                   # j = 0, 1, ...
        if p[idx] > alpha / (m - j):
            break
        reject[idx] = True
    return reject

def benjamini_hochberg(p, q):
    """Reject the k smallest p-values, k the largest index with p_(k) <= k q / m."""
    m = len(p)
    order = np.argsort(p)
    below = np.nonzero(p[order] <= q * np.arange(1, m + 1) / m)[0]
    reject = np.zeros(m, dtype=bool)
    if below.size:
        reject[order[: below[-1] + 1]] = True
    return reject

print("Holm rejects", holm(pval, alpha).sum(), "  BH (q = 0.05) rejects",
      benjamini_hochberg(pval, alpha).sum())
# <</stepwise>>

counts = {name: int(np.sum(np.abs(t) > c)) for name, c in mult.items()}
counts["Holm"] = int(holm(pval, alpha).sum())
counts["BH"] = int(benjamini_hochberg(pval, alpha).sum())
# agree with statsmodels
assert counts["Holm"] == multipletests(pval, alpha, "holm")[0].sum()
assert counts["BH"] == multipletests(pval, alpha, "fdr_bh")[0].sum()
assert np.all(multipletests(pval, alpha, "holm")[0] == holm(pval, alpha))
assert np.all(multipletests(pval, alpha, "fdr_bh")[0] == benjamini_hochberg(pval, alpha))
# ordering of the multipliers and of the rejection sets
c = mult
assert c["unadjusted"] < c["Tukey-Kramer"] < c["Sidak"] < c["Bonferroni"] < c["Scheffe"]
assert counts["Scheffe"] <= counts["Bonferroni"] <= counts["Sidak"] <= counts["Tukey-Kramer"]
assert counts["Bonferroni"] <= counts["Holm"] <= counts["BH"] <= counts["unadjusted"]
# Holm's rejections contain Bonferroni's
assert np.all(holm(pval, alpha) >= (pval <= alpha / m))

# Holm steps and BH thresholds, for the text
order = np.argsort(pval)
p_sorted = pval[order]
k_bh = counts["BH"]
k_holm = counts["Holm"]

gen = Generated("ch13", "education", prefix="edu")
gen.int("n", n)
gen.int("nu", nu)
gen.int("m", m)
gen.text("sizes", ",\\ ".join(str(int(v)) for v in n_k))
gen.text("means", ",\\ ".join(f"{v:.2f}" for v in means))
gen.num("s", s, 3)
gen.num("F", F, 2)
gen.num("Fp", stats.f.sf(F, g - 1, nu), 1, sci=True)
gen.num("between", between, 2)
gen.num("csch", c_sch, 3)
gen.num("estlh", est_lh, 3)
gen.num("selh", se_lh, 3)
gen.num("tlh", est_lh / se_lh, 2)
gen.num("lolh", est_lh - c_sch * se_lh, 3)
gen.num("hilh", est_lh + c_sch * se_lh, 3)
gen.num("tmax", est_max / se_max, 2)
for name, key in (("unadjusted", "t"), ("Bonferroni", "bon"), ("Sidak", "sid"),
                  ("Tukey-Kramer", "tk"), ("Scheffe", "sch")):
    gen.num(f"c{key}", mult[name], 3)
    gen.int(f"n{key}", counts[name])
gen.int("nholm", counts["Holm"])
gen.int("nbh", counts["BH"])
gen.num("qtk", stats.studentized_range.ppf(1 - alpha, g, nu), 3)
# the two pairs for the width table: smallest and largest groups
small = pairs.index((0, 1))                           # levels 1 and 2 (13 and 52)
large = pairs.index((2, 5))                           # levels 3 and 6 (248 and 227)
gen.num("sesmall", se_d[small], 3)
gen.num("selarge", se_d[large], 3)
for name, key in (("unadjusted", "t"), ("Bonferroni", "bon"), ("Sidak", "sid"),
                  ("Tukey-Kramer", "tk"), ("Scheffe", "sch")):
    gen.num(f"hwsmall{key}", mult[name] * se_d[small], 3)
    gen.num(f"hwlarge{key}", mult[name] * se_d[large], 3)
    gen.num(f"ratio{key}", mult[name] / mult["unadjusted"], 2)
gen.num("pmin", p_sorted[0], 1, sci=True)
gen.num("pholmlast", p_sorted[k_holm - 1], 4)
gen.num("pholmnext", p_sorted[k_holm], 4)
gen.num("holmnextcut", alpha / (m - k_holm), 4)
gen.num("pbhlast", p_sorted[k_bh - 1], 4)
gen.num("bhlastcut", alpha * k_bh / m, 4)
gen.num("pbhnext", p_sorted[k_bh], 4)
gen.num("bhnextcut", alpha * (k_bh + 1) / m, 4)
gen.num("bonfcut", alpha / m, 5)
gen.write()

# ---- Tukey-Kramer intervals for all 21 differences ------------------------------
use_book_style()
fig, ax = plt.subplots(figsize=(4.8, 4.0))
hw = mult["Tukey-Kramer"] * se_d
hw_t = mult["unadjusted"] * se_d
ypos = np.arange(m)[::-1]
for i in range(m):
    excl = abs(diff[i]) > hw[i]
    col = COLORS["accent"] if excl else COLORS["muted"]
    ax.plot([diff[i] - hw_t[i], diff[i] + hw_t[i]], [ypos[i], ypos[i]], color=col,
            linewidth=3.0, alpha=0.25, solid_capstyle="butt")
    ax.plot([diff[i] - hw[i], diff[i] + hw[i]], [ypos[i], ypos[i]], color=col, linewidth=1.0)
    ax.plot(diff[i], ypos[i], "o", color=col, markersize=2.5)
ax.axvline(0, color=COLORS["ink"], linewidth=0.6)
ax.set_yticks(ypos)
ax.set_yticklabels([f"{k + 1} − {l + 1}" for k, l in pairs], fontsize=7)
ax.set_xlabel("difference in mean placement of Clinton")
ax.set_ylabel("education levels compared")
fig.savefig(figure_path("ch13", "tukey_education"))

# claims of the Tukey-Kramer example: the significant pairs
sig = {pairs[i] for i in range(m) if abs(t[i]) > mult["Tukey-Kramer"]}
assert sig == {(0, 4), (0, 5), (0, 6), (1, 3), (1, 4), (1, 5), (1, 6), (2, 5), (2, 6)}
assert np.argmax(means) == 0 and n_k.argmin() == 0   # smallest group, highest mean
