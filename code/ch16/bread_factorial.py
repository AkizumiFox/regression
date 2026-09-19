"""Chapter 16, Sections 2-3: a balanced 3 x 4 factorial with three replicates.

Synthetic data: loaf volume (ml) of bread from three flours (white, wholemeal, rye blend)
proofed at four temperatures (24, 28, 32, 36 C), three loaves per combination, 36 loaves
in all. The script
(1) computes the orthogonal decomposition from the Kronecker projections and checks it
    against statsmodels (all three types of sums of squares agree in a balanced layout);
(2) tests interaction and main effects, and splits the interaction into six product
    contrasts (flour contrast x orthogonal polynomial in temperature);
(3) gives Tukey intervals for the flour means averaged over temperature, simple-effect
    comparisons at each temperature, and Scheffe intervals for interaction contrasts;
(4) draws the interaction plot.
"""
import itertools

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from scipy import stats
from statsmodels.stats.anova import anova_lm

from regbook import COLORS, Generated, figure_path, use_book_style

# ---- data --------------------------------------------------------------------------
flours = ["white", "wholemeal", "rye blend"]
temps = [24, 28, 32, 36]
true_means = np.array([[1880, 2010, 2130, 2190],
                       [1660, 1760, 1810, 1750],
                       [1500, 1540, 1580, 1590]], dtype=float)
a, b, m = 3, 4, 3
data_rng = np.random.default_rng(16_03)
volume = np.round(true_means[:, :, None] + data_rng.normal(scale=55, size=(a, b, m)))

# <<decompose>>
def Jbar(k):
    return np.full((k, k), 1.0 / k)


def Cen(k):
    return np.eye(k) - Jbar(k)


y = volume.ravel()                                    # replicate fastest, then temperature, then flour
P = {"flour": np.kron(np.kron(Cen(a), Jbar(b)), Jbar(m)),
     "temperature": np.kron(np.kron(Jbar(a), Cen(b)), Jbar(m)),
     "interaction": np.kron(np.kron(Cen(a), Cen(b)), Jbar(m)),
     "error": np.kron(np.kron(np.eye(a), np.eye(b)), Cen(m))}
ss = {k: y @ Pk @ y for k, Pk in P.items()}
df = {k: int(round(np.trace(Pk))) for k, Pk in P.items()}
mse = ss["error"] / df["error"]
for k in ["flour", "temperature", "interaction", "error"]:
    F = ss[k] / df[k] / mse
    tail = f"F = {F:7.2f}  p = {stats.f.sf(F, df[k], df['error']):.2g}" if k != "error" else ""
    print(f"{k:12s} df {df[k]:2d}  SS {ss[k]:10.0f}  MS {ss[k] / df[k]:9.0f}  {tail}")
# <</decompose>>

n = a * b * m
F_ = {k: ss[k] / df[k] / mse for k in ["flour", "temperature", "interaction"]}
p_ = {k: stats.f.sf(F_[k], df[k], df["error"]) for k in F_}
assert df == {"flour": 2, "temperature": 3, "interaction": 6, "error": 24}
assert p_["interaction"] < 0.01 and p_["flour"] < 1e-10

# statsmodels, sum-to-zero coding: Types I, II and III coincide in a balanced layout
frame = pd.DataFrame({"v": y, "flour": np.repeat(flours, b * m),
                      "temp": np.tile(np.repeat(temps, m), a)})
fit = smf.ols("v ~ C(flour, Sum) * C(temp, Sum)", frame).fit()
for typ in (1, 2, 3):
    tab = anova_lm(fit, typ=typ)
    assert np.isclose(tab.loc["C(flour, Sum)", "sum_sq"], ss["flour"])
    assert np.isclose(tab.loc["C(flour, Sum):C(temp, Sum)", "sum_sq"], ss["interaction"])
assert np.isclose(fit.ssr, ss["error"])

cell = volume.mean(axis=2)
fmean, tmean, grand = cell.mean(1), cell.mean(0), cell.mean()

# ---- product contrasts ------------------------------------------------------------
# <<contrasts>>
flour_c = {"white vs others": np.array([2.0, -1, -1]), "wholemeal vs rye": np.array([0.0, 1, -1])}
temp_c = {"linear": np.array([-3.0, -1, 1, 3]), "quadratic": np.array([1.0, -1, -1, 1]),
          "cubic": np.array([-1.0, 3, -3, 1])}
rows = []
for (fn, c), (tn, d) in itertools.product(flour_c.items(), temp_c.items()):
    D = np.outer(c, d)                                # an interaction contrast in the cell means
    est = (D * cell).sum()
    ss_c = est ** 2 / ((D ** 2).sum() / m)            # its one-degree-of-freedom sum of squares
    rows.append((fn, tn, est, ss_c, ss_c / mse))
    print(f"{fn:17s} x {tn:9s}  estimate {est:8.1f}  SS {ss_c:8.0f}  F {ss_c / mse:6.2f}")
print(f"sum of the six: {sum(r[3] for r in rows):.0f}   interaction SS: {ss['interaction']:.0f}")
# <</contrasts>>
assert np.isclose(sum(r[3] for r in rows), ss["interaction"])

# temperature polynomial main-effect contrasts (averaged over flour)
temp_ss = {tn: (d @ tmean) ** 2 / ((d ** 2).sum() / (a * m)) for tn, d in temp_c.items()}
assert np.isclose(sum(temp_ss.values()), ss["temperature"])

# ---- intervals ---------------------------------------------------------------------
# <<intervals>>
s = np.sqrt(mse)
nu = df["error"]
q = stats.studentized_range.ppf(0.95, a, nu)
half_tukey = q * s / np.sqrt(b * m)                   # flour means average b*m loaves
for i, k in itertools.combinations(range(a), 2):
    d = fmean[i] - fmean[k]
    print(f"{flours[i]:9s} - {flours[k]:9s}: {d:7.1f}  Tukey [{d - half_tukey:7.1f}, {d + half_tukey:7.1f}]")
t_bonf = stats.t.ppf(1 - 0.05 / (2 * b), nu)          # white - wholemeal at each temperature
half_simple = t_bonf * s * np.sqrt(2 / m)
for j in range(b):
    d = cell[0, j] - cell[1, j]
    print(f"white - wholemeal at {temps[j]} C: {d:6.1f}  Bonferroni [{d - half_simple:6.1f}, {d + half_simple:6.1f}]")
c_scheffe = np.sqrt(df["interaction"] * stats.f.ppf(0.95, df["interaction"], nu))
# <</intervals>>

white_lin = rows[0]
se_lin = s * np.sqrt((np.outer(flour_c["white vs others"], temp_c["linear"]) ** 2).sum() / m)
scheffe_lo = white_lin[2] - c_scheffe * se_lin
scheffe_hi = white_lin[2] + c_scheffe * se_lin
assert scheffe_lo > 0                                 # survives Scheffe's correction
share_lin = white_lin[3] / ss["interaction"]
# the linear contrast (-3, -1, 1, 3) has sum e_j t_j = 40 over the temperatures, and the flour
# contrast (2, -1, -1) is twice (white - mean of the others): estimate / 80 is a slope difference
assert np.isclose(temp_c["linear"] @ np.array(temps), 40)
slopes = np.polyfit(temps, cell.T, 1)[0]
assert np.isclose(white_lin[2] / 80, slopes[0] - slopes[1:].mean())
# the most significant interaction contrast is the table of estimated interaction effects,
# and its F statistic is (a-1)(b-1) times F_AB
gam = cell - fmean[:, None] - tmean[None, :] + grand
ss_best = (gam * cell).sum() ** 2 / ((gam ** 2).sum() / m)
assert np.isclose(ss_best, ss["interaction"])
simple = cell[0] - cell[1]
# weights of the temperatures in the flour main effect: equal weights (balanced)
# versus weights proportional to (1, 1, 2, 4), emphasizing warm proofing
w = np.array([1.0, 1, 2, 4]) / 8

gen = Generated("ch16", "bread_factorial")
for i in range(a):
    for j in range(b):
        gen.num(f"cell{i + 1}{j + 1}", cell[i, j], 1)
    gen.num(f"fmean{i + 1}", fmean[i], 1)
for j in range(b):
    gen.num(f"tmean{j + 1}", tmean[j], 1)
    gen.num(f"simple{j + 1}", simple[j], 1)
gen.num("grand", grand, 1)
for k in ss:
    gen.num(f"ss_{k}", ss[k], 0)
    gen.num(f"ms_{k}", ss[k] / df[k], 0)
for k in F_:
    gen.num(f"F_{k}", F_[k], 2)
    mant, exp = f"{p_[k]:.1e}".split("e")
    gen.text(f"p_{k}", f"{mant}\\times 10^{{{int(exp)}}}")
gen.num("s", s, 1)
gen.num("q", q, 3)
gen.num("half_tukey", half_tukey, 1)
gen.num("t_bonf", t_bonf, 3)
gen.num("half_simple", half_simple, 1)
gen.num("c_scheffe", c_scheffe, 3)
gen.num("se_lin", se_lin, 1)
gen.num("scheffe_lo", scheffe_lo, 1)
gen.num("scheffe_hi", scheffe_hi, 1)
gen.num("share_lin", share_lin, 3)
for idx, r in enumerate(rows):
    gen.num(f"c{idx}_est", r[2], 1)
    gen.num(f"c{idx}_ss", r[3], 0)
    gen.num(f"c{idx}_F", r[4], 2)
for tn, v in temp_ss.items():
    gen.num(f"temp_{tn}", v, 0)
gen.num("equal_weight_diff", simple.mean(), 1)
for i, k in itertools.combinations(range(a), 2):
    d = fmean[i] - fmean[k]
    gen.num(f"diff{i + 1}{k + 1}", d, 1)
    gen.num(f"diff{i + 1}{k + 1}_lo", d - half_tukey, 1)
    gen.num(f"diff{i + 1}{k + 1}_hi", d + half_tukey, 1)
for j in range(b):
    gen.num(f"simple{j + 1}_lo", simple[j] - half_simple, 1)
    gen.num(f"simple{j + 1}_hi", simple[j] + half_simple, 1)
gen.num("slope_diff", white_lin[2] / 80, 2)
gen.num("slope_lo", scheffe_lo / 80, 2)
gen.num("slope_hi", scheffe_hi / 80, 2)
gen.num("se_cell", s / np.sqrt(m), 1)
gen.num("warm_weight_diff", w @ simple, 1)
gen.write()

# ---- interaction plot ----------------------------------------------------------------
use_book_style()
fig, ax = plt.subplots(figsize=(4.2, 2.7))
colours = [COLORS["accent"], COLORS["second"], COLORS["third"]]
marks = ["o", "s", "^"]
for i in range(a):
    ax.plot(temps, cell[i], color=colours[i], marker=marks[i], markersize=4, label=flours[i])
    for j in range(b):
        ax.scatter([temps[j] + (i - 1) * 0.3] * m, volume[i, j], s=6, color=colours[i], alpha=0.45,
                   linewidths=0)
ax.set_xticks(temps)
ax.set_xlabel("proofing temperature (°C)")
ax.set_ylabel("loaf volume (ml)")
ax.legend(frameon=False, loc="upper left", ncol=3)
ax.set_ylim(1400, 2380)
fig.tight_layout()
fig.savefig(figure_path("ch16", "bread_interaction"))
