"""Chapter 21, Section 2: diagnosing heteroscedasticity in Engel's food expenditure data.

Engel's 1857 survey of 235 Belgian working-class households (public domain,
statsmodels.datasets.engel): annual food expenditure against annual income. The script
draws the residual plots, computes the Breusch-Pagan score statistic, Koenker's studentized
version, White's test and the Goldfeld-Quandt test by hand, and checks them against
statsmodels.
"""
import matplotlib.pyplot as plt
import numpy as np
import statsmodels.api as sm
from scipy import stats
from statsmodels.stats.diagnostic import het_breuschpagan, het_goldfeldquandt, het_white

from regbook import COLORS, Generated, figure_path, use_book_style

# <<fit>>
df = sm.datasets.engel.load_pandas().data
income, food = df["income"].to_numpy(), df["foodexp"].to_numpy()
n = len(food)
X = np.column_stack([np.ones(n), income])
beta, *_ = np.linalg.lstsq(X, food, rcond=None)
fitted = X @ beta
e = food - fitted                                     # OLS residuals
Q, _ = np.linalg.qr(X)
h = np.sum(Q**2, axis=1)
s2 = e @ e / (n - 2)
# internally studentized residuals
r = e / np.sqrt(s2 * (1 - h))
# <</fit>>

# <<tests>>
def explained_ss(v, Z):
    """Explained sum of squares of v on [1, Z], about the mean of v."""
    Z1 = np.column_stack([np.ones(len(v)), Z])
    coef, *_ = np.linalg.lstsq(Z1, v, rcond=None)
    return np.sum((Z1 @ coef - v.mean()) ** 2)


# maximum likelihood estimate under H0
sig2 = e @ e / n
u = e**2 / sig2
# Breusch-Pagan score statistic
bp = explained_ss(u, income) / 2
tss = np.sum((e**2 - sig2) ** 2)
# n R^2 from regressing e^2 on [1, income]
koenker = n * explained_ss(e**2, income) / tss
white = n * explained_ss(e**2, np.column_stack([income, income**2])) / tss
# sample kurtosis of the residuals
kappa = n * np.sum(e**4) / np.sum(e**2) ** 2
print(f"Breusch-Pagan {bp:.1f}, Koenker {koenker:.1f}, "
      f"White {white:.1f} (2 df)")
print(f"ratio BP/Koenker = {bp / koenker:.3f}, "
      f"(kurtosis - 1)/2 = {(kappa - 1) / 2:.3f}")

# Goldfeld-Quandt: sort, drop the middle fifth
order = np.argsort(income)
m = (n - 47) // 2
low, high = order[:m], order[-m:]
sse = []
for g in (low, high):
    b_g = np.linalg.lstsq(X[g], food[g], rcond=None)[0]
    sse.append(np.sum((food[g] - X[g] @ b_g) ** 2))
gq = (sse[1] / (m - 2)) / (sse[0] / (m - 2))
print(f"Goldfeld-Quandt F = {gq:.2f} on ({m - 2}, {m - 2}) df, "
      f"p = {stats.f.sf(gq, m - 2, m - 2):.1e}")
# <</tests>>

# checks against statsmodels
assert np.isclose(bp, het_breuschpagan(e, X, robust=False)[0])
assert np.isclose(koenker, het_breuschpagan(e, X, robust=True)[0])
assert np.isclose(white, het_white(e, X)[0])
assert np.isclose(bp / koenker, (kappa - 1) / 2)
gq_sm = het_goldfeldquandt(food, X, idx=1, split=0.4, drop=0.2, alternative="increasing", result_object=False)
assert np.isclose(gq, gq_sm[0]), (gq, gq_sm)

# slope of log |r_i| on log income: sd proportional to income^slope
lcoef = np.polyfit(np.log(income), np.log(np.abs(r)), 1)

gen = Generated("ch21", "engel_hetero", prefix="engh")
gen.int("n", n)
gen.num("slope", beta[1], 4)
gen.num("bp", bp, 1)
gen.num("koenker", koenker, 1)
gen.num("white", white, 1)
gen.num("kappa", kappa, 2)
gen.num("ratio", bp / koenker, 3)
gen.num("gq", gq, 2)
gen.int("gq_m", m)
gen.int("gq_df", m - 2)
gen.num("gq_p", stats.f.sf(gq, m - 2, m - 2), 1, sci=True)
gen.num("koenker_p", stats.chi2.sf(koenker, 1), 1, sci=True)
gen.num("white_p", stats.chi2.sf(white, 2), 1, sci=True)
gen.num("sd_power", lcoef[0], 2)
gen.num("hmax", h.max(), 3)
gen.write()

# ---- figure ---------------------------------------------------------------------
use_book_style()
fig, axes = plt.subplots(1, 2, figsize=(5.8, 2.4))
ax = axes[0]
ax.scatter(fitted, e, s=6, color=COLORS["accent"], alpha=0.8, linewidths=0)
ax.axhline(0, color=COLORS["grid"], linewidth=0.8, zorder=0)
ax.set_xlabel("fitted food expenditure")
ax.set_ylabel("residual")
ax.set_title("(a) residuals against fitted values")
ax = axes[1]
ax.scatter(income, np.abs(r), s=6, color=COLORS["accent"], alpha=0.8, linewidths=0)
xs = np.linspace(income.min(), income.max(), 50)
ax.plot(xs, np.exp(np.polyval(lcoef, np.log(xs))), color=COLORS["second"])
ax.set_xscale("log")
ax.set_yscale("log")
ax.xaxis.set_minor_formatter(plt.NullFormatter())
ax.set_xticks([500, 1000, 2000, 4000], ["500", "1000", "2000", "4000"])
ax.set_xlabel("income (log scale)")
ax.set_ylabel(r"$|r_i|$ (log scale)")
ax.set_title(f"(b) slope {lcoef[0]:.2f} on the log scale")
fig.tight_layout()
fig.savefig(figure_path("ch21", "engel_residuals"))
