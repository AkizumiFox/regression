"""Chapter 22, Section 1: Engel's food expenditure data on the raw and the log-log scale.

235 Belgian working-class households (Engel 1857), annual income and food expenditure.
Public-domain data shipped with statsmodels (statsmodels.datasets.engel).
On the raw scale the residual spread grows with income and the residuals are heavy-tailed; on the
log-log scale both defects largely disappear and the slope is an income elasticity.
"""
import matplotlib.pyplot as plt
import numpy as np
import statsmodels.api as sm
from scipy import stats

from regbook import COLORS, Generated, figure_path, use_book_style

# <<fits>>
df = sm.datasets.engel.load_pandas().data
income, food = df["income"].to_numpy(), df["foodexp"].to_numpy()
n = len(food)

raw = sm.OLS(food, sm.add_constant(income)).fit()                   # food on income
loglog = sm.OLS(np.log(food), sm.add_constant(np.log(income))).fit() # log food on log income

def spread_ratio(fit):
    """Residual SD among the richest third over that among the poorest third."""
    order = np.argsort(income)
    low, high = order[: n // 3], order[-(n // 3):]
    return fit.resid[high].std(ddof=1) / fit.resid[low].std(ddof=1)

for name, fit in [("raw", raw), ("log-log", loglog)]:
    print(f"{name:8s} slope {fit.params[1]:.4f}  spread ratio {spread_ratio(fit):.2f}"
          f"  excess kurtosis of residuals {stats.kurtosis(fit.resid):.2f}")
# <</fits>>

ratio_raw, ratio_log = spread_ratio(raw), spread_ratio(loglog)
kurt_raw, kurt_log = stats.kurtosis(raw.resid), stats.kurtosis(loglog.resid)
assert ratio_raw > 2.5 and ratio_log < 1.6          # the fan largely closes on the log scale
assert kurt_raw > 5 and abs(kurt_log) < 0.5         # heavy tails on the raw scale, near-normal after
assert 0 < loglog.params[1] < 1                     # Engel's law: elasticity below one

# elasticity of the linear fit at the means, for comparison with the log-log slope
elast_raw = raw.params[1] * income.mean() / food.mean()
ci = loglog.conf_int(0.05)[1]
assert ci[1] < 1

# ---- exercise: R^2 on two scales is not comparable (x = 0, 1, 2 and y = 1, 8, 7) ----------------
x3, y3 = np.array([0.0, 1.0, 2.0]), np.array([1.0, 8.0, 7.0])
X3 = sm.add_constant(x3)
fit_raw3, fit_log3 = sm.OLS(y3, X3).fit(), sm.OLS(np.log(y3), X3).fit()
sse_back3 = np.sum((y3 - np.exp(fit_log3.fittedvalues)) ** 2)
assert fit_log3.rsquared > fit_raw3.rsquared and sse_back3 > fit_raw3.ssr
assert np.isclose(fit_raw3.ssr, 32 / 3)

gen = Generated("ch22", "engel_scale", prefix="scale")
gen.num("r2_raw3", fit_raw3.rsquared, 3)
gen.num("r2_log3", fit_log3.rsquared, 3)
gen.num("sse_raw3", fit_raw3.ssr, 2)
gen.num("sse_back3", sse_back3, 2)
gen.int("n", n)
gen.num("slope_raw", raw.params[1], 4)
gen.num("int_raw", raw.params[0], 1)
gen.num("slope_log", loglog.params[1], 4)
gen.num("int_log", loglog.params[0], 4)
gen.num("se_log", loglog.bse[1], 4)
gen.num("ci_lo", ci[0], 3)
gen.num("ci_hi", ci[1], 3)
gen.num("s_log", np.sqrt(loglog.scale), 4)
gen.num("ratio_raw", ratio_raw, 2)
gen.num("ratio_log", ratio_log, 2)
gen.num("kurt_raw", kurt_raw, 2)
gen.num("kurt_log", kurt_log, 2)
gen.num("elast_raw", elast_raw, 2)
gen.num("inc_min", income.min(), 0)
gen.num("inc_max", income.max(), 0)
gen.write()

# ---- figure: fit and residuals on each scale --------------------------------------------------
use_book_style()
fig, axes = plt.subplots(2, 2, figsize=(5.6, 4.3))
panels = [(income, food, raw, "income", "food expenditure", "(a) raw scale"),
          (np.log(income), np.log(food), loglog, "log income", "log food expenditure", "(c) log-log scale")]
for row, (xv, yv, fit, xl, yl, title) in enumerate(panels):
    ax = axes[row, 0]
    ax.scatter(xv, yv, s=6, color=COLORS["accent"], alpha=0.7, linewidths=0)
    xs = np.linspace(xv.min(), xv.max(), 2)
    ax.plot(xs, fit.params[0] + fit.params[1] * xs, color=COLORS["second"])
    ax.set_xlabel(xl)
    ax.set_ylabel(yl)
    ax.set_title(title)
    ax = axes[row, 1]
    ax.scatter(fit.fittedvalues, fit.resid, s=6, color=COLORS["accent"], alpha=0.7, linewidths=0)
    ax.axhline(0, color=COLORS["muted"], linewidth=0.6)
    ax.set_xlabel("fitted value")
    ax.set_ylabel("residual")
    ax.set_title("(b) raw residuals" if row == 0 else "(d) log-log residuals")
fig.tight_layout()
fig.savefig(figure_path("ch22", "engel_scales"))
