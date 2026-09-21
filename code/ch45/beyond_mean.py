"""Chapter 45, Section 1: what a fitted mean leaves undescribed.

Engel's 1857 budget survey (235 Belgian working-class households; public domain,
shipped with statsmodels as statsmodels.datasets.engel). A homoscedastic normal
model fitted to these data produces one prediction interval whose width does not
depend on income; the data say the width should roughly triple across the range.
"""
import matplotlib.pyplot as plt
import numpy as np
import statsmodels.api as sm
from scipy import stats

from regbook import COLORS, Generated, figure_path, use_book_style

data = sm.datasets.engel.load_pandas().data

# <<bands>>
income = data["income"].to_numpy()
food = data["foodexp"].to_numpy()
n = len(income)
X = np.column_stack([np.ones(n), income])

beta, *_ = np.linalg.lstsq(X, food, rcond=None)
resid = food - X @ beta
sigma = np.sqrt(resid @ resid / (n - 2))          # one number for every household

# a 95% prediction interval from the homoscedastic normal model
t = stats.t.ppf(0.975, n - 2)
lev = np.einsum("ij,jk,ik->i", X, np.linalg.inv(X.T @ X), X)
half = t * sigma * np.sqrt(1 + lev)
inside = (food >= X @ beta - half) & (food <= X @ beta + half)

# split the households into income thirds and look at the spread in each
edges = np.quantile(income, [1 / 3, 2 / 3])
third = np.digitize(income, edges)
spread = [resid[third == k].std(ddof=1) for k in range(3)]
covered = [inside[third == k].mean() for k in range(3)]

print(f"residual sd by income third: {spread[0]:.1f}, {spread[1]:.1f}, {spread[2]:.1f}")
print(f"inside the 95% band:         {covered[0]:.3f}, {covered[1]:.3f}, {covered[2]:.3f}")
# <</bands>>

assert n == 235
assert spread[2] > 2.5 * spread[0]                # the spread grows strongly with income
assert covered[0] > covered[2]                    # the band is too wide low, too narrow high

# five equal-count income groups, for the figure and for a robust look at shape
groups = np.digitize(income, np.quantile(income, [0.2, 0.4, 0.6, 0.8]))
by_group = [resid[groups == k] for k in range(5)]
centres = [np.median(income[groups == k]) for k in range(5)]


def bowley(r):
    """The quartile (Bowley) coefficient of skewness: 0 for a symmetric sample."""
    q1, q2, q3 = np.quantile(r, [0.25, 0.5, 0.75])
    return (q3 + q1 - 2 * q2) / (q3 - q1)


skew_low, skew_high = bowley(by_group[0]), bowley(by_group[4])
iqr_ratio = np.subtract(*np.quantile(by_group[4], [0.75, 0.25])) / np.subtract(
    *np.quantile(by_group[0], [0.75, 0.25]))

gen = Generated("ch45", "beyond_mean")
gen.int("n", n)
gen.num("slope", beta[1], 4)
gen.num("intercept", beta[0], 2)
gen.num("sigma", sigma, 1)
gen.num("sd_low", spread[0], 1)
gen.num("sd_mid", spread[1], 1)
gen.num("sd_high", spread[2], 1)
gen.num("sd_ratio", spread[2] / spread[0], 2)
gen.num("cov_low", covered[0], 3)
gen.num("cov_mid", covered[1], 3)
gen.num("cov_high", covered[2], 3)
gen.num("cov_all", inside.mean(), 3)
gen.num("iqr_ratio", iqr_ratio, 2)
gen.num("skew_low", skew_low, 2)
gen.num("skew_high", skew_high, 2)
gen.write()

# ---- figure ----------------------------------------------------------------
use_book_style()
fig, axes = plt.subplots(1, 2, figsize=(5.6, 2.5))

ax = axes[0]
order = np.argsort(income)
ax.scatter(income, food, s=8, color=COLORS["accent"], alpha=0.7, linewidths=0)
ax.plot(income[order], (X @ beta)[order], color=COLORS["second"])
ax.fill_between(income[order], (X @ beta - half)[order], (X @ beta + half)[order],
                color=COLORS["second"], alpha=0.12, linewidth=0)
ax.set_xlabel("household income (Belgian francs)")
ax.set_ylabel("food expenditure")
ax.set_title("(a) one line, one interval")

ax = axes[1]
ax.axhline(0, color=COLORS["grid"], linewidth=0.6, zorder=0)
box = ax.boxplot(by_group, positions=range(5), widths=[0.55] * 5, manage_ticks=False,
                 showfliers=False, patch_artist=True)
for patch in box["boxes"]:
    patch.set(facecolor=COLORS["accent"], alpha=0.25, edgecolor=COLORS["ink"], linewidth=0.6)
for key in ("whiskers", "caps", "medians"):
    for line in box[key]:
        line.set(color=COLORS["ink"], linewidth=0.7)
ax.set_xticks(range(5))
ax.set_xticklabels([f"{c:.0f}" for c in centres])
ax.set_xlim(-0.6, 4.6)
ax.set_xlabel("median income of the group")
ax.set_ylabel("residual")
ax.set_title("(b) residuals in five income groups")
fig.tight_layout()
fig.savefig(figure_path("ch45", "mean_is_not_enough"))
