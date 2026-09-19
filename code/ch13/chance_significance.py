"""Chapter 13, Section 1: how often something is 'significant' by chance.

All g group means are equal. Every pair of groups is compared with an unadjusted
two-sided t test at level 0.05 using the pooled variance estimate. The expected
number of false rejections is exactly 0.05 times the number of pairs; the chance of
at least one is far larger than 0.05, though smaller than for independent tests.
"""
import itertools

import matplotlib.pyplot as plt
import numpy as np
from scipy import stats

from regbook import COLORS, Generated, figure_path, use_book_style

alpha, reps, m_rep = 0.05, 20_000, 5                  # 5 observations per group

# <<simulate>>
def any_significant(g, m_rep, reps, alpha):
    """Share of data sets with at least one unadjusted pairwise rejection, and mean count."""
    rng = np.random.default_rng([20260919, g])        # one fixed stream per g
    nu = g * (m_rep - 1)                              # error degrees of freedom
    y = rng.normal(size=(reps, g, m_rep))             # all group means equal
    means = y.mean(axis=2)
    s2 = ((y - means[:, :, None]) ** 2).sum(axis=(1, 2)) / nu
    se = np.sqrt(s2 * 2 / m_rep)                      # standard error of a difference
    crit = stats.t.ppf(1 - alpha / 2, nu)
    count = np.zeros(reps)
    for k, l in itertools.combinations(range(g), 2):
        count += np.abs(means[:, k] - means[:, l]) / se > crit
    return np.mean(count > 0), np.mean(count)

for g in (3, 5, 10):
    pairs = g * (g - 1) // 2
    fwer, mean_count = any_significant(g, m_rep, reps, alpha)
    print(f"g = {g:2d}, {pairs:2d} pairs: P(some rejection) = {fwer:.3f}, "
          f"mean number = {mean_count:.2f} (alpha x pairs = {alpha * pairs:.2f}), "
          f"independent tests would give {1 - (1 - alpha) ** pairs:.3f}")
# <</simulate>>

gen = Generated("ch13", "chance_significance", prefix="chance")
gen.int("reps", reps)
gen.int("mrep", m_rep)
gs = np.arange(2, 13)
fwers, counts = [], []
for g in gs:
    f, c = any_significant(g, m_rep, reps, alpha)
    pairs = g * (g - 1) // 2
    fwers.append(f)
    counts.append(c)
    se = np.sqrt(f * (1 - f) / reps)
    # a single test (g = 2) has exact level alpha
    if g == 2:
        assert abs(f - alpha) < 4 * se
    # the mean number of false rejections is exactly alpha * pairs (each test is exact)
    assert abs(c - alpha * pairs) < 0.05 * alpha * pairs + 0.01
    # dependence through shared means and s makes the rate below the independent-test value
    if g >= 4:
        assert f < 1 - (1 - alpha) ** pairs
    # and it is always at most the Bonferroni bound
    assert f <= min(1.0, alpha * pairs) + 4 * se
    if g in (3, 5, 7, 10):
        gen.num(f"fwer{g}", f, 3)
        gen.num(f"indep{g}", 1 - (1 - alpha) ** pairs, 3)
        gen.num(f"count{g}", c, 2)
        gen.int(f"pairs{g}", pairs)
gen.num("indep10tests", 1 - (1 - alpha) ** 10, 3)
gen.num("indep100tests", 1 - (1 - alpha) ** 100, 3)
gen.write()

use_book_style()
fig, ax = plt.subplots(figsize=(4.6, 2.6))
pairs = gs * (gs - 1) / 2
grid = np.linspace(1, pairs.max(), 200)
ax.plot(grid, np.minimum(1, alpha * grid), color=COLORS["muted"], linestyle="--",
        label=r"Bonferroni bound $m\alpha$")
ax.plot(grid, 1 - (1 - alpha) ** grid, color=COLORS["second"],
        label=r"independent tests $1-(1-\alpha)^m$")
ax.plot(pairs, fwers, "o", color=COLORS["accent"], markersize=4,
        label="all pairs of $g$ groups (simulated)")
for g, m, f in zip(gs, pairs, fwers):
    if g in (3, 5, 8, 12):
        ax.annotate(f"$g={g}$", (m, f), textcoords="offset points", xytext=(4, -10), fontsize=7,
                    color=COLORS["accent"])
ax.axhline(alpha, color=COLORS["grid"], linewidth=0.8, zorder=0)
ax.set_xlabel("number of tests $m$")
ax.set_ylabel("P(at least one rejection)")
ax.set_ylim(0, 1.02)
ax.legend(frameon=False, loc="lower right")
fig.savefig(figure_path("ch13", "chance_significance"))
