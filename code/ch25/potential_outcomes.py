"""Chapter 25, Section 2: the randomization distribution of the difference in means.

A synthetic finite population of 20 units with heterogeneous unit effects. Every one of
the C(20, 8) completely randomized assignments of 8 units to treatment is enumerated,
so the mean and variance of the difference in means are computed exactly and compared
with Neyman's formulas (Section 25.2).
"""
import matplotlib.pyplot as plt
import numpy as np

from regbook import COLORS, Generated, figure_path, use_book_style

# <<population>>
import itertools
import numpy as np

rng = np.random.default_rng(2502)
n, n1 = 20, 8
n0 = n - n1
y0 = np.round(rng.normal(10, 3, size=n), 1)                      # potential outcome, control
y1 = np.round(y0 + 2 + 0.5 * (y0 - 10) + rng.normal(0, 1, n), 1)  # potential outcome, treated
tau = np.mean(y1 - y0)                                            # the sample average effect
S1, S0, St = np.var(y1, ddof=1), np.var(y0, ddof=1), np.var(y1 - y0, ddof=1)
print(f"tau = {tau:.3f};  S1^2 = {S1:.3f}, S0^2 = {S0:.3f}, S_tau^2 = {St:.3f}")
# <</population>>

# <<enumerate>>
idx = np.array(list(itertools.combinations(range(n), n1)))       # every possible treated set
treated = np.zeros((len(idx), n), dtype=bool)
treated[np.arange(len(idx))[:, None], idx] = True                 # one row per assignment
m1, m0 = (treated @ y1) / n1, (~treated @ y0) / n0                # arm means, one row per assignment
diff = m1 - m0                                                    # difference in means
s1 = (treated @ y1 ** 2 - n1 * m1 ** 2) / (n1 - 1)                # within-arm variances
s0 = (~treated @ y0 ** 2 - n0 * m0 ** 2) / (n0 - 1)
V_hat = s1 / n1 + s0 / n0                                         # Neyman's variance estimate
V_neyman = S1 / n1 + S0 / n0 - St / n
print(f"{len(diff)} assignments")
print(f"mean of difference in means {diff.mean():.3f}   (tau = {tau:.3f})")
print(f"variance over assignments   {diff.var():.4f}  (Neyman's formula {V_neyman:.4f})")
print(f"mean of variance estimate   {V_hat.mean():.4f}  (excess S_tau^2/n = {St / n:.4f})")
# <</enumerate>>

from math import comb
assert len(diff) == comb(n, n1)
assert np.isclose(diff.mean(), tau, atol=1e-12)
assert np.isclose(diff.var(), V_neyman, rtol=1e-10)
assert np.isclose(V_hat.mean() - V_neyman, St / n, rtol=1e-10)
# the least squares coefficient of T equals the difference in means (one assignment)
t = treated[0]
yobs = np.where(t, y1, y0)
coef = np.linalg.lstsq(np.column_stack([np.ones(n), t]), yobs, rcond=None)[0]
assert np.isclose(coef[1], diff[0])
assert np.isclose(s1[0], np.var(y1[t], ddof=1)) and np.isclose(s0[0], np.var(y0[~t], ddof=1))
# pooled-variance standard error of ordinary least squares, averaged over assignments
sp2 = ((n1 - 1) * s1 + (n0 - 1) * s0) / (n - 2)
V_ols = sp2 * (1 / n1 + 1 / n0)
assert V_ols.mean() < V_neyman                     # pooled OLS variance understates on average
coverage = np.mean(np.abs(diff - tau) <= 1.959964 * np.sqrt(V_hat))

gen = Generated("ch25", "potential_outcomes")
gen.num("tau", tau, 3)
gen.num("Sone", S1, 3)
gen.num("Szero", S0, 3)
gen.num("Stau", St, 3)
gen.int("count", len(diff))
gen.num("meandiff", diff.mean(), 3)
gen.num("vardiff", diff.var(), 4)
gen.num("vneyman", V_neyman, 4)
gen.num("meanvhat", V_hat.mean(), 4)
gen.num("excess", St / n, 4)
gen.num("meanvols", V_ols.mean(), 4)
gen.num("coverage", coverage, 3)
gen.write()

use_book_style()
fig, ax = plt.subplots(figsize=(4.6, 2.5))
ax.hist(diff, bins=60, density=True, color=COLORS["accent"], alpha=0.75)
grid = np.linspace(diff.min(), diff.max(), 300)
sd = np.sqrt(V_neyman)
ax.plot(grid, np.exp(-0.5 * ((grid - tau) / sd) ** 2) / (sd * np.sqrt(2 * np.pi)),
        color=COLORS["second"], label="normal, Neyman variance")
ax.axvline(tau, color=COLORS["ink"], linewidth=0.8, linestyle="--")
ax.set_xlabel("difference in means over all assignments")
ax.set_ylabel("density")
ax.set_ylim(0, 1.3 * ax.get_ylim()[1])
ax.legend(frameon=False, loc="upper left")
fig.tight_layout()
fig.savefig(figure_path("ch25", "randomization_distribution"))
