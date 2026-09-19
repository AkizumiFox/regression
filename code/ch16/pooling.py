"""Chapter 16, Section 4: the cost of pooling interaction sums of squares into error.

(1) Pooling after a preliminary test. In a 2 x 3 x 4 layout with m = 2 replicates the error
    has 24 df and the three-factor interaction ABC has 6. The rule "pool ABC into error if
    its F test is not significant at level 0.25, then test A at level 0.05" is compared with
    never and always pooling, when A has no effect, as the ABC noncentrality grows. The
    sums of squares are independent with the distributions of the balanced decomposition,
    so they are drawn directly: SS_A ~ chi2(1), SS_ABC ~ chi2(6, gamma), SSE ~ chi2(24).
(2) Pooling the smallest effects. In an unreplicated 2^4 design with no active effects,
    the five smallest of the 15 absolute effects are pooled as "error" and the largest is
    tested against them with the F(1, 5) critical value.
"""
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats

from regbook import COLORS, Generated, figure_path, use_book_style

alpha, alpha_pre = 0.05, 0.25
df_A, df_ABC, df_E = 1, 6, 24

# <<pretest>>
def rejection_rates(gamma, reps, rng):
    """Size of the level-0.05 test of A under three pooling rules, when A has no effect."""
    ssA = rng.chisquare(df_A, reps)
    ssABC = rng.noncentral_chisquare(df_ABC, gamma, reps) if gamma > 0 else rng.chisquare(df_ABC, reps)
    ssE = rng.chisquare(df_E, reps)
    never = ssA / (ssE / df_E) > stats.f.ppf(1 - alpha, df_A, df_E)
    pooled_F = ssA / ((ssE + ssABC) / (df_E + df_ABC))
    always = pooled_F > stats.f.ppf(1 - alpha, df_A, df_E + df_ABC)
    pool = ssABC / df_ABC / (ssE / df_E) <= stats.f.ppf(1 - alpha_pre, df_ABC, df_E)
    pretest = np.where(pool, always, never)
    return never.mean(), always.mean(), pretest.mean(), pool.mean()


rng = np.random.default_rng(20260919)
gammas = np.array([0, 2, 4, 6, 8, 10, 12, 16, 20, 25, 30])
rates = np.array([rejection_rates(g, 400_000, rng) for g in gammas])
for g, (nv, al, pt, pl) in zip(gammas, rates):
    print(f"gamma {g:3d}: never {nv:.4f}  always {al:.4f}  pre-test {pt:.4f}  (pooled in {pl:.2f})")
# <</pretest>>
se = np.sqrt(alpha * (1 - alpha) / 400_000)
assert np.all(np.abs(rates[:, 0] - alpha) < 5 * se)          # never pooling: exact size
assert abs(rates[0, 1] - alpha) < 5 * se                      # always pooling is exact when gamma = 0
assert rates[0, 2] > alpha + 5 * se                           # the pre-test inflates the size at gamma = 0
assert rates[-1, 1] < 0.02                                    # always pooling: conservative when gamma is large
peak = rates[:, 2].max()

# <<smallest>>
rng2 = np.random.default_rng(16_0404)
reps = 200_000
effects = rng2.normal(size=(reps, 15))                        # 15 null effect estimates, common variance
sq = np.sort(effects ** 2, axis=1)
F_largest = sq[:, -1] / sq[:, :5].mean(axis=1)                # largest against the five smallest
size_largest = np.mean(F_largest > stats.f.ppf(1 - alpha, 1, 5))
pooled_ms = sq[:, :5].mean(axis=1)                           # the "error" mean square; true variance 1
underestimate = 1 / pooled_ms.mean()                          # factor by which it falls short on average
print(f"largest effect declared significant: {size_largest:.3f};"
      f"  pooled mean square averages {pooled_ms.mean():.3f}, an underestimate by a factor {underestimate:.1f}")
# <</smallest>>
assert size_largest > 0.9
assert underestimate > 10

gen = Generated("ch16", "pooling")
gen.num("pre0", rates[0, 2], 4)
gen.num("peak", peak, 4)
gen.int("peak_gamma", int(gammas[np.argmax(rates[:, 2])]))
gen.num("pool0", rates[0, 3], 2)
gen.num("always30", rates[-1, 1], 4)
gen.num("pre30", rates[-1, 2], 4)
gen.num("size_largest", size_largest, 3)
gen.num("underestimate", underestimate, 1)
gen.write()

use_book_style()
fig, ax = plt.subplots(figsize=(4.0, 2.5))
ax.plot(gammas, rates[:, 0], color=COLORS["muted"], marker="o", markersize=3, label="never pool")
ax.plot(gammas, rates[:, 1], color=COLORS["accent"], marker="s", markersize=3, label="always pool")
ax.plot(gammas, rates[:, 2], color=COLORS["second"], marker="^", markersize=3, label="pool if $p > 0.25$")
ax.axhline(alpha, color=COLORS["grid"], linewidth=0.6, zorder=0)
ax.set_xlabel(r"noncentrality $\gamma_{ABC}$")
ax.set_ylabel("rejection rate of true $H_A$")
ax.legend(frameon=False)
fig.tight_layout()
fig.savefig(figure_path("ch16", "pooling_size"))
