"""Chapter 15, Sections 1 and 3: the nitrogen trial, a balanced one-way layout.

Synthetic data (fixed seed): five nitrogen rates, six plots per rate, completely randomized.
The true means follow a Mitscherlich response curve 6.8 - 2.7 exp(-0.018 N) (t/ha), and the
plot-to-plot standard deviation is 0.5 t/ha. Yields are rounded to two decimals, as recorded.
"""
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats

from regbook import COLORS, Generated, figure_path, use_book_style

# <<data>>
rates = np.array([0, 40, 80, 120, 160])            # kg N per hectare
g, m = len(rates), 6                                # 5 rates, 6 plots per rate
true_means = 6.8 - 2.7 * np.exp(-0.018 * rates)     # unknown to the analyst
rng = np.random.default_rng(1515)
group = np.repeat(np.arange(g), m)                  # plot i received rate group[i]
y = np.round(true_means[group] + rng.normal(0, 0.5, g * m), 2)   # yield, t/ha
n = len(y)
# <</data>>

# <<fit>>
Z = (group[:, None] == np.arange(g)).astype(float)  # cell-means model matrix
n_k = Z.sum(axis=0)
means = Z.T @ y / n_k                               # least squares: the group means
fitted = Z @ means                                  # M y
resid = y - fitted                                  # (I - M) y
s2 = resid @ resid / (n - g)                        # pooled within-group variance
print("group means:", means.round(3))
print(f"s = {np.sqrt(s2):.4f} on {n - g} degrees of freedom")
# <</fit>>

# projection facts: M y is the group means, the residual is orthogonal to C(Z)
M = Z @ np.diag(1 / n_k) @ Z.T
assert np.allclose(M @ y, fitted)
assert np.allclose(Z.T @ resid, 0)
# the effects form [1, Z] has the same column space and fitted values
X = np.column_stack([np.ones(n), Z])
assert np.linalg.matrix_rank(X) == g
b, *_ = np.linalg.lstsq(X, y, rcond=None)
assert np.allclose(X @ b, fitted)
# group standard deviations
sds = np.array([y[group == k].std(ddof=1) for k in range(g)])
assert np.isclose(np.mean(sds ** 2), s2)            # balanced: s^2 is the average of s_k^2

# <<anova>>
P0 = np.full((n, n), 1 / n)                         # projection onto span(1)
ss_mean = y @ P0 @ y                                # n * ybar^2
ss_between = y @ (M - P0) @ y                       # sum_k n_k (ybar_k - ybar)^2
ss_within = y @ (np.eye(n) - M) @ y                 # sum_kj (y_kj - ybar_k)^2
df_between, df_within = g - 1, n - g
F = (ss_between / df_between) / (ss_within / df_within)
p_value = stats.f.sf(F, df_between, df_within)
print(f"{'source':10s} {'df':>3s} {'SS':>9s} {'MS':>8s}")
print(f"{'rates':10s} {df_between:3d} {ss_between:9.4f} {ss_between / df_between:8.4f}")
print(f"{'plots':10s} {df_within:3d} {ss_within:9.4f} {ss_within / df_within:8.4f}")
print(f"{'corrected':10s} {n - 1:3d} {ss_between + ss_within:9.4f}")
print(f"F = {F:.2f} on ({df_between}, {df_within}) df, p = {p_value:.2g}")
# <</anova>>

ybar = y.mean()
assert np.isclose(ss_between, np.sum(n_k * (means - ybar) ** 2))
assert np.isclose(ss_within, np.sum((y - fitted) ** 2))
assert np.isclose(ss_mean + ss_between + ss_within, y @ y)
for A, r in ((P0, 1), (M - P0, g - 1), (np.eye(n) - M, n - g)):
    assert np.allclose(A @ A, A) and np.isclose(np.trace(A), r)
eta2 = ss_between / (ss_between + ss_within)
assert np.isclose(F, eta2 / (1 - eta2) * df_within / df_between)
F_crit = stats.f.ppf(0.95, df_between, df_within)
# a routine one-way ANOVA gives the same statistic
assert np.isclose(stats.f_oneway(*[y[group == k] for k in range(g)]).statistic, F)

# the true noncentrality, known because the data are simulated
gamma_true = m * np.sum((true_means - true_means.mean()) ** 2) / 0.5 ** 2

gen = Generated("ch15", "nitrogen_trial")
for k in range(g):
    gen.text(f"data{k}", ", ".join(f"{v:.2f}" for v in y[group == k]))
    gen.num(f"mean{k}", means[k], 3)
    gen.num(f"sd{k}", sds[k], 3)
gen.num("s", np.sqrt(s2), 4)
gen.num("s2", s2, 4)
gen.num("ybar", ybar, 4)
gen.num("ss_mean", ss_mean, 3)
gen.num("ss_between", ss_between, 4)
gen.num("ss_within", ss_within, 4)
gen.num("ss_total", ss_between + ss_within, 4)
gen.num("ms_between", ss_between / df_between, 4)
gen.num("ms_within", ss_within / df_within, 4)
gen.num("F", F, 2)
gen.text("p", f"{p_value:.1e}")
gen.num("F_crit", F_crit, 3)
gen.num("eta2", eta2, 3)
gen.num("gamma_true", gamma_true, 1)
gen.write()

# ---- figure: the data, group means and the grand mean ---------------------------------
use_book_style()
fig, ax = plt.subplots(figsize=(4.6, 2.7))
jitter = np.tile(np.linspace(-4, 4, m), g)
ax.scatter(rates[group] + jitter, y, s=12, color=COLORS["accent"], alpha=0.85,
           linewidths=0, label="plot yields", zorder=3)
for k in range(g):
    ax.plot([rates[k] - 9, rates[k] + 9], [means[k]] * 2, color=COLORS["second"], lw=1.6,
            label="rate mean" if k == 0 else None)
ax.axhline(ybar, color=COLORS["muted"], lw=0.8, ls="--", label="grand mean")
ax.set_xticks(rates)
ax.set_xlabel("nitrogen rate (kg N/ha)")
ax.set_ylabel("yield (t/ha)")
ax.legend(frameon=False, loc="lower right")
fig.tight_layout()
fig.savefig(figure_path("ch15", "nitrogen_data"))
