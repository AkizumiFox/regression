"""Chapter 16, Sections 1-2: the additive two-way model with one observation per cell, and
Tukey's one-degree-of-freedom test for nonadditivity.

Synthetic data: running times (seconds) of 6 solvers on 8 benchmark problems, one run
each. Times are generated multiplicatively (log time = solver + problem + noise), so the
raw scale is nonadditive and the log scale is additive. The script
(1) fits the additive model on both scales and computes Tukey's statistic;
(2) checks the suggested power 1 - theta_hat * ybar;
(3) checks by simulation that Tukey's F has the central F(1, (a-1)(b-1)-1) distribution
    under an additive normal model (its size is alpha);
(4) draws the diagnostic plot: residuals against a_i b_j / ybar.
"""
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats

from regbook import COLORS, Generated, figure_path, use_book_style

# ---- data --------------------------------------------------------------------------
gen_rng = np.random.default_rng(16_01)
a, b = 6, 8
solver_log = np.array([0.0, 0.25, -0.35, 0.55, 0.1, -0.2])            # log speed factors
problem_log = np.log([2.0, 4.5, 7.0, 12.0, 25.0, 60.0, 110.0, 240.0])  # typical seconds
noise = gen_rng.normal(scale=0.12, size=(a, b))
times = np.round(np.exp(solver_log[:, None] + problem_log[None, :] + noise), 1)


# <<tukey>>
def additive_fit(Y):
    """Row effects, column effects and residuals of the additive fit to an a x b table."""
    g = Y.mean()
    r = Y.mean(axis=1) - g                    # ybar_i. - ybar..
    c = Y.mean(axis=0) - g                    # ybar_.j - ybar..
    resid = Y - g - r[:, None] - c[None, :]
    return g, r, c, resid


def tukey_test(Y):
    """Tukey's one-degree-of-freedom test for nonadditivity in an a x b table."""
    a, b = Y.shape
    g, r, c, resid = additive_fit(Y)
    u = np.outer(r, c)                        # lies in the interaction space
    ss_n = (u * Y).sum() ** 2 / (u * u).sum()
    sse = (resid ** 2).sum()
    df = (a - 1) * (b - 1) - 1
    F = ss_n / ((sse - ss_n) / df)
    theta = (u * Y).sum() / (u * u).sum()     # coefficient of a_i b_j
    return F, stats.f.sf(F, 1, df), theta, g, sse, ss_n


for label, Y in [("seconds", times), ("log seconds", np.log(times))]:
    F, p, theta, g, sse, ss_n = tukey_test(Y)
    print(f"{label:12s} F = {F:9.2f}  p = {p:.2g}  theta = {theta:.4f}"
          f"  suggested power 1 - theta*ybar = {1 - theta * g:.2f}")
# <</tukey>>

F_raw, p_raw, theta_raw, g_raw, sse_raw, ssn_raw = tukey_test(times)
F_log, p_log, theta_log, g_log, sse_log, ssn_log = tukey_test(np.log(times))
power_raw = 1 - theta_raw * g_raw
df_n = (a - 1) * (b - 1) - 1
assert p_raw < 1e-6 and p_log > 0.05
assert abs(power_raw) < 0.3                   # the log (power 0) is suggested
# fraction of the raw residual SS that the single Tukey direction takes
share_raw = ssn_raw / sse_raw

# the ordinary additive ANOVA on each scale
def anova(Y):
    a, b = Y.shape
    g, r, c, resid = additive_fit(Y)
    return b * (r ** 2).sum(), a * (c ** 2).sum(), (resid ** 2).sum()

ssA_raw, ssB_raw, sse_raw2 = anova(times)
ssA_log, ssB_log, sse_log2 = anova(np.log(times))
assert np.isclose(sse_raw2, sse_raw)
FA_raw = ssA_raw / (a - 1) / (sse_raw / ((a - 1) * (b - 1)))
FA_log = ssA_log / (a - 1) / (sse_log / ((a - 1) * (b - 1)))
pA_raw = stats.f.sf(FA_raw, a - 1, (a - 1) * (b - 1))
pA_log = stats.f.sf(FA_log, a - 1, (a - 1) * (b - 1))
s_log = np.sqrt(sse_log / ((a - 1) * (b - 1)))
# check the additive fit against statsmodels on the log scale
import pandas as pd
import statsmodels.formula.api as smf
from statsmodels.stats.anova import anova_lm

df = pd.DataFrame({"t": np.log(times).ravel(), "solver": np.repeat(np.arange(a), b),
                   "problem": np.tile(np.arange(b), a)})
tab = anova_lm(smf.ols("t ~ C(solver) + C(problem)", df).fit())
assert np.isclose(tab.loc["C(solver)", "sum_sq"], ssA_log)
assert np.isclose(tab.loc["Residual", "sum_sq"], sse_log)

# solver differences on the log scale: ratios of running times
g, r, c, _ = additive_fit(np.log(times))
fastest, slowest = int(np.argmin(r)), int(np.argmax(r))
ratio = np.exp(r[slowest] - r[fastest])
se_diff = s_log * np.sqrt(2 / b)
q = stats.studentized_range.ppf(0.95, a, (a - 1) * (b - 1))
tukey_half = q * s_log / np.sqrt(b)
lo, hi = np.exp(r[slowest] - r[fastest] - tukey_half), np.exp(r[slowest] - r[fastest] + tukey_half)

# ---- (3) the null distribution of Tukey's F ---------------------------------------
# <<null>>
rng = np.random.default_rng(20260919)
reps = 20_000
Fs = np.empty(reps)
for t in range(reps):
    Y = 3.0 + np.arange(a)[:, None] * 0.5 + np.arange(b)[None, :] + rng.normal(size=(a, b))
    Fs[t] = tukey_test(Y)[0]
size = np.mean(Fs > stats.f.ppf(0.95, 1, (a - 1) * (b - 1) - 1))
ks = stats.kstest(Fs, stats.f(1, (a - 1) * (b - 1) - 1).cdf)
print(f"simulated size of the 5% test: {size:.4f}   Kolmogorov-Smirnov p-value: {ks.pvalue:.2f}")
# <</null>>
assert abs(size - 0.05) < 4 * np.sqrt(0.05 * 0.95 / reps)
assert ks.pvalue > 0.01

gen = Generated("ch16", "tukey_runtime")
for i in range(a):
    for j in range(b):
        gen.num(f"t{i + 1}{j + 1}", times[i, j], 1)
gen.num("F_raw", F_raw, 1)
gen.text("p_raw", f"{p_raw:.0e}".replace("e-0", "e-").replace("e-", r"\times 10^{-") + "}")
gen.num("F_log", F_log, 2)
gen.num("p_log", p_log, 2)
gen.num("theta_raw", theta_raw, 5)
gen.num("ybar_raw", g_raw, 2)
gen.num("power_raw", power_raw, 2)
gen.num("theta_log", theta_log, 3)
gen.num("share_raw", share_raw, 3)
gen.num("ssn_raw", ssn_raw, 0)
gen.num("sse_raw", sse_raw, 0)
gen.int("df_n", df_n)
gen.int("df_e", (a - 1) * (b - 1))
gen.num("FA_raw", FA_raw, 2)
gen.num("pA_raw", pA_raw, 3)
gen.num("FA_log", FA_log, 1)
gen.num("ssA_log", ssA_log, 3)
gen.num("ssB_log", ssB_log, 2)
gen.num("sse_log", sse_log, 3)
gen.num("s_log", s_log, 3)
gen.int("fastest", fastest + 1)
gen.int("slowest", slowest + 1)
gen.num("ratio", ratio, 2)
gen.num("ratio_lo", lo, 2)
gen.num("ratio_hi", hi, 2)
gen.num("q", q, 3)
gen.num("size", size, 4)
gen.num("ks_p", ks.pvalue, 2)
gen.int("reps", reps)
gen.write()

# ---- (4) the diagnostic plot --------------------------------------------------------
use_book_style()
fig, axes = plt.subplots(1, 2, figsize=(5.6, 2.4))
for ax, Y, title in [(axes[0], times, "(a) seconds"), (axes[1], np.log(times), "(b) log seconds")]:
    g, r, c, resid = additive_fit(Y)
    comp = np.outer(r, c) / g
    ax.scatter(comp.ravel(), resid.ravel(), s=10, color=COLORS["accent"], alpha=0.85, linewidths=0)
    theta = tukey_test(Y)[2]
    xs = np.linspace(comp.min(), comp.max(), 2)
    ax.plot(xs, theta * g * xs, color=COLORS["second"])
    ax.axhline(0, color=COLORS["grid"], linewidth=0.6, zorder=0)
    ax.set_xlabel(r"comparison value $r_i c_j/\bar{y}_{\cdot\cdot}$")
    ax.set_ylabel("residual of additive fit")
    ax.set_title(title + r": slope $\hat\theta\,\bar{y}_{\cdot\cdot}$ = %.2f" % (theta * g))
fig.tight_layout()
fig.savefig(figure_path("ch16", "tukey_diagnostic"))
