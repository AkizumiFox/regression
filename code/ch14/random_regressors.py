"""Chapter 14, Section 1: inference conditional on random regressors.

(a) Bivariate normal rows: the least squares slope is unconditionally a scaled t(n-1)
    variable, with variance sigma^2 / (sigma_X^2 (n-3)), not normal.
(b) Skewed random regressors, normal errors independent of them: the t statistic is still
    exactly t(n-p). With errors whose spread depends on the regressor it is not.
(c) The power of the t test depends on the design drawn; the unconditional power is the
    average of the conditional powers.
"""
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats

from regbook import COLORS, Generated, figure_path, use_book_style

gen = Generated("ch14", "random_regressors", prefix="rr")
reps = 100_000

# <<slope>>
rng = np.random.default_rng(1401)
n, beta1, sigma, sigma_x = 12, 0.5, 1.0, 2.0
x = sigma_x * rng.normal(size=(reps, n))            # a new design for every sample
y = 1.0 + beta1 * x + sigma * rng.normal(size=(reps, n))
xc = x - x.mean(axis=1, keepdims=True)
Sxx = np.sum(xc ** 2, axis=1)
b1 = np.sum(xc * y, axis=1) / Sxx                    # least squares slope, one per sample

scale = sigma / (sigma_x * np.sqrt(n - 1))           # b1 - beta1 = scale * t(n-1)
print(f"variance of the slope   {b1.var():.5f}")
print(f"sigma^2/(sigma_x^2 (n-3)) {sigma**2 / (sigma_x**2 * (n - 3)):.5f}")
print(f"sigma^2/E(Sxx)           {sigma**2 / (sigma_x**2 * (n - 1)):.5f}")
print("KS p-value against the scaled t(n-1):",
      round(stats.kstest((b1 - beta1) / scale, "t", args=(n - 1,)).pvalue, 3))
# <</slope>>
var_theory = sigma ** 2 / (sigma_x ** 2 * (n - 3))
var_naive = sigma ** 2 / (sigma_x ** 2 * (n - 1))
assert abs(b1.var() / var_theory - 1) < 0.03
ks_slope = stats.kstest((b1 - beta1) / scale, "t", args=(n - 1,)).pvalue
assert ks_slope > 1e-3
kurt = stats.kurtosis(b1)
assert abs(kurt - 6 / (n - 1 - 4)) < 0.25
gen.int("n", n)
gen.int("reps", reps)
gen.num("var_sim", b1.var(), 5)
gen.num("var_theory", var_theory, 5)
gen.num("var_naive", var_naive, 5)
gen.num("kurt_sim", kurt, 2)
gen.num("kurt_theory", 6 / (n - 5), 2)
gen.num("ks_slope", ks_slope, 2)

# <<skewed>>
rng = np.random.default_rng(1402)
n2, p2 = 10, 3
beta = np.array([1.0, 0.0, 0.5])                     # H0: beta_1 = 0 is true
tq = stats.t.ppf(0.975, n2 - p2)

def t_stats(errors_from):
    """t statistics for beta_1 = 0 over many samples with skewed random regressors."""
    out = np.empty(reps)
    for s in range(reps):
        X = np.column_stack([np.ones(n2), rng.exponential(size=n2), rng.lognormal(size=n2)])
        y = X @ beta + errors_from(X)
        XtX_inv = np.linalg.inv(X.T @ X)
        b = XtX_inv @ X.T @ y
        s2 = np.sum((y - X @ b) ** 2) / (n2 - p2)
        out[s] = b[1] / np.sqrt(s2 * XtX_inv[1, 1])
    return out

t_indep = t_stats(lambda X: rng.normal(size=n2))           # errors independent of X
t_hetero = t_stats(lambda X: X[:, 1] * rng.normal(size=n2))  # spread grows with x_1
for label, t in [("independent errors", t_indep), ("spread depends on X", t_hetero)]:
    print(f"{label:20s} rejection rate {np.mean(np.abs(t) > tq):.4f}")
# <</skewed>>
rej_indep = np.mean(np.abs(t_indep) > tq)
rej_hetero = np.mean(np.abs(t_hetero) > tq)
se_rej = np.sqrt(0.05 * 0.95 / reps)
assert abs(rej_indep - 0.05) < 4 * se_rej
assert rej_hetero > 0.05 + 10 * se_rej
ks_t = stats.kstest(t_indep, "t", args=(n2 - p2,)).pvalue
assert ks_t > 1e-3
gen.int("n2", n2)
gen.int("df2", n2 - p2)
gen.num("rej_indep", rej_indep, 4)
gen.num("rej_hetero", rej_hetero, 4)
gen.num("ks_t", ks_t, 2)

# <<power>>
rng = np.random.default_rng(1403)
n3, b3, df3 = 12, 0.6, 10                            # slope 0.6, sigma = sigma_x = 1
Sxx3 = stats.chi2.rvs(n3 - 1, size=reps, random_state=rng)   # Sxx / sigma_x^2 ~ chi^2(n-1)
gamma = b3 ** 2 * Sxx3                               # noncentrality given the design
fq = stats.f.ppf(0.95, 1, df3)
cond_power = stats.ncf.sf(fq, 1, df3, gamma)         # power of the F = t^2 test given X
print(f"conditional power: quartiles {np.percentile(cond_power, [25, 50, 75]).round(3)}")
print(f"unconditional power {cond_power.mean():.3f}")
print(f"power at the average design {stats.ncf.sf(fq, 1, df3, b3**2 * (n3 - 1)):.3f}")
# <</power>>
q25, q50, q75 = np.percentile(cond_power, [25, 50, 75])
pow_avg_design = stats.ncf.sf(fq, 1, df3, b3 ** 2 * (n3 - 1))
# direct simulation of the unconditional power as a check
xs = rng.normal(size=(20_000, n3))
ys = b3 * xs + rng.normal(size=(20_000, n3))
xcs = xs - xs.mean(axis=1, keepdims=True)
Sx = np.sum(xcs ** 2, axis=1)
bs = np.sum(xcs * ys, axis=1) / Sx
res = ys - ys.mean(axis=1, keepdims=True) - bs[:, None] * xcs
ts = bs / np.sqrt(np.sum(res ** 2, axis=1) / df3 / Sx)
direct = np.mean(ts ** 2 > fq)
assert abs(direct - cond_power.mean()) < 4 * np.sqrt(0.25 / 20_000)
assert pow_avg_design > cond_power.mean()
gen.num("pow_q25", q25, 3)
gen.num("pow_q50", q50, 3)
gen.num("pow_q75", q75, 3)
gen.num("pow_uncond", cond_power.mean(), 3)
gen.num("pow_avgdesign", pow_avg_design, 3)
gen.num("pow_direct", direct, 3)
gen.write()

# ---- figure ------------------------------------------------------------------
use_book_style()
fig, axes = plt.subplots(1, 3, figsize=(6.4, 2.2))
ax = axes[0]
u = (b1 - beta1) / np.sqrt(var_theory)               # standardized to variance one
edges = np.linspace(-5, 5, 61)
dens, _ = np.histogram(u, bins=edges, density=True)
mid = 0.5 * (edges[1:] + edges[:-1])
ok = dens > 0
ax.semilogy(mid[ok], dens[ok], "o", ms=2.2, color=COLORS["accent"], label="simulated")
g = np.linspace(-5, 5, 400)
c = np.sqrt((n - 3) / (n - 1))                       # u = t / sqrt((n-1)/(n-3))
ax.semilogy(g, stats.t.pdf(g / c, n - 1) / c, color=COLORS["ink"], label=r"scaled $t(11)$")
ax.semilogy(g, stats.norm.pdf(g), "--", color=COLORS["second"], label="normal")
ax.set_ylim(1e-4, 1)
ax.set_xlabel("standardized slope")
ax.set_ylabel("density (log scale)")
ax.set_title("(a) slope, random design")
ax.legend(frameon=False, fontsize=6.5, loc="lower center")
ax = axes[1]
ax.hist(t_indep, bins=np.linspace(-6, 6, 61), density=True, color=COLORS["accent"], alpha=0.35)
g = np.linspace(-6, 6, 400)
ax.plot(g, stats.t.pdf(g, n2 - p2), color=COLORS["ink"])
ax.set_xlabel(r"$t$ statistic")
ax.set_title(r"(b) skewed $\mathbf{X}$: $t(7)$ exactly")
ax = axes[2]
ax.hist(cond_power, bins=np.linspace(0, 1, 41), density=True, color=COLORS["third"], alpha=0.45)
ax.axvline(cond_power.mean(), color=COLORS["ink"], lw=1)
ax.set_xlabel("power given the design")
ax.set_title("(c) conditional power")
fig.tight_layout()
fig.savefig(figure_path("ch14", "random_regressors"))
