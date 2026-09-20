"""Chapter 32, Section 5: Wald inference for a fixed effect when the covariance is estimated.

A cluster-randomized trial: g clusters, half assigned to the treatment, m units in each.
The treatment contrast is estimated by the difference of arm means, and the balanced
design makes its exact distribution a t distribution on g - 2 degrees of freedom. The
simulation compares the coverage of the normal Wald interval, which treats the estimated
covariance as known, with the exact interval.
"""
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy import stats

from regbook import COLORS, Generated, figure_path, use_book_style

# <<trial>>
import numpy as np
from scipy import stats


def simulate(g, m, sigma_a, sigma, effect, B, seed):
    """B cluster-randomized trials; returns the estimate and its estimated variance."""
    rng = np.random.default_rng(seed)
    treat = np.tile([0.0, 1.0], g // 2)                 # cluster k is treated if treat[k]=1
    Y = (effect * treat[:, None] + sigma_a * rng.standard_normal((B, g, 1))
         + sigma * rng.standard_normal((B, g, m)))
    cluster_mean = Y.mean(axis=2)
    diff = (cluster_mean[:, treat == 1].mean(axis=1)
            - cluster_mean[:, treat == 0].mean(axis=1))
    # mean square between clusters within arms, on g - 2 degrees of freedom
    centred = cluster_mean - np.where(treat == 1,
                                      cluster_mean[:, treat == 1].mean(axis=1, keepdims=True),
                                      cluster_mean[:, treat == 0].mean(axis=1, keepdims=True))
    lam1 = m * np.sum(centred**2, axis=1) / (g - 2)     # estimates sigma^2 + m sigma_a^2
    return diff, 4 * lam1 / (g * m)


g, m = 8, 5
# <</trial>>

# <<demo>>
# a shorter run than the 200 000 trials quoted in the text, so the cell finishes quickly
diff, var_hat = simulate(g, m, 0.5, 1.0, 0.0, 20000, seed=5150)
t = diff / np.sqrt(var_hat)
print(f"normal interval covers {np.mean(np.abs(t) <= stats.norm.ppf(0.975)):.3f}")
print(f"t(g-2) interval covers {np.mean(np.abs(t) <= stats.t.ppf(0.975, g - 2)):.3f}")
# <</demo>>

diff, var_hat = simulate(g, m, 0.5, 1.0, 0.0, 200000, seed=5150)
t = diff / np.sqrt(var_hat)
z_cover = np.mean(np.abs(t) <= stats.norm.ppf(0.975))
t_cover = np.mean(np.abs(t) <= stats.t.ppf(0.975, g - 2))
print(f"g = {g}: normal interval covers {z_cover:.3f}, t interval covers {t_cover:.3f}")

assert abs(t_cover - 0.95) < 0.003
assert z_cover < 0.94
# the statistic really is t on g - 2 degrees of freedom
ks = stats.kstest(t, "t", args=(g - 2,))
assert ks.pvalue > 0.01, ks

# ---- coverage against the number of clusters --------------------------------------
sizes = [4, 6, 8, 12, 20, 40]
z_cov, t_cov = [], []
for gg in sizes:
    d, v = simulate(gg, m, 0.5, 1.0, 0.0, 120000, seed=1000 + gg)
    tt = d / np.sqrt(v)
    z_cov.append(np.mean(np.abs(tt) <= stats.norm.ppf(0.975)))
    t_cov.append(np.mean(np.abs(tt) <= stats.t.ppf(0.975, gg - 2)))
assert all(abs(c - 0.95) < 0.004 for c in t_cov)
# the normal interval's coverage is exactly 2 F_{t(g-2)}(1.96) - 1
exact_z = [2 * stats.t.cdf(stats.norm.ppf(0.975), gg - 2) - 1 for gg in sizes]
assert all(abs(c - e) < 0.004 for c, e in zip(z_cov, exact_z))
assert z_cov[0] < 0.82 and z_cov[-1] < 0.946

# ---- one data set: the closed forms are the REML estimates ------------------------
rng = np.random.default_rng(10)
treat = np.tile([0.0, 1.0], g // 2)
a = 0.5 * rng.standard_normal(g)
y = (0.8 * np.repeat(treat, m) + np.repeat(a, m) + rng.standard_normal(g * m))
cl = y.reshape(g, m).mean(axis=1)
arm = np.where(treat == 1, cl[treat == 1].mean(), cl[treat == 0].mean())
ms_cluster = m * np.sum((cl - arm) ** 2) / (g - 2)
ms_error = np.sum((y.reshape(g, m) - cl[:, None]) ** 2) / (g * m - g)
s2a_reml, s2_reml = (ms_cluster - ms_error) / m, ms_error

frame = pd.DataFrame({"y": y, "tr": np.repeat(treat, m), "cl": np.repeat(np.arange(g), m)})
sm_fit = sm.MixedLM.from_formula("y ~ tr", groups="cl", data=frame).fit()
assert abs(sm_fit.scale - s2_reml) < 1e-4, (sm_fit.scale, s2_reml)
assert abs(sm_fit.cov_re.iloc[0, 0] - s2a_reml) < 1e-4

est = cl[treat == 1].mean() - cl[treat == 0].mean()
se = np.sqrt(4 * ms_cluster / (g * m))
assert abs(est - sm_fit.params.iloc[1]) < 1e-8
assert abs(se - sm_fit.bse.iloc[1]) < 1e-4, (se, sm_fit.bse.iloc[1])
half_z, half_t = stats.norm.ppf(0.975) * se, stats.t.ppf(0.975, g - 2) * se
assert half_z < est < half_t      # the normal interval excludes zero, the t interval does not

gen = Generated("ch32", "fixed_effects")
gen.int("g", g)
gen.int("m", m)
gen.num("zcover", z_cover, 3)
gen.num("tcover", t_cover, 3)
gen.num("zcoverfour", z_cov[0], 3)
gen.num("zcoverforty", z_cov[-1], 3)
gen.num("est", est, 3)
gen.num("se", se, 3)
gen.num("halfz", half_z, 3)
gen.num("halft", half_t, 3)
gen.num("ratio", half_t / half_z, 2)
gen.num("s2a", s2a_reml, 3)
gen.num("s2", s2_reml, 3)
gen.num("tcrit", stats.t.ppf(0.975, g - 2), 3)

# how often MS_C < MSE, so that the REML estimate of lambda_1 sits on the boundary
lam1 = 1.0 + m * 0.5**2
p_boundary = stats.f.cdf(1.0 / lam1, g - 2, g * (m - 1))
p_boundary_null = stats.f.cdf(1.0, g - 2, g * (m - 1))

rng2 = np.random.default_rng(4242)
B2 = 200000
Ysim = (0.5 * rng2.standard_normal((B2, g, 1)) + rng2.standard_normal((B2, g, m)))
cl2 = Ysim.mean(axis=2)
arm2 = np.where(treat == 1, cl2[:, treat == 1].mean(axis=1, keepdims=True),
                cl2[:, treat == 0].mean(axis=1, keepdims=True))
msc = m * np.sum((cl2 - arm2) ** 2, axis=1) / (g - 2)
mse = np.sum((Ysim - cl2[:, :, None]) ** 2, axis=(1, 2)) / (g * (m - 1))
assert abs(np.mean(msc < mse) - p_boundary) < 0.004, (np.mean(msc < mse), p_boundary)

gen.num("pboundary", p_boundary, 3)
gen.num("pboundarynull", p_boundary_null, 3)
gen.write()

# ---- figure -----------------------------------------------------------------------
use_book_style()
fig, ax = plt.subplots(figsize=(4.4, 2.6))
ax.plot(sizes, z_cov, "o-", ms=4, color=COLORS["second"], mew=0,
        label="normal Wald interval")
ax.plot(sizes, t_cov, "s--", ms=4, color=COLORS["accent"], mew=0,
        label="$t(g-2)$ interval")
ax.axhline(0.95, color=COLORS["ink"], lw=0.6, ls=":")
ax.set_xlabel("number of clusters $g$")
ax.set_ylabel("coverage of a 95% interval")
ax.set_ylim(0.7, 1.0)
ax.legend(frameon=False, loc="lower right")
fig.tight_layout()
fig.savefig(figure_path("ch32", "wald_coverage"))
