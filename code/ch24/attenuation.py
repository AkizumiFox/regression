"""Chapter 24, Section 1: attenuation of the least squares slope under classical measurement error.

A simulated study used throughout the chapter: n = 300 patients, long-run systolic blood
pressure x (mean 130, sd 12), age (mean 55, sd 10, correlation 0.5 with x), two clinic readings
w = x + u with u ~ N(0, 9^2), and an outcome y = -4 + 0.08 x + 0 * age + e, e ~ N(0, 1).
Age has no effect once x is known.

(a) The slope on a single reading is attenuated by the reliability ratio 144/225 = 0.64.
(b) With age in the model, the slope on the reading is attenuated by the smaller conditional
    reliability 108/189, and age picks up a spurious coefficient (1 - 108/189) * 0.08 * 0.6.
(c) Repeated samples: the naive interval for the slope almost never covers 0.08, and the t test
    of "no age effect" rejects far more often than 5%.
(d) The large-sample limits (Sigma_xx + Sigma_uu)^{-1} Sigma_xx beta, checked with n = 10^6.
"""
import matplotlib.pyplot as plt

from regbook import COLORS, Generated, figure_path, use_book_style

gen = Generated("ch24", "attenuation", prefix="att")

# <<setup>>
import numpy as np
import statsmodels.api as sm

rng = np.random.default_rng(2401)
n = 300
age = rng.normal(55, 10, n)
x = 130 + 0.6 * (age - 55) + rng.normal(0, np.sqrt(108), n)   # long-run blood pressure, sd 12
w = x[:, None] + rng.normal(0, 9, (n, 2))                      # two clinic readings
y = -4 + 0.08 * x + 0.0 * age + rng.normal(0, 1.0, n)          # age has no effect given x
# <</setup>>

# <<fits>>
oracle = sm.OLS(y, sm.add_constant(x)).fit()                   # uses the true x
naive = sm.OLS(y, sm.add_constant(w[:, 0])).fit()              # uses one reading
lam = 144 / (144 + 81)                                         # reliability of one reading
print(f"slope on x:          {oracle.params[1]:.4f}")
print(f"slope on reading:    {naive.params[1]:.4f}   limit lambda*beta = {lam * 0.08:.4f}")
print("naive 95% interval:", naive.conf_int()[1].round(4))

both = sm.OLS(y, sm.add_constant(np.column_stack([w[:, 0], age]))).fit()
print(f"with age: reading {both.params[1]:.4f}, age {both.params[2]:.4f} (p = {both.pvalues[2]:.4f})")
# <</fits>>

beta, sx2, su2, sz, rho = 0.08, 144.0, 81.0, 10.0, 0.5
sxz = rho * np.sqrt(sx2) * sz                                   # Cov(x, age) = 60
pi = sxz / sz ** 2                                              # regression of x on age: 0.6
sx_z = sx2 - sxz ** 2 / sz ** 2                                 # residual variance of x given age: 108
lam_z = sx_z / (sx_z + su2)
leak = (1 - lam_z) * beta * pi
assert np.isclose(lam, 0.64) and np.isclose(sx_z, 108) and np.isclose(lam_z, 108 / 189)
ci = naive.conf_int()[1]
assert ci[1] < beta                                             # this sample's interval misses 0.08
assert abs(naive.params[1] - lam * beta) < 3 * naive.bse[1]
assert both.pvalues[2] < 0.05

gen.int("n", n)
gen.num("lam", lam, 2)
gen.num("lam_beta", lam * beta, 4)
gen.num("oracle", oracle.params[1], 4)
gen.num("oracle_se", oracle.bse[1], 4)
gen.num("naive", naive.params[1], 4)
gen.num("naive_se", naive.bse[1], 4)
gen.num("naive_lo", ci[0], 4)
gen.num("naive_hi", ci[1], 4)
gen.num("ratio", naive.params[1] / oracle.params[1], 3)
gen.num("lam_z", lam_z, 3)
gen.num("lam_z_beta", lam_z * beta, 4)
gen.num("leak", leak, 4)
gen.num("both_x", both.params[1], 4)
gen.num("both_age", both.params[2], 4)
gen.num("both_age_se", both.bse[2], 4)
gen.num("both_age_p", both.pvalues[2], 4)

# ---- (d) the limits, checked in one very large sample ------------------------------------
big = np.random.default_rng(2402)
N = 1_000_000
ageB = big.normal(55, 10, N)
xB = 130 + 0.6 * (ageB - 55) + big.normal(0, np.sqrt(108), N)
wB = xB + big.normal(0, 9, N)
yB = -4 + 0.08 * xB + big.normal(0, 1.0, N)
bB = np.linalg.lstsq(np.column_stack([np.ones(N), wB, ageB]), yB, rcond=None)[0]
Sxx = np.array([[144.0, 60.0], [60.0, 100.0]])
Suu = np.diag([81.0, 0.0])
limit = np.linalg.solve(Sxx + Suu, Sxx @ np.array([beta, 0.0]))
assert np.allclose(limit, [lam_z * beta, leak])
assert np.allclose(bB[1:], limit, atol=2e-3)
b1B = np.polyfit(wB, yB, 1)[0]
assert abs(b1B - lam * beta) < 1e-3
# residual variance of the naive simple regression: sigma^2 + lam * beta^2 * sigma_u^2
res_var = 1.0 ** 2 + lam * beta ** 2 * su2
assert abs(np.var(yB - np.polyval(np.polyfit(wB, yB, 1), wB)) - res_var) < 0.01
gen.num("res_var", res_var, 4)

# ---- (c) repeated samples ------------------------------------------------------------------
# <<repeat>>
def repeat(reps, seed):
    """Coverage of the naive 95% interval for the slope, and the rate at which age is 'significant'."""
    r = np.random.default_rng(seed)
    cover = reject = 0
    for _ in range(reps):
        a = r.normal(55, 10, n)
        xx = 130 + 0.6 * (a - 55) + r.normal(0, np.sqrt(108), n)
        ww = xx + r.normal(0, 9, n)
        yy = -4 + 0.08 * xx + r.normal(0, 1.0, n)
        lo, hi = sm.OLS(yy, sm.add_constant(ww)).fit().conf_int()[1]
        cover += lo <= 0.08 <= hi
        reject += sm.OLS(yy, sm.add_constant(np.column_stack([ww, a]))).fit().pvalues[2] < 0.05
    return cover / reps, reject / reps
# <</repeat>>

reps = 4000
cover, reject = repeat(reps, 2403)
# <<repeat-small>>
print("coverage of naive interval, rate of rejecting 'no age effect':", repeat(200, 2404))
# <</repeat-small>>
assert cover < 0.02
assert reject > 0.5
gen.int("reps", reps)
gen.num("cover", cover, 3)
gen.num("reject", reject, 3)
gen.write()

# ---- figure ------------------------------------------------------------------------------
use_book_style()
fig, axes = plt.subplots(1, 2, figsize=(5.8, 2.55))
ax = axes[0]
ax.scatter(x, y, s=6, color=COLORS["muted"], alpha=0.55, linewidths=0, label="true $x$")
ax.scatter(w[:, 0], y, s=6, color=COLORS["accent"], alpha=0.55, linewidths=0, label="reading $w$")
xs = np.linspace(95, 165, 2)
ax.plot(xs, oracle.params[0] + oracle.params[1] * xs, color=COLORS["ink"], label="fit on $x$")
ax.plot(xs, naive.params[0] + naive.params[1] * xs, color=COLORS["second"], label="fit on $w$")
ax.set_xlabel("blood pressure (mmHg)")
ax.set_ylabel("outcome $y$")
ax.set_title("(a) the same outcomes, two regressors")
leg = ax.legend(frameon=False, loc="upper left", handlelength=1.2, fontsize=7, markerscale=3)
for h in leg.legend_handles:
    h.set_alpha(1.0)
ax = axes[1]
lams = np.linspace(0.05, 1, 200)
su2s = sx2 * (1 / lams - 1)                                      # error variance giving reliability lam
lz = sx_z / (sx_z + su2s)
ax.plot(lams, lz, color=COLORS["accent"], label=r"coefficient on $w$ / $\beta_x$")
ax.plot(lams, (1 - lz) * pi / 1.0, color=COLORS["second"], label=r"coefficient on age / $\beta_x$")
ax.plot(lams, lams, color=COLORS["muted"], linestyle="--", linewidth=0.9, label=r"$\lambda$ (no age in model)")
ax.axvline(lam, color=COLORS["grid"], linewidth=0.8, zorder=0)
ax.set_xlabel(r"reliability $\lambda$ of one reading")
ax.set_ylabel("limit, in units of $\\beta_x$")
ax.set_title("(b) attenuation and leakage")
ax.set_ylim(0, 1.02)
ax.legend(frameon=False, loc="upper left", fontsize=7)
fig.tight_layout()
fig.savefig(figure_path("ch24", "attenuation"))
