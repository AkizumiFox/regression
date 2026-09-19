"""Chapter 14, Section 4: the sample correlation under bivariate normality and Fisher's z.

(a) Exact law of r: given the x values, the t statistic r sqrt(n-2)/sqrt(1-r^2) is
    t(n-2, delta sqrt(V)) with delta = rho/sqrt(1-rho^2) and V ~ chi^2(n-1); averaging over V
    gives P(r <= c). Checked by simulation.
(b) z = atanh(r): mean about atanh(rho) + rho/(2(n-1)), variance about 1/(n-3), nearly free of rho.
(c) Coverage of the Fisher z interval and of the naive interval r +- 1.96 (1-r^2)/sqrt(n).
(d) The 2009 US state data (statsmodels.datasets.statecrime, public domain): poverty and
    high-school graduation, with the z interval and the exact interval.
"""
import matplotlib.pyplot as plt
import numpy as np
import statsmodels.api as sm
from scipy import optimize, stats

from regbook import COLORS, Generated, figure_path, use_book_style

gen = Generated("ch14", "fisher_z", prefix="fz")


# <<exact>>
def r_cdf(c, rho, n, m=4000):
    """Exact P(r <= c) for n iid bivariate normal pairs with correlation rho."""
    if c < 0:                                                # law of -r is that for -rho
        return 1 - r_cdf(-c, -rho, n, m)
    t = c * np.sqrt(n - 2) / np.sqrt(1 - c ** 2)            # increasing in c
    if rho == 0:
        return stats.t.cdf(t, n - 2)
    V = stats.chi2.ppf((np.arange(m) + 0.5) / m, n - 1)     # equal-probability nodes
    delta = rho / np.sqrt(1 - rho ** 2)
    return np.mean(stats.nct.cdf(t, n - 2, delta * np.sqrt(V)))
# <</exact>>


def sample_r(rho, n, reps, rng):
    x = rng.normal(size=(reps, n))
    y = rho * x + np.sqrt(1 - rho ** 2) * rng.normal(size=(reps, n))
    xc = x - x.mean(axis=1, keepdims=True)
    yc = y - y.mean(axis=1, keepdims=True)
    return np.sum(xc * yc, axis=1) / np.sqrt(np.sum(xc ** 2, axis=1) * np.sum(yc ** 2, axis=1))


# ---- (a) and (b) -------------------------------------------------------------------
rng = np.random.default_rng(1441)
rho, n, reps = 0.8, 15, 200_000
r = sample_r(rho, n, reps, rng)
z = np.arctanh(r)
for q in (0.05, 0.25, 0.5, 0.75, 0.95):
    c = np.quantile(r, q)
    assert abs(r_cdf(c, rho, n) - q) < 0.004, (q, r_cdf(c, rho, n))
zeta = np.arctanh(rho)
assert abs(z.mean() - (zeta + rho / (2 * (n - 1)))) < 0.005
assert abs(z.var() * (n - 3) - 1) < 0.03
gen.num("rho", rho, 1)
gen.int("n", n)
gen.int("reps", reps)
gen.num("r_mean", r.mean(), 4)
gen.num("r_var", r.var(), 5)
gen.num("r_asvar", (1 - rho ** 2) ** 2 / n, 5)
gen.num("r_skew", stats.skew(r), 2)
gen.num("z_mean", z.mean(), 4)
gen.num("zeta", zeta, 4)
gen.num("z_meanapprox", zeta + rho / (2 * (n - 1)), 4)
gen.num("z_var", z.var(), 4)
gen.num("z_var_approx", 1 / (n - 3), 4)
gen.num("z_skew", stats.skew(z), 2)
gen.num("r_median", np.median(r), 3)

# variance stabilization: sd of r and of z across rho, n = 15
rhos = np.linspace(0, 0.95, 20)
sd_r, sd_z = [], []
rng = np.random.default_rng(1442)
for rh in rhos:
    rr = sample_r(rh, n, 40_000, rng)
    sd_r.append(rr.std())
    sd_z.append(np.arctanh(rr).std())
sd_r, sd_z = np.array(sd_r), np.array(sd_z)
assert np.ptp(sd_z) / sd_z.mean() < 0.1 and sd_r[0] / sd_r[-1] > 5
gen.num("sdr0", sd_r[0], 3)
gen.num("sdr95", sd_r[-1], 3)
gen.num("sdz0", sd_z[0], 3)
gen.num("sdz95", sd_z[-1], 3)

# ---- (c) coverage ----------------------------------------------------------------------
# <<coverage>>
rng = np.random.default_rng(1443)
n_c, reps_c, zq = 10, 200_000, stats.norm.ppf(0.975)
for rho_c in (0.0, 0.5, 0.9):
    r_c = sample_r(rho_c, n_c, reps_c, rng)
    z_c = np.arctanh(r_c)
    fisher = np.abs(z_c - np.arctanh(rho_c)) <= zq / np.sqrt(n_c - 3)
    naive = np.abs(r_c - rho_c) <= zq * (1 - r_c ** 2) / np.sqrt(n_c)
    print(f"rho = {rho_c}: Fisher z covers {fisher.mean():.4f}, naive covers {naive.mean():.4f}")
# <</coverage>>
    tag = f"{int(round(rho_c * 10))}"
    gen.num(f"cov_z{tag}", fisher.mean(), 3)
    gen.num(f"cov_naive{tag}", naive.mean(), 3)
    assert abs(fisher.mean() - 0.95) < 0.01
gen.int("nc", n_c)

# ---- (d) the state data -----------------------------------------------------------------
# <<states>>
data = sm.datasets.statecrime.load_pandas().data.drop(index="District of Columbia")
x, y = data["poverty"].to_numpy(), data["hs_grad"].to_numpy()
n_s = len(x)
r_s = np.corrcoef(x, y)[0, 1]
t_s = r_s * np.sqrt(n_s - 2) / np.sqrt(1 - r_s ** 2)
lo_z, hi_z = np.tanh(np.arctanh(r_s) + np.array([-1, 1]) * zq / np.sqrt(n_s - 3))
# exact interval: the rho at which the observed r sits at the 97.5% and 2.5% points
lo_x = optimize.brentq(lambda rh: r_cdf(r_s, rh, n_s) - 0.975, -0.99, 0.99)
hi_x = optimize.brentq(lambda rh: r_cdf(r_s, rh, n_s) - 0.025, -0.99, 0.99)
print(f"r = {r_s:.4f}, t = {t_s:.2f} on {n_s - 2} df")
print(f"Fisher z interval ({lo_z:.3f}, {hi_z:.3f}); exact interval ({lo_x:.3f}, {hi_x:.3f})")
# <</states>>
fit = sm.OLS(y, sm.add_constant(x)).fit()
assert np.isclose(t_s, fit.tvalues[1])
assert lo_x < r_s < hi_x and lo_z < r_s < hi_z
gen.int("ns", n_s)
gen.num("r_s", r_s, 4)
gen.num("t_s", t_s, 2)
gen.num("lo_z", lo_z, 3)
gen.num("hi_z", hi_z, 3)
gen.num("lo_x", lo_x, 3)
gen.num("hi_x", hi_x, 3)
gen.write()

# ---- figure ---------------------------------------------------------------------------------
use_book_style()
fig, axes = plt.subplots(1, 3, figsize=(6.4, 2.2))
ax = axes[0]
ax.hist(r, bins=np.linspace(0, 1, 81), density=True, color=COLORS["accent"], alpha=0.35)
g = np.linspace(0.05, 0.985, 300)
cdf = np.array([r_cdf(c, rho, n, m=800) for c in g])
ax.plot(0.5 * (g[1:] + g[:-1]), np.diff(cdf) / np.diff(g), color=COLORS["ink"])
ax.axvline(rho, color=COLORS["second"], lw=0.8, ls="--")
ax.set_xlabel(r"$r$  ($\rho=0.8$, $n=15$)")
ax.set_ylabel("density")
ax.set_title("(a) exact law of $r$")
ax = axes[1]
ax.hist(z, bins=np.linspace(0, 2.4, 81), density=True, color=COLORS["accent"], alpha=0.35)
g = np.linspace(0, 2.4, 300)
ax.plot(g, stats.norm.pdf(g, zeta + rho / (2 * (n - 1)), 1 / np.sqrt(n - 3)), color=COLORS["ink"])
ax.axvline(zeta, color=COLORS["second"], lw=0.8, ls="--")
ax.set_xlabel(r"$z=\tanh^{-1}r$")
ax.set_title("(b) Fisher's $z$")
ax = axes[2]
ax.plot(rhos, sd_r, "o-", ms=2.5, color=COLORS["second"], label=r"sd of $r$")
ax.plot(rhos, sd_z, "o-", ms=2.5, color=COLORS["accent"], label=r"sd of $z$")
ax.axhline(1 / np.sqrt(n - 3), color=COLORS["muted"], lw=0.8, ls="--")
ax.set_xlabel(r"$\rho$  ($n=15$)")
ax.set_ylim(0, 0.35)
ax.set_title("(c) variance stabilization")
ax.legend(frameon=False, loc="lower left")
fig.tight_layout()
fig.savefig(figure_path("ch14", "fisher_z"))
