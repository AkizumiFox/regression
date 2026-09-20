"""Chapter 34, Section 2: link functions, on Engel's household budget data.

Three fits with the same regressors: a normal/identity model (least squares), a
gamma/identity model, and a gamma/log model. The random component and the link are
separate choices, and the figure shows what each one changes.

Data: statsmodels.datasets.engel (public domain; 235 Belgian working-class households,
annual income and food expenditure in Belgian francs).
"""
import matplotlib.pyplot as plt
import numpy as np
import statsmodels.api as sm

from regbook import COLORS, Generated, figure_path, use_book_style

# <<fits>>
import numpy as np
import statsmodels.api as sm

engel = sm.datasets.engel.load_pandas().data
y = engel["foodexp"].to_numpy()                  # food expenditure
X = np.column_stack([np.ones(len(y)), np.log(engel["income"].to_numpy())])

gauss = sm.GLM(y, X, family=sm.families.Gaussian()).fit()            # identity link
gam_id = sm.GLM(y, X, family=sm.families.Gamma(sm.families.links.Identity())).fit()
gam_log = sm.GLM(y, X, family=sm.families.Gamma(sm.families.links.Log())).fit()

for name, fit in [("normal/identity", gauss), ("gamma/identity", gam_id), ("gamma/log", gam_log)]:
    print(f"{name:16s} coefficients {fit.params[0]:9.3f} {fit.params[1]:8.3f}"
          f"   dispersion {fit.scale:8.4f}")
# <</fits>>

# The gamma/log fit is a power law: mu = exp(b0) * income^b1, so b1 is an elasticity.
elasticity = gam_log.params[1]
# The two gamma fits differ only in the link; compare them at the largest observed income.
x_grid = np.log(np.array([500.0, 1000.0, engel["income"].max()]))
mu_id = gam_id.params[0] + gam_id.params[1] * x_grid
mu_log = np.exp(gam_log.params[0] + gam_log.params[1] * x_grid)
assert np.all(mu_id > 0) and np.all(mu_log > 0)

# Coefficient of variation: constant for the gamma (sd proportional to the mean),
# and sqrt(phi) estimates it.
cv_gamma = np.sqrt(gam_log.scale)
resid_cv = np.std((y - gam_log.fittedvalues) / gam_log.fittedvalues, ddof=2)
assert abs(cv_gamma - resid_cv) < 0.02

# ---- an identity link can leave the parameter space -------------------------
# <<identity>>
rng = np.random.default_rng(34)
x = np.arange(1.0, 9.0)
counts = rng.poisson(12.0 - 1.2 * x)                             # a mean that declines linearly
pois_id = sm.GLM(counts, np.column_stack([np.ones(8), x]),
                 family=sm.families.Poisson(sm.families.links.Identity())).fit(start_params=[10.0, -1.0])
cross = -pois_id.params[0] / pois_id.params[1]                   # where the fitted mean hits zero
print("counts        ", counts)
print(f"identity-link fitted mean reaches zero at x = {cross:.3f}")
# <</identity>>

assert 8.0 < cross < 10.5                                        # just outside the observed range
assert np.all(pois_id.fittedvalues > 0)

gen = Generated("ch34", "links")
gen.int("n", len(y))
gen.num("gauss_b0", gauss.params[0], 2)
gen.num("gauss_b1", gauss.params[1], 2)
gen.num("gam_id_b0", gam_id.params[0], 2)
gen.num("gam_id_b1", gam_id.params[1], 2)
gen.num("gam_log_b0", gam_log.params[0], 4)
gen.num("elasticity", elasticity, 4)
gen.num("phi_gamma", gam_log.scale, 4)
gen.num("phi_gamma_id", gam_id.scale, 4)
gen.num("cv_gamma_id", float(np.sqrt(gam_id.scale)), 3)
gen.num("cv_gamma", cv_gamma, 3)
gen.num("cross", cross, 2)
gen.num("mu_ratio", float(mu_log[2] / mu_id[2]), 3)
gen.write()

# ---- figure -----------------------------------------------------------------
use_book_style()
inc = engel["income"].to_numpy()
grid = np.linspace(inc.min(), inc.max(), 300)
Xg = np.column_stack([np.ones_like(grid), np.log(grid)])
fit_gauss = Xg @ gauss.params
fit_gid = Xg @ gam_id.params
fit_glog = np.exp(Xg @ gam_log.params)

use_book_style()
fig, axes = plt.subplots(1, 2, figsize=(5.8, 2.5))
ax = axes[0]
ax.scatter(inc, y, s=7, color=COLORS["muted"], alpha=0.6, linewidths=0)
ax.plot(grid, fit_gid, color=COLORS["accent"], label="identity link")
ax.plot(grid, fit_glog, color=COLORS["second"], linestyle="--", label="log link")
ax.set_xlabel("income")
ax.set_ylabel("food expenditure")
ax.set_title("(a) two links, same regressors")
ax.legend(frameon=False, fontsize=7, loc="upper left")
ax = axes[1]
ax.scatter(inc, y, s=7, color=COLORS["muted"], alpha=0.6, linewidths=0)
sd_gauss = np.sqrt(gauss.scale)
ax.plot(grid, fit_gauss, color=COLORS["accent"])
ax.fill_between(grid, fit_gauss - 2 * sd_gauss, fit_gauss + 2 * sd_gauss,
                color=COLORS["accent"], alpha=0.16, linewidth=0, label="normal")
sd_gam = np.sqrt(gam_id.scale) * fit_gid
ax.plot(grid, fit_gid, color=COLORS["third"], linestyle="--")
ax.fill_between(grid, fit_gid - 2 * sd_gam, fit_gid + 2 * sd_gam,
                color=COLORS["third"], alpha=0.16, linewidth=0, label="gamma")
ax.set_xlabel("income")
ax.set_ylabel("food expenditure")
ax.set_title("(b) two variance functions")
ax.legend(frameon=False, fontsize=7, loc="upper left")
fig.tight_layout()
fig.savefig(figure_path("ch34", "links_engel"))
