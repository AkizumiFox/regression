"""Chapter 39, Section 1: choosing the response family for a nonnegative response.

Data: statsmodels.datasets.strikes (public domain; a subset of the data used by
Kennan (1985), originally published by the US Bureau of Labor Statistics). The
response is the duration in days of 62 contract strikes; the single covariate is a
measure of unanticipated industrial production. The covariate takes only nine
distinct values, so the data come in replicate groups, and the mean-variance
relationship can be read off directly.

Every likelihood quoted here is maximized over the dispersion parameter as well as
over beta, so that the three AIC values compare maxima of the same kind.
"""
import matplotlib.pyplot as plt
import numpy as np
import statsmodels.api as sm
from scipy import optimize, special, stats

from regbook import COLORS, Generated, figure_path, use_book_style

# <<groups>>
strikes = sm.datasets.strikes.load_pandas().data       # 62 contract strikes, public domain
y = strikes["duration"].to_numpy(float)                # days
x = strikes["iprod"].to_numpy(float)                   # unanticipated industrial production

groups = strikes.groupby("iprod")["duration"].agg(["count", "mean", "var"])
groups = groups[groups["count"] >= 2]                  # nine groups, sizes 2 to 18
w = groups["count"].to_numpy() - 1.0                   # weights: m_g - 1 degrees of freedom
L = np.column_stack([np.ones(len(groups)), np.log(groups["mean"])])
zeta_fit = np.linalg.solve(L.T @ (w[:, None] * L), L.T @ (w * np.log(groups["var"])))

print(groups.round(2))
print(f"variance-function power  zeta = {zeta_fit[1]:.3f}")
# <</groups>>

# Standard error of the slope from the weighted regression.  For normal replicates
# (m-1) s^2 / sigma^2 is chi-squared on m-1 degrees of freedom, so Var(log s^2) is
# exactly trigamma((m-1)/2); the familiar 2/(m-1) is only its large-group limit and
# is far too small for the two groups of size two.
var_log_s2 = special.polygamma(1, w / 2.0)
bread = np.linalg.inv(L.T @ (w[:, None] * L))
cov_zeta = bread @ (L.T @ ((w**2 * var_log_s2)[:, None] * L)) @ bread
se_zeta = np.sqrt(np.diag(cov_zeta))[1]
assert 1.4 < zeta_fit[1] < 2.6, zeta_fit
assert se_zeta > np.sqrt(np.diag(
    bread @ (L.T @ ((w**2 * (2.0 / w))[:, None] * L)) @ bread))[1]   # exact beats 2/(m-1)

# <<families>>
X = np.column_stack([np.ones(len(y)), x])
n, p = X.shape
d = p + 1                                              # slopes plus one dispersion parameter

gam = sm.GLM(y, X, family=sm.families.Gamma(sm.families.links.Log())).fit()
igauss = sm.GLM(y, X, family=sm.families.InverseGaussian(sm.families.links.Log())).fit()
lognorm = sm.OLS(np.log(y), X).fit()                   # a normal linear model for log y


def gamma_loglik(yy, mu, nu):
    return np.sum(stats.gamma.logpdf(yy, nu, scale=mu / nu))


def gamma_shape(yy, mu):
    """Maximum likelihood shape: the root of log nu - digamma(nu) + mean log(y/mu)."""
    c = np.mean(np.log(yy / mu))
    return optimize.brentq(lambda v: np.log(v) - special.digamma(v) + c, 1e-4, 1e4)


def igauss_lambda(yy, mu):
    """Maximum likelihood precision of the inverse Gaussian, in closed form."""
    return 1.0 / np.mean((yy - mu) ** 2 / (mu**2 * yy))


nu_hat = gamma_shape(y, gam.fittedvalues)
lam_hat = igauss_lambda(y, igauss.fittedvalues)
ll_gam = gamma_loglik(y, gam.fittedvalues, nu_hat)
ll_ig = np.sum(stats.invgauss.logpdf(y, igauss.fittedvalues / lam_hat, scale=lam_hat))
sigma_ml = np.sqrt(lognorm.ssr / n)                    # the maximizing sigma, divisor n
ll_ln = np.sum(stats.norm.logpdf(np.log(y), lognorm.fittedvalues, sigma_ml) - np.log(y))

for name, ll in [("gamma, log link", ll_gam), ("inverse Gaussian", ll_ig),
                 ("log-normal (on the y scale)", ll_ln)]:
    print(f"{name:28s} loglik {ll:9.2f}   AIC {-2 * ll + 2 * d:9.2f}")
print(f"gamma:  slope {gam.params[1]:8.3f},  maximum likelihood shape {nu_hat:.3f}")
# <</families>>

aic_gam, aic_ig, aic_ln = -2 * ll_gam + 2 * d, -2 * ll_ig + 2 * d, -2 * ll_ln + 2 * d
aic_ln_logscale = -2 * (ll_ln + np.sum(np.log(y))) + 2 * d       # AIC as a model for log y
jacobian = 2 * np.sum(np.log(y))
assert np.isclose(aic_ln, aic_ln_logscale + jacobian)
assert aic_gam < aic_ln < aic_ig

# each quoted log-likelihood really is the maximum over the dispersion parameter
for v in (0.6, 0.8, 1.25, 1.6):
    assert gamma_loglik(y, gam.fittedvalues, v * nu_hat) < ll_gam
    assert np.sum(stats.invgauss.logpdf(y, igauss.fittedvalues / (v * lam_hat),
                                        scale=v * lam_hat)) < ll_ig
    assert np.sum(stats.norm.logpdf(np.log(y), lognorm.fittedvalues, v * sigma_ml)) \
        < ll_ln + np.sum(np.log(y))

# the Pearson moment estimates that software reports by default, for contrast
phi_pearson = gam.scale
nu_pearson = 1.0 / gam.scale
lam_pearson = 1.0 / igauss.scale
assert lam_pearson > 2 * lam_hat                        # the two disagree badly here

# mean of y under the three fits at the mean covariate value: the retransformation gap
x0 = np.array([1.0, x.mean()])
mean_gamma = np.exp(x0 @ gam.params)
median_lognorm = np.exp(x0 @ lognorm.params)
mean_lognorm = median_lognorm * np.exp(sigma_ml**2 / 2)


# log-score cross-validation: leave one out, predictive log density of the held-out
# point, the dispersion re-maximized on each training set
def loo_logscore(fit_predict):
    total = 0.0
    for i in range(n):
        keep = np.ones(n, bool)
        keep[i] = False
        total += fit_predict(X[keep], y[keep], X[i], y[i])
    return total / n


def gamma_ls(Xtr, ytr, xte, yte):
    f = sm.GLM(ytr, Xtr, family=sm.families.Gamma(sm.families.links.Log())).fit()
    mu, nu = np.exp(xte @ f.params), gamma_shape(ytr, f.fittedvalues)
    return stats.gamma.logpdf(yte, nu, scale=mu / nu)


def lognorm_ls(Xtr, ytr, xte, yte):
    f = sm.OLS(np.log(ytr), Xtr).fit()
    s = np.sqrt(f.ssr / len(ytr))
    return stats.norm.logpdf(np.log(yte), xte @ f.params, s) - np.log(yte)


def igauss_ls(Xtr, ytr, xte, yte):
    f = sm.GLM(ytr, Xtr, family=sm.families.InverseGaussian(sm.families.links.Log())).fit()
    mu, lam = np.exp(xte @ f.params), igauss_lambda(ytr, f.fittedvalues)
    return stats.invgauss.logpdf(yte, mu / lam, scale=lam)


cv_gamma = loo_logscore(gamma_ls)
cv_lognorm = loo_logscore(lognorm_ls)
cv_igauss = loo_logscore(igauss_ls)
print(f"leave-one-out log score: gamma {cv_gamma:.4f}, log-normal {cv_lognorm:.4f}, "
      f"inverse Gaussian {cv_igauss:.4f}")
assert cv_gamma > cv_lognorm > cv_igauss

gen = Generated("ch39", "family_choice")
gen.int("n", n)
gen.int("ngroups", len(groups))
gen.int("mmax", int(groups["count"].max()))
gen.num("ymin", y.min(), 0)
gen.num("ymax", y.max(), 0)
gen.num("zeta", zeta_fit[1], 3)
gen.num("sezeta", se_zeta, 3)
gen.num("dist2", (2.0 - zeta_fit[1]) / se_zeta, 1)
gen.num("dist1", (zeta_fit[1] - 1.0) / se_zeta, 1)
gen.num("aicgamma", aic_gam, 1)
gen.num("aiclognorm", aic_ln, 1)
gen.num("aiclognormraw", aic_ln_logscale, 1)
gen.num("jacobian", jacobian, 1)
gen.num("aicigauss", aic_ig, 1)
gen.num("aicgap", aic_gam - aic_ln_logscale, 1)
gen.num("slope", gam.params[1], 2)
gen.num("se_slope", gam.bse[1], 2)
gen.num("nu", nu_hat, 3)
gen.num("nupearson", nu_pearson, 3)
gen.num("phi", phi_pearson, 3)
gen.num("cvgamma", cv_gamma, 4)
gen.num("cvlognorm", cv_lognorm, 4)
gen.num("cvigauss", cv_igauss, 4)
gen.num("meangamma", mean_gamma, 1)
gen.num("medianlognorm", median_lognorm, 1)
gen.num("meanlognorm", mean_lognorm, 1)
gen.write()

# ---- (a) the mean-variance diagram, (b) the three fitted densities ----------
use_book_style()
fig, axes = plt.subplots(1, 2, figsize=(5.6, 2.5))

ax = axes[0]
ax.scatter(np.log(groups["mean"]), np.log(groups["var"]),
           s=8 + 4 * groups["count"], color=COLORS["accent"], alpha=0.85, linewidths=0)
gx = np.linspace(np.log(groups["mean"]).min() - 0.2, np.log(groups["mean"]).max() + 0.2, 2)
ax.plot(gx, zeta_fit[0] + zeta_fit[1] * gx, color=COLORS["second"],
        label=r"fitted, $\zeta=%.2f$" % zeta_fit[1])
cx, cy = np.log(groups["mean"]).mean(), np.log(groups["var"]).mean()
ax.plot(gx, cy + 1.0 * (gx - cx), color=COLORS["muted"], linestyle=":",
        label=r"$\zeta=1$ (Poisson)")
ax.plot(gx, cy + 2.0 * (gx - cx), color=COLORS["third"], linestyle="--",
        label=r"$\zeta=2$ (gamma)")
ax.set_xlabel(r"$\log$ group mean")
ax.set_ylabel(r"$\log$ group variance")
ax.set_title("(a) nine replicate groups")
ax.legend(loc="upper left", frameon=False)

ax = axes[1]
grid = np.linspace(1, 220, 400)
ax.hist(y, bins=np.arange(0, 230, 15), density=True, color=COLORS["grid"],
        edgecolor="white", linewidth=0.4)
ax.plot(grid, stats.gamma.pdf(grid, nu_hat, scale=mean_gamma / nu_hat),
        color=COLORS["accent"], label="gamma")
ax.plot(grid, stats.norm.pdf(np.log(grid), x0 @ lognorm.params, sigma_ml) / grid,
        color=COLORS["second"], linestyle="--", label="log-normal")
mu_ig = np.exp(x0 @ igauss.params)
ax.plot(grid, stats.invgauss.pdf(grid, mu_ig / lam_hat, scale=lam_hat),
        color=COLORS["third"], linestyle="-.", label="inverse Gaussian")
ax.set_xlabel("strike duration (days)")
ax.set_ylabel("density")
ax.set_title("(b) fitted densities at the mean")
ax.legend(frameon=False)
fig.tight_layout()
fig.savefig(figure_path("ch39", "variance_power"))
