"""Chapter 39, Section 2: diagnosing a badly fitting generalized linear model.

Data: statsmodels.datasets.anes96 (public domain; the 1996 American National
Election Study, as distributed with statsmodels). The response is the number of
days in the past week on which the respondent watched the television news, taken
as a binomial count out of seven; the covariates are age, education and income.
The model fits badly, and this script works out which part of it is wrong.
"""
import matplotlib.pyplot as plt
import numpy as np
import statsmodels.api as sm
from scipy import optimize, special, stats

from regbook import COLORS, Generated, figure_path, use_book_style

# <<fit>>
anes = sm.datasets.anes96.load_pandas().data
m = 7.0                                                   # days in a week
y = anes["TVnews"].to_numpy(float)                        # successes out of seven
X = np.column_stack([np.ones(len(y)), anes["age"], anes["educ"], anes["income"]])
n, p = X.shape

fit = sm.GLM(np.column_stack([y, m - y]), X, family=sm.families.Binomial()).fit()
mu = fit.fittedvalues * m                                 # fitted expected days
eta = fit.predict(which="linear")
pearson = (y - mu) / np.sqrt(mu * (1 - mu / m))           # r^P, with phi = 1

print(f"deviance {fit.deviance:8.1f} on {n - p} degrees of freedom")
print(f"Pearson  {np.sum(pearson**2):8.1f};  dispersion estimate "
      f"{np.sum(pearson**2) / (n - p):.3f}")
# <</fit>>

phi_hat = np.sum(pearson**2) / (n - p)
assert fit.deviance / (n - p) > 2 and phi_hat > 2

# weighted hat matrix and a Cook-type influence measure
W = m * fit.fittedvalues * (1 - fit.fittedvalues)         # working weights
root = np.sqrt(W)[:, None] * X
h = np.einsum("ij,jk,ik->i", root, np.linalg.inv(root.T @ root), root)
assert np.isclose(h.sum(), p)
cook = pearson**2 * h / (p * (1 - h) ** 2)

# <<envelope>>
B = 19                                                    # simulated data sets


def halfnormal_quantiles(k):
    """The horizontal axis: expected half-normal order statistics."""
    i = np.arange(1, k + 1)
    return stats.norm.ppf((i + k - 0.125) / (2 * k + 0.5))


def envelope(draw, refit_abs_resid, seed):
    """Smallest and largest of B ordered |residual| vectors from the fitted model."""
    rng = np.random.default_rng(seed)
    out = np.array([np.sort(refit_abs_resid(draw(rng))) for _ in range(B)])
    return out.min(axis=0), out.max(axis=0)


def binomial_abs_resid(ysim):
    f = sm.GLM(np.column_stack([ysim, m - ysim]), X, family=sm.families.Binomial()).fit()
    mu_s = f.fittedvalues * m
    return np.abs((ysim - mu_s) / np.sqrt(mu_s * (1 - mu_s / m)))


q = halfnormal_quantiles(n)
lo, hi = envelope(lambda rng: rng.binomial(int(m), fit.fittedvalues).astype(float),
                  binomial_abs_resid, seed=3901)
ordered = np.sort(np.abs(pearson))
outside, above = np.mean((ordered > hi) | (ordered < lo)), np.mean(ordered > hi)
print(f"binomial model: {100 * outside:.1f}% of the ordered residuals leave the envelope, "
      f"{100 * above:.1f}% above it")
# <</envelope>>

assert outside > 0.4 and above > 0.4

# <<checks>>
# (i) is the link right?  Add the squared linear predictor and test it.
fit_link = sm.GLM(np.column_stack([y, m - y]), np.column_stack([X, eta**2]),
                  family=sm.families.Binomial()).fit()
lr_link = fit.deviance - fit_link.deviance

# (ii) is the age term right?  A partial residual, on the link scale.
working = (y - mu) / (mu * (1 - mu / m))                  # (y - mu) dEta/dMu
partial_age = fit.params[1] * anes["age"].to_numpy() + working

# (iii) is the variance function right?  Regress log (r^P)^2 on the linear predictor.
lin = sm.OLS(np.log(pearson**2 + 1e-12), sm.add_constant(eta)).fit()

print(f"link test: deviance drop {lr_link:.2f} on 1 d.f., p = {stats.chi2.sf(lr_link, 1):.3f}")
print(f"slope of log (r^P)^2 on the linear predictor: {lin.params[1]:.3f} "
      f"(t = {lin.tvalues[1]:.2f})")
# <</checks>>

# the three findings the example turns on, asserted rather than merely printed
assert stats.chi2.sf(lr_link, 1) > 0.2                    # the link test is non-significant
assert h.max() < 0.05                                     # every leverage is tiny
assert lin.params[1] > 0.2 and lin.tvalues[1] > 3         # the dispersion trends upward

# local means of |r^P| in ten equal bins of the linear predictor: the picture that
# panel (b) draws, and the one the regression slope above summarizes
edges = np.quantile(eta, np.linspace(0, 1, 11))
which = np.clip(np.searchsorted(edges, eta, side="right") - 1, 0, 9)
bin_eta = np.array([eta[which == k].mean() for k in range(10)])
bin_abs = np.array([np.abs(pearson)[which == k].mean() for k in range(10)])
assert np.all(bin_abs > np.sqrt(2 / np.pi)) and bin_abs[-1] > bin_abs[0]

# a quadratic in age: the mean model has nothing left to give
Xq = np.column_stack([X, (anes["age"].to_numpy() - anes["age"].mean()) ** 2])
fit_q = sm.GLM(np.column_stack([y, m - y]), Xq, family=sm.families.Binomial()).fit()
lr_age2 = fit.deviance - fit_q.deviance
assert lr_age2 < 4

# <<betabinom>>
scale = np.array([1.0, X[:, 1].std(), X[:, 2].std(), X[:, 3].std()])   # keep BFGS conditioned


def betabinom_fit(ysim, start):
    """Maximize the beta-binomial likelihood: logit mean X beta, precision s."""
    def negll(theta):
        pi = special.expit((X / scale) @ theta[:p])
        s = np.exp(theta[p])
        v = -np.sum(special.betaln(ysim + s * pi, m - ysim + s * (1 - pi))
                    - special.betaln(s * pi, s * (1 - pi)))
        return v if np.isfinite(v) else 1e12

    theta = np.append(start[:p] * scale, start[p])
    for _ in range(3):                                    # restart until the gradient is flat
        opt = optimize.minimize(negll, theta, method="BFGS")
        theta = opt.x
    return np.append(theta[:p] / scale, theta[p]), opt


theta_bb, opt = betabinom_fit(y, np.append(fit.params, 0.1))
beta_bb, s_bb = theta_bb[:p], np.exp(theta_bb[p])
rho = 1.0 / (1.0 + s_bb)                                  # correlation within the week
print(f"beta-binomial: s = {s_bb:.3f}, within-respondent correlation {rho:.3f}")
for j, name in enumerate(["intercept", "age", "educ", "income"]):
    print(f"  {name:10s} binomial {fit.params[j]:8.4f}   beta-binomial {beta_bb[j]:8.4f}")
# <</betabinom>>

assert np.max(np.abs(opt.jac)) < 1e-3 and 0.2 < rho < 0.8


def betabinom_negll(theta):
    pi = special.expit(X @ theta[:p])
    s = np.exp(theta[p])
    return -np.sum(special.betaln(y + s * pi, m - y + s * (1 - pi))
                   - special.betaln(s * pi, s * (1 - pi)))


def numerical_hessian(f, theta, step=1e-4):
    k = len(theta)
    H = np.empty((k, k))
    for a in range(k):
        for b in range(k):
            ea, eb = np.zeros(k), np.zeros(k)
            ea[a], eb[b] = step, step
            H[a, b] = (f(theta + ea + eb) - f(theta + ea - eb)
                       - f(theta - ea + eb) + f(theta - ea - eb)) / (4 * step**2)
    return H


se_bb = np.sqrt(np.diag(np.linalg.inv(numerical_hessian(betabinom_negll, theta_bb))))[:p]
assert np.all(se_bb > fit.bse) and np.all(se_bb < 10 * fit.bse)
se_quasi = fit.bse * np.sqrt(phi_hat)

# the envelope again, now under the beta-binomial fit
pi_bb = special.expit(X @ beta_bb)
var_bb = m * pi_bb * (1 - pi_bb) * (m + s_bb) / (1 + s_bb)
pearson_bb = (y - m * pi_bb) / np.sqrt(var_bb)


def bb_draw(rng):
    return rng.binomial(int(m), rng.beta(s_bb * pi_bb, s_bb * (1 - pi_bb))).astype(float)


def bb_abs_resid(ysim):
    th, _ = betabinom_fit(ysim, theta_bb)
    pi = special.expit(X @ th[:p])
    s = np.exp(th[p])
    return np.abs((ysim - m * pi) / np.sqrt(m * pi * (1 - pi) * (m + s) / (1 + s)))


lo_bb, hi_bb = envelope(bb_draw, bb_abs_resid, seed=3902)
ordered_bb = np.sort(np.abs(pearson_bb))
outside_bb = np.mean((ordered_bb > hi_bb) | (ordered_bb < lo_bb))
above_bb = np.mean(ordered_bb > hi_bb)
print(f"beta-binomial: {100 * outside_bb:.1f}% leave the envelope, {100 * above_bb:.1f}% above it")
assert above_bb < 0.02 and outside_bb < 0.3

gen = Generated("ch39", "diagnostics")
gen.int("n", n)
gen.int("p", p)
gen.int("df", n - p)
gen.int("B", B)
gen.int("nzero", int(np.sum(y == 0)))
gen.int("nseven", int(np.sum(y == m)))
gen.num("deviance", fit.deviance, 1)
gen.num("chisq", float(np.sum(pearson**2)), 1)
gen.num("phi", phi_hat, 3)
gen.num("sqrtphi", float(np.sqrt(phi_hat)), 3)
gen.num("outside", 100 * outside, 1)
gen.num("above", 100 * above, 1)
gen.num("outsidebb", 100 * outside_bb, 1)
gen.num("abovebb", 100 * above_bb, 1)
gen.num("lrlink", lr_link, 2)
gen.num("plink", float(stats.chi2.sf(lr_link, 1)), 3)
gen.num("varslope", lin.params[1], 3)
gen.num("vart", lin.tvalues[1], 2)
gen.num("lrage2", lr_age2, 2)
gen.num("hmax", h.max(), 4)
gen.num("hmean", h.mean(), 4)
gen.num("cookmax", cook.max(), 4)
gen.num("s", s_bb, 3)
gen.num("rho", rho, 3)
gen.num("varratio", (m + s_bb) / (1 + s_bb), 2)
gen.num("binlo", float(bin_abs[0]), 2)
gen.num("binhi", float(bin_abs.max()), 2)
gen.num("beta_age", fit.params[1], 4)
gen.num("se_age", fit.bse[1], 4)
gen.num("se_age_quasi", se_quasi[1], 4)
gen.num("beta_age_bb", beta_bb[1], 4)
gen.num("se_age_bb", se_bb[1], 4)
gen.num("beta_educ", fit.params[2], 4)
gen.num("se_educ", fit.bse[2], 4)
gen.num("beta_educ_bb", beta_bb[2], 4)
gen.num("se_educ_bb", se_bb[2], 4)
gen.num("se_educ_quasi", se_quasi[2], 4)
gen.write()

# ---- Figure 1: half-normal plots with simulated envelopes -------------------
use_book_style()
fig, axes = plt.subplots(1, 2, figsize=(5.6, 2.6), sharex=True, sharey=True)
for ax, (r, bnd, title) in zip(axes, [(pearson, (lo, hi), "(a) binomial"),
                                      (pearson_bb, (lo_bb, hi_bb), "(b) beta-binomial")]):
    ax.fill_between(q, bnd[0], bnd[1], color=COLORS["grid"], alpha=0.8, linewidth=0)
    ax.scatter(q, np.sort(np.abs(r)), s=3, color=COLORS["accent"], linewidths=0)
    ax.set_xlabel("half-normal quantile")
    ax.set_title(title)
axes[0].set_ylabel(r"$|r^{P}|$, ordered")
fig.tight_layout()
fig.savefig(figure_path("ch39", "halfnormal_envelope"))

# ---- Figure 2: partial residuals, variance check, influence ----------------
fig, axes = plt.subplots(1, 3, figsize=(6.3, 2.3))
age = anes["age"].to_numpy(float)
ax = axes[0]
ax.scatter(age, partial_age, s=3, color=COLORS["accent"], alpha=0.45, linewidths=0)
smooth = sm.nonparametric.lowess(partial_age, age, frac=0.5, return_sorted=True)
ax.plot(smooth[:, 0], smooth[:, 1], color=COLORS["second"])
grid = np.linspace(age.min(), age.max(), 2)
ax.plot(grid, fit.params[1] * grid, color=COLORS["muted"], linestyle=":")
ax.set_xlabel("age (years)")
ax.set_ylabel("partial residual")
ax.set_title("(a) the age term")

ax = axes[1]
ax.scatter(eta, np.abs(pearson), s=3, color=COLORS["accent"], alpha=0.45, linewidths=0)
ax.plot(bin_eta, bin_abs, "o-", color=COLORS["second"], markersize=3, linewidth=0.9)
ax.axhline(np.sqrt(2 / np.pi), color=COLORS["muted"], linestyle=":")
ax.set_xlabel(r"linear predictor")
ax.set_ylabel(r"$|r^{P}|$")
ax.set_title("(b) the variance function")

ax = axes[2]
ax.vlines(np.arange(n), 0, cook, color=COLORS["accent"], linewidth=0.5)
for i in np.argsort(cook)[-2:]:
    ax.annotate(str(i), (i, cook[i]), fontsize=6, ha="center", va="bottom")
ax.set_xlabel("observation")
ax.set_ylabel("Cook-type distance")
ax.set_title("(c) influence")
fig.tight_layout()
fig.savefig(figure_path("ch39", "glm_diagnostics"))
