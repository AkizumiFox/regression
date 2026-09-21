"""Chapter 40, Section 2: a gamma random-intercept GLMM for Grunfeld's investment panel,
and the marginal moments induced by a random intercept.

Public-domain data shipped with statsmodels (statsmodels.datasets.grunfeld): 11 US firms
observed in each of the 20 years 1935-1954. Chapters 32 and 33 fitted a linear mixed model to
these data; here the mean is multiplicative and the conditional variance grows with its square.
The integrated likelihood is evaluated by adaptive Gauss-Hermite quadrature.
"""
import matplotlib.pyplot as plt
import numpy as np
import statsmodels.api as sm
from numpy.polynomial.hermite import hermgauss
from scipy.optimize import minimize
from scipy.special import gammaln, polygamma

from regbook import COLORS, Generated, figure_path, use_book_style

# <<data>>
d = sm.datasets.grunfeld.load_pandas().data.sort_values(["firm", "year"])
m, n_i = d["firm"].nunique(), 20                      # 11 firms, 20 years each
y = d["invest"].to_numpy().reshape(m, n_i)            # rows are firms
X = np.stack([np.ones((m, n_i)),
              np.log(d["value"].to_numpy()).reshape(m, n_i),
              np.log(d["capital"].to_numpy()).reshape(m, n_i)], axis=2)
# <</data>>


# <<aghq>>
def cluster_mode(beta, tau, nu):
    """Newton's method for the mode and curvature of the integrand, firm by firm."""
    eta0, u = X @ beta, np.zeros(m)
    for _ in range(100):
        s = (y * np.exp(-eta0 - u[:, None])).sum(1)
        step = (-n_i * nu + nu * s - u / tau**2) / (-nu * s - 1 / tau**2)
        u = u - step
        if np.max(np.abs(step)) < 1e-12:
            break
    s = (y * np.exp(-eta0 - u[:, None])).sum(1)
    return u, 1.0 / np.sqrt(nu * s + 1 / tau**2)       # mode, curvature^(-1/2)


def integrated_loglik(beta, tau, nu, K):
    """Adaptive Gauss-Hermite approximation with K nodes; K = 1 is Laplace."""
    node, w = hermgauss(K)
    u_hat, sigma = cluster_mode(beta, tau, nu)
    pts = u_hat[:, None] + np.sqrt(2) * sigma[:, None] * node          # (m, K)
    eta = (X @ beta)[:, :, None] + pts[:, None, :]                     # (m, n_i, K)
    cond = (nu * np.log(nu) - nu * eta + (nu - 1) * np.log(y)[:, :, None]
            - nu * y[:, :, None] * np.exp(-eta) - gammaln(nu)).sum(1)  # (m, K)
    prior = -0.5 * (pts / tau) ** 2 - np.log(tau) - 0.5 * np.log(2 * np.pi)
    terms = np.log(w) + node**2 + cond + prior
    return float(np.sum(np.log(np.sqrt(2) * sigma) + np.logaddexp.reduce(terms, axis=1)))


def fit_glmm(K, start):
    obj = lambda p: -integrated_loglik(p[:3], np.exp(p[3]), np.exp(p[4]), K)
    r = minimize(obj, start, method="Nelder-Mead",
                 options=dict(maxiter=20000, maxfev=20000, xatol=1e-9, fatol=1e-11))
    return r.x, -r.fun
# <</aghq>>


# <<fit>>
flat = np.column_stack([X[:, :, 0].ravel(), X[:, :, 1].ravel(), X[:, :, 2].ravel()])
glm = sm.GLM(y.ravel(), flat, family=sm.families.Gamma(sm.families.links.Log())).fit()
phi_glm = glm.pearson_chi2 / glm.df_resid             # the marginal (independence) fit

par, loglik = fit_glmm(9, np.r_[glm.params, np.log(0.5), np.log(1 / phi_glm)])
beta_hat, tau_hat, phi_hat = par[:3], np.exp(par[3]), np.exp(-par[4])
u_hat, _ = cluster_mode(beta_hat, tau_hat, 1 / phi_hat)

print(f"marginal  GLM  : beta = {glm.params.round(4)},  phi = {phi_glm:.4f}")
print(f"conditional GLMM: beta = {beta_hat.round(4)},  tau = {tau_hat:.4f},  phi = {phi_hat:.4f}")
# <</fit>>

# Laplace (K = 1) against 9-node adaptive quadrature: with n_i = 20 they agree to four places
par1, loglik1 = fit_glmm(1, par)
assert np.max(np.abs(par1 - par)) < 5e-4, np.max(np.abs(par1 - par))
assert abs(loglik1 - loglik) < 1e-2
par15, loglik15 = fit_glmm(15, par)
assert np.max(np.abs(par15 - par)) < 1e-6

# the fit improves on the independence fit, and the random effects have the right scale
assert loglik > -940 and tau_hat > 0.4
assert abs(np.std(u_hat, ddof=0) - tau_hat) < 0.2

# <<shrinkage>>
# the firm's own estimate: maximize firm i's conditional log-likelihood over u alone,
# with beta held at its fitted value. For the gamma with a log link this is in closed form.
a_free = np.log(np.mean(y * np.exp(-X @ beta_hat), axis=1))
shrink = float(np.sum(a_free * u_hat) / np.sum(a_free**2))     # slope through the origin
factor = tau_hat**2 / (tau_hat**2 + phi_hat / n_i)             # the linear-model factor

print(f"empirical shrinkage {shrink:.4f}, linear-model factor {factor:.4f}")
# <</shrinkage>>

# the same firms with only their first three years, at the same fitted parameters
def predict_u(yy, XX, beta, tau, nu):
    u, k = np.zeros(yy.shape[0]), yy.shape[1]
    for _ in range(100):
        s = (yy * np.exp(-(XX @ beta) - u[:, None])).sum(1)
        u = u - (-k * nu + nu * s - u / tau**2) / (-nu * s - 1 / tau**2)
    return u


a_short = np.log(np.mean(y[:, :3] * np.exp(-X[:, :3] @ beta_hat), axis=1))
u_short = predict_u(y[:, :3], X[:, :3], beta_hat, tau_hat, 1 / phi_hat)
shrink_short = float(np.sum(a_short * u_short) / np.sum(a_short**2))
factor_short = tau_hat**2 / (tau_hat**2 + phi_hat / 3)
assert abs(shrink_short - factor_short) < 0.05 and shrink_short < shrink

assert 0 < shrink < 1 and abs(shrink - factor) < 0.03
assert np.all(np.abs(u_hat) <= np.abs(a_free) + 1e-9)          # every effect is pulled in

# ---- marginal moments induced by a random intercept (Poisson, log link) ------
# <<poisson>>
rng = np.random.default_rng(40_002)
tau, eta1, eta2 = 0.6, 0.4, -0.2
u = rng.normal(0, tau, 1_000_000)
y1 = rng.poisson(np.exp(eta1 + u))
y2 = rng.poisson(np.exp(eta2 + u))
mu1, mu2 = np.exp(eta1 + tau**2 / 2), np.exp(eta2 + tau**2 / 2)

print(f"E Y1      : {y1.mean():.4f}  (theory {mu1:.4f})")
print(f"Var Y1    : {y1.var():.4f}  (theory {mu1 + mu1**2 * (np.exp(tau**2) - 1):.4f})")
print(f"Cov(Y1,Y2): {np.cov(y1, y2)[0, 1]:.4f}  (theory {mu1 * mu2 * (np.exp(tau**2) - 1):.4f})")
# <</poisson>>

assert abs(y1.mean() - mu1) < 0.008
assert abs(y1.var() - (mu1 + mu1**2 * (np.exp(tau**2) - 1))) < 0.05
assert abs(np.cov(y1, y2)[0, 1] - mu1 * mu2 * (np.exp(tau**2) - 1)) < 0.02
kappa = 1 / (np.exp(tau**2) - 1)                       # the NB2 shape with the same variance
assert abs(mu1 + mu1**2 / kappa - (mu1 + mu1**2 * (np.exp(tau**2) - 1))) < 1e-12

gen = Generated("ch40", "glmm_grunfeld", prefix="grunfeld")
gen.int("m", m)
gen.int("n_i", n_i)
gen.num("beta_value_glm", glm.params[1], 4)
gen.num("beta_capital_glm", glm.params[2], 4)
gen.num("phi_glm", phi_glm, 4)
gen.num("beta_value", beta_hat[1], 4)
gen.num("beta_capital", beta_hat[2], 4)
gen.num("intercept", beta_hat[0], 4)
gen.num("tau", tau_hat, 4)
gen.num("phi", phi_hat, 4)
gen.num("cv", np.sqrt(phi_hat), 4)
# Var(log Y | u) for a gamma with shape 1/phi is the trigamma function of the shape
var_log_cond = float(polygamma(1, 1 / phi_hat))
gen.num("var_log_cond", var_log_cond, 4)
gen.num("icc_log", tau_hat**2 / (tau_hat**2 + var_log_cond), 3)
gen.num("loglik", loglik, 3)
gen.num("shrink", shrink, 4)
gen.num("factor", factor, 4)
gen.num("shrink_short", shrink_short, 4)
gen.num("intercept_marginal", beta_hat[0] + tau_hat**2 / 2, 4)
gen.num("kappa", float(kappa), 4)
gen.write()

# ---- figure ------------------------------------------------------------------
use_book_style()
fig, axes = plt.subplots(1, 2, figsize=(5.6, 2.5))

ax = axes[0]
lim = 1.15 * max(np.max(np.abs(a_free)), np.max(np.abs(a_short)))
ax.plot([-lim, lim], [-lim, lim], color=COLORS["grid"], linewidth=0.8, zorder=0)
ax.plot([-lim, lim], [-shrink_short * lim, shrink_short * lim], color=COLORS["grid"],
        linewidth=0.8, zorder=0)
ax.scatter(a_free, u_hat, s=14, color=COLORS["accent"], zorder=3, linewidths=0,
           label=f"20 years, slope {shrink:.2f}")
ax.scatter(a_short, u_short, s=14, marker="^", color=COLORS["third"], zorder=2, linewidths=0,
           label=f"first 3 years, slope {shrink_short:.2f}")
ax.set_xlabel(r"the firm's own estimate $\hat a_i$")
ax.set_ylabel(r"predicted random effect $\hat u_i$")
ax.set_title("(a) shrinkage of the firm effects")
ax.legend(frameon=False, loc="upper left", fontsize=7)

ax = axes[1]
ns = np.arange(1, 41)
for t, colour in [(tau_hat, COLORS["accent"]), (0.25, COLORS["second"]), (1.0, COLORS["third"])]:
    ax.plot(ns, t**2 / (t**2 + phi_hat / ns), color=colour, label=rf"$\tau = {t:.2f}$")
ax.axvline(n_i, color=COLORS["grid"], linewidth=0.8, zorder=0)
ax.set_xlabel(r"observations per cluster $n_i$")
ax.set_ylabel("shrinkage factor")
ax.set_ylim(0, 1.05)
ax.set_title(r"(b) $\tau^2/(\tau^2+\phi/n_i)$")
ax.legend(frameon=False, loc="lower right", fontsize=7)

fig.tight_layout()
fig.savefig(figure_path("ch40", "shrinkage"))
