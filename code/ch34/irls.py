"""Chapter 34, Section 4: iteratively reweighted least squares, written from scratch.

Fits a Poisson/log, a binomial/logit, a binomial/probit and a gamma/log model by Fisher
scoring and checks every one against statsmodels. Then measures how fast the iteration
converges for a canonical and for a non-canonical link, and shows one failure mode.
"""
import matplotlib.pyplot as plt
import numpy as np
import statsmodels.api as sm
from scipy import stats

from regbook import COLORS, Generated, figure_path, use_book_style

# <<families>>
import numpy as np
from scipy import stats

# A family is the variance function V, the link g, its inverse h, the derivative
# h' = d(mu)/d(eta), and a rule for the starting fitted means.
def poisson_log():
    return dict(name="Poisson/log", V=lambda mu: mu, g=np.log, h=np.exp,
                dmu=lambda eta: np.exp(eta), start=lambda y, w: y + 0.1)

def binomial_logit():
    return dict(name="binomial/logit", V=lambda mu: mu * (1 - mu),
                g=lambda mu: np.log(mu / (1 - mu)), h=lambda eta: 1 / (1 + np.exp(-eta)),
                dmu=lambda eta: 1 / (2 + np.exp(eta) + np.exp(-eta)),
                start=lambda y, w: (w * y + 0.5) / (w + 1.0))

def binomial_probit():
    return dict(name="binomial/probit", V=lambda mu: mu * (1 - mu),
                g=stats.norm.ppf, h=stats.norm.cdf, dmu=stats.norm.pdf,
                start=lambda y, w: (w * y + 0.5) / (w + 1.0))

def gamma_log():
    return dict(name="gamma/log", V=lambda mu: mu ** 2, g=np.log, h=np.exp,
                dmu=lambda eta: np.exp(eta), start=lambda y, w: y)

def irls(X, y, family, w=None, tol=1e-12, max_iter=60, trace=False):
    """Fisher scoring for a generalized linear model, as weighted least squares.

    At each step the working response is z = eta + (y - mu)/h'(eta) and the working
    weight is W = w h'(eta)^2 / V(mu); beta is the weighted least squares coefficient
    of z on X. The dispersion phi never enters, because it cancels.
    """
    w = np.ones(len(y)) if w is None else np.asarray(w, float)
    mu = family["start"](y, w)
    eta = family["g"](mu)
    beta, path = np.zeros(X.shape[1]), []
    for _ in range(max_iter):
        d = family["dmu"](eta)                             # h'(eta)
        W = w * d ** 2 / family["V"](mu)                   # working weights
        z = eta + (y - mu) / d                             # working response
        XtW = X.T * W
        beta_new = np.linalg.solve(XtW @ X, XtW @ z)       # one weighted least squares fit
        path.append(beta_new)
        eta = X @ beta_new
        mu = family["h"](eta)
        if np.max(np.abs(beta_new - beta)) < tol * (1 + np.max(np.abs(beta_new))):
            beta = beta_new
            break
        beta = beta_new
    W = w * family["dmu"](eta) ** 2 / family["V"](mu)      # weights at the solution
    cov_unscaled = np.linalg.inv((X.T * W) @ X)            # (X^T W X)^{-1}
    return (beta, cov_unscaled, np.array(path)) if trace else (beta, cov_unscaled)
# <</families>>

# ---- check against statsmodels ----------------------------------------------
# <<check>>
import statsmodels.api as sm

rng = np.random.default_rng(3401)
n = 300
Xs = np.column_stack([np.ones(n), rng.normal(size=n), rng.binomial(1, 0.4, n)])
beta_true = np.array([0.7, 0.5, -0.4])
counts = rng.poisson(np.exp(Xs @ beta_true))

beta_hat, cov = irls(Xs, counts.astype(float), poisson_log())
sm_fit = sm.GLM(counts, Xs, family=sm.families.Poisson()).fit()
print("from scratch:", np.round(beta_hat, 6))
print("statsmodels :", np.round(sm_fit.params, 6))
print("standard errors:", np.round(np.sqrt(np.diag(cov)), 6), np.round(sm_fit.bse, 6))
# <</check>>

assert np.allclose(beta_hat, sm_fit.params, atol=1e-8)
assert np.allclose(np.sqrt(np.diag(cov)), sm_fit.bse, atol=1e-8)

# binomial, grouped: m trials in each of n groups, response = observed proportion
m = rng.integers(5, 30, size=n).astype(float)
p_true = 1 / (1 + np.exp(-(Xs @ np.array([0.2, 0.9, -0.6]))))
succ = rng.binomial(m.astype(int), p_true).astype(float)
prop = succ / m

b_logit, cov_logit = irls(Xs, prop, binomial_logit(), w=m)
sm_logit = sm.GLM(np.column_stack([succ, m - succ]), Xs, family=sm.families.Binomial()).fit()
assert np.allclose(b_logit, sm_logit.params, atol=1e-8)
assert np.allclose(np.sqrt(np.diag(cov_logit)), sm_logit.bse, atol=1e-8)

b_probit, cov_probit = irls(Xs, prop, binomial_probit(), w=m)
sm_probit = sm.GLM(np.column_stack([succ, m - succ]), Xs,
                   family=sm.families.Binomial(sm.families.links.Probit())).fit()
assert np.allclose(b_probit, sm_probit.params, atol=1e-8)
assert np.allclose(np.sqrt(np.diag(cov_probit)), sm_probit.bse, atol=1e-8)

# gamma with log link on Engel's data; here the dispersion is estimated afterwards
engel = sm.datasets.engel.load_pandas().data
y_e = engel["foodexp"].to_numpy()
X_e = np.column_stack([np.ones(len(y_e)), np.log(engel["income"].to_numpy())])
b_gam, cov_gam = irls(X_e, y_e, gamma_log())
sm_gam = sm.GLM(y_e, X_e, family=sm.families.Gamma(sm.families.links.Log())).fit()
assert np.allclose(b_gam, sm_gam.params, atol=1e-9)
mu_e = np.exp(X_e @ b_gam)
phi_hat = np.sum((y_e - mu_e) ** 2 / mu_e ** 2) / (len(y_e) - 2)     # Pearson estimate
assert np.isclose(phi_hat, sm_gam.scale)
se_gam = np.sqrt(phi_hat * np.diag(cov_gam))
assert np.allclose(se_gam, sm_gam.bse, atol=1e-8)

# ---- how fast does it converge? ---------------------------------------------
paths = {}
for fam, data in [(poisson_log(), (Xs, counts.astype(float), None)),
                  (binomial_logit(), (Xs, prop, m)),
                  (binomial_probit(), (Xs, prop, m)),
                  (gamma_log(), (X_e, y_e, None))]:
    XX, yy, ww = data
    bb, _, path = irls(XX, yy, fam, w=ww, trace=True, max_iter=40)
    paths[fam["name"]] = np.max(np.abs(path - bb), axis=1)

steps = {name: int(np.argmax(err < 1e-10)) + 1 for name, err in paths.items()}
print("iterations to 1e-10:", steps)

# ---- a failure mode: separated binary data ---------------------------------
# <<separation>>
x_sep = np.array([-2.4, -1.7, -1.1, -0.6, -0.2, 0.3, 0.8, 1.4, 1.9, 2.6])
y_sep = (x_sep > 0).astype(float)                        # the two groups do not overlap
X_sep = np.column_stack([np.ones(10), x_sep])

for iters in range(1, 9):
    with np.errstate(divide="ignore", invalid="ignore"):
        b, _ = irls(X_sep, y_sep, binomial_logit(), tol=0.0, max_iter=iters)
        mu = 1 / (1 + np.exp(-(X_sep @ b)))
        dev = -2 * np.sum(np.where(y_sep == 1, np.log(mu), np.log(1 - mu)))
    print(f"{iters} steps: slope {b[1]:9.4f}   deviance {dev:.3e}"
          f"   largest working weight {np.max(mu * (1 - mu)):.3e}")
# <</separation>>

slopes, devs, wmax = [], [], []
for k in range(1, 9):
    with np.errstate(divide="ignore", invalid="ignore"):
        b, _ = irls(X_sep, y_sep, binomial_logit(), tol=0.0, max_iter=k)
        mu = 1 / (1 + np.exp(-(X_sep @ b)))
        slopes.append(b[1])
        devs.append(-2 * np.sum(np.where(y_sep == 1, np.log(mu), np.log(1 - mu))))
        wmax.append(np.max(mu * (1 - mu)))
assert all(slopes[k] < slopes[k + 1] for k in range(6))   # the slope never settles
assert all(devs[k] > devs[k + 1] for k in range(6))       # the deviance keeps falling
assert np.isnan(slopes[7])                               # step 8: the weights underflow

gen = Generated("ch34", "irls")
gen.int("n_sim", n)
gen.num("pois_b1", beta_hat[1], 4)
gen.num("pois_se1", float(np.sqrt(cov[1, 1])), 4)
gen.num("gam_b0", b_gam[0], 4)
gen.num("gam_b1", b_gam[1], 4)
gen.num("gam_phi", float(phi_hat), 4)
gen.num("gam_se1", float(se_gam[1]), 4)
gen.int("steps_pois", steps["Poisson/log"])
gen.int("steps_logit", steps["binomial/logit"])
gen.int("steps_probit", steps["binomial/probit"])
gen.int("steps_gamma", steps["gamma/log"])
gen.num("sep_slope1", slopes[0], 2)
gen.num("sep_slope4", slopes[3], 2)
gen.num("sep_slope7", slopes[6], 2)
gen.num("sep_dev1", devs[0], 3)
gen.num("sep_dev4", devs[3], 3)
gen.num("sep_w1", wmax[0], 3)
gen.num("sep_w4", wmax[3], 3)
gen.num("sep_dev7", devs[6], 4)
gen.num("sep_w7", wmax[6], 4)
gen.write()

# ---- figure -----------------------------------------------------------------
use_book_style()
fig, ax = plt.subplots(figsize=(4.2, 2.6))
styles = {"Poisson/log": (COLORS["accent"], "-", "o"),
          "binomial/logit": (COLORS["third"], "-", "s"),
          "binomial/probit": (COLORS["second"], "--", "^"),
          "gamma/log": (COLORS["thread"], "--", "d")}
for name, err in paths.items():
    col, ls, mk = styles[name]
    keep = err > 0
    ax.plot(np.arange(1, len(err) + 1)[keep], err[keep], ls, color=col, marker=mk,
            markersize=3, label=name)
ax.set_yscale("log")
ax.set_ylim(1e-14, 1e1)
ax.set_xlabel("iteration $t$")
ax.set_ylabel(r"$\max_j|\beta_j^{(t)}-\hat\beta_j|$")
ax.legend(frameon=False, fontsize=7, loc="lower left")
fig.tight_layout()
fig.savefig(figure_path("ch34", "irls_convergence"))
