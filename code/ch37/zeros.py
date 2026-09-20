"""Chapter 37, Section 5: zero-inflated and hurdle models for the physician-visit counts.

Data: statsmodels.datasets.randhie (public domain), 20190 person-years, of which
31.2 percent report no visit.

Fits, in numpy and scipy: Poisson, NB2, zero-inflated Poisson and zero-inflated NB2
(both by EM), and a hurdle model (logistic for the zeros, zero-truncated NB2 for the
positives). Compares the fitted marginal distributions of the count with the
observed one, computes AIC and the Vuong statistic, and draws the comparison.
"""
import matplotlib.pyplot as plt
import numpy as np
import statsmodels.api as sm
from scipy import optimize, special, stats

from regbook import COLORS, Generated, figure_path, use_book_style

data = sm.datasets.randhie.load_pandas().data
y = data["mdvis"].to_numpy(float)
names = ["lncoins", "idp", "physlm", "disea", "hlthg", "hlthf", "hlthp"]
X = np.column_stack([np.ones(len(y))] + [data[v].to_numpy(float) for v in names])
n, p = X.shape
log_factorial = special.gammaln(y + 1)
positive = y > 0


# ---- building blocks --------------------------------------------------------
# <<pieces>>
import numpy as np
from scipy import optimize, special, stats


def poisson_irls(X, y, w=None, start=None, tol=1e-11, maxit=200):
    """Weighted Poisson regression with a log link; w are prior weights."""
    w = np.ones(len(y)) if w is None else w
    beta = np.zeros(X.shape[1]) if start is None else start.copy()
    if start is None:
        beta[0] = np.log(max((w * y).sum() / w.sum(), 1e-8))
    for _ in range(maxit):
        mu = np.exp(X @ beta)
        XW = X * (w * mu)[:, None]
        step = np.linalg.solve(X.T @ XW, XW.T @ (X @ beta + (y - mu) / mu))
        if np.max(np.abs(step - beta)) < tol:
            return step
        beta = step
    return beta


def logistic_irls(X, r, tol=1e-11, maxit=200):
    """Logistic regression with a response r in [0, 1] (a fractional response is allowed)."""
    gamma = np.zeros(X.shape[1])
    for _ in range(maxit):
        pi = 1 / (1 + np.exp(-X @ gamma))
        v = np.clip(pi * (1 - pi), 1e-10, None)
        XW = X * v[:, None]
        step = np.linalg.solve(X.T @ XW, XW.T @ (X @ gamma + (r - pi) / v))
        if np.max(np.abs(step - gamma)) < tol:
            return step
        gamma = step
    return gamma


def nb_logpmf(y, mu, kappa):
    """log P(Y = y) for the negative binomial with mean mu and shape kappa."""
    return (special.gammaln(y + kappa) - special.gammaln(kappa) - special.gammaln(y + 1)
            + kappa * np.log(kappa / (kappa + mu)) + y * np.log(mu / (mu + kappa)))


def poisson_logpmf(y, mu):
    return y * np.log(mu) - mu - special.gammaln(y + 1)
# <</pieces>>


# <<nbstep>>
def nb_step(X, y, beta, kappa, w):
    """One IRLS step in beta and one Newton step in log kappa, at prior weights w."""
    mu = np.exp(X @ beta)
    XW = X * (w * mu / (1 + mu / kappa))[:, None]        # NB2 working weights, log link
    beta = np.linalg.solve(X.T @ XW, XW.T @ (X @ beta + (y - mu) / mu))
    mu = np.exp(X @ beta)
    s = np.sum(w * (special.digamma(y + kappa) - special.digamma(kappa)
                    + np.log(kappa / (kappa + mu)) + 1 - (kappa + y) / (kappa + mu)))
    d = np.sum(w * (special.polygamma(1, y + kappa) - special.polygamma(1, kappa)
                    + 1 / kappa - 1 / (kappa + mu) - (mu - y) / (kappa + mu) ** 2))
    return beta, float(np.clip(np.exp(np.log(kappa) - s / (kappa * d)), 1e-4, 1e4))


def nb_fit(X, y, w=None, start=None, tol=1e-10, maxit=500):
    """Weighted NB2 maximum likelihood: iterate nb_step to convergence."""
    w = np.ones(len(y)) if w is None else w
    beta, kappa = (poisson_irls(X, y, w), 1.0) if start is None else start
    for _ in range(maxit):
        new_beta, new_kappa = nb_step(X, y, beta, kappa, w)
        done = max(np.abs(new_beta - beta).max(), abs(new_kappa - kappa)) < tol
        beta, kappa = new_beta, new_kappa
        if done:
            break
    return beta, kappa
# <</nbstep>>


# ---- the two baseline fits --------------------------------------------------
beta_pois = poisson_irls(X, y)
loglik_pois = poisson_logpmf(y, np.exp(X @ beta_pois)).sum()
beta_nb, kappa_nb = nb_fit(X, y)
loglik_nb = nb_logpmf(y, np.exp(X @ beta_nb), kappa_nb).sum()
print(f"Poisson {loglik_pois:.1f}   NB2 {loglik_nb:.1f} (kappa {kappa_nb:.3f})")

# ---- zero-inflated fits by EM ----------------------------------------------
# <<em>>
def zero_inflated_em(X, Z, y, negbin=False, tol=1e-10, maxit=2000):
    """EM for a zero-inflated Poisson or NB2 model; returns (gamma, beta, kappa, loglik).

    gamma indexes the logistic model for the degenerate component, beta the log-linear
    model for the count component. The latent indicator is z_i = 1 if observation i
    comes from the degenerate component at zero.
    """
    zero = y == 0
    gamma = np.zeros(Z.shape[1])
    beta, kappa = (poisson_irls(X, y), np.inf)
    if negbin:
        beta, kappa = nb_fit(X, y)
    previous = -np.inf
    for _ in range(maxit):
        mu = np.exp(X @ beta)
        pi = 1 / (1 + np.exp(-Z @ gamma))
        log_count = nb_logpmf(y, mu, kappa) if negbin else poisson_logpmf(y, mu)
        # E step: the posterior probability that a zero came from the degenerate part
        w = np.where(zero, pi / (pi + (1 - pi) * np.exp(log_count)), 0.0)
        # M step: a fractional logistic regression and a weighted count regression
        gamma = logistic_irls(Z, w)
        if negbin:
            beta, kappa = nb_step(X, y, beta, kappa, 1 - w)   # one step is enough
        else:
            beta = poisson_irls(X, y, 1 - w, start=beta)
        pi = 1 / (1 + np.exp(-Z @ gamma))
        mu = np.exp(X @ beta)
        log_count = nb_logpmf(y, mu, kappa) if negbin else poisson_logpmf(y, mu)
        loglik = np.sum(np.where(zero,
                                 np.log(pi + (1 - pi) * np.exp(log_count)),
                                 np.log1p(-pi) + log_count))
        if loglik - previous < tol * abs(loglik):
            break
        previous = loglik
    return gamma, beta, kappa, loglik


gamma_zip, beta_zip, _, loglik_zip = zero_inflated_em(X, X, y)
gamma_zinb, beta_zinb, kappa_zinb, loglik_zinb = zero_inflated_em(X, X, y, negbin=True)
print(f"ZIP  log-likelihood {loglik_zip:.2f}")
print(f"ZINB log-likelihood {loglik_zinb:.2f}   kappa {kappa_zinb:.3f}")
# <</em>>

assert loglik_zip > loglik_pois
assert loglik_zinb > max(loglik_nb, loglik_zip)

# the M step takes one weighted step rather than a full refit: a generalized EM, whose
# observed-data log-likelihood still increases at every iteration
path = [zero_inflated_em(X, X, y, negbin=True, maxit=k)[3] for k in range(1, 15)]
assert np.all(np.diff(path) > 0)

# ---- the hurdle model -------------------------------------------------------
# <<hurdle>>
def truncated_nb_loglik(theta, X, y):
    """Log-likelihood of the zero-truncated NB2 model, parameters (beta, log kappa)."""
    mu = np.exp(np.clip(X @ theta[:-1], -20.0, 20.0))
    k = np.exp(theta[-1])
    log_p0 = k * np.log(k / (k + mu))                 # log P(Y = 0) before truncation
    return np.sum(nb_logpmf(y, mu, k) - np.log(-np.expm1(log_p0)))


def hurdle_nb(X, y):
    """Logistic model for y > 0, zero-truncated NB2 for the positive counts."""
    gamma = logistic_irls(X, (y > 0).astype(float))
    pi = 1 / (1 + np.exp(-X @ gamma))
    loglik_zero = np.sum(np.where(y > 0, np.log(pi), np.log1p(-pi)))

    Xp, yp = X[y > 0], y[y > 0]
    beta0, kappa0 = nb_fit(Xp, yp)                    # an untruncated fit as the start
    out = optimize.minimize(lambda t: -truncated_nb_loglik(t, Xp, yp),
                            np.append(beta0, np.log(kappa0)), method="BFGS",
                            options={"gtol": 1e-5, "maxiter": 800})
    return gamma, out.x[:-1], np.exp(out.x[-1]), loglik_zero - out.fun


gamma_h, beta_h, kappa_h, loglik_hurdle = hurdle_nb(X, y)
print(f"hurdle log-likelihood {loglik_hurdle:.2f}   kappa {kappa_h:.3f}")
# <</hurdle>>

assert loglik_hurdle > loglik_nb

# the hurdle optimum survives a derivative-free refinement
_Xp, _yp = X[positive], y[positive]
_theta = np.append(beta_h, np.log(kappa_h))
_refined = optimize.minimize(lambda t: -truncated_nb_loglik(t, _Xp, _yp), _theta,
                             method="Nelder-Mead",
                             options={"maxiter": 20000, "fatol": 1e-8, "xatol": 1e-8})
assert -_refined.fun - truncated_nb_loglik(_theta, _Xp, _yp) < 1e-6

# ---- fitted marginal distribution of the count ------------------------------
K = 11
grid = np.arange(K)


def average_pmf(logpmf_at):
    return np.array([np.mean(np.exp(logpmf_at(k))) for k in grid])


mu_p, mu_nb = np.exp(X @ beta_pois), np.exp(X @ beta_nb)
pi_zinb = 1 / (1 + np.exp(-X @ gamma_zinb))
mu_zinb = np.exp(X @ beta_zinb)
pi_h = 1 / (1 + np.exp(-X @ gamma_h))
mu_h = np.exp(X @ beta_h)
p0_h = np.exp(kappa_h * np.log(kappa_h / (kappa_h + mu_h)))

fitted = {
    "Poisson": average_pmf(lambda k: poisson_logpmf(k, mu_p)),
    "NB2": average_pmf(lambda k: nb_logpmf(k, mu_nb, kappa_nb)),
    "ZINB": np.array([np.mean(np.where(k == 0, pi_zinb, 0.0)
                              + (1 - pi_zinb) * np.exp(nb_logpmf(k, mu_zinb, kappa_zinb)))
                      for k in grid]),
    "hurdle NB": np.array([np.mean(1 - pi_h) if k == 0 else
                           np.mean(pi_h * np.exp(nb_logpmf(k, mu_h, kappa_h)) / (1 - p0_h))
                           for k in grid]),
}
observed = np.array([np.mean(y == k) for k in grid])

pi_zip, mu_zip = 1 / (1 + np.exp(-X @ gamma_zip)), np.exp(X @ beta_zip)
zip_zero = np.mean(pi_zip + (1 - pi_zip) * np.exp(-mu_zip))
assert abs(zip_zero - observed[0]) < 0.01

assert fitted["Poisson"][0] < 0.3 * observed[0]            # the Poisson misses the zeros
assert abs(fitted["NB2"][0] - observed[0]) < 0.01          # the NB2 already matches them
assert abs(fitted["hurdle NB"][0] - observed[0]) < 1e-10   # the hurdle fits the zeros exactly
assert loglik_hurdle > loglik_zinb > loglik_nb
# above two visits the hurdle and ZINB fitted frequencies are visually indistinguishable
assert np.max(np.abs(fitted["hurdle NB"][3:] - fitted["ZINB"][3:])) < 0.005

# ---- Vuong statistics -------------------------------------------------------
# <<vuong>>
def vuong(log_f1, log_f2):
    """Vuong's statistic for model 1 against model 2 (not valid for nested models)."""
    m = log_f1 - log_f2
    return np.sqrt(len(m)) * m.mean() / m.std(ddof=1)


log_zinb = np.where(y == 0, np.log(pi_zinb + (1 - pi_zinb) * np.exp(nb_logpmf(0.0, mu_zinb, kappa_zinb))),
                    np.log1p(-pi_zinb) + nb_logpmf(y, mu_zinb, kappa_zinb))
log_hurdle = np.where(y == 0, np.log1p(-pi_h),
                      np.log(pi_h) + nb_logpmf(y, mu_h, kappa_h) - np.log1p(-p0_h))
log_nb = nb_logpmf(y, mu_nb, kappa_nb)

print(f"ZINB against NB2 (nested, so not valid):  V = {vuong(log_zinb, log_nb):6.2f}")
print(f"ZINB against the hurdle model:            V = {vuong(log_zinb, log_hurdle):6.2f}")
# <</vuong>>

assert np.isclose(log_zinb.sum(), loglik_zinb)
assert np.isclose(log_hurdle.sum(), loglik_hurdle)

v_nested = vuong(log_zinb, log_nb)
v_hurdle = vuong(log_zinb, log_hurdle)

# where the hurdle model's advantage over the ZINB comes from: split each model's
# log-likelihood into the part contributed by the zeros and the part by the positives
zero_part = {"hurdle": log_hurdle[~positive].sum(), "zinb": log_zinb[~positive].sum()}
positive_part = {"hurdle": log_hurdle[positive].sum(), "zinb": log_zinb[positive].sum()}
assert zero_part["hurdle"] < zero_part["zinb"]          # the hurdle loses on the zeros
assert positive_part["hurdle"] > positive_part["zinb"]  # and wins on the positive counts
assert (positive_part["hurdle"] - positive_part["zinb"]
        > zero_part["zinb"] - zero_part["hurdle"])      # the second gain is the larger

k_pars = {"Poisson": p, "NB2": p + 1, "ZIP": 2 * p, "ZINB": 2 * p + 1,
          "hurdle NB": 2 * p + 1}
logliks = {"Poisson": loglik_pois, "NB2": loglik_nb, "ZIP": loglik_zip,
           "ZINB": loglik_zinb, "hurdle NB": loglik_hurdle}
aic = {name: -2 * logliks[name] + 2 * k_pars[name] for name in logliks}
bic = {name: -2 * logliks[name] + np.log(n) * k_pars[name] for name in logliks}
assert min(aic, key=aic.get) == "hurdle NB"
assert min(bic, key=bic.get) == "hurdle NB"

gen = Generated("ch37", "zeros")
gen.int("n", n)
gen.num("zero_fraction", observed[0], 4)
for name in logliks:
    tag = name.replace(" ", "_").lower()
    gen.num("loglik_" + tag, logliks[name], 1)
    gen.num("aic_" + tag, aic[name], 1)
    gen.num("bic_" + tag, bic[name], 1)
    gen.int("k_" + tag, k_pars[name])
    if name in fitted:
        gen.num("zero_" + tag, fitted[name][0], 4)
gen.num("zero_zip", zip_zero, 4)
gen.num("kappa_nb", kappa_nb, 3)
gen.num("kappa_zinb", kappa_zinb, 3)
gen.num("kappa_hurdle", kappa_h, 3)
gen.num("v_nested", v_nested, 2)
gen.num("v_hurdle", v_hurdle, 2)
gen.num("lr_zinb_nb", 2 * (loglik_zinb - loglik_nb), 1)
gen.num("mean_pi_zinb", pi_zinb.mean(), 4)
for tag in ("hurdle", "zinb"):
    gen.num("zeropart_" + tag, zero_part[tag], 1)
    gen.num("pospart_" + tag, positive_part[tag], 1)
gen.num("zeropart_loss", zero_part["zinb"] - zero_part["hurdle"], 1)
gen.num("pospart_gain", positive_part["hurdle"] - positive_part["zinb"], 1)
gen.num("net_gain", loglik_hurdle - loglik_zinb, 1)
for k in range(4):
    gen.num(f"obs_{k}", observed[k], 4)
    for name in fitted:
        gen.num(f"fit_{name.replace(' ', '_').lower()}_{k}", fitted[name][k], 4)
gen.write()

# ---- the figure -------------------------------------------------------------
use_book_style()
fig, ax = plt.subplots(figsize=(5.4, 3.0))
ax.bar(grid, observed, width=0.72, color=COLORS["grid"], label="observed", zorder=1)
styles = [("Poisson", "accent", "-"), ("NB2", "third", (0, (4, 2))),
          ("ZINB", "second", (0, (1, 1.4))), ("hurdle NB", "thread", (0, (5, 1.5, 1, 1.5)))]
for name, colour, dash in styles:
    ax.plot(grid, fitted[name], linestyle=dash, color=COLORS[colour], marker="o",
            markersize=2.6, linewidth=1.0, label=name, zorder=3)
ax.set_xlabel("number of physician visits")
ax.set_ylabel("proportion of person-years")
ax.set_xticks(grid)
ax.legend(frameon=False, fontsize=8)
fig.tight_layout()
fig.savefig(figure_path("ch37", "zero_inflation"))
