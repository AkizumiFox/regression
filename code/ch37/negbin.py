"""Chapter 37, Section 4: negative binomial models for the physician-visit counts.

Data: statsmodels.datasets.randhie (public domain), 20190 person-years.

Fits the Poisson, NB2 and NB1 models by maximum likelihood in numpy/scipy, compares
their standard errors with the quasi-Poisson and sandwich alternatives, runs the
likelihood ratio and score tests for overdispersion, and draws the three variance
functions against binned sample variances.
"""
import matplotlib.pyplot as plt
import numpy as np
import statsmodels.api as sm
from scipy import optimize, special, stats

from regbook import COLORS, Generated, figure_path, use_book_style

data = sm.datasets.randhie.load_pandas().data
y = data["mdvis"].to_numpy(float)
names = ["intercept", "lncoins", "idp", "physlm", "disea", "hlthg", "hlthf", "hlthp"]
X = np.column_stack([np.ones(len(y))] + [data[v].to_numpy(float) for v in names[1:]])
n, p = X.shape

# ---- the Poisson fit, for comparison ---------------------------------------
beta_pois = np.zeros(p)
beta_pois[0] = np.log(y.mean())
for _ in range(50):
    mu = np.exp(X @ beta_pois)
    XW = X * mu[:, None]
    beta_pois = np.linalg.solve(X.T @ XW, XW.T @ (X @ beta_pois + (y - mu) / mu))
mu_pois = np.exp(X @ beta_pois)
info_pois = X.T @ (X * mu_pois[:, None])
se_pois = np.sqrt(np.diag(np.linalg.inv(info_pois)))
log_factorial = special.gammaln(y + 1)
loglik_pois = np.sum(y * np.log(mu_pois) - mu_pois - log_factorial)
print(f"Poisson log-likelihood {loglik_pois:.2f}")

# <<nb2>>
import numpy as np
from scipy import optimize, special, stats


def nb_loglik(beta, kappa, X, y):
    """Negative binomial log-likelihood: mean exp(X beta), shape kappa (NB2 if scalar)."""
    mu = np.exp(X @ beta)
    return np.sum(special.gammaln(y + kappa) - special.gammaln(kappa)
                  - special.gammaln(y + 1)
                  + kappa * np.log(kappa / (kappa + mu))
                  + y * np.log(mu / (mu + kappa)))


def kappa_score(kappa, mu, y):
    """Derivative of the log-likelihood in the shape parameter."""
    return np.sum(special.digamma(y + kappa) - special.digamma(kappa)
                  + np.log(kappa / (kappa + mu)) + 1 - (kappa + y) / (kappa + mu))


def fit_nb2(X, y, tol=1e-11, maxit=200):
    """Alternate IRLS in beta with a one-dimensional search in kappa."""
    beta = np.zeros(X.shape[1])
    beta[0] = np.log(y.mean())
    kappa = 1.0
    for _ in range(maxit):
        mu = np.exp(X @ beta)
        w = mu / (1 + mu / kappa)                  # NB2 working weights, log link
        XW = X * w[:, None]
        step = np.linalg.solve(X.T @ XW, XW.T @ (X @ beta + (y - mu) / mu))
        new_kappa = optimize.brentq(kappa_score, 1e-4, 1e4,
                                    args=(np.exp(X @ step), y), xtol=1e-12)
        done = max(np.abs(step - beta).max(), abs(new_kappa - kappa)) < tol
        beta, kappa = step, new_kappa
        if done:
            break
    return beta, kappa, nb_loglik(beta, kappa, X, y)


beta_nb2, kappa, loglik_nb2 = fit_nb2(X, y)
mu_nb2 = np.exp(X @ beta_nb2)
weights = mu_nb2 / (1 + mu_nb2 / kappa)            # IRLS weights for NB2 with a log link
se_nb2 = np.sqrt(np.diag(np.linalg.inv(X.T @ (X * weights[:, None]))))

print(f"kappa = {kappa:.4f}   log-likelihood {loglik_nb2:.2f}")
for name, b, s, sp in zip(names, beta_nb2, se_nb2, se_pois):
    print(f"{name:10s} {b:8.4f}  se {s:.4f}  (Poisson se {sp:.4f})")
# <</nb2>>

assert np.abs(X.T @ ((y - mu_nb2) / (1 + mu_nb2 / kappa))).max() < 1e-3
assert np.all(se_nb2 > se_pois)

# ---- the two tests for overdispersion --------------------------------------
# <<tests>>
lr = 2 * (loglik_nb2 - loglik_pois)
# the null puts 1/kappa on the boundary, so the null distribution of L is the mixture
# (1/2) chi2(0) + (1/2) chi2(1) and the p-value is the normal tail at sqrt(L)
log_p = stats.norm.logsf(np.sqrt(lr))

score = np.sum((y - mu_pois) ** 2 - y) / np.sqrt(2 * np.sum(mu_pois ** 2))
print(f"likelihood ratio {lr:.1f}, log p-value {log_p:.0f}")
print(f"score statistic {score:.1f}, log p-value {stats.norm.logsf(score):.0f}")
# <</tests>>

assert lr > 0 and score > 0

# ---- NB1 --------------------------------------------------------------------
def fit_nb1(X, y, start):
    """NB1: variance phi * mu, that is shape kappa_i = alpha * mu_i."""
    def objective(theta):
        beta, alpha = theta[:-1], np.exp(theta[-1])
        mu = np.exp(X @ beta)
        k = alpha * mu
        return -np.sum(special.gammaln(y + k) - special.gammaln(k) - special.gammaln(y + 1)
                       + k * np.log(k / (k + mu)) + y * np.log(mu / (mu + k)))

    out = optimize.minimize(objective, np.append(start, 0.0), method="BFGS",
                            options={"gtol": 1e-7, "maxiter": 500})
    return out.x[:-1], np.exp(out.x[-1]), -out.fun


beta_nb1, alpha1, loglik_nb1 = fit_nb1(X, y, beta_pois)
phi_nb1 = 1 + 1 / alpha1                            # variance = phi * mu
assert loglik_nb1 > loglik_pois


def nb1_loglik(beta, alpha):
    mu = np.exp(X @ beta)
    k = alpha * mu
    return np.sum(special.gammaln(y + k) - special.gammaln(k) - special.gammaln(y + 1)
                  + k * np.log(k / (k + mu)) + y * np.log(mu / (mu + k)))


# the reported NB1 fit is a local maximum
for factor in (0.99, 1.01):
    assert nb1_loglik(beta_nb1, alpha1 * factor) < loglik_nb1
    for j in range(p):
        shifted = beta_nb1.copy()
        shifted[j] *= factor
        assert nb1_loglik(shifted, alpha1) <= loglik_nb1 + 1e-9

# check the NB1 score of prp-cnt-negbin(f) against a numerical derivative, and confirm
# that the NB1 estimate differs from the Poisson one, which is also the quasi-Poisson one
mu_nb1 = np.exp(X @ beta_nb1)
nb1_score = alpha1 * X.T @ (mu_nb1 * (special.digamma(y + alpha1 * mu_nb1)
                                      - special.digamma(alpha1 * mu_nb1)
                                      + np.log(alpha1 / (1 + alpha1))))
eps = 1e-6
numeric = np.array([(nb1_loglik(beta_nb1 + eps * e, alpha1)
                     - nb1_loglik(beta_nb1 - eps * e, alpha1)) / (2 * eps)
                    for e in np.eye(p)])
assert np.abs(nb1_score - numeric).max() < 1e-3
assert np.abs(beta_nb1 - beta_pois).max() > 1e-3

# ---- quasi-Poisson and sandwich standard errors ----------------------------
pearson = np.sum((y - mu_pois) ** 2 / mu_pois)
phi_hat = pearson / (n - p)
se_quasi = se_pois * np.sqrt(phi_hat)
meat = X.T @ (X * ((y - mu_pois) ** 2)[:, None])
inv = np.linalg.inv(info_pois)
se_sandwich = np.sqrt(np.diag(inv @ meat @ inv))

aic = {"Poisson": -2 * loglik_pois + 2 * p,
       "NB2": -2 * loglik_nb2 + 2 * (p + 1),
       "NB1": -2 * loglik_nb1 + 2 * (p + 1)}
assert aic["NB1"] < aic["NB2"] < aic["Poisson"]

j = names.index("idp")
gen = Generated("ch37", "negbin")
gen.int("n", n)
gen.int("p", p)
gen.num("kappa", kappa, 4)
gen.num("alpha", 1 / kappa, 4)
gen.num("loglik_pois", loglik_pois, 2)
gen.num("loglik_nb2", loglik_nb2, 2)
gen.num("loglik_nb1", loglik_nb1, 2)
gen.num("lr", lr, 1)
gen.num("score", score, 1)
gen.num("log_p_lr", log_p, 0)
gen.num("phi", phi_hat, 3)
gen.num("phi_nb1", phi_nb1, 3)
gen.num("root_phi", np.sqrt(phi_hat), 3)
gen.num("b_idp_pois", beta_pois[j], 4)
gen.num("b_idp_nb2", beta_nb2[j], 4)
gen.num("se_idp_pois", se_pois[j], 4)
gen.num("se_idp_nb2", se_nb2[j], 4)
gen.num("se_idp_quasi", se_quasi[j], 4)
gen.num("se_idp_sandwich", se_sandwich[j], 4)
gen.num("ratio_nb2_pois", se_nb2[j] / se_pois[j], 3)
gen.num("b_disea_nb2", beta_nb2[names.index("disea")], 4)
gen.num("rr_disea_nb2", np.exp(beta_nb2[names.index("disea")]), 4)
gen.num("aic_pois", aic["Poisson"], 1)
gen.num("aic_nb2", aic["NB2"], 1)
gen.num("aic_nb1", aic["NB1"], 1)
gen.num("var_at_mean", y.mean() + y.mean() ** 2 / kappa, 3)
gen.num("mean", y.mean(), 3)
gen.write()

# ---- the variance functions against binned sample variances ----------------
order = np.argsort(mu_pois)
bins = np.array_split(order, 25)
bin_mean = np.array([y[b].mean() for b in bins])
bin_var = np.array([y[b].var(ddof=1) for b in bins])
bin_fitted = np.array([mu_pois[b].mean() for b in bins])
assert np.all(bin_var > bin_mean * 0.9)

# the text's reading of the figure: NB1 is the better description of the crowded low
# range of fitted means, NB2 of the upper end
nb1_curve, nb2_curve = phi_nb1 * bin_fitted, bin_fitted + bin_fitted ** 2 / kappa
low = bin_fitted < np.median(bin_fitted)
assert np.abs(bin_var - nb1_curve)[low].mean() < np.abs(bin_var - nb2_curve)[low].mean()
assert np.abs(bin_var - nb2_curve)[~low].mean() < np.abs(bin_var - nb1_curve)[~low].mean()

use_book_style()
fig, ax = plt.subplots(figsize=(4.4, 3.1))
grid = np.linspace(bin_fitted.min() * 0.9, bin_fitted.max() * 1.05, 200)
ax.scatter(bin_fitted, bin_var, s=14, color=COLORS["ink"], zorder=3,
           label="binned sample variance")
ax.plot(grid, grid, color=COLORS["accent"], label=r"Poisson: $V(\mu)=\mu$")
ax.plot(grid, phi_nb1 * grid, color=COLORS["third"], linestyle=(0, (4, 2)),
        label=r"NB1: $V(\mu)=\phi\mu$")
ax.plot(grid, grid + grid ** 2 / kappa, color=COLORS["second"], linestyle=(0, (1, 1.6)),
        label=r"NB2: $V(\mu)=\mu+\mu^2/\kappa$")
ax.set_xlabel(r"fitted mean $\hat\mu$ (bin average)")
ax.set_ylabel("sample variance in the bin")
ax.legend(frameon=False, fontsize=7.5, loc="upper left")
fig.tight_layout()
fig.savefig(figure_path("ch37", "variance_functions"))
