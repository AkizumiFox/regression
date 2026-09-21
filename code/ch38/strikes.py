"""Chapter 38: the worked example. Strike durations against unanticipated output.

62 US manufacturing strikes (Kennan 1985, from Bureau of Labor Statistics records;
public-domain data shipped with statsmodels as statsmodels.datasets.strikes). The
response is the duration of the strike in days, the regressor is a measure of
unanticipated industrial production. A Poisson log-linear mean model is fitted by the
quasi-score equations, and the four standard covariance estimates are compared.
"""
import numpy as np
import statsmodels.api as sm
from scipy import stats

from regbook import Generated

# <<fit>>
data = sm.datasets.strikes.load_pandas().data
y = data["duration"].to_numpy(float)          # length of the strike, in days
x = data["iprod"].to_numpy(float)             # unanticipated industrial production
X = np.column_stack([np.ones(len(y)), x])
n, p = X.shape


def quasi_irls(X, y, V, steps=60):
    """Solve the quasi-score equations for a log link and variance function V."""
    beta = np.zeros(X.shape[1])
    beta[0] = np.log(y.mean())
    for _ in range(steps):
        eta = X @ beta
        mu = np.exp(eta)
        W = mu ** 2 / V(mu)                   # working weight (dmu/deta)^2 / V(mu)
        z = eta + (y - mu) / mu               # working response
        beta = np.linalg.solve((X.T * W) @ X, (X.T * W) @ z)
    return beta


beta_lin = quasi_irls(X, y, lambda m: m)              # V(mu) = mu
beta_scaled = quasi_irls(X, y, lambda m: 7.5 * m)     # V(mu) = phi mu: same equations
beta_sq = quasi_irls(X, y, lambda m: m ** 2)          # V(mu) = mu^2: different equations
print(f"V = mu       {beta_lin}")
print(f"V = 7.5 mu   {beta_scaled}")
print(f"V = mu^2     {beta_sq}")
# <</fit>>

assert np.allclose(beta_lin, beta_scaled, atol=1e-10)
assert not np.allclose(beta_lin, beta_sq, atol=1e-3)
check = sm.GLM(y, X, family=sm.families.Poisson()).fit()
assert np.allclose(check.params, beta_lin, atol=1e-8)

# <<dispersion>>
mu = np.exp(X @ beta_lin)
pearson = np.sum((y - mu) ** 2 / mu)                          # X^2 with V(mu) = mu
deviance = 2 * np.sum(y * np.log(y / mu) - (y - mu))          # quasi-deviance, V(mu) = mu
phi_pearson = pearson / (n - p)
phi_deviance = deviance / (n - p)
score_stat = np.sum((y - mu) ** 2 - y) / np.sqrt(2 * np.sum(mu ** 2))
print(f"Pearson X^2 {pearson:.1f} on {n - p} df, phi_hat = {phi_pearson:.3f}")
print(f"deviance    {deviance:.1f} on {n - p} df, phi_hat = {phi_deviance:.3f}")
print(f"score statistic for overdispersion: {score_stat:.1f}")
# <</dispersion>>

# <<covariances>>
W = mu                                        # log link with V(mu) = mu: W_ii = mu_i
bread = np.linalg.inv((X.T * W) @ X)          # (X' W X)^{-1}
root = np.sqrt(W)[:, None] * X
hat = np.sum(root * (root @ bread), axis=1)   # leverages of the weighted hat matrix
resid = y - mu

cov_model = bread                             # phi = 1: the Poisson standard errors
cov_quasi = phi_pearson * bread               # phi estimated from the same fit
cov_sand = bread @ ((X * resid[:, None] ** 2).T @ X) @ bread
cov_corr = bread @ ((X * (resid ** 2 / (1 - hat) ** 2)[:, None]).T @ X) @ bread

for name, C in [("Poisson", cov_model), ("quasi", cov_quasi),
                ("sandwich", cov_sand), ("corrected", cov_corr)]:
    se = np.sqrt(np.diag(C))
    print(f"{name:10s} slope {beta_lin[1]:8.3f}   se {se[1]:6.3f}   t {beta_lin[1] / se[1]:7.2f}")
print(f"largest leverage {hat.max():.4f}, sum {hat.sum():.4f}")
# <</covariances>>

assert np.isclose(hat.sum(), p)
se_model, se_quasi = np.sqrt(np.diag(cov_model)), np.sqrt(np.diag(cov_quasi))
se_sand, se_corr = np.sqrt(np.diag(cov_sand)), np.sqrt(np.diag(cov_corr))
assert np.isclose(se_quasi[1] / se_model[1], np.sqrt(phi_pearson))
assert se_corr[1] > se_sand[1]
# the sandwich agrees with the model-based error: evidence for V(mu) = mu (Section 38.3)
assert abs(se_sand[1] / se_quasi[1] - 1) < 0.05

# the largest count and the raw residual it contributes to the quasi-score (Section 38.4)
imax = int(np.argmax(y))
ymax, mu_at_ymax = y[imax], mu[imax]
print(f"largest count {ymax:.0f} at fitted mean {mu_at_ymax:.2f},"
      f" raw residual {ymax - mu_at_ymax:.1f}")
assert mu_at_ymax > mu.min() and ymax - mu_at_ymax > 100

# <<tests>>
beta_null = quasi_irls(X[:, :1], y, lambda m: m)
mu_null = np.exp(X[:, :1] @ beta_null)
deviance_null = 2 * np.sum(y * np.log(y / mu_null) - (y - mu_null))
naive_chisq = deviance_null - deviance                    # the Poisson likelihood ratio
quasi_f = (deviance_null - deviance) / phi_pearson        # referred to F(1, n - p)
print(f"naive chi-squared {naive_chisq:.1f}, p = {stats.chi2.sf(naive_chisq, 1):.2e}")
print(f"quasi F {quasi_f:.2f} on (1, {n - p}) df, p = {stats.f.sf(quasi_f, 1, n - p):.4f}")
# <</tests>>

assert naive_chisq / quasi_f > 30

# ---- negative binomial fits, for the comparison in Section 38.4 -------------
with np.errstate(divide="ignore", invalid="ignore"):
    nb2 = sm.NegativeBinomial(y, X, loglike_method="nb2").fit(disp=0, maxiter=500)
    nb1 = sm.NegativeBinomial(y, X, loglike_method="nb1").fit(disp=0, maxiter=500)
alpha_nb2, alpha_nb1 = nb2.params[-1], nb1.params[-1]
assert alpha_nb2 > 0 and alpha_nb1 > 0
# NB1 has variance (1 + alpha) mu, so its dispersion should be near the Pearson phi
assert abs((1 + alpha_nb1) - phi_pearson) < 0.25 * phi_pearson
# at these fitted means the NB2 variance is dominated by mu^2, so the NB2 slope nearly
# coincides with the quasi fit at V(mu) = mu^2 (Section 38.4)
assert abs(nb2.params[1] - beta_sq[1]) < 0.1
print(f"NB2 slope {nb2.params[1]:.3f} (se {nb2.bse[1]:.3f}), 1/kappa = {alpha_nb2:.3f}")
print(f"NB1 slope {nb1.params[1]:.3f} (se {nb1.bse[1]:.3f}), phi = {1 + alpha_nb1:.3f}")

gen = Generated("ch38", "strikes")
gen.int("n", n)
gen.int("df", n - p)
gen.num("ybar", y.mean(), 2)
gen.num("yvar", y.var(ddof=1), 1)
gen.num("beta0", beta_lin[0], 3)
gen.num("beta1", beta_lin[1], 3)
gen.num("beta1_sq", beta_sq[1], 3)
gen.num("pearson", pearson, 1)
gen.num("deviance", deviance, 1)
gen.num("phi_pearson", phi_pearson, 2)
gen.num("phi_deviance", phi_deviance, 2)
gen.num("score", score_stat, 1)
gen.num("se_model", se_model[1], 3)
gen.num("se_quasi", se_quasi[1], 3)
gen.num("se_sand", se_sand[1], 3)
gen.num("se_corr", se_corr[1], 3)
gen.num("t_model", beta_lin[1] / se_model[1], 2)
gen.num("t_quasi", beta_lin[1] / se_quasi[1], 2)
gen.num("t_sand", beta_lin[1] / se_sand[1], 2)
gen.num("t_corr", beta_lin[1] / se_corr[1], 2)
gen.num("ratio_sand_quasi", se_sand[1] / se_quasi[1], 3)
gen.num("hmax", hat.max(), 3)
gen.num("naive_chisq", naive_chisq, 1)
gen.num("quasi_f", quasi_f, 2)
gen.num("quasi_p", stats.f.sf(quasi_f, 1, n - p), 4)
gen.num("nb2_beta1", nb2.params[1], 3)
gen.num("nb2_se", nb2.bse[1], 3)
gen.num("nb2_alpha", alpha_nb2, 3)
gen.num("nb1_beta1", nb1.params[1], 3)
gen.num("nb1_se", nb1.bse[1], 3)
gen.num("nb1_phi", 1 + alpha_nb1, 2)
gen.num("ymax", ymax, 0)
gen.num("mu_at_ymax", mu_at_ymax, 2)
gen.num("resid_at_ymax", ymax - mu_at_ymax, 1)
gen.write()
