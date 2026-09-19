"""Chapter 7, Section 3: maximum likelihood in the normal linear model, and the profile likelihood.

Brownlee's stack loss data (public domain, shipped with statsmodels): 21 days of a plant
oxidizing ammonia; response = stack loss, regressors = air flow, cooling water temperature,
acid concentration. The MLE of beta is least squares; the MLE of sigma^2 divides by n.
The profile log-likelihood of one coefficient is an increasing function of its t statistic.
"""
import matplotlib.pyplot as plt
import numpy as np
from scipy import optimize, stats

from regbook import COLORS, Generated, figure_path, use_book_style

# <<fit>>
import numpy as np
import statsmodels.api as sm
from scipy import optimize, stats

data = sm.datasets.stackloss.load_pandas().data
y = data["STACKLOSS"].to_numpy()
X = np.column_stack([np.ones(len(y)), data[["AIRFLOW", "WATERTEMP", "ACIDCONC"]]])
n, p = X.shape

beta_hat, *_ = np.linalg.lstsq(X, y, rcond=None)
sse = np.sum((y - X @ beta_hat) ** 2)
sigma2_mle, s2 = sse / n, sse / (n - p)
loglik_max = -n / 2 * (np.log(2 * np.pi * sigma2_mle) + 1)
print("beta_hat   ", np.round(beta_hat, 4))
print(f"SSE = {sse:.3f}   MLE of sigma^2 = {sigma2_mle:.3f}   s^2 = {s2:.3f}")
print(f"maximized log-likelihood = {loglik_max:.3f}")
# <</fit>>


# <<numerical>>
def negloglik(theta):
    """Minus the normal log-likelihood; theta = (beta, log sigma^2)."""
    b, log_s2 = theta[:p], theta[p]
    r = y - X @ b
    return n / 2 * (np.log(2 * np.pi) + log_s2) + r @ r / (2 * np.exp(log_s2))


start = np.r_[np.zeros(p), np.log(np.var(y))]
opt = optimize.minimize(negloglik, start, method="BFGS", options={"gtol": 1e-8})
print("numerical MLE of beta  ", np.round(opt.x[:p], 4))
print("numerical MLE sigma^2  ", round(float(np.exp(opt.x[p])), 3))
# <</numerical>>

assert np.allclose(opt.x[:p], beta_hat, atol=1e-3)
assert np.isclose(np.exp(opt.x[p]), sigma2_mle, rtol=1e-4)
assert np.isclose(-opt.fun, loglik_max, atol=1e-6)


# <<profile>>
def profile_loglik(b1):
    """Profile log-likelihood of the air-flow coefficient: maximize over everything else."""
    others = np.delete(X, 1, axis=1)
    z = y - b1 * X[:, 1]                                 # move the fixed term to the left side
    coef, *_ = np.linalg.lstsq(others, z, rcond=None)
    sse_b = np.sum((z - others @ coef) ** 2)             # SSE(b1): smallest SSE with beta_1 = b1
    return -n / 2 * (np.log(2 * np.pi * sse_b / n) + 1)


se1 = np.sqrt(s2 * np.linalg.inv(X.T @ X)[1, 1])
grid = np.linspace(beta_hat[1] - 4 * se1, beta_hat[1] + 4 * se1, 201)
lp = np.array([profile_loglik(b) for b in grid])
t = (beta_hat[1] - grid) / se1
print("max |2(l_max - l_p) - n log(1 + t^2/(n-p))| =",
      np.abs(2 * (loglik_max - lp) - n * np.log1p(t ** 2 / (n - p))).max())
# <</profile>>

assert np.allclose(2 * (loglik_max - lp), n * np.log1p(t ** 2 / (n - p)), atol=1e-9)
assert np.isclose(lp.max(), loglik_max, atol=1e-3)

# likelihood-ratio interval {b : 2(l_max - l_p(b)) <= chi^2_{1, 0.95}} versus the t interval
chi = stats.chi2.ppf(0.95, 1)
t_lr = np.sqrt((n - p) * (np.exp(chi / n) - 1))          # the |t| at which the cutoff is reached
tq = stats.t.ppf(0.975, n - p)
lr_int = beta_hat[1] + np.array([-1, 1]) * t_lr * se1
t_int = beta_hat[1] + np.array([-1, 1]) * tq * se1
for b in lr_int:
    assert np.isclose(2 * (loglik_max - profile_loglik(b)), chi, atol=1e-8)
assert t_lr < tq                                         # the LR interval is a little too short

# mean squared error of c * SSE as an estimator of sigma^2, relative to sigma^4
k = n - p


def rel_mse(divisor):
    return 2 * k / divisor ** 2 + (k / divisor - 1) ** 2


mse = {"unbiased": rel_mse(k), "mle": rel_mse(n), "best": rel_mse(k + 2)}
assert mse["best"] < min(mse["unbiased"], mse["mle"])
assert np.isclose(mse["best"], 2 / (k + 2))

gen = Generated("ch07", "mle_profile", prefix="mle")
gen.int("n", n)
gen.int("p", p)
for j in range(p):
    gen.num(f"b{j}", beta_hat[j], 4)
gen.num("sse", sse, 3)
gen.num("s2mle", sigma2_mle, 3)
gen.num("s2", s2, 3)
gen.num("loglik", loglik_max, 3)
gen.num("se1", se1, 4)
gen.num("tlr", t_lr, 3)
gen.num("tq", tq, 3)
gen.num("shortfall_pct", 100 * (1 - t_lr / tq), 0)
gen.num("chi", chi, 3)
gen.num("lr_lo", lr_int[0], 3)
gen.num("lr_hi", lr_int[1], 3)
gen.num("t_lo", t_int[0], 3)
gen.num("t_hi", t_int[1], 3)
gen.int("k", k)
for key, v in mse.items():
    gen.num(f"mse:{key}", v, 4)
gen.write()

use_book_style()
fig, ax = plt.subplots(figsize=(4.2, 2.6))
ax.plot(grid, lp - loglik_max, color=COLORS["accent"], label="profile log-likelihood")
ax.axhline(-chi / 2, color=COLORS["muted"], lw=0.7, ls="--")
ax.text(grid[0], -chi / 2 + 0.15, r"cutoff $-\chi^2_{1,0.95}/2$", fontsize=7, color=COLORS["muted"])
ax.plot(lr_int, [-chi / 2] * 2, "|", color=COLORS["accent"], ms=9, mew=1.4)
ax.plot(t_int, [-2.9, -2.9], "-|", color=COLORS["second"], lw=1.2, ms=7, label=r"95% $t$ interval")
ax.axvline(beta_hat[1], color=COLORS["grid"], lw=0.6, zorder=0)
ax.set_xlabel(r"air-flow coefficient $\beta_1$")
ax.set_ylabel(r"$\ell_p(\beta_1)-\ell(\hat{\boldsymbol{\theta}})$")
ax.set_ylim(-4.2, 0.3)
ax.legend(frameon=False, loc="lower center", fontsize=7)
fig.savefig(figure_path("ch07", "profile_likelihood"))
