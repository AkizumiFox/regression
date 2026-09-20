"""Chapter 37, Section 1: the Poisson log-linear model fitted to physician visits.

Data: the RAND Health Insurance Experiment sample shipped with statsmodels
(statsmodels.datasets.randhie, public domain), 20190 person-years. Response: the
number of outpatient visits to a physician in the year.

Computes the maximum likelihood fit by Fisher scoring, checks the likelihood
equations and the fitted-marginal property, reports rate ratios and the deviance,
and draws the log-likelihood surface of a two-parameter submodel with the scoring
path on it.
"""
import matplotlib.pyplot as plt
import numpy as np
import statsmodels.api as sm

from regbook import COLORS, Generated, figure_path, use_book_style

# <<fit>>
import numpy as np
import statsmodels.api as sm

data = sm.datasets.randhie.load_pandas().data        # 20190 person-years, public domain
y = data["mdvis"].to_numpy(float)                    # outpatient physician visits in the year
names = ["intercept", "lncoins", "idp", "physlm", "disea", "hlthg", "hlthf", "hlthp"]
X = np.column_stack([np.ones(len(y))] + [data[v].to_numpy(float) for v in names[1:]])


def poisson_fit(X, y, offset=0.0, tol=1e-11, maxit=50):
    """Fisher scoring for log mu = X beta + offset; returns beta and its iterates."""
    beta = np.zeros(X.shape[1])
    beta[0] = np.log(y.mean())
    path = [beta]
    for _ in range(maxit):
        mu = np.exp(X @ beta + offset)               # current fitted means
        z = X @ beta + (y - mu) / mu                 # working response
        XW = X * mu[:, None]                         # weights w_i = mu_i
        step = np.linalg.solve(X.T @ XW, XW.T @ z)
        converged = np.max(np.abs(step - beta)) < tol
        beta = step
        path.append(beta)
        if converged:
            break
    return beta, np.array(path)


beta, path = poisson_fit(X, y)
mu = np.exp(X @ beta)
se = np.sqrt(np.diag(np.linalg.inv(X.T @ (X * mu[:, None]))))

print("largest entry of X^T (y - mu):", np.abs(X.T @ (y - mu)).max())
for name, b, s in zip(names, beta, se):
    print(f"{name:10s} {b:8.4f} ({s:.4f})   rate ratio {np.exp(b):.4f}")
# <</fit>>

assert np.abs(X.T @ (y - mu)).max() < 1e-6
assert len(path) <= 8

# the fitted marginal totals of each health-status group match the observed totals
for column in ["hlthg", "hlthf", "hlthp"]:
    g = data[column].to_numpy(float) > 0.5
    assert np.isclose(mu[g].sum(), y[g].sum())
assert np.isclose(mu.sum(), y.sum())

# <<deviance>>
dev = 2 * np.sum(np.where(y > 0, y * np.log(np.where(y > 0, y, 1) / mu), 0.0) - (y - mu))
pearson = np.sum((y - mu) ** 2 / mu)
df = len(y) - X.shape[1]
print(f"deviance {dev:.1f} on {df} df   (ratio {dev / df:.3f})")
print(f"Pearson  {pearson:.1f} on {df} df   (ratio {pearson / df:.3f})")
print("fitted total", mu.sum(), " observed total", y.sum())
# <</deviance>>

assert np.isclose(dev, 2 * np.sum(np.where(y > 0, y * np.log(np.where(y > 0, y, 1) / mu), 0.0)))

# log-likelihood of the fit (dropping the constant -sum log y!)
log_factorial = np.array([np.sum(np.log(np.arange(1, int(v) + 1))) for v in y])
loglik = np.sum(y * np.log(mu) - mu - log_factorial)

# a Wald interval for the coinsurance rate ratio
j = names.index("lncoins")
lo, hi = np.exp(beta[j] - 1.96 * se[j]), np.exp(beta[j] + 1.96 * se[j])
assert lo < np.exp(beta[j]) < hi

# the full range of lncoins, from free care to the 95 percent coinsurance plan
span = data["lncoins"].max() - data["lncoins"].min()
ratio_span = np.exp(beta[j] * span)

# disea is a continuous chronic-disease score, not a count, so quote the effect of moving
# across its interquartile range as well as the effect of one unit
q1, q3 = np.percentile(data["disea"].to_numpy(float), [25, 75])
rr_disea_iqr = np.exp(beta[names.index("disea")] * (q3 - q1))
assert data["disea"].max() > 50 and np.mean(data["disea"] % 1 != 0) > 0.9

gen = Generated("ch37", "poisson_visits")
gen.int("n", len(y))
gen.int("p", X.shape[1])
gen.int("steps", len(path) - 1)
gen.num("mean", y.mean(), 3)
gen.num("var_over_mean", y.var(ddof=1) / y.mean(), 3)
for name, b, s in zip(names, beta, se):
    gen.num("b_" + name, b, 4)
    gen.num("se_" + name, s, 4)
    gen.num("rr_" + name, np.exp(b), 4)
gen.num("rr_lncoins_lo", lo, 4)
gen.num("rr_lncoins_hi", hi, 4)
gen.num("disea_q1", q1, 2)
gen.num("disea_q3", q3, 2)
gen.num("rr_disea_iqr", rr_disea_iqr, 3)
gen.num("span", span, 3)
gen.num("rr_span", ratio_span, 3)
gen.num("deviance", dev, 1)
gen.num("pearson", pearson, 1)
gen.int("df", df)
gen.num("dev_ratio", dev / df, 3)
gen.num("pearson_ratio", pearson / df, 3)
gen.num("loglik", loglik, 2)
gen.write()

# ---- the log-likelihood surface of a two-parameter submodel ------------------
x1 = data["disea"].to_numpy(float)
X2 = np.column_stack([np.ones(len(y)), x1])
beta2, path2 = poisson_fit(X2, y)
log_fact_sum = log_factorial.sum()


def loglik2(b0, b1):
    eta = b0 + b1 * x1
    return np.sum(y * eta - np.exp(eta)) - log_fact_sum


grid0 = np.linspace(0.42, 1.12, 120)
grid1 = np.linspace(-0.004, 0.060, 120)
Z = np.array([[loglik2(b0, b1) for b0 in grid0] for b1 in grid1])
assert Z.max() <= loglik2(*beta2) + 1e-6          # the fit maximizes the surface

use_book_style()
fig, ax = plt.subplots(figsize=(4.4, 3.1))
levels = loglik2(*beta2) - np.array([6000, 2400, 900, 340, 120, 40, 10])
ax.contour(grid0, grid1, Z, levels=levels, colors=COLORS["grid"], linewidths=0.7,
           linestyles="solid")
ax.plot(path2[:, 0], path2[:, 1], "-o", color=COLORS["accent"], markersize=3.4,
        linewidth=1.0, label="Fisher scoring")
ax.plot([beta2[0]], [beta2[1]], "*", color=COLORS["second"], markersize=11,
        label="maximum", zorder=5)
ax.set_xlabel(r"intercept $\beta_0$")
ax.set_ylabel(r"slope $\beta_1$ (chronic-disease score)")
ax.legend(loc="lower left", frameon=False)
fig.tight_layout()
fig.savefig(figure_path("ch37", "scoring_path"))

gen2 = Generated("ch37", "scoring_path")
gen2.int("steps", len(path2) - 1)
gen2.num("b0", beta2[0], 4)
gen2.num("b1", beta2[1], 4)
gen2.write()
