"""Chapter 41, Section 4: the EM algorithm.

Part one is the linear model with missing responses, where EM has a closed form:
the fixed point reproduces the complete-case coefficients exactly and divides the
residual sum of squares by n - m, not by n. It is run on the CO2 record of
Section 41.1.

Part two is the harder case, a missing covariate whose missingness depends on the
(observed) response. The E step needs the conditional distribution of the
covariate given the response, the observed-data log-likelihood is recorded at
every sweep, and the linear rate of convergence is compared with the largest
eigenvalue of the matrix of missing information.
"""
import matplotlib.pyplot as plt
import numpy as np
import statsmodels.api as sm

from regbook import COLORS, Generated, figure_path, use_book_style

# ---- part one: missing responses --------------------------------------------
# <<responses>>
co2 = sm.datasets.co2.load_pandas().data["co2"]
t = (co2.index.year + (co2.index.dayofyear - 0.5) / 365.25).to_numpy()
tc = t - 1980.0
X = np.column_stack([np.ones(len(t)), tc, tc**2,
                     np.cos(2 * np.pi * t), np.sin(2 * np.pi * t),
                     np.cos(4 * np.pi * t), np.sin(4 * np.pi * t)])
y = co2.to_numpy()
obs = ~np.isnan(y)
n, m, p = len(y), int((~obs).sum()), X.shape[1]

beta, sigma2 = np.zeros(p), 1.0                       # any starting value will do
for sweep in range(60):
    filled = np.where(obs, y, X @ beta)               # E step: gaps get their fitted values
    beta, *_ = np.linalg.lstsq(X, filled, rcond=None)  # M step for the coefficients
    sse = float(np.sum((filled - X @ beta) ** 2))
    sigma2 = (sse + m * sigma2) / n                   # M step for the variance: add m*sigma^2
print(f"EM after {sweep + 1} sweeps: sigma-hat = {np.sqrt(sigma2):.4f} ppm")

beta_cc, *_ = np.linalg.lstsq(X[obs], y[obs], rcond=None)
sse_o = float(np.sum((y[obs] - X[obs] @ beta_cc) ** 2))
print(f"observed-data MLE:          {np.sqrt(sse_o / (n - m)):.4f} ppm")
print(f"treating the fills as data: {np.sqrt(sse_o / n):.4f} ppm")
# <</responses>>

assert np.allclose(beta, beta_cc, atol=1e-8)
assert abs(sigma2 - sse_o / (n - m)) < 1e-9
naive_ratio = np.sqrt(n / (n - m))

# ---- part two: a missing covariate ------------------------------------------
# <<covariate>>
rng = np.random.default_rng(4104)
N = 400
theta_true = dict(b0=1.0, b1=0.8, s2=1.0, mu=0.5, t2=1.0)
x = theta_true["mu"] + np.sqrt(theta_true["t2"]) * rng.normal(size=N)
y2 = theta_true["b0"] + theta_true["b1"] * x + np.sqrt(theta_true["s2"]) * rng.normal(size=N)
seen = rng.uniform(size=N) > 1 / (1 + np.exp(-(-0.4 + 2.0 * (y2 - y2.mean()))))
print(f"{N - seen.sum()} of {N} covariates missing, more often when y is large")


def estep(th, y, seen):
    """Conditional mean and second moment of the covariate given the response."""
    b0, b1, s2, mu, t2 = th
    vy = b1**2 * t2 + s2                              # marginal variance of y
    cond_mean = mu + (b1 * t2 / vy) * (y - b0 - b1 * mu)
    cond_var = t2 - (b1 * t2) ** 2 / vy
    ex = np.where(seen, x, cond_mean)
    ex2 = np.where(seen, x**2, cond_mean**2 + cond_var)
    return ex, ex2


def mstep(ex, ex2, y):
    mu, n_ = ex.mean(), len(y)
    t2 = ex2.mean() - mu**2
    sxx = ex2.sum() - n_ * mu**2
    sxy = (ex * y).sum() - n_ * mu * y.mean()
    b1 = sxy / sxx
    b0 = y.mean() - b1 * mu
    s2 = ((y**2).sum() - 2 * b0 * y.sum() - 2 * b1 * (ex * y).sum()
          + n_ * b0**2 + 2 * b0 * b1 * ex.sum() + b1**2 * ex2.sum()) / n_
    return np.array([b0, b1, s2, mu, t2])


def loglik(th, y, seen):
    """Observed-data log-likelihood: joint for complete rows, marginal for the rest."""
    b0, b1, s2, mu, t2 = th
    ll = -0.5 * np.sum(np.log(2 * np.pi * s2) + (y[seen] - b0 - b1 * x[seen]) ** 2 / s2)
    ll += -0.5 * np.sum(np.log(2 * np.pi * t2) + (x[seen] - mu) ** 2 / t2)
    vy = b1**2 * t2 + s2
    ll += -0.5 * np.sum(np.log(2 * np.pi * vy) + (y[~seen] - b0 - b1 * mu) ** 2 / vy)
    return ll


th = np.array([0.0, 0.0, 1.0, 0.0, 1.0])
path, iterates = [loglik(th, y2, seen)], [th]
for sweep in range(200):
    ex, ex2 = estep(th, y2, seen)
    th = mstep(ex, ex2, y2)
    path.append(loglik(th, y2, seen))
    iterates.append(th)
print("EM   b0, b1 =", np.round(th[:2], 4))

cc = np.polyfit(x[seen], y2[seen], 1)[::-1]           # complete cases only
print("complete cases:", np.round(cc, 4))
# <</covariate>>

path, iterates = np.array(path), np.array(iterates)
assert np.all(np.diff(path) > -1e-9)                  # monotone, by the theorem
assert abs(th[1] - theta_true["b1"]) < 0.12
assert cc[1] < th[1] - 0.15                           # the complete-case slope is attenuated


def numerical_hessian(f, th, h=1e-4):
    k = len(th)
    H = np.zeros((k, k))
    for i in range(k):
        for j in range(k):
            ei, ej = np.zeros(k), np.zeros(k)
            ei[i] = ej[j] = h
            H[i, j] = (f(th + ei + ej) - f(th + ei - ej)
                       - f(th - ei + ej) + f(th - ei - ej)) / (4 * h * h)
    return H


info_obs = -numerical_hessian(lambda u: loglik(u, y2, seen), th)
ex, ex2 = estep(th, y2, seen)


def q_function(u):
    """E[complete-data log-likelihood at u | observed data, EM solution]."""
    b0, b1, s2, mu, t2 = u
    q = -0.5 * (N * np.log(2 * np.pi * s2)
                + ((y2**2).sum() - 2 * b0 * y2.sum() - 2 * b1 * (ex * y2).sum()
                   + N * b0**2 + 2 * b0 * b1 * ex.sum() + b1**2 * ex2.sum()) / s2)
    q += -0.5 * (N * np.log(2 * np.pi * t2) + (ex2.sum() - 2 * mu * ex.sum() + N * mu**2) / t2)
    return q


info_com = -numerical_hessian(q_function, th)
rate_theory = float(np.max(np.abs(np.linalg.eigvals(
    np.eye(5) - np.linalg.solve(info_com, info_obs)))))
dist = np.linalg.norm(iterates - th, axis=1)          # distance to the fixed point
good = np.flatnonzero(dist > 1e-11)[-40:]             # before it underflows
rate_empirical = float(np.exp(np.polyfit(good, np.log(dist[good]), 1)[0]))
gap = path[-1] - path[:-1]
good_l = np.flatnonzero(gap > 1e-9)[-30:]
rate_loglik = float(np.exp(np.polyfit(good_l, np.log(gap[good_l]), 1)[0]))
print(f"rate: empirical {rate_empirical:.4f}, largest eigenvalue {rate_theory:.4f}")
assert abs(rate_empirical - rate_theory) < 0.01
assert abs(rate_loglik - rate_theory**2) < 0.01       # the likelihood gap halves twice as fast

se_obs = np.sqrt(np.diag(np.linalg.inv(info_obs)))
se_com = np.sqrt(np.diag(np.linalg.inv(info_com)))
boot = []
rng_b = np.random.default_rng(4105)
x_full, x_orig = x.copy(), x
for _ in range(300):
    idx = rng_b.integers(0, N, N)
    x = x_full[idx]                                   # estep/loglik read the module-level x
    yb, sb = y2[idx], seen[idx]
    tb = np.array([0.0, 0.0, 1.0, 0.0, 1.0])
    for _ in range(300):
        e1, e2 = estep(tb, yb, sb)
        tb = mstep(e1, e2, yb)
    boot.append(tb[1])
x = x_orig
se_boot = float(np.std(boot, ddof=1))
print(f"se(b1): observed information {se_obs[1]:.4f}, bootstrap {se_boot:.4f},"
      f" complete-data information {se_com[1]:.4f}")
assert se_com[1] < 0.85 * se_obs[1]                   # the complete-data information is too large
assert abs(se_boot - se_obs[1]) < 0.25 * se_obs[1]

gen = Generated("ch41", "em")
gen.int("n", n)
gen.int("m", m)
gen.num("sigma_em", np.sqrt(sigma2), 4)
gen.num("sigma_naive", np.sqrt(sse_o / n), 4)
gen.num("naive_ratio", naive_ratio, 4)
gen.int("N", N)
gen.int("nmis", int(N - seen.sum()))
gen.num("b1_em", th[1], 3)
gen.num("b0_em", th[0], 3)
gen.num("b1_cc", cc[1], 3)
gen.num("rate_emp", rate_empirical, 3)
gen.num("rate_eig", rate_theory, 3)
gen.num("rate_loglik", rate_loglik, 3)
gen.num("se_obs", se_obs[1], 4)
gen.num("se_boot", se_boot, 4)
gen.num("se_com", se_com[1], 4)
gen.int("sweeps_em", int(np.argmax(dist < 1e-8)))
gen.write()

# ---- figure ------------------------------------------------------------------
use_book_style()
fig, axes = plt.subplots(1, 2, figsize=(5.6, 2.3))
ax = axes[0]
ax.plot(np.arange(3, 61), path[3:61], color=COLORS["accent"], marker="o",
        markersize=1.8, linewidth=0.8)
ax.set_xlabel("sweep (the first three are off the scale)")
ax.set_ylabel("observed-data log-likelihood")
ax.set_title("(a) monotone increase")
ax = axes[1]
k = np.arange(good[-1] + 1)
ax.semilogy(k, dist[: good[-1] + 1], color=COLORS["accent"], label="EM")
ax.semilogy(k, dist[20] * rate_theory ** (k - 20.0), color=COLORS["second"],
            linestyle="--", label=f"rate {rate_theory:.3f}")
ax.set_xlabel("sweep")
ax.set_ylabel(r"distance to $\hat\theta$")
ax.set_title("(b) linear rate")
ax.legend(frameon=False)
fig.tight_layout()
fig.savefig(figure_path("ch41", "em_convergence"))
