"""Chapter 12, Section 6: credible intervals, a credible ellipse and the posterior predictive
distribution in the conjugate normal-inverse-gamma model, and the frequentist coverage of a
credible interval.

Stack loss data (statsmodels.datasets.stackloss, public domain) with the prior of Chapter 7.
"""
import matplotlib.pyplot as plt
import numpy as np
import statsmodels.api as sm
from scipy import stats

from regbook import COLORS, Generated, figure_path, use_book_style

# <<posterior>>
data = sm.datasets.stackloss.load_pandas().data
y = data["STACKLOSS"].to_numpy()
X = np.column_stack([np.ones(len(y)), data[["AIRFLOW", "WATERTEMP", "ACIDCONC"]]])
n, p = X.shape
XtX, Xty = X.T @ X, X.T @ y


def nig_posterior(m0, P0, a0, b0, y=y):
    """Conjugate update; P0 is the prior precision V0^{-1} (per unit sigma^2)."""
    Vn = np.linalg.inv(P0 + XtX)
    mn = Vn @ (P0 @ m0 + X.T @ y)
    an = a0 + n / 2
    bn = b0 + (y @ y + m0 @ P0 @ m0 - mn @ (P0 + XtX) @ mn) / 2
    return mn, Vn, an, bn


m0 = np.array([0.0, 1.0, 1.0, 0.0])                     # the prior of Chapter 7
P0 = np.linalg.inv(np.diag([100.0**2, 0.1**2, 0.2**2, 0.1**2]))
a0, b0 = 3.0, 18.0
mn, Vn, an, bn = nig_posterior(m0, P0, a0, b0)
# <</posterior>>

# <<predictive>>
x0 = np.array([1.0, 65.0, 23.0, 87.0])                  # a new operating condition
qb = stats.t.ppf(0.975, 2 * an)
centre = x0 @ mn
half_mean = qb * np.sqrt(bn / an * (x0 @ Vn @ x0))     # credible interval for x0^T beta
half_new = qb * np.sqrt(bn / an * (1 + x0 @ Vn @ x0))  # posterior predictive interval for Y0

beta_hat = np.linalg.solve(XtX, Xty)                   # frequentist intervals, for comparison
s2 = np.sum((y - X @ beta_hat) ** 2) / (n - p)
qf = stats.t.ppf(0.975, n - p)
h0 = x0 @ np.linalg.solve(XtX, x0)
print(f"Bayes: mean response {centre:.2f} +- {half_mean:.2f}; new day {centre:.2f} +- {half_new:.2f}")
print(f"least squares: mean response {x0 @ beta_hat:.2f} +- {qf * np.sqrt(s2 * h0):.2f}; "
      f"new day +- {qf * np.sqrt(s2 * (1 + h0)):.2f}")
# <</predictive>>

# ---- checks: the flat prior reproduces the frequentist intervals exactly ------------------------
flat_m, flat_V, flat_a, flat_b = beta_hat, np.linalg.inv(XtX), (n - p) / 2, s2 * (n - p) / 2
assert np.isclose(flat_b / flat_a, s2) and 2 * flat_a == n - p
assert np.isclose(stats.t.ppf(0.975, 2 * flat_a) * np.sqrt(flat_b / flat_a * (1 + x0 @ flat_V @ x0)),
                  qf * np.sqrt(s2 * (1 + h0)))
# the proper prior with P0 -> 0 and a0 -> -p/2, b0 -> 0 approaches the flat-prior posterior
mt, Vt, at, bt = nig_posterior(m0, 1e-10 * np.eye(p), -p / 2, 0.0)
assert np.allclose(mt, beta_hat, atol=1e-5) and np.isclose(bt, flat_b, rtol=1e-6) and at == flat_a

# ---- Monte Carlo from the posterior: marginal and predictive quantiles, ellipse content ----------
rng = np.random.default_rng(1209)
draws = 400_000
sig2 = 1 / rng.gamma(an, 1 / bn, size=draws)           # sigma^2 | y ~ IG(a_n, b_n)
beta_draws = mn + np.sqrt(sig2)[:, None] * rng.multivariate_normal(np.zeros(p), Vn, size=draws)
y0_draws = beta_draws @ x0 + np.sqrt(sig2) * rng.normal(size=draws)
mc_new = np.quantile(y0_draws, [0.025, 0.975])
mc_mean = np.quantile(beta_draws @ x0, [0.025, 0.975])
assert np.allclose(mc_new, [centre - half_new, centre + half_new], atol=0.08)
assert np.allclose(mc_mean, [centre - half_mean, centre + half_mean], atol=0.05)
L = np.eye(p)[:, [1, 2]]                                # air flow and water temperature
Wb = L.T @ Vn @ L
dev = beta_draws @ L - mn @ L
Q = np.einsum("ij,jk,ik->i", dev, np.linalg.inv(Wb), dev) / (2 * bn / an)
content = np.mean(Q <= stats.f.ppf(0.95, 2, 2 * an))
assert abs(content - 0.95) < 0.003

# ---- frequentist coverage of the 95% credible interval for the air-flow coefficient ------------
def credible_coverage(beta_true, sigma_true, P0=P0, reps=4000, seed=0):
    rs = np.random.default_rng(seed)
    Y = X @ beta_true + sigma_true * rs.normal(size=(reps, n))
    Vn_ = np.linalg.inv(P0 + XtX)
    Mn = (Y @ X + P0 @ m0) @ Vn_                        # rows: posterior means (Vn symmetric)
    quad = np.einsum("ij,jk,ik->i", Mn, P0 + XtX, Mn)
    Bn = b0 + (np.sum(Y**2, axis=1) + m0 @ P0 @ m0 - quad) / 2
    half = qb * np.sqrt(Bn / an * Vn_[1, 1])
    return np.mean(np.abs(Mn[:, 1] - beta_true[1]) <= half)


sigma_true = np.sqrt(s2)
P0_strong = P0.copy()
P0_strong[1:, 1:] *= 25                                  # slope prior sd divided by 5
grid = np.linspace(0.3, 1.7, 57)
cov_curve, cov_strong = [], []
for k, b_air in enumerate(grid):
    bt_ = beta_hat.copy()
    bt_[1] = b_air
    cov_curve.append(credible_coverage(bt_, sigma_true, seed=100 + k))
    cov_strong.append(credible_coverage(bt_, sigma_true, P0=P0_strong, seed=200 + k))
cov_curve, cov_strong = np.array(cov_curve), np.array(cov_strong)
cov_at = {v: credible_coverage(np.r_[beta_hat[0], v, beta_hat[2:]], sigma_true, reps=20_000, seed=7)
          for v in (0.7, 1.0, 1.5)}
cov_at_strong = {v: credible_coverage(np.r_[beta_hat[0], v, beta_hat[2:]], sigma_true, P0=P0_strong,
                                      reps=20_000, seed=8) for v in (0.7, 1.0, 1.5)}

# averaged over the prior, the coverage is exactly 1 - alpha
reps_prior = 40_000
rs = np.random.default_rng(1210)
s2_prior = 1 / rs.gamma(a0, 1 / b0, size=reps_prior)
beta_prior = m0 + np.sqrt(s2_prior)[:, None] * rs.normal(size=(reps_prior, p)) * np.sqrt(np.diag(np.linalg.inv(P0)))
Yp = np.einsum("ij,kj->ki", X, beta_prior) + np.sqrt(s2_prior)[:, None] * rs.normal(size=(reps_prior, n))
Mn = (Yp @ X + P0 @ m0) @ Vn
Bn = b0 + (np.sum(Yp**2, axis=1) + m0 @ P0 @ m0 - np.einsum("ij,jk,ik->i", Mn, P0 + XtX, Mn)) / 2
cov_prior = np.mean(np.abs(Mn[:, 1] - beta_prior[:, 1]) <= qb * np.sqrt(Bn / an * Vn[1, 1]))
print(f"coverage of the credible interval: at the prior mean {cov_at[1.0]:.3f}, at 0.7 {cov_at[0.7]:.3f}, "
      f"at 1.5 {cov_at[1.5]:.3f}; averaged over the prior {cov_prior:.4f}")
assert abs(cov_prior - 0.95) < 4 * np.sqrt(0.95 * 0.05 / reps_prior)
assert cov_at[1.0] > 0.95 and cov_curve.min() > 0.85
assert cov_at_strong[1.0] > 0.99 and cov_at_strong[0.7] < 0.05 and cov_at_strong[1.5] < 0.01

gen = Generated("ch12", "bayes_intervals", prefix="bay")
gen.num("an", an, 1)
gen.int("df", 2 * an)
gen.num("qb", qb, 3)
gen.num("qf", qf, 3)
gen.num("centre", centre, 2)
gen.num("hmean", half_mean, 2)
gen.num("hnew", half_new, 2)
gen.num("fit", x0 @ beta_hat, 2)
gen.num("fhmean", qf * np.sqrt(s2 * h0), 2)
gen.num("fhnew", qf * np.sqrt(s2 * (1 + h0)), 2)
gen.num("mcnew:lo", mc_new[0], 2)
gen.num("mcnew:hi", mc_new[1], 2)
gen.num("content", content, 4)
gen.int("draws", draws)
gen.num("cov:prior", cov_prior, 4)
for v, c in cov_at.items():
    gen.num(f"cov:{v}", c, 3)
for v, c in cov_at_strong.items():
    gen.num(f"covs:{v}", c, 3)
gen.num("covmin", cov_curve.min(), 3)
gen.num("covmax", cov_curve.max(), 3)
gen.int("repsprior", reps_prior)
gen.write()

# ---- figure ----------------------------------------------------------------------------------
use_book_style()
fig, ax = plt.subplots(figsize=(4.4, 3.5))
ax.plot(grid, cov_curve, color=COLORS["thread"], label="prior of Chapter 7")
ax.plot(grid, cov_strong, color=COLORS["second"], linestyle="--", label="slope prior sd divided by 5")
ax.axhline(0.95, color=COLORS["accent"], linewidth=1.0, label="flat prior (the $t$ interval)")
ax.axvline(1.0, color=COLORS["muted"], linewidth=0.6, linestyle=":", label="prior mean of the coefficient")
ax.set_xlabel("true air-flow coefficient")
ax.set_ylabel("coverage of 95% interval")
ax.set_ylim(0, 1.02)
ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.22), ncol=1, frameon=False)
fig.tight_layout()
fig.savefig(figure_path("ch12", "credible_coverage"))
