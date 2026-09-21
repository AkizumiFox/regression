"""Chapter 38, Section 4: what ignoring overdispersion does to a test.

Null data with a constant mean and variance phi * mu are fitted by a Poisson
log-linear model with an intercept and a slope. The drop in deviance is referred to
chi-squared(1) (the naive test), to F(1, n - 2) after division by the Pearson
dispersion (the quasi-F test), and a Wald test uses the sandwich covariance.
"""
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats

from regbook import COLORS, Generated, figure_path, use_book_style

N, PHI, MU0 = 100, 4.0, 6.0
REPS = 5000


# <<size>>
def naive_size(phi, q):
    """Size of a nominal 5% chi-squared test when the statistic is inflated by phi."""
    return stats.chi2.sf(stats.chi2.ppf(0.95, q) / phi, q)


print("phi   q = 1    q = 3    q = 6")
for phi in (1.5, 2.0, 3.0, 5.0, 10.0):
    print(f"{phi:4.1f} " + "  ".join(f"{naive_size(phi, q):7.3f}" for q in (1, 3, 6)))
# <</size>>

assert naive_size(1.0, 1) == 0.05 or abs(naive_size(1.0, 1) - 0.05) < 1e-12
assert naive_size(3.0, 6) > naive_size(3.0, 1)   # the inflation bites harder for large q


def poisson_fit(X, y, steps=20):
    beta = np.zeros(X.shape[1])
    beta[0] = np.log(max(y.mean(), 0.5))
    for _ in range(steps):
        mu = np.exp(X @ beta)
        z = X @ beta + (y - mu) / mu
        beta = np.linalg.solve((X.T * mu) @ X, (X.T * mu) @ z)
    return beta, np.exp(X @ beta)


def poisson_deviance(y, mu):
    return 2 * np.sum(np.where(y > 0, y * np.log(np.maximum(y, 1e-300) / mu), 0.0) - (y - mu))


# <<simulation>>
def null_statistics(n, phi, mu0, reps, rng):
    """Three statistics for the slope, on data whose true slope is zero."""
    x = np.linspace(-1.0, 1.0, n)
    X = np.column_stack([np.ones(n), x])
    shape = np.full(n, mu0 / (phi - 1.0))            # negative binomial with variance phi*mu
    out = np.empty((reps, 3))
    for r in range(reps):
        y = rng.negative_binomial(shape, shape / (shape + mu0)).astype(float)
        beta, mu = poisson_fit(X, y)
        beta0, mu_null = poisson_fit(X[:, :1], y)
        drop = poisson_deviance(y, mu_null) - poisson_deviance(y, mu)
        phi_hat = np.sum((y - mu) ** 2 / mu) / (n - 2)
        bread = np.linalg.inv((X.T * mu) @ X)
        sand = bread @ ((X * ((y - mu) ** 2)[:, None]).T @ X) @ bread
        out[r] = [drop, drop / phi_hat, beta[1] ** 2 / sand[1, 1]]
    return out


rng = np.random.default_rng(3841)
stat = null_statistics(N, PHI, MU0, 400, rng)
rates = [np.mean(stat[:, 0] > stats.chi2.ppf(0.95, 1)),
         np.mean(stat[:, 1] > stats.f.ppf(0.95, 1, N - 2)),
         np.mean(stat[:, 2] > stats.chi2.ppf(0.95, 1))]
print(f"rejection rates at the 5% level: naive {rates[0]:.3f},"
      f" quasi-F {rates[1]:.3f}, sandwich Wald {rates[2]:.3f}")
# <</simulation>>

rng = np.random.default_rng(38411)
stat = null_statistics(N, PHI, MU0, REPS, rng)
naive, quasi_f, wald = stat[:, 0], stat[:, 1], stat[:, 2]
rate_naive = np.mean(naive > stats.chi2.ppf(0.95, 1))
rate_quasi = np.mean(quasi_f > stats.f.ppf(0.95, 1, N - 2))
rate_wald = np.mean(wald > stats.chi2.ppf(0.95, 1))

assert rate_naive > 0.25                      # the naive test is far too liberal
assert abs(rate_quasi - 0.05) < 0.02          # the quasi-F test is about right
assert abs(naive.mean() / quasi_f.mean() - PHI) < 1.0
print(f"mean of the naive statistic {naive.mean():.2f} (chi-squared(1) has mean 1)")
print(f"{REPS} replicates: naive {rate_naive:.3f}, quasi-F {rate_quasi:.3f},"
      f" sandwich {rate_wald:.3f}")

gen = Generated("ch38", "overdispersion_tests")
gen.int("n", N)
gen.int("reps", REPS)
gen.num("phi", PHI, 1)
gen.num("mu0", MU0, 1)
gen.num("rate_naive", rate_naive, 3)
gen.num("rate_quasi", rate_quasi, 3)
gen.num("rate_wald", rate_wald, 3)
gen.num("mean_naive", naive.mean(), 2)
gen.num("mean_quasi", quasi_f.mean(), 2)
for phi in (1.5, 2.0, 3.0, 4.0, 5.0, 10.0):
    tag = str(phi).replace(".", "")
    gen.num(f"size1_{tag}", naive_size(phi, 1), 3)
    gen.num(f"size6_{tag}", naive_size(phi, 6), 3)
gen.write()

# ---- the figure ------------------------------------------------------------
use_book_style()
fig, axes = plt.subplots(1, 2, figsize=(5.8, 2.4))
grid = np.linspace(0.02, 14, 400)

bins = np.linspace(0.0, 14.0, 43)
weights = np.full(REPS, 1.0 / REPS / (bins[1] - bins[0]))

ax = axes[0]
ax.hist(naive, bins=bins, weights=weights, color=COLORS["accent"],
        alpha=0.45, linewidth=0)
ax.plot(grid, stats.chi2.pdf(grid, 1), color=COLORS["ink"])
ax.axvline(stats.chi2.ppf(0.95, 1), color=COLORS["second"], linestyle="--", linewidth=0.9)
ax.set_ylim(0, 1.0)
ax.set_xlabel("drop in deviance")
ax.set_ylabel("density")
ax.set_title(r"(a) naive, against $\chi^2(1)$")

ax = axes[1]
ax.hist(quasi_f, bins=bins, weights=weights, color=COLORS["third"],
        alpha=0.45, linewidth=0)
ax.plot(grid, stats.f.pdf(grid, 1, N - 2), color=COLORS["ink"])
ax.axvline(stats.f.ppf(0.95, 1, N - 2), color=COLORS["second"], linestyle="--", linewidth=0.9)
ax.set_ylim(0, 1.0)
ax.set_xlabel("drop in deviance $\\div\\ \\hat\\phi$")
ax.set_title(r"(b) quasi-$F$, against $F(1,98)$")

fig.tight_layout()
fig.savefig(figure_path("ch38", "overdispersion_tests"))
