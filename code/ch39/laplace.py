"""Chapter 39, Section 4: the Laplace approximation to a posterior and to a marginal likelihood.

A two-parameter logistic regression with a normal prior. The marginal likelihood is
computed to machine accuracy by adaptive Gauss-Hermite quadrature on a tensor grid
centred at the posterior mode, and compared with the Laplace value as the sample
size grows: the relative error falls like 1/n.
"""
import matplotlib.pyplot as plt
import numpy as np
from scipy import special, stats

from regbook import COLORS, Generated, figure_path, use_book_style

TAU = 5.0                                       # prior standard deviation, N(0, TAU^2 I)
BETA0 = np.array([-0.4, 1.2])                   # the truth used to simulate


def simulate(n, seed):
    rng = np.random.default_rng(seed)
    X = np.column_stack([np.ones(n), rng.uniform(-2, 2, n)])
    y = rng.binomial(1, special.expit(X @ BETA0)).astype(float)
    return X, y


# <<laplace>>
def log_post(beta, X, y):
    """Log likelihood plus log prior, up to a constant free of beta."""
    eta = X @ beta
    return np.sum(y * eta - np.logaddexp(0.0, eta)) - np.sum(beta**2) / (2 * TAU**2)


def mode_and_hessian(X, y):
    """Posterior mode and the negative Hessian there, by Newton's method."""
    beta = np.zeros(X.shape[1])
    for _ in range(100):
        pi = special.expit(X @ beta)
        H = X.T @ ((pi * (1 - pi))[:, None] * X) + np.eye(len(beta)) / TAU**2
        g = X.T @ (y - pi) - beta / TAU**2
        beta = beta + np.linalg.solve(H, g)
    return beta, H


def log_marginal_laplace(X, y):
    """log of the Laplace value of the integral of exp(log_post)."""
    beta, H = mode_and_hessian(X, y)
    p = len(beta)
    return (log_post(beta, X, y) + p / 2 * np.log(2 * np.pi)
            - 0.5 * np.linalg.slogdet(H)[1] - p / 2 * np.log(2 * np.pi * TAU**2))
# <</laplace>>


def log_marginal_quadrature(X, y, nodes=60):
    """The same integral by adaptive Gauss-Hermite quadrature: the reference value."""
    beta, H = mode_and_hessian(X, y)
    p = len(beta)
    L = np.linalg.cholesky(np.linalg.inv(H))
    z, w = np.polynomial.hermite_e.hermegauss(nodes)            # weight exp(-z^2/2)
    Z = np.stack(np.meshgrid(z, z, indexing="ij"), -1).reshape(-1, p)
    logw = np.log(w)[:, None] + np.log(w)[None, :]
    pts = beta + Z @ L.T
    vals = np.array([log_post(b, X, y) for b in pts]) + (Z**2).sum(1) / 2 + logw.ravel()
    return (special.logsumexp(vals) + np.linalg.slogdet(L)[1]
            - p / 2 * np.log(2 * np.pi * TAU**2))


sizes = [25, 50, 100, 200, 400, 800, 1600]
rel_err = []
for k, n in enumerate(sizes):
    X, y = simulate(n, 3910 + k)
    a, b = log_marginal_laplace(X, y), log_marginal_quadrature(X, y)
    rel_err.append(abs(np.expm1(a - b)))
rel_err = np.array(rel_err)
slope = np.polyfit(np.log(sizes), np.log(rel_err), 1)[0]
print("n, relative error of the Laplace marginal likelihood:")
for n, e in zip(sizes, rel_err):
    print(f"  {n:5d}  {e:.3e}   n * error = {n * e:.3f}")
print(f"log-log slope {slope:.3f}")
products = np.array(sizes, float) * rel_err
assert -1.3 < slope < -0.7
assert rel_err[0] < 0.2
# two-sided: the product n x error is bounded away from zero and from infinity, which
# is what an O(1/n) rate with a nonzero constant means
assert 1.0 < products.min() and products.max() < 4.0

# the BIC comparison, on the same sequence
bic_gap = []
for k, n in enumerate(sizes):
    X, y = simulate(n, 3910 + k)
    beta, H = mode_and_hessian(X, y)
    pi = special.expit(X @ beta)
    loglik = np.sum(y * np.log(pi) + (1 - y) * np.log1p(-pi))
    bic = -2 * loglik + X.shape[1] * np.log(n)
    bic_gap.append(bic + 2 * log_marginal_quadrature(X, y))
print("BIC minus -2 log m(y):", np.round(bic_gap, 3))
assert np.ptp(bic_gap) < 4.0                 # an O(1) gap, not one that grows with n

# the posterior and its normal approximation at the smallest sample size
X25, y25 = simulate(25, 3910)
beta25, H25 = mode_and_hessian(X25, y25)
V25 = np.linalg.inv(H25)
grid = np.linspace(beta25[1] - 4 * np.sqrt(V25[1, 1]), beta25[1] + 5 * np.sqrt(V25[1, 1]), 220)
prof = np.array([
    np.trapezoid(np.exp([log_post(np.array([b0, b1]), X25, y25)
                         for b0 in np.linspace(beta25[0] - 6 * np.sqrt(V25[0, 0]),
                                               beta25[0] + 6 * np.sqrt(V25[0, 0]), 220)]),
                 np.linspace(beta25[0] - 6 * np.sqrt(V25[0, 0]),
                             beta25[0] + 6 * np.sqrt(V25[0, 0]), 220))
    for b1 in grid])
prof = prof / np.trapezoid(prof, grid)
approx = stats.norm.pdf(grid, beta25[1], np.sqrt(V25[1, 1]))
tv = 0.5 * np.trapezoid(np.abs(prof - approx), grid)
mean_exact = np.trapezoid(grid * prof, grid)
print(f"n = 25: posterior mean of the slope {mean_exact:.4f}, mode {beta25[1]:.4f}, "
      f"total variation distance to the normal approximation {tv:.4f}")
assert 0.1 < tv < 0.3

gen = Generated("ch39", "laplace")
gen.int("nsmall", 25)
gen.int("nbig", sizes[-1])
gen.num("tau", TAU, 0)
gen.num("err25", rel_err[0], 5)
gen.num("errbig", rel_err[-1], 6)
gen.num("slope", slope, 2)
gen.num("prodmin", float(products.min()), 3)
gen.num("prodmax", float(products.max()), 3)
gen.num("bicgapmin", float(np.min(bic_gap)), 2)
gen.num("bicgapmax", float(np.max(bic_gap)), 2)
gen.num("mode25", beta25[1], 4)
gen.num("mean25", mean_exact, 4)
gen.num("sd25", float(np.sqrt(V25[1, 1])), 4)
gen.num("tv25", tv, 4)
gen.write()

# ---- (a) exact posterior vs normal approximation, (b) the 1/n rate ----------
use_book_style()
fig, axes = plt.subplots(1, 2, figsize=(5.6, 2.4))
ax = axes[0]
ax.plot(grid, prof, color=COLORS["accent"], label="marginal posterior")
ax.plot(grid, approx, color=COLORS["second"], linestyle="--", label="Laplace normal")
ax.axvline(beta25[1], color=COLORS["muted"], linewidth=0.6)
ax.set_xlabel(r"slope $\beta_1$")
ax.set_ylabel("density")
ax.set_title("(a) 25 observations")
ax.legend(frameon=False, loc="upper left")

ax = axes[1]
ax.loglog(sizes, rel_err, "o-", color=COLORS["accent"], markersize=3,
          label="Laplace relative error")
ref = rel_err[0] * sizes[0] / np.array(sizes, float)
ax.loglog(sizes, ref, color=COLORS["muted"], linestyle=":", label=r"slope $-1$")
ax.set_xlabel("sample size")
ax.set_ylabel("relative error")
ax.set_title("(b) error of the marginal likelihood")
ax.legend(frameon=False)
fig.tight_layout()
fig.savefig(figure_path("ch39", "laplace_error"))
