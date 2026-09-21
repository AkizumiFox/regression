"""Chapter 38, Section 3: model-based against sandwich covariance for a quasi-Poisson fit.

The mean model log mu = b0 + b1 x is always correct. The working variance function is
always V(mu) = phi mu. The data are generated twice: with a truly quadratic variance
(so the working variance function is wrong in shape), and with a truly linear one (so it
is right, and the model-based covariance is valid). Coverage of nominal 95% intervals
for b1 is recorded for four standard errors.
"""
import matplotlib.pyplot as plt
import numpy as np

from regbook import COLORS, Generated, figure_path, use_book_style

BETA = np.array([1.5, 1.0])
HALF = 2.0            # the regressor runs over [-HALF, HALF]


# <<coverage>>
def quasi_poisson(X, y, steps=25):
    """IRLS for the quasi-score equations with a log link and V(mu) = mu."""
    beta = np.zeros(X.shape[1])
    beta[0] = np.log(max(y.mean(), 0.5))
    for _ in range(steps):
        mu = np.exp(X @ beta)
        z = X @ beta + (y - mu) / mu
        beta = np.linalg.solve((X.T * mu) @ X, (X.T * mu) @ z)
    return beta, np.exp(X @ beta)


def slope_variances(X, y):
    """The four estimates of Var(b1): naive, phi-scaled, sandwich, leverage-corrected."""
    beta, mu = quasi_poisson(X, y)
    n, p = X.shape
    bread = np.linalg.inv((X.T * mu) @ X)
    root = np.sqrt(mu)[:, None] * X
    hat = np.sum(root * (root @ bread), axis=1)
    r2 = (y - mu) ** 2
    phi = np.sum(r2 / mu) / (n - p)
    meat = lambda w: bread @ ((X * w[:, None]).T @ X) @ bread
    return beta[1], np.array([bread[1, 1], phi * bread[1, 1],
                              meat(r2)[1, 1], meat(r2 / (1 - hat) ** 2)[1, 1]])


def coverage(n, quadratic, reps, rng, kappa=0.5, phi=9.0):
    """Coverage of nominal 95% normal intervals for b1 over `reps` simulated data sets."""
    x = np.linspace(-HALF, HALF, n)
    X = np.column_stack([np.ones(n), x])
    mu = np.exp(X @ BETA)
    shape = np.full(n, kappa) if quadratic else mu / (phi - 1.0)
    hit = np.zeros(4)
    for _ in range(reps):
        y = rng.negative_binomial(shape, shape / (shape + mu)).astype(float)
        b1, var = slope_variances(X, y)
        hit += np.abs(b1 - BETA[1]) <= 1.96 * np.sqrt(var)
    return hit / reps


rng = np.random.default_rng(3821)
names = ["Poisson", "phi-scaled", "sandwich", "corrected"]
for n in (25, 100, 400):
    quad = coverage(n, True, 300, rng)
    lin = coverage(n, False, 300, rng)
    print(f"n = {n:4d}  quadratic variance " + "  ".join(
        f"{k} {v:.3f}" for k, v in zip(names, quad)))
    print(f"          linear variance    " + "  ".join(
        f"{k} {v:.3f}" for k, v in zip(names, lin)))
# <</coverage>>

sizes = np.array([20, 40, 80, 160, 320, 640, 1280])
REPS = 3000
rng = np.random.default_rng(38211)
quad = np.array([coverage(n, True, REPS, rng) for n in sizes])
lin = np.array([coverage(n, False, REPS, rng) for n in sizes])

# the sandwich improves with n in both panels; the naive Poisson interval never does
assert quad[-1, 2] > quad[0, 2] and lin[-1, 2] > lin[0, 2]
assert quad[:, 0].max() < 0.55 and lin[:, 0].max() < 0.55
# with a correctly specified variance function the phi-scaled interval is already
# near nominal at the smallest sample size, and the sandwich is not
assert lin[0, 1] > lin[0, 2] + 0.02
assert abs(lin[:, 1] - 0.95).max() < 0.05
# with a quadratic true variance the model-based interval does not reach nominal coverage
assert quad[-1, 1] < quad[-1, 2] - 0.02

gen = Generated("ch38", "sandwich_coverage")
gen.int("reps", REPS)
for j, key in enumerate(["naive", "phi", "sand", "corr"]):
    for i, n in [(0, 20), (len(sizes) - 1, 1280)]:
        gen.num(f"quad_{key}_{n}", quad[i, j], 3)
        gen.num(f"lin_{key}_{n}", lin[i, j], 3)
gen.num("quad_sand_80", quad[2, 2], 3)
gen.num("quad_corr_80", quad[2, 3], 3)
gen.write()

# ---- the figure ------------------------------------------------------------
use_book_style()
fig, axes = plt.subplots(1, 2, figsize=(5.8, 2.5), sharey=True)
styles = [(COLORS["muted"], ":", "Poisson"), (COLORS["accent"], "-", r"$\hat\phi$-scaled"),
          (COLORS["second"], "--", "sandwich"), (COLORS["third"], "-.", "corrected")]
for ax, table, title in [(axes[0], quad, "(a) true variance quadratic"),
                         (axes[1], lin, "(b) true variance linear")]:
    for j, (colour, dash, label) in enumerate(styles):
        ax.plot(sizes, table[:, j], color=colour, linestyle=dash, marker="o",
                markersize=2.5, label=label)
    ax.axhline(0.95, color=COLORS["grid"], linewidth=0.8, zorder=0)
    ax.set_xscale("log")
    ax.set_xticks(sizes)
    ax.set_xticklabels([str(s) for s in sizes])
    ax.set_xlabel("$n$")
    ax.set_title(title)
axes[0].set_ylabel("coverage")
axes[0].set_ylim(0.25, 1.0)
axes[1].legend(frameon=False, loc="center right", ncol=2, handlelength=1.6,
               columnspacing=1.0, borderaxespad=0.2)
fig.tight_layout()
fig.savefig(figure_path("ch38", "sandwich_coverage"))
