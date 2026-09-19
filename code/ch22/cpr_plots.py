"""Chapter 22, Section 3: component-plus-residual plots, and when another regressor distorts them.

Simulated data (fixed seed): E(Y) = 3 log(x1) + x2 with x1 uniform on [1, 10].
(a) x2 linearly related to x1 (correlation near 0.96): the plot still shows the logarithm.
(b) x2 a nonlinear function of x1 (square root, plus a little noise): part of the curvature is
    absorbed by x2. The expected partial residuals are g - M2(g - b x1) exactly.
(c) Mallows' augmented partial residuals (x1^2 added to the fit) recover most of the curvature.
"""
import matplotlib.pyplot as plt
import numpy as np

from regbook import COLORS, Generated, figure_path, use_book_style

# <<cpr>>
rng = np.random.default_rng(2204)
n = 120
x1 = rng.uniform(1, 10, n)
g = 3 * np.log(x1)                                        # the true shape in x1
standardize = lambda v: (v - v.mean()) / v.std()
x2_lin = standardize(x1) + rng.normal(0, 0.3, n)          # linear in x1, plus noise
x2_sqrt = standardize(np.sqrt(x1)) + rng.normal(0, 0.1, n)  # nonlinear in x1, plus noise


def lstsq(X, v):
    coef, *_ = np.linalg.lstsq(X, v, rcond=None)
    return coef


def partial_residuals(y, x2, augmented=False):
    """e + b1 x1 (+ b11 x1^2 if augmented), from the fit of y on 1, x1 (, x1^2), x2."""
    X = np.column_stack([np.ones(len(y)), x1] + ([x1**2] if augmented else []) + [x2])
    coef = lstsq(X, y)
    e = y - X @ coef
    return e + coef[1] * x1 + (coef[2] * x1**2 if augmented else 0)


def curvature_retained(v):
    """Coefficient of the nonlinear part of g in the nonlinear part of v (1 = shape kept)."""
    Z = np.column_stack([np.ones(n), x1])
    nonlin = lambda u: u - Z @ lstsq(Z, u)
    return nonlin(v) @ nonlin(g) / (nonlin(g) @ nonlin(g))


noise = rng.normal(0, 0.6, n)
results = {}
for name, x2, aug in [("linear", x2_lin, False), ("square root", x2_sqrt, False),
                      ("square root, augmented", x2_sqrt, True)]:
    mean = g + x2                                         # E(Y); partial residuals are linear in y
    results[name] = (partial_residuals(mean + noise, x2, aug), partial_residuals(mean, x2, aug))
    print(f"{name:23s} corr(x1, x2) {np.corrcoef(x1, x2)[0, 1]:.2f}   curvature retained "
          f"{curvature_retained(results[name][1]):.2f}")
# <</cpr>>

# the exact formula E(pr) = g - M2 (g - b x1), with M2 projecting onto C([1, x2])
slopes = []
for x2 in (x2_lin, x2_sqrt):
    X2 = np.column_stack([np.ones(n), x2])
    M2 = lambda v: X2 @ lstsq(X2, v)
    x1t = x1 - M2(x1)
    b = x1t @ g / (x1t @ x1t)
    assert np.allclose(partial_residuals(g + x2, x2), g - M2(g - b * x1))
    slopes.append(b)
assert slopes[0] > 0 > slopes[1]                         # with x2 close to sqrt(x1) the expected slope is negative
# a least squares line in the plot, with or without intercept, has slope b1 and leaves the fit's residuals
X = np.column_stack([np.ones(n), x1, x2_sqrt])
y = g + x2_sqrt + noise
coef = lstsq(X, y)
pr = partial_residuals(y, x2_sqrt)
assert np.isclose(x1 @ pr / (x1 @ x1), coef[1])
Z = np.column_stack([np.ones(n), x1])
c = lstsq(Z, pr)
assert np.allclose(c, [0, coef[1]], atol=1e-10) and np.allclose(pr - Z @ c, y - X @ coef)

ret = {k: curvature_retained(v[1]) for k, v in results.items()}
assert ret["linear"] > 0.95 and ret["square root"] < 0.6 and ret["square root, augmented"] > 0.8

gen = Generated("ch22", "cpr_plots", prefix="cpr")
gen.int("n", n)
gen.num("ret_lin", ret["linear"], 2)
gen.num("ret_sqrt", ret["square root"], 2)
gen.num("ret_aug", ret["square root, augmented"], 2)
gen.num("corr_lin", np.corrcoef(x1, x2_lin)[0, 1], 2)
gen.num("corr_sqrt", np.corrcoef(x1, x2_sqrt)[0, 1], 2)
gen.write()

# ---- figure ----------------------------------------------------------------------------------
use_book_style()
fig, axes = plt.subplots(1, 3, figsize=(5.8, 2.2), sharey=True)
order = np.argsort(x1)
titles = {"linear": "(a) $x_2$ linear in $x_1$", "square root": r"(b) $x_2\approx\sqrt{x_1}$",
          "square root, augmented": "(c) (b), augmented"}
for ax, (name, (pr, pr_mean)) in zip(axes, results.items()):
    ax.scatter(x1, pr - pr.mean(), s=5, color=COLORS["accent"], alpha=0.6, linewidths=0)
    # expected partial residuals: markers, since eq. 22.3.1 depends on the simulated x2 and is not a curve in x1
    ax.scatter(x1, pr_mean - pr_mean.mean(), s=4, marker="s", color=COLORS["second"], linewidths=0)
    ax.plot(x1[order], (g - g.mean())[order], color=COLORS["muted"],
            linewidth=1.0, linestyle="--")
    ax.set_title(titles[name])
    ax.set_xlabel("$x_1$")
axes[0].set_ylabel("partial residual (centred)")
fig.tight_layout()
fig.savefig(figure_path("ch22", "cpr_plots"))
