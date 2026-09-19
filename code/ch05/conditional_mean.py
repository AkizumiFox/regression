"""Chapter 5, Section 1: what a regression model claims.

(i) Figure: two data-generating processes with the same regression function
    E(Y | x) = 1 + 0.5 x. In (a) the errors have constant variance and are symmetric;
    in (b) they are skewed and their spread grows with x. Both satisfy the mean
    assumption; only (a) satisfies the second-moment assumptions of the linear model.
(ii) Simulation: the least squares slope under a fixed design and under a random
    design. It is unbiased in both; its variance is sigma^2 / S_xx given the design,
    and E(sigma^2 / S_xx) when the design is redrawn with every sample.
"""
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats

from regbook import COLORS, Generated, figure_path, use_book_style

gen = Generated("ch05", "conditional_mean", prefix="cm")

# <<fixed-random>>
rng = np.random.default_rng(2005)
n, beta0, beta1, sigma = 15, 1.0, 0.5, 1.0
reps = 50_000


def slope(x, y):
    """Least squares slope of y on x (with an intercept)."""
    xc = x - x.mean()
    return xc @ (y - y.mean()) / (xc @ xc)


x_fixed = rng.uniform(0, 10, size=n)            # drawn once, then held fixed
Sxx_fixed = np.sum((x_fixed - x_fixed.mean()) ** 2)
b_fixed = np.empty(reps)
b_random = np.empty(reps)
Sxx_random = np.empty(reps)
for r in range(reps):
    e = sigma * rng.standard_normal(n)
    b_fixed[r] = slope(x_fixed, beta0 + beta1 * x_fixed + e)
    x = rng.uniform(0, 10, size=n)              # a new design for every sample
    Sxx_random[r] = np.sum((x - x.mean()) ** 2)
    b_random[r] = slope(x, beta0 + beta1 * x + e)

print(f"fixed design:  mean {b_fixed.mean():.4f}  var {b_fixed.var():.5f}"
      f"  sigma^2/Sxx = {sigma**2 / Sxx_fixed:.5f}")
print(f"random design: mean {b_random.mean():.4f}  var {b_random.var():.5f}"
      f"  E(sigma^2/Sxx) = {np.mean(sigma**2 / Sxx_random):.5f}")
# <</fixed-random>>

se_fixed = b_fixed.std() / np.sqrt(reps)
se_random = b_random.std() / np.sqrt(reps)
assert abs(b_fixed.mean() - beta1) < 4 * se_fixed
assert abs(b_random.mean() - beta1) < 4 * se_random
assert abs(b_fixed.var() / (sigma**2 / Sxx_fixed) - 1) < 0.02
e_inv = np.mean(sigma**2 / Sxx_random)
assert abs(b_random.var() / e_inv - 1) < 0.02
# Jensen: E(1/Sxx) >= 1/E(Sxx), and E(Sxx) = (n-1) * 100/12 for U(0,10)
ESxx = (n - 1) * 100 / 12
assert e_inv > sigma**2 / ESxx

gen.int("n", n)
gen.text("reps", f"{reps:,}")
gen.num("Sxxfixed", Sxx_fixed, 2)
gen.num("meanfixed", b_fixed.mean(), 4)
gen.num("meanrandom", b_random.mean(), 4)
gen.num("varfixed", b_fixed.var(), 5)
gen.num("varrandom", b_random.var(), 5)
gen.num("theoryfixed", sigma**2 / Sxx_fixed, 5)
gen.num("theoryrandom", e_inv, 5)
gen.num("ESxx", ESxx, 2)
gen.num("jensen", sigma**2 / ESxx, 5)
gen.write()

# ---- figure: two processes with the same regression function -----------------
use_book_style()
rng_fig = np.random.default_rng(7)
m = lambda x: 1 + 0.5 * x
fig, axes = plt.subplots(1, 2, figsize=(5.6, 2.5), sharey=True)
xs_pts = rng_fig.uniform(0, 10, size=150)
# (a) normal, constant variance
ya = m(xs_pts) + rng_fig.normal(0, 0.8, size=150)
# (b) centred gamma errors, shape 2, standard deviation 0.15 + 0.25 x
k = 2.0
sd_b = 0.15 + 0.25 * xs_pts
yb = m(xs_pts) + sd_b * (rng_fig.gamma(k, 1 / np.sqrt(k), size=150) - np.sqrt(k))
for ax, yy, title, kind in ((axes[0], ya, "(a) linear model assumptions hold", "a"),
                            (axes[1], yb, "(b) same mean, other distribution", "b")):
    ax.scatter(xs_pts, yy, s=6, color=COLORS["accent"], alpha=0.55, linewidths=0)
    grid = np.linspace(0, 10, 2)
    ax.plot(grid, m(grid), color=COLORS["second"])
    for x0 in (2.0, 5.0, 8.0):
        if kind == "a":
            t = np.linspace(-2.6, 2.6, 200)
            dens = stats.norm.pdf(t, scale=0.8)
        else:
            s0 = 0.15 + 0.25 * x0
            t = np.linspace(-1.4 * s0, 3.2 * s0, 200)
            u = t / s0 * 1.0 + np.sqrt(k)          # back to the gamma(k, 1/sqrt(k)) scale
            dens = stats.gamma.pdf(u, k, scale=1 / np.sqrt(k)) / s0
        ax.plot(x0 - 1.4 * dens / dens.max(), m(x0) + t, color=COLORS["third"], linewidth=0.9)
        ax.plot([x0, x0], [m(x0) + t.min(), m(x0) + t.max()], color=COLORS["grid"],
                linewidth=0.6, zorder=0)
    ax.set_xlabel("x")
    ax.set_title(title)
axes[0].set_ylabel("y")
fig.tight_layout()
fig.savefig(figure_path("ch05", "conditional_mean"))
