"""Chapter 6, Section 11: least squares estimates the population projection,
even when the conditional mean is not linear."""
import matplotlib.pyplot as plt
import numpy as np

from regbook import COLORS, Generated, figure_path, use_book_style

rng = np.random.default_rng(611)

# <<population>>
# X ~ Uniform(0, 2), E(Y | X) = exp(X), Var(Y | X) = 0.25.
# Population projection of Y onto span{1, X}: slope Cov(X, Y)/Var(X), intercept E Y - slope * E X.
EX, VarX = 1.0, 4.0 / 12.0
EY = (np.exp(2) - 1) / 2                             # E exp(X)
EXY = (np.exp(2) + 1) / 2                            # E X exp(X)
slope_pop = (EXY - EX * EY) / VarX
intercept_pop = EY - slope_pop * EX

for n in [20, 200, 2000, 20000]:
    x = rng.uniform(0, 2, n)
    y = np.exp(x) + rng.normal(0, 0.5, n)
    b, *_ = np.linalg.lstsq(np.column_stack([np.ones(n), x]), y, rcond=None)
    print(f"n = {n:6d}   intercept {b[0]:.3f}   slope {b[1]:.3f}")
print(f"population    intercept {intercept_pop:.3f}   slope {slope_pop:.3f}")
# <</population>>

# Monte Carlo check of the population values
x = rng.uniform(0, 2, 2_000_000)
yy = np.exp(x) + rng.normal(0, 0.5, x.size)
bb, *_ = np.linalg.lstsq(np.column_stack([np.ones(x.size), x]), yy, rcond=None)
assert np.allclose(bb, [intercept_pop, slope_pop], atol=0.01)

gen = Generated("ch06", "best_linear_predictor", prefix="blp")
gen.num("slope", slope_pop, 3)
gen.num("intercept", intercept_pop, 3)
gen.write()

use_book_style()
fig, ax = plt.subplots(figsize=(3.8, 2.6))
xs = np.sort(rng.uniform(0, 2, 400))
ys = np.exp(xs) + rng.normal(0, 0.5, xs.size)
ax.scatter(xs, ys, s=5, color=COLORS["muted"], alpha=0.5, linewidths=0)
g = np.linspace(0, 2, 100)
ax.plot(g, np.exp(g), color=COLORS["ink"], label=r"$\mathrm{E}(Y\mid X)=e^X$")
ax.plot(g, intercept_pop + slope_pop * g, color=COLORS["second"], label="population projection")
b, *_ = np.linalg.lstsq(np.column_stack([np.ones(xs.size), xs]), ys, rcond=None)
ax.plot(g, b[0] + b[1] * g, color=COLORS["accent"], ls="--", label="least squares, $n=400$")
ax.set_xlabel("$x$"); ax.set_ylabel("$y$")
ax.legend(frameon=False, loc="upper left")
fig.savefig(figure_path("ch06", "best_linear_predictor"))
