"""Chapter 5, Section 8: where the linear model's assumptions visibly fail.

RAND Health Insurance Experiment (statsmodels.datasets.randhie). The number of physician
visits is a count: it cannot be negative, it has many zeros, and its variance grows with
its mean. A linear model and a Poisson regression with a log link (Chapter 37) are fitted
to the chronic-disease score.
"""
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import statsmodels.api as sm

from regbook import COLORS, Generated, figure_path, use_book_style

# <<counts>>
rand = sm.datasets.randhie.load_pandas().data
y = rand["mdvis"].to_numpy(dtype=float)
x = rand["disea"].to_numpy()
X = sm.add_constant(x)

bins = pd.qcut(x, 10, duplicates="drop")        # ten groups by disease score
groups = pd.DataFrame({"y": y, "x": x}).groupby(bins, observed=True)
table = groups.agg(x=("x", "mean"), mean=("y", "mean"), var=("y", "var"))
print(table.round(2).to_string(index=False))
print(f"share of zeros {np.mean(y == 0):.3f}; variance/mean overall {y.var() / y.mean():.2f}")

linear = sm.OLS(y, X).fit()
poisson = sm.GLM(y, X, family=sm.families.Poisson()).fit()
print("linear model:    ", np.round(linear.params, 4))
print("Poisson, log link:", np.round(poisson.params, 4))
mu = poisson.fittedvalues
print("X^T (y - mu_hat) =", X.T @ (y - mu))
# <</counts>>

ratio = table["var"] / table["mean"]
assert np.all(ratio > 3)                              # far more variable than Poisson
assert np.corrcoef(table["mean"], table["var"])[0, 1] > 0.9
# the Poisson score equations: X^T (y - mu_hat) = 0
assert np.allclose(X.T @ (y - mu), 0, atol=1e-6 * len(y))

gen = Generated("ch05", "panorama_counts", prefix="pan")
gen.num("zeros", 100 * np.mean(y == 0), 1)
gen.num("dispersion", y.var() / y.mean(), 2)
gen.num("ratiomin", ratio.min(), 1)
gen.num("ratiomax", ratio.max(), 1)
gen.num("lin0", linear.params[0], 3)
gen.num("lin1", linear.params[1], 4)
gen.num("pois0", poisson.params[0], 3)
gen.num("pois1", poisson.params[1], 4)
gen.num("poispct", 100 * (np.exp(poisson.params[1]) - 1), 2)
gen.num("maxmdvis", y.max(), 0)
gen.write()

# ---- figure -------------------------------------------------------------------
use_book_style()
fig, axes = plt.subplots(1, 2, figsize=(5.6, 2.5))
ax = axes[0]
ax.plot(table["x"], table["mean"], "o", markersize=3.5, color=COLORS["accent"], label="decile means")
xs = np.linspace(x.min(), np.percentile(x, 99.5), 200)
ax.plot(xs, linear.params[0] + linear.params[1] * xs, color=COLORS["second"], label="linear model")
ax.plot(xs, np.exp(poisson.params[0] + poisson.params[1] * xs), color=COLORS["third"],
        linestyle="--", label="Poisson, log link")
ax.set_xlabel("chronic-disease score")
ax.set_ylabel("mean visits per year")
ax.set_title("(a) models for the mean")
ax.legend(frameon=False, loc="upper left")
ax = axes[1]
ax.plot(table["mean"], table["var"], "o", markersize=3.5, color=COLORS["accent"])
mm = np.linspace(0, table["mean"].max() * 1.05, 2)
ax.plot(mm, mm, color=COLORS["muted"], linestyle=":", label="variance = mean")
ax.set_xlabel("mean visits")
ax.set_ylabel("variance of visits")
ax.set_title("(b) variance against mean")
ax.set_xlim(0, None)
ax.set_ylim(0, None)
ax.legend(frameon=False, loc="upper left")
fig.tight_layout()
fig.savefig(figure_path("ch05", "panorama_counts"))
