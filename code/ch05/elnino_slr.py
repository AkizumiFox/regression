"""Chapter 5, Section 2: simple linear regression by the formulas S_xy / S_xx.

Sea surface temperature (degrees C) in the Nino 1+2 region off Peru and Ecuador,
1950-2010 (NOAA, public domain, statsmodels.datasets.elnino). Response: the December
mean; regressor: the August mean of the same year, four months earlier.
"""
import matplotlib.pyplot as plt
import numpy as np
import statsmodels.api as sm

from regbook import COLORS, Generated, figure_path, use_book_style

# <<fit>>
sst = sm.datasets.elnino.load_pandas().data
x = sst["AUG"].to_numpy()          # August temperature
y = sst["DEC"].to_numpy()          # December temperature, same year
n = len(y)

xbar, ybar = x.mean(), y.mean()
Sxx = np.sum((x - xbar) ** 2)
Sxy = np.sum((x - xbar) * (y - ybar))
Syy = np.sum((y - ybar) ** 2)
b1 = Sxy / Sxx                     # slope
b0 = ybar - b1 * xbar              # intercept
fitted = b0 + b1 * x
resid = y - fitted
print(f"n = {n}, slope = {b1:.4f}, intercept = {b0:.4f}")
# <</fit>>

# <<checks>>
print("sum of residuals      ", resid.sum())
print("sum of x * residuals  ", (x * resid).sum())
SSE = resid @ resid
R2 = 1 - SSE / Syy
r = Sxy / np.sqrt(Sxx * Syy)
print(f"R^2 = {R2:.4f}, r^2 = {r**2:.4f}")
print("numpy.polyfit agrees: ", np.polyfit(x, y, 1))
# <</checks>>

assert abs(resid.sum()) < 1e-9 and abs((x * resid).sum()) < 1e-8
assert np.allclose(np.polyfit(x, y, 1), [b1, b0])
assert np.isclose(R2, r**2)
assert np.isclose(SSE, Syy - Sxy**2 / Sxx)
# the fitted line passes through (xbar, ybar)
assert np.isclose(b0 + b1 * xbar, ybar)
# slope = r * s_y / s_x
assert np.isclose(b1, r * y.std() / x.std())
# least squares really is the minimum: perturbations increase the SSE
rng = np.random.default_rng(1)
for _ in range(100):
    d0, d1 = rng.normal(scale=0.1, size=2)
    assert np.sum((y - (b0 + d0) - (b1 + d1) * x) ** 2) > SSE

# through the origin: slope sum(xy)/sum(x^2)
b_origin = (x @ y) / (x @ x)
sse_origin = np.sum((y - b_origin * x) ** 2)
assert sse_origin > SSE

years = sst["YEAR"].to_numpy().astype(int)
i97 = int(np.argmax(x))
# without the largest August value
mask = np.arange(n) != i97
b1_drop = np.polyfit(x[mask], y[mask], 1)[0]

gen = Generated("ch05", "elnino_slr", prefix="sst")
gen.int("n", n)
gen.int("first", years[0])
gen.int("last", years[-1])
gen.num("xbar", xbar, 3)
gen.num("ybar", ybar, 3)
gen.num("Sxx", Sxx, 3)
gen.num("Sxy", Sxy, 3)
gen.num("Syy", Syy, 3)
gen.num("b1", b1, 4)
gen.num("b0", b0, 4)
gen.num("SSE", SSE, 3)
gen.num("R2", R2, 4)
gen.num("r", r, 4)
gen.num("sy", y.std(ddof=1), 3)
gen.num("sx", x.std(ddof=1), 3)
gen.num("s", np.sqrt(SSE / (n - 2)), 4)
gen.int("maxyear", years[i97])
gen.num("maxx", x[i97], 2)
gen.num("maxy", y[i97], 2)
gen.num("maxresid", resid[i97], 3)
gen.num("b1drop", b1_drop, 4)
gen.num("borigin", b_origin, 4)
gen.num("pred22", b0 + b1 * 22.0, 2)
gen.write()

# ---- figure -------------------------------------------------------------------
use_book_style()
fig, axes = plt.subplots(1, 2, figsize=(5.6, 2.5))
ax = axes[0]
ax.scatter(x, y, s=10, color=COLORS["accent"], alpha=0.8, linewidths=0)
xs = np.linspace(x.min() - 0.2, x.max() + 0.2, 2)
ax.plot(xs, b0 + b1 * xs, color=COLORS["second"])
ax.plot([xbar], [ybar], marker="+", color=COLORS["ink"], markersize=8)
ax.annotate(str(years[i97]), (x[i97], y[i97]), textcoords="offset points", xytext=(-26, -3),
            fontsize=7, color=COLORS["muted"])
ax.set_xlabel("August temperature (°C)")
ax.set_ylabel("December temperature (°C)")
ax.set_title("(a) data and least squares line")
ax = axes[1]
ax.scatter(x, resid, s=10, color=COLORS["accent"], alpha=0.8, linewidths=0)
ax.axhline(0, color=COLORS["grid"], linewidth=0.8, zorder=0)
ax.set_xlabel("August temperature (°C)")
ax.set_ylabel("residual (°C)")
ax.set_title("(b) residuals")
fig.tight_layout()
fig.savefig(figure_path("ch05", "elnino_slr"))
