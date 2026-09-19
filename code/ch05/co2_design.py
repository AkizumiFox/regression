"""Chapter 5, Section 3: a model that is curved in time but linear in its parameters.

Monthly mean atmospheric CO2 (ppm) at Mauna Loa, March 1958 to December 2001
(Keeling and Whorf; public domain, statsmodels.datasets.co2; weekly values averaged by
month). Model matrix: intercept, a quadratic trend in centred time, and two pairs of
annual harmonics.
"""
import matplotlib.pyplot as plt
import numpy as np
import statsmodels.api as sm

from regbook import COLORS, Generated, figure_path, use_book_style

# <<design>>
co2 = sm.datasets.co2.load_pandas().data["co2"]
monthly = co2.resample("MS").mean().dropna()        # weekly values -> monthly means
y = monthly.to_numpy()
t = (monthly.index.year + (monthly.index.month - 0.5) / 12).to_numpy()
tc = t - 1980.0                                      # years since 1980

X = np.column_stack([
    np.ones_like(tc),                                # intercept
    tc, tc**2,                                       # quadratic trend
    np.cos(2 * np.pi * t), np.sin(2 * np.pi * t),    # annual cycle
    np.cos(4 * np.pi * t), np.sin(4 * np.pi * t),    # half-year harmonic
])
n, p = X.shape
print("n =", n, " p =", p, " rank =", np.linalg.matrix_rank(X))
# <</design>>

# <<fit>>
beta_hat = np.linalg.solve(X.T @ X, X.T @ y)         # normal equations (fine here)
fitted = X @ beta_hat
resid = y - fitted
print("coefficients:", np.round(beta_hat, 4))
print(f"residual standard deviation {np.sqrt(resid @ resid / (n - p)):.3f} ppm")

beta_line = np.linalg.solve(X[:, :2].T @ X[:, :2], X[:, :2].T @ y)
print(f"straight line only: slope {beta_line[1]:.3f} ppm per year")
# <</fit>>

assert np.linalg.matrix_rank(X) == p
lstsq, *_ = np.linalg.lstsq(X, y, rcond=None)
assert np.allclose(beta_hat, lstsq, rtol=1e-8)
assert np.allclose(X.T @ resid, 0, atol=1e-6 * np.abs(X.T @ y).max())
sse_line = np.sum((y - X[:, :2] @ beta_line) ** 2)
sse_full = resid @ resid
assert sse_full < sse_line

amp1 = np.hypot(beta_hat[3], beta_hat[4])
amp2 = np.hypot(beta_hat[5], beta_hat[6])
growth_1960 = beta_hat[1] + 2 * beta_hat[2] * (1960 - 1980)
growth_2000 = beta_hat[1] + 2 * beta_hat[2] * (2000 - 1980)
seasonal = X[:, 3:] @ beta_hat[3:]
peak_to_trough = seasonal.max() - seasonal.min()
R2 = 1 - sse_full / np.sum((y - y.mean()) ** 2)

gen = Generated("ch05", "co2_design", prefix="co2")
gen.int("n", n)
gen.int("p", p)
for j, name in enumerate(["b0", "b1", "b2", "c1", "s1", "c2", "s2"]):
    gen.num(name, beta_hat[j], 4)
gen.num("amp1", amp1, 2)
gen.num("amp2", amp2, 2)
gen.num("g1960", growth_1960, 2)
gen.num("g2000", growth_2000, 2)
gen.num("ptt", peak_to_trough, 1)
gen.num("sdres", np.sqrt(sse_full / (n - p)), 3)
gen.num("sdline", np.sqrt(sse_line / (n - 2)), 3)
gen.num("lineslope", beta_line[1], 3)
gen.num("R2", R2, 4)
r1 = np.corrcoef(resid[:-1], resid[1:])[0, 1]
assert r1 > 0.5
gen.num("r1", r1, 2)
gen.write()

# ---- figure -------------------------------------------------------------------
use_book_style()
fig, axes = plt.subplots(1, 2, figsize=(5.6, 2.5), gridspec_kw={"width_ratios": [1.6, 1]})
ax = axes[0]
ax.plot(t, y, color=COLORS["accent"], linewidth=0.6, label="monthly mean")
ax.plot(t, X[:, :3] @ beta_hat[:3], color=COLORS["second"], linewidth=1.0, label="quadratic trend")
ax.set_xlabel("year")
ax.set_ylabel(r"CO$_2$ (ppm)")
ax.set_title("(a) data and fitted trend")
ax.legend(frameon=False, loc="upper left")
ax = axes[1]
month = np.arange(1, 13)
tt = 1980 + (month - 0.5) / 12
season_curve = (beta_hat[3] * np.cos(2 * np.pi * tt) + beta_hat[4] * np.sin(2 * np.pi * tt)
                + beta_hat[5] * np.cos(4 * np.pi * tt) + beta_hat[6] * np.sin(4 * np.pi * tt))
detr = y - X[:, :3] @ beta_hat[:3]
ax.scatter(monthly.index.month + np.random.default_rng(0).uniform(-0.2, 0.2, n), detr, s=3,
           color=COLORS["accent"], alpha=0.35, linewidths=0)
ax.plot(month, season_curve, color=COLORS["second"], marker="o", markersize=2.5)
ax.axhline(0, color=COLORS["grid"], linewidth=0.6, zorder=0)
ax.set_xticks([1, 4, 7, 10])
ax.set_xticklabels(["Jan", "Apr", "Jul", "Oct"])
ax.set_ylabel("detrended (ppm)")
ax.set_title("(b) fitted seasonal cycle")
fig.tight_layout()
fig.savefig(figure_path("ch05", "co2_design"))
