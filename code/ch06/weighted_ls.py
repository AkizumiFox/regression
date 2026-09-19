"""Chapter 6, Section 9: weighted least squares as projection in a weighted inner product.

Engel's 1857 data on household income and food expenditure (Belgian working-class
households), public domain, statsmodels.datasets.engel. Spread grows with income.
"""
import matplotlib.pyplot as plt
import numpy as np
import statsmodels.api as sm

from regbook import COLORS, Generated, figure_path, use_book_style

df = sm.datasets.engel.load_pandas().data
x, y = df["income"].to_numpy(), df["foodexp"].to_numpy()
X = np.column_stack([np.ones(len(x)), x])

# <<wls>>
w = 1.0 / x**2                                   # working model: sd proportional to income
sw = np.sqrt(w)
beta_ols, *_ = np.linalg.lstsq(X, y, rcond=None)
beta_wls, *_ = np.linalg.lstsq(X * sw[:, None], y * sw, rcond=None)   # OLS on (W^1/2 X, W^1/2 y)

# leverages in the weighted geometry: squared row norms of an orthonormal basis of C(W^1/2 X)
Qw, _ = np.linalg.qr(X * sw[:, None])
h_wls = np.sum(Qw**2, axis=1)
Q, _ = np.linalg.qr(X)
h_ols = np.sum(Q**2, axis=1)
print("OLS:", beta_ols, " WLS:", beta_wls)
print("max leverage OLS %.3f, WLS %.3f" % (h_ols.max(), h_wls.max()))
# <</wls>>

W = np.diag(w)
P = X @ np.linalg.solve(X.T @ W @ X, X.T @ W)
assert np.allclose(P @ P, P) and np.allclose(W @ P, (W @ P).T)
assert np.allclose(np.diag(P), h_wls) and np.isclose(h_wls.sum(), 2)
assert np.allclose(X @ beta_wls, P @ y)

gen = Generated("ch06", "weighted_ls", prefix="wls")
gen.int("n", len(y))
gen.num("ols_slope", beta_ols[1], 4)
gen.num("wls_slope", beta_wls[1], 4)
gen.num("ols_int", beta_ols[0], 1)
gen.num("wls_int", beta_wls[0], 1)
gen.num("hmax_ols", h_ols.max(), 3)
gen.num("hmax_wls", h_wls.max(), 3)
gen.num("xmax", x.max(), 0)
gen.num("x_at_hmax_wls", x[np.argmax(h_wls)], 0)
gen.write()

use_book_style()
fig, axes = plt.subplots(1, 2, figsize=(5.6, 2.4))
ax = axes[0]
ax.scatter(x, y, s=6, color=COLORS["muted"], alpha=0.7, linewidths=0)
xs = np.linspace(x.min(), x.max(), 2)
ax.plot(xs, beta_ols[0] + beta_ols[1] * xs, color=COLORS["accent"], label="OLS")
ax.plot(xs, beta_wls[0] + beta_wls[1] * xs, color=COLORS["second"], ls="--", label="WLS")
ax.set_xlabel("household income"); ax.set_ylabel("food expenditure")
ax.legend(frameon=False); ax.set_title("(a) fits")
ax = axes[1]
ax.scatter(x, h_ols, s=6, color=COLORS["accent"], linewidths=0, label="OLS")
ax.scatter(x, h_wls, s=6, color=COLORS["second"], linewidths=0, label="WLS")
ax.set_xlabel("household income"); ax.set_ylabel("leverage")
ax.legend(frameon=False); ax.set_title("(b) leverage in each geometry")
fig.tight_layout()
fig.savefig(figure_path("ch06", "weighted_ls"))
