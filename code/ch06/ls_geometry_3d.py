"""Chapter 6, Section 1: least squares with n = 3 observations, drawn in R^3."""
import numpy as np
import matplotlib.pyplot as plt

from regbook import COLORS, Generated, figure_path, use_book_style

# <<setup>>
y = np.array([1.0, 4.0, 2.0])
X = np.column_stack([np.ones(3), [0.0, 1.0, 3.0]])   # intercept and one regressor

beta_hat, *_ = np.linalg.lstsq(X, y, rcond=None)
y_hat = X @ beta_hat              # the point of C(X) nearest to y
e_hat = y - y_hat                 # the residual vector

print("fitted  ", y_hat)
print("residual", e_hat)
print("X^T e   ", X.T @ e_hat)    # zero: the residual is orthogonal to C(X)
# <</setup>>

assert np.allclose(X.T @ e_hat, 0, atol=1e-12)
assert np.isclose(y @ y, y_hat @ y_hat + e_hat @ e_hat, atol=1e-12)
# any other point of the plane is farther away
rng = np.random.default_rng(1)
for b in rng.normal(size=(1000, 2)):
    assert np.sum((y - X @ (beta_hat + b)) ** 2) > e_hat @ e_hat

gen = Generated("ch06", "ls_geometry_3d", prefix="geo")
for i in range(3):
    gen.num(f"yhat{i+1}", y_hat[i], 3)
    gen.num(f"ehat{i+1}", e_hat[i], 3)
gen.num("b0", beta_hat[0], 3)
gen.num("b1", beta_hat[1], 3)
gen.num("yy", y @ y, 3)
gen.num("yhyh", y_hat @ y_hat, 3)
gen.num("ee", e_hat @ e_hat, 3)
gen.write()

# ---- figure ------------------------------------------------------------------
use_book_style()
fig = plt.figure(figsize=(4.2, 3.6))
ax = fig.add_subplot(projection="3d")
s, t = np.meshgrid(np.linspace(-0.4, 2.6, 2), np.linspace(-0.3, 1.25, 2))
P = s[..., None] * X[:, 0] + t[..., None] * X[:, 1]
ax.plot_surface(P[..., 0], P[..., 1], P[..., 2], color=COLORS["accent"], alpha=0.12, linewidth=0)
o = np.zeros(3)


def arrow(a, b, **kw):
    ax.quiver(*a, *(b - a), arrow_length_ratio=0.08, **kw)


arrow(o, y, color=COLORS["ink"], linewidth=1.3)
arrow(o, y_hat, color=COLORS["accent"], linewidth=1.3)
arrow(y_hat, y, color=COLORS["second"], linewidth=1.3)
arrow(o, X[:, 0], color=COLORS["muted"], linewidth=0.8)
arrow(o, X[:, 1], color=COLORS["muted"], linewidth=0.8)
# right-angle marker at y_hat
u = -y_hat / np.linalg.norm(y_hat) * 0.28
v = e_hat / np.linalg.norm(e_hat) * 0.28
ax.plot(*np.array([y_hat + u, y_hat + u + v, y_hat + v]).T, color=COLORS["ink"], linewidth=0.6)
ax.text(*(y * 1.06), r"$\mathbf{y}$", color=COLORS["ink"])
ax.text(*(y_hat + np.array([0.15, 0.1, -0.35])), r"$\hat{\mathbf{y}}=\mathbf{M}\mathbf{y}$", color=COLORS["accent"])
ax.text(*((y + y_hat) / 2 + np.array([0.1, 0.05, 0])), r"$\hat{\boldsymbol{\varepsilon}}$", color=COLORS["second"])
ax.text(*(X[:, 0] * 1.08), r"$\mathbf{1}$", color=COLORS["muted"])
ax.text(*(X[:, 1] * 1.05), r"$\mathbf{x}$", color=COLORS["muted"])
ax.text(2.3, 2.9, 5.6, r"$\mathcal{C}(\mathbf{X})$", color=COLORS["accent"])
ax.set_xlabel("observation 1", labelpad=-12)
ax.set_ylabel("observation 2", labelpad=-12)
ax.set_zlabel("observation 3", labelpad=-12)
ax.set_xticklabels([]); ax.set_yticklabels([]); ax.set_zticklabels([])
for axis in (ax.xaxis, ax.yaxis, ax.zaxis):
    axis.set_pane_color((1, 1, 1, 0))
    axis._axinfo["grid"]["color"] = (0.85, 0.87, 0.9, 1)
ax.view_init(elev=22, azim=-35)
ax.set_box_aspect(None, zoom=0.9)
fig.savefig(figure_path("ch06", "ls_geometry_3d"))
