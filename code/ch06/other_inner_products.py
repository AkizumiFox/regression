"""Chapter 6, Section 9: projections for other inner products, and oblique projections."""
import matplotlib.pyplot as plt
import numpy as np

from regbook import COLORS, Generated, figure_path, use_book_style

rng = np.random.default_rng(9)
n = 8
X = np.column_stack([np.ones(n), np.arange(n, dtype=float)])
y = X @ np.array([1.0, 0.5]) + rng.normal(size=n)
t = np.arange(n)
V = 0.7 ** np.abs(t[:, None] - t[None, :])             # AR(1)-type covariance

# <<vproj>>
Vinv = np.linalg.inv(V)
P_V = X @ np.linalg.solve(X.T @ Vinv @ X, X.T @ Vinv)    # V^{-1}-orthogonal projection onto C(X)

print("idempotent:", np.allclose(P_V @ P_V, P_V))
print("symmetric: ", np.allclose(P_V, P_V.T))                   # False: oblique in the usual geometry
print("self-adjoint for <u,v> = u'V^{-1}v:", np.allclose(Vinv @ P_V, (Vinv @ P_V).T))

# GLS = ordinary least squares after whitening with a square root of V^{-1}
L = np.linalg.cholesky(V)                               # V = L L'
Xw, yw = np.linalg.solve(L, X), np.linalg.solve(L, y)
beta_whitened, *_ = np.linalg.lstsq(Xw, yw, rcond=None)
beta_gls = np.linalg.solve(X.T @ Vinv @ X, X.T @ Vinv @ y)
print("GLS via P_V equals OLS on whitened data:", np.allclose(beta_gls, beta_whitened))
# <</vproj>>

assert np.allclose(P_V @ P_V, P_V) and not np.allclose(P_V, P_V.T)
assert np.allclose(Vinv @ P_V, (Vinv @ P_V).T)
assert np.allclose(X @ beta_gls, P_V @ y)
assert np.allclose(beta_gls, beta_whitened)
# The V^{-1} residual is orthogonal to C(X) in the V^{-1} inner product ...
r = y - P_V @ y
assert np.allclose(X.T @ Vinv @ r, 0)
# ... and P_V y is the nearest point of C(X) in the V^{-1} norm, not the Euclidean one.
M = X @ np.linalg.solve(X.T @ X, X.T)
dist_V = lambda z: (y - z) @ Vinv @ (y - z)
dist_E = lambda z: (y - z) @ (y - z)
assert dist_V(P_V @ y) <= dist_V(M @ y) + 1e-12
assert dist_E(M @ y) <= dist_E(P_V @ y) + 1e-12

gen = Generated("ch06", "other_inner_products", prefix="ip")
gen.num("asym", np.abs(P_V - P_V.T).max(), 3)
gen.num("ols_slope", np.linalg.lstsq(X, y, rcond=None)[0][1], 3)
gen.num("gls_slope", beta_gls[1], 3)
gen.write()

# ---- figure: orthogonal vs oblique projection onto a line in R^2 -------------
use_book_style()
fig, axes = plt.subplots(1, 2, figsize=(5.2, 2.3))
u = np.array([1.0, 0.35]); u /= np.linalg.norm(u)
yv = np.array([1.2, 1.6])
w = np.array([0.55, 1.0])                                # oblique direction, not orthogonal to u
orth = (yv @ u) * u
# oblique: move from y along w until hitting span(u): y - s w = c u
c, s = np.linalg.solve(np.column_stack([u, w]), yv)
obl = c * u


def draw(ax, foot, title, color):
    ax.plot([-0.3 * u[0], 2.4 * u[0]], [-0.3 * u[1], 2.4 * u[1]], color=COLORS["accent"], alpha=0.3, lw=4)
    for a, b, col in [((0, 0), yv, COLORS["ink"]), ((0, 0), foot, COLORS["accent"]), (foot, yv, color)]:
        ax.annotate("", xy=b, xytext=a, arrowprops=dict(arrowstyle="-|>", color=col, lw=1.1, shrinkA=0, shrinkB=0, mutation_scale=8))
    ax.text(*(yv + [0.05, 0.03]), r"$\mathbf{y}$")
    ax.text(*(foot + [0.05, -0.2]), r"$\mathbf{P}\mathbf{y}$", color=COLORS["accent"])
    ax.text(2.4 * u[0], 2.4 * u[1] + 0.12, r"$\mathcal{S}$", color=COLORS["accent"])
    ax.set_aspect("equal"); ax.axis("off"); ax.set_title(title)
    ax.set_xlim(-0.4, 2.7); ax.set_ylim(-0.3, 1.9)


draw(axes[0], orth, "(a) orthogonal: shortest residual", COLORS["second"])
draw(axes[1], obl, "(b) oblique: along a fixed direction", COLORS["third"])
fig.tight_layout()
fig.savefig(figure_path("ch06", "oblique_projection"))
