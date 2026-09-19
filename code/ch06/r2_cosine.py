"""Chapter 6, Section 7: R^2 is a squared cosine.

With n = 3 and one regressor, the centred vectors live in the plane 1^perp,
so the whole picture is two-dimensional.
"""
import matplotlib.pyplot as plt
import numpy as np
import statsmodels.api as sm

from regbook import COLORS, Generated, figure_path, use_book_style

# <<cosine>>
def r_squared_as_cosine(y, X):
    """X must contain the intercept column. Returns (R^2 from sums of squares, cos^2 of angle)."""
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    y_hat = X @ beta
    yc, yhc = y - y.mean(), y_hat - y.mean()          # centre both vectors
    r2 = 1 - np.sum((y - y_hat) ** 2) / np.sum(yc ** 2)
    cos = (yc @ yhc) / (np.linalg.norm(yc) * np.linalg.norm(yhc))
    return r2, cos ** 2
# <</cosine>>

data = sm.datasets.statecrime.load_pandas().data.drop(index="District of Columbia")
y = data["murder"].to_numpy()
X = np.column_stack([np.ones(len(y)), data["poverty"], data["single"], data["urban"]])
r2, cos2 = r_squared_as_cosine(y, X)
assert np.isclose(r2, cos2, rtol=1e-12)
angle = np.degrees(np.arccos(np.sqrt(cos2)))

gen = Generated("ch06", "r2_cosine", prefix="r2")
gen.num("R2", r2, 4)
gen.num("angle", angle, 1)
# for simple regression R^2 is the squared correlation
r2s, _ = r_squared_as_cosine(y, X[:, :2])
corr = np.corrcoef(y, X[:, 1])[0, 1]
assert np.isclose(r2s, corr ** 2)
gen.num("R2simple", r2s, 4)
gen.num("corr", corr, 4)

# ---- figure: n = 3, the plane orthogonal to 1 --------------------------------
y3 = np.array([2.0, 1.0, 4.0])
x3 = np.array([0.0, 1.0, 3.0])
B = np.array([[1, -1, 0], [1, 1, -2]], dtype=float).T
B /= -np.linalg.norm(B, axis=0)                         # orthonormal basis for 1^perp (sign chosen for the picture)
yc, xc = y3 - y3.mean(), x3 - x3.mean()
proj = (xc @ yc) / (xc @ xc) * xc
cy, cx, cp = B.T @ yc, B.T @ xc, B.T @ proj
r2_3, cos2_3 = r_squared_as_cosine(y3, np.column_stack([np.ones(3), x3]))
assert np.isclose(r2_3, cos2_3)
gen.num("R2small", r2_3, 4)
gen.write()

use_book_style()
fig, ax = plt.subplots(figsize=(3.2, 2.6))


def arrow(a, b, color):
    ax.annotate("", xy=b, xytext=a,
                arrowprops=dict(arrowstyle="-|>", color=color, lw=1.2, shrinkA=0, shrinkB=0, mutation_scale=9))


d = cx / np.linalg.norm(cx)
L = 2.6
ax.plot([-0.4 * L * d[0], L * d[0]], [-0.4 * L * d[1], L * d[1]], color=COLORS["accent"], alpha=0.25, lw=4, solid_capstyle="butt")
arrow((0, 0), cy, COLORS["ink"])
arrow((0, 0), cp, COLORS["accent"])
arrow(cp, cy, COLORS["second"])
# right angle at the foot of the perpendicular
u, v = -cp / np.linalg.norm(cp) * 0.18, (cy - cp) / np.linalg.norm(cy - cp) * 0.18
ax.plot(*np.array([cp + u, cp + u + v, cp + v]).T, color=COLORS["ink"], lw=0.6)
a0, a1 = np.arctan2(cp[1], cp[0]), np.arctan2(cy[1], cy[0])
th = np.linspace(a0, a1, 30)
ax.plot(0.5 * np.cos(th), 0.5 * np.sin(th), color=COLORS["ink"], lw=0.7)
ax.text(0.68 * np.cos(th[15]), 0.68 * np.sin(th[15]), r"$\theta$", ha="center", va="center")
ax.text(*(cy + np.array([-0.1, 0.0])), r"$\mathbf{y}-\bar y\mathbf{1}$", ha="right", va="center")
ax.text(*(cp + np.array([0.15, -0.12])), r"$\hat{\mathbf{y}}-\bar y\mathbf{1}$", color=COLORS["accent"], ha="left", va="top")
ax.text(*((cp + cy) / 2 + np.array([0.0, 0.12])), r"$\hat{\boldsymbol{\varepsilon}}$", color=COLORS["second"], ha="center", va="bottom")
ax.text(*(L * d + np.array([0.05, 0])), r"$\mathrm{span}(\mathbf{x}-\bar x\mathbf{1})$", color=COLORS["accent"], ha="left", va="center")
ax.plot(0, 0, "o", ms=2.5, color=COLORS["ink"])
ax.set_aspect("equal")
pts = np.array([cy, cp, L * d, -0.4 * L * d, [0, 0]])
ax.set_xlim(pts[:, 0].min() - 1.3, pts[:, 0].max() + 1.6)
ax.set_ylim(pts[:, 1].min() - 0.3, pts[:, 1].max() + 0.3)
ax.axis("off")
ax.set_title(r"inside $\mathbf{1}^\perp$:  $R^2=\cos^2\theta=%.3f$" % r2_3)
fig.savefig(figure_path("ch06", "r2_cosine"))
