"""Chapter 3, Section 5: conditional distributions.

(a) A bivariate normal: density contours, the two conditional-mean lines, the
    major axis, and a least squares line from a sample.
(b) An exact trivariate example (fractions) checked by simulation."""
from fractions import Fraction as F

import matplotlib.pyplot as plt
import numpy as np

from regbook import COLORS, Generated, figure_path, use_book_style

rng = np.random.default_rng(305)
gen = Generated("ch03", "conditional", prefix="cond")

# ---- (a) bivariate picture -------------------------------------------------
mx, my, sx, sy, rho = 1.0, 2.0, 1.0, 1.5, 0.7
Sigma = np.array([[sx**2, rho * sx * sy], [rho * sx * sy, sy**2]])
slope_y_on_x = rho * sy / sx                    # E(Y | X = x)
slope_x_on_y = rho * sx / sy                    # E(X | Y = y), as dx/dy
cond_sd = sy * np.sqrt(1 - rho**2)

n = 60
XY = rng.multivariate_normal([mx, my], Sigma, size=n)
b_ls, a_ls = np.polyfit(XY[:, 0], XY[:, 1], 1)

# vertical tangency points of the ellipse Delta^2 = c^2 lie on the line E(Y | X = x)
Linv = np.linalg.inv(np.linalg.cholesky(Sigma))
theta = np.linspace(0, 2 * np.pi, 200_001)
circle = np.vstack([np.cos(theta), np.sin(theta)])
L = np.linalg.cholesky(Sigma)
for c in [1.0, 2.0, 3.0]:
    pts = np.array([[mx], [my]]) + c * L @ circle
    i = np.argmax(pts[0])
    assert np.isclose(pts[0, i], mx + c * sx, atol=1e-6)
    assert np.isclose(pts[1, i], my + slope_y_on_x * (pts[0, i] - mx), atol=1e-3)

# conditional law by simulation: slice a large sample near x0
big = rng.multivariate_normal([mx, my], Sigma, size=4_000_000)
x0 = 2.0
sl = big[np.abs(big[:, 0] - x0) < 0.01, 1]
assert abs(sl.mean() - (my + slope_y_on_x * (x0 - mx))) < 0.02
assert abs(sl.std() - cond_sd) < 0.02
# least squares slope is within 3 standard errors of the population slope
resid = XY[:, 1] - a_ls - b_ls * XY[:, 0]
se = np.sqrt(resid @ resid / (n - 2) / np.sum((XY[:, 0] - XY[:, 0].mean())**2))
assert abs(b_ls - slope_y_on_x) < 3 * se
# the major axis lies strictly between the two regression lines
lam, U = np.linalg.eigh(Sigma)
slope_major = U[1, 1] / U[0, 1]
assert slope_y_on_x < slope_major < 1 / slope_x_on_y

gen.num("mx", mx, 0); gen.num("my", my, 0); gen.num("sx", sx, 0); gen.num("sy", sy, 1)
gen.num("rho", rho, 1)
gen.num("slope", slope_y_on_x, 2)
gen.num("slopexy", 1 / slope_x_on_y, 2)
gen.num("slopemajor", slope_major, 2)
gen.num("condsd", cond_sd, 3)
gen.int("n", n)
levels = [1, 4, 9]
gen.text("levels", ",".join(str(v) for v in levels))
gen.num("bls", b_ls, 2)
gen.num("se", se, 2)

use_book_style()
fig, ax = plt.subplots(figsize=(4.2, 3.3))
gx, gy = np.meshgrid(np.linspace(-3, 5, 300), np.linspace(-4, 8, 300))
D = np.stack([gx - mx, gy - my], axis=-1)
D2 = np.einsum("...i,ij,...j->...", D, np.linalg.inv(Sigma), D)
ax.contour(gx, gy, D2, levels=levels, colors=COLORS["muted"], linewidths=0.7)
ax.scatter(XY[:, 0], XY[:, 1], s=6, color=COLORS["muted"], alpha=0.6, linewidths=0)
xs = np.linspace(-2.2, 4.2, 2)
ax.plot(xs, my + slope_y_on_x * (xs - mx), color=COLORS["accent"], label=r"$\mathrm{E}(Y\mid X=x)$")
ys = np.linspace(-3.0, 7.0, 2)
ax.plot(mx + slope_x_on_y * (ys - my), ys, color=COLORS["third"], label=r"$\mathrm{E}(X\mid Y=y)$")
ax.plot(xs, my + slope_major * (xs - mx), color=COLORS["ink"], lw=0.8, ls=":", label="major axis")
ax.plot(xs, a_ls + b_ls * xs, color=COLORS["second"], ls="--", label=f"least squares, $n={n}$")
for c in np.sqrt(levels):
    for s in (-1, 1):
        xt = mx + s * c * sx
        ax.plot([xt], [my + slope_y_on_x * (xt - mx)], "o", ms=3, color=COLORS["accent"])
ax.set_xlim(-2.5, 4.5); ax.set_ylim(-3.5, 7.5)
ax.set_xlabel("$x$"); ax.set_ylabel("$y$")
ax.legend(frameon=False, loc="upper left", fontsize=7.5)
fig.savefig(figure_path("ch03", "conditional"))

# ---- (b) exact trivariate example --------------------------------------------
# <<exact>>
mu = [F(1), F(0), F(-1)]
S = [[F(4), F(1), F(2)],
     [F(1), F(3), F(2)],
     [F(2), F(2), F(2)]]

def inv2(M):
    det = M[0][0] * M[1][1] - M[0][1] * M[1][0]
    return [[M[1][1] / det, -M[0][1] / det], [-M[1][0] / det, M[0][0] / det]]

# Y1 given (Y2, Y3): B = sigma_12 Sigma_22^{-1}, variance sigma_11 - B sigma_21
S22inv = inv2([[S[1][1], S[1][2]], [S[2][1], S[2][2]]])
s12 = [S[0][1], S[0][2]]
B = [s12[0] * S22inv[0][j] + s12[1] * S22inv[1][j] for j in range(2)]
var_1_23 = S[0][0] - (B[0] * s12[0] + B[1] * s12[1])
const = mu[0] - B[0] * mu[1] - B[1] * mu[2]
print("E(Y1 | y2, y3) =", const, "+", B[0], "y2 +", B[1], "y3;  variance", var_1_23)

# (Y1, Y2) given Y3, and the partial correlation of Y1 and Y2 given Y3
C = [[S[i][j] - S[i][2] * S[2][j] / S[2][2] for j in range(2)] for i in range(2)]
print("Cov((Y1, Y2) | Y3) =", C)
# <</exact>>

Snp = np.array(S, dtype=float)
assert np.all(np.linalg.eigvalsh(Snp) > 0)
Bnp = Snp[0, 1:] @ np.linalg.inv(Snp[1:, 1:])
assert np.allclose(Bnp, [float(b) for b in B])
assert np.isclose(float(var_1_23), np.linalg.det(Snp) / np.linalg.det(Snp[1:, 1:]))
rho12 = float(S[0][1]) / np.sqrt(float(S[0][0] * S[1][1]))
rho12_3 = float(C[0][1]) / np.sqrt(float(C[0][0] * C[1][1]))
Om = np.linalg.inv(Snp)
assert np.isclose(rho12_3, -Om[0, 1] / np.sqrt(Om[0, 0] * Om[1, 1]))
R = Snp / np.sqrt(np.outer(np.diag(Snp), np.diag(Snp)))
r12, r13, r23 = R[0, 1], R[0, 2], R[1, 2]
assert np.isclose(rho12_3, (r12 - r13 * r23) / np.sqrt((1 - r13**2) * (1 - r23**2)))
# simulation: regress Y1 on (1, Y2, Y3); residual correlation of Y1, Y2 given Y3
Ysim = rng.multivariate_normal(np.array(mu, dtype=float), Snp, size=1_000_000)
Xd = np.column_stack([np.ones(len(Ysim)), Ysim[:, 1:]])
coef, *_ = np.linalg.lstsq(Xd, Ysim[:, 0], rcond=None)
assert np.allclose(coef, [float(const), float(B[0]), float(B[1])], atol=0.01)
assert abs(np.var(Ysim[:, 0] - Xd @ coef) - float(var_1_23)) < 0.01
Z3 = np.column_stack([np.ones(len(Ysim)), Ysim[:, 2]])
e1 = Ysim[:, 0] - Z3 @ np.linalg.lstsq(Z3, Ysim[:, 0], rcond=None)[0]
e2 = Ysim[:, 1] - Z3 @ np.linalg.lstsq(Z3, Ysim[:, 1], rcond=None)[0]
assert abs(np.corrcoef(e1, e2)[0, 1] - rho12_3) < 0.005
assert rho12 > 0 > rho12_3


def tex(fr):
    fr = F(fr)
    if fr.denominator == 1:
        return f"{fr.numerator}" if fr >= 0 else f"\\ensuremath{{-}}{-fr.numerator}"
    sgn = "-" if fr < 0 else ""
    return f"\\ensuremath{{{sgn}\\tfrac{{{abs(fr.numerator)}}}{{{fr.denominator}}}}}"


gen.text("Bone", tex(B[0])); gen.text("Btwo", tex(B[1]))


def affine(c0, coefs, names):
    """Typeset c0 + sum coef * name with integer coefficients, e.g. 3 - y_2 + 2y_3."""
    out = f"{c0}"
    for c, v in zip(coefs, names):
        assert F(c).denominator == 1
        c = int(c)
        mag = "" if abs(c) == 1 else f"{abs(c)}"
        out += f"{'-' if c < 0 else '+'}{mag}{v}"
    return out


gen.text("condmean", affine(int(const), B, ["y_2", "y_3"]))
gen.text("const", tex(const)); gen.text("var", tex(var_1_23))
gen.text("Caa", tex(C[0][0])); gen.text("Cab", tex(C[0][1])); gen.text("Cbb", tex(C[1][1]))
gen.text("marginalslope", tex(S[0][1] / S[1][1]))
gen.text("marginalvar", tex(S[0][0] - S[0][1]**2 / S[1][1]))
gen.num("rho12", rho12, 3)
gen.num("rho12three", rho12_3, 3)
gen.write()
