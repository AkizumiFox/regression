"""Chapter 31, Section 3: a regression whose response has a singular covariance matrix.

A balanced one-way layout is summarized by the estimated treatment effects under the
sum-to-zero side condition; those effects add to zero exactly, so their covariance matrix
is singular. Regressing them on the dose is a general Gauss-Markov model in which one
linear function of the coefficients is known without error.
"""
import matplotlib.pyplot as plt
import numpy as np

from regbook import COLORS, Generated, figure_path, use_book_style

# <<setup>>
import numpy as np

dose = np.array([0.0, 1.0, 2.0, 4.0, 8.0])           # a = 5 treatments
a, reps, sigma = len(dose), 6, 1.2
X = np.column_stack([np.ones(a), dose])
V = np.eye(a) - np.ones((a, a)) / a                  # Cov(effects) = (sigma^2/reps) V
xbar = dose.mean()
Sxx = float(((dose - xbar) ** 2).sum())

Vplus = V                                            # V is symmetric idempotent, so V+ = V
N = np.ones((a, 1)) / np.sqrt(a)                     # an orthonormal basis of C(V) perp
print("rank(V) =", np.linalg.matrix_rank(V), " rank([X V]) =",
      np.linalg.matrix_rank(np.column_stack([X, V])), " rank(X) =", np.linalg.matrix_rank(X))
print("the perfectly known function:", (X.T @ N).ravel())
# <</setup>>

assert np.linalg.matrix_rank(np.column_stack([X, V])) == a          # C(X : V) is everything
assert np.allclose(V @ V, V) and np.allclose(V @ np.ones(a), 0)
df = np.linalg.matrix_rank(np.column_stack([X, V])) - np.linalg.matrix_rank(X)
assert df == a - 2


def blue(y):
    """The BLUE of X beta: minimize (y - Xb)' V+ (y - Xb) over b with y - Xb in C(V)."""
    from scipy.optimize import minimize
    con = {"type": "eq", "fun": lambda b: N.T @ (y - X @ b)}
    obj = lambda b: float((y - X @ b) @ Vplus @ (y - X @ b))
    out = minimize(obj, np.zeros(2), constraints=[con], tol=1e-14)
    return out.x


def rao(y):
    """Rao's unified formula, with T = V + X X' and the Moore-Penrose inverse of T."""
    T = V + X @ X.T
    Tp = np.linalg.pinv(T)
    A = X @ np.linalg.pinv(X.T @ Tp @ X) @ X.T @ Tp
    assert np.allclose(A @ X, X) and np.allclose(A @ V @ (np.eye(a) - X @ np.linalg.pinv(X)), 0,
                                                 atol=1e-9)
    return A @ y


def two_step(y, Kginv=None):
    """The explicit two-step form: pin down b_* from the exact relations, then GLS inside U."""
    K = N.T @ X                                      # 1 x 2: the perfectly known function
    if Kginv is None:
        Kginv = np.linalg.pinv(K)
    b_star = Kginv @ (N.T @ y)
    U = (dose - xbar)[:, None]                       # a basis of C(X) cap C(V)
    r = y - X @ b_star
    return X @ b_star + U @ np.linalg.solve(U.T @ Vplus @ U, U.T @ Vplus @ r)


def pinv_gls(Xm, Vp, y):
    """What software does when it silently replaces V^{-1} by a pseudo-inverse."""
    return np.linalg.pinv(Xm.T @ Vp @ Xm) @ Xm.T @ Vp @ y


# <<simulate>>
rng = np.random.default_rng(3107)
trials = 200_000
cell = rng.normal(scale=sigma / np.sqrt(reps), size=(trials, a))     # treatment means
cell += 0.35 * (dose - xbar)                                         # the true effects
effects = cell - cell.mean(axis=1, keepdims=True)                    # sum to zero exactly
B = effects @ np.linalg.solve(X.T @ X, X.T).T                        # ordinary least squares

tau2 = sigma**2 / reps
print("Var(slope):     %.5f simulated, %.5f from tau^2/Sxx" % (B[:, 1].var(), tau2 / Sxx))
print("Var(intercept): %.5f simulated, %.5f from xbar^2 tau^2/Sxx"
      % (B[:, 0].var(), xbar**2 * tau2 / Sxx))
print("what OLS reports for the intercept: %.5f" % (tau2 * (1 / a + xbar**2 / Sxx)))
# <</simulate>>

assert np.allclose(effects.sum(axis=1), 0, atol=1e-10)
assert np.allclose(a * B[:, 0] + dose.sum() * B[:, 1], 0, atol=1e-9)   # known without error
# a single data set: the three routes to the BLUE agree with ordinary least squares
y0 = effects[0]
b_ols = np.linalg.solve(X.T @ X, X.T @ y0)
assert np.allclose(blue(y0), b_ols, atol=1e-6)
assert np.allclose(rao(y0), X @ b_ols, atol=1e-8)
assert np.allclose(two_step(y0), X @ b_ols, rtol=0, atol=1e-10)
# and the two-step form does not depend on which generalized inverse of K is used
K = N.T @ X
for Kg in (np.array([[1.0 / K[0, 0]], [0.0]]), np.array([[0.0], [1.0 / K[0, 1]]])):
    assert np.allclose(K @ Kg, 1.0)
    assert np.allclose(two_step(y0, Kg), X @ b_ols, rtol=0, atol=1e-10)

# a pseudo-inverse plugged into the GLS formula solves the wrong problem: it minimizes
# over all b instead of over the feasible set, and can return a different point estimate
b_pinv = pinv_gls(X, Vplus, y0)
assert abs(b_pinv[0]) < 1e-12                                 # it reports a zero intercept
assert abs(a * b_pinv[0] + dose.sum() * b_pinv[1]) > 1e-3     # violating the exact relation
assert abs(b_pinv[0] - b_ols[0]) > 1e-3                       # so it is not the BLUE
# the smallest example of the same failure: the second observation is measured exactly
X2, V2, y2 = np.array([[1.0], [1.0]]), np.diag([1.0, 0.0]), np.array([3.0, 5.0])
assert np.isclose(pinv_gls(X2, V2, y2)[0], 3.0)     # minimizes (3 - b)^2 over every b
assert np.isclose(y2[1], 5.0)                       # but Y_2 = beta with probability one
# Kruskal's condition holds here, which is why ordinary least squares is already best
M = X @ np.linalg.solve(X.T @ X, X.T)
assert np.allclose((np.eye(a) - M) @ V @ X, 0, atol=1e-12)

var_slope = tau2 / Sxx
var_intercept = xbar**2 * var_slope
cov_01 = -xbar * var_slope
reported = tau2 * (1 / a + xbar**2 / Sxx)        # E(s^2) = tau^2 times the usual factor
assert abs(B[:, 1].var() / var_slope - 1) < 0.01
assert abs(B[:, 0].var() / var_intercept - 1) < 0.01
assert abs(np.cov(B.T)[0, 1] / cov_01 - 1) < 0.01
assert reported > 1.5 * var_intercept

# the residual sum of squares in the V-plus geometry is the ordinary one
E = effects - B @ X.T
rss_V = np.einsum("ij,jk,ik->i", E, Vplus, E)
assert np.allclose(rss_V, np.sum(E**2, axis=1))
assert abs(rss_V.mean() / (tau2 * df) - 1) < 0.02
assert abs(rss_V.mean() / (a - 2) / tau2 - 1) < 0.02   # E(s^2) = tau^2 here

gen = Generated("ch31", "singular_effects", prefix="sing")
gen.int("a", a)
gen.int("reps", reps)
gen.int("df", df)
gen.num("sigma", sigma, 1)
gen.num("xbar", xbar, 1)
gen.num("Sxx", Sxx, 1)
gen.num("tau2", tau2, 2)
gen.num("var_slope", var_slope, 5)
gen.num("var_slope_sim", B[:, 1].var(), 5)
gen.num("var_intercept", var_intercept, 4)
gen.num("var_intercept_sim", B[:, 0].var(), 4)
gen.num("reported", reported, 4)
gen.num("reported_ratio", reported / var_intercept, 2)
gen.num("rss_mean", rss_V.mean(), 3)
gen.write()

# ---- Figure: the geometry of consistency, with three observations ------------
use_book_style()
V3 = np.diag([1.0, 1.0, 0.0])                        # C(V) is the horizontal plane
X3 = np.column_stack([[1.0, 1.0, 0.0], [0.0, 0.0, 1.0]])
y3 = np.array([2.7, 0.1, 1.5])
d = X3[:, 0] / np.sqrt(2.0)                          # spans C(X) cap C(V)
alpha = (y3[0] + y3[1]) / 2
fit3 = np.array([alpha, alpha, y3[2]])               # the BLUE of X beta

assert np.linalg.matrix_rank(np.column_stack([X3, V3])) == 3      # y is always consistent
assert np.isclose(np.dot(y3 - fit3, d), 0)           # the residual is V+-orthogonal to the line
assert np.isclose((y3 - fit3)[2], 0)                 # and lies in C(V)

fig = plt.figure(figsize=(4.2, 3.0))
ax = fig.add_subplot(projection="3d")
s, t = np.meshgrid(np.linspace(-0.6, 2.6, 2), np.linspace(-0.4, 2.2, 2))
P = s[..., None] * X3[:, 0] + t[..., None] * X3[:, 1]
ax.plot_surface(P[..., 0], P[..., 1], P[..., 2], color=COLORS["accent"], alpha=0.16,
                linewidth=0)
u, v = np.meshgrid(np.linspace(-0.6, 3.0, 2), np.linspace(-0.6, 3.0, 2))
ax.plot_surface(u, v, np.zeros_like(u), color=COLORS["muted"], alpha=0.10, linewidth=0)
ends = np.array([[-0.5, -0.5, 0.0], [2.6, 2.6, 0.0]])
ax.plot(*ends.T, color=COLORS["muted"], linewidth=1.0)
ax.plot(*(ends + [0, 0, y3[2]]).T, color=COLORS["ink"], linewidth=1.2, linestyle="--")
ax.quiver(0, 0, 0, *y3, color=COLORS["ink"], linewidth=1.1, arrow_length_ratio=0.08)
ax.scatter(*fit3, color=COLORS["accent"], s=14, depthshade=False)
ax.quiver(*fit3, *(y3 - fit3), color=COLORS["second"], linewidth=1.2, arrow_length_ratio=0.16)
ax.text(*(y3 + [0.10, -0.15, 0.10]), r"$\mathbf{y}$", color=COLORS["ink"])
ax.text(*(fit3 + [-0.95, 0.0, 0.16]), r"$\mathbf{X}\tilde{\boldsymbol{\beta}}$",
        color=COLORS["accent"], fontsize=8)
ax.text(-0.5, -0.5, y3[2] + 0.22, r"$\mathcal{C}(\mathbf{X})\cap(\mathbf{y}+\mathcal{C}(\mathbf{V}))$",
        color=COLORS["ink"], fontsize=7, ha="left")
ax.text(2.6, 2.6, -0.50, r"$\mathcal{C}(\mathbf{X})\cap\mathcal{C}(\mathbf{V})$",
        color=COLORS["muted"], fontsize=7, ha="right")
ax.text(2.3, 2.3, 2.05, r"$\mathcal{C}(\mathbf{X})$", color=COLORS["accent"], fontsize=8,
        ha="right")
ax.text(2.9, -0.55, 0.08, r"$\mathcal{C}(\mathbf{V})$", color=COLORS["muted"], fontsize=7,
        ha="right")
ax.set_xlim(-0.5, 3.0)
ax.set_ylim(-0.5, 3.0)
ax.set_zlim(-0.3, 2.3)
ax.set_box_aspect((1, 1, 0.8))
ax.view_init(elev=20, azim=-20)
ax.set_axis_off()
fig.tight_layout()
fig.savefig(figure_path("ch31", "singular_geometry"))
