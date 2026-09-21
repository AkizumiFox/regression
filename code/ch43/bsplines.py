"""Chapter 43, Section 2: the B-spline basis, its properties, and its conditioning.

Implements the de Boor recurrence in numpy, checks every claim of Proposition 43.2
numerically, and compares the conditioning of the B-spline and truncated power bases
for the same knots.
"""
import matplotlib.pyplot as plt
import numpy as np

from regbook import COLORS, Generated, figure_path, use_book_style

# <<basis>>
def knot_vector(a, b, K, d):
    """d+1 copies of a, K equally spaced interior knots, d+1 copies of b."""
    interior = a + (b - a) * np.arange(1, K + 1) / (K + 1)
    return np.concatenate([np.full(d + 1, a), interior, np.full(d + 1, b)])


def uniform_knots(a, b, K, d):
    """K interior knots equally spaced in (a, b), continued d further steps at each end."""
    delta = (b - a) / (K + 1)
    return a + delta * np.arange(-d, K + 2 + d)


def bspline_basis(x, kn, d):
    """The B-splines of degree d on the knot vector kn, by the de Boor recurrence."""
    x = np.atleast_1d(np.asarray(x, float))
    B = np.array([(x >= kn[j]) & (x < kn[j + 1]) for j in range(len(kn) - 1)], float).T
    last = np.max(np.nonzero(kn[:-1] < kn[1:])[0])       # close the rightmost interval
    B[x == kn[-1], last] = 1.0
    for deg in range(1, d + 1):
        C = np.zeros((len(x), B.shape[1] - 1))
        for j in range(C.shape[1]):
            if kn[j + deg] > kn[j]:                       # a term with a zero denominator
                C[:, j] += (x - kn[j]) / (kn[j + deg] - kn[j]) * B[:, j]
            if kn[j + deg + 1] > kn[j + 1]:               # contributes nothing
                C[:, j] += (kn[j + deg + 1] - x) / (kn[j + deg + 1] - kn[j + 1]) * B[:, j + 1]
        B = C
    return B


D, K, A, Bnd = 3, 6, 0.0, 1.0                            # cubic, six interior knots on [0, 1]
kn = knot_vector(A, Bnd, K, D)
# <</basis>>

grid = np.linspace(A, Bnd, 1001)
Bg = bspline_basis(grid, kn, D)

# --- the properties of Proposition 43.2 --------------------------------------
assert Bg.shape[1] == D + 1 + K                                  # the dimension count
assert np.all(Bg >= 0)                                           # nonnegativity
assert np.allclose(Bg.sum(axis=1), 1.0)                          # partition of unity
for j in range(Bg.shape[1]):                                     # local support
    inside = (grid >= kn[j]) & (grid <= kn[j + D + 1])
    assert np.all(Bg[~inside, j] == 0.0)
    assert np.any(Bg[inside, j] > 0.0)

# the derivative formula, checked against a central difference
eps = 1e-6
Bm, Bp = bspline_basis(grid[1:-1] - eps, kn, D), bspline_basis(grid[1:-1] + eps, kn, D)
numeric = (Bp - Bm) / (2 * eps)
lower = bspline_basis(grid[1:-1], kn, D - 1)                     # degree d-1, one more of them
formula = np.zeros_like(numeric)
for j in range(numeric.shape[1]):
    if kn[j + D] > kn[j]:
        formula[:, j] += D * lower[:, j] / (kn[j + D] - kn[j])
    if kn[j + D + 1] > kn[j + 1]:
        formula[:, j] -= D * lower[:, j + 1] / (kn[j + D + 1] - kn[j + 1])
assert np.max(np.abs(numeric - formula)) < 1e-5

# the Greville abscissae reproduce the identity function
greville = np.array([kn[j + 1:j + D + 1].mean() for j in range(D + 1 + K)])
assert np.allclose(Bg @ greville, grid, atol=1e-12)

# --- the same space as the truncated power basis ------------------------------
# <<equivalence>>
def truncated_power(x, interior, d):
    """The basis 1, x, ..., x^d, (x - kappa_1)_+^d, ..., (x - kappa_K)_+^d."""
    x = np.asarray(x, float)
    powers = [x**j for j in range(d + 1)]
    return np.column_stack(powers + [np.maximum(x - k, 0.0) ** d for k in interior])


rng = np.random.default_rng(43043)
n = 300
x = np.sort(rng.uniform(0, 1, n))
y = np.sin(7 * x) + 0.4 * x**2 + 0.3 * rng.normal(size=n)

interior = kn[D + 1:D + 1 + K]
Bx, Tx = bspline_basis(x, kn, D), truncated_power(x, interior, D)
fit_b = Bx @ np.linalg.lstsq(Bx, y, rcond=None)[0]
fit_t = Tx @ np.linalg.lstsq(Tx, y, rcond=None)[0]
print("largest difference between the two fits:", np.max(np.abs(fit_b - fit_t)))
print("condition numbers:", np.linalg.cond(Bx), np.linalg.cond(Tx))
# <</equivalence>>

assert np.max(np.abs(fit_b - fit_t)) < 1e-8
assert Bx.shape == Tx.shape

# --- the uniform layout used for P-splines spans the same functions on [a, b] --
Ux = bspline_basis(x, uniform_knots(A, Bnd, K, D), D)
fit_u = Ux @ np.linalg.lstsq(Ux, y, rcond=None)[0]
assert Ux.shape == Bx.shape and np.linalg.matrix_rank(Ux) == D + 1 + K
assert np.max(np.abs(fit_b - fit_u)) < 1e-8
kn_u = uniform_knots(A, Bnd, K, D)
delta = (Bnd - A) / (K + 1)
greville_u = np.array([kn_u[j + 1:j + D + 1].mean() for j in range(D + 1 + K)])
assert np.allclose(np.diff(greville_u), delta)                      # equally spaced
level = A + (np.arange(1, D + 2 + K) - (D + 1) / 2) * delta         # and at this level
assert np.allclose(greville_u, level)
assert np.allclose(bspline_basis(grid, kn_u, D) @ greville_u, grid, atol=1e-12)

# --- conditioning as the number of knots grows --------------------------------
rows = []
for Kk in (5, 10, 20, 40):
    knk = knot_vector(A, Bnd, Kk, D)
    cb = np.linalg.cond(bspline_basis(x, knk, D))
    ct = np.linalg.cond(truncated_power(x, knk[D + 1:D + 1 + Kk], D))
    rows.append((Kk, cb, ct))
    assert ct > cb
assert rows[-1][2] / rows[-1][1] > 1e4

# banded structure: at most d+1 nonzero entries in a row
assert np.max((Bx > 0).sum(axis=1)) <= D + 1

gen = Generated("ch43", "bsplines")
gen.int("d", D)
gen.int("K", K)
gen.int("dim", D + 1 + K)
gen.int("n", n)
gen.num("fit_gap", np.max(np.abs(fit_b - fit_t)), 1, sci=True)
gen.num("fit_gap_u", np.max(np.abs(fit_b - fit_u)), 1, sci=True)
for Kk, cb, ct in rows:
    gen.num(f"condb{Kk}", cb, 1)
    gen.num(f"condt{Kk}", ct, 2, sci=True)
gen.num("condb", np.linalg.cond(Bx), 1)                             # at the K of the example
gen.num("condt", np.linalg.cond(Tx), 1, sci=True)
gen.write()

# ---- figure: the basis, and a fit as a weighted sum of it --------------------
use_book_style()
fig, axes = plt.subplots(1, 2, figsize=(5.8, 2.4))
ax = axes[0]
for j in range(Bg.shape[1]):
    ax.plot(grid, Bg[:, j], color=COLORS["accent"], linewidth=1.0)
ax.plot(grid, Bg.sum(axis=1), color=COLORS["second"], linestyle="--", linewidth=1.0)
for k in interior:
    ax.axvline(k, color=COLORS["grid"], linewidth=0.6, zorder=0)
ax.set_ylim(0, 1.15)
ax.set_xlabel("$x$")
ax.set_title("(a) cubic B-splines, 6 interior knots")
ax.text(0.5, 1.04, "sum $=1$", color=COLORS["second"], ha="center", fontsize=8)
ax = axes[1]
coef = np.linalg.lstsq(Bx, y, rcond=None)[0]
ax.scatter(x, y, s=5, color=COLORS["muted"], alpha=0.6, linewidths=0)
for j in range(Bg.shape[1]):
    ax.plot(grid, coef[j] * Bg[:, j], color=COLORS["third"], linewidth=0.7)
ax.plot(grid, Bg @ coef, color=COLORS["accent"], linewidth=1.3)
ax.scatter(greville, coef, s=9, color=COLORS["second"], zorder=5)
for k in interior:
    ax.axvline(k, color=COLORS["grid"], linewidth=0.6, zorder=0)
ax.set_xlabel("$x$")
ax.set_title("(b) the fit and its scaled basis")
fig.tight_layout()
fig.savefig(figure_path("ch43", "bspline_basis"))
