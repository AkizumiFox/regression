"""Chapter 43, Section 4: the smoothing spline, the Reinsch algorithm, and its rows.

Builds the Reinsch matrices Q and R from the design, checks the penalty against an
independently computed natural cubic interpolant, checks that minimizing over the whole
cubic spline space with the exact integral penalty gives the same fit (Theorem 43.4),
and draws the rows of the smoother matrix.
"""
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import CubicSpline
from scipy.linalg import solveh_banded

from regbook import COLORS, Generated, figure_path, use_book_style


def f_true(x):
    """The mean function of Section 43.1."""
    return 2 * x + np.exp(-25 * (x - 0.35) ** 2) - 0.75 * np.exp(-50 * (x - 0.80) ** 2)


# <<reinsch>>
def reinsch(t):
    """The band matrices Q (N x (N-2)) and R ((N-2) x (N-2)) of the Reinsch form."""
    N = len(t)
    h = np.diff(t)
    Q = np.zeros((N, N - 2))
    R = np.zeros((N - 2, N - 2))
    for j in range(N - 2):                       # column j belongs to interior knot t[j+1]
        Q[j, j], Q[j + 1, j], Q[j + 2, j] = 1 / h[j], -1 / h[j] - 1 / h[j + 1], 1 / h[j + 1]
        R[j, j] = (h[j] + h[j + 1]) / 3
        if j + 1 < N - 2:
            R[j, j + 1] = R[j + 1, j] = h[j + 1] / 6
    return Q, R


def smoothing_spline(t, y, lam):
    """Fitted values of the smoothing spline, by solving a pentadiagonal system."""
    Q, R = reinsch(t)
    gamma = np.linalg.solve(R + lam * Q.T @ Q, Q.T @ y)   # the second derivatives
    return y - lam * Q @ gamma, gamma


def smoother_matrix(t, lam):
    """The matrix S_lambda with fitted values S_lambda y (formed only to look at it)."""
    Q, R = reinsch(t)
    return np.eye(len(t)) - lam * Q @ np.linalg.solve(R + lam * Q.T @ Q, Q.T)
# <</reinsch>>


n, sigma = 120, 0.25
rng = np.random.default_rng(430430)
t = np.sort(rng.uniform(0, 1, n))
y = f_true(t) + sigma * rng.normal(size=n)

# --- the penalty matrix reproduces the roughness of the natural interpolant ----
Q, R = reinsch(t)
z = f_true(t)          # interpolate the true curve, whose roughness is finite
gamma = np.linalg.solve(R, Q.T @ z)
penalty_reinsch = gamma @ R @ gamma

spline = CubicSpline(t, z, bc_type="natural")            # an independent implementation
second = spline.derivative(2)
pieces = [np.polynomial.polynomial.Polynomial([second(t[i]), (second(t[i + 1]) - second(t[i])) / (t[i + 1] - t[i])])
          for i in range(n - 1)]
penalty_direct = sum(np.polynomial.polynomial.Polynomial.integ((p * p))(t[i + 1] - t[i])
                     for i, p in enumerate(pieces))
assert abs(penalty_reinsch - penalty_direct) < 1e-8 * penalty_direct

# --- the fit, and the smoother matrix ----------------------------------------
lam = 1e-4
fit, _ = smoothing_spline(t, y, lam)
S = smoother_matrix(t, lam)
assert np.allclose(fit, S @ y, atol=1e-9)
assert np.allclose(S, S.T, atol=1e-9)
df = np.trace(S)
print("effective degrees of freedom at lambda = 1e-4:", round(df, 2))
print("rows of S sum to one, largest departure:", np.max(np.abs(S.sum(axis=1) - 1)))
assert np.max(np.abs(S.sum(axis=1) - 1)) < 1e-8      # constants are reproduced
assert np.max(np.abs(S @ t - t)) < 1e-8              # and so are straight lines

# a pentadiagonal system: nothing outside the five central diagonals
M = R + lam * Q.T @ Q
i, j = np.indices(M.shape)
assert np.all(M[np.abs(i - j) > 2] == 0.0)
ab = np.zeros((3, M.shape[0]))
for k in range(3):
    ab[2 - k, k:] = np.diag(M, k)
gamma_band = solveh_banded(ab, Q.T @ y)
fit_band = y - lam * Q @ gamma_band
band_gap = np.max(np.abs(fit_band - fit))
assert band_gap < 1e-9

# --- Theorem 43.4 checked against the whole cubic spline space ----------------
def bspline_basis(x, kn, d):
    """The B-splines of degree d on the knot vector kn, by the de Boor recurrence."""
    x = np.atleast_1d(np.asarray(x, float))
    B = np.array([(x >= kn[j]) & (x < kn[j + 1]) for j in range(len(kn) - 1)], float).T
    last = np.max(np.nonzero(kn[:-1] < kn[1:])[0])
    B[x == kn[-1], last] = 1.0
    for deg in range(1, d + 1):
        C = np.zeros((len(x), B.shape[1] - 1))
        for jj in range(C.shape[1]):
            if kn[jj + deg] > kn[jj]:
                C[:, jj] += (x - kn[jj]) / (kn[jj + deg] - kn[jj]) * B[:, jj]
            if kn[jj + deg + 1] > kn[jj + 1]:
                C[:, jj] += (kn[jj + deg + 1] - x) / (kn[jj + deg + 1] - kn[jj + 1]) * B[:, jj + 1]
        B = C
    return B


def deriv_matrix(kn, d):
    """G with (B^d)' = B^{d-1} G: the derivative formula of Proposition 43.2 in matrix form."""
    m = len(kn) - d - 1
    G = np.zeros((m + 1, m))
    for j in range(m):
        if kn[j + d] > kn[j]:
            G[j, j] = d / (kn[j + d] - kn[j])
        if kn[j + d + 1] > kn[j + 1]:
            G[j + 1, j] = -d / (kn[j + d + 1] - kn[j + 1])
    return G


kn = np.concatenate([np.full(4, t[0]), t[1:-1], np.full(4, t[-1])])   # a knot at every design point
Bt = bspline_basis(t, kn, 3)
assert Bt.shape[1] == n + 2                                # cubic splines: dimension N + 2
nodes, weights = np.polynomial.legendre.leggauss(3)        # exact for the quadratic integrand
Omega = np.zeros((n + 2, n + 2))
for i in range(n - 1):
    mid, half = 0.5 * (t[i] + t[i + 1]), 0.5 * (t[i + 1] - t[i])
    u = mid + half * nodes
    Bpp = bspline_basis(u, kn, 1) @ deriv_matrix(kn, 2) @ deriv_matrix(kn, 3)
    Omega += half * (Bpp * (weights[:, None])).T @ Bpp
coef = np.linalg.solve(Bt.T @ Bt + lam * Omega, Bt.T @ y)
fit_basis = Bt @ coef
basis_gap = np.max(np.abs(fit_basis - fit))
assert basis_gap < 1e-4 * np.ptp(y)

# --- three amounts of smoothing ----------------------------------------------
def lambda_for_df(t, target):
    lo, hi = -12.0, 6.0
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        if np.trace(smoother_matrix(t, 10.0**mid)) > target:
            lo = mid
        else:
            hi = mid
    return 10.0 ** (0.5 * (lo + hi))


fits = {}
for target in (4.0, 9.0, 30.0):
    lamt = lambda_for_df(t, target)
    fits[target] = (lamt, smoothing_spline(t, y, lamt)[0])

# --- the rows of the smoother matrix, at three target points -----------------
S9 = smoother_matrix(t, fits[9.0][0])
rows = [int(np.argmin(np.abs(t - c))) for c in (0.0, 0.3, 0.6)]
assert np.sum(S9[rows[2]] < 0) >= 2                    # the weights change sign
assert np.argmax(S9[rows[2]]) == rows[2]               # and peak at the target point
assert np.isclose(S9[rows[2]].sum(), 1.0, atol=1e-8)   # the rows sum to one

gen = Generated("ch43", "smoothing_spline")
gen.int("n", n)
gen.num("sigma", sigma, 2)
gen.num("lam", lam, 0, sci=True)
gen.num("df", df, 2)
gen.num("penalty", penalty_reinsch, 1)
gen.num("penalty_rel", abs(penalty_reinsch - penalty_direct) / penalty_direct, 1, sci=True)
gen.num("penalty_gap", abs(penalty_reinsch - penalty_direct), 1, sci=True)
gen.num("band_gap", band_gap, 1, sci=True)
gen.num("basis_gap", basis_gap, 1, sci=True)
gen.int("dim_cubic", n + 2)
gen.num("lam4", fits[4.0][0], 1, sci=True)
gen.num("lam9", fits[9.0][0], 1, sci=True)
gen.num("lam30", fits[30.0][0], 1, sci=True)
gen.num("minweight", S9[rows[2]].min(), 4)
gen.num("maxweight", S9[rows[2]].max(), 4)
gen.int("nneg", int(np.sum(S9[rows[2]] < 0)))
gen.write()

# ---- figure -----------------------------------------------------------------
use_book_style()
fig, axes = plt.subplots(1, 2, figsize=(5.8, 2.4))
ax = axes[0]
ax.scatter(t, y, s=6, color=COLORS["muted"], alpha=0.7, linewidths=0)
ax.plot(t, f_true(t), color=COLORS["ink"], linewidth=1.0, label="true $f$")
for target, col, ls in ((4.0, COLORS["third"], "--"), (9.0, COLORS["accent"], "-"), (30.0, COLORS["second"], "-")):
    ax.plot(t, fits[target][1], color=col, linestyle=ls, linewidth=1.1, label=f"df $={target:g}$")
ax.set_xlabel("$x$")
ax.set_ylabel("$y$")
ax.set_ylim(-0.9, 3.2)
ax.set_title("(a) smoothing splines")
ax.legend(frameon=False, loc="upper left", ncol=2, handlelength=1.4, columnspacing=1.0)
ax = axes[1]
ax.axhline(0.0, color=COLORS["grid"], linewidth=0.6, zorder=0)
for r, col, ls in zip(rows, (COLORS["second"], COLORS["accent"], COLORS["third"]), ("--", "-", ":")):
    ax.plot(t, S9[r], color=col, linestyle=ls, linewidth=1.0, label=f"row at $x={t[r]:.2f}$")
ax.set_xlabel("$x$")
ax.set_ylabel("weight")
ax.set_title("(b) rows of $S_\\lambda$, df $=9$")
ax.legend(frameon=False, loc="upper right", handlelength=1.6)
fig.tight_layout()
fig.savefig(figure_path("ch43", "smoothing_spline"))
