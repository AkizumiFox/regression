"""Chapter 10, Section 7: rank-revealing QR with column pivoting, its guarantees and its
failure on Kahan's matrix, and what different conventions do with aliased columns.

Data: statsmodels.datasets.grunfeld (public domain): investment of 11 US firms, 1935-1954.
"""
import warnings

import matplotlib.pyplot as plt
import numpy as np
import scipy.linalg as sla
import statsmodels.api as sm

from regbook import COLORS, Generated, figure_path, use_book_style

gen = Generated("ch10", "pivoted_qr", prefix="piv")


# <<pivoted>>
def pivoted_qr(X):
    """Householder QR with column pivoting: X[:, perm] = Q R, |r_11| >= |r_22| >= ..."""
    A = np.array(X, dtype=float)
    n, p = A.shape
    perm = np.arange(p)
    for k in range(min(n, p)):
        norms = np.sum(A[k:, k:] ** 2, axis=0)          # squared norms of the trailing columns
        j = k + int(np.argmax(norms))
        A[:, [k, j]], perm[[k, j]] = A[:, [j, k]], perm[[j, k]]
        x = A[k:, k]
        alpha = -np.copysign(np.linalg.norm(x), x[0])
        if alpha == 0.0:
            break                                        # the trailing block is exactly zero
        v = x.copy()
        v[0] -= alpha
        v /= np.linalg.norm(v)
        A[k:, k:] -= 2.0 * np.outer(v, v @ A[k:, k:])
    return np.triu(A[:p, :]), perm


def numerical_rank(R, tol):
    """Number of diagonal entries of the pivoted triangle above tol * |r_11|."""
    d = np.abs(np.diag(R))
    return int(np.sum(d > tol * d[0]))
# <</pivoted>>


rng = np.random.default_rng(107)
# ---- agreement with LAPACK and the guarantees of the theorem --------------------------------
for trial in range(20):
    Xt = rng.normal(size=(30, 8)) @ np.diag(10.0 ** rng.uniform(-6, 0, 8)) @ rng.normal(size=(8, 8))
    R, perm = pivoted_qr(Xt)
    _, Rs, Ps = sla.qr(Xt, pivoting=True, mode="economic")
    assert np.array_equal(perm, Ps)
    assert np.allclose(np.abs(R), np.abs(Rs), rtol=1e-8, atol=1e-12 * np.abs(Rs[0, 0]))
    d = np.abs(np.diag(R))
    assert np.all(np.diff(d) <= 1e-12 * d[0])                           # nonincreasing diagonal
    s = np.linalg.svd(Xt, compute_uv=False)
    for k in range(1, 8):
        R22 = R[k:, k:]
        assert d[k] ** 2 >= np.max(np.sum(R22 ** 2, axis=0)) * (1 - 1e-10)   # pivot dominates the rest
        assert s[k] <= np.linalg.norm(R22, 2) * (1 + 1e-10)                  # sigma_{k+1} <= ||R22||
        assert np.linalg.norm(R22, 2) <= np.sqrt(8 - k) * d[k] * (1 + 1e-10)
        assert np.linalg.svd(R[:k, :k], compute_uv=False)[-1] <= s[k - 1] * (1 + 1e-10)

# ---- Kahan's matrix: no small diagonal, yet nearly singular ---------------------------------
# <<kahan>>
m, c = 100, 0.2
sn = np.sqrt(1 - c ** 2)
K = np.diag(sn ** np.arange(m)) @ (np.eye(m) - c * np.triu(np.ones((m, m)), 1))
K = K @ np.diag(1 - 1e-10 * np.arange(m))           # tiny tilt so that pivoting keeps the order
RK, permK = pivoted_qr(K)
sK = np.linalg.svd(K, compute_uv=False)
print("pivot order unchanged:", np.array_equal(permK, np.arange(m)))
print(f"smallest |r_kk| = {np.abs(RK[-1, -1]):.3f},  smallest singular value = {sK[-1]:.2e}")
# <</kahan>>

assert np.array_equal(permK, np.arange(m))
assert np.abs(RK[-1, -1]) > 0.1 and sK[-1] < 1e-8
gen.num("kahan_rmin", np.abs(RK[-1, -1]), 3)
gen.num("kahan_smin", sK[-1], 1, sci=True)
gen.int("kahan_m", m)
# the bound sigma_min >= |r_mm| / (something exponential) still holds: sigma_m <= |r_mm|
assert sK[-1] <= np.abs(RK[-1, -1])

# a matrix of numerical rank 6 with noise: the diagonal reveals the gap
Xg = rng.normal(size=(50, 6)) @ rng.normal(size=(6, 10)) + 1e-9 * rng.normal(size=(50, 10))
Rg, _ = pivoted_qr(Xg)
sg = np.linalg.svd(Xg, compute_uv=False)
assert numerical_rank(Rg, 1e-7) == 6 and np.sum(sg > 1e-7 * sg[0]) == 6

# ---- aliased columns: the Grunfeld data with firm indicators --------------------------------
# <<aliased>>
gr = sm.datasets.grunfeld.load_pandas().data
firms = list(dict.fromkeys(gr["firm"]))
D = np.column_stack([(gr["firm"] == f).to_numpy(float) for f in firms])   # 11 indicators
y = gr["invest"].to_numpy()
X = np.column_stack([np.ones(len(y)), gr[["value", "capital"]], D])      # 14 columns, rank 13
n, p = X.shape
R, perm = pivoted_qr(X)
r = numerical_rank(R, 1e-10)
print("p =", p, " numerical rank =", r, " column left out:", perm[r:])
# <</aliased>>

assert r == 13 and np.linalg.matrix_rank(X) == 13


# <<conventions>>
def basic_solution(X, y, tol=1e-10):
    """Drop the columns that pivoted QR puts last; fit the rest (zeros for dropped columns)."""
    R, perm = pivoted_qr(X)
    r = numerical_rank(R, tol)
    b = np.zeros(X.shape[1])
    b[perm[:r]] = np.linalg.lstsq(X[:, perm[:r]], y, rcond=None)[0]
    return b


def in_order_solution(X, y, tol=1e-7):
    """Keep a column only if it is not (nearly) a combination of the columns kept before it."""
    kept = []
    for j in range(X.shape[1]):
        Z = X[:, kept + [j]]
        Q, R = np.linalg.qr(Z)
        if abs(R[-1, -1]) > tol * np.linalg.norm(X[:, j]):
            kept.append(j)
    b = np.full(X.shape[1], np.nan)                  # aliased coefficients reported as missing
    b[kept] = np.linalg.lstsq(X[:, kept], y, rcond=None)[0]
    return b


solutions = {
    "minimum norm (pinv)": np.linalg.pinv(X) @ y,
    "pivoted QR, basic": basic_solution(X, y),
    "in order, aliased = NA": in_order_solution(X, y),
}
for name, b in solutions.items():
    fitted = X @ np.nan_to_num(b)
    print(f"{name:24s} value {b[1]:.5f}  capital {b[2]:.5f}  const {b[0]:9.3f}"
          f"  GM - US Steel {b[3] - b[4]:9.3f}  SSE {np.sum((y - fitted) ** 2):.1f}")
# <</conventions>>

ref = solutions["minimum norm (pinv)"]
fits = [X @ np.nan_to_num(b) for b in solutions.values()]
for b, f in zip(solutions.values(), fits):
    assert np.allclose(f, fits[0], rtol=1e-8, atol=1e-6)
    assert np.allclose(b[1:3], ref[1:3], rtol=1e-8)                        # slopes are estimable
    assert np.isclose(b[3] - b[4], ref[3] - ref[4], rtol=1e-8)             # firm contrasts too
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    sm_fit = sm.OLS(y, X).fit()
assert np.allclose(sm_fit.params, ref, rtol=1e-8, atol=1e-8)            # statsmodels: minimum norm
b_io = solutions["in order, aliased = NA"]
assert np.isnan(b_io[-1]) and np.sum(np.isnan(b_io)) == 1                # the last indicator is aliased
b_basic = solutions["pivoted QR, basic"]
assert np.sum(b_basic == 0) == 1
dropped = firms[int(np.where(b_basic == 0)[0][0]) - 3] if np.where(b_basic == 0)[0][0] >= 3 else "constant"
gen.text("dropped_basic", dropped)
gen.num("value", ref[1], 4)
gen.num("capital", ref[2], 4)
gen.num("contrast", ref[3] - ref[4], 2)
gen.num("sse", np.sum((y - fits[0]) ** 2), 0)
gen.num("const_mn", ref[0], 2)
gen.num("const_io", b_io[0], 2)
gen.num("const_basic", b_basic[0], 2)
gen.int("n", n)
gen.int("p", p)

# ---- a nearly aliased column: the tolerance decides ------------------------------------------
# <<tolerance>>
rng = np.random.default_rng(108)
z = X[:, 1] + X[:, 2] + 1e-9 * np.linalg.norm(X[:, 1]) * rng.normal(size=n) / np.sqrt(n)
Xz = np.column_stack([X[:, :3], z])              # constant, value, capital, value + capital + tiny
s = np.linalg.svd(Xz, compute_uv=False)
print("relative singular values:", np.array2string(s / s[0], precision=2))
for tol in [None, 1e-8]:
    b, _, rank, _ = np.linalg.lstsq(Xz, y, rcond=tol)
    print(f"rcond = {tol}: rank {rank}, largest |b_j| = {np.abs(b).max():.3g},"
          f" SSE = {np.sum((y - Xz @ b) ** 2):.1f}")
# <</tolerance>>

b_def, _, rank_def, _ = np.linalg.lstsq(Xz, y, rcond=None)
b_tol, _, rank_tol, _ = np.linalg.lstsq(Xz, y, rcond=1e-8)
assert rank_def == 4 and rank_tol == 3
assert np.abs(b_def).max() > 1e3 * np.abs(b_tol).max()
sse_def, sse_tol = np.sum((y - Xz @ b_def) ** 2), np.sum((y - Xz @ b_tol) ** 2)
assert abs(sse_def - sse_tol) / sse_tol < 1e-2
gen.num("near_smin", s[-1] / s[0], 1, sci=True)
gen.num("near_bmax_def", np.abs(b_def).max(), 1, sci=True)
gen.num("near_bmax_tol", np.abs(b_tol).max(), 2)
gen.num("near_sse_def", sse_def, 1)
gen.num("near_sse_tol", sse_tol, 1)
gen.num("near_sse_pct", 100 * (sse_tol - sse_def) / sse_tol, 2)
gen.write()

use_book_style()
fig, axes = plt.subplots(1, 2, figsize=(5.8, 2.5))
ax = axes[0]
ax.semilogy(np.arange(1, m + 1), np.abs(np.diag(RK)), color=COLORS["accent"], label=r"$|r_{kk}|$")
ax.semilogy(np.arange(1, m + 1), sK, color=COLORS["second"], ls="--", label=r"$\sigma_k$")
ax.set_xlabel(r"$k$")
ax.set_title("(a) Kahan's matrix")
ax.legend(frameon=False, loc="lower left")
ax = axes[1]
k10 = np.arange(1, 11)
ax.semilogy(k10, np.abs(np.diag(Rg)), "o-", ms=3, color=COLORS["accent"], label=r"$|r_{kk}|$")
ax.semilogy(k10, sg, "s--", ms=3, color=COLORS["second"], label=r"$\sigma_k$")
ax.axhline(1e-7 * sg[0], color=COLORS["muted"], ls=":", lw=0.8)
ax.text(10, 1e-7 * sg[0] * 3, "tolerance", ha="right", fontsize=7, color=COLORS["muted"])
ax.set_xlabel(r"$k$")
ax.set_title("(b) numerical rank 6 plus noise")
ax.legend(frameon=False, loc="lower left")
fig.tight_layout()
fig.savefig(figure_path("ch10", "pivoted_qr"))
