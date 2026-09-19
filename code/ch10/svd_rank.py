"""Chapter 10, Section 4: least squares through the SVD, the minimum-norm solution, the
discontinuity of the pseudoinverse, Weyl's inequality, and truncated SVD as a
bias-variance trade-off (checked by simulation).
"""
import matplotlib.pyplot as plt
import numpy as np

from regbook import COLORS, Generated, figure_path, use_book_style

gen = Generated("ch10", "svd_rank", prefix="svd")


# <<svd_ls>>
def svd_ls(X, y, tol=None):
    """Minimum-norm least squares solution, treating singular values <= tol as zero."""
    U, s, Vt = np.linalg.svd(X, full_matrices=False)
    if tol is None:
        tol = max(X.shape) * np.finfo(float).eps * s[0]   # the usual default
    r = int(np.sum(s > tol))                               # numerical rank
    c = U[:, :r].T @ y                                     # coordinates of the fit
    b = Vt[:r].T @ (c / s[:r])
    sse = y @ y - c @ c
    return b, r, sse, s
# <</svd_ls>>


# ---- a rank-deficient design: SVD solution is the minimum-norm one -------------------
rng = np.random.default_rng(104)
n = 30
g = np.repeat([0, 1, 2], 10)
Xr = np.column_stack([np.ones(n), np.eye(3)[g], rng.normal(size=n)])   # intercept + 3 indicators: rank 4
yr = Xr @ np.array([1.0, 0.5, -0.5, 0.0, 2.0]) + rng.normal(scale=0.3, size=n)
b_svd, r_svd, sse_svd, s_r = svd_ls(Xr, yr)
assert r_svd == 4
assert np.allclose(b_svd, np.linalg.pinv(Xr) @ yr)
assert np.allclose(b_svd, np.linalg.lstsq(Xr, yr, rcond=None)[0])
null = np.array([1.0, -1.0, -1.0, -1.0, 0.0])                      # X null = 0
assert np.allclose(Xr @ null, 0)
assert abs(null @ b_svd) < 1e-10                                  # b_svd is orthogonal to N(X)
other = b_svd + 0.7 * null
assert np.allclose(Xr @ other, Xr @ b_svd) and np.linalg.norm(other) > np.linalg.norm(b_svd)
assert np.isclose(sse_svd, np.sum((yr - Xr @ b_svd) ** 2))
assert s_r[-1] < 1e-14 * s_r[0]
gen.num("rd_smin", s_r[-1] / s_r[0], 0, sci=True)

# ---- the pseudoinverse is discontinuous ---------------------------------------------------
# <<discontinuity>>
x1 = np.array([1.0, 2.0, 3.0, 4.0])
d = np.array([1.0, -1.0, -1.0, 1.0])            # orthogonal to x1
yd = np.array([1.0, 3.0, 2.0, 5.0])
for t in [1e-2, 1e-6, 1e-10, 0.0]:
    Xt = np.column_stack([x1, x1 + t * d])        # second column tends to the first
    b = np.linalg.pinv(Xt, rcond=1e-15) @ yd
    print(f"t = {t:7.0e}   b = {np.round(b, 3)}   fit = {np.round(Xt @ b, 3)}")
# <</discontinuity>>

b0 = np.linalg.pinv(np.column_stack([x1, x1])) @ yd
b6 = np.linalg.pinv(np.column_stack([x1, x1 + 1e-6 * d]), rcond=1e-15) @ yd
assert np.allclose(b0[0], b0[1])
assert np.abs(b6).max() > 1e5                     # coefficients of order 1/t
gen.num("disc_b0", b0[0], 3)
gen.num("disc_b6", b6[0], 0, sci=False)
gen.num("disc_b6_2", b6[1], 0, sci=False)
fit_lim = np.column_stack([x1, x1]) @ b0
fit6 = np.column_stack([x1, x1 + 1e-6 * d]) @ b6
assert np.linalg.norm(fit6 - fit_lim) > 0.5      # fitted values jump too: the limit drops the d-direction

# ---- Weyl: singular values move by at most ||E|| -----------------------------------------
worst = 0.0
for _ in range(200):
    A = rng.normal(size=(12, 5))
    E = rng.normal(size=(12, 5)) * 10.0 ** rng.uniform(-8, 0)
    ds = np.abs(np.linalg.svd(A + E, compute_uv=False) - np.linalg.svd(A, compute_uv=False))
    worst = max(worst, ds.max() / np.linalg.norm(E, 2))
assert worst <= 1 + 1e-12
gen.num("weyl_worst", worst, 3)

# ---- R-SVD: an SVD of R gives an SVD of X ---------------------------------------------
Xb = rng.normal(size=(500, 6))
Q, R = np.linalg.qr(Xb)
UR, sR, VtR = np.linalg.svd(R)
assert np.allclose(np.linalg.svd(Xb, compute_uv=False), sR)
assert np.allclose(Q @ UR @ np.diag(sR) @ VtR, Xb)

# ---- truncated SVD: bias, variance, mean squared error ------------------------------------
# <<tsvd_setup>>
n, p, sigma = 60, 8, 0.1
t = np.linspace(0, 1, n)
X = np.vander(t, p, increasing=True)            # 1, t, ..., t^7: badly conditioned
U, s, Vt = np.linalg.svd(X, full_matrices=False)
a = np.array([2.0, -1.5, 1.0, -0.6, 0.3, -0.1, 0.05, -0.02])   # decaying components v_i' beta
beta = Vt.T @ a
print("condition number:", s[0] / s[-1])
print("singular values:", np.array2string(s, precision=2))
# <</tsvd_setup>>

# <<tsvd_mse>>
bias2 = np.array([np.sum((Vt[k:] @ beta) ** 2) for k in range(1, p + 1)])   # ||(I - V_k V_k') beta||^2
var = np.array([sigma ** 2 * np.sum(1 / s[:k] ** 2) for k in range(1, p + 1)])
mse = bias2 + var
print("k    bias^2      variance     MSE")
for k in range(p):
    print(f"{k + 1}  {bias2[k]:10.3e}  {var[k]:10.3e}  {mse[k]:10.3e}")
# <</tsvd_mse>>

# simulation check of the MSE formula
reps = 4000
Y = X @ beta + sigma * rng.normal(size=(reps, n))
C = Y @ U                                            # rows: U'y for each replicate
sim = []
for k in range(1, p + 1):
    B = (C[:, :k] / s[:k]) @ Vt[:k]                  # truncated SVD estimates
    sim.append(np.mean(np.sum((B - beta) ** 2, axis=1)))
sim = np.array(sim)
assert np.allclose(sim, mse, rtol=0.08)
kbest = int(np.argmin(mse)) + 1
assert kbest < p                                     # truncation beats least squares here
gen.num("kappa_poly", s[0] / s[-1], 1, sci=True)
gen.int("kbest", kbest)
gen.num("mse_best", mse[kbest - 1], 2)
gen.num("mse_ls", mse[-1], 0)
gen.num("mse_ratio", mse[-1] / mse[kbest - 1], 0)
gen.num("smin", s[-1], 1, sci=True)
gen.num("bias_best", bias2[kbest - 1], 2)
gen.num("var_best", var[kbest - 1], 2)
gen.int("reps", reps)
gen.write()

use_book_style()
fig, axes = plt.subplots(1, 2, figsize=(5.8, 2.5))
ax = axes[0]
ax.semilogy(np.arange(1, p + 1), s, "o-", ms=3, color=COLORS["accent"])
ax.axhline(max(X.shape) * np.finfo(float).eps * s[0], color=COLORS["muted"], ls=":", lw=0.8)
ax.axhline(sigma, color=COLORS["second"], ls="--", lw=0.8)
ax.text(p, sigma * 1.6, r"noise level $\sigma$", ha="right", fontsize=7, color=COLORS["second"])
ax.text(p, max(X.shape) * np.finfo(float).eps * s[0] * 2, "default tolerance", ha="right", fontsize=7, color=COLORS["muted"])
ax.set_xlabel(r"index $i$")
ax.set_ylabel(r"singular value $\sigma_i$")
ax.set_title("(a) singular values")
ax = axes[1]
k = np.arange(1, p + 1)
ax.semilogy(k[:-1], bias2[:-1], "s-", ms=3, color=COLORS["third"], label=r"bias$^2$")
ax.semilogy(k, var, "^-", ms=3, color=COLORS["second"], label="variance")
ax.semilogy(k, mse, "o-", ms=3, color=COLORS["ink"], label="MSE")
ax.semilogy(k, sim, "x", ms=5, color=COLORS["accent"], label="simulated MSE")
ax.set_xlabel(r"number of components $k$")
ax.set_title("(b) truncated SVD estimate")
ax.legend(frameon=False, fontsize=7, loc="upper left")
fig.tight_layout()
fig.savefig(figure_path("ch10", "truncated_svd"))
