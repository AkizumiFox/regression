"""Chapter 10, Section 2: the Cholesky factorization, the augmented cross-product matrix,
the sweep operator, and accurate formation of centred cross-products.

Data: the 2009 US state data shipped with statsmodels (statsmodels.datasets.statecrime,
public domain), as in Chapter 6: murder rate on poverty, single-parent households and
urbanization, District of Columbia removed.
"""
import matplotlib.pyplot as plt
import numpy as np
import statsmodels.api as sm
from fractions import Fraction

from regbook import COLORS, Generated, figure_path, use_book_style

gen = Generated("ch10", "cholesky", prefix="chol")


# <<cholesky>>
def cholesky_upper(A):
    """Upper triangular R with positive diagonal and R'R = A (A positive definite)."""
    A = np.array(A, dtype=float)
    p = A.shape[0]
    R = np.zeros_like(A)
    for j in range(p):
        d = A[j, j] - R[:j, j] @ R[:j, j]
        if d <= 0:
            raise np.linalg.LinAlgError(f"not positive definite at step {j + 1}")
        R[j, j] = np.sqrt(d)
        R[j, j + 1:] = (A[j, j + 1:] - R[:j, j] @ R[:j, j + 1:]) / R[j, j]
    return R


def solve_upper(R, b):
    """Back substitution for R x = b, R upper triangular."""
    x = np.zeros(np.shape(b))
    for i in range(len(b) - 1, -1, -1):
        x[i] = (b[i] - R[i, i + 1:] @ x[i + 1:]) / R[i, i]
    return x


def solve_lower(L, b):
    """Forward substitution for L x = b, L lower triangular."""
    x = np.zeros(np.shape(b))
    for i in range(len(b)):
        x[i] = (b[i] - L[i, :i] @ x[:i]) / L[i, i]
    return x
# <</cholesky>>


# <<data>>
data = sm.datasets.statecrime.load_pandas().data.drop(index="District of Columbia")
y = data["murder"].to_numpy()
X = np.column_stack([np.ones(len(y)), data[["poverty", "single", "urban"]]])
n, p = X.shape
# <</data>>

# <<augmented>>
Z = np.column_stack([X, y])              # augmented matrix [X, y]
S = Z.T @ Z                              # (p+1) x (p+1) cross-product matrix
T = cholesky_upper(S)                    # T = [[R, z], [0, d]]
R, z, d = T[:p, :p], T[:p, p], T[p, p]

beta = solve_upper(R, z)                 # R b = z
sse = d ** 2                             # the last diagonal entry is sqrt(SSE)
K = solve_upper(R, np.eye(p))            # R^{-1}
XtX_inv = K @ K.T                        # (X'X)^{-1} = R^{-1} R^{-T}
h = np.array([np.sum(solve_lower(R.T, x) ** 2) for x in X])   # leverages
print("coefficients:", np.round(beta, 4))
print(f"SSE = {sse:.4f},  d = {d:.4f}")
# <</augmented>>

fit = sm.OLS(y, X).fit()
assert np.allclose(cholesky_upper(S), np.linalg.cholesky(S).T, rtol=1e-12, atol=1e-12)
assert np.allclose(beta, fit.params, rtol=1e-9)
assert np.isclose(sse, fit.ssr, rtol=1e-9)
assert np.allclose(XtX_inv, np.linalg.inv(X.T @ X), rtol=1e-8)
assert np.allclose(h, fit.get_influence().hat_matrix_diag, rtol=1e-8)
assert np.isclose(np.prod(np.diag(R)) ** 2, np.linalg.det(X.T @ X), rtol=1e-9)
se = np.sqrt(sse / (n - p) * np.diag(XtX_inv))
assert np.allclose(se, fit.bse, rtol=1e-8)
qr_beta = np.linalg.lstsq(X, y, rcond=None)[0]
rel_ne_qr = np.linalg.norm(beta - qr_beta) / np.linalg.norm(qr_beta)
assert rel_ne_qr < 1e-10
gen.int("n", n)
gen.int("p", p)
for k, b in enumerate(beta):
    gen.num(f"b{k}", b, 4)
gen.num("sse", sse, 2)
gen.num("d", d, 4)
gen.num("kappa", np.linalg.cond(X), 0)
gen.num("rel_ne_qr", rel_ne_qr, 0, sci=True)
gen.num("hmax", h.max(), 3)
gen.num("det_log10", np.log10(np.prod(np.diag(R)) ** 2), 2)


# ---- the sweep operator -----------------------------------------------------------------
# <<sweep>>
def sweep(A, k):
    """Sweep the symmetric matrix A on pivot k."""
    A = np.array(A, dtype=float)
    a = A[k, k]
    B = A - np.outer(A[:, k], A[k, :]) / a
    B[k, :] = A[k, :] / a
    B[:, k] = A[:, k] / a
    B[k, k] = -1 / a
    return B


W = S.copy()
for k in range(p):                       # sweep the p columns of X, one at a time
    W = sweep(W, k)
print("coefficients from the sweep:", np.round(W[:p, p], 4))
print(f"SSE from the sweep: {W[p, p]:.4f}")
# <</sweep>>

assert np.allclose(W[:p, p], beta, rtol=1e-9)
assert np.isclose(W[p, p], sse, rtol=1e-9)
assert np.allclose(-W[:p, :p], XtX_inv, rtol=1e-8)
# sweeps on different pivots commute
assert np.allclose(sweep(sweep(S, 0), 2), sweep(sweep(S, 2), 0))
# after sweeping pivots 0 and 2 the corner holds the SSE of the submodel [1, single]
W02 = sweep(sweep(S, 0), 2)
sub = sm.OLS(y, X[:, [0, 2]]).fit()
assert np.isclose(W02[p, p], sub.ssr, rtol=1e-9)
gen.num("sse_sub", sub.ssr, 2)


# ---- centred cross-products: textbook, two-pass and updating (Welford) -----------------
# <<centring>>
def ss_textbook(x):
    return np.sum(x * x) - len(x) * np.mean(x) ** 2        # one pass, then subtract


def ss_two_pass(x):
    m = np.mean(x)                                         # first pass: the mean
    return np.sum((x - m) ** 2)                            # second pass: deviations


def ss_updating(x):
    m, s = 0.0, 0.0
    for k, xk in enumerate(x, start=1):                    # one pass, no cancellation
        d = xk - m
        m = m + d / k
        s = s + (k - 1) / k * d * d
    return s
# <</centring>>


def ss_exact(x):
    xs = [Fraction(v) for v in x]
    m = sum(xs) / len(xs)
    return sum((v - m) ** 2 for v in xs)


rng = np.random.default_rng(102)
zc = rng.normal(size=1000)
offsets = 10.0 ** np.arange(0, 9)
rows = []
for c in offsets:
    x = c + zc                                             # sd about 1, mean c
    exact = ss_exact(x)
    errs = [abs(float((Fraction(f(x)) - exact) / exact)) for f in (ss_textbook, ss_two_pass, ss_updating)]
    kappa_col = np.linalg.norm(x) / np.sqrt(float(exact))   # uncentred / centred length
    rows.append((c, kappa_col, *errs))
rows = np.array(rows)
u = np.finfo(float).eps / 2
# textbook error grows like kappa^2 u; the other two stay near u (times modest factors)
assert np.all(rows[:, 2] <= 10 * len(zc) * rows[:, 1] ** 2 * u)
assert rows[-1, 2] > 1e-2 and rows[4, 2] > 1e-9
# two-pass stays near u; updating grows like kappa u; textbook like kappa^2 u
assert np.all(rows[:, 3] < 1e-14)
assert np.all(rows[:, 4] <= len(zc) * rows[:, 1] * u)
assert np.all(rows[4:, 4] < rows[4:, 2] / 100)

# the vector form of prp-cmp-welford: running means and centred cross-product matrices
W = rng.normal(size=(200, 4)) @ rng.normal(size=(4, 4)) + 50.0
m_k, C_k = np.zeros(4), np.zeros((4, 4))
for k, wk in enumerate(W, start=1):
    d = wk - m_k
    m_k = m_k + d / k
    C_k = C_k + (k - 1) / k * np.outer(d, d)
    if k >= 2:
        assert np.allclose(m_k, W[:k].mean(axis=0), rtol=1e-12, atol=1e-10)
        assert np.allclose(C_k, np.cov(W[:k].T) * (k - 1), rtol=1e-10, atol=1e-8)

r6 = rows[offsets == 1e6][0]
gen.num("kc6", r6[1], 1, sci=True)
gen.num("err_text6", r6[2], 1, sci=True)
gen.num("err_two6", r6[3], 1, sci=True)
gen.num("err_upd6", r6[4], 1, sci=True)
r8 = rows[offsets == 1e8][0]
gen.num("err_text8", r8[2], 1)
gen.num("err_upd8", r8[4], 1, sci=True)
gen.num("err_two8", r8[3], 1, sci=True)
gen.write()

use_book_style()
fig, ax = plt.subplots(figsize=(3.9, 2.7))
ax.loglog(rows[:, 1], np.maximum(rows[:, 2], 1e-18), "o-", ms=3, color=COLORS["second"], label="textbook (one pass)")
ax.loglog(rows[:, 1], np.maximum(rows[:, 3], 1e-18), "s-", ms=3, color=COLORS["accent"], label="two pass")
ax.loglog(rows[:, 1], np.maximum(rows[:, 4], 1e-18), "^-", ms=3, color=COLORS["third"], label="updating (one pass)")
kk = np.logspace(0, np.log10(rows[:, 1].max()), 50)
ax.loglog(kk, kk ** 2 * u, ":", color=COLORS["second"], lw=0.8, label=r"$\kappa^2 u$")
ax.loglog(kk, kk * u, ":", color=COLORS["third"], lw=0.8, label=r"$\kappa u$")
ax.set_xlabel(r"column condition $\|\mathbf{x}\|/\|\mathbf{x}-\bar{x}\mathbf{1}\|$")
ax.set_ylabel("relative error of centred SS")
ax.set_ylim(1e-18, 1e2)
ax.legend(frameon=False, loc="upper left")
fig.savefig(figure_path("ch10", "centring_accuracy"))
