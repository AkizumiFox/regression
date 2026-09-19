"""Chapter 10, Section 3: Householder and Givens QR, implemented directly and checked
against LAPACK (numpy), with the regression quantities read off the factorization, and
the loss of orthogonality of classical and modified Gram-Schmidt.

Data for the worked example: statsmodels.datasets.statecrime (public domain), as in
Chapter 6 and Section 10.2.
"""
import matplotlib.pyplot as plt
import numpy as np
import statsmodels.api as sm

from regbook import COLORS, Generated, figure_path, use_book_style

gen = Generated("ch10", "householder_givens", prefix="hg")


# <<householder>>
def house(x):
    """Unit vector v with (I - 2 v v') x = alpha e_1, alpha = -sign(x_1) ||x||."""
    alpha = -np.copysign(np.linalg.norm(x), x[0])
    v = x.astype(float).copy()
    v[0] -= alpha                        # x_1 - alpha adds two numbers of the same sign
    return v / np.linalg.norm(v), alpha


def householder_qr(X):
    """Return the reflection vectors and the p x p triangle R of X = Q [R; 0]."""
    A = np.array(X, dtype=float)
    n, p = A.shape
    vs = []
    for j in range(p):
        v, alpha = house(A[j:, j])
        A[j:, j:] -= 2.0 * np.outer(v, v @ A[j:, j:])    # apply I - 2vv' to the trailing block
        vs.append(v)
    return vs, np.triu(A[:p, :])


def apply_Qt(vs, y):
    """Q'y, applying the reflections H_1, ..., H_p in turn (Q itself is never formed)."""
    z = np.array(y, dtype=float)
    for j, v in enumerate(vs):
        z[j:] -= 2.0 * v * (v @ z[j:])
    return z


def apply_Q(vs, z):
    """Q z, applying the reflections in reverse order."""
    w = np.array(z, dtype=float)
    for j in range(len(vs) - 1, -1, -1):
        v = vs[j]
        w[j:] -= 2.0 * v * (v @ w[j:])
    return w
# <</householder>>


# <<givens>>
def givens(a, b):
    """c, s, r with [[c, s], [-s, c]] @ [a, b] = [r, 0]."""
    if b == 0.0:
        return 1.0, 0.0, a
    r = np.hypot(a, b)                   # sqrt(a^2 + b^2) without overflow
    return a / r, b / r, r


def givens_qr_R(X):
    """R of X = Q [R; 0], zeroing the subdiagonal column by column with plane rotations."""
    A = np.array(X, dtype=float)
    n, p = A.shape
    for j in range(p):
        for i in range(n - 1, j, -1):    # zero A[i, j] against the row above it
            c, s, r = givens(A[i - 1, j], A[i, j])
            top, bot = A[i - 1, j:].copy(), A[i, j:].copy()
            A[i - 1, j:] = c * top + s * bot
            A[i, j:] = -s * top + c * bot
    return np.triu(A[:p, :])
# <</givens>>


# <<rowwise>>
def add_row(R, z, x, y):
    """Rotate one observation (x, y) into the triangle R and the vector z = (Q'y)[:p].

    Returns the new R, z and the entry that leaves: its square is the increase in SSE.
    """
    R, z, x, y = R.copy(), z.copy(), np.array(x, dtype=float), float(y)
    for k in range(len(x)):
        c, s, r = givens(R[k, k], x[k])
        Rk, xk = R[k, k:].copy(), x[k:].copy()
        R[k, k:], x[k:] = c * Rk + s * xk, -s * Rk + c * xk
        z[k], y = c * z[k] + s * y, -s * z[k] + c * y
    return R, z, y


def rowwise_ls(X, y):
    """Least squares with one pass over the rows and O(p^2) memory."""
    p = X.shape[1]
    R, z, sse = np.zeros((p, p)), np.zeros(p), 0.0
    for xi, yi in zip(X, y):
        R, z, leftover = add_row(R, z, xi, yi)
        sse += leftover ** 2
    return R, z, sse
# <</rowwise>>


# ---- the regression quantities from the factorization (state data) -----------------------
# <<statecrime>>
data = sm.datasets.statecrime.load_pandas().data.drop(index="District of Columbia")
y = data["murder"].to_numpy()
X = np.column_stack([np.ones(len(y)), data[["poverty", "single", "urban"]]])
n, p = X.shape

vs, R = householder_qr(X)
c = apply_Qt(vs, y)                      # c = Q'y = (c1, c2)
c1, c2 = c[:p], c[p:]
beta = np.linalg.solve(R, c1)            # triangular system R b = c1
sse = c2 @ c2                            # SSE = ||c2||^2
resid = apply_Q(vs, np.concatenate([np.zeros(p), c2]))
seq_ss = c1 ** 2                         # sequential sums of squares, in column order
print("R diagonal:", np.round(np.diag(R), 4))
print("coefficients:", np.round(beta, 4), f"  SSE = {sse:.4f}")
print("sequential SS (intercept, poverty, single, urban):", np.round(seq_ss, 3))
# <</statecrime>>

fit = sm.OLS(y, X).fit()
Qnp, Rnp = np.linalg.qr(X)
# same triangle as LAPACK up to the signs of its rows
signs = np.sign(np.diag(R)) * np.sign(np.diag(Rnp))
assert np.allclose(R, signs[:, None] * Rnp, rtol=1e-12, atol=1e-12)
assert np.allclose(beta, fit.params, rtol=1e-10)
assert np.isclose(sse, fit.ssr, rtol=1e-10)
assert np.allclose(resid, fit.resid, atol=1e-10)
assert np.allclose(np.abs(givens_qr_R(X)), np.abs(R), rtol=1e-10, atol=1e-10)
R_row, z_row, sse_row = rowwise_ls(X, y)
print("coefficients (row-wise Givens):", np.round(np.linalg.solve(R_row, z_row), 4))
print(f"SSE (row-wise Givens): {sse_row:.4f}")
assert np.allclose(np.abs(R_row), np.abs(R), rtol=1e-10, atol=1e-10)
assert np.allclose(np.linalg.solve(R_row, z_row), beta, rtol=1e-10)
assert np.isclose(sse_row, sse, rtol=1e-10)
# sequential SS = drops in SSE along the nested sequence of models
ssr_chain = [y @ y] + [sm.OLS(y, X[:, :k]).fit().ssr for k in range(1, p + 1)]
assert np.allclose(seq_ss, -np.diff(ssr_chain), rtol=1e-9)
assert np.isclose(seq_ss[0], n * y.mean() ** 2)
# leverages are squared row lengths of Q1
Q1 = np.column_stack([apply_Q(vs, e) for e in np.eye(n)[:, :p].T])
assert np.allclose(Q1.T @ Q1, np.eye(p), atol=1e-13)
assert np.allclose(np.sum(Q1 ** 2, axis=1), fit.get_influence().hat_matrix_diag, rtol=1e-10)
names = ["int", "pov", "sgl", "urb"]
for k in range(p):
    gen.num(f"seq_{names[k]}", seq_ss[k], 2)
    gen.num(f"r{k}", R[k, k], 3)
gen.num("sse", sse, 2)
gen.int("n", n)

# ---- orthogonality of computed Q: CGS, MGS and Householder ----------------------------
# <<gram_schmidt>>
def cgs(X):
    """Classical Gram-Schmidt: inner products with the original column."""
    n, p = X.shape
    Q = np.zeros((n, p))
    for j in range(p):
        w = X[:, j] - Q[:, :j] @ (Q[:, :j].T @ X[:, j])
        Q[:, j] = w / np.linalg.norm(w)
    return Q


def mgs(X):
    """Modified Gram-Schmidt: subtract each direction from the working vector at once."""
    Q = np.array(X, dtype=float)
    n, p = Q.shape
    for j in range(p):
        Q[:, j] /= np.linalg.norm(Q[:, j])
        Q[:, j + 1:] -= np.outer(Q[:, j], Q[:, j] @ Q[:, j + 1:])
    return Q
# <</gram_schmidt>>


rng = np.random.default_rng(103)
n_t, p_t = 100, 12
U0, _ = np.linalg.qr(rng.normal(size=(n_t, p_t)))
V0, _ = np.linalg.qr(rng.normal(size=(p_t, p_t)))
rows = []
for logk in np.arange(1, 15.5, 1.0):
    s = np.logspace(0, -logk, p_t)                  # singular values from 1 to 10^-logk
    Xt = U0 @ np.diag(s) @ V0.T
    vs_t, _ = householder_qr(Xt)
    Qh = np.column_stack([apply_Q(vs_t, e) for e in np.eye(n_t)[:, :p_t].T])
    loss = [np.linalg.norm(Qm.T @ Qm - np.eye(p_t), 2) for Qm in (cgs(Xt), mgs(Xt), Qh)]
    rows.append((10.0 ** logk, *loss))
rows = np.array(rows)
u = np.finfo(float).eps / 2
assert np.all(rows[:, 3] < 1e-13)                            # Householder: orthogonal to working accuracy
assert np.all(rows[:, 2] < 100 * rows[:, 0] * u * p_t)       # MGS: loss of order kappa u
assert rows[rows[:, 0] == 1e8][0, 1] > 1e-2                  # CGS: far from orthogonal already at 1e8
k8 = rows[rows[:, 0] == 1e8][0]
gen.num("cgs8", k8[1], 1, sci=True)
gen.num("mgs8", k8[2], 1, sci=True)
gen.num("hh8", k8[3], 1, sci=True)
gen.write()

use_book_style()
fig, ax = plt.subplots(figsize=(3.9, 2.7))
ax.loglog(rows[:, 0], rows[:, 1], "o-", ms=3, color=COLORS["second"], label="classical Gram–Schmidt")
ax.loglog(rows[:, 0], rows[:, 2], "s-", ms=3, color=COLORS["accent"], label="modified Gram–Schmidt")
ax.loglog(rows[:, 0], rows[:, 3], "^-", ms=3, color=COLORS["third"], label="Householder")
ax.loglog(rows[:, 0], rows[:, 0] * u, ":", color=COLORS["accent"], lw=0.8, label=r"$\kappa u$")
ax.set_xlabel(r"condition number $\kappa(\mathbf{X})$")
ax.set_ylabel(r"$\|\mathbf{Q}^\top\mathbf{Q}-\mathbf{I}\|_2$")
ax.set_ylim(1e-17, 1e2)
ax.legend(frameon=False, loc="upper left")
fig.savefig(figure_path("ch10", "orthogonality_loss"))
