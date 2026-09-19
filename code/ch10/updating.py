"""Chapter 10, Section 6: updating and downdating.

(1) Deleting an observation from a triangular factor (LINPACK-style downdating), checked
    against refitting; all leave-one-out quantities from one fit.
(2) Adding and deleting a column of a QR factorization.
(3) A moving-window regression maintained three ways: rank-one updates of (X'X)^{-1},
    Givens updating and downdating of the triangle, and refitting from scratch.

Data: statsmodels.datasets.statecrime and statsmodels.datasets.macrodata (both public domain).
"""
import matplotlib.pyplot as plt
import numpy as np
import statsmodels.api as sm

from regbook import COLORS, Generated, figure_path, use_book_style

gen = Generated("ch10", "updating", prefix="upd")


# <<rotations>>
def givens(a, b):
    """c, s, r with [[c, s], [-s, c]] @ [a, b] = [r, 0]."""
    if b == 0.0:
        return 1.0, 0.0, a
    r = np.hypot(a, b)
    return a / r, b / r, r


def update(T, z):
    """Triangle of [A; z'] from the triangle T of A: rotate the new row into T."""
    T, z = T.copy(), np.array(z, dtype=float)
    for k in range(len(z)):
        c, s, r = givens(T[k, k], z[k])
        Tk, zk = T[k, k:].copy(), z[k:].copy()
        T[k, k:], z[k:] = c * Tk + s * zk, -s * Tk + c * zk
    return T


def downdate(T, z):
    """Triangle of A with the row z' removed, from the triangle T of A (T'T = A'A)."""
    m = len(z)
    a = np.linalg.solve(T.T, z)                   # T'a = z; ||a||^2 is the leverage of z
    alpha2 = 1.0 - a @ a
    if alpha2 <= 0:
        raise np.linalg.LinAlgError("the row cannot be removed (leverage >= 1)")
    alpha = np.sqrt(alpha2)
    T, w = T.copy(), np.zeros(m)                  # w: the extra row, initially zero
    for k in range(m - 1, -1, -1):                # rotate (a_k, alpha) into (0, new alpha)
        c, s, r = givens(alpha, a[k])
        alpha = r
        Tk = T[k, k:].copy()
        T[k, k:], w[k:] = c * Tk - s * w[k:], s * Tk + c * w[k:]
    return T                                      # on exit w equals z
# <</rotations>>


def to_positive(T):
    s = np.sign(np.diag(T))
    s[s == 0] = 1
    return s[:, None] * T


# ---- (1) deletion: downdating and leave-one-out quantities --------------------------------
# <<loo>>
data = sm.datasets.statecrime.load_pandas().data           # all 50 states and DC
y = data["murder"].to_numpy()
X = np.column_stack([np.ones(len(y)), data[["poverty", "single", "urban"]]])
n, p = X.shape

Q, R = np.linalg.qr(X)
beta = np.linalg.solve(R, Q.T @ y)
e = y - X @ beta
sse = e @ e
h = np.sum(Q ** 2, axis=1)                          # leverages
press_resid = e / (1 - h)                           # y_i minus its prediction without case i
sse_loo = sse - e ** 2 / (1 - h)                    # SSE after deleting case i
s_loo = np.sqrt(sse_loo / (n - 1 - p))
t_ext = e / (s_loo * np.sqrt(1 - h))                # externally studentized residuals
i = int(np.argmax(np.abs(press_resid)))
print(f"largest PRESS residual: {data.index[i]}, e = {e[i]:.3f}, h = {h[i]:.3f},"
      f" e/(1-h) = {press_resid[i]:.3f}, t = {t_ext[i]:.2f}")
print(f"PRESS = {np.sum(press_resid ** 2):.2f}, SSE = {sse:.2f}")
# <</loo>>

# check against brute-force refitting and against downdating the augmented triangle
Z = np.column_stack([X, y])
T = np.linalg.qr(Z, mode="r")
for j in range(n):
    keep = np.arange(n) != j
    fit_j = sm.OLS(y[keep], X[keep]).fit()
    assert np.isclose(y[j] - X[j] @ fit_j.params, press_resid[j], rtol=1e-9, atol=1e-10)
    assert np.isclose(fit_j.ssr, sse_loo[j], rtol=1e-9)
    Tj = to_positive(downdate(T, Z[j]))
    assert np.allclose(np.abs(Tj), np.abs(np.linalg.qr(Z[keep], mode="r")), rtol=1e-8, atol=1e-8)
    assert np.isclose(Tj[p, p] ** 2, fit_j.ssr, rtol=1e-8)                 # corner = sqrt(SSE_(j))
    assert np.allclose(np.linalg.solve(Tj[:p, :p], Tj[:p, p]), fit_j.params, rtol=1e-8)
infl = sm.OLS(y, X).fit().get_influence()
assert np.allclose(t_ext, infl.resid_studentized_external, rtol=1e-9)
assert np.allclose(to_positive(update(downdate(T, Z[5]), Z[5])), to_positive(T), atol=1e-9)
assert data.index[i] == "District of Columbia"
gen.text("loo_name", data.index[i])
gen.num("loo_e", e[i], 3)
gen.num("loo_h", h[i], 3)
gen.num("loo_press", press_resid[i], 3)
gen.num("loo_t", t_ext[i], 2)
gen.num("loo_sse", sse, 2)
gen.num("loo_sse_i", sse_loo[i], 2)
gen.num("press", np.sum(press_resid ** 2), 2)
gen.num("press_ratio", np.sum(press_resid ** 2) / sse, 2)
gen.num("press_balanced", (n / (n - p)) ** 2, 2)
gen.int("n", n)


# ---- (2) adding and deleting columns --------------------------------------------------------
# <<add_column>>
def add_column(Q1, R, x):
    """Thin QR of [X, x] from the thin QR of X (one Gram-Schmidt step, repeated once)."""
    r = Q1.T @ x
    w = x - Q1 @ r
    r2 = Q1.T @ w                                 # reorthogonalize: guards against cancellation
    w, r = w - Q1 @ r2, r + r2
    rho = np.linalg.norm(w)
    Rn = np.block([[R, r[:, None]], [np.zeros((1, R.shape[1])), np.array([[rho]])]])
    return np.column_stack([Q1, w / rho]), Rn
# <</add_column>>


def delete_column(R, j):
    """Triangle of X with column j removed, from the triangle R of X: restore with rotations."""
    H = np.delete(R, j, axis=1)                   # upper Hessenberg from column j on
    for k in range(j, H.shape[1]):
        c, s, r = givens(H[k, k], H[k + 1, k])
        top, bot = H[k, k:].copy(), H[k + 1, k:].copy()
        H[k, k:], H[k + 1, k:] = c * top + s * bot, -s * top + c * bot
    return H[:-1]


keep_rows = data.index != "District of Columbia"
yk, Xk = y[keep_rows], X[keep_rows]
Q1, R1 = np.linalg.qr(Xk[:, :3])                  # intercept, poverty, single
Q2, R2 = add_column(Q1, R1, Xk[:, 3])             # add urban
assert np.allclose(Q2 @ R2, Xk) and np.allclose(Q2.T @ Q2, np.eye(4), atol=1e-13)
drop = (Q2[:, 3] @ yk) ** 2                       # SSE falls by the square of the new coordinate
assert np.isclose(drop, sm.OLS(yk, Xk[:, :3]).fit().ssr - sm.OLS(yk, Xk).fit().ssr, rtol=1e-8)
Rd = delete_column(np.linalg.qr(Xk, mode="r"), 1)  # delete poverty
assert np.allclose(np.abs(Rd), np.abs(np.linalg.qr(Xk[:, [0, 2, 3]], mode="r")), rtol=1e-10, atol=1e-10)

# ---- (3) moving-window regression -----------------------------------------------------------
# <<window_setup>>
macro = sm.datasets.macrodata.load_pandas().data
yw = macro["realcons"].to_numpy()
Xw = np.column_stack([np.ones(len(yw)), macro[["realdpi", "realinv", "pop", "unemp"]]])
w = 40                                            # ten years of quarterly data
print("kappa of the full design:", f"{np.linalg.cond(Xw):.2e}")
# <</window_setup>>

# <<window>>
K = np.linalg.inv(Xw[:w].T @ Xw[:w])              # (a) rank-one updates of the inverse
g = Xw[:w].T @ yw[:w]
Zw = np.column_stack([Xw, yw])
Tw = np.linalg.qr(Zw[:w], mode="r")               # (b) triangle of [X, y]
p1 = Xw.shape[1]
err_sm, err_giv = [], []
for t in range(w, len(yw)):
    new, old = Xw[t], Xw[t - w]
    Kx = K @ new                                  # add the new quarter (Sherman-Morrison)
    K = K - np.outer(Kx, Kx) / (1 + new @ Kx)
    Kx = K @ old                                  # delete the oldest quarter
    K = K + np.outer(Kx, Kx) / (1 - old @ Kx)
    g = g + new * yw[t] - old * yw[t - w]
    b_sm = K @ g
    Tw = downdate(update(Tw, Zw[t]), Zw[t - w])
    b_giv = np.linalg.solve(Tw[:p1, :p1], Tw[:p1, p1])
    b_ref = np.linalg.lstsq(Xw[t - w + 1:t + 1], yw[t - w + 1:t + 1], rcond=None)[0]   # (c) refit
    err_sm.append(np.linalg.norm(b_sm - b_ref) / np.linalg.norm(b_ref))
    err_giv.append(np.linalg.norm(b_giv - b_ref) / np.linalg.norm(b_ref))
print(f"final relative error: inverse updating {err_sm[-1]:.1e}, Givens {err_giv[-1]:.1e}")
# <</window>>

err_sm, err_giv = np.array(err_sm), np.array(err_giv)
steps = len(err_sm)
kw = np.max([np.linalg.cond(Xw[t - w + 1:t + 1]) for t in range(w, len(yw))])
assert err_giv.max() < 1e-6
assert np.median(err_sm) > 100 * np.median(err_giv)
gen.int("steps", steps)
gen.num("kw_max", kw, 1, sci=True)
gen.num("err_sm_final", err_sm[-1], 1, sci=True)
gen.num("err_giv_final", err_giv[-1], 1, sci=True)
gen.num("err_sm_max", err_sm.max(), 1, sci=True)
gen.num("err_giv_max", err_giv.max(), 1, sci=True)
gen.write()

use_book_style()
fig, ax = plt.subplots(figsize=(3.9, 2.6))
idx = np.arange(1, steps + 1)
ax.semilogy(idx, err_sm, color=COLORS["second"], label=r"updating $(\mathbf{X}^\top\mathbf{X})^{-1}$")
ax.semilogy(idx, err_giv, color=COLORS["accent"], label="updating the triangle")
ax.set_xlabel("window step (quarters)")
ax.set_ylabel("relative error vs. refitting")
ax.set_ylim(1e-17, 1e0)
ax.legend(frameon=False, loc="upper left")
fig.savefig(figure_path("ch10", "moving_window"))
