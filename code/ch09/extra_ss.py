"""Chapter 9, Section 2: extra sums of squares on the 2009 US state data.

Adding regressors to a fitted model, the sum of squares of a linear hypothesis,
restricted least squares, and the extra sums of squares SS(x | S) for every subset S.
Public-domain data shipped with statsmodels (statsmodels.datasets.statecrime).
"""
from itertools import combinations

import numpy as np
import statsmodels.api as sm

from regbook import Generated

data = sm.datasets.statecrime.load_pandas().data.drop(index="District of Columbia")
y = data["murder"].to_numpy()
n = len(y)
one = np.ones(n)
cols = {"poverty": data["poverty"].to_numpy(), "single": data["single"].to_numpy(),
        "urban": data["urban"].to_numpy()}


def sse(Z):
    b, *_ = np.linalg.lstsq(Z, y, rcond=None)
    e = y - Z @ b
    return e @ e


# <<adding>>
X = np.column_stack([one, cols["single"]])                 # the model already fitted
Z = np.column_stack([cols["poverty"], cols["urban"]])       # the regressors to add
XtX_inv = np.linalg.inv(X.T @ X)
beta_X = XtX_inv @ X.T @ y
R = np.eye(n) - X @ XtX_inv @ X.T                           # residual projection of the old model

gamma = np.linalg.solve(Z.T @ R @ Z, Z.T @ R @ y)          # new coefficients
beta_new = beta_X - XtX_inv @ X.T @ Z @ gamma               # old coefficients, corrected
extra_ss = gamma @ Z.T @ R @ y                              # SS(Z | X) = SSE_X - SSE_W
print("gamma =", gamma.round(4), " corrected beta =", beta_new.round(4))
print(f"SSE_X = {y @ R @ y:.2f},  extra SS = {extra_ss:.2f}")
# <</adding>>

W = np.column_stack([X, Z])
fit_W = sm.OLS(y, W).fit()
assert np.allclose(fit_W.params, np.concatenate([beta_new, gamma]))
assert np.isclose(extra_ss, sse(X) - sse(W))
# covariance formulas of the adding-regressors theorem (sigma^2 = 1)
L = XtX_inv @ X.T @ Z
Minv = np.linalg.inv(Z.T @ R @ Z)
cov = np.linalg.inv(W.T @ W)
assert np.allclose(cov[2:, 2:], Minv)
assert np.allclose(cov[:2, :2], XtX_inv + L @ Minv @ L.T)
assert np.allclose(cov[:2, 2:], -L @ Minv)

# ---- a single column: partial SS = beta_j^2 / c_jj = t^2 s^2 -----------------
# <<single>>
Xf = np.column_stack([one, cols["poverty"], cols["single"], cols["urban"]])
fit = sm.OLS(y, Xf).fit()
C = np.linalg.inv(Xf.T @ Xf)
j = 1                                                       # poverty
partial_ss = fit.params[j] ** 2 / C[j, j]
print(f"beta = {fit.params[j]:.4f}, c_jj = {C[j, j]:.6f}, beta^2/c_jj = {partial_ss:.2f}")
print(f"t = {fit.tvalues[j]:.3f},  t^2 * s^2 = {fit.tvalues[j] ** 2 * fit.scale:.2f}")
# <</single>>
X_drop = np.delete(Xf, j, axis=1)
assert np.isclose(partial_ss, sse(X_drop) - sse(Xf))
assert np.isclose(partial_ss, fit.tvalues[j] ** 2 * fit.scale)


# ---- the sum of squares of a linear hypothesis ------------------------------
# <<hypothesis>>
def hypothesis_ss(X, y, Lam, d=None):
    """(L'b - d)' [L' (X'X)^- L]^{-1} (L'b - d) for estimable L'beta = d."""
    G = np.linalg.pinv(X.T @ X)
    b = G @ X.T @ y
    u = Lam.T @ b - (0 if d is None else d)
    return u @ np.linalg.solve(Lam.T @ G @ Lam, u)


Lam_equal = np.array([[0.0, 1.0, -1.0, 0.0]]).T             # beta_poverty = beta_single
ss_equal = hypothesis_ss(Xf, y, Lam_equal)
X_equal = np.column_stack([one, cols["poverty"] + cols["single"], cols["urban"]])
print(f"H: equal slopes   formula {ss_equal:.3f}   SSE_0 - SSE {sse(X_equal) - sse(Xf):.3f}")
# <</hypothesis>>
assert np.isclose(ss_equal, sse(X_equal) - sse(Xf))
Lam_two = np.zeros((4, 2))
Lam_two[1, 0] = Lam_two[3, 1] = 1.0                         # poverty and urban both zero
ss_two = hypothesis_ss(Xf, y, Lam_two)
assert np.isclose(ss_two, sse(np.column_stack([one, cols["single"]])) - sse(Xf))
assert np.isclose(ss_two, extra_ss)

# restricted least squares with a nonzero right-hand side: beta_single = 0.5
lam = np.array([0.0, 0.0, 1.0, 0.0])
d = 0.5
b = fit.params
b_H = b - C @ lam * (lam @ b - d) / (lam @ C @ lam)
assert np.isclose(lam @ b_H, d)
b_direct, *_ = np.linalg.lstsq(np.delete(Xf, 2, axis=1), y - d * cols["single"], rcond=None)
assert np.allclose(np.delete(b_H, 2), b_direct)
e_H = y - Xf @ b_H
ss_restricted = e_H @ e_H - fit.ssr
assert np.isclose(ss_restricted, hypothesis_ss(Xf, y, lam[:, None], np.array([d])))
assert np.isclose(ss_restricted, (Xf @ (b - b_H)) @ (Xf @ (b - b_H)))

# ---- SS(x | S) for every subset S of the other two regressors ----------------
# <<subsets>>
table = {}
for name, x in cols.items():
    others = [k for k in cols if k != name]
    for size in range(3):
        for S in combinations(others, size):
            base = np.column_stack([one] + [cols[k] for k in S])
            table[name, S] = sse(base) - sse(np.column_stack([base, x]))
            print(f"SS({name} | 1{''.join(', ' + k for k in S)}) = {table[name, S]:7.2f}")
# <</subsets>>

gen = Generated("ch09", "extra_ss", prefix="xss")
for k, v in zip(["b0", "b1"], beta_X):
    gen.num(f"betaX:{k}", v, 4)
for k, v in zip(["b0", "b1"], beta_new):
    gen.num(f"betanew:{k}", v, 4)
gen.num("gamma:pov", gamma[0], 4)
gen.num("gamma:urb", gamma[1], 4)
gen.num("sseX", sse(X), 2)
gen.num("sseW", sse(W), 2)
gen.num("extra", extra_ss, 2)
gen.num("beta_pov", fit.params[1], 4)
gen.num("cjj", C[1, 1], 6)
gen.num("partial", partial_ss, 2)
gen.num("t", fit.tvalues[1], 3)
gen.num("t2", fit.tvalues[1] ** 2, 3)
gen.num("s2", fit.scale, 3)
gen.num("beta_single", fit.params[2], 4)
gen.num("ss_equal", ss_equal, 3)
gen.num("F_equal", ss_equal / fit.scale, 3)
gen.num("ss_restricted", ss_restricted, 3)
gen.num("b_H_pov", b_H[1], 4)
gen.num("b_H_urb", b_H[3], 4)
short = {"poverty": "pov", "single": "sin", "urban": "urb"}
for (name, S), v in table.items():
    gen.num(f"R:{short[name]}:{'-'.join(short[k] for k in S) or 'none'}", v, 2)
gen.write()
