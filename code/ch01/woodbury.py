"""Chapter 1, Section 6: the Sherman-Morrison-Woodbury identity and the determinant lemma.

(1) Deleting one observation from a regression is a rank-one downdate of X'X. The stack
    loss data (Brownlee; public domain, statsmodels.datasets.stackloss) has n = 21.
(2) A covariance matrix 'diagonal + low rank' (a factor model) is inverted by solving a
    k x k system instead of an n x n one.
"""
import numpy as np
import statsmodels.api as sm

from regbook import Generated

data = sm.datasets.stackloss.load_pandas().data
X = np.column_stack([np.ones(len(data)), data[["AIRFLOW", "WATERTEMP", "ACIDCONC"]]])
n, p = X.shape

# <<deletion>>
A_inv = np.linalg.inv(X.T @ X)
h = np.einsum("ij,jk,ik->i", X, A_inv, X)             # leverages x_i' (X'X)^{-1} x_i

worst = 0.0
for i in range(n):
    x = X[i]
    # Sherman-Morrison with u = -x, v = x:  (X'X - x x')^{-1}
    B_update = A_inv + np.outer(A_inv @ x, x @ A_inv) / (1 - h[i])
    X_del = np.delete(X, i, axis=0)
    B_direct = np.linalg.inv(X_del.T @ X_del)
    worst = max(worst, np.abs(B_update - B_direct).max() / np.abs(B_direct).max())
print(f"largest relative discrepancy over all {n} deletions: {worst:.1e}")
# <</deletion>>
assert worst < 1e-10
# determinant lemma for the deletion: det(X'X - x x') = det(X'X) (1 - h_i)
i_max = int(np.argmax(h))
X_del = np.delete(X, i_max, axis=0)
ratio = np.linalg.det(X_del.T @ X_del) / np.linalg.det(X.T @ X)
assert np.isclose(ratio, 1 - h[i_max])

# <<factor>>
rng = np.random.default_rng(2024)
N, k = 1500, 4
L = rng.normal(size=(N, k))                          # loadings
d = rng.uniform(0.5, 2.0, size=N)                    # specific variances
Sigma = np.diag(d) + L @ L.T
v = rng.normal(size=N)

# Woodbury: (D + L L')^{-1} = D^{-1} - D^{-1} L (I + L' D^{-1} L)^{-1} L' D^{-1}
DinvL = L / d[:, None]
core = np.eye(k) + L.T @ DinvL                       # only k x k
w_woodbury = v / d - DinvL @ np.linalg.solve(core, DinvL.T @ v)
w_direct = np.linalg.solve(Sigma, v)                 # N x N solve

# determinant lemma: log det(D + L L') = log det D + log det(I + L' D^{-1} L)
logdet_lemma = np.sum(np.log(d)) + np.linalg.slogdet(core)[1]
logdet_direct = np.linalg.slogdet(Sigma)[1]
print("solve discrepancy:", np.abs(w_woodbury - w_direct).max() / np.abs(w_direct).max())
print("log-det discrepancy:", abs(logdet_lemma - logdet_direct))
# <</factor>>
rel_solve = np.abs(w_woodbury - w_direct).max() / np.abs(w_direct).max()
assert rel_solve < 1e-8
assert abs(logdet_lemma - logdet_direct) < 1e-8 * abs(logdet_direct)

gen = Generated("ch01", "woodbury", prefix="smw")
gen.int("n", n)
gen.int("p", p)
gen.num("worst", worst, 0, sci=True)
gen.int("imax", i_max + 1)
gen.num("hmax", h[i_max], 3)
gen.num("detratio", ratio, 3)
gen.int("N", N)
gen.int("k", k)
gen.num("relsolve", rel_solve, 0, sci=True)
gen.num("logdet", logdet_direct, 2)
gen.write()
