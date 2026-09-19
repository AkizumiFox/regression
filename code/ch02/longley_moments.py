"""Chapter 2, Section 6: sample covariance and correlation matrices of the Longley data.

Sixteen annual observations (1947-1962) on seven US macroeconomic series, public domain,
shipped with statsmodels (statsmodels.datasets.longley). The eigenvalues show how scale
dominates a covariance matrix and how a nearly degenerate linear combination shows up
as a tiny eigenvalue of the correlation matrix.
"""
import numpy as np
import statsmodels.api as sm

from regbook import Generated

data = sm.datasets.longley.load_pandas().data
names = list(data.columns)

# <<sample>>
Y = data.to_numpy(dtype=float)                    # n x p data matrix, one row per year
n, p = Y.shape
C = np.eye(n) - np.ones((n, n)) / n               # centring matrix I - J/n
S = Y.T @ C @ Y / (n - 1)                         # sample covariance matrix
d = np.sqrt(np.diag(S))
R = S / np.outer(d, d)                            # D^{-1} S D^{-1}, D = diag(d)

lam_S = np.linalg.eigvalsh(S)
lam_R, V = np.linalg.eigh(R)                      # ascending eigenvalues
print("eigenvalues of S:", lam_S)
print("eigenvalues of R:", lam_R)
v = V[:, 0]                                       # direction of least variance
Z = (Y - Y.mean(0)) / d                           # standardized variables
print("sample variance of Z v:", np.var(Z @ v, ddof=1))
# <</sample>>

assert np.allclose(S, np.cov(Y, rowvar=False))
assert np.allclose(R, np.corrcoef(Y, rowvar=False))
assert np.all(lam_S > 0) and np.all(lam_R > 0)                 # nonnegative definite (here pd)
assert np.isclose(lam_R.sum(), p)                              # tr R = p
assert np.allclose(np.diag(d) @ R @ np.diag(d), S)
assert np.isclose(np.var(Z @ v, ddof=1), lam_R[0])             # variance of v'Z is v'Rv
assert np.all(np.abs(R) <= 1 + 1e-12)
top = np.argsort(-np.abs(v))[:2]
assert {names[i] for i in top} == {"GNP", "YEAR"}

# with n <= p the sample covariance must be singular: rank at most n - 1
k = 5
Sk = np.cov(Y[:k], rowvar=False)
rank_k = np.linalg.matrix_rank(Sk, tol=1e-8 * np.linalg.norm(Sk))
assert rank_k == k - 1

gen = Generated("ch02", "longley_moments", prefix="lon")
gen.int("n", n)
gen.int("p", p)
gen.num("S_max", lam_S[-1], 2, sci=True)
gen.num("S_min", lam_S[0], 2, sci=True)
gen.num("S_ratio", lam_S[-1] / lam_S[0], 1, sci=True)
gen.num("R_max", lam_R[-1], 3)
gen.num("R_second", lam_R[-2], 3)
gen.num("R_min", lam_R[0], 5)
gen.num("R_min_sci", lam_R[0], 2, sci=True)
gen.num("R_ratio", lam_R[-1] / lam_R[0], 0)
gen.num("R_share", lam_R[-1] / p * 100, 0)
gen.num("sd_min", np.sqrt(lam_R[0]), 3)
sgn = np.sign(v[names.index("YEAR")])
gen.num("v_year", sgn * v[names.index("YEAR")], 3)
gen.num("v_gnp", sgn * v[names.index("GNP")], 3)
gen.num("r_gnp_year", R[names.index("GNP"), names.index("YEAR")], 3)
gen.num("sd_gnp", d[names.index("GNP")], 0)
gen.num("sd_defl", d[names.index("GNPDEFL")], 2)
gen.int("k", k)
gen.int("rank_k", rank_k)
gen.write()
