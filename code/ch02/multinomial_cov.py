"""Chapter 2, Section 2: a singular covariance matrix (multinomial counts).

Counts in k = 3 categories from m independent trials always sum to m, so the covariance
matrix m(diag(pi) - pi pi') is singular with null vector 1, and every draw lies in the
plane {y : 1'y = m} = E(Y) + C(Sigma).
"""
import numpy as np

from regbook import Generated

rng = np.random.default_rng(11)

# <<multinomial>>
m, pi = 20, np.array([0.5, 0.3, 0.2])
mean = m * pi
Sigma = m * (np.diag(pi) - np.outer(pi, pi))

lam, U = np.linalg.eigh(Sigma)
print("eigenvalues:", np.round(lam, 6) + 0.0)     # one of them is zero
print("Sigma @ 1 =", np.round(Sigma @ np.ones(3), 12) + 0.0)

Y = rng.multinomial(m, pi, size=5000).astype(float)
print("1'Y for the first draws:", Y[:5].sum(axis=1))
# <</multinomial>>

one = np.ones(3)
assert np.allclose(Sigma @ one, 0) and np.linalg.matrix_rank(Sigma) == 2
assert np.isclose(lam[0], 0, atol=1e-12) and np.all(lam[1:] > 0)
assert np.allclose(Y.sum(axis=1), m)
# Y - E(Y) lies in C(Sigma): projecting onto the null space gives zero
P0 = np.outer(one, one) / 3
assert np.allclose((Y - mean) @ P0, 0)
S = np.cov(Y, rowvar=False)
assert np.allclose(S @ one, 0, atol=1e-10) and np.allclose(S, Sigma, atol=0.2)
# whitening on the support: Sigma^{+1/2} (Y - mu) has covariance P = I - J/3
root_pinv = U[:, 1:] @ np.diag(lam[1:] ** -0.5) @ U[:, 1:].T
assert np.allclose(root_pinv @ Sigma @ root_pinv, np.eye(3) - P0)

gen = Generated("ch02", "multinomial_cov", prefix="mn")
gen.int("m", m)
gen.num("lam1", lam[2], 3)
gen.num("lam2", lam[1], 3)
gen.int("draws", len(Y))
gen.write()
