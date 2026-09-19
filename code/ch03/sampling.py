"""Chapter 3, Sections 1-2: generating multivariate normal vectors from A Z + mu,
evaluating the log density through a Cholesky factor, and the chi-squared law of
the Mahalanobis distance."""
import numpy as np
from scipy import stats

from regbook import Generated

rng = np.random.default_rng(301)

mu = np.array([1.0, -2.0, 0.5])
Sigma = np.array([[4.0, 1.2, -0.8],
                  [1.2, 2.0, 0.6],
                  [-0.8, 0.6, 1.5]])
N = 400_000

# <<draw>>
def draw_normal(mu, A, size, rng):
    """Rows are independent copies of mu + A Z with Z ~ N(0, I)."""
    Z = rng.standard_normal((size, A.shape[1]))
    return mu + Z @ A.T

L = np.linalg.cholesky(Sigma)                  # lower triangular, L L^T = Sigma
lam, U = np.linalg.eigh(Sigma)
root = U @ np.diag(np.sqrt(lam)) @ U.T         # symmetric square root

Y_chol = draw_normal(mu, L, N, rng)
Y_root = draw_normal(mu, root, N, rng)
for Y in (Y_chol, Y_root):
    print(np.round(np.cov(Y, rowvar=False) - Sigma, 3))
# <</draw>>

assert np.allclose(L @ L.T, Sigma) and np.allclose(root @ root, Sigma)
for Y in (Y_chol, Y_root):
    assert np.allclose(Y.mean(axis=0), mu, atol=0.02)
    assert np.allclose(np.cov(Y, rowvar=False), Sigma, atol=0.05)
# same law: a linear combination has matching quantiles under both factors
a = np.array([0.3, -1.0, 2.0])
q = [0.05, 0.25, 0.5, 0.75, 0.95]
exact = stats.norm.ppf(q, loc=a @ mu, scale=np.sqrt(a @ Sigma @ a))
assert np.allclose(np.quantile(Y_chol @ a, q), exact, atol=0.05)
assert np.allclose(np.quantile(Y_root @ a, q), exact, atol=0.05)
cov_err = max(np.abs(np.cov(Y, rowvar=False) - Sigma).max() for Y in (Y_chol, Y_root))

# <<logpdf>>
def normal_logpdf(y, mu, Sigma):
    """log density of N(mu, Sigma) at the rows of y, Sigma positive definite."""
    L = np.linalg.cholesky(Sigma)
    w = np.linalg.solve(L, (y - mu).T)         # w = L^{-1}(y - mu), so w^T w = Delta^2
    D2 = np.sum(w**2, axis=0)
    logdet = 2.0 * np.sum(np.log(np.diag(L)))
    return -0.5 * (len(mu) * np.log(2 * np.pi) + logdet + D2), D2

logf, D2 = normal_logpdf(Y_chol, mu, Sigma)
inside = np.mean(D2 <= stats.chi2.ppf(0.95, df=3))
print(f"fraction inside the 95% ellipsoid: {inside:.4f}")
# <</logpdf>>

assert np.allclose(logf[:1000], stats.multivariate_normal(mu, Sigma).logpdf(Y_chol[:1000]))
assert abs(inside - 0.95) < 0.002
# Delta^2 has mean n and variance 2n
assert abs(D2.mean() - 3) < 0.02 and abs(D2.var() - 6) < 0.1

# singular normal: Sigma = A A^T with A of size 3 x 2 lives on a plane through mu
A = np.array([[1.0, 0.0], [1.0, 1.0], [0.0, 2.0]])
Ys = draw_normal(mu, A, 10_000, rng)
Q, _ = np.linalg.qr(A)
off_plane = np.abs((Ys - mu) - (Ys - mu) @ Q @ Q.T).max()
normal_dir = np.cross(A[:, 0], A[:, 1])
assert np.allclose(A @ A.T @ normal_dir, 0)
assert off_plane < 1e-12
assert np.linalg.matrix_rank(A @ A.T) == 2

# Example (three overlapping sums): U = (Z1+Z2, Z2-Z3, Z1+Z3) has rank 2 and U1 - U2 - U3 = 0
B = np.array([[1.0, 1.0, 0.0], [0.0, 1.0, -1.0], [1.0, 0.0, 1.0]])
assert np.allclose(B @ B.T, [[2, 1, 1], [1, 2, -1], [1, -1, 2]])
assert np.linalg.matrix_rank(B @ B.T) == 2 and np.allclose(np.array([1, -1, -1]) @ B, 0)
assert all(np.isclose(np.linalg.det((B @ B.T)[np.ix_(p, p)]), 3) for p in ([0, 1], [0, 2], [1, 2]))

gen = Generated("ch03", "sampling", prefix="samp")
gen.text("N", f"{N:,}".replace(",", "{,}"))
gen.num("coverr", cov_err, 3)
gen.num("inside", inside, 4)
gen.num("D2mean", D2.mean(), 3)
gen.num("D2var", D2.var(), 3)
gen.num("offplane", off_plane, 1, sci=True)
gen.write()
