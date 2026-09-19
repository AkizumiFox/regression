"""Chapter 4, Section 2: mean and variance of two quadratic estimators of sigma^2.

S^2 = Y^T C Y/(n-1) with C = I - J/n, and the mean-square successive difference
Q = sum (Y_{i+1} - Y_i)^2 / (2(n-1)) = Y^T D^T D Y/(2(n-1)). With a linear trend in the
mean, both are biased; the formulas E = tr(A Sigma) + mu^T A mu and
Var = 2 tr((A Sigma)^2) + 4 mu^T A Sigma A mu are checked by simulation.
"""
import numpy as np

from regbook import Generated

rng = np.random.default_rng(4402)

n, sigma, theta = 20, 1.0, 0.1
t = np.arange(1, n + 1)
mu = 5.0 + theta * t                                   # mean with a linear trend
Sigma = sigma**2 * np.eye(n)

# <<moments>>
def qf_moments(A, mu, Sigma):
    """Mean and variance of Y^T A Y for Y ~ N(mu, Sigma), A symmetric."""
    AS = A @ Sigma
    mean = np.trace(AS) + mu @ A @ mu
    var = 2 * np.trace(AS @ AS) + 4 * mu @ AS @ A @ mu
    return mean, var


C = np.eye(n) - np.ones((n, n)) / n                    # centring matrix
D = np.diff(np.eye(n), axis=0)                         # (n-1) x n first differences
A_s2 = C / (n - 1)
A_q = D.T @ D / (2 * (n - 1))
for name, A in [("S2", A_s2), ("Q", A_q)]:
    print(name, qf_moments(A, mu, Sigma), qf_moments(A, 0 * mu, Sigma))
# <</moments>>

# closed forms quoted in the text
m_s2, v_s2 = qf_moments(A_s2, mu, Sigma)
m_q, v_q = qf_moments(A_q, mu, Sigma)
m_s20, v_s20 = qf_moments(A_s2, 0 * mu, Sigma)
m_q0, v_q0 = qf_moments(A_q, 0 * mu, Sigma)
assert np.isclose(m_s20, sigma**2) and np.isclose(m_q0, sigma**2)
assert np.isclose(v_s20, 2 * sigma**4 / (n - 1))
assert np.isclose(np.trace((D.T @ D) @ (D.T @ D)), 6 * n - 8)
assert np.isclose(v_q0, sigma**4 * (3 * n - 4) / (n - 1) ** 2)
assert np.isclose(m_s2 - sigma**2, theta**2 * n * (n + 1) / 12)
assert np.isclose(m_q - sigma**2, theta**2 / 2)

# simulation check of both moments, with and without the trend
reps = 400_000
E = sigma * rng.standard_normal((reps, n))
for mean_vec, targets in [(mu, [(A_s2, m_s2, v_s2), (A_q, m_q, v_q)]),
                          (0 * mu, [(A_s2, m_s20, v_s20), (A_q, m_q0, v_q0)])]:
    Y = mean_vec + E
    for A, m, v in targets:
        vals = np.einsum("ij,jk,ik->i", Y, A, Y)
        assert abs(vals.mean() - m) < 5 * np.sqrt(v / reps)
        assert abs(vals.var() / v - 1) < 0.03

# covariance of the linear and quadratic forms: Cov(Y, Y^T A Y) = 2 Sigma A mu
Y = mu + E
vals = np.einsum("ij,jk,ik->i", Y, A_s2, Y)
cov_emp = ((Y - Y.mean(0)) * (vals - vals.mean())[:, None]).mean(0)
assert np.allclose(cov_emp, 2 * Sigma @ A_s2 @ mu, atol=0.01)

gen = Generated("ch04", "successive_differences", prefix="sd")
gen.int("n", n)
gen.num("theta", theta, 1)
gen.num("ms2", m_s2, 3)
gen.num("mq", m_q, 3)
gen.num("vs2zero", v_s20, 4)
gen.num("vqzero", v_q0, 4)
gen.num("vs2", v_s2, 4)
gen.num("vq", v_q, 4)
gen.num("effratio", v_q0 / v_s20, 3)
gen.num("rmses2", np.sqrt(v_s2 + (m_s2 - sigma**2) ** 2), 3)
gen.num("rmseq", np.sqrt(v_q + (m_q - sigma**2) ** 2), 3)
gen.text("reps", f"{reps:,}".replace(",", "{,}"))
gen.write()
