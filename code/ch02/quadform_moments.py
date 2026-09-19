"""Chapter 2, Section 3: mean and variance of a quadratic form without normality.

Part 1 checks E(Y'AY) = tr(A Sigma) + mu'A mu for a correlated, skewed random vector,
exactly (by enumerating a discrete distribution) and by simulation.
Part 2 checks the variance formula for independent components with common central
moments (Theorem thm:rv-quadform-variance) on the residual sum of squares of a
straight-line fit when the true mean is curved and the errors are skewed.
"""
import itertools

import numpy as np

from regbook import Generated

rng = np.random.default_rng(20260917)

# ---- Part 1: the mean of a quadratic form -----------------------------------
# <<mean>>
mu = np.array([1.0, -1.0, 2.0])
B = np.array([[1.0, 0.0, 0.0],
              [0.5, 1.0, 0.0],
              [-0.5, 0.8, 0.6]])
Sigma = B @ B.T                                   # Cov(mu + B U) when Cov(U) = I
A = np.array([[2.0, 1.0, 0.0],
              [1.0, 1.0, -1.0],
              [0.0, -1.0, 3.0]])

theory = np.trace(A @ Sigma) + mu @ A @ mu        # Theorem: no normality needed

N = 1_000_000
U = rng.exponential(size=(N, 3)) - 1.0            # skewed, mean 0, variance 1
Y = mu + U @ B.T
Q = np.einsum("ij,jk,ik->i", Y, A, Y)             # Y_i' A Y_i for every draw
print(f"theory {theory:.4f}   simulation {Q.mean():.4f} ± {Q.std() / np.sqrt(N):.4f}")
# <</mean>>

se = Q.std() / np.sqrt(N)
assert abs(Q.mean() - theory) < 5 * se
# the sample covariance of the draws is close to Sigma
assert np.allclose(np.cov(Y, rowvar=False), Sigma, atol=0.02)

# exact check: U_i iid two-point, values -1 (prob 2/3) and 2 (prob 1/3): mean 0, variance 2
vals, probs = np.array([-1.0, 2.0]), np.array([2 / 3, 1 / 3])
assert np.isclose(probs @ vals, 0) and np.isclose(probs @ vals**2, 2)
exact = 0.0
for idx in itertools.product(range(2), repeat=3):
    u = vals[list(idx)]
    y = mu + B @ u
    exact += np.prod(probs[list(idx)]) * (y @ A @ y)
assert np.isclose(exact, np.trace(A @ (2 * Sigma)) + mu @ A @ mu, rtol=1e-13)

# ---- Part 2: the variance needs third and fourth moments --------------------
x = np.arange(1.0, 7.0)                           # n = 6 design points
X = np.column_stack([np.ones_like(x), x])
M = X @ np.linalg.solve(X.T @ X, X.T)
A2 = np.eye(6) - M                                # RSS = Y'(I - M)Y for a straight line
theta = 0.25 * (x - 3.5) ** 2                     # true mean is curved: the line is wrong
a = np.diag(A2)


# <<variance>>
def var_quadform(A, theta, m2, m3, m4):
    """Var(Y'AY) for independent Y_i with means theta_i and common central moments."""
    a = np.diag(A)
    return ((m4 - 3 * m2**2) * a @ a + 2 * m2**2 * np.trace(A @ A)
            + 4 * m2 * theta @ A @ A @ theta + 4 * m3 * theta @ A @ a)


# A2 = I - M for a straight line at x = 1..6; theta = true (curved) mean
mean_rss = np.trace(A2) + theta @ A2 @ theta      # sigma^2 = 1
var_expo = var_quadform(A2, theta, 1.0, 2.0, 9.0) # centred exponential errors
var_norm = var_quadform(A2, theta, 1.0, 0.0, 3.0) # what normal errors would give

E = rng.exponential(size=(N, 6)) - 1.0
R = np.einsum("ij,jk,ik->i", theta + E, A2, theta + E)
print(f"mean  {mean_rss:.3f} vs {R.mean():.3f}")
print(f"var   {var_expo:.3f} vs {R.var():.3f}   (normal-theory value {var_norm:.3f})")
# <</variance>>

assert abs(R.mean() - mean_rss) < 5 * R.std() / np.sqrt(N)
assert abs(R.var() / var_expo - 1) < 0.04

# exact check of the variance formula by enumeration (two-point errors, 2^6 atoms)
m2, m3, m4 = (probs @ vals**k for k in (2, 3, 4))
e1 = e2 = 0.0
for idx in itertools.product(range(2), repeat=6):
    y = theta + vals[list(idx)]
    q = y @ A2 @ y
    w = np.prod(probs[list(idx)])
    e1 += w * q
    e2 += w * q * q
assert np.isclose(e1, m2 * np.trace(A2) + theta @ A2 @ theta, rtol=1e-12)
assert np.isclose(e2 - e1**2, var_quadform(A2, theta, m2, m3, m4), rtol=1e-10)
# the normal case reduces to 2 sigma^4 tr(A^2) + 4 sigma^2 theta'A^2 theta
assert np.isclose(var_norm, 2 * np.trace(A2 @ A2) + 4 * theta @ A2 @ A2 @ theta)

gen = Generated("ch02", "quadform_moments", prefix="qfm")
gen.num("trASigma", np.trace(A @ Sigma), 1)
gen.num("muAmu", mu @ A @ mu, 1)
gen.num("theory", theory, 1)
gen.num("sim", Q.mean(), 3)
gen.num("simse", se, 3)
gen.text("N", "10^6")
gen.num("trA", np.trace(A2), 0)
gen.num("bias", theta @ A2 @ theta, 3)
gen.num("mean_rss", mean_rss, 3)
gen.num("sim_rss", R.mean(), 3)
gen.num("var_expo", var_expo, 2)
gen.num("var_norm", var_norm, 2)
gen.num("sim_var", R.var(), 2)
gen.num("kurt_term", (9.0 - 3.0) * a @ a, 2)
gen.num("skew_term", 4 * 2.0 * theta @ A2 @ a, 2)
gen.num("ratio", var_expo / var_norm, 2)
gen.write()
