"""Chapter 27, Section 1: a biased estimator that beats least squares, in canonical coordinates.

Two regressors with sample correlation 0.95 and n = 50, scaled so that X^T X = n R. Errors have
sigma = 1 and the true coefficients are (1, 0.5). The script computes the exact mean squared
error of least squares, of dropping the minor principal direction, of halving it, and of the
oracle shrinkage factors, and checks the matrix comparison of Theorem (bias and variance).
"""
import numpy as np

from regbook import Generated

# <<example>>
n, rho, sigma = 50, 0.95, 1.0
XtX = n * np.array([[1.0, rho], [rho, 1.0]])
beta = np.array([1.0, 0.5])

evals, V = np.linalg.eigh(XtX)
order = np.argsort(evals)[::-1]
d2, V = evals[order], V[:, order]          # d_1^2 >= d_2^2, eigenvectors v_1, v_2
alpha = V.T @ beta                         # canonical coefficients
tau2 = d2 * alpha ** 2 / sigma ** 2        # squared signal-to-noise ratio of each direction


def mse(f):
    """Total MSE of the estimator that multiplies alpha_hat_j by f_j."""
    f = np.asarray(f, dtype=float)
    return np.sum(f ** 2 * sigma ** 2 / d2 + (1 - f) ** 2 * alpha ** 2)


f_oracle = tau2 / (1 + tau2)
for label, f in [("least squares", [1, 1]), ("drop v_2", [1, 0]),
                 ("halve v_2", [1, 0.5]), ("oracle", f_oracle)]:
    print(f"{label:14s} f = {np.round(f, 3)}   MSE = {mse(f):.4f}")
print("tau^2 =", np.round(tau2, 3))
# <</example>>

assert np.isclose(d2[0], n * (1 + rho)) and np.isclose(d2[1], n * (1 - rho))
assert mse([1, 0]) < mse([1, 1]) and mse(f_oracle) <= min(mse([1, 0]), mse([1, 0.5]))
# component 2: dropping beats keeping iff tau^2 < 1; halving iff tau^2 < 3
assert (tau2[1] < 1) and mse([1, 0.5]) < mse([1, 1])

# matrix comparison for "drop v_2": D = sigma^2 v2 v2^T / d2^2 (singular), b = -alpha_2 v_2
v2 = V[:, 1]
D = sigma ** 2 * np.outer(v2, v2) / d2[1]
b = -alpha[1] * v2
crit = b @ np.linalg.pinv(D) @ b
assert np.isclose(crit, tau2[1])
mse_ls = sigma ** 2 * np.linalg.inv(XtX)
A = np.outer(V[:, 0], V[:, 0])             # the estimator V diag(1, 0) V^T beta_hat
mse_drop = sigma ** 2 * A @ np.linalg.inv(XtX) @ A.T + np.outer(b, b)
assert np.linalg.eigvalsh(mse_ls - mse_drop).min() > -1e-12   # dominance in the matrix sense

gen = Generated("ch27", "mse_example")
gen.num("d2_1", d2[0], 1)
gen.num("d2_2", d2[1], 1)
gen.num("alpha_1", alpha[0], 4)
gen.num("alpha_2", abs(alpha[1]), 4)
gen.num("tau2_1", tau2[0], 1)
gen.num("tau2_2", tau2[1], 4)
gen.num("mse_ls", mse([1, 1]), 4)
gen.num("mse_drop", mse([1, 0]), 4)
gen.num("mse_half", mse([1, 0.5]), 4)
gen.num("mse_oracle", mse(f_oracle), 4)
gen.num("f_oracle_1", f_oracle[0], 3)
gen.num("f_oracle_2", f_oracle[1], 3)
gen.write()
