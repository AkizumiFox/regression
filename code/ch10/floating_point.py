"""Chapter 10, Section 1: floating-point arithmetic, rounding in inner products, and the
information that is lost when X^T X is formed.

Everything here is exact or deterministic: the rounding errors are compared with bounds
computed in exact rational arithmetic (fractions.Fraction).
"""
from fractions import Fraction

import numpy as np

from regbook import Generated

gen = Generated("ch10", "floating_point", prefix="fp")

# <<unit>>
eps = np.finfo(float).eps          # spacing of the doubles just above 1: 2^-52
u = eps / 2                        # unit roundoff: 2^-53
print("eps =", eps, "  u =", u)
print("1 + u == 1:", 1.0 + u == 1.0, "   1 + eps == 1:", 1.0 + eps == 1.0)
print("0.1 + 0.2 == 0.3:", 0.1 + 0.2 == 0.3)
# <</unit>>

assert eps == 2.0 ** -52 and u == 2.0 ** -53
assert 1.0 + u == 1.0 and 1.0 + eps > 1.0
assert Fraction(0.1) != Fraction(1, 10)          # 0.1 is not a double
gen.num("eps", eps, 2, sci=True)
gen.num("u", u, 2, sci=True)


# ---- rounding error of inner products against the bound gamma_n |x|'|y| -------------
# <<inner>>
def gamma(k, u=2.0 ** -53):
    return k * u / (1 - k * u)


def dot_recursive(x, y):
    """Inner product by recursive summation, one rounding per operation."""
    s = 0.0
    for a, b in zip(x, y):
        s = s + a * b
    return s


rng = np.random.default_rng(101)
worst = 0.0
for trial in range(200):
    n = int(rng.integers(2, 400))
    x = rng.normal(size=n) * 10.0 ** rng.integers(-3, 4, size=n)
    y = rng.normal(size=n)
    exact = sum(Fraction(a) * Fraction(b) for a, b in zip(x, y))       # exact rational value
    error = abs(Fraction(dot_recursive(x, y)) - exact)
    bound = Fraction(gamma(n)) * sum(abs(Fraction(a) * Fraction(b)) for a, b in zip(x, y))
    assert error <= bound
    worst = max(worst, float(error / bound))
print(f"largest error / bound over 200 inner products: {worst:.3f}")
# <</inner>>

assert worst < 0.5               # the bound is a worst case; typical errors are smaller
gen.num("dot_worst_ratio", worst, 3)
gen.num("gamma_1e6", gamma(10 ** 6), 1, sci=True)

# ---- the Lauchli matrix: X^T X rounds to a singular matrix ------------------------------
# <<lauchli>>
delta = 1e-9
X = np.array([[1.0, 1.0],
              [delta, 0.0],
              [0.0, delta]])
y = X @ np.array([1.0, 1.0])            # exact data: the solution is b = (1, 1), residual 0

A = X.T @ X                             # 1 + delta^2 rounds to 1
print("computed X'X =\n", A)
print("cond(X) =", np.linalg.cond(X))
try:
    np.linalg.cholesky(A)
except np.linalg.LinAlgError as err:
    print("Cholesky of the computed X'X fails:", err)

Q, R = np.linalg.qr(X)
b_qr = np.linalg.solve(R, Q.T @ y)
print("QR solution:", b_qr)
# <</lauchli>>

assert np.all(A == 1.0)                                  # every entry rounded to exactly 1
assert np.linalg.matrix_rank(A) == 1
exact_gram = [[Fraction(1) + Fraction(delta) ** 2, Fraction(1)], [Fraction(1), Fraction(1) + Fraction(delta) ** 2]]
assert exact_gram[0][0] != 1                             # the exact Gram matrix is nonsingular
assert np.allclose(b_qr, [1.0, 1.0], rtol=1e-6)
kappa_l = np.linalg.cond(X)
assert abs(kappa_l - np.sqrt(2 + delta ** 2) / delta) / kappa_l < 1e-8
gen.num("lauchli_kappa", kappa_l, 2, sci=True)
gen.num("lauchli_kappa_sq", kappa_l ** 2, 1, sci=True)
gen.num("lauchli_err_qr", np.linalg.norm(b_qr - 1) / np.sqrt(2), 1, sci=True)

# ---- the Gram-rounding bound: when is fl(X^T X) guaranteed positive definite? --------
# prp-cmp-gram-rounding: if gamma_n * p * kappa(X)^2 < 1 the computed Gram matrix is positive definite.
kappa_threshold = 1 / np.sqrt(u)
gen.num("kappa_threshold", kappa_threshold, 1, sci=True)
assert 6e7 < kappa_threshold < 1e8

gen.write()
