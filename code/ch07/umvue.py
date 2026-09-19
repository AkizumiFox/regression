"""Chapter 7, Section 4: sufficiency, Rao-Blackwell, UMVUE and the Cramer-Rao bound.

(i) Rao-Blackwellizing any linear unbiased estimator in the normal model gives least squares:
    the reflection y -> (2M - I) y keeps the sufficient statistic and the distribution fixed.
(ii) The UMVUE of sigma is a constant multiple of s; its constant and its variance.
(iii) Efficiency of s^2 relative to the Cramer-Rao bound.
"""
import numpy as np
from scipy import special

from regbook import Generated

# <<reflection>>
import numpy as np

x = np.array([1.0, 2.0, 2.5, 3.0, 4.0, 5.0, 5.5, 7.0, 8.0, 9.0, 10.5, 12.0])
n = len(x)
X = np.column_stack([np.ones(n), x])
M = X @ np.linalg.solve(X.T @ X, X.T)
R = 2 * M - np.eye(n)                            # reflection: keeps C(X), flips C(X)-perp

rng = np.random.default_rng(71)
y = X @ np.array([2.0, 0.5]) + rng.normal(size=n)
y_ref = R @ y
print("X^T y unchanged:", np.allclose(X.T @ y_ref, X.T @ y))
print("y^T y unchanged:", np.isclose(y_ref @ y_ref, y @ y))

a_end = np.zeros(n)
a_end[[0, -1]] = [-1.0, 1.0]
a_end /= x[-1] - x[0]                            # the end-point slope estimator
a_ls = X @ np.linalg.solve(X.T @ X, [0.0, 1.0])
print("average of a^T y over the two points:", (a_end @ y + a_end @ y_ref) / 2)
print("least squares slope:                 ", a_ls @ y)
# <</reflection>>

assert np.allclose(R @ R, np.eye(n)) and np.allclose(R, R.T)
assert np.allclose(R @ X, X)                     # the mean X beta is fixed, so the law N(X beta, s^2 I) is too
assert np.isclose((a_end @ y + a_end @ y_ref) / 2, a_ls @ y)
for _ in range(20):                              # the identity holds for every y and every LUE
    yy = rng.normal(size=n) * 3
    a = a_end + (np.eye(n) - M) @ rng.normal(size=n)
    a = a + a_ls - M @ a                         # force M a = a_ls, so a is unbiased for the slope
    assert np.allclose(X.T @ a, [0.0, 1.0])
    assert np.isclose((a @ yy + a @ (R @ yy)) / 2, a_ls @ yy)

# <<sigma>>
from scipy import special


def mean_s_factor(k):
    """b_k = E(s)/sigma when k s^2 / sigma^2 ~ chi^2(k)."""
    return np.sqrt(2 / k) * np.exp(special.gammaln((k + 1) / 2) - special.gammaln(k / 2))


for k in [2, 5, 10, 17, 50, 200]:
    b = mean_s_factor(k)
    print(f"k = {k:3d}   E(s)/sigma = {b:.4f}   UMVUE of sigma = s / {b:.4f}")
# <</sigma>>

k, nn = 17, 21                                   # the stack loss sizes of Section 7.3
b17 = mean_s_factor(k)
sims = np.sqrt(rng.chisquare(k, size=1_000_000) / k)   # s / sigma
assert abs(sims.mean() / b17 - 1) < 1e-3
assert abs((sims / b17).mean() - 1) < 1e-3
var_umvue_sigma = 1 / b17 ** 2 - 1               # Var(s / b_k) / sigma^2
assert abs((sims / b17).var() / var_umvue_sigma - 1) < 0.01
crb_sigma = 1 / (2 * nn)                         # Cramer-Rao bound for sigma, in units of sigma^2
eff_s2 = k / nn                                  # (2 sigma^4 / n) / (2 sigma^4 / (n - p))
assert var_umvue_sigma > crb_sigma

gen = Generated("ch07", "umvue", prefix="umv")
gen.num("b17", b17, 4)
gen.num("inv_b17", 1 / b17, 4)
gen.num("corr_pct", 100 * (1 / b17 - 1), 1)
gen.num("var_sigma", var_umvue_sigma, 4)
gen.num("crb_sigma", crb_sigma, 4)
gen.num("eff_sigma", crb_sigma / var_umvue_sigma, 3)
gen.num("eff_s2", eff_s2, 3)
gen.num("b2", mean_s_factor(2), 4)
gen.num("b200", mean_s_factor(200), 4)
gen.write()
