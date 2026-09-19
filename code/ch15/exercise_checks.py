"""Chapter 15: checks of the numbers used in exercise solutions."""
import itertools

import numpy as np
from scipy import stats

from regbook import Generated

gen = Generated("ch15", "exercise_checks")

# Section 15.2, B1: plots per rate for power 0.8 for mu_5 - mu_3 = sigma (delta^2 = m / 2)
def power_contrast(m):
    nu = 5 * (m - 1)
    return stats.ncf.sf(stats.f.ppf(0.95, 1, nu), 1, nu, m / 2)

m_needed = next(m for m in range(2, 100) if power_contrast(m) >= 0.8)
assert m_needed == 17 and power_contrast(16) < 0.8
z = stats.norm.ppf(0.975) + stats.norm.ppf(0.8)
assert abs(2 * z ** 2 - 15.7) < 0.05
gen.int("m_needed", m_needed)

# Section 15.2, A1: estimate and standard error from the rounded means
means = np.array([4.005, 5.237, 6.057, 6.515, 6.652])
c = np.array([-0.5, -0.5, 0, 0.5, 0.5])
assert np.isclose(c @ means, 1.9625)
assert np.isclose(0.5211 * np.sqrt(np.sum(c ** 2) / 6), 0.2127, atol=5e-4)

# Section 15.4, B1: completing the two planned contrasts
C = np.array([[-1, 0.25, 0.25, 0.25, 0.25], [0, 0, -1, 0, 1], [0, 1, 0, -1, 0], [0, 1, -1, 1, -1]])
G = C @ C.T
assert np.allclose(G, np.diag(np.diag(G))) and np.allclose(C.sum(axis=1), 0)

# Section 15.4, B2: orthogonal polynomials for rates 0, 1, 2, 4
s = np.array([0.0, 1, 2, 4])
Q, _ = np.linalg.qr(np.vander(s, 3, increasing=True))
lin, quad = np.array([-7, -3, 1, 9.0]), np.array([7, -4, -8, 5.0])
assert np.isclose(abs(Q[:, 1] @ lin) / np.linalg.norm(lin), 1)
assert np.isclose(abs(Q[:, 2] @ quad) / np.linalg.norm(quad), 1)

# Section 15.4, B3: sizes 2, 4, 4, 2
nk = np.array([2, 4, 4, 2.0])
L, Qd, Cu = np.array([-3, -1, 1, 3.0]), np.array([1, -1, -1, 1.0]), np.array([-1, 3, -3, 1.0])
assert np.isclose(np.sum(L * Qd / nk), 0) and np.isclose(np.sum(L * Cu / nk), 1.5)

# Section 15.4, B4: tree contrast from rounded means
assert np.isclose(24 * 6 / 30 * (means[1:].mean() - means[0]) ** 2, 21.375, atol=1e-3)

# Section 15.5, A1: expected mean squares with variances 1, 2, 4, 8
var = np.array([1, 2, 4, 8.0])
def ems(sizes):
    nk = np.array(sizes, float)
    n, g = nk.sum(), nk.size
    return np.sum((1 - nk / n) * var) / (g - 1), np.sum((nk - 1) * var) / (n - g)
assert np.allclose(ems([10] * 4), [3.75, 3.75])
assert np.allclose(ems([4, 6, 10, 20]), [3.2, 201 / 36])
assert np.allclose(ems([20, 10, 6, 4]), [4.2, 2.25])

# Section 15.5, B4 and the allocation example: guaranteed h
h = lambda sizes: min(a * b / (a + b) for a, b in itertools.combinations(sizes, 2))
assert np.isclose(h([6] * 5), 3) and np.isclose(h([10, 5, 5, 5, 5]), 2.5)
assert np.isclose(8 * h([7, 9, 11, 13]), 31.5)

# Section 15.5, C1: level of the F test with intraclass correlation 0.1, m = 10, g = 4
factor = 1 + 10 * 0.1 / 0.9
level = stats.f.sf(stats.f.ppf(0.95, 3, 36) / factor, 3, 36)
rng = np.random.default_rng(3)
sims = 20000
a = rng.normal(0, np.sqrt(0.1), (sims, 4, 1))
e = rng.normal(0, np.sqrt(0.9), (sims, 4, 10))
yk = a + e
msb = 10 * np.sum((yk.mean(axis=2) - yk.mean(axis=(1, 2))[:, None]) ** 2, axis=1) / 3
mse = np.sum((yk - yk.mean(axis=2, keepdims=True)) ** 2, axis=(1, 2)) / 36
assert abs(np.mean(msb / mse > stats.f.ppf(0.95, 3, 36)) - level) < 0.01
gen.num("intraclass_level", level, 3)

# Section 15.3, B2: null probability that F < 1 (omega-squared numerator negative)
p_below_one = stats.f.cdf(1, 4, 25)
assert 0.55 < p_below_one < 0.6
assert 0.65 < stats.f.cdf(1, 1, 25) < 0.7 and stats.f.cdf(1, 20, 5) < 0.5
gen.num("p_f_below_one", p_below_one, 2)

# Section 15.5, C2
assert np.isclose(12 / 7.5, 1.6)
gen.write()
