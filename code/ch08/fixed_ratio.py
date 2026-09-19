"""Chapter 8, Section 2: a fertilizer trial in which nitrogen and phosphorus always came in the ratio 2:1.

The model E(y) = b0 + b1 N + b2 P has a model matrix of rank 2. The mean response is
estimable at every dose on the line P = N/2 and at no dose off it. Simulated yields with a
fixed seed.
"""
import numpy as np

from regbook import Generated

# <<design>>
import numpy as np

N = np.array([0, 0, 20, 20, 40, 40, 60, 60, 80, 80], dtype=float)   # kg/ha nitrogen
P = N / 2                                                           # phosphorus, always half
X = np.column_stack([np.ones(N.size), N, P])
rng = np.random.default_rng(80)
y = 3.0 + 0.020 * N + 0.030 * P + rng.normal(scale=0.25, size=N.size)   # yield, t/ha
print("rank(X) =", np.linalg.matrix_rank(X), "  null vector (0, 1, -2):", X @ [0, 1, -2])
# <</design>>

# <<solutions>>
b_min = np.linalg.pinv(X) @ y                                  # minimum-norm solution
b_noP = np.append(np.linalg.lstsq(X[:, :2], y, rcond=None)[0], 0.0)     # drop P
b_noN = np.insert(np.linalg.lstsq(X[:, [0, 2]], y, rcond=None)[0], 1, 0.0)  # drop N
new_points = {"N=40, P=20 (on the line)": [1, 40, 20],
              "N=40, P=0  (off the line)": [1, 40, 0]}
for name, x0 in new_points.items():
    print(name, [round(float(np.dot(x0, b)), 3) for b in (b_min, b_noP, b_noN)])
# <</solutions>>

null = np.array([0.0, 1.0, -2.0])
assert np.linalg.matrix_rank(X) == 2 and np.allclose(X @ null, 0)
for b in (b_min, b_noP, b_noN):
    assert np.allclose(X.T @ X @ b, X.T @ y)                    # all solve the normal equations
on, off = np.array([1, 40, 20.0]), np.array([1, 40, 0.0])
assert abs(on @ null) < 1e-12 and abs(off @ null) > 1        # on is in C(X^T), off is not
preds_on = [on @ b for b in (b_min, b_noP, b_noN)]
preds_off = [off @ b for b in (b_min, b_noP, b_noN)]
assert np.allclose(preds_on, preds_on[0])
assert np.ptp(preds_off) > 0.5

# two unbiased estimators of the same estimable function lambda = (1, 40, 20)
# (i) the average of the two plots that received exactly that dose
a1 = np.zeros(N.size)
a1[N == 40] = 0.5
# (ii) the least squares estimator rho^T M y, written as a^T y with a = M a1
M = X @ np.linalg.pinv(X)
a2 = M @ a1
assert np.allclose(a1 @ X, on) and np.allclose(a2 @ X, on)
v1, v2 = a1 @ a1, a2 @ a2                                     # variances in units of sigma^2
assert v2 < v1
assert np.isclose(a2 @ y, preds_on[0])

gen = Generated("ch08", "fixed_ratio", prefix="fr")
gen.int("n", N.size)
for key, b in zip(["min", "noP", "noN"], (b_min, b_noP, b_noN)):
    for i in range(3):
        gen.num(f"{key}:b{i}", b[i], 4)
    gen.num(f"{key}:on", on @ b, 3)
    gen.num(f"{key}:off", off @ b, 3)
gen.num("v1", v1, 3)
gen.num("v2", v2, 3)
gen.write()
