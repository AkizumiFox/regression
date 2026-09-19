"""Chapter 8, Section 4: side conditions pick one least squares solution; estimable functions do not care.

A one-way layout with four groups and an intercept. Five side conditions, each a single
nonestimable linear restriction, give five different coefficient vectors, all solving the
normal equations. A restriction on an estimable function is not a side condition: it changes
the fit. Simulated data with a fixed seed.
"""
import numpy as np
import scipy.linalg
import statsmodels.api as sm

from regbook import Generated

# <<layout>>
import numpy as np

sizes = [5, 3, 6, 4]
groups = np.repeat(np.arange(4), sizes)
rng = np.random.default_rng(84)
y = np.array([12.0, 15.0, 11.0, 14.0])[groups] + rng.normal(scale=1.0, size=groups.size)
X = np.column_stack([np.ones(groups.size), np.eye(4)[groups]])   # (mu, alpha_1..alpha_4)
A, c = X.T @ X, X.T @ y
# <</layout>>

# <<sideconditions>>
n_k = np.array(sizes, dtype=float)
side = {"alpha_1 = 0":         [0, 1, 0, 0, 0],
        "alpha_4 = 0":         [0, 0, 0, 0, 1],
        "sum alpha_k = 0":     [0, 1, 1, 1, 1],
        "sum n_k alpha_k = 0": [0, *n_k],
        "mu = 0":              [1, 0, 0, 0, 0]}
solutions = {}
for name, t in side.items():
    T = np.array([t], dtype=float)
    G = np.linalg.inv(A + T.T @ T)          # (X'X + T'T)^{-1}: a g-inverse of X'X
    b = G @ c
    solutions[name] = b
    print(f"{name:20s} b = {np.round(b, 3)}   T b = {(T @ b)[0]:.1e}")
b_mn = np.linalg.pinv(X) @ y
print(f"{'minimum norm':20s} b = {np.round(b_mn, 3)}",
      f"  mu - sum(alpha) = {b_mn[0] - b_mn[1:].sum():.1e}")
# <</sideconditions>>

means = np.array([y[groups == k].mean() for k in range(4)])
null = np.array([1.0, -1, -1, -1, -1])
assert np.linalg.matrix_rank(X) == 4 and np.allclose(X @ null, 0)
for name, t in side.items():
    T = np.array([t], dtype=float)
    assert abs(T @ null) > 1e-8                                # nonestimable
    assert np.linalg.matrix_rank(np.vstack([X, T])) == 5       # right rank
    G = np.linalg.inv(A + T.T @ T)
    assert np.allclose(A @ G @ A, A)                           # a generalized inverse
    b = solutions[name]
    assert np.allclose(A @ b, c) and np.allclose(T @ b, 0)
    assert np.allclose(b[0] + b[1:], means)                    # mu + alpha_k = group mean
    assert np.isclose(b[1] - b[2], means[0] - means[1])
assert np.allclose(A @ b_mn, c) and np.isclose(b_mn[0], b_mn[1:].sum())
# the minimum-norm solution is the side-condition solution with T = null^T
assert np.allclose(np.linalg.inv(A + np.outer(null, null)) @ c, b_mn)
# interpretations
assert np.isclose(solutions["sum alpha_k = 0"][0], means.mean())
assert np.isclose(solutions["sum n_k alpha_k = 0"][0], y.mean())
assert np.isclose(solutions["alpha_1 = 0"][0], means[0])

# <<estimablerestriction>>
# "alpha_1 = alpha_2" is a restriction on an estimable function: it changes the model
T_bad = np.array([[0, 1, -1, 0, 0]], dtype=float)
K = np.block([[A, T_bad.T], [T_bad, np.zeros((1, 1))]])        # Lagrange (bordered) system
b_bad = np.linalg.lstsq(K, np.append(c, 0.0), rcond=None)[0][:5]
sse_free = np.sum((y - X @ solutions["alpha_1 = 0"]) ** 2)
sse_bad = np.sum((y - X @ b_bad) ** 2)
print(f"SSE without restriction {sse_free:.3f};  with alpha_1 = alpha_2 {sse_bad:.3f}")
# <</estimablerestriction>>
assert abs(T_bad @ null) < 1e-12
assert sse_bad > sse_free + 1

# <<software>>
b_sm = sm.OLS(y, X).fit().params                 # statsmodels default: method="pinv"
b_np = np.linalg.lstsq(X, y, rcond=None)[0]      # LAPACK gelsd: minimum-norm solution
_, R, piv = scipy.linalg.qr(X, mode="economic", pivoting=True)
print("statsmodels:", np.round(b_sm, 3))
print("numpy lstsq:", np.round(b_np, 3))
print("pivoted QR column order:", piv, " |R_kk| =", np.round(np.abs(np.diag(R)), 3))
# <</software>>
assert np.allclose(b_sm, b_mn) and np.allclose(b_np, b_mn)
assert np.abs(np.diag(R))[-1] < 1e-10
dropped = int(piv[-1])
assert list(piv) == [0, 3, 1, 4, 2] and piv[-1] == 1 + int(np.argmin(sizes))

gen = Generated("ch08", "side_conditions", prefix="sc")
keys = {"alpha_1 = 0": "ref1", "alpha_4 = 0": "ref4", "sum alpha_k = 0": "sum",
        "sum n_k alpha_k = 0": "wsum", "mu = 0": "cell"}
for name, key in keys.items():
    for i, v in enumerate(solutions[name]):
        gen.num(f"{key}:b{i}", 0.0 if abs(v) < 1e-9 else v, 3)
for i, v in enumerate(b_mn):
    gen.num(f"mn:b{i}", v, 3)
for k in range(4):
    gen.num(f"mean{k + 1}", means[k], 3)
gen.num("ybar", y.mean(), 3)
gen.num("meanofmeans", means.mean(), 3)
gen.num("sse", sse_free, 3)
gen.num("ssebad", sse_bad, 3)
gen.int("dropped", dropped)
gen.write()
