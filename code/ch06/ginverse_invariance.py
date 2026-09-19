"""Chapter 6, Section 4: one projection, many coefficient vectors.

A one-way layout with an intercept is rank deficient. Different generalized
inverses of X^T X give different least squares coefficient vectors, but all of
them give the same projection M and therefore the same fitted values.
"""
import numpy as np

from regbook import Generated

# <<layout>>
groups = np.array([0, 0, 0, 1, 1, 2, 2, 2, 2])     # three groups, sizes 3, 2, 4
y = np.array([5.1, 6.3, 5.7, 8.2, 7.4, 4.0, 4.9, 3.8, 4.5])
n, g = len(y), 3
X = np.column_stack([np.ones(n), np.eye(g)[groups]])  # [1, indicator columns]
A = X.T @ X
print("n =", n, " p =", X.shape[1], " rank(X) =", np.linalg.matrix_rank(X))
# <</layout>>


# <<ginverses>>
def g_drop(j):
    """Generalized inverse obtained by deleting row and column j of A."""
    keep = [k for k in range(A.shape[0]) if k != j]
    G = np.zeros_like(A)
    G[np.ix_(keep, keep)] = np.linalg.inv(A[np.ix_(keep, keep)])
    return G


rng = np.random.default_rng(6)
G0 = np.linalg.pinv(A)                             # Moore-Penrose inverse
U, V = rng.normal(size=A.shape), rng.normal(size=A.shape)
G_random = G0 + (np.eye(4) - G0 @ A) @ U + V @ (np.eye(4) - A @ G0)

ginverses = {
    "Moore-Penrose": G0,
    "drop intercept": g_drop(0),
    "drop group 1": g_drop(1),
    "random": G_random,
}
# <</ginverses>>

# <<compare>>
M_ref = X @ G0 @ X.T
for name, G in ginverses.items():
    assert np.allclose(A @ G @ A, A)               # G really is a g-inverse
    beta = G @ X.T @ y
    M = X @ G @ X.T
    print(f"{name:15s} beta = {np.round(beta, 3)}",
          f" max|M - M_ref| = {np.abs(M - M_ref).max():.1e}",
          f" group1-group2 = {beta[1] - beta[2]:.3f}")
# <</compare>>

gen = Generated("ch06", "ginverse_invariance", prefix="gi")
betas = {}
for key, (name, G) in zip(["mp", "dropint", "dropgone", "random"], ginverses.items()):
    beta = G @ X.T @ y
    M = X @ G @ X.T
    assert np.allclose(M, M_ref, atol=1e-10)
    assert np.allclose(M, M.T, atol=1e-10) and np.allclose(M @ M, M, atol=1e-10)
    betas[key] = beta
    for i, b in enumerate(beta):
        gen.num(f"{key}:b{i}", b, 3)
    gen.num(f"{key}:contrast", beta[1] - beta[2], 3)
    d = np.abs(M - M_ref).max()
    gen.num(f"{key}:mdiff", d, 1, sci=True) if d > 0 else gen.text(f"{key}:mdiff", "0")
means = [y[groups == k].mean() for k in range(g)]
for k, m in enumerate(means):
    gen.num(f"mean{k+1}", m, 3)
assert np.allclose(M_ref @ y, np.array(means)[groups])
assert not np.allclose(betas["mp"], betas["dropint"])
gen.int("rank", np.linalg.matrix_rank(X))
gen.num("trM", np.trace(M_ref), 3)
gen.write()
