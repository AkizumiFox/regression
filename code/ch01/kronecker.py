"""Chapter 1, Section 10: Kronecker products and vec, checked numerically.

The balanced two-way layout (a rows, b columns, one observation per cell) has a model
matrix built from Kronecker products, and so do its projections.
"""
import numpy as np

from regbook import Generated

rng = np.random.default_rng(11)


def vec(M):
    return M.reshape(-1, order="F")                     # stack the columns


# <<kron>>
A, B, C = rng.normal(size=(3, 4)), rng.normal(size=(4, 2)), rng.normal(size=(2, 5))
assert np.allclose(vec(A @ B @ C), np.kron(C.T, A) @ vec(B))          # vec(ABC)

# balanced two-way layout: a = 3 rows, b = 4 columns, one observation per cell
a, b = 3, 4
one = lambda k: np.ones((k, 1))
X = np.hstack([np.kron(one(a), one(b)),        # intercept
               np.kron(np.eye(a), one(b)),     # row indicators
               np.kron(one(a), np.eye(b))])    # column indicators
print("rank(X) =", np.linalg.matrix_rank(X), "= a + b - 1 =", a + b - 1)
# <</kron>>
assert np.linalg.matrix_rank(X) == a + b - 1

# mixed-product rule, transpose, inverse, trace, determinant, eigenvalues
A2, P2 = rng.normal(size=(3, 3)), rng.normal(size=(2, 2))
B2, R2 = rng.normal(size=(3, 3)), rng.normal(size=(2, 2))
assert np.allclose(np.kron(A2, P2) @ np.kron(B2, R2), np.kron(A2 @ B2, P2 @ R2))
assert np.allclose(np.kron(A2, P2).T, np.kron(A2.T, P2.T))
assert np.allclose(np.linalg.inv(np.kron(A2, P2)), np.kron(np.linalg.inv(A2), np.linalg.inv(P2)))
assert np.isclose(np.trace(np.kron(A2, P2)), np.trace(A2) * np.trace(P2))
assert np.isclose(np.linalg.det(np.kron(A2, P2)), np.linalg.det(A2) ** 2 * np.linalg.det(P2) ** 3)
S1, S2 = A2 @ A2.T, P2 @ P2.T
ev = np.sort(np.linalg.eigvalsh(np.kron(S1, S2)))
assert np.allclose(ev, np.sort(np.outer(np.linalg.eigvalsh(S1), np.linalg.eigvalsh(S2)).ravel()))
assert np.isclose(np.trace(A2.T @ B2), vec(A2) @ vec(B2))
G = rng.normal(size=(4, 6)); G[:, 3:] = G[:, :3] @ rng.normal(size=(3, 3))          # rank 3
H = rng.normal(size=(2, 2)); H[:, 1] = 2 * H[:, 0]                              # rank 1
assert np.linalg.matrix_rank(np.kron(G, H)) == 3 * 1
assert np.allclose(np.kron(G, H) @ np.kron(np.linalg.pinv(G), np.linalg.pinv(H)) @ np.kron(G, H),
                   np.kron(G, H))
# (b) with a g-inverse that is not Moore-Penrose: G^- = G^+ + (I - G^+ G) Z
Gg = np.linalg.pinv(G) + (np.eye(6) - np.linalg.pinv(G) @ G) @ rng.normal(size=(6, 4))
Hg = np.linalg.pinv(H) + (np.eye(2) - np.linalg.pinv(H) @ H) @ rng.normal(size=(2, 2))
assert np.allclose(G @ Gg @ G, G) and not np.allclose(Gg, np.linalg.pinv(G))
assert np.allclose(np.kron(G, H) @ np.kron(Gg, Hg) @ np.kron(G, H), np.kron(G, H))
# (e) orthogonal, symmetric, idempotent and nonnegative definite factors
Q1, Q2 = np.linalg.qr(A2)[0], np.linalg.qr(P2)[0]
assert np.allclose(np.kron(Q1, Q2).T @ np.kron(Q1, Q2), np.eye(6))
Sym = np.kron(A2 + A2.T, P2 + P2.T)
assert np.allclose(Sym, Sym.T)
Jb3, Jb2 = np.ones((3, 3)) / 3, np.ones((2, 2)) / 2
Idem = np.kron(Jb3, np.eye(2) - Jb2)
assert np.allclose(Idem @ Idem, Idem)
N1 = G.T @ G                                             # singular, nonnegative definite
assert np.linalg.eigvalsh(np.kron(N1, S2)).min() > -1e-9

# projections of the two-way layout are Kronecker products of averaging matrices
Jbar = lambda k: np.ones((k, k)) / k
M = X @ np.linalg.pinv(X.T @ X) @ X.T
M_kron = (np.kron(np.eye(a), Jbar(b)) + np.kron(Jbar(a), np.eye(b)) - np.kron(Jbar(a), Jbar(b)))
assert np.allclose(M, M_kron)
assert np.isclose(np.trace(M), a + b - 1)

gen = Generated("ch01", "kronecker", prefix="kron")
gen.int("a", a)
gen.int("b", b)
gen.int("rank", a + b - 1)
gen.write()
