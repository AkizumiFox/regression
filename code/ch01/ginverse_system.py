"""Chapter 1, Sections 3-4: generalized inverses solve consistent systems.

The unknowns are the six cell entries of a 2 x 3 table; the equations prescribe the two
row totals and the three column totals. The system is consistent exactly when the row
totals and the column totals have the same grand total.
"""
from fractions import Fraction
from itertools import combinations

import numpy as np

from regbook import Generated

# <<system>>
# cells ordered (1,1), (1,2), (1,3), (2,1), (2,2), (2,3)
A = np.array([[1, 1, 1, 0, 0, 0],      # row 1 total
              [0, 0, 0, 1, 1, 1],      # row 2 total
              [1, 0, 0, 1, 0, 0],      # column 1 total
              [0, 1, 0, 0, 1, 0],      # column 2 total
              [0, 0, 1, 0, 0, 1]],     # column 3 total
             dtype=float)
b = np.array([9, 6, 4, 5, 6], dtype=float)       # both sets of totals add to 15
r = np.linalg.matrix_rank(A)
print("rank(A) =", r, " rank([A, b]) =", np.linalg.matrix_rank(np.column_stack([A, b])))
# <</system>>
m, n = A.shape
assert r == 4 and n - r == 2


# <<ginverses>>
def ginverse_from_submatrix(A, rows, cols):
    """Place the inverse of the nonsingular block A[rows, cols] at [cols, rows]."""
    G = np.zeros(A.T.shape)
    G[np.ix_(cols, rows)] = np.linalg.inv(A[np.ix_(rows, cols)])
    return G


rows, cols = [0, 1, 2, 3], [0, 1, 2, 3]                  # r independent rows and columns
G1 = ginverse_from_submatrix(A, rows, cols)
G2 = np.linalg.pinv(A)                                   # Moore-Penrose inverse
rng = np.random.default_rng(1)
Z = rng.integers(-2, 3, size=A.T.shape).astype(float)
G3 = G1 + (np.eye(n) - G1 @ A) @ Z                       # another, non-reflexive one

for G in (G1, G2, G3):
    assert np.allclose(A @ G @ A, A)                     # each is a g-inverse
    assert np.allclose(A @ G @ b, b)                     # consistency test passes
    x = G @ b
    print(np.round(x, 4), " A x - b =", np.abs(A @ x - b).max())
# <</ginverses>>

# every r x r submatrix built from r independent rows and r independent columns is nonsingular
indep = lambda S: np.linalg.matrix_rank(S) == r
n_checked = 0
for R in combinations(range(m), r):
    if not indep(A[list(R), :]):
        continue
    for Cc in combinations(range(n), r):
        if indep(A[:, list(Cc)]):
            assert abs(np.linalg.det(A[np.ix_(R, Cc)])) > 1e-9
            n_checked += 1

# <<inconsistent>>
b_bad = b.copy()
b_bad[4] = 7                                             # column totals now add to 16
for G in (G1, G2, G3):
    assert not np.allclose(A @ G @ b_bad, b_bad)         # the test detects inconsistency
print("rank([A, b_bad]) =", np.linalg.matrix_rank(np.column_stack([A, b_bad])))
# <</inconsistent>>
assert np.linalg.matrix_rank(np.column_stack([A, b_bad])) == r + 1

# general solution G b + (I - G A) z sweeps out x0 + N(A)
x1, x2, x3 = G1 @ b, G2 @ b, G3 @ b
K = np.eye(n) - G1 @ A
assert np.linalg.matrix_rank(K) == n - r and np.allclose(A @ K, 0)
for x in (x2, x3):
    z = np.linalg.lstsq(K, x - x1, rcond=None)[0]
    assert np.allclose(x1 + K @ z, x)
# minimum norm: the Moore-Penrose solution lies in C(A^T) and is shortest
assert np.allclose(A.T @ np.linalg.lstsq(A.T, x2, rcond=None)[0], x2)
for _ in range(200):
    x = x1 + K @ rng.normal(size=n)
    assert np.linalg.norm(x) >= np.linalg.norm(x2) - 1e-12
# invariance: q'x is the same for all solutions iff q in C(A^T)
q_bad = np.array([1, 0, 0, 0, 0, 0.])        # a single cell
in_row_space = lambda q: np.allclose(q @ G1 @ A, q)
assert not in_row_space(q_bad)
assert not np.isclose(q_bad @ x1, q_bad @ x2)
q_row = A[0] - A[2]                          # row-1 total minus column-1 total = x12 + x13 - x21
assert in_row_space(q_row)
assert np.isclose(q_row @ x1, q_row @ x2) and np.isclose(q_row @ x1, q_row @ x3)
# additive form of the minimum-norm solution: x_ij = r_i/3 + c_j/2 - 15/6
row_tot, col_tot = b[:2], b[2:]
additive = (row_tot[:, None] / 3 + col_tot[None, :] / 2 - 15 / 6).ravel()
assert np.allclose(x2, additive)
# G3 is a g-inverse but not reflexive: rank(G3) > r, so G3 A G3 != G3
assert np.linalg.matrix_rank(G3) > r and not np.allclose(G3 @ A @ G3, G3)
assert np.isclose(np.trace(A @ G1), r) and np.isclose(np.trace(G3 @ A), r)


def frac(v):
    f = Fraction(v).limit_denominator(100)
    s = f"{abs(f.numerator)}" if f.denominator == 1 else f"{abs(f.numerator)}/{f.denominator}"
    return ("\\ensuremath{-}" if f < 0 else "") + s


gen = Generated("ch01", "ginverse_system", prefix="gi")
for name, x in (("a", x1), ("b", x2), ("c", x3)):
    gen.text(f"x{name}", ",\\ ".join(frac(v) for v in x))
    gen.text(f"x{name}11", frac(x[0]))
    gen.num(f"norm{name}", np.linalg.norm(x), 3)
gen.text("qrow", frac(q_row @ x1))
gen.int("nsub", n_checked)
gen.write()
