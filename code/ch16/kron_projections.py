"""Chapter 16: projections of balanced layouts as Kronecker products of averaging and centring matrices.

Checks, by assertion, every structural claim of the chapter about balanced layouts:
(1) the two-way additive model (Section 16.1): P0 + PA + PB is the projection onto
    C([1, Z_A, Z_B]), and the residual projection splits into interaction and pure error;
(2) the balanced two-way layout with interaction (Section 16.3): five mutually orthogonal
    projections summing to I, with ranks 1, a-1, b-1, (a-1)(b-1), ab(m-1);
(3) the general balanced factorial (Section 16.4): 2^k projections plus error, ranks as
    products, and each P_S y is the inclusion-exclusion combination of marginal means;
(4) the nested layout (Section 16.5): I_a (x) C_b (x) Jbar_m = P_B + P_AB, and the nested
    projections are those of the model with indicators of B within A.
"""
import itertools

import numpy as np

from regbook import Generated

rng = np.random.default_rng(1601)


# <<kron>>
def Jbar(k):
    return np.full((k, k), 1.0 / k)          # averaging: replaces a k-vector by its mean


def Cen(k):
    return np.eye(k) - Jbar(k)               # centring: subtracts the mean


def kron_all(mats):
    out = np.ones((1, 1))
    for A in mats:
        out = np.kron(out, A)
    return out


def factorial_projections(levels, m):
    """P_S for every subset S of the factors (a 0/1 tuple), and the error projection."""
    P = {}
    for S in itertools.product([0, 1], repeat=len(levels)):
        P[S] = kron_all([Cen(l) if s else Jbar(l) for l, s in zip(levels, S)] + [Jbar(m)])
    P["error"] = kron_all([np.eye(l) for l in levels] + [Cen(m)])
    return P


def rank(A):
    return int(round(np.trace(A)))           # rank = trace for a projection
# <</kron>>


def proj(Z):
    """Orthogonal projection onto C(Z), from an SVD (Z need not have full rank)."""
    U, s, _ = np.linalg.svd(Z, full_matrices=False)
    U = U[:, s > 1e-10 * s[0]]
    return U @ U.T


def is_projection(P):
    return np.allclose(P, P.T) and np.allclose(P @ P, P)


# ---- (1) and (2): the two-way layout, a = 3, b = 4, m = 2 ---------------------------
# <<twoway>>
a, b, m = 3, 4, 2
P = factorial_projections([a, b], m)
names = {(0, 0): "mean", (1, 0): "A", (0, 1): "B", (1, 1): "AB", "error": "error"}
for key, name in names.items():
    print(f"{name:6s} rank {rank(P[key])}")
keys = list(names)
worst = max(np.abs(P[s] @ P[t]).max() for s, t in itertools.combinations(keys, 2))
print("largest entry of any product P_s P_t, s != t:", f"{worst:.1e}")
print("sum is the identity:", np.allclose(sum(P[k] for k in keys), np.eye(a * b * m)))
# <</twoway>>

n = a * b * m
for key in keys:
    assert is_projection(P[key])
assert worst < 1e-12
assert [rank(P[k]) for k in keys] == [1, a - 1, b - 1, (a - 1) * (b - 1), a * b * (m - 1)]

# indicator matrices, observations ordered with k fastest, then j, then i
row = np.repeat(np.arange(a), b * m)
col = np.tile(np.repeat(np.arange(b), m), a)
ZA, ZB = np.eye(a)[row], np.eye(b)[col]
ZAB = np.eye(a * b)[row * b + col]
one = np.ones((n, 1))
M_add = proj(np.column_stack([one, ZA, ZB]))
M_int = proj(np.column_stack([one, ZA, ZB, ZAB]))
assert np.allclose(M_add, P[(0, 0)] + P[(1, 0)] + P[(0, 1)])
assert np.allclose(M_int, M_add + P[(1, 1)])
assert np.allclose(np.eye(n) - M_int, P["error"])
# the additive residual projection = interaction + pure error
assert np.allclose(np.eye(n) - M_add, P[(1, 1)] + P["error"])
# action on data: cell, row, column, grand means
y = rng.normal(size=n)
Y = y.reshape(a, b, m)
cell, rmean, cmean, g = Y.mean(2), Y.mean((1, 2)), Y.mean((0, 2)), Y.mean()
ab_expected = (cell - rmean[:, None] - cmean[None, :] + g)
assert np.allclose((P[(1, 1)] @ y).reshape(a, b, m), ab_expected[:, :, None])
assert np.allclose((P[(1, 0)] @ y).reshape(a, b, m), (rmean - g)[:, None, None] * np.ones((a, b, m)))
assert np.allclose((P["error"] @ y).reshape(a, b, m), Y - cell[:, :, None])
# with one observation per cell the additive residual space is the interaction space
P1 = factorial_projections([a, b], 1)
assert np.allclose(np.eye(a * b) - proj(np.column_stack([np.ones(a * b), np.eye(a)[np.repeat(np.arange(a), b)],
                                                          np.eye(b)[np.tile(np.arange(b), a)]])), P1[(1, 1)])
assert rank(P1["error"]) == 0

# ---- (3): three factors, levels 2, 3, 4, m = 2 ----------------------------------------
# <<threeway>>
levels, m3 = [2, 3, 4], 2
P3 = factorial_projections(levels, m3)
for S in P3:
    if S == "error":
        continue
    label = "".join(f for f, s in zip("ABC", S) if s) or "mean"
    print(f"{label:5s} rank {rank(P3[S]):2d}   predicted {int(np.prod([l - 1 for l, s in zip(levels, S) if s]))}")
print("error rank", rank(P3["error"]), "  n =", int(np.prod(levels)) * m3)
# <</threeway>>

keys3 = list(P3)
n3 = int(np.prod(levels)) * m3
assert np.allclose(sum(P3[k] for k in keys3), np.eye(n3))
assert max(np.abs(P3[s] @ P3[t]).max() for s, t in itertools.combinations(keys3, 2)) < 1e-12
for S in keys3:
    if S != "error":
        assert rank(P3[S]) == int(np.prod([l - 1 for l, s in zip(levels, S) if s]))
assert rank(P3["error"]) == int(np.prod(levels)) * (m3 - 1)

# P_S y by inclusion-exclusion of marginal means: P_S = sum_{T subset S} (-1)^{|S|-|T|} A_T,
# where A_T averages over every factor outside T (and over replicates)
def averaging(T):
    return kron_all([np.eye(l) if t else Jbar(l) for l, t in zip(levels, T)] + [Jbar(m3)])

for S in keys3:
    if S == "error":
        continue
    total = np.zeros((n3, n3))
    for T in itertools.product([0, 1], repeat=3):
        if all(t <= s for t, s in zip(T, S)):
            total += (-1) ** (sum(S) - sum(T)) * averaging(T)
    assert np.allclose(total, P3[S])
# the space of vectors depending only on the factors in S has projection sum_{T subset S} P_T
for S in itertools.product([0, 1], repeat=3):
    total = sum(P3[T] for T in P3 if T != "error" and all(t <= s for t, s in zip(T, S)))
    assert np.allclose(total, averaging(S))

# a hierarchical model (A, B, C, AB) has projection P0 + PA + PB + PC + PAB,
# the same as the projection onto its indicator columns
y3 = rng.normal(size=n3)
idx = np.array(list(itertools.product(*[range(l) for l in levels], range(m3))))
iA, iB, iC = idx[:, 0], idx[:, 1], idx[:, 2]
X_h = np.column_stack([np.ones(n3), np.eye(2)[iA], np.eye(3)[iB], np.eye(4)[iC], np.eye(6)[iA * 3 + iB]])
M_h = proj(X_h)
assert np.allclose(M_h, P3[(0, 0, 0)] + P3[(1, 0, 0)] + P3[(0, 1, 0)] + P3[(0, 0, 1)] + P3[(1, 1, 0)])

# ---- (4): nested layout, a = 3 levels of A, b = 4 of B within each, m = 2 -----------
# <<nested>>
a, b, m = 3, 4, 2
PB_in_A = kron_all([np.eye(a), Cen(b), Jbar(m)])        # B within A
P2 = factorial_projections([a, b], m)
print("B(A) = B + AB of the crossed layout:", np.allclose(PB_in_A, P2[(0, 1)] + P2[(1, 1)]))
print("rank of B(A):", rank(PB_in_A), " = a(b-1) =", a * (b - 1))
# <</nested>>
assert np.allclose(PB_in_A, P2[(0, 1)] + P2[(1, 1)])
assert rank(PB_in_A) == a * (b - 1)
n = a * b * m
row = np.repeat(np.arange(a), b * m)
cls = np.repeat(np.arange(a * b), m)                    # b levels of B within each A level: a*b labels
M_nest = proj(np.column_stack([np.ones(n), np.eye(a)[row], np.eye(a * b)[cls]]))
assert np.allclose(M_nest, P2[(0, 0)] + P2[(1, 0)] + PB_in_A)
assert np.allclose(np.eye(n) - M_nest, P2["error"])

gen = Generated("ch16", "kron_projections")
gen.int("rank_ab", (3 - 1) * (4 - 1))
gen.int("n3", n3)
gen.int("err3", rank(P3["error"]))
gen.write()
