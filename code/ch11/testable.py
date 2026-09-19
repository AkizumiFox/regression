"""Chapter 11, Section 3: testable hypotheses in a rank-deficient model.

The state murder data of Section 11.1 with one indicator for *every* Census region as well as
an intercept, so that X is 50 x 8 of rank 7:
    beta = (mu, alpha_N, alpha_M, alpha_S, alpha_W, b_poverty, b_single, b_urban).
Hypotheses Lambda' beta = d are tested with the generalized-inverse formula and, for comparison,
by fitting the reduced model with an offset. Two different generalized inverses give the same answers.
"""
import numpy as np
import statsmodels.api as sm
from scipy import linalg, stats

from regbook import Generated

# <<setup>>
data = sm.datasets.statecrime.load_pandas().data
data.index = data.index.str.strip()
data = data.drop(index="District of Columbia")
codes = "SWWSWWNSSSWWMMMMSSNSNMMSMWMWNNWNSMMSWNNSMSSWNSWSMW"
region = np.array(list(codes))
y = data["murder"].to_numpy()
n = len(y)
D = np.column_stack([(region == g).astype(float) for g in "NMSW"])   # all four regions
X = np.column_stack([np.ones(n), D, data["poverty"], data["single"], data["urban"]])
r = np.linalg.matrix_rank(X)                                          # 7, not 8
# <</setup>>

assert X.shape == (n, 8) and r == 7

# <<general>>
G = np.linalg.pinv(X.T @ X)                     # one generalized inverse of X'X
b = G @ X.T @ y                                 # one least squares solution
s2 = np.sum((y - X @ b) ** 2) / (n - r)

def glh_test(Lam, d):
    """F test of the testable hypothesis Lam' beta = d, via any g-inverse of Lam' G Lam."""
    u = Lam.T @ b - d
    W = Lam.T @ G @ Lam
    q = np.linalg.matrix_rank(Lam)
    ss_h = u @ np.linalg.pinv(W) @ u
    F = ss_h / q / s2
    return ss_h, q, F, stats.f.sf(F, q, n - r)

def contrast(pairs):
    """Columns e_i - e_j on the region coefficients (positions 1..4 = N, M, S, W)."""
    L = np.zeros((8, len(pairs)))
    for k, (i, j) in enumerate(pairs):
        L[1 + i, k], L[1 + j, k] = 1.0, -1.0
    return L

all_pairs = contrast([(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)])   # rank 3, redundant
three_equal = contrast([(1, 2), (2, 3)])                                 # M = S = W
for name, L in (("four regions equal", all_pairs), ("M = S = W", three_equal)):
    ss_h, q, F, p = glh_test(L, np.zeros(L.shape[1]))
    print(f"{name:20s} SS_H = {ss_h:.3f}  q = {q}  F = {F:.3f}  p = {p:.4f}")
# <</general>>

ss_all, q_all, F_all, p_all = glh_test(all_pairs, np.zeros(6))
ss_3, q_3, F_3, p_3 = glh_test(three_equal, np.zeros(2))
assert q_all == 3 and np.isclose(F_all, 2.440, atol=5e-4)        # the test of Section 11.1

# estimability: every column of Lambda lies in C(X') = row space of X
def estimable(L):
    return np.allclose(X.T @ np.linalg.lstsq(X.T, L, rcond=None)[0], L)

assert estimable(all_pairs) and estimable(three_equal)

# a second generalized inverse: invert X'X with the alpha_W row and column deleted
keep = [0, 1, 2, 3, 5, 6, 7]
G2 = np.zeros((8, 8))
G2[np.ix_(keep, keep)] = np.linalg.inv((X.T @ X)[np.ix_(keep, keep)])
assert np.allclose(X.T @ X @ G2 @ X.T @ X, X.T @ X)
b2 = G2 @ X.T @ y
assert not np.allclose(b, b2) and np.allclose(X @ b, X @ b2)
for L in (all_pairs, three_equal):
    assert np.allclose(L.T @ b, L.T @ b2)
    assert np.allclose(L.T @ G @ L, L.T @ G2 @ L)


# reduced model with an offset: Lambda' beta = d  <=>  beta = b0 + U gamma, C(U) = N(Lambda')
def sse(Z, target):
    coef, *_ = np.linalg.lstsq(Z, target, rcond=None)
    return np.sum((target - Z @ coef) ** 2)


sse_full = sse(X, y)
for L, ss_h in ((all_pairs, ss_all), (three_equal, ss_3)):
    U = linalg.null_space(L.T)
    assert np.isclose(sse(X @ U, y) - sse_full, ss_h)

# <<offset>>
lam = np.array([0, 1, -1 / 3, -1 / 3, -1 / 3, 0, 0, 0])   # Northeast minus the other three
d = -1.0                                                  # hypothesized difference
b0 = np.array([0, -1.0, 0, 0, 0, 0, 0, 0])                 # one solution of lam' b0 = d
U = linalg.null_space(lam[None, :])                       # beta - b0 ranges over C(U)
sse_H = sse(X @ U, y - X @ b0)                            # reduced model, offset X b0
est = lam @ b
se = np.sqrt(s2 * lam @ G @ lam)
print(f"estimate {est:.3f} (se {se:.3f}); SSE_H - SSE = {sse_H - sse_full:.3f}")
print(f"formula: {(est - d) ** 2 / (lam @ G @ lam):.3f}")
# <</offset>>

assert np.isclose(sse_H - sse_full, (est - d) ** 2 / (lam @ G @ lam))
ss_shift, q_shift, F_shift, p_shift = glh_test(lam[:, None], np.array([d]))
assert q_shift == 1 and np.isclose(F_shift, ((est - d) / se) ** 2)
# the offset may be any solution of lam' b0 = d
b0_other = np.array([0, 0, 1.0, 1.0, 1.0, 0, 0, 0])
assert np.isclose(lam @ b0_other, d)
assert np.isclose(sse(X @ U, y - X @ b0_other), sse_H)

# a nonestimable constraint, alpha_N = 0, does not restrict the mean at all
e_N = np.eye(8)[:, [1]]
assert not estimable(e_N)
U_N = linalg.null_space(e_N.T)
assert np.isclose(sse(X @ U_N, y), sse_full)
assert np.linalg.matrix_rank(X @ U_N) == r

# <<restricted>>
def restricted(Lam, d, Ginv):
    """Least squares under the testable constraint Lam' b = d (Lam of full column rank)."""
    bh = Ginv @ X.T @ y
    W = Lam.T @ Ginv @ Lam
    return bh - Ginv @ Lam @ np.linalg.solve(W, Lam.T @ bh - d)

bH = restricted(three_equal, np.zeros(2), G)
bH2 = restricted(three_equal, np.zeros(2), G2)
print("constraint holds:", np.allclose(three_equal.T @ bH, 0))
print("fitted values agree:", np.allclose(X @ bH, X @ bH2))
# <</restricted>>

assert np.allclose(three_equal.T @ bH, 0) and np.allclose(three_equal.T @ bH2, 0)
assert np.allclose(X @ bH, X @ bH2)
# the restricted fit is the fit of the model with Midwest, South and West merged
X_merged = np.column_stack([np.ones(n), D[:, 0], data["poverty"], data["single"], data["urban"]])
coef_m, *_ = np.linalg.lstsq(X_merged, y, rcond=None)
assert np.allclose(X @ bH, X_merged @ coef_m)
assert np.isclose(np.sum((y - X @ bH) ** 2) - sse_full, ss_3)
north_gap_H = bH[1] - bH[2]                    # estimable under the constraint
assert np.isclose(north_gap_H, coef_m[1])

gen = Generated("ch11", "testable")
gen.int("p", X.shape[1])
gen.int("r", r)
gen.num("s2", s2, 3)
gen.num("ss_all", ss_all, 2)
gen.num("F_all", F_all, 3)
gen.num("ss_3", ss_3, 3)
gen.num("F_3", F_3, 3)
gen.num("p_3", p_3, 3)
gen.num("est", est, 3)
gen.num("se", se, 3)
t0 = est / se
gen.num("t0", t0, 3)
gen.num("p0", 2 * stats.t.sf(abs(t0), n - r), 3)
gen.num("ss_shift", ss_shift, 3)
gen.num("F_shift", F_shift, 3)
gen.num("p_shift", p_shift, 3)
gen.num("north_gap_H", north_gap_H, 3)
gen.num("b_single_H", bH[6], 3)
gen.num("b_single", b[6], 3)
gen.write()
