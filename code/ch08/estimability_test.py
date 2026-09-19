"""Chapter 8, Section 3: testing estimability, exactly and in floating point.

Part 1 uses a one-way layout with an intercept (four groups) and checks the
characterizations of estimability for several coefficient vectors, with three different
generalized inverses. Part 2 shows that in a nearly rank-deficient design, estimability
becomes a matter of degree: the variance of a "barely estimable" function grows like 1/delta^2.
"""
import matplotlib.pyplot as plt
import numpy as np

from regbook import COLORS, Generated, figure_path, use_book_style

# <<layout>>
import numpy as np

groups = np.repeat([0, 1, 2, 3], [4, 6, 3, 5])        # four groups, sizes 4, 6, 3, 5
X = np.column_stack([np.ones(groups.size), np.eye(4)[groups]])  # [1, z1, z2, z3, z4]
A = X.T @ X
# <</layout>>


# <<tests>>
def is_estimable(lam, X, tol=None):
    """lambda is estimable iff it is orthogonal to N(X); use the SVD of X."""
    U, s, Vt = np.linalg.svd(X)
    if tol is None:
        tol = max(X.shape) * np.finfo(float).eps * s[0]
    r = int(np.sum(s > tol))
    V0 = Vt[r:].T                                  # orthonormal basis of N(X)
    return np.linalg.norm(V0.T @ lam) <= 1e-8 * max(1.0, np.linalg.norm(lam))

G_mp = np.linalg.pinv(A)                            # Moore-Penrose inverse of X^T X
G_ref = np.zeros((5, 5))
G_ref[1:, 1:] = np.linalg.inv(A[1:, 1:])            # delete the intercept row and column
rng = np.random.default_rng(3)
G_odd = G_mp + (np.eye(5) - G_mp @ A) @ rng.normal(size=(5, 5))   # non-symmetric g-inverse

candidates = {"mu + alpha_1": [1, 1, 0, 0, 0],
              "alpha_1 - alpha_2": [0, 1, -1, 0, 0],
              "alpha_1 + alpha_2 - 2 alpha_3": [0, 1, 1, -2, 0],
              "alpha_1": [0, 1, 0, 0, 0],
              "alpha_1 + ... + alpha_4": [0, 1, 1, 1, 1]}
for name, lam in candidates.items():
    lam = np.array(lam, dtype=float)
    H_checks = [np.allclose(lam @ G @ A, lam) for G in (G_mp, G_ref, G_odd)]
    variances = [lam @ G @ lam for G in (G_mp, G_ref, G_odd)]
    print(f"{name:30s} SVD test: {is_estimable(lam, X)!s:5s}",
          f" lambda'GX'X = lambda': {H_checks}",
          f" lambda'G lambda: {np.round(variances, 3)}")
# <</tests>>

est = {k: is_estimable(np.array(v, float), X) for k, v in candidates.items()}
assert est == {"mu + alpha_1": True, "alpha_1 - alpha_2": True,
               "alpha_1 + alpha_2 - 2 alpha_3": True, "alpha_1": False,
               "alpha_1 + ... + alpha_4": False}
for G in (G_mp, G_ref, G_odd):
    assert np.allclose(A @ G @ A, A)
var_table = {}
for name, lam in candidates.items():
    lam = np.array(lam, float)
    checks = [np.allclose(lam @ G @ A, lam) for G in (G_mp, G_ref, G_odd)]
    assert all(checks) == est[name] and any(checks) == est[name]   # "some" iff "every"
    v = [lam @ G @ lam for G in (G_mp, G_ref, G_odd)]
    if est[name]:
        assert np.allclose(v, v[0])                 # invariant for estimable functions
    else:
        assert np.ptp(v) > 1e-3                     # meaningless otherwise
    var_table[name] = v
# variance of alpha_1 - alpha_2 is 1/4 + 1/6
assert np.isclose(var_table["alpha_1 - alpha_2"][0], 1 / 4 + 1 / 6)
assert np.isclose(var_table["alpha_1 + alpha_2 - 2 alpha_3"][0], 1 / 4 + 1 / 6 + 4 / 3)
# the estimable space has dimension rank(X) = 4: H = G X^T X has rank 4 and trace 4
H = G_odd @ A
assert np.linalg.matrix_rank(H) == 4 and np.isclose(np.trace(H), 4) and np.allclose(H @ H, H)

gen = Generated("ch08", "estimability_test", prefix="et")
gen.num("v12", var_table["alpha_1 - alpha_2"][0], 3)
gen.num("v123", var_table["alpha_1 + alpha_2 - 2 alpha_3"][0], 3)
gen.num("vmu1", var_table["mu + alpha_1"][0], 3)
for key, i in zip(["mp", "ref", "odd"], range(3)):
    gen.num(f"a1:{key}", var_table["alpha_1"][i], 3)
    gen.num(f"asum:{key}", var_table["alpha_1 + ... + alpha_4"][i], 3)

# ---- Part 2: near rank deficiency ----------------------------------------------
# <<near>>
t = np.linspace(-1, 1, 12)
w = np.cos(np.pi * t)                                # a fixed direction to perturb along
def near_design(delta):
    """Columns 1, t and t + delta*w: exactly rank deficient at delta = 0."""
    return np.column_stack([np.ones_like(t), t, t + delta * w])

for delta in [1.0, 1e-2, 1e-4, 1e-8, 0.0]:
    Xd = near_design(delta)
    s = np.linalg.svd(Xd, compute_uv=False)
    G = np.linalg.pinv(Xd.T @ Xd)
    v_sum = np.array([0, 1, 1]) @ G @ [0, 1, 1]     # beta_1 + beta_2
    v_one = np.array([0, 1, 0]) @ G @ [0, 1, 0]     # beta_1 alone
    print(f"delta={delta:7.0e}  smallest singular value {s[-1]:8.1e}",
          f" Var(b1+b2)/s2 = {v_sum:8.3f}  Var(b1)/s2 = {v_one:9.2e}",
          f" beta_1 estimable? {is_estimable(np.array([0., 1, 0]), Xd)}")
# <</near>>

deltas = np.logspace(-6, 0, 40)
v_sum, v_one, v_diff = [], [], []
for d in deltas:
    Xd = near_design(d)
    Ci = np.linalg.inv(Xd.T @ Xd)
    v_sum.append(np.array([0, 1, 1]) @ Ci @ [0, 1, 1])
    v_one.append(Ci[1, 1])
    v_diff.append(np.array([0, 1, -1]) @ Ci @ [0, 1, -1])
v_sum, v_one, v_diff = map(np.array, (v_sum, v_one, v_diff))
# Var(b1) grows like 1/delta^2; Var(b1 + b2) stays bounded
slope = np.polyfit(np.log10(deltas[:20]), np.log10(v_one[:20]), 1)[0]
assert abs(slope + 2) < 0.05
assert v_sum.max() < 2 * v_sum.min() + 1
assert is_estimable(np.array([0., 1, 0]), near_design(1e-4))
assert not is_estimable(np.array([0., 1, 0]), near_design(0.0))
assert is_estimable(np.array([0., 1, 1]), near_design(0.0))
# at delta = 1e-8 the SVD of X still sees rank 3, but pinv(X^T X) sees rank 2:
# forming X^T X squares the smallest singular value, 1.8e-8 -> 3e-16, below its tolerance
X8 = near_design(1e-8)
assert is_estimable(np.array([0., 1, 0]), X8)
assert np.linalg.matrix_rank(X8) == 3 and np.linalg.matrix_rank(X8.T @ X8) == 2
v8 = np.linalg.pinv(X8.T @ X8)[1, 1]
assert v8 < 1                                        # the huge true variance has vanished
gen.num("near:v8", v8, 4)
_, s8, Vt8 = np.linalg.svd(X8, full_matrices=False)
true8 = np.sum((Vt8[:, 1] / s8) ** 2)                 # Var(b1)/sigma^2 from the SVD of X
assert true8 > 1e14
gen.int("near:ratio_exp", int(np.floor(np.log10(true8 / v8))))
X4 = near_design(1e-4)
C4 = np.linalg.inv(X4.T @ X4)
gen.num("near:vsum", C4[1:, 1:].sum(), 3)
gen.num("near:vone", C4[1, 1], 2, sci=True)
gen.num("near:smin", np.linalg.svd(X4, compute_uv=False)[-1], 1, sci=True)
gen.num("near:slope", slope, 2)
gen.write()

use_book_style()
fig, ax = plt.subplots(figsize=(3.9, 2.6))
ax.loglog(deltas, v_one, color=COLORS["second"], label=r"$\hat\beta_1$ alone")
ax.loglog(deltas, v_diff, "--", color=COLORS["thread"], label=r"$\hat\beta_1-\hat\beta_2$")
ax.loglog(deltas, v_sum, color=COLORS["accent"], label=r"$\hat\beta_1+\hat\beta_2$")
ax.set_xlabel(r"perturbation $\delta$")
ax.set_ylabel(r"variance $/\,\sigma^2$")
ax.legend(frameon=False, loc="upper right")
fig.savefig(figure_path("ch08", "near_estimable"))
