"""Chapter 26, Section 1: what collinearity does to the variances of linear functions.

(a) Two regressors with sample correlation about 0.99: the canonical coordinates of the
    slope pair, the variances of the slopes, their sum and their difference, and a simulation
    of the sampling distribution of (beta_1_hat, beta_2_hat) (figure).
(b) A family of designs [X1, X1 c + delta z] approaching rank deficiency: a linear function
    keeps a bounded variance as delta -> 0 iff it is estimable in the limiting design.
"""
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import norm

from regbook import COLORS, Generated, figure_path, use_book_style

gen = Generated("ch26", "variance")

# <<design>>
import numpy as np
rng = np.random.default_rng(2601)
n, sigma = 30, 1.0
x1 = rng.normal(size=n)
x2 = x1 + 0.15 * rng.normal(size=n)                  # x2 nearly repeats x1
x1c, x2c = x1 - x1.mean(), x2 - x2.mean()             # centring makes the intercept orthogonal
X = np.column_stack([np.ones(n), x1c, x2c])
G = np.linalg.inv(X.T @ X)                            # Cov(beta_hat) / sigma^2

lam, V = np.linalg.eigh(X[:, 1:].T @ X[:, 1:])        # canonical coordinates of the slopes
print("correlation of x1 and x2:", round(np.corrcoef(x1, x2)[0, 1], 4))
print("eigenvalues:", np.round(lam, 3), " weak direction:", np.round(V[:, 0], 3))

def sd(a):
    """Standard deviation of a' beta_hat, in units of sigma."""
    return np.sqrt(a @ G @ a)

print("sd of b1, b2   :", round(sd(np.array([0, 1, 0])), 3), round(sd(np.array([0, 0, 1])), 3))
print("sd of b1 + b2  :", round(sd(np.array([0, 1, 1])), 3))
print("sd of b1 - b2  :", round(sd(np.array([0, 1, -1])), 3))
# <</design>>

r12 = np.corrcoef(x1, x2)[0, 1]
lam_small, lam_large = lam[0], lam[1]
v_small = V[:, 0] * np.sign(V[0, 0])
# canonical coordinates: v' beta_hat has variance sigma^2 / lambda, and they are uncorrelated
C = V.T @ G[1:, 1:] @ V
assert np.allclose(C, np.diag(1 / lam))
# extreme variances over unit vectors, and the kappa^2 ratio
assert np.isclose(sd(np.r_[0, v_small]) ** 2, 1 / lam_small)
kappa2 = lam_large / lam_small
assert np.isclose(np.linalg.cond(X[:, 1:]) ** 2, kappa2)
# FWL form: Var(b_j) = sigma^2 / ||(I - M_(j)) x_j||^2 >= sigma^2 / lambda_min
res = x1c - x2c * (x2c @ x1c) / (x2c @ x2c)
assert np.isclose(G[1, 1], 1 / (res @ res)) and res @ res >= lam_small
sd_b1, sd_b2 = sd(np.array([0, 1, 0])), sd(np.array([0, 0, 1]))
sd_sum, sd_diff = sd(np.array([0, 1, 1])), sd(np.array([0, 1, -1]))

# <<simulate>>
beta = np.array([0.0, 1.0, 1.0])                      # equal true slopes
reps = 2000
Y = X @ beta + sigma * rng.normal(size=(reps, n))     # 2000 replicate responses, same design
B = np.linalg.solve(X.T @ X, X.T @ Y.T).T             # one row of estimates per replicate
print("share of replicates with b2 < 0:", np.mean(B[:, 2] < 0))
print("share with b1 < 0 or b2 < 0:", np.mean((B[:, 1:] < 0).any(axis=1)))
print("range of b1 + b2:", np.round([B[:, 1:].sum(1).min(), B[:, 1:].sum(1).max()], 3))
# <</simulate>>

share_neg = np.mean(B[:, 2] < 0)
sums = B[:, 1:].sum(1)
assert 0.1 < share_neg < 0.5
share_either = np.mean((B[:, 1:] < 0).any(axis=1))
assert share_either > 0.5 and share_either > 1.9 * share_neg   # the two events are nearly disjoint
assert abs(np.std(B[:, 1]) / sd_b1 - 1) < 0.05
assert abs(np.std(sums) / sd_sum - 1) < 0.05

# ---- (b) near-estimability: X_delta = [X1, X1 c + delta z] ------------------------------
t = np.linspace(-1, 1, 12)
X1 = np.column_stack([np.ones_like(t), t])
c = np.array([0.0, 1.0])                              # limiting design [1, t, t]
w = np.cos(np.pi * t)
z = w - X1 @ np.linalg.lstsq(X1, w, rcond=None)[0]
z /= np.linalg.norm(z)                                # unit vector orthogonal to C(X1)
G1 = np.linalg.inv(X1.T @ X1)
for delta in [1.0, 1e-1, 1e-3]:
    Xd = np.column_stack([X1, X1 @ c + delta * z])
    Gd = np.linalg.inv(Xd.T @ Xd)
    for a in [np.array([0, 1, 1.0]), np.array([0, 1, 0.0]), np.array([1, 0.5, -0.5])]:
        a1, a2 = a[:2], a[2]
        formula = a1 @ G1 @ a1 + (a2 - c @ a1) ** 2 / delta**2
        assert np.isclose(a @ Gd @ a, formula, rtol=1e-8)
# a in C(X0'), X0 = [X1, X1 c]: variance independent of delta and equal to the limiting BLUE's
a = np.array([0.0, 1.0, 1.0]) * 0.5                   # (beta_1 + beta_2)/2, estimable in the limit
X0 = np.column_stack([X1, X1 @ c])
w0 = np.linalg.lstsq(X0.T, a, rcond=None)[0]
assert np.allclose(X0.T @ w0, a)
M0 = X0 @ np.linalg.pinv(X0)
var_limit = w0 @ M0 @ w0
assert np.isclose(var_limit, a[:2] @ G1 @ a[:2])

gen.int("n", n)
gen.num("r12", r12, 4)
gen.num("lam_small", lam_small, 3)
gen.num("lam_small4", lam_small, 4)
gen.num("lam_large", lam_large, 2)
gen.num("v1", abs(v_small[0]), 3)
gen.num("kappa2", kappa2, 1)
gen.num("sd_b1", sd_b1, 3)
gen.num("sd_b2", sd_b2, 3)
gen.num("sd_sum", sd_sum, 3)
gen.num("sd_diff", sd_diff, 3)
gen.num("sd_orth", 1 / np.sqrt(x1c @ x1c), 3)
gen.num("share_neg", 100 * share_neg, 1)
gen.num("share_either", 100 * share_either, 1)
var_diff, up_diff, low_diff = sd_diff**2, 2 / lam_small, 2 / lam_large
assert low_diff < var_diff < up_diff
gen.num("var_diff", var_diff, 3)
gen.num("up_diff", up_diff, 3)
gen.num("low_diff", low_diff, 4)
gen.num("vp_weight", (np.array([1, -1]) @ V[:, 0]) ** 2, 4)
gen.num("sum_min", sums.min(), 3)
gen.num("sum_max", sums.max(), 3)
gen.num("b2_min", B[:, 2].min(), 3)
gen.num("b2_max", B[:, 2].max(), 3)
gen.int("reps", reps)
p_neg = norm.cdf(-1 / sd_b2)                          # P(b2 < 0) under normal errors, beta_2 = 1
assert abs(share_neg - p_neg) < 3 * np.sqrt(p_neg * (1 - p_neg) / reps)
gen.num("p_neg", p_neg, 3)
gen.num("z_neg", 1 / sd_b2, 3)
gen.num("vif1", 1 / (1 - r12**2), 1)
gen.num("sd_ratio", sd_b1 / (1 / np.sqrt(x1c @ x1c)), 1)
gen.write()

# ---- figure ------------------------------------------------------------------------------
use_book_style()
fig, axes = plt.subplots(1, 2, figsize=(5.8, 2.7))
ax = axes[0]
ax.scatter(x1c, x2c, s=10, color=COLORS["accent"], linewidths=0)
ax.set_xlabel(r"$x_1$ (centred)")
ax.set_ylabel(r"$x_2$ (centred)")
ax.set_title("(a) the design")
ax.set_aspect("equal")
ax = axes[1]
ax.scatter(B[:, 1], B[:, 2], s=2, color=COLORS["accent"], alpha=0.35, linewidths=0)
ax.plot([1], [1], "o", color=COLORS["second"], ms=5, zorder=3)
g = np.linspace(-5.5, 7.5, 2)
ax.plot(g, 2 - g, color=COLORS["second"], linewidth=0.7, linestyle="--")
ax.axhline(0, color=COLORS["grid"], linewidth=0.6, zorder=0)
ax.axvline(0, color=COLORS["grid"], linewidth=0.6, zorder=0)
ax.set_xlim(-5.5, 7.5); ax.set_ylim(-5.5, 7.5)
ax.set_aspect("equal")
ax.set_xlabel(r"$\hat\beta_1$")
ax.set_ylabel(r"$\hat\beta_2$")
ax.set_title(r"(b) %d replicate estimates" % reps)
fig.tight_layout()
fig.savefig(figure_path("ch26", "canonical_directions"))
