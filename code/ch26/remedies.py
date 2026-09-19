"""Chapter 26, Section 5: remedies, on the two-regressor design of Section 1.

(a) One well-placed new observation along the weak direction, against more observations that
    repeat the pattern of the old ones (Sherman-Morrison and the eigenvalue bound).
(b) Dropping x2: the short regression's slope and its bias (the omitted-variable formula).
"""
import numpy as np

from regbook import Generated

gen = Generated("ch26", "remedies")

# <<design>>
import numpy as np
rng = np.random.default_rng(2601)                     # the design of Section 26.1
n = 30
x1 = rng.normal(size=n)
x2 = x1 + 0.15 * rng.normal(size=n)
X = np.column_stack([np.ones(n), x1 - x1.mean(), x2 - x2.mean()])
G = np.linalg.inv(X.T @ X)
lam, V = np.linalg.eigh(X.T @ X)                      # ascending eigenvalues
v_weak = V[:, 0]                                      # the weak direction (it has no intercept part)
# <</design>>

# <<new-point>>
c = 1.5
x_new = np.r_[1.0, c * v_weak[1:] / np.linalg.norm(v_weak[1:])]   # breaks the pattern
print("new point (centred x1, x2):", np.round(x_new[1:], 3))
G_new = np.linalg.inv(X.T @ X + np.outer(x_new, x_new))
print("sd of b1 before / after:", round(np.sqrt(G[1, 1]), 3), round(np.sqrt(G_new[1, 1]), 3))
S_old = X[:, 1:].T @ X[:, 1:]                         # slope block: the intercept is orthogonal
xs = x_new[1:]
S_new = S_old + n / (n + 1) * np.outer(xs, xs)        # centred slope block after recentring
print("smallest slope eigenvalue before / after:", round(np.linalg.eigvalsh(S_old)[0], 3),
      round(np.linalg.eigvalsh(S_new)[0], 3))
# repeating the old pattern k times multiplies X'X by k: Var(b1) falls by the factor 1/k
k_equiv = G[1, 1] / G_new[1, 1]
print("equivalent number of copies of the whole data set:", round(k_equiv, 1))
# <</new-point>>

assert abs(v_weak[0]) < 1e-12
# Sherman-Morrison: Var reduction (a'Gx)^2/(1 + x'Gx)
a = np.array([0, 1.0, 0])
assert np.isclose(G[1, 1] - G_new[1, 1], (a @ G @ x_new) ** 2 / (1 + x_new @ G @ x_new))
# The slope estimates are governed by the centred slope block (the Schur complement of the
# intercept), whose inverse is the slope block of G_new. Adding the case (1, xs) turns it into
# S + n/(n+1) xs xs': recentring about the new mean costs the factor n/(n+1).
assert np.allclose(np.linalg.inv(S_new), G_new[1:, 1:])
lam_s = np.linalg.eigvalsh(S_old)
lam_s_new = np.linalg.eigvalsh(S_new)
# xs is along the weak eigenvector of S_old, so the proposition applies with c^2 n/(n+1)
assert np.isclose(lam_s_new[0], min(lam_s[0] + c**2 * n / (n + 1), lam_s[1]))
# random directions with the same length never do better on the centred block
for _ in range(200):
    u = rng.normal(size=2); u *= c / np.linalg.norm(u)
    assert np.linalg.eigvalsh(S_old + n / (n + 1) * np.outer(u, u))[0] <= lam_s_new[0] + 1e-12
# the new point is inside the range of each regressor separately
assert np.all(np.abs(xs) < np.abs(X[:, 1:]).max(axis=0))

# <<drop>>
x1c, x2c = X[:, 1], X[:, 2]
pi = (x1c @ x2c) / (x1c @ x1c)                        # slope of x2 on x1
sd_short = 1 / np.sqrt(x1c @ x1c)                     # sd of the short slope, units of sigma
print("slope of x2 on x1:", round(pi, 4), "  sd of the short slope:", round(sd_short, 3))
# <</drop>>
assert 0.95 < pi < 1.05

gen.num("c", c, 1)
gen.num("xnew1", x_new[1], 3)
gen.num("xnew2", x_new[2], 3)
gen.num("sd_b1_before", np.sqrt(G[1, 1]), 3)
gen.num("sd_b1_after", np.sqrt(G_new[1, 1]), 3)
gen.num("lam_before", lam_s[0], 3)
gen.num("lam_after", lam_s_new[0], 3)
gen.num("k_equiv", k_equiv, 1)
gen.num("pi", pi, 4)
gen.num("sd_short", sd_short, 3)
gen.write()
