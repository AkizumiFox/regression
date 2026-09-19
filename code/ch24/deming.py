"""Chapter 24, Section 3: Deming and orthogonal regression in a method-comparison study.

Sixty specimens with true concentration x ~ Uniform(2, 20). Method A reads w = x + u,
u ~ N(0, 0.8^2); method B reads y = 0.3 + 1.05 x + e, e ~ N(0, 1.0^2). Repeatability studies give
the error-variance ratio delta = 1.0^2 / 0.8^2 = 1.5625.

Checks: the Deming slope lies between the direct and reverse regression slopes; for delta = 1 it
equals the slope of the first principal axis (orthogonal regression); orthogonal regression is
not equivariant under a change of units of y; over repeated samples the Deming slope is close to
unbiased while least squares is attenuated.
"""
from regbook import Generated

gen = Generated("ch24", "deming", prefix="dem")

# <<deming>>
import numpy as np


def deming(w, y, delta):
    """Slope and intercept of Deming regression; delta = Var(error in y) / Var(error in w)."""
    S = np.cov(w, y)
    s_ww, s_wy, s_yy = S[0, 0], S[0, 1], S[1, 1]
    d = s_yy - delta * s_ww
    b = (d + np.sqrt(d ** 2 + 4 * delta * s_wy ** 2)) / (2 * s_wy)
    return b, y.mean() - b * w.mean()


rng = np.random.default_rng(2430)
x = rng.uniform(2, 20, 60)                        # true concentrations
w = x + rng.normal(0, 0.8, 60)                    # method A
y = 0.3 + 1.05 * x + rng.normal(0, 1.0, 60)       # method B
delta = 1.0 ** 2 / 0.8 ** 2

b_ls = np.polyfit(w, y, 1)[0]                     # regress y on w
b_rev = 1 / np.polyfit(y, w, 1)[0]                # regress w on y, inverted
b_dem, a_dem = deming(w, y, delta)
b_orth, _ = deming(w, y, 1.0)                     # orthogonal regression
print(f"least squares {b_ls:.4f}, Deming {b_dem:.4f} (intercept {a_dem:.3f}), "
      f"orthogonal {b_orth:.4f}, reverse {b_rev:.4f}")
# <</deming>>

assert b_ls < b_orth < b_dem < b_rev or b_ls < b_dem < b_orth < b_rev
assert b_ls < b_dem < b_rev
# orthogonal regression = first principal axis of the sample covariance matrix
vals, vecs = np.linalg.eigh(np.cov(w, y))
v = vecs[:, 1]                                    # eigenvector of the largest eigenvalue
assert np.isclose(v[1] / v[0], b_orth)
# Deming = orthogonal regression after dividing y by sqrt(delta)
b_scaled, _ = deming(w, y / np.sqrt(delta), 1.0)
assert np.isclose(b_scaled * np.sqrt(delta), b_dem)
# orthogonal regression is not equivariant under a change of units of y; Deming (with delta
# rescaled accordingly) and least squares are
b_orth10, _ = deming(w, 10 * y, 1.0)
b_dem10, _ = deming(w, 10 * y, 100 * delta)
assert np.isclose(b_dem10, 10 * b_dem)
assert abs(b_orth10 / 10 - b_orth) > 0.005
gen.num("delta", delta, 4)
gen.num("b_ls", b_ls, 4)
gen.num("b_rev", b_rev, 4)
gen.num("b_dem", b_dem, 4)
gen.num("a_dem", a_dem, 3)
gen.num("b_orth", b_orth, 4)
gen.num("b_orth10", b_orth10 / 10, 4)

# ---- repeated samples ------------------------------------------------------------------------
r = np.random.default_rng(2431)
reps = 20_000
X = r.uniform(2, 20, (reps, 60))
Wr = X + r.normal(0, 0.8, (reps, 60))
Yr = 0.3 + 1.05 * X + r.normal(0, 1.0, (reps, 60))
bl = np.array([np.polyfit(a, b, 1)[0] for a, b in zip(Wr, Yr)])
bd = np.array([deming(a, b, delta)[0] for a, b in zip(Wr, Yr)])
lam_pop = 27 / (27 + 0.64)                        # Var(Uniform(2, 20)) = 18^2 / 12 = 27
assert abs(bl.mean() - 1.05 * lam_pop) < 0.002
assert abs(bd.mean() - 1.05) < 0.002
gen.int("reps", reps)
gen.num("mean_ls", bl.mean(), 4)
gen.num("mean_dem", bd.mean(), 4)
gen.num("lam_pop", lam_pop, 4)
gen.write()
