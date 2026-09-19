"""Chapter 1, Section 11: derivative formulas checked against finite differences."""
import numpy as np

from regbook import Generated

rng = np.random.default_rng(5)
p = 4


def num_grad(f, x, h=1e-6):
    """Central-difference gradient of a scalar function of a vector or matrix."""
    g = np.zeros_like(x)
    for idx in np.ndindex(x.shape):
        e = np.zeros_like(x)
        e[idx] = h
        g[idx] = (f(x + e) - f(x - e)) / (2 * h)
    return g


# <<gradients>>
a = rng.normal(size=p)
A = rng.normal(size=(p, p))                         # not symmetric
x = rng.normal(size=p)
checks = {
    "a'x":        (num_grad(lambda v: a @ v, x), a),
    "x'Ax":       (num_grad(lambda v: v @ A @ v, x), (A + A.T) @ x),
}
X = rng.normal(size=(20, p))
y = rng.normal(size=20)
b = rng.normal(size=p)
grad_ls = -2 * X.T @ (y - X @ b)
checks["||y - Xb||^2"] = (num_grad(lambda v: np.sum((y - X @ v) ** 2), b), grad_ls)

S = rng.normal(size=(p, p)); S = S @ S.T + p * np.eye(p)   # positive definite
checks["log det S"] = (num_grad(lambda W: np.linalg.slogdet(W)[1], S), np.linalg.inv(S).T)
checks["tr(AS)"] = (num_grad(lambda W: np.trace(A @ W), S), A.T)
for name, (numeric, formula) in checks.items():
    print(f"{name:14s} max error {np.abs(numeric - formula).max():.1e}")
# <</gradients>>
worst = max(np.abs(nm - fm).max() for nm, fm in checks.values())
assert worst < 1e-6

# derivatives with respect to a scalar: A(t) = S + t E
E = rng.normal(size=(p, p)); E = E + E.T
h = 1e-6
dinv = (np.linalg.inv(S + h * E) - np.linalg.inv(S - h * E)) / (2 * h)
assert np.allclose(dinv, -np.linalg.inv(S) @ E @ np.linalg.inv(S), atol=1e-6)
dlogdet = (np.linalg.slogdet(S + h * E)[1] - np.linalg.slogdet(S - h * E)[1]) / (2 * h)
assert np.isclose(dlogdet, np.trace(np.linalg.solve(S, E)), atol=1e-6)

# completing the square: minimizers of x'Ax - 2 b'x for nonnegative definite A
Bm = rng.normal(size=(3, p))
Ann = Bm.T @ Bm                                        # rank 3, nonnegative definite
bb = Ann @ rng.normal(size=p)                          # b in C(A): a minimum exists
G = np.linalg.pinv(Ann)
xstar = G @ bb
fmin = xstar @ Ann @ xstar - 2 * bb @ xstar
assert np.isclose(fmin, -bb @ G @ bb)
for _ in range(500):
    z = rng.normal(size=p) * 3
    assert z @ Ann @ z - 2 * bb @ z >= fmin - 1e-9

gen = Generated("ch01", "calculus_check", prefix="calc")
gen.num("worst", worst, 0, sci=True)
gen.write()
