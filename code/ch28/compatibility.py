"""Chapter 28, Section 2: compatibility constants of Gaussian designs with p > n.

phi^2(S) = min { s ||X v||^2 / n : ||v_S||_1 = 1, ||v_{S^c}||_1 <= 3 } for S = {1, ..., s}.
For each sign pattern of v_S (up to a global sign) the minimization is convex; it is solved by
accelerated projected gradient and certified by the Frank-Wolfe duality gap. Designs: p = 1000,
independent N(0, 1) entries with n = 200 (the design of Section 3) and n = 400, and an
equicorrelated design (correlation 0.5) with n = 400; columns scaled to ||x_j||^2 = n. For every
design the population compatibility constant is at least 1 - rho, and the coherence bound
1 - 16 s mu of the text is also computed.
"""
import numpy as np

from regbook import Generated

gen = Generated("ch28", "compatibility")


# <<compat>>
def proj_simplex(u):
    """Euclidean projection onto {u >= 0, sum(u) = 1}."""
    v = np.sort(u)[::-1]
    css = np.cumsum(v) - 1
    k = np.flatnonzero(v - css / np.arange(1, len(u) + 1) > 0)[-1]
    return np.maximum(u - css[k] / (k + 1), 0)


def proj_l1(w, radius):
    """Euclidean projection onto the l1 ball of the given radius."""
    if np.abs(w).sum() <= radius:
        return w
    return np.sign(w) * proj_simplex(np.abs(w) / radius) * radius


def compatibility(X, S, radius=3.0, iters=4000):
    """Return phi^2(S) and a certified lower bound (value minus duality gap)."""
    n, p = X.shape
    s = len(S)
    Sc = np.setdiff1d(np.arange(p), S)
    A, C = X[:, S], X[:, Sc]
    L = 2 * s / n * np.linalg.norm(X, 2) ** 2         # Lipschitz constant of the gradient
    best, lower = np.inf, np.inf
    for code in range(2 ** (s - 1)):                  # sign patterns of v_S up to a global sign
        sg = np.array([1.0] + [-1.0 if (code >> k) & 1 else 1.0 for k in range(s - 1)])
        u, w = np.full(s, 1 / s), np.zeros(len(Sc))
        yu, yw, t = u.copy(), w.copy(), 1.0
        for _ in range(iters):                        # accelerated projected gradient
            r = A @ (sg * yu) + C @ yw
            gu, gw = 2 * s / n * sg * (A.T @ r), 2 * s / n * (C.T @ r)
            u_new, w_new = proj_simplex(yu - gu / L), proj_l1(yw - gw / L, radius)
            t_new = (1 + np.sqrt(1 + 4 * t * t)) / 2
            yu = u_new + (t - 1) / t_new * (u_new - u)
            yw = w_new + (t - 1) / t_new * (w_new - w)
            u, w, t = u_new, w_new, t_new
        r = A @ (sg * u) + C @ w
        value = s / n * r @ r
        gu, gw = 2 * s / n * sg * (A.T @ r), 2 * s / n * (C.T @ r)
        gap = gu @ u - gu.min() + gw @ w + radius * np.abs(gw).max()
        best, lower = min(best, value), min(lower, value - gap)
    return best, lower
# <</compat>>


def coherence(X):
    G = X.T @ X / X.shape[0]
    return np.max(np.abs(G - np.diag(np.diag(G))))


# <<small>>
rng = np.random.default_rng(2820)
n, p = 100, 300
X = rng.normal(size=(n, p))
X /= np.sqrt(np.sum(X ** 2, axis=0) / n)            # ||x_j||^2 = n
print("smallest eigenvalue of X^T X / n:", np.linalg.eigvalsh(X.T @ X / n)[0])
for s in [1, 2, 3]:
    value, lower = compatibility(X, np.arange(s), iters=1500)
    print(f"s = {s}: phi^2 = {value:.3f}")
# <</small>>
assert abs(np.linalg.eigvalsh(X.T @ X / n)[0]) < 1e-10

designs = {}
rng = np.random.default_rng(2830)                   # the design of Section 3
Z = rng.normal(size=(200, 1000))
designs["indep200"] = (Z / np.sqrt(np.sum(Z ** 2, axis=0) / 200), 0.0)
rng = np.random.default_rng(2821)
Z = rng.normal(size=(400, 1000))
designs["indep400"] = (Z / np.sqrt(np.sum(Z ** 2, axis=0) / 400), 0.0)
rho = 0.5
Z = np.sqrt(1 - rho) * rng.normal(size=(400, 1000)) + np.sqrt(rho) * rng.normal(size=(400, 1))
designs["equi400"] = (Z / np.sqrt(np.sum(Z ** 2, axis=0) / 400), rho)

table = {}
for name, (X, r) in designs.items():
    mu = coherence(X)
    gen.num(f"{name}_mu", mu, 3)
    for s in [1, 2, 5]:
        value, lower = compatibility(X, np.arange(s))
        assert value - lower < 0.02 * value + 1e-4, (name, s, value, lower)
        table[name, s] = lower
        gen.num(f"{name}_s{s}", lower, 4 if lower < 0.01 else 3)
        gen.num(f"{name}_coh{s}", 1 - 16 * s * mu, 2)
        print(f"{name}  s = {s}: phi^2 >= {lower:.4f}  (value {value:.4f}); coherence bound {1 - 16 * s * mu:.2f}")
    # the sample constant is below the population lower bound 1 - rho once s is moderate
    assert table[name, 5] < 1 - r
    # and it falls as s grows
    assert table[name, 1] > table[name, 2] > table[name, 5]
    # the coherence bound is vacuous for all of these
    assert 1 - 16 * mu < 0
# doubling n raises the constant
assert all(table["indep400", s] > table["indep200", s] for s in [1, 2, 5])
gen.write()
