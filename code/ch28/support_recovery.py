"""Chapter 28, Section 4: when the lasso selects the right variables.

Two Gaussian designs with p = 500. The response depends on x_1 and x_2 (coefficients 0.5 and 0.5);
x_3 has correlation r with each of them; x_4, ..., x_p are independent of everything. The
irrepresentable quantity for x_3 is 2r in the population: r = 0.3 satisfies the condition
(0.6 < 1), r = 0.6 violates it (1.2 > 1). For each n one design is drawn and held fixed, and
normal errors (sigma = 1) are drawn 200 times; the lasso path is computed on a grid of 60
penalties from 1.5 down to 0.1, and we record whether some penalty on the grid gives exactly the right signs.
"""
import matplotlib.pyplot as plt
import numpy as np

from regbook import COLORS, Generated, figure_path, use_book_style

gen = Generated("ch28", "support_recovery")

# <<lasso>>
def soft(z, t):
    """Soft thresholding S(z, t) = sign(z) max(|z| - t, 0)."""
    return np.sign(z) * np.maximum(np.abs(z) - t, 0.0)


def lasso_cd(X, Y, lam, B=None, tol=1e-9, max_sweeps=2000):
    """Minimize ||y - X b||^2 / (2n) + lam ||b||_1 for every column y of Y, by coordinate descent."""
    n, p = X.shape
    Y = Y.reshape(n, -1)
    B = np.zeros((p, Y.shape[1])) if B is None else B.copy()
    R = Y - X @ B                                    # residuals, one column per response
    scale = np.sum(X ** 2, axis=0) / n
    active = np.arange(p)
    for sweep in range(max_sweeps):
        change = 0.0
        for j in active:
            z = X[:, j] @ R / n + scale[j] * B[j]
            new = soft(z, lam) / scale[j]
            d = new - B[j]
            if np.any(d != 0):
                R -= np.outer(X[:, j], d)
                B[j] = new
                change = max(change, np.max(np.abs(d)))
        if change < tol:
            if len(active) == p:                     # converged over all coordinates
                break
            active = np.arange(p)                    # check every coordinate once more
        else:
            active = np.flatnonzero(np.any(B != 0, axis=1)) if sweep % 5 else np.arange(p)
    return B


def kkt_violation(X, Y, B, lam):
    """Largest violation of the lasso optimality conditions (zero at a solution)."""
    n = X.shape[0]
    G = X.T @ (Y.reshape(n, -1) - X @ B) / n         # gradient part, should lie in lam * subdifferential
    off = np.where(B == 0, np.maximum(np.abs(G) - lam, 0), np.abs(G - lam * np.sign(B)))
    return off.max()
# <</lasso>>


# <<design>>
def make_design(rng, n, p, r):
    """x_1, x_2 independent; x_3 correlated r with each; x_4, ..., x_p independent noise."""
    Z = rng.normal(size=(n, p))
    Z[:, 2] = r * Z[:, 0] + r * Z[:, 1] + np.sqrt(1 - 2 * r ** 2) * Z[:, 2]
    return Z / np.sqrt(np.sum(Z ** 2, axis=0) / n)  # ||x_j||^2 = n


def irrepresentable(X, S, signs):
    """max over j outside S of |x_j^T X_S (X_S^T X_S)^{-1} sign(beta_S)|."""
    A = X[:, S]
    w = A @ np.linalg.solve(A.T @ A, signs)
    out = np.delete(X.T @ w, S)
    return np.max(np.abs(out))


def recovers_signs(X, Y, beta, lams):
    """For each column of Y: does some penalty in lams give exactly sign(beta)?"""
    ok = np.zeros(Y.shape[1], dtype=bool)
    B = None
    for lam in lams:                                 # from large to small, warm starts
        B = lasso_cd(X, Y, lam, B)
        ok |= np.all(np.sign(B) == np.sign(beta)[:, None], axis=0)
    return ok
# <</design>>


# <<small>>
rng = np.random.default_rng(2840)
n, p, sigma = 100, 200, 1.0
beta = np.zeros(p)
beta[:2] = 0.5
S = np.array([0, 1])
lams = np.geomspace(1.0, 0.1, 30)
for r in [0.3, 0.6]:
    X = make_design(rng, n, p, r)
    Y = (X @ beta)[:, None] + sigma * rng.normal(size=(n, 100))
    print(f"r = {r}: irrepresentable quantity {irrepresentable(X, S, np.ones(2)):.2f},"
          f" sign recovery in {recovers_signs(X, Y, beta, lams).mean():.2f} of 100 data sets")
# <</small>>

p, sigma, reps = 500, 1.0, 200
beta = np.zeros(p)
beta[:2] = 0.5
S = np.array([0, 1])
lams = np.geomspace(1.5, 0.1, 60)
ns = [50, 100, 200, 400, 800]
prob = {}
irr = {}
for r in [0.3, 0.6]:
    rng = np.random.default_rng(2841 if r == 0.3 else 2842)
    prob[r], irr[r] = [], []
    for n in ns:
        X = make_design(rng, n, p, r)
        Y = (X @ beta)[:, None] + sigma * rng.normal(size=(n, reps))
        ok = recovers_signs(X, Y, beta, lams)
        prob[r].append(ok.mean())
        irr[r].append(irrepresentable(X, S, np.ones(2)))
        print(f"r = {r}, n = {n}: irrepresentable {irr[r][-1]:.3f}, P(sign recovery) = {ok.mean():.3f}")
    prob[r], irr[r] = np.array(prob[r]), np.array(irr[r])

se = np.sqrt(0.25 / reps)
# condition violated in the sample: recovery probability at most 1/2 (proposition of the text)
for k, n in enumerate(ns):
    if irr[0.6][k] > 1:
        assert prob[0.6][k] <= 0.5 + 3 * se
# condition holds (with margin) for the larger samples: recovery probability tends to one
big = np.array(ns) >= 200
assert np.all(irr[0.3][big] < 1)
assert prob[0.3][-1] > 0.95 and prob[0.3][0] < 0.6
# violated: the probability rises slowly with n, largest at n = 800, and stays far below 1/2
assert prob[0.6][-1] == prob[0.6].max() and prob[0.6][-1] > prob[0.6][ns.index(100)]
assert prob[0.6].max() < 0.2
gen.int("p", p)
gen.int("reps", reps)
gen.num("irr_ok_800", irr[0.3][-1], 2)
gen.num("irr_bad_800", irr[0.6][-1], 2)
gen.num("prob_ok_50", prob[0.3][0], 2)
gen.num("prob_ok_100", prob[0.3][ns.index(100)], 2)
gen.num("prob_ok_800", prob[0.3][-1], 3)
gen.num("prob_bad_100", prob[0.6][ns.index(100)], 2)
gen.num("prob_bad_800", prob[0.6][-1], 3)

# a single path in the violated design with n = 800: x_3 enters first and stays
rng = np.random.default_rng(2843)
n = 800
X = make_design(rng, n, p, 0.6)
y = X @ beta + sigma * rng.normal(size=n)
grid = np.geomspace(1.6, 0.02, 100)
path, B = [], None
for lam in grid:
    B = lasso_cd(X, y, lam, B)
    path.append(B[:, 0].copy())
path = np.array(path)
first = [np.flatnonzero(path[:, j] != 0)[0] if np.any(path[:, j] != 0) else len(grid) for j in range(3)]
assert first[2] < min(first[0], first[1])            # x_3 enters before x_1 and x_2
nz3 = path[first[2]:, 2] != 0
assert np.all(nz3)                                    # and never leaves
gen.num("inner3", X[:, 2] @ (X @ beta) / n, 2)       # x_3^T X beta / n
gen.write()

# ---- figure ---------------------------------------------------------------------------------
use_book_style()
fig, axes = plt.subplots(1, 2, figsize=(5.8, 2.5))
ax = axes[0]
ax.plot(ns, prob[0.3], "o-", ms=3, color=COLORS["accent"], label="$r=0.3$ (condition holds)")
ax.plot(ns, prob[0.6], "s-", ms=3, color=COLORS["second"], label="$r=0.6$ (condition fails)")
ax.axhline(0.5, color=COLORS["grid"], lw=0.8, zorder=0)
ax.set_xscale("log")
ax.set_xticks(ns)
ax.set_xticklabels([str(k) for k in ns])
ax.xaxis.set_minor_formatter(plt.NullFormatter())
ax.set_ylim(-0.03, 1.03)
ax.set_xlabel("sample size $n$ ($p=500$)")
ax.set_ylabel("P(some $\\lambda$ recovers the signs)")
ax.set_title("(a) sign recovery")
ax.text(800, 0.92, "$r=0.3$ (condition holds)", color=COLORS["accent"], ha="right", va="top", fontsize=8)
ax.text(800, 0.24, "$r=0.6$ (condition fails)", color=COLORS["second"], ha="right", va="bottom", fontsize=8)
ax = axes[1]
for j in range(3, 60):
    ax.plot(grid, path[:, j], color=COLORS["grid"], lw=0.6)
ax.plot(grid, path[:, 0], color=COLORS["accent"], label="$x_1$")
ax.plot(grid, path[:, 1], color=COLORS["accent"], ls="--", label="$x_2$")
ax.plot(grid, path[:, 2], color=COLORS["second"], label="$x_3$ (inactive)")
ax.set_xscale("log")
ax.invert_xaxis()
ax.set_xlabel("$\\lambda$")
ax.set_ylabel("coefficient")
ax.set_title("(b) one path, $r=0.6$, $n=800$")
ax.legend(frameon=False, loc="upper left")
fig.tight_layout()
fig.savefig(figure_path("ch28", "support_recovery"))
