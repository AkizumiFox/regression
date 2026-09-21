"""Chapter 43, Section 6: the mixed-model form of a penalized spline, REML, and the bands.

Checks the explicit transformation that turns the difference penalty into a normal prior,
estimates the smoothing parameter by REML, and measures the coverage of the Wahba-Nychka
bands, pointwise and averaged across the function.
"""
import matplotlib.pyplot as plt
import numpy as np

from regbook import COLORS, Generated, figure_path, use_book_style


def f_true(x):
    """The mean function of Section 43.1."""
    return 2 * x + np.exp(-25 * (x - 0.35) ** 2) - 0.75 * np.exp(-50 * (x - 0.80) ** 2)


def uniform_knots(a, b, K, d):
    """K interior knots equally spaced in (a, b), continued d further steps at each end."""
    delta = (b - a) / (K + 1)
    return a + delta * np.arange(-d, K + 2 + d)


def bspline_basis(x, kn, d):
    """The B-splines of degree d on the knot vector kn, by the de Boor recurrence."""
    x = np.atleast_1d(np.asarray(x, float))
    B = np.array([(x >= kn[j]) & (x < kn[j + 1]) for j in range(len(kn) - 1)], float).T
    last = np.max(np.nonzero(kn[:-1] < kn[1:])[0])
    B[x == kn[-1], last] = 1.0
    for deg in range(1, d + 1):
        C = np.zeros((len(x), B.shape[1] - 1))
        for j in range(C.shape[1]):
            if kn[j + deg] > kn[j]:
                C[:, j] += (x - kn[j]) / (kn[j + deg] - kn[j]) * B[:, j]
            if kn[j + deg + 1] > kn[j + 1]:
                C[:, j] += (kn[j + deg + 1] - x) / (kn[j + deg + 1] - kn[j + 1]) * B[:, j + 1]
        B = C
    return B


# <<transform>>
def split_penalty(P, tol=1e-9):
    """Bases X_g of N(P) and Z_g of its complement with Z_g' P Z_g = I and X_g' P = 0."""
    e, U = np.linalg.eigh(P)
    pos = e > tol * e.max()
    L = U[:, pos] * np.sqrt(e[pos])                    # P = L L' with L of full column rank
    Xg = U[:, ~pos]                                    # the null space of the penalty
    Zg = L @ np.linalg.inv(L.T @ L)
    return Xg, Zg
# <</transform>>


n, sigma, K, d, k = 120, 0.25, 20, 3, 2
x = (np.arange(1, n + 1) - 0.5) / n
kn = uniform_knots(0.0, 1.0, K, d)
B = bspline_basis(x, kn, d)
m = B.shape[1]
Dk = np.diff(np.eye(m), n=k, axis=0)
P = Dk.T @ Dk
Xg, Zg = split_penalty(P)
X, Z = B @ Xg, B @ Zg
r = Z.shape[1]

# --- the transformation does what Theorem 43.6 says ---------------------------
assert Xg.shape[1] == k and r == m - k
assert np.allclose(P @ Xg, 0.0, atol=1e-8)                       # X_g spans the null space
assert np.allclose(Zg.T @ P @ Zg, np.eye(r), atol=1e-8)          # and the penalty becomes u'u
assert abs(np.linalg.det(np.column_stack([Xg, Zg]))) > 1e-8      # a change of basis
rng = np.random.default_rng(431)
g = rng.normal(size=m)                                           # any coefficient vector
beta_u = np.linalg.solve(np.column_stack([Xg, Zg]), g)
print("the penalty, in the old and the new coordinates:", g @ P @ g, beta_u[k:] @ beta_u[k:])
assert np.isclose(g @ P @ g, beta_u[k:] @ beta_u[k:])            # the penalty is ||u||^2

# --- the penalized fit is the mixed-model fit ---------------------------------
y = f_true(x) + sigma * rng.normal(size=n)
lam = 1.0
fit_pen = B @ np.linalg.solve(B.T @ B + lam * P, B.T @ y)
C = np.block([[X.T @ X, X.T @ Z], [Z.T @ X, Z.T @ Z + lam * np.eye(r)]])   # Henderson's equations
sol = np.linalg.solve(C, np.concatenate([X.T @ y, Z.T @ y]))
fit_mix = X @ sol[:k] + Z @ sol[k:]
print("penalized fit versus mixed-model fit:", np.max(np.abs(fit_pen - fit_mix)))
assert np.max(np.abs(fit_pen - fit_mix)) < 1e-8

# --- REML for the smoothing parameter ----------------------------------------
# <<reml>>
Uz, sz, _ = np.linalg.svd(Z, full_matrices=False)     # the design is fixed, so do this once


def reml(y, lam):
    """Profiled REML criterion (to be minimized) for V = I + Z Z' / lam."""
    w = sz**2 / (lam + sz**2)
    def Vi(M):                                         # multiplication by V^{-1}
        M2 = M if M.ndim == 2 else M[:, None]
        out = M2 - Uz @ (w[:, None] * (Uz.T @ M2))
        return out if M.ndim == 2 else out[:, 0]

    logdetV = np.sum(np.log1p(sz**2 / lam))
    ViX = Vi(X)
    A = X.T @ ViX
    resid = Vi(y) - ViX @ np.linalg.solve(A, ViX.T @ y)
    s2 = y @ resid / (n - k)
    return 0.5 * ((n - k) * np.log(s2) + logdetV + np.linalg.slogdet(A)[1]), s2


def reml_lambda(y, grid):
    return grid[int(np.argmin([reml(y, l)[0] for l in grid]))]
# <</reml>>


grid = np.geomspace(1e-6, 1e6, 121)
lam_reml = reml_lambda(y, grid)
S_reml = B @ np.linalg.solve(B.T @ B + lam_reml * P, B.T)
df_reml = np.trace(S_reml)
s2_reml = reml(y, lam_reml)[1]
assert 3 < df_reml < 25

# GCV on the same data set, for comparison
gcvs = []
for l in grid:
    S = B @ np.linalg.solve(B.T @ B + l * P, B.T)
    gcvs.append(np.mean((y - S @ y) ** 2) / (1 - np.trace(S) / n) ** 2)
lam_gcv = grid[int(np.argmin(gcvs))]
df_gcv = np.trace(B @ np.linalg.solve(B.T @ B + lam_gcv * P, B.T))
print(f"REML: lambda = {lam_reml:.3g}, df = {df_reml:.2f}; GCV: df = {df_gcv:.2f}")

# --- coverage of the Bayesian bands ------------------------------------------
reps = 400
z975 = 1.959963984540054
covered = np.zeros(n)
covered_naive = np.zeros(n)
widths, dfs = [], []
for _ in range(reps):
    yy = f_true(x) + sigma * rng.normal(size=n)
    l = reml_lambda(yy, grid)
    A = np.linalg.solve(B.T @ B + l * P, B.T)
    S = B @ A
    fit = S @ yy
    dfs.append(np.trace(S))
    s2 = np.sum((yy - fit) ** 2) / (n - np.trace(S))
    half = z975 * np.sqrt(s2 * np.diag(S))            # the posterior standard deviation
    covered += np.abs(fit - f_true(x)) <= half
    half_naive = z975 * np.sqrt(s2 * np.diag(S @ S.T))   # the variance-only interval
    covered_naive += np.abs(fit - f_true(x)) <= half_naive
    widths.append(np.mean(half))
cover = covered / reps
cover_naive = covered_naive / reps
mean_cover = cover.mean()
mean_cover_naive = cover_naive.mean()
assert mean_cover > mean_cover_naive
assert 0.90 < mean_cover < 0.99
assert cover.min() < 0.93 and cover.max() > 0.97        # pointwise coverage is not uniform

# the same bands read as frequentist intervals for E(f-hat), ignoring the bias
S = B @ np.linalg.solve(B.T @ B + lam_reml * P, B.T)
freq_sd = np.sqrt(s2_reml * np.diag(S @ S.T))
bayes_sd = np.sqrt(s2_reml * np.diag(S))
assert np.all(bayes_sd > freq_sd)
ratio = float(np.mean(bayes_sd / freq_sd))

gen = Generated("ch43", "mixed_bayes")
gen.int("n", n)
gen.int("K", K)
gen.int("m", m)
gen.int("r", r)
gen.int("reps", reps)
gen.num("sigma", sigma, 2)
gen.num("lam_reml", lam_reml, 2, sci=True)
gen.num("df_reml", df_reml, 2)
gen.num("sigma_reml", np.sqrt(s2_reml), 3)
gen.num("df_gcv", df_gcv, 2)
gen.num("mean_cover", mean_cover, 3)
gen.num("mean_cover_naive", mean_cover_naive, 3)
gen.num("min_cover", cover.min(), 3)
gen.num("max_cover", cover.max(), 3)
gen.num("ratio", ratio, 2)
gen.num("mean_half", float(np.mean(widths)), 3)
gen.write()

# ---- figure -----------------------------------------------------------------
use_book_style()
fig, axes = plt.subplots(1, 2, figsize=(5.8, 2.4))
ax = axes[0]
fit = S @ y
half = z975 * bayes_sd
ax.scatter(x, y, s=6, color=COLORS["muted"], alpha=0.7, linewidths=0)
ax.fill_between(x, fit - half, fit + half, color=COLORS["accent"], alpha=0.18, linewidth=0)
ax.plot(x, f_true(x), color=COLORS["ink"], linewidth=1.0, label="true $f$")
ax.plot(x, fit, color=COLORS["accent"], linewidth=1.1, label=f"REML fit, df $={df_reml:.1f}$")
ax.set_xlabel("$x$")
ax.set_ylabel("$y$")
ax.set_ylim(-0.9, 3.1)
ax.set_title("(a) fit and Bayesian band")
ax.legend(frameon=False, loc="upper left", handlelength=1.4)
ax = axes[1]
ax.axhline(0.95, color=COLORS["ink"], linewidth=0.8, linestyle="--")
ax.plot(x, cover, color=COLORS["accent"], linewidth=1.0)
ax.axhline(mean_cover, color=COLORS["second"], linewidth=0.9, linestyle=":")
ax.text(0.02, 0.873, f"average {mean_cover:.3f}", color=COLORS["second"], fontsize=7)
ax.set_xlabel("$x$")
ax.set_ylabel("coverage")
ax.set_ylim(0.86, 1.005)
ax.set_title("(b) pointwise coverage, nominal 0.95")
fig.tight_layout()
fig.savefig(figure_path("ch43", "bayes_bands"))
