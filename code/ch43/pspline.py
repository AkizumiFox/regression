"""Chapter 43, Section 3: penalized splines (P-splines) on the annual Nile flow.

Data: statsmodels.datasets.nile, annual flow of the Nile at Aswan 1871-1970, public domain.
Shows the difference penalty, the effect of the smoothing parameter, the polynomial limit
as lambda grows, and how little the number of knots matters once the penalty is there.
"""
import matplotlib.pyplot as plt
import numpy as np
import statsmodels.api as sm

from regbook import COLORS, Generated, figure_path, use_book_style


def uniform_knots(a, b, K, d):
    """K interior knots equally spaced in (a, b), continued d further steps at each end."""
    delta = (b - a) / (K + 1)
    return a + delta * np.arange(-d, K + 2 + d)


def knot_vector(a, b, K, d):
    """The clamped sequence: d+1 copies of a, K equally spaced interior knots, d+1 copies of b."""
    interior = a + (b - a) * np.arange(1, K + 1) / (K + 1)
    return np.concatenate([np.full(d + 1, a), interior, np.full(d + 1, b)])


def bspline_basis(x, kn, d):
    """The B-splines of degree d on the knot vector kn, by the de Boor recurrence."""
    x = np.atleast_1d(np.asarray(x, float))
    B = np.array([(x >= kn[j]) & (x < kn[j + 1]) for j in range(len(kn) - 1)], float).T
    last = np.max(np.nonzero(kn[:-1] < kn[1:])[0])       # close the rightmost interval
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


# <<pspline>>
def difference_matrix(m, k):
    """The k-th difference operator on m coefficients, an (m-k) x m matrix."""
    return np.diff(np.eye(m), n=k, axis=0)


def pspline(x, y, lam, K=20, d=3, k=2, a=None, b=None):
    """Penalized spline fit: minimize ||y - B gamma||^2 + lam ||Delta^k gamma||^2."""
    a = x.min() if a is None else a
    b = x.max() if b is None else b
    kn = uniform_knots(a, b, K, d)
    B = bspline_basis(x, kn, d)
    Dk = difference_matrix(B.shape[1], k)
    A = B.T @ B + lam * Dk.T @ Dk
    gamma = np.linalg.solve(A, B.T @ y)
    df = np.trace(np.linalg.solve(A, B.T @ B))          # trace of the smoother matrix
    return gamma, df, kn, d


def evaluate(gamma, kn, d, grid):
    return bspline_basis(grid, kn, d) @ gamma
# <</pspline>>


def lambda_for_df(x, y, target, **kw):
    """Bisect on log(lambda) until the effective degrees of freedom hit the target."""
    lo, hi = -8.0, 16.0
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        if pspline(x, y, 10.0**mid, **kw)[1] > target:
            lo = mid
        else:
            hi = mid
    return 10.0 ** (0.5 * (lo + hi))


data = sm.datasets.nile.load_pandas().data
year = data["year"].to_numpy()
y = data["volume"].to_numpy()
n = len(y)
x = (year - year.min()) / (year.max() - year.min())      # the covariate, rescaled to [0, 1]
grid = np.linspace(0, 1, 401)
grid_year = year.min() + grid * (year.max() - year.min())

# --- three amounts of smoothing ----------------------------------------------
curves = {}
for target in (3.0, 8.0, 20.0):
    lam = lambda_for_df(x, y, target)
    gamma, df, kn, d = pspline(x, y, lam)
    curves[target] = (lam, df, evaluate(gamma, kn, d, grid))
    assert abs(df - target) < 1e-6

# --- the limit lambda -> infinity is the least squares line -------------------
# <<limit>>
gamma_big, df_big, kn, d = pspline(x, y, 1e12)
line = np.polyval(np.polyfit(x, y, 1), grid)
print("df at lambda = 1e12:", round(df_big, 6))
print("largest gap from the least squares line:", np.max(np.abs(evaluate(gamma_big, kn, d, grid) - line)))
# <</limit>>
assert abs(df_big - 2.0) < 1e-3
gap_line = np.max(np.abs(evaluate(gamma_big, kn, d, grid) - line))
assert gap_line < 1e-3 * (y.max() - y.min())

# the same limit on the CLAMPED sequence, whose Greville abscissae are not equally spaced:
# the limiting fit is then a spline close to, but not equal to, the least squares line
kn_c = knot_vector(0.0, 1.0, 20, 3)
Bc = bspline_basis(x, kn_c, 3)
Dc = difference_matrix(Bc.shape[1], 2)
gamma_c = np.linalg.solve(Bc.T @ Bc + 1e12 * Dc.T @ Dc, Bc.T @ y)
gap_line_clamped = np.max(np.abs(bspline_basis(grid, kn_c, 3) @ gamma_c - line))
assert gap_line_clamped > 100 * gap_line

gamma_flat, df_flat, kn1, d1 = pspline(x, y, 1e12, k=1)
assert abs(df_flat - 1.0) < 1e-3
assert np.max(np.abs(evaluate(gamma_flat, kn1, d1, grid) - y.mean())) < 1e-3 * (y.max() - y.min())

# --- the number of knots hardly matters once the penalty is there -------------
knot_fits = {}
for K in (10, 20, 40, 80):
    lam = lambda_for_df(x, y, 8.0, K=K)
    gamma, df, kn, d = pspline(x, y, lam, K=K)
    knot_fits[K] = (lam, evaluate(gamma, kn, d, grid))
    assert abs(df - 8.0) < 1e-6
gaps = {K: np.max(np.abs(knot_fits[K][1] - knot_fits[80][1])) for K in (10, 20, 40)}
spread = np.max(np.abs(np.ptp(np.array([knot_fits[K][1] for K in (10, 20, 40, 80)]), axis=0)))
assert spread < 0.03 * (y.max() - y.min())

# the unpenalized fits, by contrast, are all over the place
raw = {}
for K in (10, 20, 40):
    gamma, df, kn, d = pspline(x, y, 0.0, K=K)
    raw[K] = evaluate(gamma, kn, d, grid)
spread_raw = np.max(np.abs(np.ptp(np.array([raw[K] for K in (10, 20, 40)]), axis=0)))
assert spread_raw > 8 * spread

gen = Generated("ch43", "pspline")
gen.int("n", n)
gen.int("y0", int(year.min()))
gen.int("y1", int(year.max()))
gen.num("lam3", curves[3.0][0], 1, sci=True)
gen.num("lam8", curves[8.0][0], 1, sci=True)
gen.num("lam20", curves[20.0][0], 1, sci=True)
gen.num("spread", spread, 2)
gen.num("spread_raw", spread_raw, 1)
gen.num("gap10", gaps[10], 2)
gen.num("gap40", gaps[40], 3)
gen.num("range", y.max() - y.min(), 0)
gen.num("spread_pct", 100 * spread / (y.max() - y.min()), 2)
gen.num("gap_line", gap_line, 1, sci=True)
gen.num("gap_clamped", gap_line_clamped, 1)
gen.num("gap_clamped_pct", 100 * gap_line_clamped / (y.max() - y.min()), 1)
gen.write()

# ---- figure -----------------------------------------------------------------
use_book_style()
fig, axes = plt.subplots(1, 2, figsize=(5.8, 2.4))
ax = axes[0]
ax.scatter(year, y, s=6, color=COLORS["muted"], alpha=0.8, linewidths=0)
for target, col, ls in ((3.0, COLORS["third"], "--"), (8.0, COLORS["accent"], "-"), (20.0, COLORS["second"], "-")):
    ax.plot(grid_year, curves[target][2], color=col, linestyle=ls, linewidth=1.1,
            label=f"df $= {target:g}$")
ax.set_xlabel("year")
ax.set_ylabel("flow")
ax.set_title("(a) three amounts of smoothing")
ax.legend(frameon=False, loc="upper right", handlelength=1.4)
ax = axes[1]
ax.scatter(year, y, s=6, color=COLORS["muted"], alpha=0.5, linewidths=0)
styles = {10: ("-", COLORS["accent"]), 20: ("--", COLORS["second"]), 40: (":", COLORS["third"]),
          80: ("-.", COLORS["thread"])}
for K, (ls, col) in styles.items():
    ax.plot(grid_year, knot_fits[K][1], linestyle=ls, color=col, linewidth=1.0, label=f"$K={K}$")
ax.set_xlabel("year")
ax.set_title("(b) four knot counts, df $=8$")
ax.legend(frameon=False, loc="upper right", ncol=2, handlelength=1.6, columnspacing=1.0)
fig.tight_layout()
fig.savefig(figure_path("ch43", "pspline_nile"))
