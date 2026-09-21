"""Chapter 42, Section 3: continuity constraints on a piecewise cubic, and the dimension count.

Annual flow of the Nile at Aswan, 1871-1970 (10^8 m^3; public domain,
statsmodels.datasets.nile). Two interior knots, 1899 and 1935. Fitting the three
intervals separately uses 12 parameters; imposing continuity of the function and of its
first and second derivatives removes 2, 4 and 6 of them, leaving the cubic spline space
of dimension d + 1 + K = 6.
"""
import matplotlib.pyplot as plt
import numpy as np
import statsmodels.api as sm
from scipy import stats

from regbook import COLORS, Generated, figure_path, use_book_style

# <<bases>>
import numpy as np
import statsmodels.api as sm

nile = sm.datasets.nile.load_pandas().data
year = nile["year"].to_numpy().astype(float)
y = nile["volume"].to_numpy()
n = len(y)
knots = np.array([1899.0, 1935.0])
d = 3
u = lambda s: 2 * (s - year.min()) / (year.max() - year.min()) - 1


def separate(s):
    """One free cubic per interval: (K+1)(d+1) columns, no constraints at the knots."""
    piece = np.digitize(s, knots)
    return np.column_stack([np.where(piece == k, u(s) ** j, 0.0)
                            for k in range(len(knots) + 1) for j in range(d + 1)])


def truncated(s, m):
    """Continuity of the function and of its first m derivatives at every knot."""
    cols = [u(s) ** j for j in range(d + 1)]
    for k in knots:
        cols += [np.clip(u(s) - u(k), 0.0, None) ** r for r in range(m + 1, d + 1)]
    return np.column_stack(cols)


def fit(X):
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    resid = y - X @ beta
    return np.linalg.matrix_rank(X), resid @ resid, X @ beta


for name, X in [("separate cubics", separate(year)), ("continuous", truncated(year, 0)),
                ("continuous derivative", truncated(year, 1)), ("cubic spline", truncated(year, 2))]:
    p, rss, _ = fit(X)
    print(f"{name:24s} p = {p:2d}   residual sum of squares {rss:10.0f}")
# <</bases>>

# <<ftests>>
p_full, rss_full, _ = fit(separate(year))
for m, label in ((0, "continuity"), (1, "first derivative"), (2, "second derivative")):
    p_m, rss_m, _ = fit(truncated(year, m))
    q = p_full - p_m
    F = ((rss_m - rss_full) / q) / (rss_full / (n - p_full))
    print(f"up to the {label:17s}: q = {q}, F = {F:5.2f}, p-value {1 - stats.f.cdf(F, q, n - p_full):.4f}")
# <</ftests>>

# <<step>>
step = np.column_stack([np.ones(n), (year >= 1899.0).astype(float)])
beta, *_ = np.linalg.lstsq(step, y, rcond=None)
resid = y - step @ beta
se = np.sqrt(resid @ resid / (n - 2) * np.linalg.inv(step.T @ step)[1, 1])
print(f"mean before 1899 {beta[0]:.1f}, after {beta[0] + beta[1]:.1f}, "
      f"difference {beta[1]:.1f} (se {se:.1f}, t = {beta[1] / se:.2f})")
print(f"two parameters, residual sum of squares {resid @ resid:.0f}")
# <</step>>

# ---- assertions -------------------------------------------------------------
results = {}
for m in (-1, 0, 1, 2):
    X = separate(year) if m < 0 else truncated(year, m)
    results[m] = fit(X)
    # the dimension count: (K+1)(d+1) constraints removed one at a time
    assert results[m][0] == (len(knots) + 1) * (d + 1) - (m + 1) * len(knots)
assert results[2][0] == d + 1 + len(knots)
# nested spaces: more constraints, larger residual sum of squares
assert results[-1][1] < results[0][1] < results[1][1] < results[2][1]
# the truncated power columns span a subspace of the unconstrained space
joint = np.column_stack([separate(year), truncated(year, 2)])
assert np.linalg.matrix_rank(joint) == results[-1][0]
# the fitted curves have exactly the continuity claimed: recover each interval's cubic
# and compare its derivatives with the next one's at the knot
def piece_polynomials(m):
    beta, *_ = np.linalg.lstsq(separate(year) if m < 0 else truncated(year, m), y, rcond=None)
    design = separate if m < 0 else (lambda s: truncated(s, m))
    edges = np.r_[year.min(), knots, year.max()]
    out = []
    for lo_k, hi_k in zip(edges[:-1], edges[1:]):
        s = np.linspace(lo_k + 1e-6, hi_k - 1e-6, 20)
        out.append(np.polyfit(u(s), design(s) @ beta, d))
    return out


for m in (-1, 0, 1, 2):
    pieces = piece_polynomials(m)
    for k, v in enumerate(u(knots)):
        jumps = [abs(np.polyval(np.polyder(pieces[k + 1], r), v)
                     - np.polyval(np.polyder(pieces[k], r), v)) for r in range(d + 1)]
        scale = max(abs(np.polyval(np.polyder(pieces[k], r), v)) for r in range(d + 1)) + 1.0
        for r in range(d + 1):
            if r <= m:
                assert jumps[r] < 1e-6 * scale          # continuous up to order m
            elif r == m + 1:
                assert jumps[r] > 1e-3 * scale          # and no further

F_stats, p_values = {}, {}
for m in (0, 1, 2):
    q = results[-1][0] - results[m][0]
    F_stats[m] = ((results[m][1] - results[-1][1]) / q) / (results[-1][1] / (n - results[-1][0]))
    p_values[m] = 1 - stats.f.cdf(F_stats[m], q, n - results[-1][0])
assert p_values[0] < 0.01                     # continuity at the knots is rejected

step_rss = resid @ resid
assert step_rss < results[2][1]               # two parameters beat the six-parameter spline
cubic_poly = np.vander(u(year), 4, increasing=True)
_, rss_cubic, _ = fit(cubic_poly)
assert step_rss < rss_cubic

gen = Generated("ch42", "piecewise_nile", prefix="pw")
gen.int("n", n)
names = {-1: "sep", 0: "c0", 1: "c1", 2: "c2"}
for m, key in names.items():
    gen.int(f"p{key}", results[m][0])
    gen.num(f"rss{key}", results[m][1], 0)
for m in (0, 1, 2):
    gen.num(f"F{m}", F_stats[m], 2)
    gen.num(f"pval{m}", p_values[m], 4)
gen.num("mean1", beta[0], 1)
gen.num("mean2", beta[0] + beta[1], 1)
gen.num("diff", beta[1], 1)
gen.num("sediff", se, 1)
gen.num("tdiff", beta[1] / se, 2)
gen.num("rssstep", step_rss, 0)
gen.num("rsscubic", rss_cubic, 0)
gen.write()

# ---- figure -----------------------------------------------------------------
use_book_style()
fig, axes = plt.subplots(2, 2, figsize=(5.8, 4.0), sharex=True, sharey=True)
titles = {-1: "(a) separate cubics", 0: "(b) continuous",
          1: "(c) continuous first derivative", 2: "(d) cubic spline"}
for ax, m in zip(axes.ravel(), (-1, 0, 1, 2)):
    ax.scatter(year, y, s=5, color=COLORS["muted"], alpha=0.55, linewidths=0)
    for k in knots:
        ax.axvline(k, color=COLORS["grid"], linewidth=0.6, zorder=0)
    if m < 0:
        for lo_k, hi_k in zip(np.r_[year.min(), knots], np.r_[knots, year.max()]):
            seg = np.linspace(lo_k, hi_k, 200)
            inside = (seg >= lo_k) & (seg <= hi_k)
            ax.plot(seg[inside], separate(np.clip(seg[inside], lo_k, hi_k - 1e-9)) @
                    np.linalg.lstsq(separate(year), y, rcond=None)[0],
                    color=COLORS["accent"])
    else:
        g = np.linspace(year.min(), year.max(), 600)
        ax.plot(g, truncated(g, m) @ np.linalg.lstsq(truncated(year, m), y, rcond=None)[0],
                color=COLORS["accent"])
    ax.set_title(f"{titles[m]}: {results[m][0]} parameters")
    ax.set_ylim(400, 1500)
for ax in axes[1]:
    ax.set_xlabel("year")
for ax in axes[:, 0]:
    ax.set_ylabel("flow ($10^8$ m$^3$)")
fig.tight_layout()
fig.savefig(figure_path("ch42", "piecewise_nile"))
