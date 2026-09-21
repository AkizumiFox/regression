"""Chapter 42, Section 1: what a higher polynomial degree buys and what it costs.

Monthly mean atmospheric CO2 (ppm) at Mauna Loa, March 1958 to December 2001
(public domain, statsmodels.datasets.co2; weekly values averaged by month), the series
of Example 5.3 (exm-lm-co2). The four harmonic columns are held fixed and only the
trend basis changes: powers of the time variable rescaled to [-1, 1].
"""
import matplotlib.pyplot as plt
import numpy as np
import statsmodels.api as sm

from regbook import COLORS, Generated, figure_path, use_book_style

# <<setup>>
import numpy as np
import statsmodels.api as sm

co2 = sm.datasets.co2.load_pandas().data["co2"]
monthly = co2.resample("MS").mean().dropna()         # weekly values -> monthly means
y = monthly.to_numpy()
t = (monthly.index.year + (monthly.index.month - 0.5) / 12).to_numpy()
n = len(y)


def harmonics(s):
    """The seasonal columns of Example 5.3, held fixed throughout the chapter."""
    return np.column_stack([np.cos(2 * np.pi * s), np.sin(2 * np.pi * s),
                            np.cos(4 * np.pi * s), np.sin(4 * np.pi * s)])


lo, hi = t.min(), t.max()
u = lambda s: 2 * (s - lo) / (hi - lo) - 1           # time rescaled to [-1, 1]

print(f"n = {n}, time from {lo:.2f} to {hi:.2f}")
# <</setup>>

# <<conditioning>>
for d in (2, 4, 8, 12, 16):
    raw = np.linalg.cond(np.vander(t, d + 1, increasing=True))
    centred = np.linalg.cond(np.vander(t - 1980.0, d + 1, increasing=True))
    scaled = np.linalg.cond(np.vander(u(t), d + 1, increasing=True))
    print(f"degree {d:2d}:  raw {raw:.2e}   centred {centred:.2e}   scaled {scaled:.2e}")
# <</conditioning>>

# <<degree>>
def trend_design(s, d):
    """Model matrix: powers 1, u, ..., u^d of rescaled time, then the four harmonics."""
    return np.column_stack([np.vander(u(s), d + 1, increasing=True), harmonics(s)])


for d in (2, 8, 16):
    X = trend_design(t, d)
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    resid = y - X @ beta
    p = X.shape[1]
    lev = np.einsum("ij,ij->i", X @ np.linalg.pinv(X.T @ X), X)   # leverages h_ii
    print(f"degree {d:2d}:  p = {p:2d}   sigma_hat = {np.sqrt(resid @ resid / (n - p)):.4f}"
          f"   mean h = {lev.mean():.4f}   max h = {lev.max():.4f}")
# <</degree>>

# <<influence>>
def trend_weight(d, i, s):
    """Change in the fitted trend at times s per unit change in the response at i."""
    X = trend_design(t, d)
    coef = np.linalg.pinv(X.T @ X) @ X[i]            # the row of (X'X)^{-1}X' for case i
    return np.vander(u(s), d + 1, increasing=True) @ coef[:d + 1]


grid = np.linspace(lo, hi, 600)
for d in (2, 12):
    w = 10 * trend_weight(d, n - 1, grid)            # move December 2001 by 10 ppm
    print(f"degree {d:2d}: the trend moves by {w[-1]:6.3f} ppm at the end, {w[0]:6.3f} ppm in 1958,"
          f"   and changes sign {np.sum(np.diff(np.sign(w)) != 0)} times")
# <</influence>>

# ---- assertions -------------------------------------------------------------
records = {}
for d in (2, 4, 8, 12, 16):
    records[d] = (np.linalg.cond(np.vander(t, d + 1, increasing=True)),
                  np.linalg.cond(np.vander(t - 1980.0, d + 1, increasing=True)),
                  np.linalg.cond(np.vander(u(t), d + 1, increasing=True)))
    # centring and rescaling do not change the column space, only the conditioning
    assert records[d][2] < records[d][1] < records[d][0]

fits = {}
for d in (2, 8, 16):
    X = trend_design(t, d)
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    resid = y - X @ beta
    p = X.shape[1]
    lev = np.einsum("ij,ij->i", X @ np.linalg.pinv(X.T @ X), X)
    assert np.isclose(lev.mean(), p / n)                       # mean leverage is p/n
    assert np.isclose(lev.sum(), np.trace(X @ np.linalg.pinv(X.T @ X) @ X.T))
    fits[d] = (p, beta, np.sqrt(resid @ resid / (n - p)), lev, resid)

# the residual standard deviation falls with the degree, the maximum leverage rises
assert fits[2][2] > fits[8][2] > fits[16][2]
assert fits[2][3].max() < fits[8][3].max() < fits[16][3].max()
# and the most influential month is the first one, at the left-hand edge of the range
for d in (2, 8, 16):
    assert fits[d][3].argmax() == 0
# the seasonal coefficients are almost unchanged by the trend basis
harm_shift = max(np.max(np.abs(fits[d][1][-4:] - fits[2][1][-4:])) for d in (8, 16))
assert harm_shift < 0.02

grid = np.linspace(lo, hi, 600)
se = {}
for d in (2, 8, 16):
    X = trend_design(t, d)
    Xg = trend_design(grid, d)
    G = np.linalg.pinv(X.T @ X)
    se[d] = fits[d][2] * np.sqrt(np.einsum("ij,jk,ik->i", Xg, G, Xg))
    # the standard error of the fitted mean is largest at an end of the range
    assert se[d].argmax() in (0, len(grid) - 1)
ratio = {d: se[d][0] / np.interp(1980.0, grid, se[d]) for d in (2, 8, 16)}
assert ratio[2] < ratio[8] < ratio[16]

w = {d: 10 * trend_weight(d, n - 1, grid) for d in (2, 12)}
signs = {d: int(np.sum(np.diff(np.sign(w[d])) != 0)) for d in (2, 12)}
for d in (2, 12):
    assert signs[d] <= d                                        # at most d sign changes
    assert np.min(np.abs(w[d])) > 0 or signs[d] > 0             # never identically zero
assert abs(w[12][0]) > abs(w[2][0])

gen = Generated("ch42", "polynomial_degree", prefix="deg")
gen.int("n", n)
for d, (raw, centred, scaled) in records.items():
    gen.num(f"kraw{d}", raw, 2, sci=True)
    gen.num(f"kcen{d}", centred, 2, sci=True)
    gen.num(f"ksca{d}", scaled, 2, sci=True)
for d in (2, 8, 16):
    gen.int(f"p{d}", fits[d][0])
    gen.num(f"sigma{d}", fits[d][2], 4)
    gen.num(f"hmax{d}", fits[d][3].max(), 4)
    gen.num(f"hbar{d}", fits[d][3].mean(), 4)
    gen.num(f"se0{d}", se[d][0], 4)
    gen.num(f"semid{d}", np.interp(1980.0, grid, se[d]), 4)
    gen.num(f"ratio{d}", ratio[d], 2)
gen.num("harmshift", harm_shift, 4)
for d in (2, 12):
    gen.num(f"wlast{d}", w[d][-1], 3)
    gen.num(f"wfirst{d}", w[d][0], 3)
    gen.int(f"signs{d}", signs[d])
gen.write()

# ---- figure -----------------------------------------------------------------
use_book_style()
fig, axes = plt.subplots(1, 2, figsize=(5.8, 2.5))
ax = axes[0]
ax.scatter(t, fits[2][4], s=2, color=COLORS["muted"], alpha=0.35, linewidths=0,
           label="residual, degree 2")
trend2 = np.vander(u(grid), 3, increasing=True) @ fits[2][1][:3]
for d, colour in ((8, COLORS["accent"]), (16, COLORS["second"])):
    trend = np.vander(u(grid), d + 1, increasing=True) @ fits[d][1][:d + 1]
    ax.plot(grid, trend - trend2, color=colour, label=f"degree {d}")
ax.axhline(0, color=COLORS["grid"], linewidth=0.6, zorder=0)
ax.set_xlabel("year")
ax.set_ylabel("ppm")
ax.set_title("(a) what the extra degrees fit")
ax.legend(frameon=False, loc="lower left", fontsize=6.5, ncol=1)
ax = axes[1]
for d, colour in ((2, COLORS["third"]), (8, COLORS["accent"]), (16, COLORS["second"])):
    ax.plot(grid, se[d], color=colour, label=f"degree {d}")
ax.set_xlabel("year")
ax.set_ylabel("standard error (ppm)")
ax.set_title("(b) standard error of the fitted mean")
ax.legend(frameon=False, loc="upper center", fontsize=6.5)
fig.tight_layout()
fig.savefig(figure_path("ch42", "polynomial_degree"))
