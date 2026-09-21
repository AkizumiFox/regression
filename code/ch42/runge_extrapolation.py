"""Chapter 42, Section 1: two failures of a high-degree polynomial.

(a) Runge's phenomenon: interpolating 1/(1+25x^2) at equally spaced nodes on [-1,1]
    diverges as the number of nodes grows, while the same degree at Chebyshev nodes
    converges. (b) Extrapolation: polynomial trends fitted to the Mauna Loa CO2 series
    up to the end of 1990 and used to predict 1991-2001.
"""
import matplotlib.pyplot as plt
import numpy as np
import statsmodels.api as sm

from regbook import COLORS, Generated, figure_path, use_book_style

# <<runge>>
import numpy as np

runge = lambda x: 1.0 / (1.0 + 25.0 * x ** 2)
grid = np.linspace(-1.0, 1.0, 2001)

for m in (5, 11, 21, 31):
    equal = np.linspace(-1.0, 1.0, m)                          # m equally spaced nodes
    cheb = np.cos((2 * np.arange(m) + 1) * np.pi / (2 * m))     # m Chebyshev nodes
    err = [np.max(np.abs(np.polyval(np.polyfit(z, runge(z), m - 1), grid) - runge(grid)))
           for z in (equal, cheb)]
    print(f"{m:2d} nodes, degree {m - 1:2d}:  equally spaced {err[0]:10.4f}   Chebyshev {err[1]:8.5f}")
# <</runge>>

# <<extrapolate>>
import statsmodels.api as sm

co2 = sm.datasets.co2.load_pandas().data["co2"]
monthly = co2.resample("MS").mean().dropna()
y = monthly.to_numpy()
t = (monthly.index.year + (monthly.index.month - 0.5) / 12).to_numpy()

train = t < 1991.0                                    # fit on 1958-1990, predict 1991-2001
lo, hi = t[train].min(), t[train].max()
u = lambda s: 2 * (s - lo) / (hi - lo) - 1


def design(s, d):
    return np.column_stack([np.vander(u(s), d + 1, increasing=True),
                            np.cos(2 * np.pi * s), np.sin(2 * np.pi * s),
                            np.cos(4 * np.pi * s), np.sin(4 * np.pi * s)])


for d in (2, 3, 4, 6, 8):
    beta, *_ = np.linalg.lstsq(design(t[train], d), y[train], rcond=None)
    inside = y[train] - design(t[train], d) @ beta
    outside = design(t[~train], d) @ beta - y[~train]
    print(f"degree {d}:  residual sd {np.sqrt(np.mean(inside ** 2)):.3f} ppm inside,"
          f"   root mean squared error {np.sqrt(np.mean(outside ** 2)):9.3f} ppm outside,"
          f"   error in Dec 2001 {outside[-1]:10.2f} ppm")
# <</extrapolate>>

# ---- assertions -------------------------------------------------------------
eq_err, ch_err = {}, {}
for m in (5, 11, 21, 31):
    equal = np.linspace(-1.0, 1.0, m)
    cheb = np.cos((2 * np.arange(m) + 1) * np.pi / (2 * m))
    eq_err[m] = np.max(np.abs(np.polyval(np.polyfit(equal, runge(equal), m - 1), grid) - runge(grid)))
    ch_err[m] = np.max(np.abs(np.polyval(np.polyfit(cheb, runge(cheb), m - 1), grid) - runge(grid)))
    # the interpolant reproduces the function at its own nodes (up to degree 20; by
    # degree 30 the monomial basis is too ill conditioned for the coefficients to be
    # computed accurately, which is the point of the section)
    if m <= 21:
        assert np.allclose(np.polyval(np.polyfit(equal, runge(equal), m - 1), equal), runge(equal), atol=1e-6)
# equally spaced nodes get worse, Chebyshev nodes get better
assert eq_err[11] < eq_err[21] < eq_err[31]
assert ch_err[11] > ch_err[21] > ch_err[31]

inside_sd, outside_rmse, dec2001 = {}, {}, {}
for d in (2, 3, 4, 6, 8):
    beta, *_ = np.linalg.lstsq(design(t[train], d), y[train], rcond=None)
    inside_sd[d] = np.sqrt(np.mean((y[train] - design(t[train], d) @ beta) ** 2))
    out = design(t[~train], d) @ beta - y[~train]
    outside_rmse[d] = np.sqrt(np.mean(out ** 2))
    dec2001[d] = out[-1]
# the fit inside improves monotonically with the degree, the fit outside does not
assert inside_sd[2] > inside_sd[3] > inside_sd[4] > inside_sd[6] > inside_sd[8]
assert outside_rmse[8] > outside_rmse[6] > 3 * outside_rmse[4]

gen = Generated("ch42", "runge_extrapolation", prefix="run")
for m in (11, 21, 31):
    gen.num(f"eq{m}", eq_err[m], 2)
    gen.num(f"ch{m}", ch_err[m], 5)
gen.int("ntrain", int(train.sum()))
gen.int("ntest", int((~train).sum()))
for d in (2, 3, 4, 6, 8):
    gen.num(f"in{d}", inside_sd[d], 3)
    gen.num(f"out{d}", outside_rmse[d], 2)
    gen.num(f"dec{d}", dec2001[d], 1)
gen.write()

# ---- figure -----------------------------------------------------------------
use_book_style()
fig, axes = plt.subplots(1, 2, figsize=(5.8, 2.5))
ax = axes[0]
ax.plot(grid, runge(grid), color=COLORS["ink"], linewidth=1.0, label=r"$1/(1+25x^2)$")
for m, colour in ((11, COLORS["accent"]), (21, COLORS["second"])):
    equal = np.linspace(-1.0, 1.0, m)
    ax.plot(grid, np.polyval(np.polyfit(equal, runge(equal), m - 1), grid), color=colour,
            linewidth=1.0, label=f"{m} equally spaced nodes")
    ax.scatter(equal, runge(equal), s=6, color=colour, linewidths=0, zorder=3)
ax.set_ylim(-1.2, 1.6)
ax.set_xlabel("x")
ax.set_title("(a) Runge's phenomenon")
ax.legend(frameon=False, loc="upper left", fontsize=6.5)
ax = axes[1]
ax.plot(t, y, color=COLORS["muted"], linewidth=0.5, label="monthly mean")
future = np.linspace(1991.0, 2002.0, 200)
for d, colour in ((2, COLORS["third"]), (4, COLORS["accent"]), (6, COLORS["second"])):
    beta, *_ = np.linalg.lstsq(design(t[train], d), y[train], rcond=None)
    ax.plot(future, design(future, d) @ beta, color=colour, linewidth=1.0, label=f"degree {d}")
ax.axvline(1991.0, color=COLORS["grid"], linewidth=0.6, zorder=0)
ax.set_ylim(305, 395)
ax.set_xlabel("year")
ax.set_ylabel(r"CO$_2$ (ppm)")
ax.set_title("(b) extrapolating past 1990")
ax.legend(frameon=False, loc="upper left", fontsize=6.5)
fig.tight_layout()
fig.savefig(figure_path("ch42", "runge_extrapolation"))
