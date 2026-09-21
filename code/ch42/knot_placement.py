"""Chapter 42, Section 4: where the knots go.

A simulated skewed design (fixed seed) with n = 200 points on [0, 1], most of them in
the left quarter. Five interior knots are placed equally spaced and at the quintiles of
the design, and the exact pointwise standard deviation of the fitted cubic spline,
sigma * sqrt(x(t)' (X'X)^{-1} x(t)), is compared. The second block shows a rank failure
when an interval holds fewer design points than the piece has parameters.
"""
import matplotlib.pyplot as plt
import numpy as np

from regbook import COLORS, Generated, figure_path, use_book_style

# <<design>>
import numpy as np

rng = np.random.default_rng(42)
n = 200
x = np.sort(rng.beta(1.5, 6.0, n))                      # a skewed design on [0, 1]
K, d = 5, 3

equal = np.linspace(x.min(), x.max(), K + 2)[1:-1]      # equally spaced knots
quantile = np.quantile(x, np.arange(1, K + 1) / (K + 1))


def spline_design(s, knots):
    """Truncated power basis for degree d with the given interior knots."""
    return np.column_stack([s ** j for j in range(d + 1)]
                           + [np.clip(s - k, 0.0, None) ** d for k in knots])


def counts(knots):
    edges = np.r_[-np.inf, knots, np.inf]
    return [int(((x >= a) & (x < b)).sum()) for a, b in zip(edges[:-1], edges[1:])]


print("equally spaced knots:", np.round(equal, 3), "  points per interval", counts(equal))
print("quantile knots:      ", np.round(quantile, 3), "  points per interval", counts(quantile))
# <</design>>

# <<variance>>
grid = np.linspace(x.min(), x.max(), 500)
for name, knots in (("equally spaced", equal), ("quantile", quantile)):
    X = spline_design(x, knots)
    G = np.linalg.inv(X.T @ X)
    Xg = spline_design(grid, knots)
    sd = np.sqrt(np.einsum("ij,jk,ik->i", Xg, G, Xg))   # in units of sigma
    print(f"{name:15s} knots: standard deviation of the fit at the median {np.interp(np.median(x), grid, sd):.3f},"
          f" largest {sd.max():.3f} at x = {grid[sd.argmax()]:.3f}")
# <</variance>>

# <<rank>>
beyond = np.r_[equal, 0.90]                             # one knot above every design point
for knots, label in ((equal, "equally spaced"), (quantile, "quantile"), (beyond, "one knot too far")):
    X = spline_design(x, knots)
    print(f"{label:17s}: {X.shape[1]} columns, rank {np.linalg.matrix_rank(X)},"
          f" condition number {np.linalg.cond(X):.2e}")
# <</rank>>

# ---- assertions -------------------------------------------------------------
sds = {}
for name, knots in (("equal", equal), ("quantile", quantile)):
    X = spline_design(x, knots)
    assert np.linalg.matrix_rank(X) == d + 1 + K        # both designs have full rank
    Xg = spline_design(grid, knots)
    sds[name] = np.sqrt(np.einsum("ij,jk,ik->i", Xg, np.linalg.inv(X.T @ X), Xg))
# quantile knots equalize the counts and cut the worst-case variance by a large factor
assert max(counts(quantile)) - min(counts(quantile)) <= 1
assert min(counts(equal)) < d + 1
assert sds["equal"].max() > 5 * sds["quantile"].max()
# but the quantile knots pay for it where the data are dense
assert sds["quantile"][np.argmin(np.abs(grid - np.median(x)))] > sds["equal"][np.argmin(np.abs(grid - np.median(x)))]

# a knot above every design point kills a column outright
X_dead = spline_design(x, beyond)
assert np.all(X_dead[:, -1] == 0.0)
assert np.linalg.matrix_rank(X_dead) == X_dead.shape[1] - 1
# both knot sets keep full rank: the price of a near-empty interval is paid in variance,
# not in rank, because the truncated power columns are not locally supported
cond_equal = np.linalg.cond(spline_design(x, equal))
cond_quantile = np.linalg.cond(spline_design(x, quantile))
assert np.linalg.matrix_rank(spline_design(x, equal)) == d + 1 + K

gen = Generated("ch42", "knot_placement", prefix="knt")
gen.int("n", n)
gen.int("K", K)
gen.text("counteq", ", ".join(str(c) for c in counts(equal)))
gen.text("countqt", ", ".join(str(c) for c in counts(quantile)))
for name in ("equal", "quantile"):
    gen.num(f"sdmax{name}", sds[name].max(), 3)
    gen.num(f"sdmed{name}", np.interp(np.median(x), grid, sds[name]), 3)
gen.num("ratio", sds["equal"].max() / sds["quantile"].max(), 1)
gen.num("condeq", cond_equal, 2, sci=True)
gen.num("condqt", cond_quantile, 2, sci=True)
gen.int("rankdead", int(np.linalg.matrix_rank(X_dead)))
gen.int("coldead", X_dead.shape[1])
gen.write()

# ---- figure -----------------------------------------------------------------
use_book_style()
fig, axes = plt.subplots(1, 2, figsize=(5.8, 2.5))
ax = axes[0]
ax.hist(x, bins=30, color=COLORS["grid"], edgecolor="none")
for k, colour, off in ((equal, COLORS["second"], -0.9), (quantile, COLORS["accent"], -2.1)):
    ax.scatter(k, np.full(K, off), s=14, color=colour, marker="|", linewidths=1.2, clip_on=False)
ax.text(0.75, -0.9, "equally spaced", color=COLORS["second"], fontsize=6.5, va="center")
ax.text(0.75, -2.1, "quantile", color=COLORS["accent"], fontsize=6.5, va="center")
ax.set_ylim(-3, 30)
ax.set_xlabel("x")
ax.set_ylabel("count")
ax.set_title("(a) a skewed design and two knot sets")
ax = axes[1]
for name, knots, colour in (("equally spaced", equal, COLORS["second"]),
                            ("quantile", quantile, COLORS["accent"])):
    key = "equal" if name.startswith("equal") else "quantile"
    ax.semilogy(grid, sds[key], color=colour, label=f"{name} knots")
ax.set_xlabel("x")
ax.set_ylabel(r"standard deviation / $\sigma$")
ax.set_title("(b) pointwise variability of the fit")
ax.legend(frameon=False, loc="upper left", fontsize=6.5)
fig.tight_layout()
fig.savefig(figure_path("ch42", "knot_placement"))
