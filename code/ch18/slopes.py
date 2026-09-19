"""Chapter 18, Section 3: unequal slopes and the Johnson-Neyman region.

Synthetic randomized comparison of two versions of a reading programme. The covariate is
the child's reading fluency before the programme (words per minute); the response is the
gain in fluency. The new version helps weak readers more, so the lines are not parallel.
"""
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats

from regbook import COLORS, Generated, figure_path, use_book_style

# <<data>>
rng = np.random.default_rng(31)
n_per = 24
arm = rng.permutation(np.repeat([0, 1], n_per))      # 0 = current version, 1 = new version
x = np.round(rng.uniform(40, 120, 2 * n_per))        # baseline fluency
y = np.where(arm == 0, 14 + 0.06 * x, 36 - 0.16 * x) + rng.normal(0, 4, 2 * n_per)
y = np.round(y, 1)
n = len(y)
# <</data>>

# <<fits>>
def line(k):
    """Least squares line within arm k, and the arm's summaries."""
    xs, ys = x[arm == k], y[arm == k]
    Sxx = np.sum((xs - xs.mean()) ** 2)
    b = np.sum((xs - xs.mean()) * (ys - ys.mean())) / Sxx
    a = ys.mean() - b * xs.mean()
    sse = np.sum((ys - a - b * xs) ** 2)
    return a, b, xs.mean(), Sxx, len(xs), sse

(a0, b0, xb0, S0, n0, e0), (a1, b1, xb1, S1, n1, e1) = line(0), line(1)
nu = n - 4
s2 = (e0 + e1) / nu                                   # separate-lines model
b_pool = (b0 * S0 + b1 * S1) / (S0 + S1)
ss_slopes = S0 * S1 / (S0 + S1) * (b1 - b0) ** 2      # sum of squares for parallelism
F_slopes = ss_slopes / s2
print(f"slopes: current {b0:.4f}, new {b1:.4f};  F = {F_slopes:.2f} on 1 and {nu} df,"
      f" p = {stats.f.sf(F_slopes, 1, nu):.2g}")
# <</fits>>

# the general formula of the theorem: sum_k b_k^2 S_k - (sum_k b_k S_k)^2 / sum_k S_k
ss_general = b0 ** 2 * S0 + b1 ** 2 * S1 - (b0 * S0 + b1 * S1) ** 2 / (S0 + S1)
assert np.isclose(ss_general, ss_slopes)
# check against a direct fit of the parallel and separate lines models
D = np.eye(2)[arm]
X_par = np.column_stack([D, x])
X_sep = np.column_stack([D, D * x[:, None]])
rss = lambda X: np.sum((y - X @ np.linalg.lstsq(X, y, rcond=None)[0]) ** 2)
assert np.isclose(rss(X_par) - rss(X_sep), ss_slopes)
assert np.isclose(rss(X_sep), s2 * nu)
p_slopes = stats.f.sf(F_slopes, 1, nu)
assert p_slopes < 0.01

# <<johnson-neyman>>
def difference(x0):
    """Estimated new-minus-current difference at x0 and its variance divided by sigma^2."""
    d = (a1 - a0) + (b1 - b0) * x0
    v = 1 / n0 + (x0 - xb0) ** 2 / S0 + 1 / n1 + (x0 - xb1) ** 2 / S1
    return d, v

def significant_region(crit2):
    """Roots of d(x)^2 - crit2 * s2 * v(x) = 0, a quadratic in x."""
    A = (b1 - b0) ** 2 - crit2 * s2 * (1 / S0 + 1 / S1)
    B = 2 * (a1 - a0) * (b1 - b0) + 2 * crit2 * s2 * (xb0 / S0 + xb1 / S1)
    C = (a1 - a0) ** 2 - crit2 * s2 * (1 / n0 + 1 / n1 + xb0 ** 2 / S0 + xb1 ** 2 / S1)
    return np.sort(np.roots([A, B, C]).real)

t2 = stats.t.ppf(0.975, nu) ** 2                      # pointwise, 5% level
w2 = 2 * stats.f.ppf(0.95, 2, nu)                     # simultaneous over all x (Scheffe, q = 2)
x_cross = -(a1 - a0) / (b1 - b0)
print(f"lines cross at x = {x_cross:.1f}")
print("pointwise 5% region: difference significant outside", significant_region(t2).round(1))
print("simultaneous region: difference significant outside", significant_region(w2).round(1))
# <</johnson-neyman>>

lo_t, hi_t = significant_region(t2)
lo_w, hi_w = significant_region(w2)
for r in (lo_t, hi_t):
    d, v = difference(r)
    assert np.isclose(d ** 2, t2 * s2 * v)
assert lo_w < lo_t < x_cross < hi_t < hi_w
assert x.min() < lo_t < x.max() < hi_t     # the upper region lies beyond the data
# only the lower region contains data beyond the simultaneous boundary
n_low_t = int(np.sum(x < lo_t))
n_high_t = int(np.sum(x > hi_t))
n_low_w = int(np.sum(x < lo_w))
n_high_w = int(np.sum(x > hi_w))
d60, v60 = difference(60.0)
d100, v100 = difference(100.0)

gen = Generated("ch18", "slopes")
gen.int("n", n)
gen.int("nu", nu)
gen.num("b0", b0, 4)
gen.num("b1", b1, 4)
gen.num("a0", a0, 2)
gen.num("a1", a1, 2)
gen.num("bpool", b_pool, 4)
gen.num("ssslopes", ss_slopes, 1)
gen.num("s2", s2, 2)
gen.num("s", np.sqrt(s2), 2)
gen.num("F", F_slopes, 1)
gen.num("p", p_slopes, 4)
gen.num("xcross", x_cross, 1)
gen.num("lot", lo_t, 1)
gen.num("hit", hi_t, 1)
gen.num("low", lo_w, 1)
gen.num("hiw", hi_w, 1)
gen.num("t2", t2, 3)
gen.num("w2", w2, 3)
gen.int("nlowt", n_low_t)
gen.int("nhight", n_high_t)
gen.int("nloww", n_low_w)
gen.int("nhighw", n_high_w)
gen.num("d60", d60, 2)
gen.num("se60", np.sqrt(s2 * v60), 2)
gen.num("d100", d100, 2)
gen.num("se100", np.sqrt(s2 * v100), 2)
gen.num("xmin", x.min(), 0)
gen.num("xmax", x.max(), 0)
gen.write()

# ---- figure: the two lines, and the estimated difference with its bands -------
use_book_style()
fig, axes = plt.subplots(1, 2, figsize=(5.8, 2.6))
ax = axes[0]
grid = np.linspace(x.min() - 3, x.max() + 3, 200)
for k, (a, b, lab, c, mk) in enumerate([(a0, b0, "current", COLORS["accent"], "o"),
                                         (a1, b1, "new", COLORS["second"], "s")]):
    ax.scatter(x[arm == k], y[arm == k], s=12, marker=mk, color=c, alpha=0.85, linewidths=0, label=lab)
    ax.plot(grid, a + b * grid, color=c)
ax.set_xlabel("baseline fluency $x$")
ax.set_ylabel("gain $y$")
ax.set_title("(a) separate lines")
ax.legend(frameon=False, fontsize=7)
ax = axes[1]
d, v = difference(grid)
se = np.sqrt(s2 * v)
ax.fill_between(grid, d - np.sqrt(w2) * se, d + np.sqrt(w2) * se, color=COLORS["grid"], alpha=0.6,
                linewidth=0, label="simultaneous")
ax.plot(grid, d - np.sqrt(t2) * se, color=COLORS["accent"], linewidth=0.8, linestyle="--",
        label="pointwise")
ax.plot(grid, d + np.sqrt(t2) * se, color=COLORS["accent"], linewidth=0.8, linestyle="--")
ax.plot(grid, d, color=COLORS["ink"])
ax.axhline(0, color=COLORS["muted"], linewidth=0.7)
ax.axvline(lo_t, color=COLORS["accent"], linewidth=0.6, linestyle=":")   # hi_t lies beyond the data
ax.set_xlabel("baseline fluency $x$")
ax.set_ylabel("new $-$ current")
ax.set_title("(b) difference and 95% bands")
ax.legend(frameon=False, fontsize=7, loc="upper right")
fig.tight_layout()
fig.savefig(figure_path("ch18", "johnson_neyman"))
