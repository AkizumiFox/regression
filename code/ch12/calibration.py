"""Chapter 12, Section 5: calibration of a simulated assay, and the three shapes of a Fieller set.

The data are simulated with a fixed seed: nine standard concentrations 0, 1, ..., 8, each measured
three times, and an unknown sample measured m = 3 times.
"""
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats

from regbook import COLORS, Generated, figure_path, use_book_style

# <<assay>>
rng = np.random.default_rng(1206)
x = np.repeat(np.arange(9.0), 3)                   # standards, in triplicate
n = len(x)
y = 0.04 + 0.115 * x + 0.015 * rng.normal(size=n)  # absorbance readings
m = 3
y0 = 0.04 + 0.115 * 5.3 + 0.015 * rng.normal(size=m)   # the unknown sample, true x0 = 5.3


def calibrate(x, y, y0bar, m, level=0.95):
    """Fieller set {x0 : |y0bar - b0 - b1 x0| <= t s sqrt(1/m + 1/n + (x0 - xbar)^2 / Sxx)}."""
    n = len(x)
    xbar, Sxx = x.mean(), np.sum((x - x.mean()) ** 2)
    b1 = np.sum((x - xbar) * (y - y.mean())) / Sxx
    s2 = np.sum((y - y.mean() - b1 * (x - xbar)) ** 2) / (n - 2)
    t2 = stats.t.ppf(0.5 + level / 2, n - 2) ** 2
    g = y0bar - y.mean()
    # quadratic in d = x0 - xbar:  A d^2 - 2 B d + K <= 0
    A = b1**2 - t2 * s2 / Sxx
    B = b1 * g
    K = g**2 - t2 * s2 * (1 / m + 1 / n)
    disc = B**2 - A * K
    if disc < 0:
        return "whole line", (-np.inf, np.inf), xbar + g / b1
    r1, r2 = sorted([(B - np.sqrt(disc)) / A, (B + np.sqrt(disc)) / A])
    if A > 0:
        return "interval", (xbar + r1, xbar + r2), xbar + g / b1
    return "two rays", (xbar + r1, xbar + r2), xbar + g / b1   # (-inf, first] and [second, inf)


shape, (lo, hi), x0_hat = calibrate(x, y, y0.mean(), m)
print(f"estimate {x0_hat:.3f}; 95% calibration set: {shape} ({lo:.3f}, {hi:.3f})")
# <</assay>>

xbar, Sxx = x.mean(), np.sum((x - x.mean()) ** 2)
b1 = np.sum((x - xbar) * (y - y.mean())) / Sxx
b0 = y.mean() - b1 * xbar
s = np.sqrt(np.sum((y - b0 - b1 * x) ** 2) / (n - 2))
tq = stats.t.ppf(0.975, n - 2)
t_slope = b1 * np.sqrt(Sxx) / s
assert shape == "interval" and lo < x0_hat < hi and lo < 5.3 < hi
# the set is the inverse image of the prediction band for the mean of m readings
def inside(x0, x, y, y0bar, m):
    n = len(x)
    xb, Sx = x.mean(), np.sum((x - x.mean()) ** 2)
    c1 = np.sum((x - xb) * (y - y.mean())) / Sx
    ss = np.sqrt(np.sum((y - y.mean() - c1 * (x - xb)) ** 2) / (n - 2))
    q = stats.t.ppf(0.975, n - 2)
    return np.abs(y0bar - y.mean() - c1 * (x0 - xb)) <= q * ss * np.sqrt(1 / m + 1 / n + (x0 - xb) ** 2 / Sx)


xs = np.linspace(-50, 60, 22001)
assert np.array_equal(inside(xs, x, y, y0.mean(), m), (xs >= lo) & (xs <= hi))
# the delta-method interval, for comparison
se_delta = s / abs(b1) * np.sqrt(1 / m + 1 / n + (x0_hat - xbar) ** 2 / Sxx)
delta = (x0_hat - tq * se_delta, x0_hat + tq * se_delta)

# ---- a weak assay: slope not significant, so the set is unbounded -------------------------------
rng_w = np.random.default_rng(1210)
y_w = 0.04 + 0.002 * x + 0.015 * rng_w.normal(size=n)
b1_w = np.sum((x - xbar) * (y_w - y_w.mean())) / Sxx
s_w = np.sqrt(np.sum((y_w - y_w.mean() - b1_w * (x - xbar)) ** 2) / (n - 2))
t_slope_w = b1_w * np.sqrt(Sxx) / s_w
assert abs(t_slope_w) < tq
readings_w = {"near": y_w.mean() + 0.005, "far": y_w.mean() + 0.03}
res_w = {k: calibrate(x, y_w, v, m) for k, v in readings_w.items()}
assert res_w["near"][0] == "whole line" and res_w["far"][0] == "two rays"
for k, v in readings_w.items():
    sh, (l_, h_), _ = res_w[k]
    ins = inside(xs, x, y_w, v, m)
    if sh == "whole line":
        assert ins.all()
    else:
        assert np.array_equal(ins, (xs <= l_) | (xs >= h_))
print("weak assay:", {k: (v[0], np.round(v[1], 2)) for k, v in res_w.items()}, f"slope t = {t_slope_w:.2f}")

# ---- coverage: the Fieller set is exact, the delta-method interval is not -----------------------
def coverage(beta1, x0, reps=20_000, seed=1208):
    rs = np.random.default_rng(seed)
    Y = 0.04 + beta1 * x + 0.015 * rs.normal(size=(reps, n))
    Y0 = 0.04 + beta1 * x0 + 0.015 * rs.normal(size=(reps, m)).mean(axis=1)
    b1s = (Y - Y.mean(axis=1, keepdims=True)) @ (x - xbar) / Sxx
    ss = np.sqrt(np.sum((Y - Y.mean(axis=1, keepdims=True) - b1s[:, None] * (x - xbar)) ** 2, axis=1) / (n - 2))
    resid0 = Y0 - Y.mean(axis=1) - b1s * (x0 - xbar)
    fieller = np.abs(resid0) <= tq * ss * np.sqrt(1 / m + 1 / n + (x0 - xbar) ** 2 / Sxx)
    xh = xbar + (Y0 - Y.mean(axis=1)) / b1s
    sed = ss / np.abs(b1s) * np.sqrt(1 / m + 1 / n + (xh - xbar) ** 2 / Sxx)
    delta_cov = np.abs(xh - x0) <= tq * sed
    bounded = b1s**2 * Sxx / ss**2 > tq**2
    return fieller.mean(), delta_cov.mean(), bounded.mean()


strong = coverage(0.115, 5.3)
weak_slope = 3 * 0.015 / np.sqrt(Sxx)            # slope three standard errors from zero
weak = coverage(weak_slope, 5.3)
weaker = coverage(2 * 0.015 / np.sqrt(Sxx), 12.0)  # two standard errors, and x0 beyond the standards
for name, c in [("strong", strong), ("weak", weak), ("weaker, x0 = 12", weaker)]:
    print(f"{name:16s}: Fieller {c[0]:.4f}, delta method {c[1]:.4f}, bounded {c[2]:.4f}")
tol = 4 * np.sqrt(0.95 * 0.05 / 20_000)
assert all(abs(c[0] - 0.95) < tol for c in (strong, weak, weaker))
assert abs(strong[1] - 0.95) < 0.01 and weak[1] > 0.97 and weaker[1] < 0.945
# probability of a bounded set = power of the slope t test = P(|t(n-2, 3)| > tq)
power = stats.nct.sf(tq, n - 2, 3) + stats.nct.cdf(-tq, n - 2, 3)
assert abs(weak[2] - power) < 0.01

gen = Generated("ch12", "calibration", prefix="cal")
gen.int("n", n)
gen.int("m", m)
gen.num("b0", b0, 4)
gen.num("b1", b1, 4)
gen.num("s", s, 4)
gen.num("tslope", t_slope, 1)
gen.num("y0bar", y0.mean(), 4)
gen.num("x0hat", x0_hat, 3)
gen.num("lo", lo, 3)
gen.num("hi", hi, 3)
gen.num("dlo", delta[0], 3)
gen.num("dhi", delta[1], 3)
gen.num("tq", tq, 3)
gen.num("w:b1", b1_w, 4)
gen.num("w:tslope", t_slope_w, 2)
gen.num("w:ybar", y_w.mean(), 4)
gen.num("w:near", readings_w["near"], 4)
gen.num("w:far", readings_w["far"], 4)
gen.num("w:far:lo", res_w["far"][1][0], 1)
gen.num("w:far:hi", res_w["far"][1][1], 1)
gen.num("cov:strong:f", strong[0], 4)
gen.num("cov:strong:d", strong[1], 4)
gen.num("cov:weak:f", weak[0], 4)
gen.num("cov:weak:d", weak[1], 4)
gen.num("cov:weak:bounded", weak[2], 3)
gen.num("cov:weaker:f", weaker[0], 4)
gen.num("cov:weaker:d", weaker[1], 4)
gen.num("cov:weaker:bounded", weaker[2], 3)
gen.num("power", power, 3)
gen.num("weakslope", weak_slope, 5)
gen.write()

# ---- figure ----------------------------------------------------------------------------------
use_book_style()
fig, axes = plt.subplots(1, 2, figsize=(5.8, 2.6))


def panel(ax, x, y, levels, xlim, title, bar0=0.03):
    xb, Sx = x.mean(), np.sum((x - x.mean()) ** 2)
    c1 = np.sum((x - xb) * (y - y.mean())) / Sx
    c0 = y.mean() - c1 * xb
    ss = np.sqrt(np.sum((y - c0 - c1 * x) ** 2) / (n - 2))
    g = np.linspace(*xlim, 3000)
    w = tq * ss * np.sqrt(1 / m + 1 / n + (g - xb) ** 2 / Sx)
    ax.fill_between(g, c0 + c1 * g - w, c0 + c1 * g + w, color=COLORS["grid"], alpha=0.8, linewidth=0)
    ax.plot(g, c0 + c1 * g, color=COLORS["ink"], linewidth=0.8)
    ax.scatter(x, y, s=6, color=COLORS["ink"], linewidths=0, zorder=3)
    for k, (level, col) in enumerate(levels):
        ax.axhline(level, color=col, linewidth=0.9, linestyle=(0, (4, 2)))
        ins = inside(g, x, y, level, m)
        ax.plot(np.where(ins, g, np.nan), np.full_like(g, bar0 + 0.05 * k), color=col, linewidth=3,
                solid_capstyle="butt", transform=ax.get_xaxis_transform())
    ax.set_xlim(*xlim)
    ax.set_xlabel("concentration")
    ax.set_title(title)
    return c0, c1, ss


panel(axes[0], x, y, [(y0.mean(), COLORS["accent"])], (-1, 9), "(a) a good assay")
axes[0].set_ylabel("reading")
inset = axes[0].inset_axes([0.58, 0.2, 0.38, 0.34])
panel(inset, x, y, [(y0.mean(), COLORS["accent"])], (5.0, 5.65), "", bar0=0.14)
inset.set_ylim(y0.mean() - 0.035, y0.mean() + 0.035)
inset.set_xlabel("")
inset.tick_params(labelsize=6)
axes[0].indicate_inset_zoom(inset, edgecolor=COLORS["muted"])
panel(axes[1], x, y_w, [(readings_w["near"], COLORS["accent"]), (readings_w["far"], COLORS["second"])],
      (-75, 45), "(b) a weak assay")
fig.tight_layout()
fig.savefig(figure_path("ch12", "calibration"))
