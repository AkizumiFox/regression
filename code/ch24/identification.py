"""Chapter 24, Section 2: what the data can and cannot identify, and Berkson error.

(a) In the normal structural model the distribution of (w, y) identifies the slope only up to
    the interval between the direct regression slope s_wy / s_ww and the reverse regression
    slope s_yy / s_wy. Computed for the blood-pressure study (one reading).
(b) Berkson error: concentrations set on a dial (w fixed), achieved concentration x = w + u.
    Least squares on w is unbiased, its t intervals keep their coverage, and s^2 estimates
    sigma^2 + beta^2 sigma_u^2.
"""
import matplotlib.pyplot as plt
from scipy import stats

from regbook import COLORS, Generated, figure_path, use_book_style

gen = Generated("ch24", "identification", prefix="idf")

# <<setup>>
import numpy as np

rng = np.random.default_rng(2401)
n = 300
age = rng.normal(55, 10, n)
x = 130 + 0.6 * (age - 55) + rng.normal(0, np.sqrt(108), n)   # long-run blood pressure, sd 12
w = x[:, None] + rng.normal(0, 9, (n, 2))                      # two clinic readings
y = -4 + 0.08 * x + 0.0 * age + rng.normal(0, 1.0, n)          # age has no effect given x
# <</setup>>

# <<interval>>
S = np.cov(w[:, 0], y)                       # sample covariance matrix of (reading, outcome)
s_ww, s_wy, s_yy = S[0, 0], S[0, 1], S[1, 1]
b_direct = s_wy / s_ww                       # assumes no measurement error (sigma_u^2 = 0)
b_reverse = s_yy / s_wy                      # assumes no equation error (sigma^2 = 0)
print(f"slopes compatible with the data: [{b_direct:.4f}, {b_reverse:.4f}]")

for b in [b_direct, 0.08, b_reverse]:        # the variances each slope would imply
    sx2 = s_wy / b
    print(f"b = {b:.4f}: sigma_x^2 = {sx2:7.2f}, sigma_u^2 = {s_ww - sx2:7.2f}, "
          f"sigma^2 = {s_yy - b * s_wy:.4f}")
# <</interval>>

r2 = s_wy ** 2 / (s_ww * s_yy)
assert b_direct < 0.08 < b_reverse
assert np.isclose(b_direct / b_reverse, r2)
# population values of the endpoints
pop_direct = 0.08 * 144 / 225
pop_reverse = (0.08 ** 2 * 144 + 1.0) / (0.08 * 144)
assert abs(b_direct - pop_direct) < 0.02 and abs(b_reverse - pop_reverse) < 0.05
gen.num("s_ww", s_ww, 2)
gen.num("s_wy", s_wy, 3)
gen.num("s_yy", s_yy, 3)
gen.num("b_direct", b_direct, 4)
gen.num("b_reverse", b_reverse, 4)
gen.num("r2", r2, 3)
sx2_true = s_wy / 0.08
gen.num("sx2_at_true", sx2_true, 1)
gen.num("su2_at_true", s_ww - sx2_true, 1)
gen.num("s2_at_true", s_yy - 0.08 * s_wy, 3)

# ---- (b) Berkson error ----------------------------------------------------------------------
# <<berkson>>
rng_b = np.random.default_rng(2410)
settings = np.repeat([10.0, 20.0, 30.0, 40.0, 50.0, 60.0], 5)      # dial settings, n = 30
Xd = np.column_stack([np.ones(30), settings])
reps = 20_000
achieved = settings + rng_b.normal(0, 4, (reps, 30))               # x = w + u, u independent of w
Y = 2 + 0.5 * achieved + rng_b.normal(0, 3, (reps, 30))
B = np.linalg.solve(Xd.T @ Xd, Xd.T @ Y.T)                          # least squares on the settings
s2 = np.sum((Y.T - Xd @ B) ** 2, axis=0) / 28
se = np.sqrt(s2 * np.linalg.inv(Xd.T @ Xd)[1, 1])
covered = np.abs(B[1] - 0.5) <= stats.t.ppf(0.975, 28) * se
print(f"mean slope {B[1].mean():.4f}, mean s^2 {s2.mean():.2f} (sigma^2 + beta^2 sigma_u^2 = 13), "
      f"coverage {covered.mean():.3f}")
# <</berkson>>
assert abs(B[1].mean() - 0.5) < 0.002
assert abs(s2.mean() - 13) < 0.1
assert abs(covered.mean() - 0.95) < 0.006
gen.num("berk_slope", B[1].mean(), 4)
gen.num("berk_s2", s2.mean(), 2)
gen.num("berk_cover", covered.mean(), 3)
gen.int("berk_reps", reps)
gen.write()

# ---- figure ----------------------------------------------------------------------------------
use_book_style()
fig, axes = plt.subplots(1, 2, figsize=(5.8, 2.55))
ax = axes[0]
ax.scatter(w[:, 0], y, s=6, color=COLORS["accent"], alpha=0.5, linewidths=0)
wbar, ybar = w[:, 0].mean(), y.mean()
xs = np.linspace(95, 165, 2)
for b in np.linspace(b_direct, b_reverse, 7):
    ax.plot(xs, ybar + b * (xs - wbar), color=COLORS["grid"], linewidth=0.7, zorder=0)
ax.plot(xs, ybar + b_direct * (xs - wbar), color=COLORS["second"], label="regress $y$ on $w$")
ax.plot(xs, ybar + b_reverse * (xs - wbar), color=COLORS["third"], label="regress $w$ on $y$, inverted")
ax.plot(xs, ybar + 0.08 * (xs - wbar), color=COLORS["ink"], linestyle="--", label="true slope")
ax.set_ylim(y.min() - 0.5, y.max() + 0.5)
ax.set_xlabel("reading $w$ (mmHg)")
ax.set_ylabel("outcome $y$")
ax.set_title("(a) lines the data cannot tell apart")
ax.legend(frameon=False, loc="upper left", fontsize=7)
ax = axes[1]
bs = np.linspace(b_direct, b_reverse, 200)
ax.plot(bs, b_direct / bs, color=COLORS["accent"], label=r"implied reliability $\lambda$")
ax.plot(bs, bs / b_reverse, color=COLORS["second"], label=r"implied $\rho^2_{xy}$")
ax.axvline(0.08, color=COLORS["ink"], linestyle="--", linewidth=0.9)
ax.set_xlabel("candidate slope $b$")
ax.set_ylabel("implied value")
ax.set_ylim(-0.02, 1.05)
ax.set_title("(b) what each slope would require")
ax.legend(frameon=False, loc="lower center", fontsize=7)
fig.tight_layout()
fig.savefig(figure_path("ch24", "identified_set"))
