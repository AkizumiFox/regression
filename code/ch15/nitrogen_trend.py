"""Chapter 15, Section 4: orthogonal contrasts and trend analysis in the nitrogen trial.

Two complete sets of orthogonal contrasts decompose the between-rates sum of squares: the
orthogonal polynomials in the rate, and "control against fertilized" followed by polynomial
contrasts among the four fertilized rates. The cubic and quartic pieces together are the
lack of fit of a quadratic response curve. Data as in nitrogen_trial.py.
"""
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats

from regbook import COLORS, Generated, figure_path, use_book_style

# <<data>>
rates = np.array([0, 40, 80, 120, 160])            # kg N per hectare
g, m = len(rates), 6                                # 5 rates, 6 plots per rate
true_means = 6.8 - 2.7 * np.exp(-0.018 * rates)     # unknown to the analyst
rng = np.random.default_rng(1515)
group = np.repeat(np.arange(g), m)                  # plot i received rate group[i]
y = np.round(true_means[group] + rng.normal(0, 0.5, g * m), 2)   # yield, t/ha
n = len(y)
# <</data>>

# <<trend>>
means = np.array([y[group == k].mean() for k in range(g)])
nu = n - g
s2 = np.sum((y - means[group]) ** 2) / nu
ss_between = m * np.sum((means - means.mean()) ** 2)

C_poly = np.array([[-2, -1, 0, 1, 2],               # linear
                   [2, -1, -2, -1, 2],              # quadratic
                   [-1, 2, 0, -2, 1],               # cubic
                   [1, -4, 6, -4, 1]], float)       # quartic
assert np.allclose(C_poly @ C_poly.T, np.diag(np.diag(C_poly @ C_poly.T)))   # orthogonal

def contrast_ss(c):
    return (c @ means) ** 2 / (c @ c / m)           # balanced: sum c_k^2 / n_k = c'c / m

ss_poly = np.array([contrast_ss(c) for c in C_poly])
F_poly = ss_poly / s2
p_poly = stats.f.sf(F_poly, 1, nu)
for name, ss, F, p in zip(["linear", "quadratic", "cubic", "quartic"], ss_poly, F_poly, p_poly):
    print(f"{name:9s} SS = {ss:8.4f}  F = {F:7.2f}  p = {p:.2g}")
print(f"sum = {ss_poly.sum():.4f}   between rates = {ss_between:.4f}")
# <</trend>>

assert np.isclose(ss_poly.sum(), ss_between)
assert np.isclose(F_poly.mean(), ss_between / (g - 1) / s2)   # F is the average of the four

# <<lack>>
ss_lof = ss_poly[2] + ss_poly[3]                    # cubic + quartic, 2 df
F_lof = (ss_lof / 2) / s2
X2 = np.column_stack([np.ones(n), rates[group], rates[group] ** 2])
coef, *_ = np.linalg.lstsq(X2, y, rcond=None)       # quadratic regression on the rate
sse_quad = np.sum((y - X2 @ coef) ** 2)
print(f"lack of fit of the quadratic: SS = {ss_lof:.4f}, F = {F_lof:.3f},"
      f" p = {stats.f.sf(F_lof, 2, nu):.2f}")
print(f"quadratic regression SSE = {sse_quad:.4f} = lack of fit + within = {ss_lof + s2 * nu:.4f}")
# <</lack>>

assert np.isclose(sse_quad, ss_lof + s2 * nu)
ss_reg_quad = np.sum((X2 @ coef - y.mean()) ** 2)
assert np.isclose(ss_reg_quad, ss_poly[0] + ss_poly[1])
p_lof = stats.f.sf(F_lof, 2, nu)

# <<second>>
C_ctrl = np.array([[-4, 1, 1, 1, 1],                # control against the four fertilized
                   [0, -3, -1, 1, 3],               # linear among 40, ..., 160
                   [0, 1, -1, -1, 1],               # quadratic among them
                   [0, -1, 3, -3, 1]], float)       # cubic among them
ss_ctrl = np.array([contrast_ss(c) for c in C_ctrl])
print("second set:", ss_ctrl.round(4), " sum =", round(ss_ctrl.sum(), 4))
# <</second>>

assert np.allclose(C_ctrl @ C_ctrl.T, np.diag(np.diag(C_ctrl @ C_ctrl.T)))
assert np.allclose(C_ctrl.sum(axis=1), 0) and np.allclose(C_poly.sum(axis=1), 0)
assert np.isclose(ss_ctrl.sum(), ss_between)
F_ctrl = ss_ctrl / s2
p_ctrl = stats.f.sf(F_ctrl, 1, nu)
# the two bases span the same space: the change of basis is invertible
T = np.linalg.lstsq(C_poly.T, C_ctrl.T, rcond=None)[0]
assert np.allclose(C_poly.T @ T, C_ctrl.T)

# fitted quadratic and linear curves in the rate (for the figure)
X1 = X2[:, :2]
coef1, *_ = np.linalg.lstsq(X1, y, rcond=None)

gen = Generated("ch15", "nitrogen_trend")
names = ["lin", "quad", "cub", "quart"]
for k, nm in enumerate(names):
    if ss_poly[k] < 1e-3:
        gen.num(f"ss_{nm}", ss_poly[k], 1, sci=True)
    else:
        gen.num(f"ss_{nm}", ss_poly[k], 4)
    if F_poly[k] < 0.01:                            # quartic: print F in scientific form
        gen.num(f"F_{nm}", F_poly[k], 1, sci=True)
    else:
        gen.num(f"F_{nm}", F_poly[k], 2)
    if p_poly[k] > 0.99:                            # never print a p-value as 1
        gen.num(f"p_{nm}", p_poly[k], 3)
    else:
        gen.text(f"p_{nm}", f"{p_poly[k]:.2g}")
    gen.num(f"share_{nm}", 100 * ss_poly[k] / ss_between, 1)
gen.num("ss_lof", ss_lof, 4)
gen.num("F_lof", F_lof, 3)
gen.num("p_lof", p_lof, 2)
gen.num("sse_quad", sse_quad, 4)
for k, nm in enumerate(["ctrl", "flin", "fquad", "fcub"]):
    gen.num(f"ss_{nm}", ss_ctrl[k], 4)
    gen.num(f"F_{nm}", F_ctrl[k], 2)
    gen.text(f"p_{nm}", f"{p_ctrl[k]:.2g}")
gen.num("b0", coef[0], 3)
gen.num("b1", coef[1], 5)
gen.num("b2", coef[2], 7, sci=False)
gen.write()

# ---- figure: means, linear and quadratic fits; the four single-df pieces -------------------
use_book_style()
fig, axes = plt.subplots(1, 2, figsize=(5.8, 2.5), gridspec_kw={"width_ratios": [1.5, 1]})
ax = axes[0]
se = np.sqrt(s2 / m)
ax.errorbar(rates, means, yerr=2 * se, fmt="o", ms=4, color=COLORS["accent"], capsize=2,
            lw=1, label=r"mean $\pm 2$ se")
xs = np.linspace(0, 160, 200)
ax.plot(xs, coef1[0] + coef1[1] * xs, color=COLORS["muted"], ls="--", lw=1, label="linear")
ax.plot(xs, coef[0] + coef[1] * xs + coef[2] * xs ** 2, color=COLORS["second"], lw=1.3,
        label="quadratic")
ax.set_xticks(rates)
ax.set_xlabel("nitrogen rate (kg N/ha)")
ax.set_ylabel("yield (t/ha)")
ax.set_title("(a) rate means and polynomial fits")
ax.legend(frameon=False, loc="lower right")
ax = axes[1]
labels = ["lin.", "quad.", "cub.", "quart."]
ax.bar(range(4), ss_poly, color=[COLORS["accent"], COLORS["second"], COLORS["muted"],
                                   COLORS["muted"]], width=0.6)
ax.axhline(s2 * stats.f.ppf(0.95, 1, nu), color=COLORS["ink"], lw=0.8, ls=":")
ax.set_yscale("log")
ax.set_xticks(range(4))
ax.set_xticklabels(labels)
ax.set_ylabel("sum of squares (log scale)")
ax.set_title("(b) single-df pieces")
fig.tight_layout()
fig.savefig(figure_path("ch15", "nitrogen_trend"))
