"""Chapter 13, Section 2: Scheffe intervals for regression slopes (2009 US state data).

Model: murder rate on poverty, single-parent percentage and urban percentage, with an
intercept, 50 states (District of Columbia excluded, as in Chapter 6). Public-domain
data shipped with statsmodels (statsmodels.datasets.statecrime).
"""
import matplotlib.pyplot as plt
import numpy as np
import statsmodels.api as sm
from matplotlib.patches import Rectangle
from scipy import stats

from regbook import COLORS, Generated, figure_path, use_book_style

data = sm.datasets.statecrime.load_pandas().data
data = data.drop(index="District of Columbia")

# <<fit>>
y = data["murder"].to_numpy()
X = np.column_stack([np.ones(len(y)), data[["poverty", "single", "urban"]]])
n, p = X.shape
nu = n - p                                            # error degrees of freedom
XtX_inv = np.linalg.inv(X.T @ X)
beta_hat = XtX_inv @ X.T @ y
s2 = np.sum((y - X @ beta_hat) ** 2) / nu
se = np.sqrt(s2 * np.diag(XtX_inv))
# <</fit>>

# <<multipliers>>
alpha, slopes = 0.05, [1, 2, 3]                       # poverty, single, urban
k = len(slopes)
mult = {
    "unadjusted": stats.t.ppf(1 - alpha / 2, nu),
    "Bonferroni": stats.t.ppf(1 - alpha / (2 * k), nu),
    "Scheffe": np.sqrt(k * stats.f.ppf(1 - alpha, k, nu)),   # q = 3 slopes
}
for name, c in mult.items():
    lo, hi = beta_hat[slopes] - c * se[slopes], beta_hat[slopes] + c * se[slopes]
    print(f"{name:10s} c = {c:.3f}  poverty [{lo[0]:.3f}, {hi[0]:.3f}]")
# <</multipliers>>

# <<snoop>>
L = np.eye(p)[:, slopes]                              # Lambda: the three slopes
W = L.T @ XtX_inv @ L
d = L.T @ beta_hat
F = d @ np.linalg.solve(W, d) / (k * s2)              # F test that all slopes are 0
a_star = np.linalg.solve(W, d)                        # the most significant combination
t_star = (a_star @ d) / np.sqrt(s2 * a_star @ W @ a_star)
print(f"F = {F:.2f};  largest |t| over all combinations = {t_star:.3f} = sqrt(3F)")
# <</snoop>>

t_pov = beta_hat[1] / se[1]
assert np.isclose(t_star ** 2, k * F)
# no combination can beat a_star (random directions)
rng = np.random.default_rng(20260919)
A = rng.normal(size=(20_000, k))
t_all = np.abs(A @ d) / np.sqrt(s2 * np.einsum("ij,jk,ik->i", A, W, A))
assert t_all.max() <= t_star + 1e-9
# poverty: significant unadjusted and with Bonferroni, not with Scheffe
assert mult["unadjusted"] < mult["Bonferroni"] < abs(t_pov) < mult["Scheffe"]
# the F test rejects iff some Scheffe interval excludes zero
assert (F > stats.f.ppf(1 - alpha, k, nu)) == (t_star > mult["Scheffe"])

# two-slope family (poverty, single) for the figure
L2 = np.eye(p)[:, [1, 2]]
W2 = L2.T @ XtX_inv @ L2
c2_sch = np.sqrt(2 * stats.f.ppf(1 - alpha, 2, nu))
c2_bon = stats.t.ppf(1 - alpha / 4, nu)
c2_t = mult["unadjusted"]
theta = np.linspace(0, 2 * np.pi, 721)
chol = np.linalg.cholesky(W2)
ell = beta_hat[[1, 2]][:, None] + c2_sch * np.sqrt(s2) * chol @ np.vstack([np.cos(theta), np.sin(theta)])
# the ellipse's shadows on the axes are the Scheffe intervals
for j in range(2):
    assert np.isclose(ell[j].max(), beta_hat[1 + j] + c2_sch * se[1 + j], rtol=1e-4)
    assert np.isclose(ell[j].min(), beta_hat[1 + j] - c2_sch * se[1 + j], rtol=1e-4)
# and along any direction a: shadow = estimate +- c * se(a'beta_hat)
a = np.array([1.0, 1.0])
sh = a @ ell
se_a = np.sqrt(s2 * a @ W2 @ a)
assert np.isclose(sh.max(), a @ beta_hat[[1, 2]] + c2_sch * se_a, rtol=1e-4)

gen = Generated("ch13", "scheffe_crime", prefix="schc")
gen.int("n", n)
gen.int("nu", nu)
gen.num("bpov", beta_hat[1], 4)
gen.num("sepov", se[1], 4)
gen.num("tpov", t_pov, 3)
gen.num("bsing", beta_hat[2], 4)
gen.num("burb", beta_hat[3], 4)
gen.num("ct", mult["unadjusted"], 3)
gen.num("cbon", mult["Bonferroni"], 3)
gen.num("csch", mult["Scheffe"], 3)
for name, key in (("unadjusted", "t"), ("Bonferroni", "bon"), ("Scheffe", "sch")):
    c = mult[name]
    gen.num(f"lo{key}", beta_hat[1] - c * se[1], 3)
    gen.num(f"hi{key}", beta_hat[1] + c * se[1], 3)
gen.num("F", F, 2)
gen.num("tstar", t_star, 3)
astar = a_star / a_star[np.argmax(np.abs(a_star))]
gen.text("astar", ",\\ ".join(f"{v:.2f}".replace("-", "\\ensuremath{-}") for v in astar))
gen.num("c2sch", c2_sch, 3)
gen.num("c2bon", c2_bon, 3)
gen.write()

use_book_style()
fig, ax = plt.subplots(figsize=(5.6, 3.0))
b = beta_hat[[1, 2]]
ax.plot(ell[0], ell[1], color=COLORS["accent"], label="95% confidence ellipse")
for c, col, ls, lab in ((c2_sch, COLORS["accent"], "--", "Scheffé intervals (shadows)"),
                        (c2_bon, COLORS["second"], "-", "Bonferroni rectangle"),
                        (c2_t, COLORS["muted"], ":", "unadjusted rectangle")):
    w = c * se[[1, 2]]
    ax.add_patch(Rectangle(b - w, 2 * w[0], 2 * w[1], fill=False, edgecolor=col,
                           linestyle=ls, linewidth=0.9, label=lab))
# tangent lines for the direction a = (1, 1): the Scheffe interval for the sum of slopes
xs = np.linspace(ell[0].min() - 0.05, ell[0].max() + 0.05, 2)
for sgn in (-1, 1):
    level = a @ b + sgn * c2_sch * se_a
    ax.plot(xs, level - xs, color=COLORS["third"], linewidth=0.7, linestyle="-.",
            label="interval for sum of slopes" if sgn > 0 else None)
ax.plot(*b, "o", color=COLORS["ink"], markersize=3)
ax.set_xlim(ell[0].min() - 0.06, ell[0].max() + 0.06)
ax.set_ylim(ell[1].min() - 0.06, ell[1].max() + 0.06)
ax.set_xlabel("coefficient of poverty")
ax.set_ylabel("coefficient of single parenthood")
ax.legend(frameon=False, fontsize=7, loc="upper left", bbox_to_anchor=(1.02, 1.0))
fig.savefig(figure_path("ch13", "scheffe_shadow"))
