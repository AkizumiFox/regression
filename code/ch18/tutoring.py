"""Chapter 18, Sections 1-3: analysis of covariance for a synthetic randomized experiment.

Forty students are randomized, ten to each of four study formats; a pretest score x is
recorded before randomization and a final-examination score y afterwards. The data are
synthetic (fixed seed) and generated from a parallel-lines model, so the analysis can be
checked against the truth.
"""
import itertools

import matplotlib.pyplot as plt
import numpy as np
import statsmodels.formula.api as smf
import pandas as pd
from scipy import stats

from regbook import COLORS, Generated, figure_path, use_book_style

# <<data>>
rng = np.random.default_rng(84)
g, m = 4, 10                                          # four formats, ten students each
n = g * m
group = rng.permutation(np.repeat(np.arange(g), m))   # the randomization
x = np.round(rng.normal(60, 10, n))                   # pretest, measured before assignment
mu_true = np.array([70.0, 74.0, 71.0, 77.0])
y = np.round(mu_true[group] + 0.8 * (x - 60) + rng.normal(0, 6, n), 1)

Z = np.eye(g)[group]                                  # indicator matrix of the formats
n_k = Z.sum(axis=0)
x_bar, y_bar = Z.T @ x / n_k, Z.T @ y / n_k           # group means
print("pretest means", x_bar.round(1), " final means", y_bar.round(1))
# <</data>>


# <<ancova>>
def within(v):
    """(I - M) v: deviations from the group means (M projects onto C(Z))."""
    return v - (Z.T @ v / n_k)[group]

E_xx, E_xy, E_yy = within(x) @ within(x), within(x) @ within(y), within(y) @ within(y)
beta = E_xy / E_xx                                    # slope from within-group variation only
sse = E_yy - E_xy ** 2 / E_xx                         # ANCOVA residual sum of squares
nu = n - g - 1

T_xx = np.sum((x - x.mean()) ** 2)                    # the same with I - M0, M0 = mean only
T_xy = np.sum((x - x.mean()) * (y - y.mean()))
T_yy = np.sum((y - y.mean()) ** 2)
sse0 = T_yy - T_xy ** 2 / T_xx                        # reduced model: one line for everyone
ss_trt = sse0 - sse                                   # adjusted treatment sum of squares
F_adj = (ss_trt / (g - 1)) / (sse / nu)
F_unadj = ((T_yy - E_yy) / (g - 1)) / (E_yy / (n - g))
print(f"slope {beta:.4f};  SSE: ANOVA {E_yy:.2f} on {n - g} df, ANCOVA {sse:.2f} on {nu} df")
print(f"treatments: unadjusted F = {F_unadj:.2f} (p = {stats.f.sf(F_unadj, g - 1, n - g):.3f}),"
      f" adjusted F = {F_adj:.2f} (p = {stats.f.sf(F_adj, g - 1, nu):.4f})")
# <</ancova>>

# Theorem checks: the estimates agree with a direct least squares fit of [Z, x]
W = np.column_stack([Z, x])
coef, *_ = np.linalg.lstsq(W, y, rcond=None)
assert np.isclose(coef[-1], beta)
assert np.allclose(coef[:g], y_bar - beta * x_bar)    # X beta-hat = M (y - x gamma-hat)
assert np.isclose(np.sum((y - W @ coef) ** 2), sse)
df = pd.DataFrame({"y": y, "x": x, "grp": group})
fit_sm = smf.ols("y ~ C(grp) + x", df).fit()
assert np.isclose(fit_sm.ssr, sse)
fit_red = smf.ols("y ~ x", df).fit()
assert np.isclose(fit_red.ssr, sse0)
assert np.isclose(fit_sm.compare_f_test(fit_red)[0], F_adj)
assert F_adj > F_unadj                                # adjustment sharpened the comparison
assert sse / nu < E_yy / (n - g)                      # and reduced the error mean square
p_unadj = stats.f.sf(F_unadj, g - 1, n - g)
p_adj = stats.f.sf(F_adj, g - 1, nu)
assert p_unadj > 0.05 and p_adj < 0.01
F_slope = (E_xy ** 2 / E_xx) / (sse / nu)
s_anova, s_ancova = np.sqrt(E_yy / (n - g)), np.sqrt(sse / nu)

# <<adjusted>>
s2 = sse / nu
x0 = x.mean()                                         # the covariate value for adjustment
adj = y_bar - beta * (x_bar - x0)                     # adjusted means
se_adj = np.sqrt(s2 * (1 / n_k + (x_bar - x0) ** 2 / E_xx))
pairs = list(itertools.combinations(range(g), 2))
q = stats.studentized_range.ppf(0.95, g, nu) / np.sqrt(2)   # Tukey-Kramer multiplier
c_s = np.sqrt((g - 1) * stats.f.ppf(0.95, g - 1, nu))       # Scheffe multiplier
for k, l in pairs:
    d = adj[k] - adj[l]
    se_d = np.sqrt(s2 * (1 / n_k[k] + 1 / n_k[l] + (x_bar[k] - x_bar[l]) ** 2 / E_xx))
    print(f"{k + 1}-{l + 1}: raw {y_bar[k] - y_bar[l]:6.2f}  adjusted {d:6.2f}  se {se_d:.2f}"
          f"  Tukey-type +-{q * se_d:.2f}  Scheffe +-{c_s * se_d:.2f}")
print("adjusted means", adj.round(2), " standard errors", se_adj.round(2))
# <</adjusted>>

# adjusted means from the fitted model at x = x0, and their covariance from (X^T X)^{-1}
C = np.linalg.inv(W.T @ W)
L = np.column_stack([np.eye(g), np.full(g, x0)])
assert np.allclose(L @ coef, adj)
assert np.allclose(np.sqrt(s2 * np.diag(L @ C @ L.T)), se_adj)
assert np.isclose(np.sum(n_k * adj) / n, y.mean())    # weighted average of adjusted means
# the adjusted F statistic is the Scheffe maximum over all contrasts of adjusted means
H = np.column_stack([np.eye(g - 1), -np.ones(g - 1)])  # contrasts mu_k - mu_g
Lh = H @ L
u = Lh @ coef
F_check = u @ np.linalg.solve(Lh @ C @ Lh.T, u) / (g - 1) / s2
assert np.isclose(F_check, F_adj)
diffs = {}
for k, l in pairs:
    d = adj[k] - adj[l]
    se_d = np.sqrt(s2 * (1 / n_k[k] + 1 / n_k[l] + (x_bar[k] - x_bar[l]) ** 2 / E_xx))
    diffs[(k, l)] = (d, se_d)
    lam = L[k] - L[l]
    assert np.isclose(se_d, np.sqrt(s2 * lam @ C @ lam))
sig_tk = [(k, l) for (k, l), (d, se) in diffs.items() if abs(d) > q * se]
sig_sch = [(k, l) for (k, l), (d, se) in diffs.items() if abs(d) > c_s * se]
print("significant (Tukey-type):", sig_tk, " (Scheffe):", sig_sch)
assert set(sig_tk) == {(1, 2), (2, 3)}
assert set(sig_sch) == {(1, 2), (2, 3)}

# <<parallel>>
X_sep = np.column_stack([Z, Z * x[:, None]])          # a separate line for each format
E_xx_k = np.array([within(x)[group == k] @ within(x)[group == k] for k in range(g)])
E_xy_k = np.array([within(x)[group == k] @ within(y)[group == k] for k in range(g)])
slopes = E_xy_k / E_xx_k
sse_sep = sse - (np.sum(E_xy_k ** 2 / E_xx_k) - E_xy ** 2 / E_xx)
F_par = ((sse - sse_sep) / (g - 1)) / (sse_sep / (n - 2 * g))
print("separate slopes", slopes.round(3),
      f" F = {F_par:.2f}, p = {stats.f.sf(F_par, g - 1, n - 2 * g):.3f}")
# <</parallel>>

coef_sep, *_ = np.linalg.lstsq(X_sep, y, rcond=None)
assert np.allclose(coef_sep[g:], slopes)
assert np.isclose(np.sum((y - X_sep @ coef_sep) ** 2), sse_sep)
assert np.isclose(beta, np.sum(E_xx_k * slopes) / np.sum(E_xx_k))  # pooled slope is a weighted mean
p_par = stats.f.sf(F_par, g - 1, n - 2 * g)
assert p_par > 0.2

gen = Generated("ch18", "tutoring")
gen.int("n", n)
gen.int("nu", nu)
for k in range(g):
    gen.num(f"xbar{k + 1}", x_bar[k], 1)
    gen.num(f"ybar{k + 1}", y_bar[k], 2)
    gen.num(f"adj{k + 1}", adj[k], 2)
    gen.num(f"seadj{k + 1}", se_adj[k], 2)
    gen.num(f"slope{k + 1}", slopes[k], 3)
gen.num("xmean", x0, 2)
gen.num("ymean", y.mean(), 2)
gen.num("Exx", E_xx, 1)
gen.num("Exy", E_xy, 1)
gen.num("Eyy", E_yy, 2)
gen.num("Txx", T_xx, 1)
gen.num("Txy", T_xy, 1)
gen.num("Tyy", T_yy, 2)
gen.num("beta", beta, 4)
gen.num("sse", sse, 2)
gen.num("sse0", sse0, 2)
gen.num("sstrt", ss_trt, 2)
gen.num("ssbetween", T_yy - E_yy, 2)
gen.num("ssslope", E_xy ** 2 / E_xx, 1)
gen.num("Fadj", F_adj, 2)
gen.num("Funadj", F_unadj, 2)
gen.num("padj", p_adj, 4)
gen.num("punadj", p_unadj, 3)
gen.num("Fslope", F_slope, 1)
gen.num("sanova", s_anova, 2)
gen.num("sancova", s_ancova, 2)
gen.num("msanova", E_yy / (n - g), 1)
gen.num("msancova", sse / nu, 1)
gen.num("q", q, 3)
gen.num("cs", c_s, 3)
for (k, l), (d, se_d) in diffs.items():
    gen.num(f"d{k + 1}{l + 1}", d, 2)
    gen.num(f"se{k + 1}{l + 1}", se_d, 2)
    gen.num(f"raw{k + 1}{l + 1}", y_bar[k] - y_bar[l], 2)
    gen.num(f"hwq{k + 1}{l + 1}", q * se_d, 2)
    gen.num(f"hws{k + 1}{l + 1}", c_s * se_d, 2)
gen.num("ssesep", sse_sep, 2)
gen.num("Fpar", F_par, 2)
gen.num("ppar", p_par, 2)
gen.num("unadjse", s_anova * np.sqrt(2 / m), 2)
gen.num("Bxx", T_xx - E_xx, 1)
gen.num("Bxy", T_xy - E_xy, 1)
gen.num("shift12", beta * (x_bar[0] - x_bar[1]), 2)
gen.num("shift2", beta * (x_bar[1] - x0), 2)
gen.num("Txy2", T_xy ** 2 / T_xx, 2)
gen.num("Exy2", E_xy ** 2 / E_xx, 2)
gen.num("sd_x", np.sqrt(E_xx / (n - g)), 2)
gen.write()

# ---- figure: parallel lines and adjusted means ------------------------------
use_book_style()
cols = [COLORS["accent"], COLORS["second"], COLORS["third"], COLORS["thread"]]
marks = ["o", "s", "^", "D"]
styles = ["-", "--", "-", "-"]
fig, ax = plt.subplots(figsize=(4.8, 3.1))
xs = np.array([x.min() - 2, x.max() + 2])
for k in range(g):
    sel = group == k
    ax.scatter(x[sel], y[sel], s=14, marker=marks[k], color=cols[k], alpha=0.85,
               linewidths=0, label=f"format {k + 1}")
    ax.plot(xs, y_bar[k] + beta * (xs - x_bar[k]), color=cols[k], linewidth=1.0, linestyle=styles[k],
            zorder=3 if k == 1 else 2)
    ax.plot([x_bar[k]], [y_bar[k]], marker=marks[k], markersize=7, markerfacecolor="white",
            markeredgecolor=cols[k], markeredgewidth=1.3, linestyle="none")
    ax.plot([x0], [adj[k]], marker=marks[k], markersize=6, color=cols[k], linestyle="none")
ax.axvline(x0, color=COLORS["muted"], linewidth=0.8, linestyle=":", zorder=0)
ax.set_xlabel("pretest score $x$")
ax.set_ylabel("final score $y$")
ax.legend(frameon=False, loc="upper left", fontsize=7)
fig.tight_layout()
fig.savefig(figure_path("ch18", "ancova_lines"))
