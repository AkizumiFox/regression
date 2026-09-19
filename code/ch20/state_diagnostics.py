"""Chapter 20, Sections 3 and 4: deletion and influence diagnostics for the 2009 state data.

All 50 states and the District of Columbia. Response: murder rate per 100,000. Regressors:
poverty rate, percentage of single-parent households, percentage urban, with an intercept.
Public-domain data shipped with statsmodels (statsmodels.datasets.statecrime).
"""
import matplotlib.pyplot as plt
import numpy as np
import statsmodels.api as sm
from scipy import stats

from regbook import COLORS, Generated, figure_path, use_book_style

data = sm.datasets.statecrime.load_pandas().data           # District of Columbia kept
y = data["murder"].to_numpy()
X = np.column_stack([np.ones(len(y)), data[["poverty", "single", "urban"]]])
names = np.array(data.index)

# <<diagnostics>>
n, p = X.shape
Q, R = np.linalg.qr(X)
beta = np.linalg.solve(R, Q.T @ y)
e = y - X @ beta                                   # raw residuals
h = np.sum(Q ** 2, axis=1)                         # leverages
s2 = e @ e / (n - p)
r = e / np.sqrt(s2 * (1 - h))                      # internally studentized
s2_del = (e @ e - e ** 2 / (1 - h)) / (n - p - 1)  # s^2 without case i
t = e / np.sqrt(s2_del * (1 - h))                  # externally studentized

p_single = 2 * stats.t.sf(np.abs(t), n - p - 1)    # one case, chosen in advance
p_bonf = np.minimum(1, n * p_single)               # Bonferroni over all n cases
crit = stats.t.isf(0.05 / (2 * n), n - p - 1)      # 5% Bonferroni critical value
for i in np.argsort(-np.abs(t))[:3]:
    print(f"{names[i]:22s} h = {h[i]:.3f}  r = {r[i]:6.3f}  t = {t[i]:6.3f}"
          f"  Bonferroni p = {p_bonf[i]:.4f}")
print(f"critical value |t| > {crit:.3f}")
# <</diagnostics>>

# ---- the deletion formulas against brute-force refits -------------------------------
for i in range(n):
    keep = np.arange(n) != i
    b_i, *_ = np.linalg.lstsq(X[keep], y[keep], rcond=None)
    K = np.linalg.inv(X.T @ X)
    assert np.allclose(beta - b_i, K @ X[i] * e[i] / (1 - h[i]))                  # coefficients
    assert np.allclose(X @ beta - X @ b_i, X @ K @ X[i] * e[i] / (1 - h[i]))      # fitted values
    assert np.isclose(y[i] - X[i] @ b_i, e[i] / (1 - h[i]))                       # prediction residual
    sse_i = np.sum((y[keep] - X[keep] @ b_i) ** 2)
    assert np.isclose(sse_i / (n - p - 1), s2_del[i])
    # mean-shift model: the t statistic of the indicator of case i is t_i
    Z = np.column_stack([X, np.eye(n)[i]])
    fit = sm.OLS(y, Z).fit()
    assert np.isclose(fit.tvalues[-1], t[i])
assert np.all(np.abs(r) <= np.sqrt(n - p))
assert np.allclose(t, r * np.sqrt((n - p - 1) / (n - p - r ** 2)))
dc = int(np.where(names == "District of Columbia")[0][0])
ms = int(np.where(names == "Mississippi")[0][0])
assert np.argmax(np.abs(t)) == dc and p_bonf[dc] < 0.01
assert np.abs(t[dc]) > crit and np.abs(t[ms]) < crit

# after deleting the District of Columbia
keep = np.arange(n) != dc
X1, y1, nm1 = X[keep], y[keep], names[keep]
n1 = n - 1
Q1, R1 = np.linalg.qr(X1)
e1 = y1 - Q1 @ (Q1.T @ y1)
h1 = np.sum(Q1 ** 2, axis=1)
t1 = e1 / np.sqrt((e1 @ e1 - e1 ** 2 / (1 - h1)) / (n1 - p - 1) * (1 - h1))
k1 = int(np.argmax(np.abs(t1)))
pb1 = min(1.0, n1 * 2 * stats.t.sf(abs(t1[k1]), n1 - p - 1))
assert nm1[k1] == "Louisiana" and pb1 > 0.05

# <<influence>>
cook = r ** 2 * h / (p * (1 - h))                  # Cook's distance
dffits = t * np.sqrt(h / (1 - h))
K = np.linalg.inv(X.T @ X)
C = K @ X.T                                        # row j: weights of y in beta_j
dfbetas = (C * e / (1 - h)).T / np.sqrt(np.outer(s2_del, np.diag(K)))
for i in np.argsort(-cook)[:3]:
    print(f"{names[i]:22s} D = {cook[i]:.3f}  DFFITS = {dffits[i]:6.3f}"
          f"  DFBETAS(poverty) = {dfbetas[i, 1]:6.3f}")
print("conventional cut-offs: D > 4/n =", round(4 / n, 3),
      " |DFFITS| > 2 sqrt(p/n) =", round(2 * np.sqrt(p / n), 3),
      " |DFBETAS| > 2/sqrt(n) =", round(2 / np.sqrt(n), 3))
# <</influence>>

# identities of the influence section
fitted_change = np.array([np.sum((X @ K @ X[i] * e[i] / (1 - h[i])) ** 2) for i in range(n)])
assert np.allclose(cook, fitted_change / (p * s2))
assert np.allclose(cook, dffits ** 2 * s2_del / (p * s2))


def resid_on(v, Z):
    coef, *_ = np.linalg.lstsq(Z, v, rcond=None)
    return v - Z @ coef


xt = resid_on(X[:, 1], X[:, [0, 2, 3]])            # poverty, residualized on the others
ell = xt ** 2 / (xt @ xt)                          # partial leverage for poverty
assert np.isclose(ell.sum(), 1)
h_without = np.sum(np.linalg.qr(X[:, [0, 2, 3]])[0] ** 2, axis=1)
assert np.allclose(h, h_without + ell)
assert np.allclose(dfbetas[:, 1], t * xt / (np.sqrt(xt @ xt) * np.sqrt(1 - h)))
assert np.allclose(C[1], xt / (xt @ xt))
# sup over directions of the standardized change is |DFFITS|
for i in [dc, ms]:
    d_i = K @ X[i] * e[i] / (1 - h[i])
    rng = np.random.default_rng(i)
    lam = rng.normal(size=(2000, p))
    ratio = (lam @ d_i) ** 2 / (s2_del[i] * np.einsum("ij,jk,ik->i", lam, K, lam))
    assert ratio.max() <= dffits[i] ** 2 + 1e-10
    assert np.isclose((X[i] @ d_i) ** 2 / (s2_del[i] * X[i] @ K @ X[i]), dffits[i] ** 2)
covratio = (s2_del / s2) ** p / (1 - h)
assert np.allclose(covratio, 1 / ((1 - h) * ((n - p - 1 + t ** 2) / (n - p)) ** p))
alaska = int(np.where(names == "Alaska")[0][0])
assert np.argmax(ell) == alaska

fit_all = sm.OLS(y, X).fit()
fit_no = sm.OLS(y1, X1).fit()
assert np.isclose(fit_all.params[1] - fit_no.params[1], (C[1, dc] * e[dc] / (1 - h[dc])))
f_med = stats.f.ppf(0.5, p, n - p)
f_dc = stats.f.cdf(cook[dc], p, n - p)            # beta_(i) sits on the boundary of this region
cook_mean_null = np.mean(h / (p * (1 - h)))       # average of E(D_i) under the model
assert cook_mean_null >= 1 / (n - p)

press_resid = e / (1 - h)
press = np.sum(press_resid ** 2)
assert np.argmax(np.abs(press_resid)) == dc
exp_press_ratio = np.mean(1 / (1 - h))            # E(PRESS) / (n sigma^2)
assert exp_press_ratio >= n / (n - p)

gen = Generated("ch20", "state_diagnostics", prefix="st")
gen.num("sse", e @ e, 2)
gen.num("press", press, 2)
gen.num("press_dc", press_resid[dc], 3)
gen.num("press_share", press_resid[dc] ** 2 / press, 2)
gen.num("exp_press_ratio", exp_press_ratio, 3)
gen.num("lower_ratio", n / (n - p), 3)
gen.num("oos_ratio", (n + p) / n, 3)
gen.int("n", n)
gen.int("p", p)
gen.num("h_dc", h[dc], 3)
gen.num("e_dc", e[dc], 3)
gen.num("s", np.sqrt(s2), 3)
gen.num("s_dc", np.sqrt(s2_del[dc]), 3)
gen.num("r_dc", r[dc], 3)
gen.num("t_dc", t[dc], 3)
gen.num("psingle_dc", p_single[dc], 1, sci=True)
gen.num("pbonf_dc", p_bonf[dc], 4)
gen.num("crit", crit, 3)
gen.num("t_ms", t[ms], 3)
gen.num("pbonf_ms", p_bonf[ms], 3)
gen.num("t_la", t1[k1], 3)
gen.num("pbonf_la", pb1, 3)
gen.num("D_dc", cook[dc], 3)
gen.num("D_ms", cook[ms], 3)
gen.num("dffits_dc", dffits[dc], 3)
gen.num("dfb_dc", dfbetas[dc, 1], 3)
gen.num("ell_dc", ell[dc], 3)
gen.num("ell_ak", ell[alaska], 3)
gen.num("dfb_ak", dfbetas[alaska, 1], 3)
gen.num("t_ak", t[alaska], 3)
gen.num("cut_D", 4 / n, 3)
gen.num("cut_dffits", 2 * np.sqrt(p / n), 3)
gen.num("cut_dfbetas", 2 / np.sqrt(n), 3)
gen.num("cut_h", 2 * p / n, 3)
gen.num("f_med", f_med, 3)
gen.num("f_dc", f_dc, 3)
gen.num("cook_mean", cook_mean_null, 4)
gen.num("cook_bound", 1 / (n - p), 4)
gen.num("covratio_dc", covratio[dc], 3)
gen.num("pov_all", fit_all.params[1], 4)
gen.num("pov_no", fit_no.params[1], 4)
gen.num("se_all", fit_all.bse[1], 4)
gen.num("se_no", fit_no.bse[1], 4)
gen.num("tpov_all", fit_all.tvalues[1], 2)
gen.num("tpov_no", fit_no.tvalues[1], 2)
gen.num("single_all", fit_all.params[2], 4)
gen.num("single_no", fit_no.params[2], 4)
gen.int("n_flag_D", int(np.sum(cook > 4 / n)))
gen.int("n_flag_dffits", int(np.sum(np.abs(dffits) > 2 * np.sqrt(p / n))))
gen.int("n_flag_dfb", int(np.sum(np.abs(dfbetas[:, 1]) > 2 / np.sqrt(n))))
gen.int("n_flag_h", int(np.sum(h > 2 * p / n)))
gen.write()

# ---- influence plot: leverage against internally studentized residual -----------------
use_book_style()
fig, ax = plt.subplots(figsize=(5.2, 3.4))
hh = np.linspace(0.005, 0.6, 400)
for Dval, ls in [(0.5, ":"), (1.0, "--")]:
    rr = np.sqrt(Dval * p * (1 - hh) / hh)
    ax.plot(hh, rr, color=COLORS["muted"], linestyle=ls, linewidth=0.8)
    ax.plot(hh, -rr, color=COLORS["muted"], linestyle=ls, linewidth=0.8)
    ax.annotate(f"$D={Dval:g}$", (0.56, np.sqrt(Dval * p * 0.44 / 0.56)), fontsize=7,
                color=COLORS["muted"], ha="center", va="bottom")
ax.axvline(2 * p / n, color=COLORS["grid"], linewidth=0.8, zorder=0)
ax.axhline(0, color=COLORS["grid"], linewidth=0.6, zorder=0)
size = 12 + 160 * cook / cook.max()
ax.scatter(h, r, s=size, color=COLORS["accent"], alpha=0.7, linewidths=0)
for i, dx, dy in [(dc, -0.012, 0.25), (ms, 0.01, -0.45), (alaska, 0.008, 0.25),
                  (int(np.where(names == "New York")[0][0]), 0.008, -0.35)]:
    ax.annotate(names[i], (h[i], r[i]), xytext=(h[i] + dx, r[i] + dy), fontsize=7,
                ha="right" if dx < 0 else "left")
ax.set_xlim(0, 0.6)
ax.set_ylim(-3.6, 4.4)
ax.set_xlabel(r"leverage $h_{ii}$")
ax.set_ylabel(r"studentized residual $r_i$")
fig.tight_layout()
fig.savefig(figure_path("ch20", "influence_plot"))
