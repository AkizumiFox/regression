"""Chapter 21, Section 3: weighted least squares with estimated weights.

(a) Engel's food expenditure data (public domain, statsmodels.datasets.engel): a power
variance function sigma_i^2 = sigma^2 income_i^gamma is estimated by regressing the log
squared OLS residuals on log income, and the regression is refitted by weighted least
squares, with model-based and sandwich standard errors.
(b) A simulation with ten design points x = 1, ..., 10, m replicates at each, and error
standard deviation proportional to x. Weights estimated from the within-group sample
variances are compared with weights from a fitted variance function.
"""
import matplotlib.pyplot as plt
import numpy as np
import statsmodels.api as sm
from scipy import stats

from regbook import COLORS, Generated, figure_path, use_book_style

# <<engel>>
df = sm.datasets.engel.load_pandas().data
income, food = df["income"].to_numpy(), df["foodexp"].to_numpy()
n = len(food)
X = np.column_stack([np.ones(n), income])


def wls(X, y, w):
    """Weighted least squares: coefficients, model-based and HC3 covariances."""
    sw = np.sqrt(w)
    Q, R = np.linalg.qr(X * sw[:, None])              # OLS on the whitened data
    b = np.linalg.solve(R, Q.T @ (y * sw))
    ew = (y - X @ b) * sw                             # whitened residuals
    # leverages in the weighted geometry
    hw = np.sum(Q**2, axis=1)
    Rinv = np.linalg.inv(R)
    bread = Rinv @ Rinv.T                             # (X'WX)^{-1}
    s2w = ew @ ew / (len(y) - X.shape[1])
    meat = (Q * (ew / (1 - hw))[:, None]).T @ (Q * (ew / (1 - hw))[:, None])
    return b, s2w * bread, Rinv @ meat @ Rinv.T


b_ols, V_ols, V_ols_hc3 = wls(X, food, np.ones(n))
e = food - X @ b_ols
Z = np.column_stack([np.ones(n), np.log(income)])
# variance ~ income^gamma
gamma = np.linalg.lstsq(Z, np.log(e**2), rcond=None)[0][1]
b_fw, V_fw, V_fw_hc3 = wls(X, food, income**-gamma)
print(f"estimated power gamma = {gamma:.2f}")
print(f"OLS  slope {b_ols[1]:.4f}  model se {np.sqrt(V_ols[1, 1]):.4f}  "
      f"HC3 se {np.sqrt(V_ols_hc3[1, 1]):.4f}")
print(f"FWLS slope {b_fw[1]:.4f}  model se {np.sqrt(V_fw[1, 1]):.4f}  "
      f"HC3 se {np.sqrt(V_fw_hc3[1, 1]):.4f}")
# <</engel>>

# iterate: re-estimate gamma from the weighted fit's residuals until it settles
g, history = gamma, [gamma]
for _ in range(50):
    b_it, _, _ = wls(X, food, income**-g)
    g_new = np.linalg.lstsq(Z, np.log((food - X @ b_it) ** 2), rcond=None)[0][1]
    history.append(g_new)
    if abs(g_new - g) < 1e-8:
        break
    g = g_new
assert abs(g_new - g) < 1e-8
b_it, V_it, _ = wls(X, food, income**-g)

# checks against statsmodels
fit = sm.WLS(food, X, weights=income**-gamma).fit()
assert np.allclose(fit.params, b_fw) and np.allclose(fit.bse, np.sqrt(np.diag(V_fw)))
assert np.allclose(np.sqrt(np.diag(fit.get_robustcov_results("HC3").cov_params())), np.sqrt(np.diag(V_fw_hc3)))
assert np.allclose(np.sqrt(np.diag(V_ols_hc3)), sm.OLS(food, X).fit(cov_type="HC3").bse)

gen = Generated("ch21", "feasible_wls", prefix="fwls")
gen.num("gamma", gamma, 2)
gen.num("gamma_iter", g, 2)
gen.int("iterations", len(history) - 1)
gen.num("ols_slope", b_ols[1], 4)
gen.num("ols_se", np.sqrt(V_ols[1, 1]), 4)
gen.num("ols_hc3", np.sqrt(V_ols_hc3[1, 1]), 4)
gen.num("fw_slope", b_fw[1], 4)
gen.num("fw_se", np.sqrt(V_fw[1, 1]), 4)
gen.num("fw_hc3", np.sqrt(V_fw_hc3[1, 1]), 4)
gen.num("it_slope", b_it[1], 4)
gen.num("fw_int", b_fw[0], 1)

# ---- (b) replicated design: estimated weights ---------------------------------------
# <<replicates>>
rng = np.random.default_rng(2106)
# ten design points, sd proportional to x
xg = np.arange(1.0, 11.0)


def replicate_study(m, reps):
    """Variance of four slope estimators relative to WLS with true weights."""
    x = np.repeat(xg, m)
    Xr = np.column_stack([np.ones(len(x)), x])
    Y = 1 + 0.5 * x + rng.normal(size=(reps, len(x))) * x

    def slopes(W):
        A = np.einsum("ri,ij,ik->rjk", W, Xr, Xr)     # X'WX for each replicate
        c = np.einsum("ri,ij,ri->rj", W, Xr, Y)
        return np.linalg.solve(A, c[..., None])[..., 0], A

    b_o, _ = slopes(np.ones_like(Y))
    b_w, _ = slopes(np.tile(1 / x**2, (reps, 1)))
    # within-group variances
    s2g = Y.reshape(reps, len(xg), m).var(axis=2, ddof=1)
    Wg = np.repeat(1 / s2g, m, axis=1)
    b_g, A_g = slopes(Wg)
    # smooth model: log s2 on log x
    Zg = np.column_stack([np.ones(len(xg)), np.log(xg)])
    coef = np.linalg.lstsq(Zg, np.log(s2g).T, rcond=None)[0]
    b_m, _ = slopes(np.repeat(np.exp(-(Zg @ coef).T), m, axis=1))
    var = np.array([np.var(b[:, 1]) for b in (b_o, b_w, b_g, b_m)])
    # nominal interval for b_g
    res = Y - b_g @ Xr.T
    se = np.sqrt((Wg * res**2).sum(axis=1) / (len(x) - 2)
                 * np.linalg.inv(A_g)[:, 1, 1])
    q = stats.t.ppf(0.975, len(x) - 2)
    cover = np.mean(np.abs(b_g[:, 1] - 0.5) <= q * se)
    return var / var[1], cover


for m in (2, 3, 5, 10, 30):
    rel, cover = replicate_study(m, 2000)
    print(f"m = {m:2d}: variance relative to WLS  OLS {rel[0]:.2f}  "
          f"group weights {rel[2]:.2f}  fitted weights {rel[3]:.2f};  "
          f"coverage {cover:.2f}")
# <</replicates>>

ms = [2, 3, 4, 5, 7, 10, 15, 20, 30]
study = {m: replicate_study(m, 20_000) for m in ms}
for m in (2, 3, 5, 10, 30):
    rel, cover = study[m]
    gen.num(f"rel_ols_{m}", rel[0], 2)
    gen.num(f"rel_group_{m}", rel[2], 2)
    gen.num(f"rel_model_{m}", rel[3], 2)
    gen.num(f"cover_{m}", cover, 2)
# group weights worse than OLS at m = 2
assert study[2][0][2] > study[2][0][0]
assert study[30][0][2] < 1.1 and study[3][0][3] < 1.2
assert study[2][1] < 0.6 and study[30][1] > 0.9
gen.write()

# ---- figure -------------------------------------------------------------------------
use_book_style()
fig, axes = plt.subplots(1, 2, figsize=(5.8, 2.4))
ax = axes[0]
for k, label, color, ls in [(0, "OLS", COLORS["muted"], "-"), (2, "group-variance weights", COLORS["second"], "-"),
                            (3, "fitted variance function", COLORS["third"], "--")]:
    ax.plot(ms, [study[m][0][k] for m in ms], color=color, ls=ls, marker="o", markersize=2.5, label=label)
ax.axhline(1, color=COLORS["accent"], linewidth=0.8, label="WLS, true weights")
ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xticks(ms[::2], [str(m) for m in ms[::2]])
ax.set_yticks([1, 2, 5], ["1", "2", "5"])
ax.minorticks_off()
ax.set_xlabel("replicates per design point, $m$")
ax.set_ylabel("variance / WLS variance")
ax.legend(frameon=False, fontsize=6.5)
ax.set_title("(a) efficiency of the slope")
ax = axes[1]
ax.plot(ms, [study[m][1] for m in ms], color=COLORS["second"], marker="o", markersize=2.5)
ax.axhline(0.95, color=COLORS["grid"], linewidth=0.8)
ax.set_xscale("log")
ax.set_xticks(ms[::2], [str(m) for m in ms[::2]])
ax.minorticks_off()
ax.set_ylim(0.4, 1.0)
ax.set_xlabel("replicates per design point, $m$")
ax.set_ylabel("coverage")
ax.set_title("(b) nominal 95% interval, group weights")
fig.tight_layout()
fig.savefig(figure_path("ch21", "feasible_wls"))
