"""Chapter 33, Section 5: growth curves for Grunfeld's panel.

Log gross investment of 11 US firms in each of the years 1935-1954, fitted with a
straight line whose intercept and slope vary from firm to firm. Two-stage
(method of moments) estimates of the variance components, best linear unbiased
predictors of the firm curves, and a prediction for a firm not in the data.
"""
import matplotlib.pyplot as plt
import numpy as np
import statsmodels.api as sm
from scipy import stats
from scipy.optimize import minimize

from regbook import COLORS, Generated, figure_path, use_book_style

# <<data>>
df = sm.datasets.grunfeld.load_pandas().data
firms = sorted(df["firm"].unique())
years = np.sort(df["year"].unique())
N, m = len(firms), len(years)
t = (years - years.mean()) / 10                      # decades from the middle of the record
X = np.column_stack([np.ones(m), t])                 # the same design for every firm
p = X.shape[1]
Y = np.column_stack([np.log(df[df["firm"] == f].sort_values("year")["invest"].to_numpy())
                     for f in firms])                # m x N: one column per firm
# <</data>>
assert Y.shape == (m, N)

# <<twostage>>
XtXinv = np.linalg.inv(X.T @ X)
B = XtXinv @ X.T @ Y                                 # p x N: each firm's own least squares line
beta_hat = B.mean(axis=1)                            # the average of the firm lines
E = Y - X @ B                                        # within-firm residuals
sigma2 = np.sum(E**2) / (N * (m - p))                # pooled within-firm variance
Sb = np.cov(B, ddof=1)                               # spread of the firm lines
D_hat = Sb - sigma2 * XtXinv                         # ... minus what estimation noise adds
print("beta_hat", beta_hat.round(4), " sigma", round(np.sqrt(sigma2), 4))
print("D_hat\n", D_hat.round(4))
print("implied sd of firm intercepts and slopes:", np.sqrt(np.diag(D_hat)).round(4))
# <</twostage>>
assert np.min(np.linalg.eigvalsh(D_hat)) > 0         # a proper covariance matrix here

Sigma = X @ D_hat @ X.T + sigma2 * np.eye(m)
# ordinary least squares is generalized least squares: C(Sigma X) is inside C(X)
Sinv = np.linalg.inv(Sigma)
assert np.allclose(np.linalg.solve(X.T @ Sinv @ X, X.T @ Sinv @ X @ XtXinv @ X.T),
                   XtXinv @ X.T)
assert np.allclose(beta_hat, XtXinv @ X.T @ Y.mean(axis=1))
cov_beta = (D_hat + sigma2 * XtXinv) / N             # exact covariance of beta_hat
se_beta = np.sqrt(np.diag(cov_beta))

# <<blup>>
W = D_hat @ np.linalg.inv(D_hat + sigma2 * XtXinv)   # the shrinkage matrix
B_blup = beta_hat[:, None] + W @ (B - beta_hat[:, None])
u_direct = D_hat @ X.T @ np.linalg.solve(Sigma, Y - (X @ beta_hat)[:, None])
print("shrinkage matrix\n", W.round(4))
print("largest difference from the direct formula:",
      np.abs(B_blup - beta_hat[:, None] - u_direct).max())
# <</blup>>
assert np.allclose(B_blup - beta_hat[:, None], u_direct)
assert 0 < np.linalg.eigvalsh(W).min() and np.linalg.eigvalsh(W).max() < 1
# the predictors are shrunk towards the population line
assert np.all(np.abs(B_blup - beta_hat[:, None]) <= np.abs(B - beta_hat[:, None]) + 1e-9)
shrink_slope = np.std(B_blup[1], ddof=1) / np.std(B[1], ddof=1)

# <<prediction>>
C = sigma2 * XtXinv                                  # noise in a firm's own line
V_known = np.linalg.inv(np.linalg.inv(D_hat) + np.linalg.inv(C)) \
    + (np.eye(p) - W) @ C / N                        # error of the predicted firm curve
V_new = D_hat + cov_beta                             # error for a firm not in the data
x_end = np.array([1.0, t[-1]])
se_curve_known = np.sqrt(x_end @ V_known @ x_end)
se_curve_new = np.sqrt(x_end @ V_new @ x_end)
pred_new = x_end @ beta_hat
print(f"1954 curve: se for a firm in the data {se_curve_known:.4f}, "
      f"for a new firm {se_curve_new:.4f}")
# <</prediction>>
assert np.allclose(V_known, W @ C + (np.eye(p) - W) @ C / N)
assert np.allclose(W @ D_hat, D_hat @ X.T @ np.linalg.solve(Sigma, X) @ D_hat)
assert np.allclose(W @ C, np.linalg.inv(np.linalg.inv(D_hat) + np.linalg.inv(C)))
assert se_curve_new > 3 * se_curve_known

# the prediction error variances checked against a simulation from the fitted model
rs = np.random.default_rng(404)
LD = np.linalg.cholesky(D_hat)
errs_known, errs_new = [], []
for _ in range(20000):
    Us = LD @ rs.standard_normal((p, N))
    Ys = X @ (beta_hat[:, None] + Us) + np.sqrt(sigma2) * rs.standard_normal((m, N))
    Bs = XtXinv @ X.T @ Ys
    bs = Bs.mean(axis=1)
    Bp = bs[:, None] + W @ (Bs - bs[:, None])
    errs_known.append(x_end @ (Bp[:, 0] - beta_hat - Us[:, 0]))
    u_new = LD @ rs.standard_normal(p)
    errs_new.append(x_end @ (bs - beta_hat - u_new))
assert abs(np.std(errs_known) / se_curve_known - 1) < 0.04, np.std(errs_known)
assert abs(np.std(errs_new) / se_curve_new - 1) < 0.04, np.std(errs_new)


# ---- the same model by restricted maximum likelihood --------------------------
def reml_objective(th):
    D = np.array([[np.exp(th[1]), th[3]], [th[3], np.exp(th[2])]])
    S = np.exp(th[0]) * np.eye(m) + X @ D @ X.T
    if np.min(np.linalg.eigvalsh(S)) <= 1e-10:
        return 1e6
    Si = np.linalg.inv(S)
    b = np.linalg.solve(N * X.T @ Si @ X, X.T @ Si @ Y.sum(axis=1))
    R = Y - (X @ b)[:, None]
    _, ld = np.linalg.slogdet(S)
    _, ldi = np.linalg.slogdet(N * X.T @ Si @ X)
    return N * ld + np.sum(R * (Si @ R)) + ldi


opt = minimize(reml_objective, [-2.0, -1.0, -2.0, 0.0], method="Nelder-Mead",
               options={"maxiter": 20000, "xatol": 1e-10, "fatol": 1e-10})
D_reml = np.array([[np.exp(opt.x[1]), opt.x[3]], [opt.x[3], np.exp(opt.x[2])]])
sigma2_reml = np.exp(opt.x[0])
print("REML: sigma", round(np.sqrt(sigma2_reml), 4), " sd of intercepts and slopes",
      np.sqrt(np.diag(D_reml)).round(4))
# the closed form of prp-cls-growth(c) IS the REML maximizer: agreement to 4 decimals
assert abs(np.sqrt(sigma2_reml) - np.sqrt(sigma2)) < 5e-5
assert np.all(np.abs(np.sqrt(np.diag(D_reml)) - np.sqrt(np.diag(D_hat))) < 5e-5)
assert round(np.sqrt(sigma2_reml), 4) == round(np.sqrt(sigma2), 4)
assert np.all(np.round(np.sqrt(np.diag(D_reml)), 4) == np.round(np.sqrt(np.diag(D_hat)), 4))

gen = Generated("ch33", "growth_curves")
gen.int("N", N)
gen.int("m", m)
gen.num("b0", beta_hat[0], 4)
gen.num("b1", beta_hat[1], 4)
gen.num("se_b0", se_beta[0], 4)
gen.num("se_b1", se_beta[1], 4)
gen.num("t_b1", beta_hat[1] / se_beta[1], 2)
gen.num("sigma", np.sqrt(sigma2), 4)
gen.num("d11", D_hat[0, 0], 4)
gen.num("d12", D_hat[0, 1], 4)
gen.num("d22", D_hat[1, 1], 4)
gen.num("sd_int", np.sqrt(D_hat[0, 0]), 4)
gen.num("sd_slope", np.sqrt(D_hat[1, 1]), 4)
gen.num("corr_int_slope", D_hat[0, 1] / np.sqrt(D_hat[0, 0] * D_hat[1, 1]), 3)
gen.num("w11", W[0, 0], 4)
gen.num("w22", W[1, 1], 4)
gen.num("w12", W[0, 1], 4)
gen.num("w21", W[1, 0], 4)
gen.num("shrink_slope", shrink_slope, 3)
gen.num("slope_min", B[1].min(), 3)
gen.num("slope_max", B[1].max(), 3)
gen.num("blup_min", B_blup[1].min(), 3)
gen.num("blup_max", B_blup[1].max(), 3)
gen.num("se_curve_known", se_curve_known, 4)
gen.num("se_curve_new", se_curve_new, 4)
gen.num("se_blup_known", float(np.std(errs_known)), 4)
gen.num("se_blup_sim_new", float(np.std(errs_new)), 4)
gen.num("pred_new", pred_new, 4)
gen.num("pred_lo", pred_new - stats.norm.ppf(0.975) * se_curve_new, 3)
gen.num("pred_hi", pred_new + stats.norm.ppf(0.975) * se_curve_new, 3)
gen.num("sigma_reml", np.sqrt(sigma2_reml), 4)
gen.num("sd_int_reml", np.sqrt(D_reml[0, 0]), 4)
gen.num("sd_slope_reml", np.sqrt(D_reml[1, 1]), 4)
j_small = int(np.argmax(np.abs(B[1] - B_blup[1])))
gen.text("firm_shrunk", firms[j_small])
gen.num("slope_shrunk_own", B[1, j_small], 3)
gen.num("slope_shrunk_blup", B_blup[1, j_small], 3)
gen.write()

# ---- figure ------------------------------------------------------------------
use_book_style()
fig, axes = plt.subplots(3, 4, figsize=(5.8, 4.1))
for j, f in enumerate(firms):
    ax = axes.ravel()[j]
    ax.plot(years, Y[:, j], "o", color=COLORS["muted"], markersize=1.8)
    ax.plot(years, X @ B[:, j], ":", color=COLORS["second"], linewidth=1.0)
    ax.plot(years, X @ B_blup[:, j], color=COLORS["accent"], linewidth=1.1)
    ax.plot(years, X @ beta_hat, "--", color=COLORS["third"], linewidth=0.8)
    ax.set_title(f, fontsize=6.5, pad=2)
    ax.tick_params(labelsize=6)
    ax.set_xticks([1935, 1945, 1955])
    ax.set_xlim(1933, 1956)
ax = axes.ravel()[-1]
ax.plot(B[1], B_blup[1], "o", color=COLORS["accent"], markersize=2.4)
lims = [B[1].min() - 0.06, B[1].max() + 0.06]
ax.plot(lims, lims, ":", color=COLORS["muted"], linewidth=0.8)
ax.axhline(beta_hat[1], color=COLORS["third"], linestyle="--", linewidth=0.8)
ax.set_xlim(lims)
ax.set_ylim(lims)
ax.set_title("slope: own vs predicted", fontsize=6.5, pad=2)
ax.tick_params(labelsize=6)
ax.set_xticks([0.2, 0.5, 0.8])
handles = [plt.Line2D([], [], marker="o", linestyle="none", color=COLORS["muted"],
                      markersize=2.4, label="data"),
           plt.Line2D([], [], linestyle=":", color=COLORS["second"], label="firm's own line"),
           plt.Line2D([], [], color=COLORS["accent"], label="predicted curve"),
           plt.Line2D([], [], linestyle="--", color=COLORS["third"], label="population curve")]
fig.supylabel("log investment", fontsize=8)
fig.tight_layout(rect=(0, 0.06, 1, 1))
fig.legend(handles=handles, fontsize=7, frameon=False, ncol=4,
           loc="lower center", bbox_to_anchor=(0.5, 0.0))
fig.savefig(figure_path("ch33", "growth_curves"))
