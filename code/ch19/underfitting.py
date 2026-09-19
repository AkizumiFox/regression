"""Chapter 19, Section 1: omitting a real regressor versus keeping an irrelevant one.

(a) Two regressors with sample correlation r. The short regression (x2 omitted) has mean
    squared error for beta_1 equal to Var(long) * (1 - r^2 + r^2 gamma), where gamma is the
    noncentrality of the F test of beta_2 = 0. The short estimator wins iff gamma < 1.
(b) Polynomial fits to a smooth mean: the risk of the fitted values at the design points is
    sigma^2 p + ||(I - M_p) mu||^2, a bias-variance trade-off with an interior minimum.
"""
import matplotlib.pyplot as plt
import numpy as np

from regbook import COLORS, Generated, figure_path, use_book_style

gen = Generated("ch19", "underfitting", prefix="uf")

# <<design>>
import numpy as np
rng = np.random.default_rng(1901)
n, r, sigma = 30, 0.9, 1.0

# two centred, unit-variance regressors with sample correlation exactly r:
# orthonormalize (1, z1, z2) and keep the last two columns, which are centred
Z = np.linalg.qr(np.column_stack([np.ones(n), rng.normal(size=(n, 2))]))[0][:, 1:]
Z *= np.sqrt(n - 1)
x1 = Z[:, 0]
x2 = r * Z[:, 0] + np.sqrt(1 - r**2) * Z[:, 1]
X1 = np.column_stack([np.ones(n), x1])                # the short model
X = np.column_stack([X1, x2])                         # the long model

M1 = X1 @ np.linalg.solve(X1.T @ X1, X1.T)
S22_1 = x2 @ (np.eye(n) - M1) @ x2                    # ||(I - M1) x2||^2
beta2 = np.sqrt(0.5 * sigma**2 / S22_1)               # chosen so that gamma = 0.5
gamma = beta2**2 * S22_1 / sigma**2

var_long = sigma**2 * np.linalg.inv(X.T @ X)[1, 1]    # Var of the long slope on x1
Pi = np.linalg.solve(X1.T @ X1, X1.T @ x2)[1]         # slope of x2 on x1
var_short = sigma**2 * np.linalg.inv(X1.T @ X1)[1, 1]
mse_short = var_short + (Pi * beta2) ** 2             # variance + squared bias
print(f"gamma = {gamma:.3f}, MSE(short)/Var(long) = {mse_short / var_long:.4f}")
print(f"formula 1 - r^2 + r^2 gamma      = {1 - r**2 + r**2 * gamma:.4f}")
# <</design>>

assert np.isclose(np.corrcoef(x1, x2)[0, 1], r)
assert np.isclose(Pi, r)
assert np.isclose(mse_short / var_long, 1 - r**2 + r**2 * gamma)

# <<simulate>>
reps = 200_000
beta = np.array([1.0, 1.0, beta2])
Y = X @ beta + sigma * rng.normal(size=(reps, n))
b_long = np.linalg.solve(X.T @ X, X.T @ Y.T)[1]
b_short = np.linalg.solve(X1.T @ X1, X1.T @ Y.T)[1]
print("simulated MSE ratio:", np.mean((b_short - 1) ** 2) / np.mean((b_long - 1) ** 2))
# <</simulate>>
mse_ratio_sim = np.mean((b_short - 1) ** 2) / np.mean((b_long - 1) ** 2)
assert abs(mse_ratio_sim / (mse_short / var_long) - 1) < 0.02
# short-model s^2 is biased upwards by gamma sigma^2/(n - p1)
s2_short = np.sum((Y - (M1 @ Y.T).T) ** 2, axis=1) / (n - 2)
Es2 = sigma**2 * (1 + gamma / (n - 2))
assert abs(s2_short.mean() / Es2 - 1) < 0.005

gen.int("n", n)
gen.num("r", r, 1)
gen.num("gamma", gamma, 1)
gen.num("beta2", beta2, 3)
gen.num("ratio", mse_short / var_long, 3)
gen.num("ratio_sim", mse_ratio_sim, 3)
gen.num("bias", Pi * beta2, 3)
gen.num("sd_long", np.sqrt(var_long), 3)
gen.num("sd_short", np.sqrt(var_short), 3)
gen.num("Es2", Es2, 4)
gen.int("reps", reps)

# ---- (b) polynomial degree: bias^2 + variance -------------------------------------------
xg = np.linspace(-1, 1, 40)
mu = np.exp(xg) * np.sin(3 * xg)                       # a smooth mean, not a polynomial
sig_b = 0.25
degrees = np.arange(0, 11)
risk, bias2 = [], []
for d in degrees:
    Q = np.linalg.qr(np.vander(xg, d + 1, increasing=True))[0]
    resid = mu - Q @ (Q.T @ mu)
    bias2.append(resid @ resid)
    risk.append(sig_b**2 * (d + 1) + resid @ resid)
risk, bias2 = np.array(risk), np.array(bias2)
d_best = int(degrees[np.argmin(risk)])
assert np.all(np.diff(bias2) <= 1e-12)                  # bias never increases with the degree
assert 2 <= d_best <= 8
gen.int("n_poly", len(xg))
gen.num("sigma_poly", sig_b, 2)
gen.int("d_best", d_best)
gen.num("risk_best", risk[d_best] / sig_b**2, 2)
gen.num("risk_max", risk[-1] / sig_b**2, 2)
gen.num("risk_line", risk[1] / sig_b**2, 2)
gen.write()

# ---- figure -----------------------------------------------------------------------------
use_book_style()
fig, axes = plt.subplots(1, 2, figsize=(5.8, 2.5))
ax = axes[0]
g = np.linspace(0, 3, 50)
for rr, col in [(0.3, COLORS["third"]), (0.6, COLORS["accent"]), (0.9, COLORS["second"])]:
    ax.plot(g, 1 - rr**2 + rr**2 * g, color=col, label=f"$r={rr}$")
ax.plot([gamma], [mse_ratio_sim], "o", color=COLORS["second"], ms=4)
ax.axhline(1, color=COLORS["grid"], lw=0.8, zorder=0)
ax.axvline(1, color=COLORS["muted"], lw=0.8, ls=":")
ax.set_xlabel(r"noncentrality $\gamma$")
ax.set_ylabel(r"MSE(short) / Var(long)")
ax.set_title(r"(a) omitting $x_2$: estimating $\beta_1$")
ax.legend(frameon=False, loc="upper left")
ax = axes[1]
ax.plot(degrees, bias2 / sig_b**2, "s-", ms=3, color=COLORS["second"], label="squared bias")
ax.plot(degrees, degrees + 1, "^-", ms=3, color=COLORS["third"], label="variance")
ax.plot(degrees, risk / sig_b**2, "o-", ms=3, color=COLORS["accent"], label="risk")
ax.set_yscale("log")
ax.set_ylim(0.3, 1500)
ax.set_xlabel("polynomial degree")
ax.set_ylabel(r"risk / $\sigma^2$")
ax.set_title("(b) fitted values at the design")
ax.legend(frameon=False)
fig.tight_layout()
fig.savefig(figure_path("ch19", "underfitting"))
