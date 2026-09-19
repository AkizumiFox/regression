"""Chapter 20, Section 1: what residual plots show.

(1) Four simulated data sets fitted by a straight line: correct model, a curved mean, a
    variance that grows with the mean, and one gross error. Residuals against fitted values.
(2) A curved effect of one regressor that is correlated with another: the added-variable
    plot and the component-plus-residual plot have the same slope and the same residuals,
    but only the second shows the curve on the regressor's own scale.
"""
import matplotlib.pyplot as plt
import numpy as np

from regbook import COLORS, Generated, figure_path, use_book_style


def ls_fit(X, y):
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    return beta, X @ beta, y - X @ beta


# <<patterns>>
rng = np.random.default_rng(7)
n = 80
x = rng.uniform(0, 10, n)
eps = rng.normal(size=n)
X = np.column_stack([np.ones(n), x])
responses = {
    "(a) model correct": 1 + 0.5 * x + eps,
    "(b) curved mean": 1 + 0.5 * x + 0.08 * (x - 5) ** 2 + 0.6 * eps,
    "(c) variance grows": 1 + 0.5 * x + (0.1 + 0.25 * x) * eps,
    "(d) one gross error": 1 + 0.5 * x + eps + 8.0 * (np.arange(n) == 17),
}
for name, y in responses.items():
    beta, fitted, resid = ls_fit(X, y)
    print(f"{name:22s} slope {beta[1]:.3f}   residuals . fitted = {resid @ fitted:.1e}")
# <</patterns>>

stats_rec = {}
for name, y in responses.items():
    beta, fitted, resid = ls_fit(X, y)
    assert abs(resid @ fitted) < 1e-8 and abs(resid @ x) < 1e-8 and abs(resid.sum()) < 1e-9
    # regressing the residuals on y gives slope 1 - R^2
    yc = y - y.mean()
    r2 = 1 - resid @ resid / (yc @ yc)
    slope_on_y = (resid @ yc) / (yc @ yc)
    assert np.isclose(slope_on_y, 1 - r2)
    stats_rec[name] = (r2, slope_on_y)

# residuals estimate (I - M) delta when the mean is misspecified
delta = 0.08 * (x - 5) ** 2
Q, _ = np.linalg.qr(X)
proj_delta = delta - Q @ (Q.T @ delta)
reps = 4000
acc = np.zeros(n)
for _ in range(reps):
    yy = 1 + 0.5 * x + delta + 0.6 * rng.normal(size=n)
    acc += yy - Q @ (Q.T @ yy)
assert np.max(np.abs(acc / reps - proj_delta)) < 0.05

# <<partial>>
rng = np.random.default_rng(11)
n2 = 120
x1 = rng.uniform(0, 4, n2)
x2 = 0.8 * x1 + rng.normal(0, 0.6, n2)             # correlated with x1
y2 = 1 + x2 + 0.7 * (x1 - 2) ** 2 + rng.normal(0, 0.4, n2)
X2 = np.column_stack([np.ones(n2), x1, x2])        # the fitted model is linear in x1
beta2, fitted2, resid2 = ls_fit(X2, y2)

component = resid2 + beta2[1] * x1                 # component-plus-residual for x1
Z = np.column_stack([np.ones(n2), x2])             # everything except x1
x1_tilde = x1 - Z @ np.linalg.lstsq(Z, x1, rcond=None)[0]
y_tilde = y2 - Z @ np.linalg.lstsq(Z, y2, rcond=None)[0]   # added-variable coordinates
print("coefficient of x1:", beta2[1])
print("slope in component-plus-residual plot:", ls_fit(np.column_stack([np.ones(n2), x1]), component)[0][1])
print("slope in added-variable plot:", (x1_tilde @ y_tilde) / (x1_tilde @ x1_tilde))
# <</partial>>

b_cr, _, res_cr = ls_fit(np.column_stack([np.ones(n2), x1]), component)
b_av = (x1_tilde @ y_tilde) / (x1_tilde @ x1_tilde)
assert np.isclose(b_cr[1], beta2[1]) and np.isclose(b_av, beta2[1])
assert np.isclose(b_cr[0], 0, atol=1e-10)
assert np.allclose(res_cr, resid2) and np.allclose(y_tilde - b_av * x1_tilde, resid2)
corr12 = np.corrcoef(x1, x2)[0, 1]
# curvature visible in the component-plus-residual plot: quadratic term in each plot
q_cr = ls_fit(np.column_stack([np.ones(n2), x1, x1 ** 2]), component)[0][2]
xa = x1_tilde
q_av = ls_fit(np.column_stack([np.ones(n2), xa, xa ** 2]), y_tilde)[0][2]
assert abs(q_cr - 0.7) < 0.1

gen = Generated("ch20", "residual_plots", prefix="rp")
gen.int("n", n)
gen.int("n2", n2)
gen.num("corr12", corr12, 2)
gen.num("beta1", beta2[1], 3)
gen.num("q_cr", q_cr, 3)
gen.num("q_av", q_av, 3)
gen.write()

# ---- figure 1: four residual plots ----------------------------------------------------
use_book_style()
fig, axes = plt.subplots(1, 4, figsize=(7.0, 1.9), sharey=False)
for ax, (name, y) in zip(axes, responses.items()):
    _, fitted, resid = ls_fit(X, y)
    ax.scatter(fitted, resid, s=6, color=COLORS["accent"], alpha=0.8, linewidths=0)
    ax.axhline(0, color=COLORS["grid"], linewidth=0.8, zorder=0)
    ax.set_title(name, fontsize=8)
    ax.set_xlabel("fitted value", fontsize=8)
    ax.tick_params(labelsize=7)
axes[0].set_ylabel("residual")
fig.tight_layout(w_pad=0.6)
fig.savefig(figure_path("ch20", "residual_patterns"))

# ---- figure 2: added-variable versus component-plus-residual -------------------------
fig, axes = plt.subplots(1, 2, figsize=(5.6, 2.4))
ax = axes[0]
ax.scatter(x1_tilde, y_tilde, s=7, color=COLORS["accent"], alpha=0.8, linewidths=0)
g = np.linspace(x1_tilde.min(), x1_tilde.max(), 2)
ax.plot(g, b_av * g, color=COLORS["second"])
ax.set_xlabel(r"$x_1$ residualized on $1, x_2$")
ax.set_ylabel(r"$y$ residualized on $1, x_2$")
ax.set_title("(a) added-variable plot")
ax = axes[1]
ax.scatter(x1, component, s=7, color=COLORS["accent"], alpha=0.8, linewidths=0)
g = np.linspace(0, 4, 2)
ax.plot(g, beta2[1] * g, color=COLORS["second"])
ax.set_xlabel(r"$x_1$")
ax.set_ylabel(r"$\hat{\varepsilon}_i+\hat\beta_1x_{i1}$")
ax.set_title("(b) component-plus-residual plot")
fig.tight_layout()
fig.savefig(figure_path("ch20", "partial_residual"))
