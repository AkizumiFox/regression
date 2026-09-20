"""Chapter 32, Section 3: BLUP, Henderson's equations and shrinkage on the proficiency study.

Three things are checked here. (i) Henderson's mixed model equations, solved as a
(p + q) x (p + q) system, return the GLS estimate and the BLUP computed from V^{-1}.
(ii) The closed form of the prediction error variance in the balanced one-way model
matches a Monte Carlo estimate. (iii) The BLUP beats the unshrunk laboratory deviation
in mean squared error, by the margin the algebra predicts.
"""
import matplotlib.pyplot as plt
import numpy as np

from regbook import COLORS, Generated, figure_path, use_book_style

# <<setup>>
import numpy as np

g, m = 12, 4
n = g * m
mu, sigma_a, sigma = 2.50, 0.20, 0.25           # the values used to simulate

rng = np.random.default_rng(20320048)
lab_effect = sigma_a * rng.standard_normal(g)
y = mu + np.repeat(lab_effect, m) + sigma * rng.standard_normal(n)

X = np.ones((n, 1))                             # one fixed effect, the grand mean
Z = np.repeat(np.eye(g), m, axis=0)             # laboratory indicators

Y = y.reshape(g, m)                             # moment estimates, Section 32.1
lab_mean = Y.mean(axis=1)
ms_between = m * np.sum((lab_mean - Y.mean()) ** 2) / (g - 1)
ms_within = np.sum((Y - lab_mean[:, None]) ** 2) / (n - g)
s2, s2a = ms_within, (ms_between - ms_within) / m

G = s2a * np.eye(g)
R = s2 * np.eye(n)
V = Z @ G @ Z.T + R
# <</setup>>

# <<blup>>
Vinv = np.linalg.inv(V)
beta_gls = np.linalg.solve(X.T @ Vinv @ X, X.T @ Vinv @ y)
u_blup = G @ Z.T @ Vinv @ (y - X @ beta_gls)

# Henderson's mixed model equations: one (p + q) x (p + q) system, no n x n inverse
Rinv = np.linalg.inv(R)
C = np.block([[X.T @ Rinv @ X, X.T @ Rinv @ Z],
              [Z.T @ Rinv @ X, Z.T @ Rinv @ Z + np.linalg.inv(G)]])
rhs = np.concatenate([X.T @ Rinv @ y, Z.T @ Rinv @ y])
solution = np.linalg.solve(C, rhs)

print("GLS estimate   ", beta_gls, solution[:1])
print("first two BLUPs", u_blup[:2], solution[1:3])
# <</blup>>

assert np.allclose(solution[:1], beta_gls)
assert np.allclose(solution[1:], u_blup)

# <<shrink>>
shrinkage = m * s2a / (s2 + m * s2a)            # the shrinkage factor B
print(f"B = {shrinkage:.4f}")
print("u_hat  ", np.round(u_blup[:4], 4))
print("B(ybar_k - ybar)", np.round(shrinkage * (lab_mean - y.mean())[:4], 4))
# <</shrink>>

assert np.allclose(u_blup, shrinkage * (lab_mean - y.mean()))

# ---- prediction error variance, closed form against Monte Carlo --------------------
lam1 = s2 + m * s2a
pev_formula = s2a / lam1 * (s2 + s2a * m / g)
pev_matrix = np.diag(G - G @ Z.T @ Vinv @ Z @ G)[0] + \
    (X.T @ Vinv @ Z @ G[:, 0]) ** 2 / (X.T @ Vinv @ X)[0, 0]
assert np.isclose(pev_formula, pev_matrix.item())

# the raw deviation ybar_k - ybar as a predictor of a_k
raw_var = s2 / m * (1 - 1 / g) + s2a / g

B = 60000
sim = np.random.default_rng(1109)
a_sim = np.sqrt(s2a) * sim.standard_normal((B, g))
Ysim = mu + a_sim[:, :, None] + np.sqrt(s2) * sim.standard_normal((B, g, m))
lm = Ysim.mean(axis=2)
dev = lm - lm.mean(axis=1, keepdims=True)
mse_blup = np.mean((shrinkage * dev - a_sim) ** 2)
mse_raw = np.mean((dev - a_sim) ** 2)
assert abs(mse_blup / pev_formula - 1) < 0.02, (mse_blup, pev_formula)
assert abs(mse_raw / raw_var - 1) < 0.02, (mse_raw, raw_var)
assert pev_formula < raw_var

gen = Generated("ch32", "blup")
gen.num("shrinkage", shrinkage, 3)
gen.num("pev", pev_formula, 5)
gen.num("pevsd", np.sqrt(pev_formula), 4)
gen.num("rawvar", raw_var, 5)
gen.num("rawsd", np.sqrt(raw_var), 4)
gen.num("ratio", raw_var / pev_formula, 2)
gen.num("umax", u_blup.max(), 3)
gen.num("umin", u_blup.min(), 3)
gen.num("devmax", (lab_mean - y.mean()).max(), 3)
gen.num("devmin", (lab_mean - y.mean()).min(), 3)
gen.num("beta", beta_gls[0], 4)
gen.write()

# ---- figure -----------------------------------------------------------------------
use_book_style()
fig, axes = plt.subplots(1, 2, figsize=(5.8, 2.6))
ax = axes[0]
dev_obs = lab_mean - y.mean()
order = np.argsort(dev_obs)
half = 1.96 * np.sqrt(pev_formula)
for j, k in enumerate(order):
    ax.plot([dev_obs[k], u_blup[k]], [j, j], color=COLORS["grid"], lw=0.8, zorder=0)
    ax.plot([u_blup[k] - half, u_blup[k] + half], [j, j], color=COLORS["accent"], lw=0.8)
ax.plot(dev_obs[order], range(g), "o", ms=4, color=COLORS["second"], mew=0, label="$\\bar y_k-\\bar y$")
ax.plot(u_blup[order], range(g), "o", ms=4, color=COLORS["accent"], mew=0, label="$\\hat a_k$")
ax.axvline(0, color=COLORS["ink"], lw=0.6, ls="--")
ax.set_yticks([])
ax.set_ylabel("laboratory, sorted")
ax.set_xlabel("deviation from the overall mean")
ax.legend(frameon=False, loc="lower right", fontsize=7)
ax.set_title("(a) shrinkage of the 12 deviations")

ax = axes[1]
iccs = np.linspace(0.01, 0.9, 200)
Bfac = m * iccs / (m * iccs + (1 - iccs))
pev = (iccs / (m * iccs + 1 - iccs)) * ((1 - iccs) + iccs * m / g)
raw = (1 - iccs) / m * (1 - 1 / g) + iccs / g
ax.plot(iccs, pev, color=COLORS["accent"], label="BLUP")
ax.plot(iccs, raw, "--", color=COLORS["second"], label=r"$\bar y_k-\bar y$")
ax.set_xlabel("intraclass correlation")
ax.set_ylabel("mean squared error")
ax.legend(frameon=False)
ax.set_title("(b) predicting $a_k$, with $\\sigma_a^2+\\sigma^2=1$")
fig.tight_layout()
fig.savefig(figure_path("ch32", "blup_shrinkage"))
