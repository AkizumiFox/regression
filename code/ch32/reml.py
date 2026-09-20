"""Chapter 32, Section 4: maximum likelihood and REML for the proficiency study.

The two log-likelihoods are coded from their general matrix definitions and checked
against (i) the closed forms of the balanced one-way model, (ii) statsmodels' MixedLM.
A simulation measures the downward bias of the maximum likelihood estimate of the
between-laboratory variance and compares it with the exact value lambda_1 / (g m).
"""
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy import optimize

from regbook import COLORS, Generated, figure_path, use_book_style

# <<likelihood>>
import numpy as np

g, m = 12, 4
n = g * m
mu, sigma_a, sigma = 2.50, 0.20, 0.25           # the values used to simulate

rng = np.random.default_rng(20320048)
lab_effect = sigma_a * rng.standard_normal(g)
y = mu + np.repeat(lab_effect, m) + sigma * rng.standard_normal(n)

X = np.ones((n, 1))
Z = np.repeat(np.eye(g), m, axis=0)


def loglik(s2, s2a, reml=True):
    """Log-likelihood (reml=False) or restricted log-likelihood of the variances."""
    V = s2a * (Z @ Z.T) + s2 * np.eye(n)
    L = np.linalg.cholesky(V)
    Vi_X = np.linalg.solve(V, X)
    W = X.T @ Vi_X                              # X' V^{-1} X
    r = y - X @ np.linalg.solve(W, Vi_X.T @ y)  # GLS residual
    quad = r @ np.linalg.solve(V, r)            # = y' P y
    logdetV = 2 * np.sum(np.log(np.diag(L)))
    out = -0.5 * (n * np.log(2 * np.pi) + logdetV + quad)
    if reml:
        out += -0.5 * (np.linalg.slogdet(W)[1] - np.linalg.slogdet(X.T @ X)[1]
                       - X.shape[1] * np.log(2 * np.pi))
    return out
# <</likelihood>>

# <<closed>>
Y = y.reshape(g, m)
lab_mean = Y.mean(axis=1)
ms_between = m * np.sum((lab_mean - Y.mean()) ** 2) / (g - 1)
ms_within = np.sum((Y - lab_mean[:, None]) ** 2) / (n - g)

reml_s2, reml_s2a = ms_within, (ms_between - ms_within) / m
ml_s2 = ms_within
ml_s2a = ((g - 1) / g * ms_between - ms_within) / m

print(f"REML  sigma^2 {reml_s2:.4f}  sigma_a^2 {reml_s2a:.4f}")
print(f"ML    sigma^2 {ml_s2:.4f}  sigma_a^2 {ml_s2a:.4f}")
print("restricted log-likelihood", round(loglik(reml_s2, reml_s2a), 4))
print("log-likelihood           ", round(loglik(ml_s2, ml_s2a, reml=False), 4))
# <</closed>>

# the closed forms really are the maximizers of the coded likelihoods
for reml, pair in [(True, (reml_s2, reml_s2a)), (False, (ml_s2, ml_s2a))]:
    best = optimize.minimize(lambda t: -loglik(np.exp(t[0]), np.exp(t[1]), reml=reml),
                             np.log(pair), method="Nelder-Mead",
                             options={"xatol": 1e-10, "fatol": 1e-12, "maxiter": 5000})
    assert np.allclose(np.exp(best.x), pair, rtol=1e-4), (reml, np.exp(best.x), pair)

# and they agree with statsmodels
frame = pd.DataFrame({"y": y, "lab": np.repeat(np.arange(g), m)})
model = sm.MixedLM.from_formula("y ~ 1", groups="lab", data=frame)
sm_reml, sm_ml = model.fit(reml=True), model.fit(reml=False)
assert abs(sm_reml.cov_re.iloc[0, 0] - reml_s2a) < 1e-5
assert abs(sm_reml.scale - reml_s2) < 1e-5
assert abs(sm_ml.cov_re.iloc[0, 0] - ml_s2a) < 1e-5
assert abs(sm_ml.scale - ml_s2) < 1e-5

# ---- the bias of maximum likelihood ------------------------------------------------
lam1 = sigma**2 + m * sigma_a**2
bias_exact = -lam1 / (g * m)
B = 40000
sim = np.random.default_rng(404)
a_sim = sigma_a * sim.standard_normal((B, g))
Ysim = mu + a_sim[:, :, None] + sigma * sim.standard_normal((B, g, m))
lm = Ysim.mean(axis=2)
msb = m * np.sum((lm - lm.mean(axis=1, keepdims=True)) ** 2, axis=1) / (g - 1)
msw = np.sum((Ysim - lm[:, :, None]) ** 2, axis=(1, 2)) / (n - g)
reml_sim = (msb - msw) / m
ml_sim = ((g - 1) / g * msb - msw) / m          # the interior (unconstrained) formula
ml_trunc = np.maximum(ml_sim, 0.0)              # the maximum likelihood estimator itself
assert abs(reml_sim.mean() - sigma_a**2) < 0.002
assert abs((ml_sim.mean() - sigma_a**2) - bias_exact) < 0.002
assert abs(np.mean(ml_sim < 0) - np.mean(reml_sim < 0)) < 0.08
# truncation at zero can only raise the estimate, so it shrinks the bias in modulus
assert 0 > ml_trunc.mean() - sigma_a**2 > ml_sim.mean() - sigma_a**2
# the ML formula is negative exactly when MSB < g MSE / (g - 1), a wider range than REML's
assert np.all((ml_sim < 0) == (msb < g * msw / (g - 1)))
assert np.all((reml_sim < 0) == (msb < msw))

# ---- EM reproduces the closed-form ML estimates ------------------------------------
# <<em>>
ZZt = Z @ Z.T
s2, s2a = 0.1, 0.1                              # a deliberately poor starting value
for _ in range(200):
    V = s2a * ZZt + s2 * np.eye(n)
    Vi = np.linalg.inv(V)
    W = X.T @ Vi @ X
    beta_hat = np.linalg.solve(W, X.T @ Vi @ y)
    u_hat = s2a * Z.T @ Vi @ (y - X @ beta_hat)
    S = np.linalg.inv(np.eye(g) / s2a + Z.T @ Z / s2)   # Cov(u | y), not (C^{-1})_22
    resid = y - X @ beta_hat - Z @ u_hat
    s2a = (u_hat @ u_hat + np.trace(S)) / g
    s2 = (resid @ resid + np.trace(Z @ S @ Z.T)) / n

print(f"EM    sigma^2 {s2:.4f}  sigma_a^2 {s2a:.4f}")
# <</em>>
assert abs(s2 - ml_s2) < 1e-9 and abs(s2a - ml_s2a) < 1e-9, (s2, s2a)

# using the trailing block of C^{-1} instead maximizes neither likelihood
h2, h2a = 0.1, 0.1
for _ in range(2000):
    Vh = h2a * ZZt + h2 * np.eye(n)
    Vih = np.linalg.inv(Vh)
    Wh = X.T @ Vih @ X
    bh = np.linalg.solve(Wh, X.T @ Vih @ y)
    uh = h2a * Z.T @ Vih @ (y - X @ bh)
    GZV = h2a * Z.T @ Vih
    Sh = (h2a * np.eye(g) - GZV @ Z * h2a
          + GZV @ X @ np.linalg.solve(Wh, X.T @ Vih @ Z * h2a))
    rh = y - X @ bh - Z @ uh
    h2a = (uh @ uh + np.trace(Sh)) / g
    h2 = (rh @ rh + np.trace(Z @ Sh @ Z.T)) / n
assert abs(h2 - ml_s2) > 1e-3 and abs(h2 - reml_s2) > 1e-3, (h2, h2a)
assert abs(h2a - ml_s2a) > 1e-3 and abs(h2a - reml_s2a) > 1e-4, (h2, h2a)

gen = Generated("ch32", "reml")
gen.num("remls2", reml_s2, 4)
gen.num("remls2a", reml_s2a, 4)
gen.num("mls2", ml_s2, 4)
gen.num("mls2a", ml_s2a, 4)
gen.num("ratio", ml_s2a / reml_s2a, 3)
gen.num("lr", loglik(reml_s2, reml_s2a), 3)
gen.num("ll", loglik(ml_s2, ml_s2a, reml=False), 3)
gen.num("biasexact", bias_exact, 4)
gen.num("biassim", ml_sim.mean() - sigma_a**2, 4)
gen.num("biastrunc", ml_trunc.mean() - sigma_a**2, 4)
gen.num("pmlneg", 100 * np.mean(ml_sim < 0), 1)
gen.num("hybs2", h2, 4)
gen.num("hybs2a", h2a, 4)
gen.num("relbias", -bias_exact / sigma_a**2, 3)
gen.num("smreml", sm_reml.cov_re.iloc[0, 0], 4)
gen.num("smml", sm_ml.cov_re.iloc[0, 0], 4)
gen.write()

# ---- figure: the two profile log-likelihoods of sigma_a^2 --------------------------
use_book_style()
grid = np.linspace(0.0, 0.16, 161)


def profile(s2a, reml):
    f = optimize.minimize_scalar(lambda t: -loglik(np.exp(t), s2a, reml=reml),
                                 bracket=(np.log(0.02), np.log(0.2)))
    return -f.fun


prof_r = np.array([profile(v, True) for v in grid])
prof_m = np.array([profile(v, False) for v in grid])
assert abs(grid[prof_r.argmax()] - reml_s2a) < 0.002
assert abs(grid[prof_m.argmax()] - ml_s2a) < 0.002

# the profile-likelihood interval for sigma_a^2: where REML drops by 1.92
def cross(lo, hi):
    return optimize.brentq(lambda v: profile(v, True) - prof_r.max() + 1.92, lo, hi)


prof_lo, prof_hi = cross(1e-6, reml_s2a), cross(reml_s2a, 0.5)
gen2 = Generated("ch32", "remlci")
gen2.num("lo", prof_lo, 4)
gen2.num("hi", prof_hi, 4)
gen2.write()

fig, ax = plt.subplots(figsize=(4.6, 2.7))
ax.plot(grid, prof_r - prof_r.max(), color=COLORS["accent"], label="REML")
ax.plot(grid, prof_m - prof_m.max(), "--", color=COLORS["second"], label="ML")
ax.plot([reml_s2a], [0], "o", ms=4, color=COLORS["accent"], mew=0)
ax.plot([ml_s2a], [0], "o", ms=4, color=COLORS["second"], mew=0)
ax.axhline(-1.92, color=COLORS["grid"], lw=0.7)
ax.set_ylim(-4.2, 0.35)
ax.set_xlabel(r"$\sigma_a^2$")
ax.set_ylabel("profile log-likelihood, relative")
ax.legend(frameon=False)
fig.tight_layout()
fig.savefig(figure_path("ch32", "reml_profile"))
