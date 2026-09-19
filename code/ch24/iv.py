"""Chapter 24, Section 4: instrumental variables.

(a) Blood-pressure study: the second reading instruments the first (age instruments itself).
    IV estimate, homoscedastic and sandwich standard errors, and the first-stage F statistic.
(b) Weak instruments: y = x + e, x = pi z + v, corr(e, v) = 0.9, n = 200, z fixed with
    sum z^2 = n. For concentration parameters mu^2 = n pi^2 in {1, 10, 100}, the distribution
    of the IV estimate, the coverage of the nominal 95% Wald interval, and a check of the exact
    representation (rho eta + sqrt(1 - rho^2) xi) / (mu + eta).
"""
import matplotlib.pyplot as plt
from scipy import stats

from regbook import COLORS, Generated, figure_path, use_book_style

gen = Generated("ch24", "iv", prefix="iv")

# <<setup>>
import numpy as np

rng = np.random.default_rng(2401)
n = 300
age = rng.normal(55, 10, n)
x = 130 + 0.6 * (age - 55) + rng.normal(0, np.sqrt(108), n)   # long-run blood pressure, sd 12
w = x[:, None] + rng.normal(0, 9, (n, 2))                      # two clinic readings
y = -4 + 0.08 * x + 0.0 * age + rng.normal(0, 1.0, n)          # age has no effect given x
# <</setup>>

# <<replicate-iv>>
X = np.column_stack([np.ones(n), w[:, 0], age])    # regressors: first reading, age
Z = np.column_stack([np.ones(n), w[:, 1], age])    # instruments: second reading, age
b_iv = np.linalg.solve(Z.T @ X, Z.T @ y)
e = y - X @ b_iv                                    # residuals use X, not the first-stage fit
ZXinv = np.linalg.inv(Z.T @ X)
cov_h = e @ e / n * ZXinv @ (Z.T @ Z) @ ZXinv.T               # homoscedastic
cov_s = ZXinv @ ((Z * e[:, None] ** 2).T @ Z) @ ZXinv.T       # heteroscedasticity-robust
print("IV (1, x, age):", b_iv.round(4))
print("se, homoscedastic:", np.sqrt(np.diag(cov_h)).round(4))
print("se, sandwich:     ", np.sqrt(np.diag(cov_s)).round(4))

# first stage: regress the first reading on the instruments; F for the excluded instrument
g, rss1 = np.linalg.lstsq(Z, w[:, 0], rcond=None)[:2]
rss0 = np.sum((w[:, 0] - np.column_stack([np.ones(n), age]) @
               np.linalg.lstsq(np.column_stack([np.ones(n), age]), w[:, 0], rcond=None)[0]) ** 2)
F1 = (rss0 - rss1[0]) / (rss1[0] / (n - 3))
print(f"first-stage F = {F1:.1f}")
# <</replicate-iv>>

se_h = np.sqrt(np.diag(cov_h))
se_s = np.sqrt(np.diag(cov_s))
# 2SLS through the projection onto C(Z) gives the same estimate (just identified)
Xhat = Z @ np.linalg.lstsq(Z, X, rcond=None)[0]
b_2sls = np.linalg.lstsq(Xhat, y, rcond=None)[0]
assert np.allclose(b_2sls, b_iv)
# the naive second-stage residuals are wrong
e_wrong = y - Xhat @ b_2sls
assert abs(e_wrong @ e_wrong - e @ e) > 1.0
assert abs(b_iv[1] - 0.08) < 2.5 * se_h[1]
assert abs(b_iv[2]) < 2.5 * se_h[2]
assert F1 > 100
gen.num("iv_x", b_iv[1], 4)
gen.num("iv_x_se", se_h[1], 4)
gen.num("iv_x_se_s", se_s[1], 4)
gen.num("iv_age", b_iv[2], 4)
gen.num("iv_age_se", se_h[2], 4)
gen.num("F1", F1, 1)
gen.num("s2_right", e @ e / n, 3)
gen.num("s2_wrong", e_wrong @ e_wrong / n, 3)

# ---- (b) weak instruments --------------------------------------------------------------------
# <<weak>>
def iv_draws(mu2, reps, seed, n=200, rho=0.9):
    """IV estimates of beta = 1 and Wald-interval coverage when the concentration is mu2."""
    r = np.random.default_rng(seed)
    z = r.normal(size=n)
    z = z * np.sqrt(n / (z @ z))                  # fixed instrument with sum z^2 = n
    pi = np.sqrt(mu2 / n)
    ev = r.multivariate_normal([0, 0], [[1, rho], [rho, 1]], size=(reps, n))
    xx = pi * z + ev[:, :, 1]
    yy = 1.0 * xx + ev[:, :, 0]
    b = (yy @ z) / (xx @ z)
    res = yy - b[:, None] * xx
    se = np.sqrt((res ** 2).mean(axis=1) * n / (xx @ z) ** 2)
    return b, np.abs(b - 1) <= 1.96 * se


for mu2 in [1, 10, 100]:
    b, hit = iv_draws(mu2, 2000, 2440)
    print(f"mu^2 = {mu2:3d}: median {np.median(b):.3f}, "
          f"quartiles {np.percentile(b, 25):.3f} to {np.percentile(b, 75):.3f}, coverage {hit.mean():.3f}")
# <</weak>>

reps = 20_000
res = {}
for mu2 in [1, 10, 100]:
    b, hit = iv_draws(mu2, reps, 2441)
    # exact representation: (b - 1) = (rho eta + sqrt(1 - rho^2) xi) / (mu + eta)
    r2 = np.random.default_rng(99)
    eta, xi = r2.normal(size=reps), r2.normal(size=reps)
    rep = (0.9 * eta + np.sqrt(1 - 0.81) * xi) / (np.sqrt(mu2) + eta)
    assert stats.ks_2samp(b - 1, rep).pvalue > 1e-3
    res[mu2] = b
    gen.num(f"med_{mu2}", np.median(b), 3)
    gen.num(f"q1_{mu2}", np.percentile(b, 25), 3)
    gen.num(f"q3_{mu2}", np.percentile(b, 75), 3)
    gen.num(f"cover_{mu2}", hit.mean(), 3)
gen.num("outside_1", np.mean((res[1] <= -1) | (res[1] >= 3)), 3)
ols_limit = 1 + 0.9                                # plim of least squares as pi -> 0
assert np.median(res[1]) - 1 > 0.25
assert float(gen.items['cover_1']) < float(gen.items['cover_10']) < float(gen.items['cover_100'])
assert abs(np.median(res[100]) - 1) < 0.02
gen.int("reps", reps)
gen.write()

# ---- figure ----------------------------------------------------------------------------------
use_book_style()
fig, ax = plt.subplots(figsize=(4.6, 2.5))
bins = np.linspace(-1, 3, 161)
for mu2, col in [(1, COLORS["second"]), (10, COLORS["thread"]), (100, COLORS["accent"])]:
    inside = res[mu2][(res[mu2] > -1) & (res[mu2] < 3)]         # area = fraction inside the window
    ax.hist(inside, bins=bins, weights=np.full(len(inside), 1 / (reps * (bins[1] - bins[0]))),
            histtype="step", color=col, linewidth=1.1, label=rf"$\mu^2={mu2}$")
ax.axvline(1.0, color=COLORS["ink"], linestyle="--", linewidth=0.9)
ax.axvline(ols_limit, color=COLORS["muted"], linestyle=":", linewidth=1.0)
ax.text(1.0, 1.01, r"$\beta$", fontsize=8, ha="center", va="bottom",
        transform=ax.get_xaxis_transform())          # above the axes, clear of the peak
ax.text(ols_limit + 0.03, ax.get_ylim()[1] * 0.92, "least squares limit", fontsize=7, color=COLORS["muted"])
ax.set_xlim(-1, 3)
ax.set_xlabel("IV estimate")
ax.set_ylabel("density")
ax.legend(frameon=False, loc="upper left", fontsize=7)
fig.tight_layout()
fig.savefig(figure_path("ch24", "iv_weak"))
