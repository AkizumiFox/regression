"""Chapter 34, Section 5: the chi-squared limits, and the three test statistics.

(i) Simulates the null distribution of the deviance for grouped and for ungrouped binary
data, to show where the chi-squared approximation holds and where it does not.
(ii) Draws the Hauck-Donner effect: the Wald statistic for a 2 x 2 table is not monotone
in the evidence, while the likelihood ratio statistic is.
(iii) Computes Wald, score and likelihood ratio statistics for one coefficient of a real fit.
"""
import matplotlib.pyplot as plt
import numpy as np
import statsmodels.api as sm
from scipy import stats

from regbook import COLORS, Generated, figure_path, use_book_style


def logit_irls(X, y, w, steps=25):
    """Maximum likelihood for a binomial model with logit link; y is a proportion."""
    mu = (w * y + 0.5) / (w + 1.0)
    eta = np.log(mu / (1 - mu))
    beta = np.zeros(X.shape[1])
    for _ in range(steps):
        d = mu * (1 - mu)
        W = w * d
        z = eta + (y - mu) / d
        beta = np.linalg.solve((X.T * W) @ X, (X.T * W) @ z)
        eta = X @ beta
        mu = 1 / (1 + np.exp(-eta))
    return beta, mu


def binomial_deviance(y, mu, w):
    t1 = np.where(y > 0, y * np.log(np.where(y > 0, y, 1.0) / mu), 0.0)
    t2 = np.where(y < 1, (1 - y) * np.log(np.where(y < 1, 1 - y, 1.0) / (1 - mu)), 0.0)
    return 2 * np.sum(w * (t1 + t2))


# ---- (i) the null distribution of the deviance -------------------------------
# <<deviance-null>>
def null_deviances(X, beta, m, reps, rng):
    """Residual deviance of the fitted model, on `reps` data sets simulated from it."""
    p = 1 / (1 + np.exp(-(X @ beta)))
    w = np.full(len(p), float(m))
    out = np.empty(reps)
    for r in range(reps):
        y = rng.binomial(m, p).astype(float) / m
        _, mu = logit_irls(X, y, w)
        out[r] = binomial_deviance(y, mu, w)
    return out

rng = np.random.default_rng(3407)
beta0 = np.array([0.3, 0.9])
G, m = 40, 20                                       # grouped: 40 groups of 20 trials
Xg = np.column_stack([np.ones(G), np.linspace(-1.5, 1.5, G)])
N = 800                                             # ungrouped: 800 single trials
Xu = np.column_stack([np.ones(N), np.linspace(-1.5, 1.5, N)])

reps = 2000
Dg = null_deviances(Xg, beta0, m, reps, rng)
Du = null_deviances(Xu, beta0, 1, reps, rng)

print(f"grouped   : mean {Dg.mean():8.2f} (df {G - 2}), sd {Dg.std():6.2f}"
      f" (sqrt(2 df) {np.sqrt(2 * (G - 2)):.2f})")
print(f"ungrouped : mean {Du.mean():8.2f} (df {N - 2}), sd {Du.std():6.2f}"
      f" (sqrt(2 df) {np.sqrt(2 * (N - 2)):.2f})")
# <</deviance-null>>

assert abs(Dg.mean() - (G - 2)) < 0.15 * (G - 2)          # grouped: chi-squared is close
assert abs(Du.mean() - (N - 2)) > 0.2 * (N - 2)           # ungrouped: it is not
assert Du.std() < 0.6 * np.sqrt(2 * (N - 2))

# The two formats are different designs, but they carry very nearly the same Fisher
# information: the standard errors of beta-hat agree to within 2 per cent.
p_g = 1 / (1 + np.exp(-(Xg @ beta0)))
p_u = 1 / (1 + np.exp(-(Xu @ beta0)))
se_grouped = np.sqrt(np.diag(np.linalg.inv((Xg.T * (m * p_g * (1 - p_g))) @ Xg)))
se_ungrouped = np.sqrt(np.diag(np.linalg.inv((Xu.T * (p_u * (1 - p_u))) @ Xu)))
print("standard errors: grouped", np.round(se_grouped, 4),
      " ungrouped", np.round(se_ungrouped, 4))
assert np.max(np.abs(se_ungrouped / se_grouped - 1)) < 0.02

# For ungrouped binary data with the canonical link, D is a function of beta-hat alone:
# D = -2 [ beta^T X^T y + sum log(1 - mu) ] and X^T y = X^T mu.
pu = 1 / (1 + np.exp(-(Xu @ beta0)))
yb = rng.binomial(1, pu).astype(float)
b_one, mu_one = logit_irls(Xu, yb, np.ones(N))
D_direct = binomial_deviance(yb, mu_one, np.ones(N))
D_from_beta = -2 * (b_one @ (Xu.T @ mu_one) + np.sum(np.log(1 - mu_one)))
assert abs(D_direct - D_from_beta) < 1e-8

# ---- (ii) the Hauck-Donner effect --------------------------------------------
# <<hauck-donner>>
m_hd, s_a = 40, 20                                  # group A: 20 successes out of 40
s_b = np.arange(21, 40)                             # group B: from 21 to 39 successes
odds_ratio = (s_b / (m_hd - s_b)) / (s_a / (m_hd - s_a))
beta1 = np.log(odds_ratio)                          # the fitted log odds ratio
se1 = np.sqrt(1 / s_a + 1 / (m_hd - s_a) + 1 / s_b + 1 / (m_hd - s_b))
wald = (beta1 / se1) ** 2

def binom_ll(s, m, p):
    return s * np.log(p) + (m - s) * np.log(1 - p)

p_pool = (s_a + s_b) / (2 * m_hd)
lr = 2 * (binom_ll(s_a, m_hd, s_a / m_hd) + binom_ll(s_b, m_hd, s_b / m_hd)
          - binom_ll(s_a, m_hd, p_pool) - binom_ll(s_b, m_hd, p_pool))

for k in (0, 10, 15, 17, 18):
    print(f"successes in B = {s_b[k]:2d}:  log odds ratio {beta1[k]:6.3f}"
          f"   Wald {wald[k]:7.3f}   likelihood ratio {lr[k]:7.3f}")
# <</hauck-donner>>

assert np.all(np.diff(lr) > 0)                      # the likelihood ratio statistic increases
assert wald[-1] < wald[-2]                          # the Wald statistic eventually falls
peak = int(np.argmax(wald))
assert peak < len(wald) - 1

# ---- (iii) the three statistics for one coefficient --------------------------
anes = sm.datasets.anes96.load_pandas().data
y_a = anes["vote"].to_numpy(float)
pid_a = anes["PID"].to_numpy(int)
X_a = np.column_stack([np.eye(7)[pid_a], (anes["age"].to_numpy(float) - 45.0) / 10.0])
w_a = np.ones(len(y_a))

b_full, mu_full = logit_irls(X_a, y_a, w_a)
b_null, mu_null = logit_irls(X_a[:, :7], y_a, w_a)
W_full = mu_full * (1 - mu_full)
info_full = (X_a.T * W_full) @ X_a
wald_age = b_full[7] ** 2 / np.linalg.inv(info_full)[7, 7]

mu0 = mu_null                                        # fitted means under the null model
U0 = X_a.T @ (y_a - mu0)                             # the full score, at the null fit
W0 = mu0 * (1 - mu0)
score_age = U0 @ np.linalg.solve((X_a.T * W0) @ X_a, U0)

ll = lambda mu: np.sum(y_a * np.log(mu) + (1 - y_a) * np.log(1 - mu))
lr_age = 2 * (ll(mu_full) - ll(mu_null))
print(f"age coefficient: Wald {wald_age:.4f}, score {score_age:.4f}, LR {lr_age:.4f}")
assert max(abs(wald_age - lr_age), abs(score_age - lr_age)) < 0.05

gen = Generated("ch34", "inference")
gen.int("reps", reps)
gen.int("G", G)
gen.int("m", m)
gen.int("N", N)
gen.num("dev_grouped_mean", float(Dg.mean()), 2)
gen.num("dev_grouped_sd", float(Dg.std()), 2)
gen.num("dev_ungrouped_mean", float(Du.mean()), 2)
gen.num("dev_ungrouped_sd", float(Du.std()), 2)
gen.num("sqrt2df_ungrouped", float(np.sqrt(2 * (N - 2))), 2)
gen.int("peak", int(s_b[peak]))
gen.num("wald_peak", float(wald[peak]), 3)
gen.num("wald_last", float(wald[-1]), 3)
gen.num("lr_last", float(lr[-1]), 3)
gen.num("beta_last", float(beta1[-1]), 3)
gen.num("wald_age", float(wald_age), 3)
gen.num("score_age", float(score_age), 3)
gen.num("lr_age", float(lr_age), 3)
gen.write()

# ---- figures -----------------------------------------------------------------
use_book_style()
fig, axes = plt.subplots(1, 2, figsize=(5.8, 2.4))
ax = axes[0]
ax.hist(Dg, bins=40, density=True, color=COLORS["accent"], alpha=0.5, linewidth=0)
xs = np.linspace(10, 80, 300)
ax.plot(xs, stats.chi2(G - 2).pdf(xs), color=COLORS["ink"])
ax.set_xlabel("deviance")
ax.set_ylabel("density")
ax.set_title(f"(a) grouped: {G} groups of {m}", fontsize=8)
ax = axes[1]
ax.hist(Du, bins=40, density=True, color=COLORS["second"], alpha=0.5, linewidth=0)
xs = np.linspace(660, 960, 400)
ax.plot(xs, stats.chi2(N - 2).pdf(xs), color=COLORS["ink"])
ax.set_xlim(660, 1100)
ax.set_xlabel("deviance")
ax.set_title(f"(b) ungrouped: {N} binary responses", fontsize=8)
fig.tight_layout()
fig.savefig(figure_path("ch34", "deviance_null"))

fig, ax = plt.subplots(figsize=(4.0, 2.5))
ax.plot(s_b, wald, color=COLORS["accent"], marker="o", markersize=3, label="Wald")
ax.plot(s_b, lr, color=COLORS["second"], marker="s", markersize=3, linestyle="--",
        label="likelihood ratio")
ax.axhline(stats.chi2(1).ppf(0.95), color=COLORS["grid"], linewidth=0.8, zorder=0)
ax.set_xlabel("successes in group B (out of 40; group A has 20)")
ax.set_ylabel("statistic")
ax.legend(frameon=False, fontsize=7, loc="upper left")
fig.tight_layout()
fig.savefig(figure_path("ch34", "hauck_donner"))

# A smaller version of the same simulation, fast enough to run in a browser cell.
reps = 200
Dg = null_deviances(Xg, beta0, m, reps, rng)
Du = null_deviances(Xu, beta0, 1, reps, rng)
print(f"grouped   : mean {Dg.mean():8.2f} (df {G - 2}), sd {Dg.std():6.2f}"
      f" (sqrt(2 df) {np.sqrt(2 * (G - 2)):.2f})")
print(f"ungrouped : mean {Du.mean():8.2f} (df {N - 2}), sd {Du.std():6.2f}"
      f" (sqrt(2 df) {np.sqrt(2 * (N - 2)):.2f})")
assert abs(Du.mean() - (N - 2)) > 0.2 * (N - 2)
