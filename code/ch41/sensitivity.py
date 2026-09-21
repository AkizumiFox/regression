"""Chapter 41, Section 6: delta-adjustment and the tipping point.

A two-arm study with a baseline covariate loses outcomes, more often in the
treated arm and more often at high baseline values. Under the assumption that the
lost outcomes look like the recorded ones with the same arm and baseline, multiple
imputation gives a clear treatment effect. The pattern-mixture delta-adjustment
shifts the imputed outcomes of the treated arm down by delta, refits, and records
where the conclusion turns over. The observed-data log-likelihood is the same for
every delta, which is the point of the exercise.
"""
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats

from regbook import COLORS, Generated, figure_path, use_book_style

# <<data>>
rng = np.random.default_rng(4108)
n, tau_true = 400, 0.6
arm = np.repeat([0, 1], n // 2)                       # 0 = control, 1 = treated
base = rng.normal(size=n)                             # a baseline measurement, never missing
y = 0.5 + tau_true * arm + 0.6 * base + rng.normal(size=n)
eta = -1.0 + 0.7 * arm + 0.8 * base                   # dropout: depends on arm and baseline
seen = rng.uniform(size=n) > 1 / (1 + np.exp(-eta))
print(f"{n - seen.sum()} of {n} outcomes lost:"
      f" {int((~seen & (arm == 1)).sum())} treated, {int((~seen & (arm == 0)).sum())} control")
# <</data>>

Z = np.column_stack([np.ones(n), arm, base])


def mi_estimate(delta, M, rng):
    """Multiple imputation under MAR, then shift the treated arm's imputations by -delta."""
    Zo, k = Z[seen], int(seen.sum())
    g, *_ = np.linalg.lstsq(Zo, y[seen], rcond=None)
    sse = float(np.sum((y[seen] - Zo @ g) ** 2))
    qs, us = [], []
    for _ in range(M):
        s2 = sse / rng.chisquare(k - Z.shape[1])
        g_star = rng.multivariate_normal(g, s2 * np.linalg.inv(Zo.T @ Zo))
        draw = Z @ g_star + np.sqrt(s2) * rng.normal(size=n) - delta * arm
        y_full = np.where(seen, y, draw)
        b, *_ = np.linalg.lstsq(Z, y_full, rcond=None)
        resid = y_full - Z @ b
        v = (resid @ resid / (n - 3)) * np.linalg.inv(Z.T @ Z)[1, 1]
        qs.append(b[1])
        us.append(v)
    M_ = len(qs)
    q_bar, u_bar, b_var = np.mean(qs), np.mean(us), np.var(qs, ddof=1)
    total = u_bar + (1 + 1 / M_) * b_var
    gamma = (1 + 1 / M_) * b_var / total
    return q_bar, total, (M_ - 1) / gamma**2


# <<tipping>>
deltas = np.linspace(0.0, 1.6, 33)
est, lo, hi = [], [], []
for d in deltas:
    q, t, nu = mi_estimate(d, 50, np.random.default_rng(4109))
    half = stats.t.ppf(0.975, nu) * np.sqrt(t)
    est.append(q)
    lo.append(q - half)
    hi.append(q + half)
est, lo, hi = map(np.array, (est, lo, hi))
tipping = float(np.interp(0.0, -lo, deltas))          # where the lower limit reaches zero
print(f"delta = 0: effect {est[0]:.3f}, 95% interval ({lo[0]:.3f}, {hi[0]:.3f})")
print(f"tipping point: delta = {tipping:.3f}")
# <</tipping>>

sd_resid = float(np.std(y[seen] - Z[seen] @ np.linalg.lstsq(Z[seen], y[seen], rcond=None)[0],
                        ddof=3))
assert lo[0] > 0 and tipping > 0.3
assert np.all(np.diff(est) < 0)                       # the effect falls steadily with delta

k = int(seen.sum())


def observed_loglik(delta):
    """Observed-data log-likelihood of the pattern-mixture model carrying this shift.

    The shift names f(y | z, R = 0); it enters the fitted mean of the unrecorded rows
    and of no other row, so the recorded outcomes never see it.
    """
    g_o, *_ = np.linalg.lstsq(Z[seen], y[seen], rcond=None)
    s2_o = float(np.sum((y[seen] - Z[seen] @ g_o) ** 2)) / k
    mu = np.where(seen, Z @ g_o, Z @ g_o - delta * arm)
    return float(-0.5 * np.sum(np.log(2 * np.pi * s2_o) + (y[seen] - mu[seen]) ** 2 / s2_o))


logliks = np.array([observed_loglik(d) for d in deltas])
assert np.allclose(logliks, logliks[0])               # every delta fits the records equally well
loglik = float(logliks[0])
print(f"observed-data log-likelihood {loglik:.3f}, the same for every delta")

gen = Generated("ch41", "sensitivity")
gen.int("n", n)
gen.int("nmis", int(n - seen.sum()))
gen.int("nmis_t", int((~seen & (arm == 1)).sum()))
gen.int("nmis_c", int((~seen & (arm == 0)).sum()))
gen.num("est0", est[0], 3)
gen.num("lo0", lo[0], 3)
gen.num("hi0", hi[0], 3)
gen.num("tipping", tipping, 3)
gen.num("sd_resid", sd_resid, 3)
gen.num("tipping_sd", tipping / sd_resid, 3)
gen.num("loglik", loglik, 2)
gen.write()

# ---- figure ------------------------------------------------------------------
use_book_style()
fig, ax = plt.subplots(figsize=(4.2, 2.4))
ax.fill_between(deltas, lo, hi, color=COLORS["accent"], alpha=0.2, linewidth=0)
ax.plot(deltas, est, color=COLORS["accent"])
ax.plot(deltas, lo, color=COLORS["accent"], linewidth=0.6)
ax.plot(deltas, hi, color=COLORS["accent"], linewidth=0.6)
ax.axhline(0.0, color=COLORS["ink"], linewidth=0.7)
ax.axvline(tipping, color=COLORS["second"], linestyle="--", linewidth=0.9)
ax.annotate(f"tipping point $\\delta={tipping:.2f}$", xy=(tipping, 0.02),
            xytext=(tipping + 0.05, 0.45), fontsize=7, ha="left", color=COLORS["second"])
ax.set_xlabel(r"shift $\delta$ applied to the treated arm's imputations")
ax.set_ylabel("estimated treatment effect")
fig.tight_layout()
fig.savefig(figure_path("ch41", "tipping_point"))
