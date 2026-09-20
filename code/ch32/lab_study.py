"""Chapter 32, Section 1: the one-way random effects model on a simulated proficiency study.

Twelve laboratories each analyse four aliquots of one homogeneous batch of soil for lead
(mg/kg). The laboratories are a sample from the population of accredited laboratories, so
the laboratory effect is random and the targets are the two variance components. The data
are simulated with a fixed seed, so the true values are known and the script can check
unbiasedness, coverage and the frequency of negative estimates by Monte Carlo.
"""
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats

from regbook import COLORS, Generated, figure_path, use_book_style

# <<data>>
import numpy as np
from scipy import stats

g, m = 12, 4                                    # laboratories, aliquots each
n = g * m
mu, sigma_a, sigma = 2.50, 0.20, 0.25           # the values used to simulate

rng = np.random.default_rng(20320048)
lab_effect = sigma_a * rng.standard_normal(g)
y = mu + np.repeat(lab_effect, m) + sigma * rng.standard_normal(n)

Y = y.reshape(g, m)                             # one row per laboratory
lab_mean, grand_mean = Y.mean(axis=1), Y.mean()

ss_between = m * np.sum((lab_mean - grand_mean) ** 2)
ss_within = np.sum((Y - lab_mean[:, None]) ** 2)
ms_between = ss_between / (g - 1)
ms_within = ss_within / (n - g)

sigma2_hat = ms_within                          # E(MS within) = sigma^2
sigma2_a_hat = (ms_between - ms_within) / m     # E(MS between) = sigma^2 + m sigma_a^2
icc_hat = sigma2_a_hat / (sigma2_a_hat + sigma2_hat)

print(f"MS between {ms_between:.4f}   MS within {ms_within:.4f}")
print(f"sigma^2 {sigma2_hat:.4f}   sigma_a^2 {sigma2_a_hat:.4f}   ICC {icc_hat:.3f}")
# <</data>>

# <<interval>>
f_obs = ms_between / ms_within                  # ~ (1 + m gamma) F(g-1, n-g)
gamma_lo = (f_obs / stats.f.ppf(0.975, g - 1, n - g) - 1) / m
gamma_hi = (f_obs / stats.f.ppf(0.025, g - 1, n - g) - 1) / m
icc_lo, icc_hi = gamma_lo / (1 + gamma_lo), gamma_hi / (1 + gamma_hi)

print(f"F = {f_obs:.2f};  gamma in ({gamma_lo:.3f}, {gamma_hi:.3f})")
print(f"intraclass correlation in ({icc_lo:.3f}, {icc_hi:.3f})")
# <</interval>>

gamma_true = sigma_a**2 / sigma**2
icc_true = gamma_true / (1 + gamma_true)

# ---- the estimators are unbiased, and the interval covers at its nominal rate ------
B = 40000
sim_rng = np.random.default_rng(7)
a_sim = sigma_a * sim_rng.standard_normal((B, g))
Ysim = mu + a_sim[:, :, None] + sigma * sim_rng.standard_normal((B, g, m))
lm = Ysim.mean(axis=2)
msb = m * np.sum((lm - lm.mean(axis=1, keepdims=True)) ** 2, axis=1) / (g - 1)
msw = np.sum((Ysim - lm[:, :, None]) ** 2, axis=(1, 2)) / (n - g)
sa2 = (msb - msw) / m
assert abs(sa2.mean() - sigma_a**2) < 0.003, sa2.mean()
assert abs(msw.mean() - sigma**2) < 0.0005, msw.mean()

fsim = msb / msw
lo = (fsim / stats.f.ppf(0.975, g - 1, n - g) - 1) / m
hi = (fsim / stats.f.ppf(0.025, g - 1, n - g) - 1) / m
coverage = np.mean((lo <= gamma_true) & (gamma_true <= hi))
assert abs(coverage - 0.95) < 0.006, coverage

# ---- negative estimates -----------------------------------------------------------
def p_negative(g, m, icc):
    """P(sigma_a^2 hat < 0) = P(F(g-1, gm-g) < 1/(1 + m gamma))."""
    gamma = icc / (1 - icc)
    return stats.f.cdf(1 / (1 + m * gamma), g - 1, g * m - g)


assert abs(np.mean(sa2 < 0) - p_negative(g, m, icc_true)) < 0.004
assert abs(p_negative(6, 4, 0.0) - 0.5) < 0.06        # near one half when sigma_a = 0

# ---- exercise: which group size estimates sigma_a^2 best, with n held fixed? -------
def var_hat(m_, n_=100, s2=1.0, s2a=0.01):
    """Var(sigma_a^2 hat) in the balanced normal model, from MSB and MSE."""
    g_ = n_ // m_
    lam1 = s2 + m_ * s2a
    return (2 / m_**2) * (lam1**2 / (g_ - 1) + s2**2 / (n_ - g_))


small = {m_: var_hat(m_) for m_ in (2, 5, 10, 25)}
large = {m_: var_hat(m_, s2a=1.0) for m_ in (2, 5, 10)}
# with a small intraclass correlation few large groups win; with a large one, many small
assert min(small, key=small.get) == 25 and min(large, key=large.get) == 2
assert small[2] / small[25] > 12

gen = Generated("ch32", "lab_study")
gen.int("g", g)
gen.int("m", m)
gen.int("n", n)
gen.num("msb", ms_between, 4)
gen.num("msw", ms_within, 4)
gen.num("sigma2", sigma2_hat, 4)
gen.num("sigma2a", sigma2_a_hat, 4)
gen.num("sda", np.sqrt(sigma2_a_hat), 3)
gen.num("sde", np.sqrt(sigma2_hat), 3)
gen.num("icc", icc_hat, 3)
gen.num("f", f_obs, 2)
gen.num("gammalo", gamma_lo, 3)
gen.num("gammahi", gamma_hi, 3)
gen.num("icclo", icc_lo, 3)
gen.num("icchi", icc_hi, 3)
gen.num("grand", grand_mean, 4)
gen.num("iccpop", icc_true, 3)
gen.num("pneg", p_negative(g, m, icc_true), 4)
gen.num("pnegsmall", p_negative(6, 4, 0.05), 3)
gen.num("pnegzero", p_negative(6, 4, 0.0), 3)
for m_, v in small.items():
    gen.num(f"varsmall{m_}", v, 4)
for m_, v in large.items():
    gen.num(f"varlarge{m_}", v, 4)
gen.write()

# ---- figure -----------------------------------------------------------------------
use_book_style()
fig, axes = plt.subplots(1, 2, figsize=(5.8, 2.5))
ax = axes[0]
for k in range(g):
    ax.plot(np.full(m, k + 1), Y[k], "o", ms=3, color=COLORS["accent"], alpha=0.75, mew=0)
    ax.plot([k + 0.65, k + 1.35], [lab_mean[k]] * 2, color=COLORS["second"], lw=1.2)
ax.axhline(grand_mean, color=COLORS["ink"], lw=0.7, ls="--")
ax.set_xticks(range(1, g + 1, 2))
ax.set_xlabel("laboratory")
ax.set_ylabel("lead (mg/kg)")
ax.set_title("(a) the proficiency study")

ax = axes[1]
iccs = np.linspace(0.0005, 0.5, 300)
for gg, style in [(6, "-"), (12, "--"), (25, ":")]:
    ax.plot(iccs, [p_negative(gg, m, r) for r in iccs], style,
            color=COLORS["accent"], label=f"$g={gg}$")
ax.set_xlabel("intraclass correlation")
ax.set_ylabel("probability of a negative estimate")
ax.set_ylim(0, 0.55)
ax.legend(frameon=False)
ax.set_title("(b) $\\hat\\sigma_a^2<0$, with $m=4$")
fig.tight_layout()
fig.savefig(figure_path("ch32", "lab_study"))
