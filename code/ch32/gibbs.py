"""Chapter 32, Section 6: a Gibbs sampler for the one-way random effects model.

The hierarchical reading of the mixed model makes every full conditional a standard
distribution, so the posterior can be sampled in a few lines. The script runs the
sampler on the proficiency study of Section 32.1 under two priors for the between-
laboratory variance, an inverse gamma with small parameters and a half-Cauchy on the
standard deviation, and compares the posterior with the REML fit and the BLUP. It then
repeats the comparison on the first five laboratories, where the prior matters.
"""
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats

from regbook import COLORS, Generated, figure_path, use_book_style

# <<sampler>>
import numpy as np

g, m = 12, 4
n = g * m
mu, sigma_a, sigma = 2.50, 0.20, 0.25           # the values used to simulate

rng = np.random.default_rng(20320048)
lab_effect = sigma_a * rng.standard_normal(g)
y = mu + np.repeat(lab_effect, m) + sigma * rng.standard_normal(n)

Y = y.reshape(g, m)                             # the REML fit of Section 32.4
ms_between = m * np.sum((Y.mean(axis=1) - Y.mean()) ** 2) / (g - 1)
ms_within = np.sum((Y - Y.mean(axis=1)[:, None]) ** 2) / (n - g)
reml_s2, reml_s2a = ms_within, (ms_between - ms_within) / m


def gibbs(y, g, m, draws=40000, burn=2000, prior="inverse-gamma", A=1.0, seed=0):
    """Posterior draws of (mu, sigma^2, sigma_a^2, a) with a flat prior on mu."""
    rng = np.random.default_rng(seed)
    Y = y.reshape(g, m)
    n = g * m
    mu_, s2, s2a, xi = Y.mean(), Y.var(), 0.1, 1.0
    out = np.empty((draws, 3 + g))
    for t in range(draws + burn):
        prec = m / s2 + 1 / s2a                                  # a_k | rest
        a = rng.normal((m / s2) * (Y.mean(axis=1) - mu_) / prec, np.sqrt(1 / prec))
        mu_ = rng.normal((y - np.repeat(a, m)).mean(), np.sqrt(s2 / n))   # mu | rest
        resid = y - mu_ - np.repeat(a, m)
        s2 = 1 / rng.gamma(0.001 + n / 2, 1 / (0.001 + resid @ resid / 2))
        if prior == "inverse-gamma":                             # sigma_a^2 | rest
            s2a = 1 / rng.gamma(0.001 + g / 2, 1 / (0.001 + a @ a / 2))
        else:                                                    # half-Cauchy(0, A)
            s2a = 1 / rng.gamma((g + 1) / 2, 1 / (1 / xi + a @ a / 2))
            xi = 1 / rng.gamma(1.0, 1 / (1 / A**2 + 1 / s2a))
        if t >= burn:
            out[t - burn] = np.concatenate([[mu_, s2, s2a], a])
    return out


# <</sampler>>

# <<demo>>
# a shorter run than the 40 000 draws quoted in the text, so the cell finishes quickly
short = gibbs(y, g, m, draws=4000, burn=500, prior="half-Cauchy", seed=12)
print(f"posterior median of sigma_a^2  {np.median(short[:, 2]):.4f}")
print(f"posterior mean of mu           {short[:, 0].mean():.4f}")
print(f"REML estimate of sigma_a^2     {reml_s2a:.4f}")
# <</demo>>

chain_ig = gibbs(y, g, m, prior="inverse-gamma", seed=11)
chain_hc = gibbs(y, g, m, prior="half-Cauchy", A=1.0, seed=12)
print(f"IG prior:  posterior median of sigma_a^2 = {np.median(chain_ig[:, 2]):.4f}")
print(f"HC prior:  posterior median of sigma_a^2 = {np.median(chain_hc[:, 2]):.4f}")
assert abs(np.median(short[:, 2]) / np.median(chain_hc[:, 2]) - 1) < 0.08

blup = (m * reml_s2a / (reml_s2 + m * reml_s2a)) * (y.reshape(g, m).mean(axis=1) - y.mean())
pev_sd = np.sqrt(reml_s2a / (reml_s2 + m * reml_s2a) * (reml_s2 + reml_s2a * m / g))

post_mean_a = chain_ig[:, 3:].mean(axis=0)
post_sd_a = chain_ig[:, 3:].std(axis=0)
assert np.max(np.abs(post_mean_a - blup)) < 0.06          # posterior mean tracks the BLUP
assert np.corrcoef(post_mean_a, blup)[0, 1] > 0.995
assert post_sd_a.mean() > pev_sd                          # and is honestly more uncertain
assert abs(chain_ig[:, 0].mean() - y.mean()) < 0.01
assert abs(np.median(chain_ig[:, 2]) / np.median(chain_hc[:, 2]) - 1) < 0.25

# ---- five laboratories: the prior now matters --------------------------------------
g5 = 5
y5 = y[: g5 * m]
small_ig = gibbs(y5, g5, m, prior="inverse-gamma", seed=21)
small_hc = gibbs(y5, g5, m, prior="half-Cauchy", A=1.0, seed=22)
ratio_small = np.median(small_ig[:, 2]) / np.median(small_hc[:, 2])
Y5 = y5.reshape(g5, m)
msb5 = m * np.sum((Y5.mean(axis=1) - Y5.mean()) ** 2) / (g5 - 1)
msw5 = np.sum((Y5 - Y5.mean(axis=1)[:, None]) ** 2) / (g5 * m - g5)
reml5 = (msb5 - msw5) / m
assert ratio_small < 0.75, ratio_small        # the inverse gamma pulls sigma_a^2 down

gen = Generated("ch32", "gibbs")
gen.int("draws", 40000)
gen.num("medig", np.median(chain_ig[:, 2]), 4)
gen.num("medhc", np.median(chain_hc[:, 2]), 4)
gen.num("meanig", chain_ig[:, 2].mean(), 4)
gen.num("mus2", chain_ig[:, 1].mean(), 4)
gen.num("mumean", chain_ig[:, 0].mean(), 4)
gen.num("musd", chain_ig[:, 0].std(), 4)
gen.num("blupmax", np.abs(blup).max(), 3)
gen.num("postmax", np.abs(post_mean_a).max(), 3)
gen.num("postsd", post_sd_a.mean(), 4)
gen.num("pevsd", pev_sd, 4)
gen.num("lam1hat", reml_s2 + m * reml_s2a, 4)
gen.num("plugse", np.sqrt((reml_s2 + m * reml_s2a) / n), 4)
gen.num("sewider", chain_ig[:, 0].std() / np.sqrt((reml_s2 + m * reml_s2a) / n), 3)
gen.num("sdwider", post_sd_a.mean() / pev_sd, 3)
gen.int("gfive", g5)
gen.num("medigfive", np.median(small_ig[:, 2]), 4)
gen.num("medhcfive", np.median(small_hc[:, 2]), 4)
gen.num("remlfive", reml5, 4)
gen.write()

# ---- figure ------------------------------------------------------------------------
use_book_style()
fig, axes = plt.subplots(1, 2, figsize=(5.8, 2.5))
for ax, (cig, chc, r, top, title) in zip(
        axes,
        [(chain_ig, chain_hc, np.sqrt(reml_s2a), 0.5, "(a) all 12 laboratories"),
         (small_ig, small_hc, np.sqrt(max(reml5, 0)), 0.8, "(b) the first 5 laboratories")]):
    grid = np.linspace(0, top, 400)
    for chain, style, label in [(cig, "-", "inverse gamma"), (chc, "--", "half-Cauchy")]:
        sd = np.sqrt(chain[:, 2])
        ax.plot(grid, stats.gaussian_kde(sd[sd < 2.0])(grid),
                style, color=COLORS["accent"] if style == "-" else COLORS["second"],
                label=label)
    ax.axvline(r, color=COLORS["ink"], lw=0.7, ls=":")
    ax.set_xlabel(r"$\sigma_a$")
    ax.set_title(title)
axes[0].set_ylabel("posterior density")
axes[1].legend(frameon=False, fontsize=7, loc="upper right")
fig.tight_layout()
fig.savefig(figure_path("ch32", "gibbs_posterior"))
