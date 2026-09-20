"""Chapter 32, Section 5: the null distribution of the REML likelihood ratio at the boundary.

In the balanced one-way model the restricted likelihood ratio statistic for
H0: sigma_a^2 = 0 is an explicit function of the analysis of variance F ratio, so its
null distribution is exact and needs no asymptotics. The script checks the formula
against a numerical maximization of the restricted likelihood, reports the point mass
at zero, and compares the exact distribution with the chi-bar-squared mixture and with
the chi-square distribution on one degree of freedom.
"""
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats

from regbook import COLORS, Generated, figure_path, use_book_style

# <<statistic>>
import numpy as np
from scipy import stats


def lrt_from_f(f, g, m):
    """Restricted likelihood ratio statistic for sigma_a^2 = 0, from the F ratio."""
    n = g * m
    a, b = (g - 1) / (n - 1), (n - g) / (n - 1)
    return np.where(f <= 1, 0.0,
                    (n - 1) * np.log(a * np.maximum(f, 1) + b) - (g - 1) * np.log(np.maximum(f, 1)))


g, m = 6, 4
n = g * m
p_zero = stats.f.cdf(1.0, g - 1, n - g)             # P(L = 0) = P(F <= 1)
crit_mix = stats.chi2.ppf(0.90, 1)                  # 5% point of the 50:50 mixture
print(f"g = {g}, m = {m}:  P(L = 0) = {p_zero:.3f},  5% point of the mixture {crit_mix:.3f}")
# <</statistic>>

# exact tail probability of L: L is increasing in F for F > 1, so invert by bisection
def tail(x, g, m):
    """P(L > x) under the null."""
    if x <= 0:
        return 1 - stats.f.cdf(1.0, g - 1, g * m - g)
    lo, hi = 1.0, 1e6
    for _ in range(200):
        mid = (lo + hi) / 2
        if lrt_from_f(mid, g, m) < x:
            lo = mid
        else:
            hi = mid
    return stats.f.sf((lo + hi) / 2, g - 1, g * m - g)


size_mix = tail(crit_mix, g, m)
size_chi1 = tail(stats.chi2.ppf(0.95, 1), g, m)
assert size_mix < 0.05 and size_chi1 < size_mix

# ---- the formula agrees with a direct maximization of the restricted likelihood ----
def reml_lrt_direct(ssb, sse, g, m):
    """-2 log ratio, from the closed-form maxima of the restricted likelihood."""
    n = g * m
    msb, mse = ssb / (g - 1), sse / (n - g)
    null = (n - 1) * np.log((ssb + sse) / (n - 1)) + (n - 1)
    if msb <= mse:
        return 0.0
    alt = (g - 1) * np.log(msb) + (n - g) * np.log(mse) + (n - 1)
    return null - alt


rng = np.random.default_rng(3141)
for _ in range(200):
    ssb = stats.chi2.rvs(g - 1, random_state=rng)
    sse = stats.chi2.rvs(n - g, random_state=rng)
    f = (ssb / (g - 1)) / (sse / (n - g))
    assert abs(reml_lrt_direct(ssb, sse, g, m) - float(lrt_from_f(f, g, m))) < 1e-9

# ---- and with a simulation of the model itself ------------------------------------
B = 200000
sim = np.random.default_rng(90210)
Ysim = sim.standard_normal((B, g, m))               # sigma_a = 0, sigma = 1
lm = Ysim.mean(axis=2)
msb = m * np.sum((lm - lm.mean(axis=1, keepdims=True)) ** 2, axis=1) / (g - 1)
msw = np.sum((Ysim - lm[:, :, None]) ** 2, axis=(1, 2)) / (n - g)
Lsim = lrt_from_f(msb / msw, g, m)
assert abs(np.mean(Lsim == 0) - p_zero) < 0.005
assert abs(np.mean(Lsim > crit_mix) - size_mix) < 0.004

# a larger design, where the mixture is close
g2, m2 = 40, 4
p_zero2 = stats.f.cdf(1.0, g2 - 1, g2 * m2 - g2)
size_mix2 = tail(crit_mix, g2, m2)
assert abs(size_mix2 - 0.05) < abs(size_mix - 0.05)

gen = Generated("ch32", "boundary")
gen.int("g", g)
gen.int("m", m)
gen.int("n", n)
gen.num("pzero", p_zero, 3)
gen.num("critmix", crit_mix, 3)
gen.num("critchi", stats.chi2.ppf(0.95, 1), 3)
gen.num("sizemix", size_mix, 4)
gen.num("sizechi", size_chi1, 4)
gen.int("g2", g2)
gen.num("pzerotwo", p_zero2, 3)
gen.num("sizemixtwo", size_mix2, 4)
gen.num("simzero", np.mean(Lsim == 0), 3)
gen.write()

# ---- figure -----------------------------------------------------------------------
use_book_style()
xs = np.linspace(0.001, 8, 300)
fig, axes = plt.subplots(1, 2, figsize=(5.8, 2.5), sharey=True)
for ax, (gg, mm, title) in zip(axes, [(g, m, "(a) $g=6$, $m=4$"),
                                      (g2, m2, "(b) $g=40$, $m=4$")]):
    ax.plot(xs, [tail(x, gg, mm) for x in xs], color=COLORS["accent"], label="exact")
    ax.plot(xs, 0.5 * stats.chi2.sf(xs, 1), "--", color=COLORS["second"],
            label=r"$\frac{1}{2}\chi^2(0)+\frac{1}{2}\chi^2(1)$")
    ax.plot(xs, stats.chi2.sf(xs, 1), ":", color=COLORS["third"], label=r"$\chi^2(1)$")
    ax.set_yscale("log")
    ax.set_ylim(2e-3, 1)
    ax.set_xlabel("$L$")
    ax.set_title(title)
axes[0].set_ylabel(r"$\Pr(L>x)$ under $\sigma_a^2=0$")
axes[1].legend(frameon=False, fontsize=7)
fig.tight_layout()
fig.savefig(figure_path("ch32", "boundary_lrt"))
