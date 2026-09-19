"""Chapter 22, Section 4: how well the classical variance-stabilizing transformations work.

Exact moments, computed by summing probability mass functions (no simulation):
(a) Poisson counts: 4 Var(sqrt Y) and 4 Var(sqrt(Y + 3/8)) against the mean.
(b) Binomial proportions with m = 20 trials: 4m Var(arcsin sqrt(S/m)) and the version with
    Anscombe's constants, 4m Var(arcsin sqrt((S + 3/8)/(m + 3/4))).
(c) Gamma data with fixed shape k (constant coefficient of variation): Var(log Y) = trigamma(k)
    exactly, whatever the mean, against the delta-method value 1/k.
"""
import matplotlib.pyplot as plt
import numpy as np
from scipy import optimize
from scipy import special
from scipy import stats

from regbook import COLORS, Generated, figure_path, use_book_style


# <<exact>>
def poisson_var(h, mu):
    """Var h(Y) for Y ~ Poisson(mu), by summing the probability mass function."""
    k = np.arange(0, int(mu + 40 * np.sqrt(mu) + 60))
    p = stats.poisson.pmf(k, mu)
    m1 = p @ h(k)
    return p @ (h(k) - m1) ** 2


def binomial_var(h, m, prob):
    """Var h(S) for S ~ Binomial(m, prob)."""
    s = np.arange(m + 1)
    p = stats.binom.pmf(s, m, prob)
    m1 = p @ h(s)
    return p @ (h(s) - m1) ** 2


for mu in (1, 2, 5, 10, 20):
    print(f"Poisson mean {mu:2d}:  4 Var sqrt(Y) = {4 * poisson_var(np.sqrt, mu):.3f}"
          f"   4 Var sqrt(Y + 3/8) = {4 * poisson_var(lambda k: np.sqrt(k + 3 / 8), mu):.3f}")
m = 20
for prob in (0.05, 0.1, 0.25, 0.5):
    plain = 4 * m * binomial_var(lambda s: np.arcsin(np.sqrt(s / m)), m, prob)
    anscombe = 4 * m * binomial_var(lambda s: np.arcsin(np.sqrt((s + 3 / 8) / (m + 3 / 4))), m, prob)
    print(f"binomial m = {m}, p = {prob:.2f}:  4m Var = {plain:.3f}   with 3/8: {anscombe:.3f}")
# <</exact>>

pv = {mu: 4 * poisson_var(np.sqrt, mu) for mu in (1, 2, 5, 10, 20)}
pa = {mu: 4 * poisson_var(lambda k: np.sqrt(k + 3 / 8), mu) for mu in (1, 2, 5, 10, 20)}
bv = {q: 4 * m * binomial_var(lambda s: np.arcsin(np.sqrt(s / m)), m, q) for q in (0.05, 0.1, 0.25, 0.5)}
ba = {q: 4 * m * binomial_var(lambda s: np.arcsin(np.sqrt((s + 3 / 8) / (m + 3 / 4))), m, q)
      for q in (0.05, 0.1, 0.25, 0.5)}
raw_ratio = poisson_var(lambda k: k, 20) / poisson_var(lambda k: k, 1)       # untransformed: 20
assert np.isclose(raw_ratio, 20)
assert all(abs(pa[mu] - 1) < abs(pv[mu] - 1) for mu in (1, 2, 5, 10, 20))  # Anscombe's shift helps
# where Anscombe's shift starts to win: the last crossing of |4 Var sqrt(Y+3/8) - 1| and |4 Var sqrt(Y) - 1|
gap = lambda u: abs(4 * poisson_var(lambda k: np.sqrt(k + 3 / 8), u) - 1) - abs(4 * poisson_var(np.sqrt, u) - 1)
grid = np.linspace(0.05, 3, 591)
signs = np.array([gap(u) for u in grid])
last_pos = np.nonzero(signs >= 0)[0].max()
crossing = optimize.brentq(gap, grid[last_pos], grid[last_pos + 1])
assert 0.6 < crossing < 0.7 and gap(crossing + 0.01) < 0
assert all(gap(u) < 0 for u in np.linspace(crossing + 0.01, 60, 400))   # wins for every larger mean
assert max(4 * poisson_var(np.sqrt, u) for u in np.linspace(1, 20, 400)) < 1.7
assert abs(pv[20] - 1) < 0.03 and abs(pa[20] - 1) < 0.002
assert abs(ba[0.25] - 1) < abs(bv[0.25] - 1) and abs(bv[0.5] - 1) < 0.1

# second order: 32 mu (Var sqrt(Y + c) - 1/4) -> 3 - 8c, which vanishes at Anscombe's c = 3/8
for c in (0, 1 / 8, 3 / 8, 1):
    assert abs(32 * 1000 * (poisson_var(lambda k: np.sqrt(k + c), 1000) - 0.25) - (3 - 8 * c)) < 0.01

# constant coefficient of variation: gamma with shape k; Var log Y does not depend on the mean
for k in (2, 5, 20):
    for scale in (0.1, 1, 30):
        lo, hi = stats.gamma.ppf([1e-12, 1 - 1e-12], k, scale=scale)
        grid = np.linspace(np.log(lo), np.log(hi), 200001)
        dens = stats.gamma.pdf(np.exp(grid), k, scale=scale) * np.exp(grid)   # density of log Y
        dx = grid[1] - grid[0]
        mean_log = np.sum(grid * dens) * dx
        var_log = np.sum((grid - mean_log) ** 2 * dens) * dx
        assert np.isclose(var_log, special.polygamma(1, k), rtol=1e-6)

gen = Generated("ch22", "variance_stabilizing", prefix="vs")
for mu in (1, 2, 5, 10, 20):
    gen.num(f"pv{mu}", pv[mu], 3)
    gen.num(f"pa{mu}", pa[mu], 3)
for q, key in [(0.05, "05"), (0.1, "10"), (0.25, "25"), (0.5, "50")]:
    gen.num(f"bv{key}", bv[q], 3)
    gen.num(f"ba{key}", ba[q], 3)
gen.num("trigamma5", special.polygamma(1, 5), 4)
gen.num("anscombe_crossing", crossing, 2)
gen.write()

# ---- figure -----------------------------------------------------------------------------------
use_book_style()
fig, axes = plt.subplots(1, 2, figsize=(5.6, 2.3))
mus = np.geomspace(0.3, 60, 120)
ax = axes[0]
ax.plot(mus, [4 * poisson_var(np.sqrt, u) for u in mus], color=COLORS["accent"], label=r"$\sqrt{y}$")
ax.plot(mus, [4 * poisson_var(lambda k: np.sqrt(k + 3 / 8), u) for u in mus], color=COLORS["second"],
        label=r"$\sqrt{y+3/8}$")
ax.axhline(1, color=COLORS["muted"], linewidth=0.6, linestyle="--")
ax.set_xscale("log")
ax.set_xlabel(r"Poisson mean $\mu$")
ax.set_ylabel(r"$4\,\mathrm{Var}$ of transformed count")
ax.set_title("(a) square root, Poisson")
ax.legend(frameon=False)
ax = axes[1]
ps = np.linspace(0.01, 0.99, 197)
ax.plot(ps, [4 * m * binomial_var(lambda s: np.arcsin(np.sqrt(s / m)), m, q) for q in ps],
        color=COLORS["accent"], label="plain")
ax.plot(ps, [4 * m * binomial_var(lambda s: np.arcsin(np.sqrt((s + 3 / 8) / (m + 3 / 4))), m, q) for q in ps],
        color=COLORS["second"], label="with 3/8, 3/4")
ax.axhline(1, color=COLORS["muted"], linewidth=0.6, linestyle="--")
ax.set_xlabel(r"success probability $p$  ($m=20$)")
ax.set_ylabel(r"$4m\,\mathrm{Var}$ of arcsine")
ax.set_title("(b) arcsine, binomial")
ax.legend(frameon=False)
fig.tight_layout()
fig.savefig(figure_path("ch22", "variance_stabilizing"))
