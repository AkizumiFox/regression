"""Chapter 4, Section 1: the noncentral chi-squared distribution.

Checks the Poisson-mixture representation of chi^2(r, gamma) against scipy's ncx2
(whose parameter nc is exactly the book's gamma = mu^T mu), the mean and variance
formulas, and the stochastic monotonicity in gamma. Draws the densities.
"""
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats

from regbook import COLORS, Generated, figure_path, use_book_style

rng = np.random.default_rng(20240401)

# <<mixture>>
def ncx2_pdf_mixture(x, r, gamma, terms=200):
    """Density of chi^2(r, gamma) as a Poisson(gamma/2) mixture of central densities."""
    k = np.arange(terms)
    weights = stats.poisson.pmf(k, gamma / 2)            # P(K = k)
    dens = stats.chi2.pdf(np.asarray(x)[:, None], r + 2 * k)
    return dens @ weights


r, gamma = 4, 6.0
x = np.linspace(0.05, 40, 400)
print(np.max(np.abs(ncx2_pdf_mixture(x, r, gamma) - stats.ncx2.pdf(x, r, gamma))))
# <</mixture>>

assert np.allclose(ncx2_pdf_mixture(x, r, gamma), stats.ncx2.pdf(x, r, gamma), atol=1e-12)

# mean r + gamma and variance 2r + 4 gamma (the book's convention: no factor 1/2)
m, v = stats.ncx2.stats(r, gamma, moments="mv")
assert np.isclose(m, r + gamma) and np.isclose(v, 2 * r + 4 * gamma)

# the definition: Z ~ N_r(mu, I) gives Z^T Z ~ chi^2(r, mu^T mu), whatever the direction of mu
mu = np.array([2.0, -1.0, 0.5, 0.0])
mu *= np.sqrt(gamma) / np.linalg.norm(mu)
Z = rng.standard_normal((200_000, r)) + mu
q = np.sum(Z**2, axis=1)
assert abs(q.mean() - (r + gamma)) < 4 * np.sqrt((2 * r + 4 * gamma) / len(q))
assert stats.kstest(q, stats.ncx2(r, gamma).cdf).pvalue > 1e-3

# stochastic monotonicity in gamma: P(chi^2(r, gamma) > c) increases with gamma
c = stats.chi2.ppf(0.95, r)
grid = np.linspace(0, 30, 301)
tail = stats.ncx2.sf(c, r, np.maximum(grid, 1e-12))
assert np.all(np.diff(tail) > 0)

# mixture weights for the text
w = stats.poisson.pmf(np.arange(4), gamma / 2)

gen = Generated("ch04", "noncentral_chisq", prefix="ncx")
gen.int("r", r)
gen.num("gamma", gamma, 0)
gen.num("mean", r + gamma, 0)
gen.num("var", 2 * r + 4 * gamma, 0)
gen.num("poismean", gamma / 2, 0)
gen.int("dfone", r + 2)
gen.int("dftwo", r + 4)
gen.int("dfthree", r + 6)
gen.num("crit", c, 3)
gen.num("tail", stats.ncx2.sf(c, r, gamma), 3)
gen.num("wzero", w[0], 3)
gen.num("wone", w[1], 3)
gen.num("wtwo", w[2], 3)
gen.num("wthree", w[3], 3)
gen.num("simmean", q.mean(), 3)
gen.write()

# ---- figure: densities of chi^2(4, gamma) ---------------------------------------
use_book_style()
fig, ax = plt.subplots(figsize=(4.4, 2.4))
xs = np.linspace(0.01, 36, 600)
styles = [(0, COLORS["ink"]), (3, COLORS["accent"]), (8, COLORS["third"]), (16, COLORS["second"])]
for g, col in styles:
    pdf = stats.chi2.pdf(xs, r) if g == 0 else stats.ncx2.pdf(xs, r, g)
    ax.plot(xs, pdf, color=col, label=rf"$\gamma={g}$")
    ax.axvline(r + g, color=col, linewidth=0.6, linestyle=":", ymax=0.25)
ax.set_xlabel("$x$")
ax.set_ylabel("density")
ax.set_xlim(0, 36)
ax.set_ylim(bottom=0)
ax.legend(frameon=False, title=r"$\chi^2(4,\gamma)$", title_fontsize=8)
fig.tight_layout()
fig.savefig(figure_path("ch04", "noncentral_chisq"))
