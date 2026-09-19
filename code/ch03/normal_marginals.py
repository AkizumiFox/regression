"""Chapter 3, Section 3: normal marginals do not make a normal vector.

An equal mixture of N_2(0, R(+rho)) and N_2(0, R(-rho)), where R(r) is the 2 x 2
correlation matrix with off-diagonal r. Both coordinates are exactly N(0, 1) and
their correlation is zero, but the pair is not bivariate normal and the
coordinates are dependent."""
import matplotlib.pyplot as plt
import numpy as np
from scipy import integrate, stats

from regbook import COLORS, Generated, figure_path, use_book_style

rng = np.random.default_rng(302)
rho, N = 0.8, 200_000

# <<mixture>>
sign = rng.choice([-1.0, 1.0], size=N)          # which component each draw uses
Z1, Z2 = rng.standard_normal((2, N))
Y1 = Z1
Y2 = sign * rho * Z1 + np.sqrt(1 - rho**2) * Z2  # corr(Y1, Y2 | sign) = sign * rho

S = (Y1 + Y2) / np.sqrt(2)                     # a direction that is not normal
print("marginal KS p-values:", stats.kstest(Y1, "norm").pvalue, stats.kstest(Y2, "norm").pvalue)
print("correlation:", np.corrcoef(Y1, Y2)[0, 1])
print("corr of squares:", np.corrcoef(Y1**2, Y2**2)[0, 1])
print("kurtosis of S:", stats.kurtosis(S, fisher=False), " theory:", 3 * (1 + rho**2))
print("KS p-value of S:", stats.kstest(S, "norm").pvalue)
# <</mixture>>

# exact check that the marginal of Y2 is N(0,1): integrate the mixture density over y1
def mix_pdf(y1, y2):
    f = lambda r: stats.multivariate_normal([0, 0], [[1, r], [r, 1]]).pdf([y1, y2])
    return 0.5 * (f(rho) + f(-rho))

for y2 in [-2.0, -0.3, 0.0, 1.1]:
    val, _ = integrate.quad(lambda t: mix_pdf(t, y2), -12, 12)
    assert np.isclose(val, stats.norm.pdf(y2), atol=1e-8)

corr = np.corrcoef(Y1, Y2)[0, 1]
corr_sq = np.corrcoef(Y1**2, Y2**2)[0, 1]
kurt = stats.kurtosis(S, fisher=False)
kurt_theory = 3 * (1 + rho**2)
p1, p2 = stats.kstest(Y1, "norm").pvalue, stats.kstest(Y2, "norm").pvalue
ks_S = stats.kstest(S, "norm")
pS = ks_S.pvalue
assert abs(corr) < 0.01                      # uncorrelated
assert abs(corr_sq - rho**2) < 0.02          # but dependent: corr(Y1^2, Y2^2) = rho^2
assert p1 > 0.01 and p2 > 0.01               # marginals look normal
assert abs(kurt - kurt_theory) < 0.15        # S is heavy tailed
assert pS < 1e-10                            # and clearly not normal
assert np.isclose(np.var(S), 1.0, atol=0.01)

gen = Generated("ch03", "normal_marginals", prefix="marg")
gen.num("rho", rho, 1)
gen.text("N", f"{N:,}".replace(",", "{,}"))
gen.num("corr", corr, 4)
gen.num("corrsq", corr_sq, 3)
gen.num("kurt", kurt, 2)
gen.num("kurttheory", kurt_theory, 2)
gen.num("p1", p1, 2)
gen.num("p2", p2, 2)
gen.num("ksS", ks_S.statistic, 3)
gen.num("corrsqtheory", rho**2, 2)
gen.write()

use_book_style()
fig, axes = plt.subplots(1, 2, figsize=(5.6, 2.55))
ax = axes[0]
k = 2500
ax.scatter(Y1[:k], Y2[:k], s=3, color=COLORS["accent"], alpha=0.35, linewidths=0)
ax.set_xlim(-3.5, 3.5); ax.set_ylim(-3.5, 3.5); ax.set_aspect("equal")
ax.set_xlabel("$Y_1$"); ax.set_ylabel("$Y_2$")
ax.set_title("(a) normal marginals, zero correlation")
ax = axes[1]
g = np.linspace(-4, 4, 400)
ax.hist(S, bins=120, range=(-4, 4), density=True, color=COLORS["grid"], edgecolor="none")
ax.plot(g, stats.norm.pdf(g), color=COLORS["second"], label=r"$\mathrm{N}(0,1)$")
dens = 0.5 * (stats.norm.pdf(g, scale=np.sqrt(1 + rho)) + stats.norm.pdf(g, scale=np.sqrt(1 - rho)))
ax.plot(g, dens, color=COLORS["ink"], ls="--", label="exact density")
ax.set_xlabel(r"$(Y_1+Y_2)/\sqrt{2}$")
ax.set_title("(b) a linear combination")
ax.set_ylim(0, 0.9)
ax.legend(frameon=False, loc="upper center", ncol=2, handlelength=1.6)
fig.tight_layout()
fig.savefig(figure_path("ch03", "normal_marginals"))
