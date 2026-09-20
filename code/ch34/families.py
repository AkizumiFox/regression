"""Chapter 34, Section 1: the exponential dispersion family.

Checks mu = b'(theta) and Var = phi b''(theta) / w against the exact moments of six
standard distributions, checks the third cumulant against (phi/w)^2 b'''(theta), and
draws the variance functions.
"""
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats

from regbook import COLORS, Generated, figure_path, use_book_style

# <<table>>
import numpy as np
from scipy import stats

# Each family: theta(mu), b(theta) and its first three derivatives, and the exact
# distribution of Y, so that the exponential-family formulas can be checked.
FAMILIES = {
    "normal": dict(
        mu=2.0, phi=1.5, w=1.0,
        theta=lambda mu: mu,
        b=lambda t: t ** 2 / 2, b1=lambda t: t, b2=lambda t: 1.0 + 0 * t, b3=lambda t: 0 * t,
        dist=lambda mu, phi, w: stats.norm(mu, np.sqrt(phi / w)),
    ),
    "binomial": dict(                                  # Y = S/m, the sample proportion
        mu=0.3, phi=1.0, w=10.0,
        theta=lambda mu: np.log(mu / (1 - mu)),
        b=lambda t: np.log1p(np.exp(t)),
        b1=lambda t: 1 / (1 + np.exp(-t)),
        b2=lambda t: np.exp(-t) / (1 + np.exp(-t)) ** 2,
        b3=lambda t: np.exp(-t) * (np.exp(-t) - 1) / (1 + np.exp(-t)) ** 3,
        dist=None,                                     # handled separately: Y is m^-1 Binomial(m, mu)
    ),
    "Poisson": dict(
        mu=3.5, phi=1.0, w=1.0,
        theta=lambda mu: np.log(mu),
        b=lambda t: np.exp(t), b1=np.exp, b2=np.exp, b3=np.exp,
        dist=lambda mu, phi, w: stats.poisson(mu),
    ),
    "gamma": dict(                                     # shape a = 1/phi
        mu=4.0, phi=0.4, w=1.0,
        theta=lambda mu: -1 / mu,
        b=lambda t: -np.log(-t), b1=lambda t: -1 / t, b2=lambda t: 1 / t ** 2,
        b3=lambda t: -2 / t ** 3,
        dist=lambda mu, phi, w: stats.gamma(a=w / phi, scale=mu * phi / w),
    ),
    "inverse Gaussian": dict(
        mu=2.0, phi=0.4, w=1.0,
        theta=lambda mu: -1 / (2 * mu ** 2),
        b=lambda t: -np.sqrt(-2 * t), b1=lambda t: (-2 * t) ** -0.5,
        b2=lambda t: (-2 * t) ** -1.5, b3=lambda t: 3 * (-2 * t) ** -2.5,
        dist=lambda mu, phi, w: stats.invgauss(mu=mu * phi / w, scale=w / phi),
    ),
}

def check(name, spec, k=None):
    """Compare b'(theta), phi b''(theta)/w and (phi/w)^2 b'''(theta) with the true moments."""
    mu, phi, w = spec["mu"], spec["phi"], spec["w"]
    t = spec["theta"](mu)
    mean_ef = spec["b1"](t)
    var_ef = phi * spec["b2"](t) / w
    cum3_ef = (phi / w) ** 2 * spec["b3"](t)
    if name == "binomial":                             # Y = S/m with S ~ Binomial(m, mu)
        m = int(w)
        s = np.arange(m + 1)
        p = stats.binom(m, mu).pmf(s)
        y = s / m
        mean, var = y @ p, (y - mu) ** 2 @ p
        cum3 = (y - mu) ** 3 @ p
    else:
        d = spec["dist"](mu, phi, w)
        mean, var = d.mean(), d.var()
        cum3 = d.stats("s") * var ** 1.5               # third cumulant = skewness * sd^3
    return mean_ef, mean, var_ef, var, cum3_ef, float(cum3)

for name, spec in FAMILIES.items():
    m_ef, m, v_ef, v, c_ef, c = check(name, spec)
    print(f"{name:17s} mean {m_ef:8.4f} {m:8.4f}   var {v_ef:8.4f} {v:8.4f}   cum3 {c_ef:8.4f} {c:8.4f}")
# <</table>>

for name, spec in FAMILIES.items():
    m_ef, m, v_ef, v, c_ef, c = check(name, spec)
    assert np.isclose(m_ef, m, rtol=1e-9, atol=1e-12), name
    assert np.isclose(v_ef, v, rtol=1e-9, atol=1e-12), name
    assert np.isclose(c_ef, c, rtol=1e-6, atol=1e-9), name

# The negative binomial with known shape k: theta = log(mu/(mu+k)), b = -k log(1 - e^theta).
k_nb, mu_nb = 2.0, 3.0
theta_nb = np.log(mu_nb / (mu_nb + k_nb))
b1_nb = k_nb * np.exp(theta_nb) / (1 - np.exp(theta_nb))
b2_nb = k_nb * np.exp(theta_nb) / (1 - np.exp(theta_nb)) ** 2
nb = stats.nbinom(n=k_nb, p=k_nb / (k_nb + mu_nb))
assert np.isclose(b1_nb, nb.mean())
assert np.isclose(b2_nb, nb.var())
assert np.isclose(b2_nb, mu_nb * (1 + mu_nb / k_nb))

# The variance of Y in the binomial case is mu(1-mu)/m: the weight w = m divides it.
assert np.isclose(FAMILIES["binomial"]["phi"] * 0.3 * 0.7 / 10.0, 0.3 * 0.7 / 10.0)

gen = Generated("ch34", "families")
gen.num("ig_var", 0.4 * 2.0 ** 3, 3)
gen.num("gamma_var", 0.4 * 4.0 ** 2, 3)
gen.num("gamma_cum3", 2 * 0.4 ** 2 * 4.0 ** 3, 3)
gen.num("nb_var", float(nb.var()), 3)
gen.num("pois_cum3", 3.5, 3)
gen.write()

# ---- variance functions -----------------------------------------------------
use_book_style()
fig, axes = plt.subplots(1, 2, figsize=(5.6, 2.4))
mu = np.linspace(0.05, 6, 400)
ax = axes[0]
for label, v, col in [("normal, $V=1$", np.ones_like(mu), COLORS["ink"]),
                      ("Poisson, $V=\\mu$", mu, COLORS["accent"]),
                      ("gamma, $V=\\mu^2$", mu ** 2, COLORS["second"]),
                      ("inv. Gaussian, $V=\\mu^3$", mu ** 3, COLORS["third"])]:
    ax.plot(mu, v, color=col, label=label)
ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlabel(r"$\mu$")
ax.set_ylabel(r"$V(\mu)$")
ax.set_title("(a) the power family")
ax.legend(frameon=False, fontsize=6.5, loc="upper left")
ax = axes[1]
p = np.linspace(0.002, 0.998, 400)
ax.plot(p, p * (1 - p), color=COLORS["accent"], label=r"binomial, $V=\mu(1-\mu)$")
mu2 = np.linspace(0.02, 1.0, 400)
ax.plot(mu2, mu2, color=COLORS["ink"], linestyle="--", label=r"Poisson, $V=\mu$")
for kk, col in [(1.0, COLORS["second"]), (4.0, COLORS["third"])]:
    ax.plot(mu2, mu2 * (1 + mu2 / kk), color=col, label=fr"neg. binomial, $k={kk:.0f}$")
ax.set_xlabel(r"$\mu$")
ax.set_ylabel(r"$V(\mu)$")
ax.set_ylim(0, 2.7)
ax.set_title("(b) bounded and overdispersed")
ax.legend(frameon=False, fontsize=6.5, loc="upper left")
fig.tight_layout()
fig.savefig(figure_path("ch34", "variance_functions"))
