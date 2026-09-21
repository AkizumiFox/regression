"""Chapter 38, Section 1: what a variance function fixes, and what it leaves open.

Two figures. The first contrasts the two ways a binomial variance can be inflated as
the group size grows. The second draws the quasi-deviance contribution
2 * integral_mu^y (y - t)/V(t) dt as a function of the fitted mean, for the four power
variance functions.
"""
import matplotlib.pyplot as plt
import numpy as np
from scipy import integrate

from regbook import COLORS, Generated, figure_path, use_book_style


# <<qdev>>
def quasi_deviance_term(y, mu, V):
    """2 * integral from mu to y of (y - t)/V(t) dt: one observation's quasi-deviance."""
    value, _ = integrate.quad(lambda t: (y - t) / V(t), mu, y)
    return 2 * value


y0 = 2.0
closed_form = {
    "V = 1": (lambda m: (y0 - m) ** 2, lambda t: 1.0),
    "V = mu": (lambda m: 2 * (y0 * np.log(y0 / m) - (y0 - m)), lambda t: t),
    "V = mu^2": (lambda m: 2 * (-np.log(y0 / m) + (y0 - m) / m), lambda t: t ** 2),
    "V = mu^3": (lambda m: (y0 - m) ** 2 / (y0 * m ** 2), lambda t: t ** 3),
}
for name, (formula, V) in closed_form.items():
    print(f"{name:9s}  at mu = 0.8: {quasi_deviance_term(y0, 0.8, V):.6f}"
          f"   closed form {formula(0.8):.6f}")
# <</qdev>>

for name, (formula, V) in closed_form.items():
    for m in (0.5, 0.8, 1.5, 3.0, 5.0):
        assert abs(quasi_deviance_term(y0, m, V) - formula(m)) < 1e-9, name
    assert abs(quasi_deviance_term(y0, y0, V)) < 1e-12

# The quasi-deviance is asymmetric for every V except the constant one: the ratio of
# the penalties for a mean too small and a mean too large by the same amount.
ratios = {}
for name, (formula, V) in closed_form.items():
    ratios[name] = formula(y0 - 1.0) / formula(y0 + 1.0)

# <<binomial>>
m_grid = np.arange(1, 41)
phi, rho = 1.8, 0.05
inflated = np.full_like(m_grid, phi, dtype=float)      # V(pi) = phi pi(1-pi)/m
exchangeable = 1 + rho * (m_grid - 1)                  # V(pi) = {1+rho(m-1)} pi(1-pi)/m
crossing = 1 + (phi - 1) / rho                         # the one group size where they agree
print(f"the two variance functions agree at m = {crossing:.1f}")
print(f"at m = 40 they differ by a factor {exchangeable[-1] / phi:.3f}")
# <</binomial>>

assert np.isclose(1 + rho * (crossing - 1), phi)
assert exchangeable[0] == 1.0                          # m = 1: no overdispersion possible

# ---- exp Q normalized over (0, 1) for the quasi-binomial specification ------
# Exercise C1 of Section 38.2: exp Q integrates, but the density it defines has
# neither of the two moments the specification asks for.


def quasi_binomial_moments(mu, phi):
    """Mean and variance of the density proportional to exp Q(mu; y) on (0, 1)."""
    def kernel(y):
        q = (y * np.log(mu) + (1 - y) * np.log(1 - mu)
             - y * np.log(y) - (1 - y) * np.log(1 - y))
        return np.exp(q / phi)

    norm, _ = integrate.quad(kernel, 0.0, 1.0)
    m1, _ = integrate.quad(lambda y: y * kernel(y), 0.0, 1.0)
    m2, _ = integrate.quad(lambda y: y * y * kernel(y), 0.0, 1.0)
    mean = m1 / norm
    return mean, m2 / norm - mean ** 2


qbin_mu, qbin_phi = 0.3, 1.0
qbin_mean, qbin_var = quasi_binomial_moments(qbin_mu, qbin_phi)
print(f"quasi-binomial exp Q at mu = {qbin_mu}, phi = {qbin_phi}: "
      f"mean {qbin_mean:.4f}, variance {qbin_var:.4f}, "
      f"target variance {qbin_phi * qbin_mu * (1 - qbin_mu):.4f}")

assert abs(qbin_mean - qbin_mu) > 0.1                  # the mean is not mu
assert qbin_var < qbin_mean * (1 - qbin_mean)          # a law on (0,1) cannot reach m(1-m)
for phi_test in (0.5, 1.0, 2.0):                       # symmetric at mu = 1/2, hence mean 1/2
    sym_mean, sym_var = quasi_binomial_moments(0.5, phi_test)
    assert abs(sym_mean - 0.5) < 1e-8
    assert sym_var < 0.25

gen = Generated("ch38", "variance_shapes")
gen.num("qbin_mean", qbin_mean, 4)
gen.num("qbin_var", qbin_var, 4)
gen.num("ratio_poisson", ratios["V = mu"], 3)
gen.num("ratio_gamma", ratios["V = mu^2"], 3)
gen.num("ratio_invgauss", ratios["V = mu^3"], 3)
gen.num("phi_example", phi, 1)
gen.num("rho_example", rho, 2)
gen.num("crossing", crossing, 1)
gen.num("factor_40", exchangeable[-1] / phi, 3)
gen.write()

# ---- two figures -----------------------------------------------------------
use_book_style()
fig, ax = plt.subplots(figsize=(3.4, 2.4))
m_plot = np.arange(1, 41)
ax.plot(m_plot, inflated, color=COLORS["accent"], label=r"$\phi\,\pi(1-\pi)/m$")
ax.plot(m_plot, exchangeable, color=COLORS["second"], linestyle="--",
        label=r"$\{1+\rho(m-1)\}\pi(1-\pi)/m$")
ax.axhline(1.0, color=COLORS["grid"], linewidth=0.6, zorder=0)
ax.plot([crossing], [phi], "o", color=COLORS["ink"], markersize=3)
ax.set_xlabel("group size $m$")
ax.set_ylabel(r"variance $\div$ binomial variance")
ax.legend(frameon=False, loc="upper left")
fig.tight_layout()
fig.savefig(figure_path("ch38", "binomial_inflation"))

fig, ax = plt.subplots(figsize=(3.4, 2.4))
mu_grid = np.linspace(0.35, 6.0, 400)
styles = [("V = 1", COLORS["muted"], "-"), ("V = mu", COLORS["accent"], "-"),
          ("V = mu^2", COLORS["second"], "--"), ("V = mu^3", COLORS["third"], "-.")]
labels = {"V = 1": r"$V=1$", "V = mu": r"$V=\mu$",
          "V = mu^2": r"$V=\mu^2$", "V = mu^3": r"$V=\mu^3$"}
for name, colour, dash in styles:
    formula = closed_form[name][0]
    ax.plot(mu_grid, formula(mu_grid), color=colour, linestyle=dash, label=labels[name])
ax.axvline(y0, color=COLORS["grid"], linewidth=0.6, zorder=0)
ax.set_ylim(0, 8.2)
ax.set_xlabel(r"fitted mean $\mu$")
ax.set_ylabel("quasi-deviance term")
ax.legend(frameon=False, loc="upper center", ncol=2, handlelength=1.5,
          columnspacing=1.0, borderaxespad=0.1)
fig.tight_layout()
fig.savefig(figure_path("ch38", "quasi_deviance_shapes"))
