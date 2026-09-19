"""Chapter 18, Section 1: how much a covariate buys in a randomized experiment.

Repeated synthetic experiments with the design of the tutoring example (four groups of
ten). The covariate is normal and independent of the assignment, and y follows the
parallel-lines model. Compares the variance of the unadjusted and the adjusted estimate
of mu_1 - mu_2 with the formula (1 - rho^2) (1 + 1 / (N - g - 2)).
"""
import matplotlib.pyplot as plt
import numpy as np

from regbook import COLORS, Generated, figure_path, use_book_style

# <<simulate>>
g, m = 4, 10
N = g * m
group = np.repeat(np.arange(g), m)
Z = np.eye(g)[group]

def one_experiment(rng, rho, reps):
    """Unadjusted and adjusted estimates of mu_1 - mu_2 (true value 0) in reps experiments."""
    x = rng.normal(size=(reps, N))                       # covariate, unit variance
    e = rng.normal(size=(reps, N))
    y = rho * x + np.sqrt(1 - rho ** 2) * e              # within-group corr(x, y) = rho, var(y) = 1
    xb, yb = x @ Z / m, y @ Z / m                        # group means, reps x g
    xw, yw = x - xb[:, group], y - yb[:, group]          # within-group deviations
    beta = np.sum(xw * yw, axis=1) / np.sum(xw ** 2, axis=1)
    unadj = yb[:, 0] - yb[:, 1]
    adj = unadj - beta * (xb[:, 0] - xb[:, 1])
    return unadj, adj

rng = np.random.default_rng(1812)
rhos = np.array([0.0, 0.2, 0.4, 0.6, 0.8, 0.9])
ratio_sim = []
for rho in rhos:
    unadj, adj = one_experiment(rng, rho, 100_000)
    ratio_sim.append(adj.var() / unadj.var())
ratio_formula = (1 - rhos ** 2) * (1 + 1 / (N - g - 2))
for rho, r_s, r_f in zip(rhos, ratio_sim, ratio_formula):
    print(f"rho = {rho:.1f}: simulated variance ratio {r_s:.3f}, formula {r_f:.3f}")
# <</simulate>>

ratio_sim = np.array(ratio_sim)
assert np.allclose(ratio_sim, ratio_formula, rtol=0.03)
rho_break = np.sqrt(1 / (N - g - 1))                     # (1 - rho^2)(1 + 1/(N-g-2)) = 1
assert np.isclose((1 - rho_break ** 2) * (1 + 1 / (N - g - 2)), 1)

gen = Generated("ch18", "precision")
gen.num("infl", 1 + 1 / (N - g - 2), 4)
gen.num("rhobreak", rho_break, 3)
for rho, r_s, r_f in zip(rhos, ratio_sim, ratio_formula):
    tag = f"{int(round(rho * 10))}"
    gen.num(f"sim{tag}", r_s, 3)
    gen.num(f"form{tag}", r_f, 3)
gen.write()

use_book_style()
fig, ax = plt.subplots(figsize=(4.4, 2.7))
grid = np.linspace(0, 0.97, 200)
ax.plot(grid, (1 - grid ** 2) * (1 + 1 / (N - g - 2)), color=COLORS["accent"],
        label=r"$(1-\rho^2)\,(1+1/(N-g-2))$")
ax.plot(grid, 1 - grid ** 2, color=COLORS["muted"], linestyle="--", linewidth=0.9,
        label=r"$1-\rho^2$ (slope known)")
ax.scatter(rhos, ratio_sim, s=18, color=COLORS["second"], zorder=3, label="simulation")
ax.axhline(1, color=COLORS["grid"], linewidth=0.7, zorder=0)
ax.set_xlabel(r"within-group correlation $\rho$ of covariate and response")
ax.set_ylabel("variance, adjusted / unadjusted")
ax.set_ylim(0, 1.12)
ax.legend(frameon=False, fontsize=7, loc="lower left")
fig.tight_layout()
fig.savefig(figure_path("ch18", "precision_gain"))
