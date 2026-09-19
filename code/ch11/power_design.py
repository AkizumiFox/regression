"""Chapter 11, Section 5: power, noncentrality and sample size.

(1) One-way layout with g = 4 groups of m observations, testing equal means at level 0.05. For a
    given range Delta of the means, the noncentrality is smallest when two means sit at the
    extremes and the others at the midpoint, and largest when the means split evenly between
    the extremes. Power curves in m for these and for equally spaced means; smallest m for power 0.9.
(2) The region test of Section 11.1: power against a Northeast shift of delta, sigma set at s.
(3) 'Observed power' is a decreasing function of the p-value.
"""
import matplotlib.pyplot as plt
import numpy as np
import statsmodels.api as sm
from scipy import stats

from regbook import COLORS, Generated, figure_path, use_book_style

rng = np.random.default_rng(1105)

# <<oneway>>
def power(gamma, q, nu, alpha=0.05):
    """P{F(q, nu, gamma) > F_alpha(q, nu)}."""
    crit = stats.f.ppf(1 - alpha, q, nu)
    return stats.ncf.sf(crit, q, nu, gamma) if gamma > 0 else alpha

g, delta_over_sigma = 4, 1.0
patterns = {                                   # means in units of the range Delta
    "least favourable": np.array([0, 0.5, 0.5, 1]),
    "equally spaced": np.array([0, 1 / 3, 2 / 3, 1]),
    "most favourable": np.array([0, 0, 1, 1]),
}

def gamma_oneway(m, means):
    """Noncentrality m * sum (mu_i - mu_bar)^2 / sigma^2 for m observations per group."""
    dev = (means - means.mean()) * delta_over_sigma
    return m * np.sum(dev ** 2)

for name, means in patterns.items():
    m = 2
    while power(gamma_oneway(m, means), g - 1, g * (m - 1)) < 0.9:
        m += 1
    print(f"{name:17s} m = {m} per group, power {power(gamma_oneway(m, means), g - 1, g * (m - 1)):.3f}")
# <</oneway>>

m_needed = {}
for name, means in patterns.items():
    m = 2
    while power(gamma_oneway(m, means), g - 1, g * (m - 1)) < 0.9:
        m += 1
    m_needed[name] = m
    # the bounds of the least-favourable proposition
    assert gamma_oneway(1, means) >= 0.5 - 1e-12 and gamma_oneway(1, means) <= g / 4 + 1e-12
assert np.isclose(gamma_oneway(1, patterns["least favourable"]), 0.5)
assert np.isclose(gamma_oneway(1, patterns["most favourable"]), g / 4)
assert m_needed["least favourable"] > m_needed["equally spaced"] > m_needed["most favourable"]
# random configurations with range 1 never beat the bounds
for _ in range(2000):
    mid = rng.uniform(0, 1, g - 2)
    gm = gamma_oneway(1, np.concatenate([[0, 1], mid]))
    assert 0.5 - 1e-12 <= gm <= g / 4 + 1e-12

# simulation check of one power value: least favourable pattern at the chosen m
m = m_needed["least favourable"]
mu = np.repeat(patterns["least favourable"] * delta_over_sigma, m)
groups = np.repeat(np.arange(g), m)
reps = 20000
Y = mu + rng.standard_normal((reps, g * m))
means = np.stack([Y[:, groups == k].mean(axis=1) for k in range(g)], axis=1)
ss_between = m * np.sum((means - means.mean(axis=1, keepdims=True)) ** 2, axis=1)
ss_within = np.sum((Y - means[:, groups]) ** 2, axis=1)
Fsim = (ss_between / (g - 1)) / (ss_within / (g * (m - 1)))
rate = np.mean(Fsim > stats.f.ppf(0.95, g - 1, g * (m - 1)))
pw_theory = power(gamma_oneway(m, patterns["least favourable"]), g - 1, g * (m - 1))
assert abs(rate - pw_theory) < 4 * np.sqrt(pw_theory * (1 - pw_theory) / reps)

# ---- (2) the region test: sensitivity to a Northeast shift ---------------------------------
data = sm.datasets.statecrime.load_pandas().data
data.index = data.index.str.strip()
data = data.drop(index="District of Columbia")
codes = "SWWSWWNSSSWWMMMMSSNSNMMSMWMWNNWNSMMSWNNSMSSWNSWSMW"
region = np.array(list(codes))
y = data["murder"].to_numpy()
n = len(y)
X0 = np.column_stack([np.ones(n), data["poverty"], data["single"], data["urban"]])
D = np.column_stack([(region == c).astype(float) for c in "NSW"])
X = np.column_stack([X0, D])

# <<sensitivity>>
b0, *_ = np.linalg.lstsq(X0, D[:, 0], rcond=None)
k2 = np.sum((D[:, 0] - X0 @ b0) ** 2)      # ||(I - M0) d_N||^2: gamma = delta^2 k2 / sigma^2
bf, *_ = np.linalg.lstsq(X, y, rcond=None)
sigma = np.sqrt(np.sum((y - X @ bf) ** 2) / (n - 7))   # planning value for sigma
for delta in (0.5, 1.0, 1.5, 2.0):
    print(f"Northeast lower by {delta}: gamma = {delta**2 * k2 / sigma**2:.2f}, "
          f"power = {power(delta**2 * k2 / sigma**2, 3, n - 7):.3f}")
# <</sensitivity>>

pw_region = {dl: power(dl ** 2 * k2 / sigma ** 2, 3, n - 7) for dl in (0.5, 1.0, 1.5, 2.0)}
assert all(pw_region[a] < pw_region[b_] for a, b_ in ((0.5, 1.0), (1.0, 1.5), (1.5, 2.0)))
# the shift that gives power 0.8
lo, hi = 0.1, 10.0
for _ in range(80):
    mid = (lo + hi) / 2
    lo, hi = (mid, hi) if power(mid ** 2 * k2 / sigma ** 2, 3, n - 7) < 0.8 else (lo, mid)
delta80 = hi
assert abs(power(delta80 ** 2 * k2 / sigma ** 2, 3, n - 7) - 0.8) < 1e-6
# k2 is the Northeast count shrunk by what the covariates explain
assert 0 < k2 < np.sum(D[:, 0])

# ---- (3) observed power as a function of the p-value --------------------------------------
q, nu = 3, n - 7
pvals = np.array([0.5, 0.2, 0.1, 0.05, 0.01, 0.001])
Fobs = stats.f.isf(pvals, q, nu)
obs_power = np.array([power(q * f, q, nu) for f in Fobs])
assert np.all(np.diff(obs_power) > 0)          # smaller p-value, larger 'observed power'

gen = Generated("ch11", "power_design")
for name, key in (("least favourable", "lf"), ("equally spaced", "es"), ("most favourable", "mf")):
    gen.int(f"m_{key}", m_needed[name])
    mm = m_needed[name]
    gen.num(f"pow_{key}", power(gamma_oneway(mm, patterns[name]), g - 1, g * (mm - 1)), 3)
    gen.num(f"gamma1_{key}", gamma_oneway(1, patterns[name]), 3)
gen.num("sim_rate", rate, 3)
gen.num("sim_theory", pw_theory, 3)
gen.int("reps", reps)
gen.num("k2", k2, 2)
gen.int("n_north", int(np.sum(D[:, 0])))
gen.num("sigma", sigma, 3)
for dl, v in pw_region.items():
    gen.num(f"pow_region_{int(dl * 10)}", v, 3)
gen.num("delta80", delta80, 2)
gen.num("obs_power_05", obs_power[3], 3)
gen.num("obs_power_01", obs_power[4], 3)
gen.num("obs_power_20", obs_power[1], 3)
gen.write()

# ---- figure: power against group size ------------------------------------------------------
use_book_style()
fig, ax = plt.subplots(figsize=(4.6, 3.0))
ms = np.arange(2, 41)
colors = {"least favourable": COLORS["second"], "equally spaced": COLORS["third"],
          "most favourable": COLORS["accent"]}
for name, means in patterns.items():
    pw = [power(gamma_oneway(mm, means), g - 1, g * (mm - 1)) for mm in ms]
    ax.plot(ms, pw, color=colors[name], label=f"{name} (m = {m_needed[name]})")
    ax.plot(m_needed[name], power(gamma_oneway(m_needed[name], means), g - 1, g * (m_needed[name] - 1)),
            "o", color=colors[name], markersize=3.5)
ax.axhline(0.9, color=COLORS["grid"], linewidth=0.8, zorder=0)
ax.axhline(0.05, color=COLORS["grid"], linewidth=0.8, zorder=0)
ax.set_xlim(2, 40)
ax.set_ylim(0, 1)
ax.set_xlabel("observations per group, m")
ax.set_ylabel("power at level 0.05")
ax.legend(loc="lower right", bbox_to_anchor=(1.0, 0.1), frameon=False)
fig.tight_layout()
fig.savefig(figure_path("ch11", "power_group_size"))
