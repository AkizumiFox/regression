"""Chapter 15, Section 2: estimable contrasts in the nitrogen trial.

Planned contrasts with t intervals, the comparisons with the control (Bonferroni, and
Dunnett's exact method for reference), Tukey's intervals for all pairs, a Scheffe interval for
a contrast suggested by the data, and the power of a planned single contrast against the
omnibus F test. Data as in nitrogen_trial.py.
"""
import itertools

import numpy as np
from scipy import stats

from regbook import Generated

# <<data>>
rates = np.array([0, 40, 80, 120, 160])            # kg N per hectare
g, m = len(rates), 6                                # 5 rates, 6 plots per rate
true_means = 6.8 - 2.7 * np.exp(-0.018 * rates)     # unknown to the analyst
rng = np.random.default_rng(1515)
group = np.repeat(np.arange(g), m)                  # plot i received rate group[i]
y = np.round(true_means[group] + rng.normal(0, 0.5, g * m), 2)   # yield, t/ha
n = len(y)
# <</data>>

# <<contrast>>
n_k = np.bincount(group).astype(float)
means = np.bincount(group, weights=y) / n_k
nu = n - g
s = np.sqrt(np.sum((y - means[group]) ** 2) / nu)

def contrast(c, psi0=0.0, mult=None):
    """Estimate, standard error, t statistic, SS and interval for sum_k c_k mu_k."""
    est = c @ means
    v = np.sum(c ** 2 / n_k)                        # Var(estimate) / sigma^2
    se = s * np.sqrt(v)
    mult = stats.t.ppf(0.975, nu) if mult is None else mult
    return est, se, (est - psi0) / se, est ** 2 / v, (est - mult * se, est + mult * se)

c_fert = np.array([-1, 0.25, 0.25, 0.25, 0.25])     # fertilized (average) minus control
c_top = np.array([0, 0, -1, 0, 1])                  # 160 against 80 kg N/ha
for name, c in (("fertilized - control", c_fert), ("160 - 80", c_top)):
    est, se, t, ss, (lo, hi) = contrast(c)
    print(f"{name:21s} {est:6.3f}  se {se:.3f}  t = {t:5.2f}  SS = {ss:6.3f}"
          f"  95% interval [{lo:.3f}, {hi:.3f}]")
# <</contrast>>

# checks: SS equals the squared length of the projection onto u_c; F = t^2
Z = (group[:, None] == np.arange(g)).astype(float)
for c in (c_fert, c_top):
    u = Z @ (c / n_k)
    est, se, t, ss, _ = contrast(c)
    assert np.isclose(u @ y, est)
    assert np.isclose((u @ y) ** 2 / (u @ u), ss)
    assert np.isclose(t ** 2, ss / s ** 2)
    assert np.isclose(abs(u @ np.ones(n)), 0)
    # scaling a contrast changes neither t nor SS
    assert np.isclose(contrast(3 * c)[3], ss) and np.isclose(contrast(3 * c)[2], t)
est_f, se_f, t_f, ss_f, ci_f = contrast(c_fert)
est_t, se_t, t_t, ss_t, ci_t = contrast(c_top)
p_top = 2 * stats.t.sf(abs(t_t), nu)

# <<control>>
k_ctrl = g - 1                                      # four comparisons with the control
bonf = stats.t.ppf(1 - 0.05 / (2 * k_ctrl), nu)     # Bonferroni multiplier
se_pair = s * np.sqrt(2 / m)
for k in range(1, g):
    d = means[k] - means[0]
    print(f"{rates[k]:3d} - 0: {d:6.3f}  Bonferroni interval"
          f" [{d - bonf * se_pair:.3f}, {d + bonf * se_pair:.3f}]")
# <</control>>

dun = stats.dunnett(*[y[group == k] for k in range(1, g)], control=y[group == 0], random_state=15)
ci_dun = dun.confidence_interval(0.95)
dunnett_mult = (ci_dun.high[0] - (means[1] - means[0])) / se_pair
assert dunnett_mult < bonf
t_unadj = stats.t.ppf(0.975, nu)

# <<tukey>>
q = stats.studentized_range.ppf(0.95, g, nu)        # q_{0.05}(g, nu)
hsd = q * s / np.sqrt(m)                            # Tukey half-width, all pairs
pairs = list(itertools.combinations(range(g), 2))
different = [(int(rates[k]), int(rates[l])) for k, l in pairs if abs(means[k] - means[l]) > hsd]
print(f"Tukey half-width {hsd:.3f}; pairs declared different: {len(different)} of {len(pairs)}")
print("not declared different:",
      [(int(rates[k]), int(rates[l])) for k, l in pairs if abs(means[k] - means[l]) <= hsd])
# <</tukey>>

from statsmodels.stats.multicomp import pairwise_tukeyhsd  # noqa: E402
res = pairwise_tukeyhsd(y, rates[group])
assert res.reject.sum() == len(different)
same = [(rates[k], rates[l]) for k, l in pairs if abs(means[k] - means[l]) <= hsd]
assert len(same) == 4
assert np.isclose(hsd / (s * np.sqrt(2 / m)), q / np.sqrt(2))

# Scheffe interval for a contrast chosen after looking: the two highest rates vs 80
c_post = np.array([0, 0, -1, 0.5, 0.5])
sch = np.sqrt((g - 1) * stats.f.ppf(0.95, g - 1, nu))
est_p, se_p, t_p, ss_p, ci_p = contrast(c_post, mult=sch)
_, _, _, _, ci_p_naive = contrast(c_post)

# <<power>>
def power_f(gamma, q, nu, alpha=0.05):
    return stats.ncf.sf(stats.f.ppf(1 - alpha, q, nu), q, nu, gamma)

for gamma in (4.0, 8.0, 12.0):                      # all of it along the planned contrast
    print(f"gamma = {gamma:4.1f}: contrast test {power_f(gamma, 1, nu):.3f},"
          f" omnibus F test {power_f(gamma, g - 1, nu):.3f}")
# <</power>>

gen = Generated("ch15", "nitrogen_contrasts")
gen.num("est_fert", est_f, 3)
gen.num("se_fert", se_f, 3)
gen.num("t_fert", t_f, 2)
gen.num("lo_fert", ci_f[0], 3)
gen.num("hi_fert", ci_f[1], 3)
gen.num("ss_fert", ss_f, 3)
gen.num("est_top", est_t, 3)
gen.num("se_top", se_t, 3)
gen.num("t_top", t_t, 2)
gen.num("p_top", p_top, 3)
gen.num("lo_top", ci_t[0], 3)
gen.num("hi_top", ci_t[1], 3)
gen.num("ss_top", ss_t, 3)
gen.num("t_unadj", t_unadj, 3)
gen.num("bonf", bonf, 3)
gen.num("dunnett", dunnett_mult, 3)
gen.num("se_pair", se_pair, 4)
for k in range(1, g):
    d = means[k] - means[0]
    gen.num(f"diff{k}", d, 3)
    gen.num(f"diff{k}_lo", d - bonf * se_pair, 3)
    gen.num(f"diff{k}_hi", d + bonf * se_pair, 3)
gen.num("q", q, 3)
gen.num("hsd", hsd, 3)
gen.int("n_different", len(different))
gen.num("sch", sch, 3)
gen.num("est_post", est_p, 3)
gen.num("se_post", se_p, 3)
gen.num("lo_post", ci_p[0], 3)
gen.num("hi_post", ci_p[1], 3)
gen.num("lo_post_naive", ci_p_naive[0], 3)
gen.num("hi_post_naive", ci_p_naive[1], 3)
for gamma in (4.0, 8.0, 12.0):
    gen.num(f"pow1_{int(gamma)}", power_f(gamma, 1, nu), 3)
    gen.num(f"pow4_{int(gamma)}", power_f(gamma, g - 1, nu), 3)
    assert power_f(gamma, 1, nu) > power_f(gamma, g - 1, nu)
gen.write()
