"""Chapter 13, Section 4: Bonferroni, Sidak and Holm for ten regression coefficients.

A regression with an intercept and ten mutually orthogonal regressors, n = 40, so the
ten t statistics are T_j = Z_j / sqrt(V / nu) with independent Z_j ~ N(delta_j, 1) and
V ~ chi^2(nu), nu = 29, independent of the Z_j (Chapter 7, thm-opt-sampling). k of the ten
coefficients are nonzero, each with delta = 3.5. All three procedures control the
familywise error rate; Holm finds more of the nonzero coefficients, and the more there
are, the larger its advantage.
"""
import matplotlib.pyplot as plt
import numpy as np
from scipy import integrate
from scipy import stats

from regbook import COLORS, Generated, figure_path, use_book_style

alpha, m, nu, delta, reps = 0.05, 10, 29, 3.5, 20_000

# <<simulate>>
def simulate(k, reps, seed):
    """FWER and average power of three procedures when k of m coefficients are nonzero."""
    rng = np.random.default_rng(seed)
    shift = np.r_[np.full(k, delta), np.zeros(m - k)]
    T = (rng.normal(size=(reps, m)) + shift) / np.sqrt(rng.chisquare(nu, (reps, 1)) / nu)
    p = 2 * stats.t.sf(np.abs(T), nu)
    rules = {
        "Bonferroni": p <= alpha / m,
        "Sidak": p <= 1 - (1 - alpha) ** (1 / m),
    }
    # Holm: sort each row, compare p_(j) with alpha / (m - j + 1), stop at the first failure
    order = np.argsort(p, axis=1)
    ps = np.take_along_axis(p, order, axis=1)
    passed = np.cumprod(ps <= alpha / (m - np.arange(m)), axis=1).astype(bool)
    holm = np.zeros_like(passed)
    np.put_along_axis(holm, order, passed, axis=1)
    rules["Holm"] = holm
    out = {}
    for name, rej in rules.items():
        fwer = np.mean(rej[:, k:].any(axis=1)) if k < m else np.nan
        power = rej[:, :k].mean() if k > 0 else np.nan
        out[name] = (fwer, power, rej)
    return out

res = simulate(5, reps, seed=5)
for name, (fwer, power, _) in res.items():
    print(f"{name:10s} FWER = {fwer:.3f}   average power = {power:.3f}")
# <</simulate>>

gen = Generated("ch13", "holm_power", prefix="hp")
gen.int("reps", reps)
gen.num("delta", delta, 1)
ks = np.arange(0, m + 1)
curves = {name: ([], []) for name in ("Bonferroni", "Sidak", "Holm")}
se = np.sqrt(alpha * (1 - alpha) / reps)
for k in ks:
    out = simulate(k, reps, seed=int(k))
    for name, (fwer, power, rej) in out.items():
        curves[name][0].append(fwer)
        curves[name][1].append(power)
        if k < m:
            assert fwer <= alpha + 4 * se             # strong control
    # Holm rejects everything Bonferroni rejects, data set by data set
    assert np.all(out["Holm"][2] >= out["Bonferroni"][2])
    if k in (1, 5, 9):
        for name, key in (("Bonferroni", "bon"), ("Sidak", "sid"), ("Holm", "holm")):
            gen.num(f"pow{key}{k}", out[name][1], 3)
    if k == 0:
        for name, key in (("Bonferroni", "bon"), ("Sidak", "sid"), ("Holm", "holm")):
            gen.num(f"fwer{key}0", out[name][0], 3)
# Holm's advantage grows with k
gain = np.array(curves["Holm"][1][1:]) - np.array(curves["Bonferroni"][1][1:])
assert gain[-1] > gain[0] and np.all(gain >= -0.002)

# Sidak's inequality with a common s: exact familywise coverage of the Sidak intervals
c_sid = stats.t.ppf((1 + (1 - alpha) ** (1 / m)) / 2, nu)
c_bon = stats.t.ppf(1 - alpha / (2 * m), nu)
chi = stats.chi(nu)                                   # S sqrt(nu) / sigma ~ chi(nu)


def coverage(c, power):
    f = lambda u: (2 * stats.norm.cdf(c * u / np.sqrt(nu)) - 1) ** power * chi.pdf(u)
    return integrate.quad(f, 0, np.inf)[0]


joint = coverage(c_sid, m)                            # P(all |T_j| <= c)
single = coverage(c_sid, 1)                           # P(|T_j| <= c)
assert np.isclose(single ** m, 1 - alpha)
assert joint > single ** m                            # strict: s is shared
joint_bon = coverage(c_bon, m)
assert joint_bon > joint > 1 - alpha
gen.num("csid", c_sid, 4)
gen.num("cbon", c_bon, 4)
gen.num("jointsid", joint, 4)
gen.num("jointbon", joint_bon, 4)
gen.write()

use_book_style()
fig, axes = plt.subplots(1, 2, figsize=(5.8, 2.4))
styles = {"Bonferroni": (COLORS["second"], "-", "o"), "Sidak": (COLORS["third"], "--", "s"),
          "Holm": (COLORS["accent"], "-", "^")}
labels = {"Bonferroni": "Bonferroni", "Sidak": "Sidak", "Holm": "Holm"}
for name, (fw, pw) in curves.items():
    col, ls, mk = styles[name]
    axes[0].plot(ks[:-1], fw[:-1], color=col, linestyle=ls, marker=mk, markersize=3,
                 label=labels[name])
    axes[1].plot(ks[1:], pw[1:], color=col, linestyle=ls, marker=mk, markersize=3,
                 label=labels[name])
axes[0].axhline(alpha, color=COLORS["grid"], linewidth=0.8, zorder=0)
axes[0].set_ylim(0, 0.08)
axes[0].set_xlabel("number of nonzero coefficients $k$")
axes[0].set_ylabel("familywise error rate")
axes[0].set_title("(a) error rate")
axes[1].set_xlabel("number of nonzero coefficients $k$")
axes[1].set_ylabel("average power")
axes[1].set_title("(b) power")
axes[1].legend(frameon=False, loc="upper left")
fig.tight_layout()
fig.savefig(figure_path("ch13", "holm_power"))
