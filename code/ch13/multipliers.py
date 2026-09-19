"""Chapter 13, Section 3: half-width multipliers for all pairwise comparisons of g means.

Balanced one-way layout, 95% simultaneous intervals for the g(g-1)/2 differences,
half-width = multiplier x standard error of a difference. Also: Scheffe's s F(s, nu)
increases with s, and Tukey's bound loses to Scheffe's for contrasts spread over many means.
"""
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats

from regbook import COLORS, Generated, figure_path, use_book_style

alpha, nu = 0.05, 30

# <<multipliers>>
def pairwise_multipliers(g, nu, alpha=0.05):
    m = g * (g - 1) // 2                              # number of pairs
    return {
        "unadjusted": stats.t.ppf(1 - alpha / 2, nu),
        "Tukey": stats.studentized_range.ppf(1 - alpha, g, nu) / np.sqrt(2),
        "Sidak": stats.t.ppf((1 + (1 - alpha) ** (1 / m)) / 2, nu),
        "Bonferroni": stats.t.ppf(1 - alpha / (2 * m), nu),
        "Scheffe": np.sqrt((g - 1) * stats.f.ppf(1 - alpha, g - 1, nu)),
    }

for g in (2, 3, 5, 10):
    row = pairwise_multipliers(g, nu)
    print(f"g = {g:2d}: " + "  ".join(f"{k} {v:.3f}" for k, v in row.items()))
# <</multipliers>>

gs = np.arange(2, 16)
table = {k: [] for k in pairwise_multipliers(3, nu)}
for g in gs:
    row = pairwise_multipliers(g, nu)
    for k, v in row.items():
        table[k].append(v)
    if g == 2:                                        # one comparison: all coincide
        assert np.allclose(list(row.values()), row["unadjusted"])
    else:
        assert (row["unadjusted"] < row["Tukey"] < row["Sidak"] < row["Bonferroni"]
                < row["Scheffe"])
for g in gs:                                          # also with nu = infinity (large)
    r = pairwise_multipliers(g, 10 ** 6)
    if g > 2:
        assert r["Tukey"] < r["Sidak"] < r["Bonferroni"] < r["Scheffe"]

# s F_alpha(s, nu) increases with s
sF = [s * stats.f.ppf(1 - alpha, s, nu) for s in range(1, 30)]
assert np.all(np.diff(sF) > 0)

# Tukey versus Scheffe for a contrast spread over g = 6 means (balanced, per unit s/sqrt(m))
g6 = 6
c = np.array([1, 1, 1, -1, -1, -1]) / 3
tukey_hw = stats.studentized_range.ppf(1 - alpha, g6, nu) * 0.5 * np.abs(c).sum()
scheffe_hw = np.sqrt((g6 - 1) * stats.f.ppf(1 - alpha, g6 - 1, nu)) * np.sqrt(np.sum(c ** 2))
pair = np.array([1, -1, 0, 0, 0, 0])
tukey_pair = stats.studentized_range.ppf(1 - alpha, g6, nu) * 0.5 * np.abs(pair).sum()
scheffe_pair = np.sqrt((g6 - 1) * stats.f.ppf(1 - alpha, g6 - 1, nu)) * np.sqrt(2)
assert scheffe_hw < tukey_hw and tukey_pair < scheffe_pair

gen = Generated("ch13", "multipliers", prefix="mult")
gen.int("nu", nu)
for g in (3, 5, 10):
    row = pairwise_multipliers(g, nu)
    for k, key in (("unadjusted", "t"), ("Tukey", "tuk"), ("Sidak", "sid"),
                   ("Bonferroni", "bon"), ("Scheffe", "sch")):
        gen.num(f"{key}{g}", row[k], 3)
gen.num("tukspread", tukey_hw, 3)
gen.num("schspread", scheffe_hw, 3)
gen.num("tukpair", tukey_pair, 3)
gen.num("schpair", scheffe_pair, 3)
gen.write()

use_book_style()
fig, ax = plt.subplots(figsize=(4.6, 2.7))
styles = {"unadjusted": (COLORS["muted"], ":"), "Tukey": (COLORS["accent"], "-"),
          "Sidak": (COLORS["third"], "--"), "Bonferroni": (COLORS["second"], "-"),
          "Scheffe": (COLORS["thread"], "-.")}
labels = {"unadjusted": "unadjusted $t$", "Tukey": "Tukey", "Sidak": "Sidak",
          "Bonferroni": "Bonferroni", "Scheffe": "Scheffé"}
for k, v in table.items():
    col, ls = styles[k]
    ax.plot(gs, v, color=col, linestyle=ls, marker="o", markersize=2.5, label=labels[k])
ax.set_xlabel("number of groups $g$")
ax.set_ylabel("multiplier of the standard error")
ax.set_xticks(gs[::2])
ax.legend(frameon=False, loc="upper left")
fig.savefig(figure_path("ch13", "multipliers"))
