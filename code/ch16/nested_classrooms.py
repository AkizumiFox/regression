"""Chapter 16, Section 5: a balanced nested layout.

Synthetic data: test scores of pupils taught under three reading curricula (factor A), in
four classrooms per curriculum (factor B, nested in A: each classroom follows one curriculum),
with six pupils per classroom, 72 pupils in all. Classroom effects are generated at random,
as they would be in a real study. The script
(1) computes the nested decomposition from Kronecker projections and checks it against
    statsmodels and against the crossed decomposition (B(A) = B + AB);
(2) tests classrooms within curricula and curricula against the pupil-level error, as the
    fixed-effects theory of Section 16.5 prescribes, and against the classroom mean square,
    the ratio used when classrooms are regarded as random (previewed only);
(3) checks by simulation, with random classroom effects and no curriculum effect, the size
    of both tests;
(4) plots the pupils' scores by classroom.
"""
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from scipy import stats
from statsmodels.stats.anova import anova_lm

from regbook import COLORS, Generated, figure_path, use_book_style

a, b, m = 3, 4, 6
curricula = ["phonics", "whole language", "mixed"]
data_rng = np.random.default_rng(16_05)
curriculum_effect = np.array([3.0, -2.0, 0.0])
classroom_effect = data_rng.normal(scale=4.0, size=(a, b))
scores = np.round(62 + curriculum_effect[:, None, None] + classroom_effect[:, :, None]
                  + data_rng.normal(scale=8.0, size=(a, b, m)))


# <<nested>>
def Jbar(k):
    return np.full((k, k), 1.0 / k)


def Cen(k):
    return np.eye(k) - Jbar(k)


y = scores.ravel()                                    # pupils fastest, then classroom, then curriculum
P = {"curriculum": np.kron(np.kron(Cen(a), Jbar(b)), Jbar(m)),
     "classroom(curriculum)": np.kron(np.kron(np.eye(a), Cen(b)), Jbar(m)),
     "error": np.kron(np.kron(np.eye(a), np.eye(b)), Cen(m))}
ss = {k: y @ Pk @ y for k, Pk in P.items()}
df = {k: int(round(np.trace(Pk))) for k, Pk in P.items()}
ms = {k: ss[k] / df[k] for k in P}
for k in P:
    print(f"{k:22s} df {df[k]:2d}  SS {ss[k]:8.1f}  MS {ms[k]:7.1f}")
F_cls = ms["classroom(curriculum)"] / ms["error"]
F_cur_fixed = ms["curriculum"] / ms["error"]
F_cur_random = ms["curriculum"] / ms["classroom(curriculum)"]
print(f"classrooms within curricula: F = {F_cls:.2f}, p = {stats.f.sf(F_cls, df['classroom(curriculum)'], df['error']):.3f}")
print(f"curricula against pupils:    F = {F_cur_fixed:.2f}, p = {stats.f.sf(F_cur_fixed, 2, df['error']):.3f}")
print(f"curricula against classrooms: F = {F_cur_random:.2f}, p = {stats.f.sf(F_cur_random, 2, df['classroom(curriculum)']):.3f}")
# <</nested>>

p_cls = stats.f.sf(F_cls, df["classroom(curriculum)"], df["error"])
p_fixed = stats.f.sf(F_cur_fixed, 2, df["error"])
p_random = stats.f.sf(F_cur_random, 2, df["classroom(curriculum)"])
assert df == {"curriculum": 2, "classroom(curriculum)": 9, "error": 60}
assert p_fixed < 0.05 < p_random

# statsmodels with classroom labels 1..4 reused within curricula: the nested term is the
# classroom-within-curriculum interaction; the crossed analysis puts the same SS in B + AB
frame = pd.DataFrame({"score": y, "cur": np.repeat(np.arange(a), b * m),
                      "cls": np.tile(np.repeat(np.arange(b), m), a)})
tab = anova_lm(smf.ols("score ~ C(cur) + C(cur):C(cls)", frame).fit())
assert np.isclose(tab.loc["C(cur):C(cls)", "sum_sq"], ss["classroom(curriculum)"])
crossed = anova_lm(smf.ols("score ~ C(cur) * C(cls)", frame).fit())
assert np.isclose(crossed.loc["C(cls)", "sum_sq"] + crossed.loc["C(cur):C(cls)", "sum_sq"],
                  ss["classroom(curriculum)"])
cls_means = scores.mean(axis=2)
cur_means = cls_means.mean(axis=1)

# ---- size of the two curriculum tests when classroom effects are random ------------
# <<size>>
rng = np.random.default_rng(20260919)
reps = 20_000
Pc, Pb, Pe = P["curriculum"], P["classroom(curriculum)"], P["error"]
rej_fixed = rej_random = 0
for _ in range(reps):
    sim = (rng.normal(scale=4.0, size=(a, b))[:, :, None] + rng.normal(scale=8.0, size=(a, b, m))).ravel()
    msc, msb, mse = sim @ Pc @ sim / 2, sim @ Pb @ sim / 9, sim @ Pe @ sim / 60
    rej_fixed += msc / mse > stats.f.ppf(0.95, 2, 60)
    rej_random += msc / msb > stats.f.ppf(0.95, 2, 9)
print(f"no curriculum effect, random classrooms: rejection rate {rej_fixed / reps:.3f} "
      f"against pupils, {rej_random / reps:.3f} against classrooms")
# <</size>>
size_fixed, size_random = rej_fixed / reps, rej_random / reps
assert size_fixed > 0.2 and abs(size_random - 0.05) < 4 * np.sqrt(0.05 * 0.95 / reps)
# expected mean squares with random classrooms (preview): E MS_A = sigma^2 + m sigma_B^2 + ...,
# E MS_B(A) = sigma^2 + m sigma_B^2, so the ratio of expectations under H_A is
ems_ratio = (64 + m * 16) / 64

gen = Generated("ch16", "nested_classrooms")
for k in P:
    key = {"curriculum": "A", "classroom(curriculum)": "BA", "error": "E"}[k]
    gen.num(f"ss_{key}", ss[k], 1)
    gen.num(f"ms_{key}", ms[k], 1)
gen.num("F_cls", F_cls, 2)
gen.num("p_cls", p_cls, 3)
gen.num("F_fixed", F_cur_fixed, 2)
gen.num("p_fixed", p_fixed, 4)
gen.num("F_random", F_cur_random, 2)
gen.num("p_random", p_random, 3)
for i in range(a):
    gen.num(f"cur{i + 1}", cur_means[i], 1)
    for j in range(b):
        gen.num(f"cls{i + 1}{j + 1}", cls_means[i, j], 1)
gen.num("size_fixed", size_fixed, 3)
gen.num("size_random", size_random, 3)
gen.num("ems_ratio", ems_ratio, 1)
gen.int("reps", reps)
gen.num("f_crit_fixed", stats.f.ppf(0.95, 2, 60), 2)
gen.num("f_crit_random", stats.f.ppf(0.95, 2, 9), 2)
gen.write()

use_book_style()
fig, ax = plt.subplots(figsize=(5.2, 2.5))
colours = [COLORS["accent"], COLORS["second"], COLORS["third"]]
pos = 0
ticks, labels = [], []
for i in range(a):
    for j in range(b):
        ax.scatter(np.full(m, pos) + np.linspace(-0.18, 0.18, m), scores[i, j], s=8, color=colours[i],
                   alpha=0.6, linewidths=0)
        ax.plot([pos - 0.3, pos + 0.3], [cls_means[i, j]] * 2, color=colours[i], linewidth=1.4)
        ticks.append(pos)
        labels.append(f"{j + 1}")
        pos += 1
    ax.plot([pos - b - 0.4, pos - 0.6], [cur_means[i]] * 2, color=colours[i], linestyle="--", linewidth=0.8)
    ax.text(pos - b / 2 - 0.5, scores.max() + 3, curricula[i], ha="center", fontsize=8, color=colours[i])
    pos += 0.8
ax.set_xticks(ticks, labels)
ax.set_xlabel("classroom within curriculum")
ax.set_ylabel("score")
ax.set_ylim(scores.min() - 4, scores.max() + 8)
fig.tight_layout()
fig.savefig(figure_path("ch16", "nested_classrooms"))
