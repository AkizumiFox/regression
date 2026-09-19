"""Chapter 16, Section 4: an unreplicated 2^4 factorial.

Synthetic data: peel strength (N) of heat seals made at two levels of each of four factors,
A = jaw temperature, B = jaw pressure, C = dwell time, D = film gauge, one seal per run
(16 runs). The script
(1) builds the +-1 contrast columns of all 15 effects and checks that they are orthogonal
    and that the effect estimates reproduce the Kronecker decomposition of Section 16.4;
(2) estimates the effects and their sums of squares, n * effect^2 / 4;
(3) tests main effects and two-factor interactions against the five higher-order
    interactions, pooled as error by a decision made before looking at the data;
(4) draws the half-normal plot of the 15 absolute effects.
"""
import itertools

import matplotlib.pyplot as plt
import numpy as np
from scipy import stats
from statsmodels.stats.multitest import multipletests

from regbook import COLORS, Generated, figure_path, use_book_style

# <<design>>
k = 4
runs = np.array(list(itertools.product([-1, 1], repeat=k)))[:, ::-1]   # standard order: A changes fastest
names = []
columns = []
for size in range(1, k + 1):
    for S in itertools.combinations(range(k), size):
        names.append("".join("ABCD"[i] for i in S))
        columns.append(np.prod(runs[:, list(S)], axis=1))
Xc = np.column_stack(columns)                   # 16 x 15 matrix of +-1 contrast columns
print("columns orthogonal:", np.allclose(Xc.T @ Xc, 16 * np.eye(15)))
# <</design>>
assert np.allclose(Xc.T @ Xc, 16 * np.eye(15)) and np.allclose(Xc.sum(0), 0)

# synthetic responses
true = 30 + 3.0 * runs[:, 0] + 2.0 * runs[:, 2] + 1.5 * runs[:, 0] * runs[:, 2] - 1.25 * runs[:, 3]
rng = np.random.default_rng(16_04)
y = np.round(true + rng.normal(scale=1.0, size=16), 1)

# <<effects>>
n = len(y)
effects = Xc.T @ y / (n / 2)                    # mean at + minus mean at -
ss = n * effects ** 2 / 4                       # one degree of freedom each
order = np.argsort(-np.abs(effects))
for i in order[:6]:
    print(f"{names[i]:4s} effect {effects[i]:6.2f}   SS {ss[i]:7.2f}")
high = [i for i, nm in enumerate(names) if len(nm) >= 3]   # ABC, ABD, ACD, BCD, ABCD
ss_pooled = ss[high].sum()
s2 = ss_pooled / len(high)
F = ss / s2
p = stats.f.sf(F, 1, len(high))
for i in order[:6]:
    print(f"{names[i]:4s} F = {F[i]:7.2f}  p = {p[i]:.3f}")
# <</effects>>

# the effects are the coordinates of y on the orthogonal contrast columns, and each SS is the
# squared length of the projection of y on its column (the 2^k case of the Kronecker decomposition)
total = sum(np.outer(c, c) / 16 for c in columns)
assert np.allclose(total, np.eye(16) - np.full((16, 16), 1 / 16))
for c, s_ in zip(columns, ss):
    assert np.isclose((c @ y) ** 2 / 16, s_)
assert np.isclose(ss.sum(), ((y - y.mean()) ** 2).sum())
# the Kronecker projection for the AC interaction equals the projection on its contrast column
C2, J2 = np.eye(2) - np.full((2, 2), 0.5), np.full((2, 2), 0.5)
P_AC = np.ones((1, 1))
for f in "DCBA":                                # A changes fastest, so it is the last Kronecker factor
    P_AC = np.kron(P_AC, C2 if f in "AC" else J2)
c_AC = columns[names.index("AC")]
assert np.allclose(P_AC, np.outer(c_AC, c_AC) / 16)

active = [names[i] for i in range(15) if p[i] < 0.01]
assert set(active) == {"A", "C", "AC", "D"}
# Holm's step-down at familywise level 0.05 over the ten effects tested against the pooled error
tested = [i for i in range(15) if i not in high]
holm_reject, _, _, _ = multipletests(p[tested], alpha=0.05, method="holm")
assert {names[i] for i, r in zip(tested, holm_reject) if r} == {"A", "C", "AC", "D"}
assert len(tested) == 10 and p[names.index("BD")] > 0.05 / 6

gen = Generated("ch16", "two_level")
for i, nm in enumerate(names):
    gen.num(f"eff_{nm}", effects[i], 2)
    gen.num(f"ss_{nm}", ss[i], 2)
    gen.num(f"F_{nm}", F[i], 1)
    gen.num(f"p_{nm}", p[i], 3)
for r in range(16):
    gen.num(f"y{r + 1}", y[r], 1)
gen.num("s2", s2, 3)
gen.num("mean", y.mean(), 2)
gen.num("f_crit", stats.f.ppf(0.95, 1, 5), 2)
gen.write()

# ---- half-normal plot -----------------------------------------------------------------
use_book_style()
fig, ax = plt.subplots(figsize=(3.6, 2.6))
abseff = np.abs(effects)
idx = np.argsort(abseff)
quant = stats.norm.ppf(0.5 + 0.5 * (np.arange(1, 16) - 0.5) / 15)
ax.scatter(quant, abseff[idx], s=12, color=COLORS["accent"], linewidths=0, zorder=3)
for qv, i in zip(quant, idx):
    if abseff[i] > 1.2:
        ax.annotate(names[i], (qv, abseff[i]), xytext=(-12, 2), textcoords="offset points", fontsize=7)
small = idx[:10]
slope = np.median(abseff[small] / quant[:10])
xs = np.linspace(0, quant.max(), 2)
ax.plot(xs, slope * xs, color=COLORS["muted"], linewidth=0.8, linestyle="--")
ax.set_xlabel("half-normal quantile")
ax.set_ylabel("absolute effect")
fig.tight_layout()
fig.savefig(figure_path("ch16", "half_normal"))
