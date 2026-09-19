"""Chapter 8, Section 6: two factors, empty cells, connectedness, and a factor with a covariate.

Part 1 (synthetic layout): in the additive two-way model the rank is a + b minus the number
of connected components of the design, and a row difference is estimable iff the two rows are
connected. Part 2 (1996 ANES, public domain): party identification by education has empty
cells, so the interaction model has rank equal to the number of occupied cells and some
interaction contrasts are not estimable. Part 3 (1996 ANES): TV news watching by education,
unadjusted and adjusted for age, and with separate age slopes for three education groups.
"""
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import statsmodels.api as sm

from regbook import COLORS, Generated, figure_path, use_book_style


# <<additive>>
import numpy as np

def additive_matrix(cells, a, b):
    """Model matrix [1, row indicators, column indicators] for occupied (i, j) cells."""
    rows = np.array([i for i, j in cells])
    cols = np.array([j for i, j in cells])
    return np.column_stack([np.ones(len(cells)), np.eye(a)[rows], np.eye(b)[cols]])

def estimable(lam, X):
    """lambda is estimable iff appending it to the rows of X does not raise the rank."""
    return np.linalg.matrix_rank(np.vstack([X, lam])) == np.linalg.matrix_rank(X)

a, b = 4, 5
cells = [(0, 0), (0, 1), (1, 1), (1, 2), (2, 3), (3, 3), (3, 4)]   # two separate blocks
X = additive_matrix(cells, a, b)
row_diff = lambda i, k: np.r_[0, np.eye(a)[i] - np.eye(a)[k], np.zeros(b)]
print("rank:", np.linalg.matrix_rank(X), " (a + b - 1 =", a + b - 1, ")")
print("alpha_1 - alpha_2 estimable:", estimable(row_diff(0, 1), X))
print("alpha_1 - alpha_3 estimable:", estimable(row_diff(0, 2), X))
X_linked = additive_matrix(cells + [(2, 2)], a, b)                # one more occupied cell
print("after adding cell (3, 3): rank", np.linalg.matrix_rank(X_linked),
      " alpha_1 - alpha_3 estimable:", estimable(row_diff(0, 2), X_linked))
# <</additive>>


def components(cells, a, b):
    """Connected components of the bipartite row-column graph (rows 0..a-1, columns a..a+b-1)."""
    parent = list(range(a + b))
    def find(u):
        while parent[u] != u:
            parent[u] = parent[parent[u]]
            u = parent[u]
        return u
    for i, j in cells:
        parent[find(i)] = find(a + j)
    return len({find(u) for u in range(a + b)})


assert components(cells, a, b) == 2
assert np.linalg.matrix_rank(X) == a + b - components(cells, a, b) == 7
assert np.linalg.matrix_rank(X_linked) == a + b - 1 == 8
assert estimable(row_diff(0, 1), X) and not estimable(row_diff(0, 2), X)
assert estimable(row_diff(0, 2), X_linked)
# random designs: rank = a + b - (number of components) always
rng = np.random.default_rng(86)
for _ in range(200):
    occ = [(i, j) for i in range(a) for j in range(b) if rng.random() < 0.3]
    used_r = {i for i, _ in occ}
    used_c = {j for _, j in occ}
    if len(used_r) < a or len(used_c) < b:
        continue
    assert np.linalg.matrix_rank(additive_matrix(occ, a, b)) == a + b - components(occ, a, b)

gen = Generated("ch08", "two_factors", prefix="tf")

# ---- Part 2: empty cells in real data -------------------------------------------
# <<emptycells>>
import statsmodels.api as sm

anes = sm.datasets.anes96.load_pandas().data
pid = anes["PID"].to_numpy().astype(int)              # 7 levels
educ = anes["educ"].to_numpy().astype(int) - 1        # 7 levels, now 0..6
counts = np.zeros((7, 7), dtype=int)
np.add.at(counts, (pid, educ), 1)
print("empty (PID, educ) cells:",
      [(int(i), int(j) + 1) for i, j in zip(*np.nonzero(counts == 0))])

n = len(pid)
cell = pid * 7 + educ
X_int = np.column_stack([np.ones(n), np.eye(7)[pid], np.eye(7)[educ], np.eye(49)[cell]])
print("interaction model: p =", X_int.shape[1], " rank =", np.linalg.matrix_rank(X_int),
      " occupied cells =", int((counts > 0).sum()))

def interaction_contrast(i, k, j, l):
    """gamma_ij - gamma_il - gamma_kj + gamma_kl, as a coefficient vector on X_int."""
    lam = np.zeros(X_int.shape[1])
    for (r, c), s in [((i, j), 1), ((i, l), -1), ((k, j), -1), ((k, l), 1)]:
        lam[15 + r * 7 + c] += s
    return lam

print("PID 0 vs 3, educ 2 vs 3:", estimable(interaction_contrast(0, 3, 1, 2), X_int))
print("PID 0 vs 3, educ 1 vs 3:", estimable(interaction_contrast(0, 3, 0, 2), X_int))
# <</emptycells>>

empty = list(zip(*np.nonzero(counts == 0)))
assert empty == [(3, 0), (5, 0)]
assert np.linalg.matrix_rank(X_int) == (counts > 0).sum() == 47
assert estimable(interaction_contrast(0, 3, 1, 2), X_int)
assert not estimable(interaction_contrast(0, 3, 0, 2), X_int)
X_add = X_int[:, :15]
assert np.linalg.matrix_rank(X_add) == 13
# in the additive model every row difference is estimable (the design is connected)
assert estimable(np.r_[0, np.eye(7)[3] - np.eye(7)[0], np.zeros(56)][:15], X_add)
gen.int("anes:p", X_int.shape[1])
gen.int("anes:rank", np.linalg.matrix_rank(X_int))
gen.int("anes:n30", counts[3, 1])
gen.int("anes:n00", counts[0, 0])

# ---- Part 3: a factor with a covariate -----------------------------------------
# <<covariate>>
tv = anes["TVnews"].to_numpy()                       # times per week watching TV news
age = anes["age"].to_numpy()
E = np.eye(7)[educ][:, [0, 1, 3, 4, 5, 6]]          # reference: high-school graduate
unadjusted = np.linalg.lstsq(np.column_stack([np.ones(n), E]), tv, rcond=None)[0]
adjusted = np.linalg.lstsq(np.column_stack([np.ones(n), E, age]), tv, rcond=None)[0]
print("some high school minus high-school graduate:",
      f"unadjusted {unadjusted[2]:.3f}, adjusted for age {adjusted[2]:.3f}")
print("mean age: some high school", age[educ == 1].mean().round(1),
      " high-school graduate", age[educ == 2].mean().round(1),
      " common age slope", adjusted[-1].round(4))
# <</covariate>>

m1, m2 = tv[educ == 1].mean(), tv[educ == 2].mean()
assert np.isclose(unadjusted[2], m1 - m2)
# adjusted difference = unadjusted difference - slope * (difference in mean ages)
slope = adjusted[-1]
assert np.isclose(adjusted[2], (m1 - m2) - slope * (age[educ == 1].mean() - age[educ == 2].mean()))
assert adjusted[2] < 0.5 * unadjusted[2]
gen.num("cov:unadj", unadjusted[2], 3)
gen.num("cov:adj", adjusted[2], 3)
gen.num("cov:slope", slope, 4)
gen.num("cov:age1", age[educ == 1].mean(), 1)
gen.num("cov:age2", age[educ == 2].mean(), 1)

# <<slopes>>
grp = np.select([educ <= 2, educ <= 4], [0, 1], 2)  # education codes 1-3 / 4-5 / 6-7
D = np.eye(3)[grp]
X_par = np.column_stack([D, age])                   # parallel lines: 3 intercepts, one slope
X_sep = np.column_stack([D, D * age[:, None]])      # separate lines: 3 intercepts, 3 slopes
b_par = np.linalg.lstsq(X_par, tv, rcond=None)[0]
b_sep = np.linalg.lstsq(X_sep, tv, rcond=None)[0]
sse_par = np.sum((tv - X_par @ b_par) ** 2)
sse_sep = np.sum((tv - X_sep @ b_sep) ** 2)
print("separate slopes:", b_sep[3:].round(4),
      f" SSE parallel {sse_par:.1f}, separate {sse_sep:.1f}")
# <</slopes>>

assert np.linalg.matrix_rank(X_sep) == 6 and sse_sep <= sse_par
# the separate-slopes fit is the same as three simple regressions
for k in range(3):
    fk = np.polyfit(age[grp == k], tv[grp == k], 1)
    assert np.isclose(fk[0], b_sep[3 + k]) and np.isclose(fk[1], b_sep[k])
# same column space as the reference-coded interaction model
X_ref = np.column_stack([np.ones(n), D[:, 1:], age, D[:, 1:] * age[:, None]])
assert np.linalg.matrix_rank(np.column_stack([X_ref, X_sep])) == 6
b_ref = np.linalg.lstsq(X_ref, tv, rcond=None)[0]
assert np.allclose(b_ref[4:], b_sep[4:] - b_sep[3])
for k in range(3):
    gen.num(f"sl:b{k}", b_sep[3 + k], 4)
    gen.num(f"sl:a{k}", b_sep[k], 3)
    gen.int(f"sl:n{k}", (grp == k).sum())
gen.num("sl:ssepar", sse_par, 1)
gen.num("sl:ssesep", sse_sep, 1)
gen.num("sl:par", b_par[3], 4)
gen.write()

# ---- figure: separate slopes with binned means ---------------------------------------
use_book_style()
fig, ax = plt.subplots(figsize=(4.4, 2.7))
colors = [COLORS["accent"], COLORS["second"], COLORS["third"]]
names = ["codes 1–3 (high school or less)", "codes 4–5", "codes 6–7 (two highest)"]
bins = np.arange(15, 95, 10)
ages = np.linspace(18, 91, 2)
for k in range(3):
    sel = grp == k
    idx = np.digitize(age[sel], bins)
    centres = [age[sel][idx == j].mean() for j in np.unique(idx) if (idx == j).sum() >= 8]
    means = [tv[sel][idx == j].mean() for j in np.unique(idx) if (idx == j).sum() >= 8]
    ax.plot(centres, means, "o", ms=3, color=colors[k])
    ax.plot(ages, b_sep[k] + b_sep[3 + k] * ages, color=colors[k], label=names[k])
ax.set_xlabel("age (years)")
ax.set_ylabel("times per week of TV news")
ax.legend(frameon=False, loc="upper left", fontsize=7)
fig.savefig(figure_path("ch08", "separate_slopes"))
