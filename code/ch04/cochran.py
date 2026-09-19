"""Chapter 4, Section 5: a Cochran decomposition checked by simulation.

One-way layout with three groups of sizes 3, 5, 4. The identity
I = J/n + (M - J/n) + (I - M) splits ||Y||^2 into the grand-mean, between-group and
within-group sums of squares. Ranks add to n, so by Cochran's theorem the three pieces
(divided by sigma^2) are independent noncentral chi-squared variables.
"""
import numpy as np
from scipy import stats

from regbook import Generated

rng = np.random.default_rng(5150)

# <<cochran>>
sizes, means, sigma = [3, 5, 4], [2.0, 3.0, 1.5], 1.0
groups = np.repeat(np.arange(3), sizes)
n = len(groups)
Z = (groups[:, None] == np.arange(3)).astype(float)    # group indicators
M = Z @ np.diag(1 / np.array(sizes)) @ Z.T             # group-mean projection
J = np.ones((n, n)) / n
pieces = [J, M - J, np.eye(n) - M]
ranks = [np.linalg.matrix_rank(P) for P in pieces]
print("ranks", ranks, "sum", sum(ranks), "n", n)

mu = np.array(means)[groups]
gammas = [mu @ P @ mu / sigma**2 for P in pieces]
Y = mu + sigma * rng.standard_normal((200_000, n))
Q = np.stack([np.einsum("ij,jk,ik->i", Y, P, Y) for P in pieces], axis=1) / sigma**2
print(np.corrcoef(Q, rowvar=False).round(3))
# <</cochran>>

assert sum(ranks) == n
assert all(np.allclose(P @ P, P) for P in pieces)
assert all(np.allclose(pieces[i] @ pieces[j], 0) for i in range(3) for j in range(3) if i != j)
assert np.isclose(gammas[2], 0)
theory = [(r + g, 2 * r + 4 * g) for r, g in zip(ranks, gammas)]
for k in range(3):
    m, v = theory[k]
    assert abs(Q[:, k].mean() - m) < 5 * np.sqrt(v / len(Q))
    assert stats.kstest(Q[:, k], stats.ncx2(ranks[k], max(gammas[k], 1e-12)).cdf).pvalue > 1e-3
R = np.corrcoef(Q, rowvar=False)
assert np.all(np.abs(R[np.triu_indices(3, 1)]) < 0.01)
# a sharper check of independence than zero correlation: a joint tail probability factorizes
c = np.median(Q, axis=0)
joint = np.mean((Q[:, 1] > c[1]) & (Q[:, 2] > c[2]))
assert abs(joint - 0.25) < 0.005

# between-group noncentrality equals sum n_k (mu_k - mu_bar)^2 / sigma^2
mbar = mu.mean()
assert np.isclose(gammas[1], sum(s * (m - mbar) ** 2 for s, m in zip(sizes, means)) / sigma**2)

gen = Generated("ch04", "cochran", prefix="coch")
gen.int("n", n)
gen.int("g", len(sizes))
gen.text("sizes", ", ".join(str(k) for k in sizes))
gen.text("means", ", ".join(f"{m:g}" for m in means))
gen.num("sigma", sigma, 0)
gen.text("reps", f"{len(Y):,}".replace(",", "{,}"))
names = ["mean", "between", "within"]
for k, name in enumerate(names):
    gen.int(f"r{name}", ranks[k])
    gen.num(f"g{name}", 0.0 if abs(gammas[k]) < 1e-9 else gammas[k], 3)
    gen.num(f"m{name}", theory[k][0], 3)
    gen.num(f"v{name}", theory[k][1], 2)
    gen.num(f"sm{name}", Q[:, k].mean(), 3)
    gen.num(f"sv{name}", Q[:, k].var(), 2)
gen.num("corrmax", np.max(np.abs(R[np.triu_indices(3, 1)])), 4)
gen.num("joint", joint, 4)
gen.write()
