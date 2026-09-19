"""Chapter 29, Section 4: searching the subsets.

(a) Branch and bound for the C_p-best subset among 2^15 subsets: the answer of exhaustive
    search, with a fraction of the fits. The bound: deleting columns never lowers SSE.
(b) A design where forward selection misses the best pair: two strongly correlated
    regressors whose difference carries the signal.
(c) Forward and backward selection by AIC on the state data.
"""
import itertools

import numpy as np
import statsmodels.api as sm

from regbook import Generated

gen = Generated("ch29", "subset_search", prefix="srch")

# ---- (a) branch and bound -----------------------------------------------------------------------
# <<bb>>
import itertools
import numpy as np
rng = np.random.default_rng(2908)
n, m = 100, 15
X = rng.normal(size=(n, m))
beta = np.r_[1.0, -0.8, 0.6, 0.4, 0.3, np.zeros(m - 5)]
y = 2 + X @ beta + rng.normal(size=n)
Xc, yc = X - X.mean(axis=0), y - y.mean()                  # centring accounts for the intercept
G, g = Xc.T @ Xc, Xc.T @ yc


def sse(S):
    S = list(S)
    if not S:
        return yc @ yc
    return yc @ yc - g[S] @ np.linalg.solve(G[np.ix_(S, S)], g[S])


s2 = sse(range(m)) / (n - m - 1)


def cp(S):
    return sse(S) / s2 - n + 2 * (len(S) + 1)


def branch_and_bound():
    """Children delete a column with a larger index than the last one deleted, so every
    subset is visited at most once; a node is not expanded when no descendant can win."""
    best = [np.inf, None]
    count = [0]

    def visit(S, last):
        e = sse(S)
        count[0] += 1
        value = e / s2 - n + 2 * (len(S) + 1)
        if value < best[0]:
            best[:] = [value, S]
        forced = [j for j in S if j <= last]              # kept by every descendant
        if e / s2 - n + 2 * (len(forced) + 1) >= best[0]:
            return                                          # SSE only grows as columns go
        for j in S:
            if j > last:
                visit(tuple(k for k in S if k != j), j)

    visit(tuple(range(m)), -1)
    return best[1], best[0], count[0]


S_bb, cp_bb, fits = branch_and_bound()
print("branch and bound:", S_bb, round(cp_bb, 3), "after", fits, "fits of", 2**m)
# <</bb>>

S_ex = min((S for k in range(m + 1) for S in itertools.combinations(range(m), k)), key=cp)
assert S_ex == S_bb and np.isclose(cp(S_ex), cp_bb)
assert S_bb == (0, 1, 2, 3, 6, 7)                        # four largest slopes + two noise columns
assert fits < 2**m / 4
gen.int("m", m)
gen.int("nsub", 2**m)
gen.int("fits", fits)
gen.num("frac", fits / 2**m, 3)
gen.int("size_bb", len(S_bb))

# ---- (b) forward selection misses the best pair --------------------------------------------------
# <<forward-miss>>
rng = np.random.default_rng(2909)
n2 = 60
z = rng.normal(size=(n2, 3))
x1 = z[:, 0]
x2 = 0.95 * z[:, 0] + np.sqrt(1 - 0.95**2) * z[:, 1]      # corr(x1, x2) about 0.95
x3 = z[:, 2]
Xf = np.column_stack([x1, x2, x3])
yf = 3 * (x1 - x2) + 0.4 * x3 + rng.normal(size=n2)       # the signal is in x1 - x2


def rss(cols):
    A = np.column_stack([np.ones(n2)] + [Xf[:, j] for j in cols])
    b, *_ = np.linalg.lstsq(A, yf, rcond=None)
    return np.sum((yf - A @ b) ** 2)


first = min(range(3), key=lambda j: rss([j]))
second = min((j for j in range(3) if j != first), key=lambda j: rss([first, j]))
best_pair = min(itertools.combinations(range(3), 2), key=rss)
print("forward chooses", sorted([first, second]), "RSS", round(rss([first, second]), 1))
print("best pair      ", list(best_pair), "RSS", round(rss(best_pair), 1))
# <</forward-miss>>
assert sorted([first, second]) != list(best_pair) and best_pair == (0, 1)
tss = np.sum((yf - yf.mean()) ** 2)
gen.num("r_x12", np.corrcoef(x1, x2)[0, 1], 2)
gen.int("first", first + 1)
gen.num("r2_fwd", 1 - rss([first, second]) / tss, 2)
gen.num("r2_best", 1 - rss(best_pair) / tss, 2)
gen.num("r_y1", np.corrcoef(yf, x1)[0, 1], 2)
gen.num("r_y2", np.corrcoef(yf, x2)[0, 1], 2)
gen.num("r_y3", np.corrcoef(yf, x3)[0, 1], 2)

# ---- (c) stepwise on the state data --------------------------------------------------------------
data = sm.datasets.statecrime.load_pandas().data.drop(index="District of Columbia")
names = ["hs_grad", "poverty", "single", "white", "urban"]
ys = np.log(data["violent"].to_numpy())
ns = len(ys)


evaluated = set()


def aic(cols):
    evaluated.add(frozenset(cols))
    A = np.column_stack([np.ones(ns)] + [data[c].to_numpy() for c in cols])
    b, *_ = np.linalg.lstsq(A, ys, rcond=None)
    return ns * np.log(np.sum((ys - A @ b) ** 2) / ns) + 2 * (A.shape[1] + 1)


def forward():
    S, path = [], []
    while True:
        cand = [c for c in names if c not in S]
        if not cand:
            return S, path
        c = min(cand, key=lambda c: aic(S + [c]))
        if aic(S + [c]) >= aic(S):
            return S, path
        S = S + [c]
        path.append(c)


def backward():
    S, path = list(names), []
    while S:
        c = min(S, key=lambda c: aic([d for d in S if d != c]))
        if aic([d for d in S if d != c]) >= aic(S):
            break
        S = [d for d in S if d != c]
        path.append(c)
    return S, path


f_model, f_path = forward()
n_forward = len(evaluated)
evaluated.clear()
b_model, b_path = backward()
n_backward = len(evaluated)
print("models evaluated: forward", n_forward, "backward", n_backward)
print("forward:", f_path, " backward removes", b_path, "->", b_model)
assert sorted(f_model) == sorted(b_model) == ["single", "urban"]
gen.text("b_path", ", ".join(b_path))
gen.int("n_forward", n_forward)
gen.int("n_backward", n_backward)
gen.write()
