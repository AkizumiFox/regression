"""Chapter 29, Section 4: comparing selection procedures by simulation.

Eight equicorrelated normal regressors (correlation 0.5), n = 60, sigma = 1, three nonzero
slopes, weak (1, 0.5, 0.25) or strong (1, 0.8, 0.6), intercept always in. Procedures: least squares on everything, the
true subset (an oracle), all subsets by AIC, BIC and leave-one-out CV, and forward and
backward stepwise by AIC. Loss: ||X(beta_hat - beta)||^2 / sigma^2 at the design points.
"""
import itertools

import numpy as np

from regbook import Generated

gen = Generated("ch29", "search_compare", prefix="cmp")

rng = np.random.default_rng(2910)
n, m, rho, reps = 60, 8, 0.5, 2000
true_set = (0, 1, 2)
L = np.linalg.cholesky(rho * np.ones((m, m)) + (1 - rho) * np.eye(m))
subsets = [S for k in range(m + 1) for S in itertools.combinations(range(m), k)]


def analyse(X, y):
    Xc, yc = X - X.mean(axis=0), y - y.mean()
    G, g = Xc.T @ Xc, Xc.T @ yc
    info = {}
    for S in subsets:
        Sl = list(S)
        if Sl:
            b = np.linalg.solve(G[np.ix_(Sl, Sl)], g[Sl])
            fitc = Xc[:, Sl] @ b
            h = 1 / n + np.einsum("ij,ij->i", Xc[:, Sl] @ np.linalg.inv(G[np.ix_(Sl, Sl)]), Xc[:, Sl])
        else:
            b, fitc, h = np.zeros(0), np.zeros(n), np.full(n, 1 / n)
        e = yc - fitc
        rss = e @ e
        info[S] = dict(b=b, rss=rss, aic=n * np.log(rss / n) + 2 * (len(S) + 2),
                       bic=n * np.log(rss / n) + np.log(n) * (len(S) + 2),
                       loo=np.mean((e / (1 - h)) ** 2))
    return info


def stepwise(info, start, direction):
    S = start
    while True:
        if direction == "forward":
            moves = [tuple(sorted(S + (j,))) for j in range(m) if j not in S]
        else:
            moves = [tuple(k for k in S if k != j) for j in S]
        if not moves:
            return S
        nxt = min(moves, key=lambda T: info[T]["aic"])
        if info[nxt]["aic"] >= info[S]["aic"]:
            return S
        S = nxt


def loss(X, S, b, beta):
    """||X_c (beta_hat - beta)||^2 on the slopes; the intercept adds sigma^2 exactly in expectation."""
    full = np.zeros(m)
    full[list(S)] = b
    Xc = X - X.mean(axis=0)
    return np.sum((Xc @ (full - beta)) ** 2)


methods = ["full", "oracle", "AIC", "BIC", "LOO", "forward", "backward"]


def run(beta, tag):
    losses = {k: np.empty(reps) for k in methods}
    hits = {k: np.zeros(reps, bool) for k in methods}
    sizes = {k: np.empty(reps) for k in methods}
    same = np.zeros(reps, bool)
    for r in range(reps):
        X = rng.normal(size=(n, m)) @ L.T
        y = 1 + X @ beta + rng.normal(size=n)
        info = analyse(X, y)
        choice = {"full": tuple(range(m)), "oracle": true_set,
                  "AIC": min(subsets, key=lambda S: info[S]["aic"]),
                  "BIC": min(subsets, key=lambda S: info[S]["bic"]),
                  "LOO": min(subsets, key=lambda S: info[S]["loo"]),
                  "forward": stepwise(info, (), "forward"),
                  "backward": stepwise(info, tuple(range(m)), "backward")}
        same[r] = choice["forward"] == choice["AIC"]
        for k, S in choice.items():
            losses[k][r] = loss(X, S, info[S]["b"], beta) + 1.0     # + 1 for the intercept
            hits[k][r] = S == true_set
            sizes[k][r] = len(S)
    for k in methods:
        print(f"{tag} {k:9s} loss {losses[k].mean():5.2f}  Pr(true) {hits[k].mean():.2f}  size {sizes[k].mean():.2f}")
        gen.num(f"{tag}_loss_{k}", losses[k].mean(), 2)
        gen.num(f"{tag}_hit_{k}", hits[k].mean(), 2)
        gen.num(f"{tag}_size_{k}", sizes[k].mean(), 2)
    gen.num(f"{tag}_same", same.mean(), 2)
    assert abs(losses["full"].mean() - (m + 1)) < 0.3      # sigma^2 p for least squares
    assert abs(losses["oracle"].mean() - 4) < 0.2
    return losses, hits


weak_l, weak_h = run(np.r_[1.0, 0.5, 0.25, np.zeros(m - 3)], "weak")
strong_l, strong_h = run(np.r_[1.0, 0.8, 0.6, np.zeros(m - 3)], "strong")
assert strong_h["BIC"].mean() > 0.6 > strong_h["AIC"].mean()
assert strong_l["BIC"].mean() < strong_l["AIC"].mean() < strong_l["full"].mean()
assert weak_l["AIC"].mean() < weak_l["full"].mean()
gen.int("n", n)
gen.int("m", m)
gen.int("reps", reps)
gen.write()
