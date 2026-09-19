"""Chapter 29, Section 5: three ways to report intervals after forward selection.

n = 100, ten equicorrelated normal regressors (correlation 0.5), sigma = 1, slopes
(0.3, 0.3, 0.3, 0, ..., 0) or all zero. Forward selection by AIC. The target of an interval
for x_j in the selected model M is the projection coefficient
beta_{M,j} = [(X_M' X_M)^{-1} X_M' mu]_j, which is what the fitted coefficient estimates.
  naive:  usual t intervals from the selected model, same data
  split:  select on 50 cases, t intervals on the other 50 (for their own design)
  Scheffe: sqrt(10 F) multiplier with the full-model s, same data (valid after any selection)
"""
import numpy as np
from scipy import stats

from regbook import Generated

gen = Generated("ch29", "splitting", prefix="spl")

# <<split>>
import numpy as np
from scipy import stats
rng = np.random.default_rng(2915)
n, m, rho = 100, 10, 0.5
L = np.linalg.cholesky(rho * np.ones((m, m)) + (1 - rho) * np.eye(m))


def design(cols, X):
    return np.column_stack([np.ones(len(X)), X[:, cols]])


def forward_aic(X, y):
    S = []
    def aic(cols):
        A = design(cols, X)
        e = y - A @ np.linalg.lstsq(A, y, rcond=None)[0]
        return len(y) * np.log(e @ e / len(y)) + 2 * A.shape[1]
    while len(S) < m:
        j = min((j for j in range(m) if j not in S), key=lambda j: aic(S + [j]))
        if aic(S + [j]) >= aic(S):
            break
        S.append(j)
    return S


def intervals(X, y, S, mult, s=None):
    """Estimates, half-widths and projection targets for the slopes of model S."""
    A = design(S, X)
    G = np.linalg.inv(A.T @ A)
    b = G @ A.T @ y
    if s is None:
        s = np.sqrt(np.sum((y - A @ b) ** 2) / (len(y) - A.shape[1]))
    return b[1:], mult(len(y) - A.shape[1]) * s * np.sqrt(np.diag(G)[1:])


def one_run(beta):
    X = rng.normal(size=(n, m)) @ L.T
    mu = X @ beta
    y = mu + rng.normal(size=n)
    out = {}
    # naive
    S = forward_aic(X, y)
    if S:
        b, h = intervals(X, y, S, lambda df: stats.t.ppf(0.975, df))
        target = np.linalg.lstsq(design(S, X), mu, rcond=None)[0][1:]
        out["naive"] = (np.abs(b - target) <= h, h)
        # Scheffe, with s from the full model and q = m
        A = design(list(range(m)), X)
        s_full = np.sqrt(np.sum((y - A @ np.linalg.lstsq(A, y, rcond=None)[0]) ** 2) / (n - m - 1))
        b, h = intervals(X, y, S, lambda df: np.sqrt(m * stats.f.ppf(0.95, m, n - m - 1)), s_full)
        out["scheffe"] = (np.abs(b - target) <= h, h)
    # split: choose on the first half, infer on the second
    half = n // 2
    S = forward_aic(X[:half], y[:half])
    if S:
        b, h = intervals(X[half:], y[half:], S, lambda df: stats.t.ppf(0.975, df))
        target = np.linalg.lstsq(design(S, X[half:]), mu[half:], rcond=None)[0][1:]
        out["split"] = (np.abs(b - target) <= h, h)
    return out


def summary(beta, reps):
    cov = {k: [] for k in ("naive", "split", "scheffe")}
    width = {k: [] for k in cov}
    for _ in range(reps):
        for k, (c, h) in one_run(beta).items():
            cov[k].extend(c)
            width[k].extend(h)
    return {k: (np.mean(cov[k]), np.mean(width[k])) for k in cov}


for label, beta in (("null", np.zeros(m)), ("signal", np.r_[0.3, 0.3, 0.3, np.zeros(m - 3)])):
    res = summary(beta, 200)
    print(label, {k: (round(c, 3), round(w, 3)) for k, (c, w) in res.items()})
# <</split>>

for label, beta in (("null", np.zeros(m)), ("signal", np.r_[0.3, 0.3, 0.3, np.zeros(m - 3)])):
    res = summary(beta, 3000)
    print(label, res)
    for k, (c, w) in res.items():
        gen.num(f"{label}_cov_{k}", c, 3)
        gen.num(f"{label}_width_{k}", w, 3)
    assert res["naive"][0] < 0.93 and abs(res["split"][0] - 0.95) < 0.015 and res["scheffe"][0] > 0.99
gen.num("sch_mult", np.sqrt(m * stats.f.ppf(0.95, m, n - m - 1)), 3)
gen.write()
