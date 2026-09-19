"""Chapter 29, Section 3: cross-validation must contain the whole procedure.

Pure noise: n = 50 cases, 200 candidate regressors, a response unrelated to all of them.
The procedure screens the 5 regressors most correlated with the response and fits least
squares on them. Ten-fold CV applied after screening on all the data is badly optimistic;
CV that repeats the screening inside every fold is honest.
"""
import numpy as np

from regbook import Generated

gen = Generated("ch29", "nested_cv", prefix="ncv")

# <<nested>>
import numpy as np
rng = np.random.default_rng(2907)
n, m, keep, K = 50, 200, 5, 10


def screen(X, y):
    """Indices of the `keep` columns most correlated with y."""
    Xc, yc = X - X.mean(axis=0), y - y.mean()
    corr = np.abs(Xc.T @ yc) / np.sqrt(np.sum(Xc**2, axis=0) * (yc @ yc))
    return np.argsort(-corr)[:keep]


def fit(X, y):
    A = np.column_stack([np.ones(len(y)), X])
    return np.linalg.lstsq(A, y, rcond=None)[0]


def predict(b, X):
    return b[0] + X @ b[1:]


def one_data_set():
    X = rng.normal(size=(n, m))
    y = rng.normal(size=n)                                   # unrelated to every column
    folds = rng.permutation(np.arange(n) % K)
    cols = screen(X, y)                                      # screening on ALL the data
    wrong = right = 0.0
    for f in range(K):
        tr, te = folds != f, folds == f
        b = fit(X[tr][:, cols], y[tr])                       # wrong: screening done outside
        wrong += np.sum((y[te] - predict(b, X[te][:, cols])) ** 2)
        cols_f = screen(X[tr], y[tr])                        # right: screen again in the fold
        b = fit(X[tr][:, cols_f], y[tr])
        right += np.sum((y[te] - predict(b, X[te][:, cols_f])) ** 2)
    b = fit(X[:, cols], y)
    truth = 1 + b[0] ** 2 + np.sum(b[1:] ** 2)               # error on a new case
    return wrong / n, right / n, truth


res = np.array([one_data_set() for _ in range(200)])
print("average CV, screening outside the folds:", res[:, 0].mean().round(3))
print("average CV, screening inside the folds: ", res[:, 1].mean().round(3))
print("average true error of the fitted model: ", res[:, 2].mean().round(3))
# <</nested>>

assert res[:, 0].mean() < 0.85 and res[:, 1].mean() > 1.1
assert abs(res[:, 1].mean() - res[:, 2].mean()) < 0.1
gen.int("n", n)
gen.int("m", m)
gen.int("keep", keep)
gen.num("wrong", res[:, 0].mean(), 3)
gen.num("right", res[:, 1].mean(), 3)
gen.num("truth", res[:, 2].mean(), 3)
gen.num("frac_below", np.mean(res[:, 0] < 1), 2)
gen.write()
