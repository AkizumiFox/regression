"""Chapter 21, Section 2: the size of the Breusch-Pagan test and of Koenker's studentized version.

Simple regression with n = 100 equally spaced x values; the variance regressor is x itself.
Homoscedastic errors with variance one from three laws of different kurtosis (normal 3,
Laplace 6, uniform 1.8). Under H0 the
Koenker statistic n R^2 tends to chi^2(1) whatever the error law, while the Breusch-Pagan
statistic tends to (kappa - 1)/2 times chi^2(1) (Theorem 21.2.2).
"""
import numpy as np
from scipy import stats

from regbook import Generated

# <<size>>
rng = np.random.default_rng(2105)
n = 100
x = np.linspace(0, 10, n)
X = np.column_stack([np.ones(n), x])
Q, _ = np.linalg.qr(X)
# centred, unit-length variance regressor
zc = (x - x.mean()) / np.linalg.norm(x - x.mean())
crit = stats.chi2.ppf(0.95, 1)
laws = {
    "normal": (lambda size: rng.normal(size=size), 3.0),
    "Laplace": (lambda size: rng.laplace(scale=np.sqrt(0.5), size=size), 6.0),
    "uniform": (lambda size: rng.uniform(-np.sqrt(3), np.sqrt(3), size=size),
                1.8),
}


def size_table(reps):
    table = {}
    for name, (draw, kappa) in laws.items():
        E = draw((reps, n))
        R = E - (E @ Q) @ Q.T                          # residuals
        R2 = R**2
        sig2 = R2.mean(axis=1, keepdims=True)
        # explained SS of e^2 on [1, x]
        ess = (R2 @ zc) ** 2
        bp = ess / (2 * sig2[:, 0] ** 2)
        koenker = n * ess / np.sum((R2 - sig2) ** 2, axis=1)
        limit_bp = stats.chi2.sf(crit / ((kappa - 1) / 2), 1)
        table[name] = (np.mean(bp > crit), np.mean(koenker > crit), limit_bp)
        print(f"{name:8s}: Breusch-Pagan rejects {table[name][0]:.3f} "
              f"(limit {limit_bp:.3f}), "
              f"Koenker rejects {table[name][1]:.3f}")
    return table


# a quick version (the book uses 40000)
quick = size_table(2000)
# <</size>>

table = size_table(40_000)
gen = Generated("ch21", "bp_size", prefix="bps")
gen.int("n", n)
tags = {"normal": "norm", "Laplace": "lap", "uniform": "unif"}
for name, (bp_rate, k_rate, limit_bp) in table.items():
    gen.num(f"bp_{tags[name]}", bp_rate, 3)
    gen.num(f"k_{tags[name]}", k_rate, 3)
    gen.num(f"limit_{tags[name]}", limit_bp, 3)
    assert abs(k_rate - 0.05) < 0.02, (name, k_rate)
assert table["Laplace"][0] > 0.15 and table["uniform"][0] < 0.01
assert abs(table["normal"][0] - 0.05) < 0.01
gen.write()
