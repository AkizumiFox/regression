"""Chapter 5, Sections 3 and 7: indicator variables for a categorical regressor.

RAND Health Insurance Experiment (public domain subset shipped with statsmodels,
statsmodels.datasets.randhie). Response mdvis = number of outpatient visits to a
physician in a year. Self-rated health has four levels: excellent (the reference),
good (hlthg), fair (hlthf), poor (hlthp).
"""
import numpy as np
import statsmodels.api as sm

from regbook import Generated

# <<dummies>>
rand = sm.datasets.randhie.load_pandas().data
y = rand["mdvis"].to_numpy(dtype=float)
D = rand[["hlthg", "hlthf", "hlthp"]].to_numpy(dtype=float)   # good, fair, poor
X = np.column_stack([np.ones(len(y)), D])                      # excellent = reference

beta_hat = np.linalg.solve(X.T @ X, X.T @ y)
level = D @ np.array([1, 2, 3])                                # 0 = excellent, ..., 3 = poor
means = np.array([y[level == k].mean() for k in range(4)])
print("coefficients:            ", np.round(beta_hat, 4))
print("group means:             ", np.round(means, 4))
print("means minus reference:   ", np.round(means - means[0], 4))
# <</dummies>>

assert np.allclose(beta_hat[0], means[0]) and np.allclose(beta_hat[1:], means[1:] - means[0])

# <<trap>>
excellent = 1.0 - D.sum(axis=1)
X_all = np.column_stack([X, excellent])          # an indicator for every level, and 1
print("columns:", X_all.shape[1], " rank:", np.linalg.matrix_rank(X_all))
print("1 - (sum of the four indicators) =", np.abs(X_all[:, 0] - X_all[:, 1:].sum(axis=1)).max())
# <</trap>>
assert np.linalg.matrix_rank(X_all) == 4

# <<adjusted>>
Xa = np.column_stack([X, rand["disea"], rand["physlm"]])
beta_a = np.linalg.solve(Xa.T @ Xa, Xa.T @ y)
print("adjusted for chronic-disease score and physical limitation:")
print(np.round(beta_a, 4))
# <</adjusted>>

counts = np.array([(level == k).sum() for k in range(4)])
gen = Generated("ch05", "rand_dummies", prefix="rand")
gen.text("n", f"{len(y):,}")
for k, name in enumerate(["exc", "good", "fair", "poor"]):
    gen.num("mean" + name, means[k], 3)
    gen.text("count" + name, f"{counts[k]:,}")
for j, name in enumerate(["b0", "bg", "bf", "bp"]):
    gen.num(name, beta_hat[j], 3)
for j, name in enumerate(["b0", "bg", "bf", "bp", "dis", "phys"]):
    gen.num("a" + name, beta_a[j], 3)
gen.num("meandisea", rand["disea"].mean(), 2)
gen.write()
