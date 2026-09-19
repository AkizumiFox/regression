"""Chapter 6, Section 8: leverages are the diagonal of the projection M."""
import numpy as np
import statsmodels.api as sm

from regbook import Generated

data = sm.datasets.statecrime.load_pandas().data          # District of Columbia kept here
X = np.column_stack([np.ones(len(data)), data["poverty"], data["single"], data["urban"]])

# <<leverage>>
Q, _ = np.linalg.qr(X)                 # orthonormal basis for C(X)
h = np.sum(Q ** 2, axis=1)             # h_ii = m_ii = ||row i of Q||^2
n, p = X.shape
print("sum of leverages:", h.sum(), " (= p =", p, ")")
print("range:", h.min(), h.max(), " bounds:", 1 / n, 1)
top = np.argsort(h)[::-1][:3]
for i in top:
    print(f"{data.index[i]:22s} h = {h[i]:.3f}   ({h[i] / (p / n):.1f} x average)")
# <</leverage>>

M = X @ np.linalg.solve(X.T @ X, X.T)
assert np.allclose(np.diag(M), h)
assert np.isclose(h.sum(), p)
assert np.all(h >= 1 / n - 1e-12) and np.all(h <= 1 + 1e-12)
gen = Generated("ch06", "leverage", prefix="lev")
gen.int("n", n)
gen.num("avg", p / n, 3)
gen.num("hmin", h.min(), 3)
for k, i in enumerate(top, 1):
    gen.text(f"name{k}", data.index[i])
    gen.num(f"h{k}", h[i], 3)
    gen.num(f"ratio{k}", h[i] / (p / n), 1)
gen.write()
