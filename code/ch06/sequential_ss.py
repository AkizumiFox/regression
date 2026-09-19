"""Chapter 6, Section 5: sequential sums of squares depend on the order of entry."""
import numpy as np
import statsmodels.api as sm

from regbook import Generated

data = sm.datasets.statecrime.load_pandas().data.drop(index="District of Columbia")
y = data["murder"].to_numpy()
cols = {"poverty": data["poverty"], "single": data["single"], "urban": data["urban"]}


# <<sequential>>
def projection(Z):
    Q, _ = np.linalg.qr(Z)
    return Q @ Q.T                       # fine for n = 50; never do this for large n


def sequential_ss(y, columns):
    """Squared lengths ||(M_j - M_{j-1}) y||^2 as columns enter one at a time."""
    n = len(y)
    Z = np.ones((n, 1))
    M_prev = projection(Z)
    pieces = {"mean": y @ M_prev @ y}
    for name, x in columns:
        Z = np.column_stack([Z, x])
        M = projection(Z)
        pieces[name] = y @ (M - M_prev) @ y
        M_prev = M
    pieces["residual"] = y @ (np.eye(n) - M_prev) @ y
    return pieces


order_a = sequential_ss(y, [(k, cols[k]) for k in ["poverty", "single", "urban"]])
order_b = sequential_ss(y, [(k, cols[k]) for k in ["single", "urban", "poverty"]])
for name in ["poverty", "single", "urban", "residual"]:
    print(f"{name:9s}  order A: {order_a[name]:8.2f}   order B: {order_b[name]:8.2f}")
# <</sequential>>

for d in (order_a, order_b):
    assert np.isclose(sum(d.values()), y @ y)
assert np.isclose(order_a["residual"], order_b["residual"])
assert np.isclose(order_a["poverty"] + order_a["single"] + order_a["urban"],
                  order_b["poverty"] + order_b["single"] + order_b["urban"])

gen = Generated("ch06", "sequential_ss", prefix="seq")
for tag, d in (("a", order_a), ("b", order_b)):
    for k, v in d.items():
        gen.num(f"{tag}:{k}", v, 2)
sst = y @ y - order_a["mean"]
gen.num("sst", sst, 2)
gen.num("ssr", sst - order_a["residual"], 2)
gen.write()
