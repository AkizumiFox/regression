"""Chapter 20, Section 2: leverage in observational data, Engel's 1857 household budgets.

Annual income and food expenditure of 235 Belgian working-class households (public domain,
statsmodels.datasets.engel). Income is strongly skewed, so a few households carry large
leverage in a straight-line fit; on the log scale they do not.
"""
import numpy as np
import statsmodels.api as sm

from regbook import Generated

engel = sm.datasets.engel.load_pandas().data

# <<leverage>>
income = engel["income"].to_numpy()
food = engel["foodexp"].to_numpy()
n = len(income)
for label, x in [("income", income), ("log income", np.log(income))]:
    X = np.column_stack([np.ones(n), x])
    Q, _ = np.linalg.qr(X)
    h = np.sum(Q ** 2, axis=1)
    top = np.argsort(h)[::-1][:3]
    print(f"{label:11s} average p/n = {2 / n:.4f};  largest:",
          ", ".join(f"row {i + 1}: h = {h[i]:.3f}" for i in top))
# <</leverage>>

X = np.column_stack([np.ones(n), income])
h = np.sum(np.linalg.qr(X)[0] ** 2, axis=1)
XL = np.column_stack([np.ones(n), np.log(income)])
hL = np.sum(np.linalg.qr(XL)[0] ** 2, axis=1)
top = int(np.argmax(h))
xbar = income.mean()
Sxx = np.sum((income - xbar) ** 2)
assert np.allclose(h, 1 / n + (income - xbar) ** 2 / Sxx)
assert top == int(np.argmax(income))
fit = sm.OLS(food, X).fit()
infl = fit.get_influence()
t = infl.resid_studentized_external
cook = infl.cooks_distance[0]
assert int(np.argmax(cook)) == top
skew = np.mean((income - xbar) ** 3) / np.mean((income - xbar) ** 2) ** 1.5
share = (income[top] - xbar) ** 2 / Sxx

gen = Generated("ch20", "engel_leverage", prefix="el")
gen.int("n", n)
gen.int("top", top + 1)
gen.num("h_top", h[top], 3)
gen.num("avg", 2 / n, 4)
gen.num("ratio", h[top] / (2 / n), 0)
gen.num("share", share, 3)
gen.num("h2", np.sort(h)[-2], 3)
gen.num("hL_top", hL[top], 3)
gen.num("hL_max", hL.max(), 3)
gen.int("hL_argmax", int(np.argmax(hL)) + 1)
gen.num("inc_top", income[top], 0)
gen.num("inc_med", np.median(income), 0)
gen.num("skew", skew, 2)
gen.num("t_top", t[top], 2)
gen.num("D_top", cook[top], 2)
gen.write()
