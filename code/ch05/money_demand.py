"""Chapter 5, Sections 4, 6 and 7: a multiple regression solved from the normal equations.

Danish money demand, quarterly 1974-1987 (Johansen and Juselius 1990; public domain,
statsmodels.datasets.danish_data). Response lrm = log real money (M2). Regressors:
lry = log real income, ibo = bond rate, ide = deposit rate (rates as fractions).
"""
import numpy as np
import statsmodels.api as sm

from regbook import Generated

# <<data>>
dk = sm.datasets.danish_data.load_pandas().data
y = dk["lrm"].to_numpy()
X = np.column_stack([np.ones(len(y)), dk["lry"], dk["ibo"], dk["ide"]])
n, p = X.shape
print("n =", n, " p =", p)
# <</data>>

# <<normal>>
XtX = X.T @ X
Xty = X.T @ y
beta_hat = np.linalg.solve(XtX, Xty)       # solve X^T X b = X^T y
fitted = X @ beta_hat
resid = y - fitted
SSE = resid @ resid

print("beta_hat =", np.round(beta_hat, 4))
print("X^T e    =", X.T @ resid)           # zero up to rounding
print("SSE two ways:", SSE, y @ y - beta_hat @ Xty)
# <</normal>>

assert np.allclose(X.T @ resid, 0, atol=1e-10)
assert np.isclose(SSE, y @ y - beta_hat @ Xty, rtol=1e-6)
lst, *_ = np.linalg.lstsq(X, y, rcond=None)
assert np.allclose(beta_hat, lst, rtol=1e-8)
H = X @ np.linalg.solve(XtX, X.T)
assert np.allclose(H, H.T) and np.allclose(H @ H, H) and np.isclose(np.trace(H), p)

# <<centred>>
Z = X[:, 1:] - X[:, 1:].mean(axis=0)       # centred regressors, no intercept column
slopes = np.linalg.solve(Z.T @ Z, Z.T @ (y - y.mean()))
intercept = y.mean() - X[:, 1:].mean(axis=0) @ slopes
print("slopes from centred data:", np.round(slopes, 4))
print("intercept recovered:     ", round(intercept, 4))
# <</centred>>

assert np.allclose(slopes, beta_hat[1:]) and np.isclose(intercept, beta_hat[0])

# <<sigma2>>
s2 = SSE / (n - p)                          # unbiased estimate of sigma^2
cov_hat = s2 * np.linalg.inv(XtX)           # estimated Cov(beta_hat)
se = np.sqrt(np.diag(cov_hat))
print(f"s^2 = {s2:.6f}, s = {np.sqrt(s2):.5f}")
print("standard errors:", np.round(se, 4))
print("SSE / n would give", SSE / n)
# <</sigma2>>

ols = sm.OLS(y, X).fit()
assert np.allclose(ols.params, beta_hat) and np.allclose(ols.bse, se)
assert np.isclose(ols.scale, s2)

# <<rescale>>
X_pct = X.copy()
X_pct[:, 2:] *= 100                         # rates in percentage points
beta_pct = np.linalg.solve(X_pct.T @ X_pct, X_pct.T @ y)
print("rates as fractions:", np.round(beta_hat, 4))
print("rates in percent:  ", np.round(beta_pct, 5))
print("same fitted values:", np.allclose(X @ beta_hat, X_pct @ beta_pct))

sd = X[:, 1:].std(axis=0, ddof=1)
standardized = beta_hat[1:] * sd / y.std(ddof=1)
print("standardized slopes:", np.round(standardized, 3))
# <</rescale>>

assert np.allclose(beta_pct[2:], beta_hat[2:] / 100) and np.allclose(beta_pct[:2], beta_hat[:2])
assert np.allclose(X @ beta_hat, X_pct @ beta_pct)
# general nonsingular reparameterization Z = X K
K = np.array([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, -1, 1.0]])
beta_K = np.linalg.solve((X @ K).T @ (X @ K), (X @ K).T @ y)
assert np.allclose(beta_K, np.linalg.solve(K, beta_hat))

# <<omitted>>
X1 = X[:, :2]                               # intercept and log income
X2 = X[:, 2:]                               # the two interest rates
b_short = np.linalg.solve(X1.T @ X1, X1.T @ y)
A = np.linalg.solve(X1.T @ X1, X1.T @ X2)   # regress each rate on (1, lry)
print("short regression slope on income:", round(b_short[1], 4))
print("long slope + A @ (rate coefs):   ", round(beta_hat[1] + A[1] @ beta_hat[2:], 4))
print("auxiliary slopes of the rates on income:", np.round(A[1], 4))
# <</omitted>>

assert np.allclose(b_short, beta_hat[:2] + A @ beta_hat[2:])

# effect of a one percentage point rise in the bond rate on money holdings
pct_bond = 100 * (np.exp(beta_hat[2] * 0.01) - 1)
# lag-one autocorrelation of residuals (the errors of a quarterly series are not independent)
r1 = np.corrcoef(resid[:-1], resid[1:])[0, 1]
simple_income = np.polyfit(X[:, 1], y, 1)[0]
corr_ibo_ide = np.corrcoef(X[:, 2], X[:, 3])[0, 1]
corr_lry_ibo = np.corrcoef(X[:, 1], X[:, 2])[0, 1]

gen = Generated("ch05", "money_demand", prefix="dk")
gen.int("n", n)
gen.int("p", p)
for j, name in enumerate(["b0", "b1", "b2", "b3"]):
    gen.num(name, beta_hat[j], 4)
    gen.num("se" + name[1], se[j], 4)
gen.num("SSE", SSE, 5)
gen.num("s2", s2, 6)
gen.num("s", np.sqrt(s2), 5)
gen.num("sigmaml", np.sqrt(SSE / n), 5)
gen.num("pctbond", pct_bond, 2)
gen.num("r1", r1, 2)
gen.num("simpleincome", simple_income, 3)
gen.num("corribo", corr_ibo_ide, 2)
gen.num("corrlryibo", corr_lry_ibo, 2)
for j, name in enumerate(["lry", "ibo", "ide"]):
    gen.num("std" + name, standardized[j], 3)
    gen.num("mean" + name, X[:, j + 1].mean(), 4)
    gen.num("sd" + name, sd[j], 4)
gen.num("ybar", y.mean(), 4)
gen.num("sdy", y.std(ddof=1), 4)
gen.num("R2", ols.rsquared, 4)
gen.num("Aibo", A[1, 0], 4)
gen.num("Aide", A[1, 1], 4)
gen.num("bias", A[1] @ beta_hat[2:], 4)
gen.write()
