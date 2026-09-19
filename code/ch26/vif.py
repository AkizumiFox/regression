"""Chapter 26, Section 2: variance inflation factors.

(a) Nonessential collinearity: x and x^2 on [10, 20] have huge VIFs that centring removes,
    while the variance of every function that means the same thing is unchanged.
(b) VIFs of the macroeconomic consumption regression (public domain,
    statsmodels.datasets.macrodata): large VIFs, yet a precisely estimated income slope.
(c) Generalized VIF for a factor: Grunfeld's investment data (public domain,
    statsmodels.datasets.grunfeld), value, capital and an 11-level firm factor. The GVIF equals
    the product of 1/(1 - rho^2) over canonical correlations and does not depend on coding.
"""
import numpy as np
import pandas as pd
import statsmodels.api as sm

from regbook import Generated

gen = Generated("ch26", "vif")

# <<macro>>
import numpy as np
import statsmodels.api as sm

def vifs(Z):
    """VIFs of the columns of Z (an intercept is always included): diagonal of R^{-1}."""
    return np.diag(np.linalg.inv(np.corrcoef(Z, rowvar=False)))

macro = sm.datasets.macrodata.load_pandas().data
names = ["realdpi", "pop", "cpi", "m1", "tbilrate", "unemp"]
Z = macro[names].to_numpy()
fit = sm.OLS(macro["realcons"].to_numpy(), sm.add_constant(Z)).fit()
for name, v, b, se in zip(names, vifs(Z), fit.params[1:], fit.bse[1:]):
    print(f"{name:9s} VIF {v:7.1f}   estimate {b:9.4f}   se {se:7.4f}   t {b / se:6.1f}")
# <</macro>>

v_mac = vifs(Z)
# VIF_j = 1/(1 - R_j^2) from the auxiliary regressions
for j in range(len(names)):
    aux = sm.OLS(Z[:, j], sm.add_constant(np.delete(Z, j, axis=1))).fit()
    assert np.isclose(v_mac[j], 1 / (1 - aux.rsquared))
# Var(b_j) = sigma^2 VIF_j / S_jj with sigma^2 estimated by s^2
S = ((Z - Z.mean(0)) ** 2).sum(0)
assert np.allclose(fit.bse[1:] ** 2, fit.scale * v_mac / S)
for j, nm in enumerate(names):
    gen.num(f"vif_{nm}", v_mac[j], 1)
gen.num("b_dpi", fit.params[1], 4)
gen.num("se_dpi", fit.bse[1], 4)
gen.num("t_dpi", fit.params[1] / fit.bse[1], 1)
gen.num("sqrt_vif_dpi", np.sqrt(v_mac[0]), 1)
gen.int("macro_n", len(Z))

# <<quadratic>>
x = np.linspace(10, 20, 21)                           # a regressor far from its origin

raw = np.column_stack([x, x**2])
cen = np.column_stack([x - x.mean(), (x - x.mean()) ** 2])
print("VIFs, raw     :", np.round(vifs(raw), 1))
print("VIFs, centred :", np.round(vifs(cen), 3))

X_raw = np.column_stack([np.ones_like(x), raw])
X_cen = np.column_stack([np.ones_like(x), cen])
G_raw, G_cen = np.linalg.inv(X_raw.T @ X_raw), np.linalg.inv(X_cen.T @ X_cen)
a_mean = np.array([0, 1, 2 * x.mean()])               # slope at x = 15 in the raw coefficients
print("sd of the slope at 0 (raw b1)       :", round(np.sqrt(G_raw[1, 1]), 3))
print("sd of the slope at 15, raw coding   :", round(np.sqrt(a_mean @ G_raw @ a_mean), 4))
print("sd of the slope at 15, centred (b1) :", round(np.sqrt(G_cen[1, 1]), 4))
# <</quadratic>>

v_raw, v_cen = vifs(raw), vifs(cen)
assert np.isclose(v_raw[0], v_raw[1]) and v_raw[0] > 100
assert np.allclose(v_cen, 1)                          # symmetric design: x - 15 and its square are uncorrelated
assert np.isclose(a_mean @ G_raw @ a_mean, G_cen[1, 1])
# the raw slope at 0 is the centred-coding function b1 - 2*15*b2
a0 = np.array([0, 1, -2 * x.mean()])
assert np.isclose(a0 @ G_cen @ a0, G_raw[1, 1])
gen.num("quad_vif_raw", v_raw[0], 1)
gen.num("quad_r2", 1 - 1 / v_raw[0], 4)
gen.num("quad_sd_slope0", np.sqrt(G_raw[1, 1]), 3)
gen.num("quad_sd_slope15", np.sqrt(G_cen[1, 1]), 4)
gen.int("quad_n", len(x))

# <<gvif>>
import pandas as pd
grun = sm.datasets.grunfeld.load_pandas().data
D = pd.get_dummies(grun["firm"], drop_first=True).to_numpy(dtype=float)   # 10 indicator columns
Zg = np.column_stack([grun["value"], grun["capital"], D])

def gvif(Z, block):
    """det(R_11) det(R_22) / det(R) for the columns in block against the rest."""
    R = np.corrcoef(Z, rowvar=False)
    rest = [k for k in range(Z.shape[1]) if k not in block]
    return (np.linalg.det(R[np.ix_(block, block)]) * np.linalg.det(R[np.ix_(rest, rest)])
            / np.linalg.det(R))

firm = list(range(2, Zg.shape[1]))
print("GVIF value  :", round(gvif(Zg, [0]), 2))
print("GVIF capital:", round(gvif(Zg, [1]), 2))
print("GVIF firm   :", round(gvif(Zg, firm), 2), " per dimension:", round(gvif(Zg, firm) ** (1 / 20), 3))
print("VIFs of the firm indicators:", np.round(vifs(Zg)[2:], 1))
# <</gvif>>

g_value, g_cap, g_firm = gvif(Zg, [0]), gvif(Zg, [1]), gvif(Zg, firm)
q = len(firm)
# q = 1: the GVIF is the ordinary VIF
assert np.isclose(g_value, vifs(Zg)[0]) and np.isclose(g_cap, vifs(Zg)[1])
# canonical correlations between the centred blocks
Zc = Zg - Zg.mean(0)
Q1 = np.linalg.qr(Zc[:, firm])[0]
Q2 = np.linalg.qr(Zc[:, :2])[0]
rho = np.linalg.svd(Q1.T @ Q2, compute_uv=False)
assert np.isclose(g_firm, np.prod(1 / (1 - rho**2)))
# determinant form: det Z1'Z1 / det Z1'(I - M2)Z1 for centred columns
M2 = Q2 @ Q2.T
Z1 = Zc[:, firm]
assert np.isclose(g_firm, np.linalg.det(Z1.T @ Z1) / np.linalg.det(Z1.T @ (Z1 - M2 @ Z1)))
# invariance to the coding of the factor: a different reference level and sum-to-zero coding
D_last = pd.get_dummies(grun["firm"]).to_numpy(dtype=float)[:, :-1]
full = pd.get_dummies(grun["firm"]).to_numpy(dtype=float)
D_sum = full[:, :-1] - full[:, [-1]]
for Dalt in (D_last, D_sum):
    Zalt = np.column_stack([grun["value"], grun["capital"], Dalt])
    assert np.isclose(gvif(Zalt, firm), g_firm)
    assert np.isclose(gvif(Zalt, [0]), g_value)
vif_dummies = vifs(Zg)[2:]
vif_dummies_last = vifs(np.column_stack([grun["value"], grun["capital"], D_last]))[2:]
assert not np.isclose(vif_dummies.max(), vif_dummies_last.max())
gen.num("g_value", g_value, 2)
gen.num("g_capital", g_cap, 2)
gen.num("g_firm", g_firm, 2)
gen.num("g_firm_adj", g_firm ** (1 / (2 * q)), 3)
gen.num("g_value_adj", np.sqrt(g_value), 3)
gen.num("rho_max", rho.max(), 4)
gen.num("rho_min", rho.min(), 4)
gen.num("dummy_vif_max", vif_dummies.max(), 1)
gen.num("dummy_vif_max_last", vif_dummies_last.max(), 1)
gen.num("dummy_vif_min", vif_dummies.min(), 1)
gen.int("grun_n", len(Zg))
gen.write()
