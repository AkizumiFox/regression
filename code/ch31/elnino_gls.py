"""Chapter 31, Sections 1 and 2: generalized least squares on the El Nino sea surface
temperatures, with a random annual shock (equicorrelation inside each year).

Monthly averaged sea surface temperature of the Pacific between 0-10 degrees South and
90-80 degrees West, 1950-2010: 61 complete years of 12 months. Public-domain data from the
US National Oceanic and Atmospheric Administration, shipped with statsmodels
(statsmodels.datasets.elnino).
"""
import matplotlib.pyplot as plt
import numpy as np
import statsmodels.api as sm
from scipy import stats

from regbook import COLORS, Generated, figure_path, use_book_style

# <<data>>
import numpy as np
import statsmodels.api as sm
from scipy import stats

frame = sm.datasets.elnino.load_pandas().data
years = frame["YEAR"].to_numpy()
temps = frame.iloc[:, 1:].to_numpy()                  # 61 years by 12 months
G, m = temps.shape
n = G * m
y = temps.ravel()                                     # the month varies fastest
year = np.repeat(years, m)
month = np.tile(np.arange(1, m + 1), G)
time = (year - 1980.0) + (month - 6.5) / 12.0         # decimal years, centred at 1980

cols = [np.ones(n), time / 10.0]                      # the trend is per decade
for k in (1, 2, 3):
    cols += [np.cos(2 * np.pi * k * month / m), np.sin(2 * np.pi * k * month / m)]
X = np.column_stack(cols)
p = X.shape[1]
# <</data>>

# <<gls>>
rho = 0.5                                             # taken as known for now


def gls_fit(X, y, rho, m):
    """GLS with equicorrelation rho inside blocks of m consecutive observations.

    Whitening is the quasi-demeaning y -> (y - c ybar_g)/sqrt(1 - rho) with
    c = 1 - sqrt((1 - rho)/(1 + (m - 1) rho)); ordinary least squares on the whitened
    data is generalized least squares on the original data.
    """
    c = 1.0 - np.sqrt((1 - rho) / (1 + (m - 1) * rho))

    def whiten(a):
        b = a.reshape(len(a) // m, m, -1) if a.ndim > 1 else a.reshape(-1, m)
        return ((b - c * b.mean(axis=1, keepdims=True)) / np.sqrt(1 - rho)).reshape(a.shape)

    Xw, yw = whiten(X), whiten(y)
    b = np.linalg.solve(Xw.T @ Xw, Xw.T @ yw)
    rss = float((yw - Xw @ b) @ (yw - Xw @ b))         # the V^{-1} residual sum of squares
    s2 = rss / (len(y) - X.shape[1])
    return b, s2 * np.linalg.inv(Xw.T @ Xw), rss, c


b_ols = np.linalg.solve(X.T @ X, X.T @ y)
rss_ols = float((y - X @ b_ols) @ (y - X @ b_ols))
V_ols = rss_ols / (n - p) * np.linalg.inv(X.T @ X)     # the standard errors usually reported
b_gls, V_gls, rss_gls, cq = gls_fit(X, y, rho, m)

print(f"trend per decade: OLS {b_ols[1]:+.4f}  GLS {b_gls[1]:+.4f}")
print(f"standard error:   OLS {np.sqrt(V_ols[1, 1]):.4f}  GLS {np.sqrt(V_gls[1, 1]):.4f}")
print(f"quasi-demeaning constant c = {cq:.4f},  s_GLS = {np.sqrt(rss_gls / (n - p)):.4f}")
# <</gls>>

# the whitened fit really is the GLS fit computed from V^{-1} directly
Z = np.zeros((n, G))
Z[np.arange(n), np.repeat(np.arange(G), m)] = 1.0      # the year indicator matrix
V = (1 - rho) * np.eye(n) + rho * Z @ Z.T
Vinv = np.linalg.inv(V)
assert np.allclose(b_gls, np.linalg.solve(X.T @ Vinv @ X, X.T @ Vinv @ y))
assert np.allclose(V_gls, rss_gls / (n - p) * np.linalg.inv(X.T @ Vinv @ X))
resid_gls = y - X @ b_gls
assert np.isclose(rss_gls, resid_gls @ Vinv @ resid_gls)
# the Woodbury form of the inverse
Vinv_wb = (np.eye(n) - rho / (1 + (m - 1) * rho) * Z @ Z.T) / (1 - rho)
assert np.allclose(Vinv, Vinv_wb)
# the GLS residual is V^{-1}-orthogonal to C(X) but not orthogonal to it
assert np.allclose(X.T @ Vinv @ resid_gls, 0, atol=1e-8)
assert np.abs(X.T @ resid_gls).max() > 1e-3
# and its ordinary sum of squares exceeds that of the OLS residual
assert resid_gls @ resid_gls > rss_ols


# <<ftest>>
def f_test(X, y, drop, rho, m):
    """F test that the coefficients listed in `drop` are zero, in the rho-geometry."""
    keep = [j for j in range(X.shape[1]) if j not in drop]
    full, null = gls_fit(X, y, rho, m), gls_fit(X[:, keep], y, rho, m)
    df1, df2 = len(drop), len(y) - X.shape[1]
    F = (null[2] - full[2]) / df1 / (full[2] / df2)
    return F, stats.f.sf(F, df1, df2)


F_trend_gls, p_trend_gls = f_test(X, y, [1], rho, m)
F_trend_ols, p_trend_ols = f_test(X, y, [1], 0.0, m)   # rho = 0 is ordinary least squares
F_harm_gls, p_harm_gls = f_test(X, y, [4, 5, 6, 7], rho, m)
F_harm_ols, p_harm_ols = f_test(X, y, [4, 5, 6, 7], 0.0, m)
print(f"trend:      F = {F_trend_ols:7.2f} (rho = 0)   F = {F_trend_gls:6.2f} (rho = 0.5)")
print(f"harmonics:  F = {F_harm_ols:7.2f} (rho = 0)   F = {F_harm_gls:6.2f} (rho = 0.5)")
# <</ftest>>

assert np.isclose(F_trend_gls, (b_gls[1] / np.sqrt(V_gls[1, 1])) ** 2)
assert np.isclose(F_trend_ols, (b_ols[1] / np.sqrt(V_ols[1, 1])) ** 2)

# ---- Section 2: what the ordinary fit loses, and what it gets wrong ----------
M = X @ np.linalg.solve(X.T @ X, X.T)
Cov_ols = np.linalg.solve(X.T @ X, X.T @ V @ X) @ np.linalg.inv(X.T @ X)   # sigma^2 = 1
Cov_gls = np.linalg.inv(X.T @ Vinv @ X)
eff = np.array([Cov_gls[j, j] / Cov_ols[j, j] for j in range(p)])
Es2 = np.trace((np.eye(n) - M) @ V) / (n - p)          # E(s^2) under the true covariance
Ginv = np.linalg.inv(X.T @ X)
se_ratio = np.sqrt(np.array([Es2 * Ginv[j, j] / Cov_ols[j, j] for j in range(p)]))
kruskal = np.abs((np.eye(n) - M) @ V @ X).max()        # zero iff C(VX) is inside C(X)
print(f"efficiency of OLS: trend {eff[1]:.4f}, first cosine {eff[2]:.6f}")
print(f"reported/true standard error: trend {se_ratio[1]:.3f}, cosine {se_ratio[2]:.3f}")

assert 0.999 < eff[1] < 1.0                    # nearly, but not exactly, best
assert abs(eff[0] - 1) < 1e-12                 # the intercept is exactly best
assert kruskal > 0.1                           # C(VX) is not inside C(X)
assert se_ratio[1] < 0.45 and se_ratio[2] > 1.4
# the intercept is exactly best because 1 is orthogonal to every other column
assert np.allclose(X[:, 0] @ X[:, 1:], 0, atol=1e-9)
assert np.allclose((np.eye(n) - M) @ V @ X @ Ginv[:, 0], 0, atol=1e-12)

# the saturated seasonal model: eleven month effects instead of three harmonics
D = np.zeros((n, m - 1))
for j in range(m - 1):
    D[month == j + 1, j] = 1.0
Xs = np.column_stack([np.ones(n), time / 10.0, D])
Ms = Xs @ np.linalg.solve(Xs.T @ Xs, Xs.T)
assert np.allclose((np.eye(n) - Ms) @ V @ Xs, 0, rtol=0, atol=1e-9)   # Kruskal's condition holds
b_sat_ols = np.linalg.solve(Xs.T @ Xs, Xs.T @ y)
b_sat_gls = np.linalg.solve(Xs.T @ Vinv @ Xs, Xs.T @ Vinv @ y)
assert np.allclose(b_sat_ols, b_sat_gls, rtol=0, atol=1e-9)         # the two fits coincide
assert np.abs(b_ols[2:] - b_gls[2:]).max() < 1e-4                   # but only nearly, with 3 harmonics
assert np.abs(b_ols[2:] - b_gls[2:]).max() > 1e-6
Cov_sat = np.linalg.solve(Xs.T @ Xs, Xs.T @ V @ Xs) @ np.linalg.inv(Xs.T @ Xs)
Gs = np.linalg.inv(Xs.T @ Xs)
Es2_sat = np.trace((np.eye(n) - Ms) @ V) / (n - Xs.shape[1])
sat_factor = np.sqrt(Es2_sat * Gs[1, 1] / Cov_sat[1, 1])      # reported / true s.e.
assert sat_factor < 0.5

# ---- the value of rho the data themselves suggest ---------------------------
e_ols = y - X @ b_ols
e_blocks = e_ols.reshape(G, m)
between = float(np.sum(e_blocks.mean(axis=1) ** 2)) * m / G
within = float(np.sum((e_blocks - e_blocks.mean(axis=1, keepdims=True)) ** 2)) / (G * (m - 1))
rho_hat = (between - within) / (between + (m - 1) * within)
print(f"moment estimate of rho from the OLS residuals: {rho_hat:.3f}")

gen = Generated("ch31", "elnino_gls", prefix="nino")
gen.int("n", n)
gen.int("G", G)
gen.int("m", m)
gen.int("p", p)
gen.num("rho", rho, 1)
gen.num("c", cq, 3)
gen.num("trend_ols", b_ols[1], 4)
gen.num("trend_gls", b_gls[1], 4)
gen.num("se_ols", np.sqrt(V_ols[1, 1]), 4)
gen.num("se_gls", np.sqrt(V_gls[1, 1]), 4)
gen.num("se_factor", np.sqrt(V_gls[1, 1] / V_ols[1, 1]), 2)
gen.num("s_gls", np.sqrt(rss_gls / (n - p)), 3)
gen.num("s_ols", np.sqrt(rss_ols / (n - p)), 3)
gen.num("rss_gls", rss_gls, 1)
gen.num("cos1_ols", b_ols[2], 4)
gen.num("cos1_gls", b_gls[2], 4)
gen.num("sin1_ols", b_ols[3], 4)
gen.num("sin1_gls", b_gls[3], 4)
gen.num("amp1", np.hypot(b_ols[2], b_ols[3]), 3)
gen.num("se_sin1_ols", np.sqrt(V_ols[3, 3]), 4)
gen.num("se_sin1_gls", np.sqrt(V_gls[3, 3]), 4)
gen.num("F_trend_ols", F_trend_ols, 1)
gen.num("F_trend_gls", F_trend_gls, 2)
gen.num("p_trend_gls", p_trend_gls, 4)
gen.num("F_harm_ols", F_harm_ols, 1)
gen.num("F_harm_gls", F_harm_gls, 1)
gen.num("eff_trend", eff[1], 4)
gen.num("eff_cos1", eff[2], 4)
gen.num("se_ratio_trend", se_ratio[1], 3)
gen.num("se_ratio_cos1", se_ratio[2], 3)
gen.num("kruskal", kruskal, 3)
gen.num("sat_factor", sat_factor, 3)
gen.num("trend_sat", b_sat_ols[1], 4)
gen.num("rho_hat", rho_hat, 3)
gen.write()

# ---- Figure 31.1.1 -----------------------------------------------------------
use_book_style()
fig, axes = plt.subplots(1, 2, figsize=(5.8, 2.4))

ax = axes[0]
ybar = e_blocks.mean(axis=1)
ax.axhline(0, color=COLORS["grid"], linewidth=0.6, zorder=0)
ax.vlines(years, 0, ybar, color=COLORS["accent"], linewidth=1.1)
sd_ind = np.sqrt(rss_ols / (n - p) / m)
for s in (-2, 2):
    ax.axhline(s * sd_ind, color=COLORS["second"], linewidth=0.8, linestyle="--")
ax.set_xlabel("year")
ax.set_ylabel(r"mean residual ($^\circ$C)")
ax.set_title("(a) annual means of the OLS residuals")

ax = axes[1]
ts = np.linspace(time.min(), time.max(), 200)
L = np.column_stack([np.ones_like(ts), ts / 10.0])     # the trend part of the fitted mean
for b, C, colour, lab in ((b_gls, V_gls, COLORS["accent"], "GLS band"),
                          (b_ols, V_ols, COLORS["second"], "OLS, usual band")):
    centre = L @ b[:2]
    half = 1.96 * np.sqrt(np.einsum("ij,jk,ik->i", L, C[:2, :2], L))
    ax.fill_between(ts + 1980, centre - half, centre + half, color=colour, alpha=0.25,
                    linewidth=0, label=lab)
ax.plot(ts + 1980, L @ b_gls[:2], color=COLORS["ink"], linewidth=1.0)
ax.set_xlabel("year")
ax.set_ylabel(r"trend component ($^\circ$C)")
ax.set_title("(b) the estimated trend")
ax.legend(loc="upper left", frameon=False)
fig.tight_layout()
fig.savefig(figure_path("ch31", "elnino_gls"))

# ---- Figure 31.2.1 -----------------------------------------------------------
grid = np.linspace(0.0, 0.9, 91)
curves = {"eff": [], "trend": [], "cos": []}
for r in grid:
    Vr = (1 - r) * np.eye(n) + r * Z @ Z.T
    cov_o = np.linalg.solve(X.T @ X, X.T @ Vr @ X) @ Ginv
    cov_g = np.linalg.inv(X.T @ np.linalg.solve(Vr, X))
    es2 = np.trace((np.eye(n) - M) @ Vr) / (n - p)
    curves["eff"].append(cov_g[1, 1] / cov_o[1, 1])
    curves["trend"].append(np.sqrt(es2 * Ginv[1, 1] / cov_o[1, 1]))
    curves["cos"].append(np.sqrt(es2 * Ginv[2, 2] / cov_o[2, 2]))

fig2, ax = plt.subplots(figsize=(3.6, 2.4))
ax.plot(grid, curves["trend"], color=COLORS["second"], label="trend, reported/true s.e.")
ax.plot(grid, curves["cos"], color=COLORS["third"], label="cosine, reported/true s.e.")
ax.plot(grid, curves["eff"], color=COLORS["accent"], label="efficiency of OLS (trend)")
ax.axhline(1, color=COLORS["grid"], linewidth=0.6, zorder=0)
ax.set_xlabel(r"$\rho$")
ax.set_ylim(0, 3.3)
ax.legend(loc="upper left", frameon=False)
fig2.tight_layout()
fig2.savefig(figure_path("ch31", "elnino_efficiency"))

assert np.isclose(curves["eff"][0], 1.0) and min(curves["eff"]) > 0.99
assert curves["trend"][-1] < 0.31 and curves["cos"][-1] > 3.0
