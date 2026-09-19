"""Chapter 18, Section 6: a lost observation in the randomized block experiment of Section 4.

The joint glued with adhesive 2 on board 4 is lost (the strip split while being cut). Exact
least squares on the 29 remaining observations is compared with the classical devices:
Yates's formula for the missing value, the covariate (indicator) method and iteration.
"""
import numpy as np
from scipy import stats

from regbook import Generated

# <<data>>
rng = np.random.default_rng(606)
t, b = 5, 6                                           # adhesives, boards
plan = np.array([rng.permutation(t) for _ in range(b)])   # which strip of each board gets which adhesive
tau = np.array([0.0, 0.40, 0.65, 0.15, 0.80])         # adhesive effects (unknown in practice)
board = rng.normal(0, 0.75, b)                        # board effects
y = np.round(9 + tau[:, None] + board[None, :] + rng.normal(0, 0.3, (t, b)), 2)
# <</data>>
y_lost = y[1, 3]

# <<exact>>
n = t * b
trt = np.repeat(np.arange(t), b)                      # y.ravel() lists adhesive 1 first
blk = np.tile(np.arange(b), t)
X = np.column_stack([np.ones(n), np.eye(t)[trt], np.eye(b)[blk]])   # complete-data design
miss = np.array([1 * b + 3])                          # adhesive 2, board 4
obs = np.setdiff1d(np.arange(n), miss)
yv = y.ravel().copy()
yv[miss] = np.nan                                     # the value is lost

def fit(Xm, yy):
    coef = np.linalg.lstsq(Xm, yy, rcond=None)[0]
    return coef, np.sum((yy - Xm @ coef) ** 2)

coef_obs, sse_obs = fit(X[obs], yv[obs])              # exact least squares, 29 observations
df_obs = len(obs) - np.linalg.matrix_rank(X[obs])     # 29 - 10 = 19
X_red = np.column_stack([np.ones(n), np.eye(b)[blk]])  # blocks only
_, sse_red = fit(X_red[obs], yv[obs])
F_exact = ((sse_red - sse_obs) / (t - 1)) / (sse_obs / df_obs)
print(f"exact: SSE {sse_obs:.3f} on {df_obs} df, adhesive SS {sse_red - sse_obs:.3f},"
      f" F = {F_exact:.2f};  fitted value for the lost cell {X[miss] @ coef_obs}")
# <</exact>>

# <<yates>>
Y_obs = np.where(np.isnan(yv), 0.0, yv).reshape(t, b)
T_i = Y_obs[1].sum()                                  # observed total of adhesive 2
B_j = Y_obs[:, 3].sum()                               # observed total of board 4
G = Y_obs.sum()
z_yates = (t * T_i + b * B_j - G) / ((t - 1) * (b - 1))
print(f"Yates's value {z_yates:.3f}")
# <</yates>>

# <<covariate>>
d = np.zeros(n)
d[miss] = 1.0                                         # one indicator column per missing cell
y_zero = np.where(np.isnan(yv), 0.0, yv)              # any number will do in the gap
coef_cov, sse_cov = fit(np.column_stack([X, d]), y_zero)
print(f"covariate method: gamma-hat = {coef_cov[-1]:.3f}, so the filled value is {-coef_cov[-1]:.3f};"
      f" SSE {sse_cov:.3f}")
# <</covariate>>

# <<iterate>>
z = np.array([np.nanmean(yv)])                        # start from the grand mean
history = [z[0]]
for step in range(30):
    y_fill = yv.copy()
    y_fill[miss] = z
    coef_fill, _ = fit(X, y_fill)                     # complete-data analysis
    z = X[miss] @ coef_fill                           # replace the gap by its fitted value
    history.append(z[0])
print("iterates", np.round(history[:6], 3))
# <</iterate>>

# ---- checks -----------------------------------------------------------------
fitted_miss = (X[miss] @ coef_obs)[0]
assert np.isclose(z_yates, fitted_miss)
assert np.isclose(-coef_cov[-1], fitted_miss)
assert np.isclose(sse_cov, sse_obs)
assert np.allclose(coef_cov[:-1] @ np.linalg.pinv(X) @ X, coef_obs @ np.linalg.pinv(X) @ X)  # same fitted means
assert np.isclose(history[-1], fitted_miss)
M = X @ np.linalg.pinv(X)
m_mm = M[miss[0], miss[0]]
assert np.isclose(m_mm, 1 / b + 1 / t - 1 / (b * t))          # the error shrinks by this factor
errs = np.abs(np.array(history[:8]) - fitted_miss)
assert np.allclose(errs[1:] / errs[:-1], m_mm)

# the filled-in data analysed as if complete
y_fill = yv.copy()
y_fill[miss] = fitted_miss
Yf = y_fill.reshape(t, b)
gm = Yf.mean()
ss_trt_fill = b * np.sum((Yf.mean(axis=1) - gm) ** 2)
sse_fill = np.sum((Yf - Yf.mean(axis=1)[:, None] - Yf.mean(axis=0)[None, :] + gm) ** 2)
assert np.isclose(sse_fill, sse_obs)                  # the residual SS is already right
df_naive = (t - 1) * (b - 1)
F_naive = (ss_trt_fill / (t - 1)) / (sse_fill / df_naive)
bias = ss_trt_fill - (sse_red - sse_obs)
assert bias > 0
assert np.isclose(bias, (B_j - (t - 1) * fitted_miss) ** 2 / (t * (t - 1)))   # Yates's bias formula
assert F_naive > F_exact

# two lost values: iteration still converges to the exact least squares values
miss2 = np.array([1 * b + 3, 4 * b + 0])
obs2 = np.setdiff1d(np.arange(n), miss2)
coef2, sse2 = fit(X[obs2], y.ravel()[obs2])
z2 = np.full(2, y.ravel()[obs2].mean())
for step in range(60):
    y_fill2 = y.ravel().copy()
    y_fill2[miss2] = z2
    z2 = X[miss2] @ fit(X, y_fill2)[0]
assert np.allclose(z2, X[miss2] @ coef2)
M_mm2 = M[np.ix_(miss2, miss2)]
rho2 = np.max(np.abs(np.linalg.eigvals(M_mm2)))
assert rho2 < 1

gen = Generated("ch18", "missing")
gen.num("lost", y_lost, 2)
gen.num("Ti", T_i, 2)
gen.num("Bj", B_j, 2)
gen.num("G", G, 2)
gen.num("z", fitted_miss, 3)
gen.num("h0", history[0], 3)
gen.num("h1", history[1], 3)
gen.num("h2", history[2], 3)
gen.num("h3", history[3], 3)
gen.num("sse", sse_obs, 4)
gen.int("df", df_obs)
gen.num("sstrt", sse_red - sse_obs, 4)
gen.num("sstrtfill", ss_trt_fill, 4)
gen.num("bias", bias, 4)
gen.num("Fexact", F_exact, 2)
gen.num("Fnaive", F_naive, 2)
gen.num("mserr", sse_obs / df_obs, 4)
gen.num("msnaive", sse_obs / df_naive, 4)
gen.num("gamma", coef_cov[-1], 3)
gen.num("rho2", rho2, 3)
gen.write()
