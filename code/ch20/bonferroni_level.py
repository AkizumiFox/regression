"""Chapter 20, Section 3: the Bonferroni outlier test at the design of the state data.

Simulated normal errors at the design of the 50 states plus the District of Columbia: how
often does max |t_i| exceed the Bonferroni critical value when there are no outliers, and how
does the power to detect one shifted case depend on its leverage?
"""
import numpy as np
import statsmodels.api as sm
from scipy import stats

from regbook import Generated

data = sm.datasets.statecrime.load_pandas().data
X = np.column_stack([np.ones(len(data)), data[["poverty", "single", "urban"]]])
names = np.array(data.index)


def max_abs_t(X, Y):
    """max_i |t_i| for each column of Y (one simulated response vector per column)."""
    n, p = X.shape
    Q, _ = np.linalg.qr(X)
    h = np.sum(Q ** 2, axis=1)[:, None]
    E = Y - Q @ (Q.T @ Y)
    sse = np.sum(E ** 2, axis=0)
    s2_del = (sse - E ** 2 / (1 - h)) / (n - p - 1)
    return np.max(np.abs(E) / np.sqrt(s2_del * (1 - h)), axis=0)


# <<level>>
n, p = X.shape
alpha = 0.05
crit = stats.t.isf(alpha / (2 * n), n - p - 1)
rng = np.random.default_rng(51)
reps = 4000
T = max_abs_t(X, rng.normal(size=(n, reps)))          # no outliers: errors N(0, 1)
print(f"Bonferroni bound {alpha}, simulated familywise level {np.mean(T > crit):.4f}")
# <</level>>

reps_big = 200_000
T_big = np.concatenate([max_abs_t(X, rng.normal(size=(n, 20_000))) for _ in range(reps_big // 20_000)])
level = np.mean(T_big > crit)
se = np.sqrt(level * (1 - level) / reps_big)
assert level <= alpha and level > alpha - 0.006

# laws of the studentized residuals of one case, here the District of Columbia
Q, _ = np.linalg.qr(X)
h = np.sum(Q ** 2, axis=1)
dc = int(np.where(names == "District of Columbia")[0][0])
Y = rng.normal(size=(n, 50_000))
E = Y - Q @ (Q.T @ Y)
sse = np.sum(E ** 2, axis=0)
r_dc = E[dc] / np.sqrt(sse / (n - p) * (1 - h[dc]))
t_dc0 = E[dc] / np.sqrt((sse - E[dc] ** 2 / (1 - h[dc])) / (n - p - 1) * (1 - h[dc]))
assert stats.kstest(r_dc ** 2 / (n - p), stats.beta(0.5, (n - p - 1) / 2).cdf).pvalue > 0.001
assert stats.kstest(t_dc0, stats.t(n - p - 1).cdf).pvalue > 0.001
assert abs(np.mean(r_dc ** 2) - 1) < 0.02

# power of the single-case test at the Bonferroni level for a shift of 4 sigma
low = int(np.argmin(h))
shift = 4.0
power = {}
for key, i in [("dc", dc), ("low", low)]:
    delta = shift * np.sqrt(1 - h[i])                 # noncentrality of t_i under the shift
    power[key] = stats.nct.sf(crit, n - p - 1, delta) + stats.nct.cdf(-crit, n - p - 1, delta)
# simulation check of the noncentral t law for the District of Columbia
Y = rng.normal(size=(n, 100_000))
Y[dc] += shift
E = Y - Q @ (Q.T @ Y)
sse = np.sum(E ** 2, axis=0)
t_dc = E[dc] / np.sqrt((sse - E[dc] ** 2 / (1 - h[dc])) / (n - p - 1) * (1 - h[dc]))
assert abs(np.mean(np.abs(t_dc) > crit) - power["dc"]) < 0.01
assert stats.kstest(t_dc, stats.nct(n - p - 1, shift * np.sqrt(1 - h[dc])).cdf).pvalue > 0.001
assert power["low"] > power["dc"]

gen = Generated("ch20", "bonferroni_level", prefix="bl")
gen.int("reps", reps_big)
gen.num("level", level, 4)
gen.num("se", se, 4)
gen.num("h_dc", h[dc], 3)
gen.num("h_low", h[low], 3)
gen.text("name_low", names[low])
gen.num("delta_dc", shift * np.sqrt(1 - h[dc]), 2)
gen.num("delta_low", shift * np.sqrt(1 - h[low]), 2)
gen.num("power_dc", power["dc"], 3)
gen.num("power_low", power["low"], 3)
gen.write()
