"""Chapter 23, Section 3: bootstrap confidence intervals for the poverty coefficient of Chapter 6.

The 2009 US state data (statsmodels.datasets.statecrime, public domain), District of Columbia
omitted as in Chapter 6: murder rate on poverty, single-parent percentage and urban percentage.
Seven 95% intervals for the poverty coefficient.
"""
import numpy as np
import statsmodels.api as sm
from scipy import stats

from regbook import Generated

data = sm.datasets.statecrime.load_pandas().data
data = data.drop(index="District of Columbia")

# <<intervals>>
rng = np.random.default_rng(2303)
B = 9999
y = data["murder"].to_numpy()
X = np.column_stack([np.ones(len(y)), data[["poverty", "single", "urban"]]])
n, p = X.shape
XtX_inv = np.linalg.inv(X.T @ X)
A = XtX_inv @ X.T
a, c = A[1], np.sqrt(XtX_inv[1, 1])               # the poverty coefficient is a @ y
h = np.sum((X @ XtX_inv) * X, axis=1)
fit = X @ (A @ y)
e = y - fit
b = a @ y
s = np.sqrt(e @ e / (n - p))
q = [0.025, 0.975]
out = {"t": b + np.array([-1, 1]) * stats.t.ppf(0.975, n - p) * s * c}

r = e / np.sqrt(1 - h)                            # residual bootstrap, modified residuals
r -= r.mean()
Ystar = fit + r[rng.integers(0, n, size=(B, n))]
bstar = Ystar @ a
sstar = np.sqrt(np.sum((Ystar - (Ystar @ A.T) @ X.T) ** 2, axis=1) / (n - p))
tstar = (bstar - b) / (sstar * c)
lo, hi = np.quantile(bstar, q)
out["percentile"] = np.array([lo, hi])
out["basic"] = np.array([2 * b - hi, 2 * b - lo])
out["studentized"] = b - np.quantile(tstar, q[::-1]) * s * c

z0 = stats.norm.ppf(np.mean(bstar < b))           # BCa
loo = b - a * e / (1 - h)                         # case-jackknife values (heuristic with residual resampling)
dev = loo.mean() - loo
acc = np.sum(dev ** 3) / (6 * np.sum(dev ** 2) ** 1.5)
zs = z0 + stats.norm.ppf(q)
out["BCa"] = np.quantile(bstar, stats.norm.cdf(z0 + zs / (1 - acc * zs)))

idx = rng.integers(0, n, size=(B, n))             # case bootstrap
bcase = np.array([np.linalg.lstsq(X[i], y[i], rcond=None)[0][1] for i in idx])
out["case"] = np.quantile(bcase, q)

se2 = np.sqrt(np.sum(a ** 2 * e ** 2 / (1 - h)))  # wild bootstrap-t with HC2
Ystar = fit + (e / np.sqrt(1 - h)) * rng.choice([-1.0, 1.0], size=(B, n))
Estar = Ystar - (Ystar @ A.T) @ X.T
wstar = (Ystar @ a - b) / np.sqrt((Estar ** 2 / (1 - h)) @ a ** 2)
out["wild-t"] = b - np.quantile(wstar, q[::-1]) * se2

for m, (l, u) in out.items():
    print(f"{m:12s} ({l:.3f}, {u:.3f})   length {u - l:.3f}")
# <</intervals>>

# ---- checks -------------------------------------------------------------------------------
assert np.isclose(b, 0.2538, atol=5e-5)            # the coefficient of Chapter 6
# deletion formula against brute force
for i in [0, 7, 23]:
    keep = np.arange(n) != i
    assert np.isclose(np.linalg.lstsq(X[keep], y[keep], rcond=None)[0][1], loo[i])
# every interval contains the estimate and excludes zero
for m, (l, u) in out.items():
    assert l < b < u and l > 0
# the studentized interval is close to the t interval (same studentization, near-normal residuals)
assert np.max(np.abs(out["studentized"] - out["t"])) < 0.02
assert se2 < s * c                                  # the sandwich se is the smaller one here
lengths = {m: u - l for m, (l, u) in out.items()}
assert lengths["percentile"] < lengths["t"]         # the percentile interval is a little short

gen = Generated("ch23", "crime_intervals", prefix="crime")
gen.int("B", B)
gen.num("b", b, 4)
gen.num("z0", z0, 3)
gen.num("acc", acc, 4)
gen.num("se_classical", s * c, 4)
gen.num("se_hc2", se2, 4)
for m, (l, u) in out.items():
    key = m.replace("-", "")
    gen.num(f"{key}:lo", l, 3)
    gen.num(f"{key}:hi", u, 3)
    gen.num(f"{key}:len", u - l, 3)
gen.write()
