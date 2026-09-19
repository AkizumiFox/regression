"""Chapter 19, Section 5: a straight line fitted to a curved mean with a random regressor.

X ~ Uniform(0, 2), E(Y | X) = X^2, Var(Y | X) = sigma^2 = 0.04 (normal errors).
The population least squares line has slope beta* = 2 and intercept -2/3. The usual standard error
estimates E(u^2)/Var(X) / n, the true asymptotic variance of the slope is E[(X-1)^2 u^2]/Var(X)^2 / n.
Coverage of the usual 95% interval, for beta* and for the design-dependent target beta_X.
"""
import numpy as np
from scipy import stats

from regbook import Generated

gen = Generated("ch19", "random_regressors", prefix="rx")

# <<theory>>
from fractions import Fraction as Fr
sigma2 = Fr(1, 25)                                   # sigma = 0.2
EU2, EU4, EU6 = Fr(1, 3), Fr(1, 5), Fr(1, 7)         # moments of U = X - 1 ~ Uniform(-1, 1)
Eu2 = EU4 - Fr(2, 3) * EU2 + Fr(1, 9) + sigma2       # u = U^2 - 1/3 + error
sandwich = (EU6 - Fr(2, 3) * EU4 + Fr(1, 9) * EU2 + sigma2 * EU2) / EU2**2
usual = Eu2 / EU2
print("n Var(slope): sandwich", float(sandwich), " usual", float(usual))
z = stats.norm.ppf(0.975)
print("limiting coverage of the usual interval:",
      2 * stats.norm.cdf(z * np.sqrt(float(usual / sandwich))) - 1)
# <</theory>>
assert sandwich == Fr(396, 945) + 3 * sigma2
limit_cov = 2 * stats.norm.cdf(z * np.sqrt(float(usual / sandwich))) - 1
gen.num("sandwich", float(sandwich), 4)
gen.num("usual", float(usual), 4)
gen.num("se_ratio", float(np.sqrt(float(sandwich / usual))), 3)
gen.num("limit_cov", limit_cov, 3)
gen.num("Eu2", float(Eu2), 4)

# <<simulate>>
def one_run(n, reps, rng):
    x = rng.uniform(0, 2, size=(reps, n))
    y = x**2 + 0.2 * rng.normal(size=(reps, n))
    xc = x - x.mean(axis=1, keepdims=True)
    Sxx = np.sum(xc**2, axis=1)
    b = np.sum(xc * y, axis=1) / Sxx                          # least squares slope
    b_X = np.sum(xc * x**2, axis=1) / Sxx                     # its conditional mean given X
    a = y.mean(axis=1) - b * x.mean(axis=1)
    s2 = np.sum((y - a[:, None] - b[:, None] * x) ** 2, axis=1) / (n - 2)
    half = stats.t.ppf(0.975, n - 2) * np.sqrt(s2 / Sxx)
    return np.mean(np.abs(b - 2.0) <= half), np.mean(np.abs(b - b_X) <= half), n * b.var()

rng = np.random.default_rng(1905)
for n in [50, 400]:
    c_star, c_X, nvar = one_run(n, 4000, rng)
    print(f"n = {n}: covers beta* {c_star:.3f}, covers beta_X {c_X:.3f}, n Var(slope) {nvar:.3f}")
# <</simulate>>
rng = np.random.default_rng(1905)
out = {}
for n in [50, 400, 1600]:
    out[n] = one_run(n, 40_000 if n < 1600 else 10_000, rng)
    gen.num(f"cov_star_{n}", out[n][0], 3)
    gen.num(f"cov_X_{n}", out[n][1], 3)
    gen.num(f"nvar_{n}", out[n][2], 3)
assert abs(out[1600][0] - limit_cov) < 0.015
assert out[400][1] > 0.97                           # conservative for the conditional target
assert abs(out[1600][2] / float(sandwich) - 1) < 0.05

# the variance decomposition Var(b) = E[sigma^2 / Sxx] + Var(b_X) at n = 50
rng = np.random.default_rng(1906)
n, reps = 50, 200_000
x = rng.uniform(0, 2, size=(reps, n))
y = x**2 + 0.2 * rng.normal(size=(reps, n))
xc = x - x.mean(axis=1, keepdims=True)
Sxx = np.sum(xc**2, axis=1)
b = np.sum(xc * y, axis=1) / Sxx
b_X = np.sum(xc * x**2, axis=1) / Sxx
within = float(sigma2) * np.mean(1 / Sxx)
between = b_X.var()
assert abs(b.var() / (within + between) - 1) < 0.02
gen.num("var_b", b.var(), 5)
gen.num("within", within, 5)
gen.num("between", between, 5)
gen.int("n_dec", n)
gen.int("reps_dec", reps)
gen.write()
