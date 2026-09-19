"""Chapter 24, Section 3: method-of-moments correction with replicate measurements.

Blood-pressure study: regress y on (1, mean of two readings, age). The per-patient matrix
D_i = diag(0, (w_i1 - w_i2)^2 / 4, 0) has expectation Cov(error of the mean reading), so
theta = (sum (w_i w_i^T - D_i))^{-1} sum w_i y_i is consistent, and the sandwich
A^{-1} B A^{-1} / n built from psi_i = w_i (y_i - w_i^T theta) + D_i theta gives its standard
errors, including the cost of estimating the error variance.

Repeated samples compare the coverage of the corrected and the naive 95% intervals.
"""
from scipy import stats

from regbook import Generated

gen = Generated("ch24", "moment_correction", prefix="mom")

# <<setup>>
import numpy as np

rng = np.random.default_rng(2401)
n = 300
age = rng.normal(55, 10, n)
x = 130 + 0.6 * (age - 55) + rng.normal(0, np.sqrt(108), n)   # long-run blood pressure, sd 12
w = x[:, None] + rng.normal(0, 9, (n, 2))                      # two clinic readings
y = -4 + 0.08 * x + 0.0 * age + rng.normal(0, 1.0, n)          # age has no effect given x
# <</setup>>


# <<correct>>
def corrected(w, age, y):
    """Moment-corrected coefficients of (1, x, age) and their sandwich standard errors."""
    n = len(y)
    W = np.column_stack([np.ones(n), w.mean(axis=1), age])   # observed regressors
    D = np.zeros((n, 3, 3))
    D[:, 1, 1] = (w[:, 0] - w[:, 1]) ** 2 / 4                  # E D_i = Var(error of the mean)
    A = (W.T @ W - D.sum(axis=0)) / n
    theta = np.linalg.solve(A, W.T @ y / n)
    psi = W * (y - W @ theta)[:, None] + D @ theta               # estimating function, one row per patient
    B = psi.T @ psi / n
    Ainv = np.linalg.inv(A)
    cov = Ainv @ B @ Ainv / n
    return theta, np.sqrt(np.diag(cov))


theta, se = corrected(w, age, y)
W = np.column_stack([np.ones(n), w.mean(axis=1), age])
naive = np.linalg.lstsq(W, y, rcond=None)[0]
s2_u = np.mean((w[:, 0] - w[:, 1]) ** 2) / 2                   # estimate of sigma_u^2 (true 81)
print(f"sigma_u^2 estimate {s2_u:.1f}")
print("naive     (1, x, age):", naive.round(4))
print("corrected (1, x, age):", theta.round(4), " se:", se.round(4))
# <</correct>>

# the naive least squares standard errors, for comparison
res = y - W @ naive
naive_se = np.sqrt(res @ res / (n - 3) * np.diag(np.linalg.inv(W.T @ W)))
# with the error variance treated as known (81), the correction uses D_i = diag(0, 81/2, 0)
Wn = W.T @ W / n
theta_known = np.linalg.solve(Wn - np.diag([0, 81 / 2, 0]), W.T @ y / n)
# the corrected slope agrees with the ratio naive / (estimated reliability given age) only
# approximately; check the exact identity for the known-variance, no-covariate case instead
Wc = np.column_stack([np.ones(n), w.mean(axis=1)])
b_naive1 = np.linalg.lstsq(Wc, y, rcond=None)[0][1]
lam_hat = 1 - (s2_u / 2) / np.var(w.mean(axis=1))
theta1 = np.linalg.solve(Wc.T @ Wc / n - np.diag([0, s2_u / 2]), Wc.T @ y / n)
assert np.isclose(theta1[1], b_naive1 / lam_hat)

assert abs(theta[1] - 0.08) < 2.5 * se[1]
assert abs(theta[2]) < 2 * se[2]
assert naive[1] < theta[1]
assert se[1] > naive_se[1]
gen.num("s2_u", s2_u, 1)
gen.num("lam_bar", 144 / (144 + 81 / 2), 3)
gen.num("lam_bar_z", 108 / (108 + 81 / 2), 3)
gen.num("naive_x", naive[1], 4)
gen.num("naive_x_se", naive_se[1], 4)
gen.num("naive_age", naive[2], 4)
gen.num("naive_age_se", naive_se[2], 4)
gen.num("corr_x", theta[1], 4)
gen.num("corr_x_se", se[1], 4)
gen.num("corr_age", theta[2], 4)
gen.num("corr_age_se", se[2], 4)
gen.num("corr_x_lo", theta[1] - 1.96 * se[1], 4)
gen.num("corr_x_hi", theta[1] + 1.96 * se[1], 4)
gen.num("known_x", theta_known[1], 4)
gen.num("se_ratio", se[1] / naive_se[1], 2)


# ---- repeated samples -----------------------------------------------------------------------
# <<coverage>>
def coverage(reps, seed):
    """Coverage of nominal 95% intervals for beta_x = 0.08: corrected (sandwich) and naive."""
    r = np.random.default_rng(seed)
    hits_c = hits_n = 0
    for _ in range(reps):
        a = r.normal(55, 10, n)
        xx = 130 + 0.6 * (a - 55) + r.normal(0, np.sqrt(108), n)
        ww = xx[:, None] + r.normal(0, 9, (n, 2))
        yy = -4 + 0.08 * xx + r.normal(0, 1.0, n)
        th, s = corrected(ww, a, yy)
        hits_c += abs(th[1] - 0.08) <= 1.96 * s[1]
        Wr = np.column_stack([np.ones(n), ww.mean(axis=1), a])
        b = np.linalg.lstsq(Wr, yy, rcond=None)[0]
        e = yy - Wr @ b
        sb = np.sqrt(e @ e / (n - 3) * np.linalg.inv(Wr.T @ Wr)[1, 1])
        hits_n += abs(b[1] - 0.08) <= 1.96 * sb
    return hits_c / reps, hits_n / reps
# <</coverage>>


reps = 10_000
cov_c, cov_n = coverage(reps, 2420)
# <<coverage-small>>
print("coverage (corrected, naive):", coverage(500, 2421))
# <</coverage-small>>
se_mc = np.sqrt(0.95 * 0.05 / reps)
assert abs(cov_c - 0.95) < 4 * se_mc + 0.005
assert cov_n < 0.3
gen.int("reps", reps)
gen.num("cov_corr", cov_c, 3)
gen.num("cov_naive", cov_n, 3)
gen.write()
