"""Chapter 21, Section 1: what residuals can and cannot reveal about the errors.

A simulated design with n = 30 observations and an intercept plus p - 1 standard normal
regressors, p = 2, 10 or 20. Errors are centred exponential (skewed, excess kurtosis 6) or
t with 5 degrees of freedom. The script computes, for each design, how much of the errors'
skewness and excess kurtosis survives in the residuals (Proposition 21.1.1), and simulates
the power of the Shapiro-Wilk test applied to the errors and to the residuals.
"""
import numpy as np
from scipy import stats

from regbook import Generated

# <<design>>
rng = np.random.default_rng(2101)
n = 30
designs = {}
for p in (2, 10, 20):
    X = np.column_stack([np.ones(n), rng.normal(size=(n, p - 1))])
    Q, _ = np.linalg.qr(X)
    M = Q @ Q.T                                      # the hat matrix
    designs[p] = (M, np.diag(M).copy())
# <</design>>

# <<shrink>>
for p, (M, h) in designs.items():
    # residual i = (row i of C) @ errors
    C = np.eye(n) - M
    # excess kurtosis of residual i / that of an error
    kurt_kept = (C**4).sum(axis=1) / (1 - h) ** 2
    # skewness of residual i / that of an error
    skew_kept = (C**3).sum(axis=1) / (1 - h) ** 1.5
    print(f"p = {p:2d}: max leverage {h.max():.2f}, "
          f"average share of skewness kept {skew_kept.mean():.3f}, "
          f"of excess kurtosis kept {kurt_kept.mean():.3f}")
# <</shrink>>

gen = Generated("ch21", "residual_errors", prefix="rerr")
gen.int("n", n)
for p, (M, h) in designs.items():
    C = np.eye(n) - M
    kk = (C**4).sum(axis=1) / (1 - h) ** 2
    sk = (C**3).sum(axis=1) / (1 - h) ** 1.5
    # Proposition 21.1.1(c): the kurtosis factor lies in [(1 - h_ii)^2, 1]
    assert np.all(kk <= 1 + 1e-12) and np.all(kk >= (1 - h) ** 2 - 1e-12)
    assert np.all(np.abs(sk) <= 1 + 1e-12)
    gen.num(f"kurt_kept_{p}", kk.mean(), 3)
    gen.num(f"skew_kept_{p}", sk.mean(), 3)
    gen.num(f"hmax_{p}", h.max(), 2)


def errors(kind, size, rng):
    """Mean-zero, variance-one errors."""
    if kind == "exp":
        return rng.exponential(size=size) - 1.0
    return rng.standard_t(5, size=size) / np.sqrt(5 / 3)


# ---- check the kurtosis factor by simulation for one residual -------------------
M, h = designs[10]
C = np.eye(n) - M
E = errors("exp", (200_000, n), rng)
R = E @ C.T
i = int(np.argmax(h))
z = R[:, i] / np.sqrt(1 - h[i])
kurt_sim = np.mean(z**4) - 3
kurt_pred = 6 * (C[i] ** 4).sum() / (1 - h[i]) ** 2
assert abs(kurt_sim - kurt_pred) < 0.15, (kurt_sim, kurt_pred)
gen.num("kurt_sim_i", kurt_sim, 2)
gen.num("kurt_pred_i", kurt_pred, 2)
gen.num("h_i", h[i], 2)

# ---- ordered residuals against ordered errors (Proposition 21.1.1(b)) -----------
E = errors("exp", (4000, n), rng)
dist2, bound = [], []
for e in E:
    r = C @ e
    dist2.append(np.sum((np.sort(r) - np.sort(e)) ** 2))
    bound.append(np.sum((M @ e) ** 2))
dist2, bound = np.array(dist2), np.array(bound)
assert np.all(dist2 <= bound + 1e-10)
assert abs(bound.mean() - 10) < 0.3             # E ||M e||^2 = sigma^2 p = 10
gen.num("ordered_mean", dist2.mean() / n, 3)
gen.num("bound_mean", bound.mean() / n, 3)

# <<power>>
def power_table(reps, rng):
    """Rejection rates of the 5% Shapiro-Wilk test on errors and on residuals."""
    table = {}
    for kind in ("exp", "t5"):
        for p, (M, h) in designs.items():
            E = errors(kind, (reps, n), rng)
            # residual vectors (M is symmetric)
            R = E - E @ M
            rej_err = np.mean([stats.shapiro(e).pvalue < 0.05 for e in E])
            rej_res = np.mean([stats.shapiro(r).pvalue < 0.05 for r in R])
            table[kind, p] = (rej_err, rej_res)
            print(f"{kind:3s} p = {p:2d}: Shapiro-Wilk rejects "
                  f"errors {rej_err:.2f}, "
                  f"residuals {rej_res:.2f}")
    return table


# a quick version (the book uses 4000 replications)
quick = power_table(400, rng)
# <</power>>

table = power_table(4000, rng)
for (kind, p), (rej_err, rej_res) in table.items():
    gen.num(f"pow_err_{kind}_{p}", rej_err, 2)
    gen.num(f"pow_res_{kind}_{p}", rej_res, 2)
    if p == 20:
        assert rej_res < rej_err - 0.1

# ---- size of Shapiro-Wilk applied to residuals under normal errors -------------
M, h = designs[20]
reps = 4000
E = rng.normal(size=(reps, n))
R = E - E @ M
S = np.sqrt((R**2).sum(axis=1) / (n - 20))
size_raw = np.mean([stats.shapiro(r).pvalue < 0.05 for r in R])
size_stud = np.mean([stats.shapiro(r / (s * np.sqrt(1 - h))).pvalue < 0.05 for r, s in zip(R, S)])
assert abs(size_raw - 0.05) < 0.02 and abs(size_stud - 0.05) < 0.02
gen.num("size_raw_20", size_raw, 3)
gen.num("size_stud_20", size_stud, 3)
gen.write()
