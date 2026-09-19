"""Chapter 11, Section 4: likelihood ratio, Wald and score statistics for the region test.

All three are increasing functions of F, so with exact calibration they are the same test.
Calibrated by chi-squared(q), as large-sample theory suggests, they are different tests with
different sizes, and the exact sizes follow from the F distribution.
"""
import numpy as np
import statsmodels.api as sm
from scipy import stats

from regbook import Generated

data = sm.datasets.statecrime.load_pandas().data
data.index = data.index.str.strip()
data = data.drop(index="District of Columbia")
codes = "SWWSWWNSSSWWMMMMSSNSNMMSMWMWNNWNSMMSWNNSMSSWNSWSMW"
region = np.array(list(codes))
y = data["murder"].to_numpy()
n = len(y)
covariates = np.column_stack([data["poverty"], data["single"], data["urban"]])
X0 = np.column_stack([np.ones(n), covariates])
D = np.column_stack([(region == g).astype(float) for g in "NSW"])
X = np.column_stack([X0, D])
r, q = 7, 3


def fit(Z):
    b, *_ = np.linalg.lstsq(Z, y, rcond=None)
    return b, np.sum((y - Z @ b) ** 2)


def loglik_max(sse):
    """Maximized normal log-likelihood of a linear model with residual sum of squares sse."""
    return -n / 2 * (np.log(2 * np.pi * sse / n) + 1)


b, sse = fit(X)
b0, sse0 = fit(X0)
F = (sse0 - sse) / q / (sse / (n - r))

# <<trinity>>
u = (sse0 - sse) / sse                          # the ratio behind every statistic
wald = n * u                                    # n (SSE0 - SSE) / SSE
lr = n * np.log1p(u)                            # n log(SSE0 / SSE)
score = n * u / (1 + u)                         # n (SSE0 - SSE) / SSE0
crit = stats.chi2.ppf(0.95, q)
for name, stat in (("Wald", wald), ("LR", lr), ("score", score)):
    print(f"{name:5s} {stat:.3f}  chi2 p = {stats.chi2.sf(stat, q):.4f}")
print(f"F     {F:.3f}  exact p = {stats.f.sf(F, q, n - r):.4f}")
# <</trinity>>

assert wald > lr > score
assert np.isclose(wald, n * q * F / (n - r))
assert np.isclose(lr, 2 * (loglik_max(sse) - loglik_max(sse0)))

# Wald from its definition: (Lambda' b)' [Lambda' (X'X)^{-1} Lambda sigma2_ML]^{-1} (Lambda' b)
Lam = np.zeros((X.shape[1], q))
Lam[4:, :] = np.eye(q)                          # the three region coefficients
V = np.linalg.inv(X.T @ X)
wald_def = b @ Lam @ np.linalg.solve(Lam.T @ V @ Lam * (sse / n), Lam.T @ b)
assert np.isclose(wald_def, wald)
# score from its definition, at the restricted MLE (b0 padded with zeros, sigma2 = SSE0/n)
bR = np.concatenate([b0, np.zeros(q)])
s2R = sse0 / n
grad = X.T @ (y - X @ bR) / s2R                 # d loglik / d beta at the restricted MLE
info = X.T @ X / s2R                            # information for beta (the sigma2 block separates)
score_def = grad @ np.linalg.solve(info, grad)
assert np.isclose(score_def, score)

# exact sizes of the chi-squared-calibrated tests, for q = 3 and r = 7 as here
def exact_size(nn, kind, rr=r, qq=q, alpha=0.05):
    c = stats.chi2.ppf(1 - alpha, qq)
    if kind == "wald":
        ucrit = c / nn
    elif kind == "lr":
        ucrit = np.expm1(c / nn)
    else:
        ucrit = c / (nn - c) if nn > c else np.inf
    return stats.f.sf(ucrit * (nn - rr) / qq, qq, nn - rr)

sizes = {nn: {k: exact_size(nn, k) for k in ("wald", "lr", "score")} for nn in (15, 30, 50, 100, 400)}
for nn, row in sizes.items():
    assert row["wald"] > row["lr"] > row["score"]
    assert row["wald"] > 0.05
assert sizes[400]["wald"] < 0.06 and sizes[400]["score"] > 0.045
assert sizes[15]["wald"] > 0.2
assert all(sizes[a]["wald"] > sizes[b_]["wald"] for a, b_ in ((15, 30), (30, 50), (50, 100), (100, 400)))

gen = Generated("ch11", "likelihood_ratio")
gen.num("wald", wald, 3)
gen.num("lr", lr, 3)
gen.num("score", score, 3)
gen.num("p_wald", stats.chi2.sf(wald, q), 4)
gen.num("p_lr", stats.chi2.sf(lr, q), 4)
gen.num("p_score", stats.chi2.sf(score, q), 4)
gen.num("p_F", stats.f.sf(F, q, n - r), 4)
gen.num("chi2crit", crit, 3)
for nn, row in sizes.items():
    for k, v in row.items():
        gen.num(f"size_{k}_{nn}", v, 3)
gen.write()
