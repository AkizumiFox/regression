"""Chapter 11, Section 2: the null distribution of F beyond normal errors.

The design is that of the region test (Section 11.1). Data are generated under the reduced
model with four error laws: independent normal; spherical multivariate t with 3 degrees of
freedom (a common random scale, so the errors are dependent); independent t with 3 degrees
of freedom; independent centred exponential. Under a spherical law, F is exactly the same
function of the error direction as under normality, so its null distribution is unchanged.
"""
import numpy as np
import statsmodels.api as sm
from scipy import stats

from regbook import Generated

rng = np.random.default_rng(1103)

data = sm.datasets.statecrime.load_pandas().data
data.index = data.index.str.strip()
data = data.drop(index="District of Columbia")
codes = "SWWSWWNSSSWWMMMMSSNSNMMSMWMWNNWNSMMSWNNSMSSWNSWSMW"
region = np.array(list(codes))
n = len(data)
covariates = np.column_stack([data["poverty"], data["single"], data["urban"]])
X0 = np.column_stack([np.ones(n), covariates])
D = np.column_stack([(region == g).astype(float) for g in "NSW"])
X = np.column_stack([X0, D])


def proj(Z):
    Q, _ = np.linalg.qr(Z)
    return Q @ Q.T


M, M0 = proj(X), proj(X0)
q, nu = 3, n - 7
crit = stats.f.ppf(0.95, q, nu)

# <<spherical>>
def F_stat(E):
    """F statistic for each row of E, used as a data vector under the reduced model."""
    num = np.sum((E @ (M - M0)) ** 2, axis=1) / q
    den = np.sum((E - E @ M) ** 2, axis=1) / nu
    return num / den

reps = 20000
Z = rng.standard_normal((reps, n))
scale = 1 / np.sqrt(rng.chisquare(3, size=(reps, 1)) / 3)       # one scale per data set
errors = {
    "normal": Z,
    "spherical t3": Z * scale,                                    # dependent, heavy-tailed
    "independent t3": rng.standard_t(3, size=(reps, n)),
    "exponential": rng.exponential(size=(reps, n)) - 1.0,         # skewed
}
for name, E in errors.items():
    print(f"{name:15s} rejection rate {np.mean(F_stat(E) > crit):.4f}")
# <</spherical>>

rates = {name: np.mean(F_stat(E) > crit) for name, E in errors.items()}
# the spherical case reproduces the normal statistics exactly, data set by data set
assert np.allclose(F_stat(errors["spherical t3"]), F_stat(Z))
se = np.sqrt(0.05 * 0.95 / reps)
assert abs(rates["normal"] - 0.05) < 3.5 * se
assert rates["spherical t3"] == rates["normal"]
assert abs(rates["independent t3"] - 0.05) < 0.01 and abs(rates["exponential"] - 0.05) < 0.01
# the null F distribution itself: Kolmogorov-Smirnov against F(3, 43)
ks = stats.kstest(F_stat(Z), stats.f(q, nu).cdf)
assert ks.pvalue > 0.001
# the spherical errors are uncorrelated but not independent: their sizes move together
Et = errors["spherical t3"]
rho_signed = stats.spearmanr(Et[:, 0], Et[:, 1]).statistic
rho_size = stats.spearmanr(np.abs(Et[:, 0]), np.abs(Et[:, 1])).statistic
assert abs(rho_signed) < 0.03 and rho_size > 0.1

gen = Generated("ch11", "null_distribution")
gen.int("reps", reps)
gen.num("rate_normal", rates["normal"], 4)
gen.num("rate_spherical", rates["spherical t3"], 4)
gen.num("rate_t3", rates["independent t3"], 4)
gen.num("rate_exp", rates["exponential"], 4)
gen.num("mc_se", se, 4)
gen.num("rho_size", rho_size, 2)
gen.write()
