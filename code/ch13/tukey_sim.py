"""Chapter 13, Section 3: the studentized range and the coverage of Tukey intervals.

(1) scipy's studentized range quantile agrees with a direct simulation of
    range(Z_1..Z_g) / sqrt(V / nu);
(2) in a balanced one-way layout the Tukey intervals for all pairs have coverage
    exactly 1 - alpha;
(3) in unbalanced layouts the Tukey-Kramer intervals cover with probability at
    least 1 - alpha (Hayter 1984): here, the design of the education example and a
    strongly unbalanced one.
"""
import itertools

import numpy as np
from scipy import stats

from regbook import Generated

alpha, reps = 0.05, 40_000

# <<quantile>>
g, nu = 5, 20
rng = np.random.default_rng(20260919)
Z = rng.normal(size=(200_000, g))
V = rng.chisquare(nu, size=200_000)
Q = (Z.max(axis=1) - Z.min(axis=1)) / np.sqrt(V / nu)   # studentized range
q_sim = np.quantile(Q, 1 - alpha)
q_exact = stats.studentized_range.ppf(1 - alpha, g, nu)
print(f"upper 5% point of the studentized range, g = {g}, nu = {nu}:"
      f" simulated {q_sim:.3f}, scipy {q_exact:.3f}")
# <</quantile>>
assert abs(q_sim - q_exact) < 0.025


# <<coverage>>
def tukey_kramer_coverage(n_k, reps, alpha, seed):
    """Share of data sets in which every Tukey-Kramer interval covers its difference."""
    rng = np.random.default_rng(seed)
    n_k = np.asarray(n_k)
    g, nu = len(n_k), n_k.sum() - len(n_k)
    q = stats.studentized_range.ppf(1 - alpha, g, nu)
    # group means (true means 0) and an independent pooled variance, as in the model
    means = rng.normal(size=(reps, g)) / np.sqrt(n_k)
    s = np.sqrt(rng.chisquare(nu, size=reps) / nu)
    ok = np.ones(reps, dtype=bool)
    for k, l in itertools.combinations(range(g), 2):
        half = q / np.sqrt(2) * s * np.sqrt(1 / n_k[k] + 1 / n_k[l])
        ok &= np.abs(means[:, k] - means[:, l]) <= half
    return ok.mean()

balanced = tukey_kramer_coverage([5] * 5, reps, alpha, seed=1)
education = tukey_kramer_coverage([13, 52, 248, 187, 90, 227, 127], reps, alpha, seed=2)
lopsided = tukey_kramer_coverage([2, 2, 2, 2, 40, 40, 40, 40], reps, alpha, seed=3)
print(f"balanced (5 x 5):           {balanced:.4f}")
print(f"education design:           {education:.4f}")
print(f"four groups of 2, four of 40: {lopsided:.4f}")
# <</coverage>>

se = np.sqrt(alpha * (1 - alpha) / reps)
assert abs(balanced - (1 - alpha)) < 4 * se           # exact in the balanced case
assert education > 1 - alpha - 2 * se                 # conservative (Hayter)
assert lopsided > 1 - alpha + 4 * se                  # visibly conservative

# the extension to all contrasts in the balanced case: coverage of Tukey's bound for
# every contrast equals that for the pairs, because the pairs are the extreme case
rng = np.random.default_rng(4)
g_b, m_b = 5, 5
nu_b = g_b * (m_b - 1)
q_b = stats.studentized_range.ppf(1 - alpha, g_b, nu_b)
means_b = rng.normal(size=(reps, g_b)) / np.sqrt(m_b)
s_b = np.sqrt(rng.chisquare(nu_b, size=reps) / nu_b)
C = rng.normal(size=(400, g_b))
C -= C.mean(axis=1, keepdims=True)
lhs = np.abs(means_b @ C.T)                           # |sum c_k (ybar_k - mu_k)|
rhs = (q_b * s_b / np.sqrt(m_b))[:, None] * 0.5 * np.abs(C).sum(axis=1)[None, :]
all_contrasts = np.all(lhs <= rhs, axis=1)
pairs_only = (means_b.max(axis=1) - means_b.min(axis=1)) <= q_b * s_b / np.sqrt(m_b)
assert np.all(all_contrasts >= pairs_only)             # pairs covered => all contrasts covered

gen = Generated("ch13", "tukey_sim", prefix="tuk")
gen.int("reps", reps)
gen.num("qsim", q_sim, 3)
gen.num("qexact", q_exact, 3)
gen.num("balanced", balanced, 4)
gen.num("education", education, 4)
gen.num("lopsided", lopsided, 4)
gen.num("mcse", se, 4)
gen.write()
