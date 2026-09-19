"""Chapter 13, Section 6: why the Newman-Keuls procedure is not developed.

Balanced one-way layout, g = 10 groups of 5 observations. The true means come in five
well separated pairs, equal within each pair. Newman-Keuls ends by testing each pair
with the two-mean critical value, so the chance of a false difference is close to that
of five unadjusted tests. Tukey's procedure keeps the rate below 0.05.
"""
import numpy as np
from scipy import stats

from regbook import Generated

alpha, g, m_rep, reps = 0.05, 10, 5, 10_000
nu = g * (m_rep - 1)
mu = np.repeat(np.arange(5) * 20.0, 2)                # pairs 0,0, 20,20, 40,40, ...

# <<newman-keuls>>
def newman_keuls(means, s, m_rep, q_crit):
    """Pairs (k, l) declared different by the Newman-Keuls step-down procedure;
    q_crit[w] is the upper-alpha point of the studentized range of w means."""
    g = len(means)
    order = np.argsort(means)
    ybar = means[order]
    se = s / np.sqrt(m_rep)
    kept = []                                         # windows found homogeneous
    different = set()
    for width in range(g, 1, -1):                     # ranges of width, width-1, ..., 2
        crit = q_crit[width] * se
        for i in range(g - width + 1):
            j = i + width - 1
            if any(a <= i and j <= b for a, b in kept):
                continue                              # inside a homogeneous window
            if ybar[j] - ybar[i] > crit:
                different.add((min(order[i], order[j]), max(order[i], order[j])))
            else:
                kept.append((i, j))
    return different

q_crit = {w: stats.studentized_range.ppf(1 - alpha, w, nu) for w in range(2, g + 1)}
rng = np.random.default_rng(20260919)
q_tukey = stats.studentized_range.ppf(1 - alpha, g, nu)
nk_false = tukey_false = 0
for r in range(reps):
    y = mu[:, None] + rng.normal(size=(g, m_rep))
    means = y.mean(axis=1)
    s = np.sqrt(np.sum((y - means[:, None]) ** 2) / nu)
    true_pairs = [(2 * i, 2 * i + 1) for i in range(5)]
    nk = newman_keuls(means, s, m_rep, q_crit)
    nk_false += any(pair in nk for pair in true_pairs)
    tukey_false += any(abs(means[a] - means[b]) > q_tukey * s / np.sqrt(m_rep)
                       for a, b in true_pairs)
print(f"P(some false difference): Newman-Keuls {nk_false / reps:.3f}, "
      f"Tukey {tukey_false / reps:.3f}")
# <</newman-keuls>>

fwer_nk, fwer_tukey = nk_false / reps, tukey_false / reps
se = np.sqrt(0.25 / reps)
assert fwer_tukey <= alpha + 3 * np.sqrt(alpha * (1 - alpha) / reps)
assert fwer_nk > 0.15                                 # far above the nominal 0.05
# with five independent level-alpha tests one would get 1 - 0.95^5
five = 1 - (1 - alpha) ** 5
assert abs(fwer_nk - five) < 0.03                     # dependence through s only
# Duncan's first step has level 1 - (1 - alpha)^(g - 1) under the complete null
duncan = 1 - (1 - alpha) ** (g - 1)

gen = Generated("ch13", "newman_keuls", prefix="nk")
gen.int("reps", reps)
gen.int("nu", nu)
gen.num("fwernk", fwer_nk, 3)
gen.num("fwertukey", fwer_tukey, 3)
gen.num("five", five, 3)
gen.num("duncan", duncan, 3)
gen.write()
