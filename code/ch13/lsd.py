"""Chapter 13, Section 1: the protected least significant difference.

Gatekeeping every pairwise t test behind a level-alpha F test controls the familywise
error rate when all means are equal (weak control), but not when one mean is far
from the others (no strong control). With three groups it does control strongly.
"""
import itertools

import numpy as np
from scipy import stats

from regbook import Generated

alpha, reps, m_rep = 0.05, 20_000, 5

# <<protected>>
def protected_lsd_fwer(mu, m_rep, reps, alpha, seed):
    """Share of data sets in which the protected LSD rejects some true equality."""
    rng = np.random.default_rng(seed)
    g = len(mu)
    nu = g * (m_rep - 1)
    y = np.asarray(mu)[None, :, None] + rng.normal(size=(reps, g, m_rep))
    means = y.mean(axis=2)
    s2 = ((y - means[:, :, None]) ** 2).sum(axis=(1, 2)) / nu
    between = m_rep * ((means - means.mean(axis=1, keepdims=True)) ** 2).sum(axis=1)
    gate = between / (g - 1) / s2 > stats.f.ppf(1 - alpha, g - 1, nu)   # the F test
    crit = stats.t.ppf(1 - alpha / 2, nu) * np.sqrt(2 * s2 / m_rep)
    false = np.zeros(reps, dtype=bool)
    for k, l in itertools.combinations(range(g), 2):
        if mu[k] == mu[l]:                            # a true null hypothesis
            false |= gate & (np.abs(means[:, k] - means[:, l]) > crit)
    return false.mean()

all_equal = protected_lsd_fwer([0.0] * 10, m_rep, reps, alpha, seed=1)
one_apart = protected_lsd_fwer([0.0] * 9 + [10.0], m_rep, reps, alpha, seed=2)
print(f"10 groups, all means equal:        FWER = {all_equal:.3f}")
print(f"10 groups, one mean far from rest: FWER = {one_apart:.3f}")
# <</protected>>

three = protected_lsd_fwer([0.0, 0.0, 1.5], m_rep, reps, alpha, seed=3)
three_far = protected_lsd_fwer([0.0, 0.0, 10.0], m_rep, reps, alpha, seed=4)
se = np.sqrt(alpha * (1 - alpha) / reps)
assert all_equal <= alpha + 4 * se                     # weak control
assert one_apart > 0.3                                 # no strong control
assert three <= alpha + 4 * se and three_far <= alpha + 4 * se   # g = 3: strong control

gen = Generated("ch13", "lsd", prefix="lsd")
gen.num("allequal", all_equal, 3)
gen.num("oneapart", one_apart, 3)
gen.num("threefar", three_far, 3)
gen.write()
