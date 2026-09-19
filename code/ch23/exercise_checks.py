"""Chapter 23: numbers quoted in exercise solutions.

(a) Section 23.4, Exercise B1: exact permutation p-value for a slope with n = 4.
(b) Section 23.6: Pr(max* = max) for the bootstrap of a maximum at n = 10 and n = 1000.
(c) Section 23.6, Exercise B1: probability that one case is drawn at least twice, n = 15.
"""
import itertools

import numpy as np
from scipy import stats

from regbook import Generated

gen = Generated("ch23", "exercise_checks", prefix="ex")

# (a)
x = np.array([1.0, 2.0, 3.0, 4.0])
y = np.array([1.2, 0.8, 2.9, 3.1])
xc = x - x.mean()
s_obs = xc @ y
s_all = np.array([xc @ np.array(p) for p in itertools.permutations(y)])
count = int(np.sum(np.abs(s_all) >= abs(s_obs) - 1e-9))
p_exact = count / len(s_all)
p_t = stats.linregress(x, y).pvalue
assert count == 6 and np.isclose(s_obs, 3.9)
top = sorted(np.round(s_all[s_all > 3.5], 6))
assert np.allclose(top, [3.7, 3.9, 4.1, 4.3])
gen.num("perm:s", s_obs, 1)
gen.int("perm:count", count)
gen.num("perm:p", p_exact, 2)
gen.num("perm:pt", p_t, 3)

# (b)
for nn in (10, 1000):
    gen.num(f"atom:{nn}", 1 - (1 - 1 / nn) ** nn, 3)
# (c) Section 23.6, Exercise B1: the far case drawn at least twice, n = 15
n15 = 15
p0, p1 = (1 - 1 / n15) ** n15, (1 - 1 / n15) ** (n15 - 1)
twice = 1 - p0 - p1
assert np.isclose(twice, stats.binom.sf(1, n15, 1 / n15))
gen.num("lev:p0", p0, 4)
gen.num("lev:p1", p1, 4)
gen.num("lev:twice", twice, 4)
gen.write()
