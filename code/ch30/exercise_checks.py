"""Chapter 30: numbers quoted in exercise solutions, and checks of the claims they make."""
import numpy as np
from scipy import stats

from regbook import Generated

gen = Generated("ch30", "exercise_checks")

# Section 30.1, exercise A1: the rules at z = 2.5 with lambda = 1
z, lam, gam, a = 2.5, 1.0, 3.0, 3.7
assert np.isclose((z - lam) / (1 - 1 / gam), 2.25)
scad = ((a - 1) * z - a * lam) / (a - 2)
assert 2 * lam < z <= a * lam
gen.num("scad_25", scad, 2)

# Section 30.3, exercise A1: entry probabilities of null groups
p1 = stats.chi2.sf(4, 1)
p10 = stats.chi2.sf(4, 10)
p10w = stats.chi2.sf(40, 10)
gen.num("chi1_4", p1, 3)
gen.num("chi10_4", p10, 3)
gen.num("chi10_40", p10w, 1, sci=True)

# Section 30.4, exercise A1: orthonormal design, nu = 0.1; coordinate 1 is chosen until 3 (0.9)^m < 1
z = np.array([3.0, -1.0, 0.5])
b = np.zeros(3)
chosen = []
for m in range(12):
    j = int(np.argmax(np.abs(z - b)))
    chosen.append(j)
    b[j] += 0.1 * (z - b)[j]
assert chosen[:11] == [0] * 11 and chosen[11] == 1
assert np.allclose(b[0], 3 * (1 - 0.9 ** 11))
gen.write()
