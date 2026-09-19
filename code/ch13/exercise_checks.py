"""Chapter 13: the small hand computations in the exercises and solutions."""
import numpy as np
from statsmodels.stats.multitest import multipletests

from regbook import Generated

alpha = 0.05

# <<five>>
p = np.array([0.001, 0.011, 0.012, 0.04, 0.2])        # the five p-values of the exercises
for method in ("bonferroni", "holm", "fdr_bh"):
    reject, adjusted, *_ = multipletests(p, alpha, method)
    print(f"{method:10s} rejects {reject.sum()}   adjusted p-values {np.round(adjusted, 4)}")
# <</five>>

assert multipletests(p, alpha, "bonferroni")[0].sum() == 1
assert multipletests(p, alpha, "holm")[0].sum() == 3
assert multipletests(p, alpha, "fdr_bh")[0].sum() == 4
holm_adj = multipletests(p, alpha, "holm")[1]
assert np.allclose(holm_adj, [0.005, 0.044, 0.044, 0.08, 0.2])
bh_adj = multipletests(p, alpha, "fdr_bh")[1]
assert np.allclose(bh_adj, [0.005, 0.02, 0.02, 0.05, 0.2])

# familywise error of m independent level-0.01 tests
fw10, fw20 = 1 - 0.99 ** 10, 1 - 0.99 ** 20
sidak20 = 1 - 0.95 ** (1 / 20)
assert 0.05 / 20 < sidak20
# range test against F test for g = 3 in the limit nu -> infinity (units of s / sqrt(m))
from scipy import stats
q3 = stats.studentized_range.ppf(1 - alpha, 3, np.inf)
chi2_2 = stats.chi2.ppf(1 - alpha, 2)                # (g - 1) F_alpha(g - 1, inf)
a_range = q3 / 2                                      # range 2a exceeds q3
a_f_mid = np.sqrt(chi2_2 / 2)                         # means (-a, 0, a): SS = 2 a^2
a_f_end = np.sqrt(3 * chi2_2 / 8)                     # means (-a, -a, a): SS = 8 a^2 / 3
assert a_f_end < a_range < a_f_mid                    # each test rejects something the other does not

gen = Generated("ch13", "exercise_checks", prefix="ex")
gen.num("q3inf", q3, 3)
gen.num("chi2two", chi2_2, 3)
gen.num("arange", a_range, 3)
gen.num("afmid", a_f_mid, 3)
gen.num("afend", a_f_end, 3)
gen.num("fw10", fw10, 4)
gen.num("fw20", fw20, 4)
gen.num("sidak20", sidak20, 5)
gen.num("bonf20", 0.05 / 20, 5)
gen.write()
