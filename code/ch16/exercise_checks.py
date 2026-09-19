"""Chapter 16: numbers quoted in exercise solutions, and checks of exercise claims.

(1) Section 16.1, the power of the additive-model F test for rows when one level of A is
    shifted by sigma, with a = 4, b = 5 and m = 1 or 2 observations per cell.
(2) Section 16.2, the 3 x 3 table with consistent orders that no increasing transformation
    makes additive (double cancellation fails), and the 2 x 2 construction.
(3) Section 16.3, the noncentralities for a single deviant cell.
(4) Section 16.5, relabelling the nested factor changes SS_B and SS_AB but not their sum.
"""
import itertools

import numpy as np
from scipy import stats

from regbook import Generated

# ---- (1) power of the row test in the additive model ---------------------------------
a, b = 4, 5
power = {}
for m in (1, 2):
    gamma = b * m * (a - 1) / a                       # b m delta^2 (a-1)/a / sigma^2 with delta = sigma
    df2 = a * b * m - a - b + 1
    crit = stats.f.ppf(0.95, a - 1, df2)
    power[m] = stats.ncf.sf(crit, a - 1, df2, gamma)
    print(f"m = {m}: gamma = {gamma:.2f}, error df = {df2}, power = {power[m]:.2f}")
assert np.isclose(b * 1 * (a - 1) / a, 3.75)

# ---- (2) the 3 x 3 table -------------------------------------------------------------
mu = np.array([[1, 3, 6], [4, 7, 9], [5, 10, 11]])
assert np.all(np.diff(mu, axis=0) > 0) and np.all(np.diff(mu, axis=1) > 0)
assert mu[1, 0] > mu[0, 1] and mu[2, 1] > mu[1, 2] and mu[2, 0] < mu[0, 2]
# no additive table u_i + v_j reproduces the order of the nine values: a linear feasibility check
from scipy.optimize import linprog

cells = [(i, j) for i in range(3) for j in range(3)]
A_ub, b_ub = [], []
for (i, j), (k, l) in itertools.permutations(cells, 2):
    if mu[i, j] > mu[k, l]:                           # need u_i + v_j >= u_k + v_l + 1
        row = np.zeros(6)
        row[i] -= 1; row[3 + j] -= 1; row[k] += 1; row[3 + l] += 1
        A_ub.append(row); b_ub.append(-1.0)
res = linprog(np.zeros(6), A_ub=np.array(A_ub), b_ub=np.array(b_ub), bounds=[(None, None)] * 6)
assert res.status == 2                                # infeasible
# a 2 x 2 table with consistent orders is always removable: the construction of the solution
mu2 = np.array([[5.0, 9.0], [1.0, 4.0]])
g_vals = {1.0: 0.0, 5.0: 2.0, 4.0: 1.0, 9.0: 3.0}     # g(mu21)=0, g(mu11)=x=2 > g(mu22)=y=1, g(mu12)=x+y
G = np.vectorize(g_vals.get)(mu2)
assert np.isclose(G[0, 0] - G[1, 0], G[0, 1] - G[1, 1])
assert all((mu2.ravel()[p] < mu2.ravel()[q]) == (G.ravel()[p] < G.ravel()[q]) for p in range(4) for q in range(4))

# ---- (3) a single deviant cell ------------------------------------------------------
a3, b3, m3, delta = 4, 6, 2, 1.0
M = np.zeros((a3, b3)); M[0, 0] = delta
alpha = M.mean(1) - M.mean(); beta = M.mean(0) - M.mean()
gam = M - M.mean(1, keepdims=True) - M.mean(0, keepdims=True) + M.mean()
assert np.isclose(b3 * m3 * (alpha ** 2).sum(), m3 * delta ** 2 * (a3 - 1) / (a3 * b3))
assert np.isclose(a3 * m3 * (beta ** 2).sum(), m3 * delta ** 2 * (b3 - 1) / (a3 * b3))
assert np.isclose(m3 * (gam ** 2).sum(), m3 * delta ** 2 * (a3 - 1) * (b3 - 1) / (a3 * b3))

# ---- (4) relabelling a nested factor ------------------------------------------------
inner = np.array([[1.0, 3.0], [1.0, 3.0]])
def crossed_ss(T):
    col = T.mean(0) - T.mean()
    inter = T - T.mean(1, keepdims=True) - T.mean(0, keepdims=True) + T.mean()
    return 2 * (col ** 2).sum(), (inter ** 2).sum()
ssB, ssAB = crossed_ss(inner)
ssB2, ssAB2 = crossed_ss(np.array([[1.0, 3.0], [3.0, 1.0]]))   # labels swapped in the second row
assert np.isclose(ssB, 4) and np.isclose(ssB2, 0) and np.isclose(ssB + ssAB, ssB2 + ssAB2)

gen = Generated("ch16", "exercise_checks")
gen.num("power_m1", power[1], 2)
gen.num("power_m2", power[2], 2)
gen.write()
