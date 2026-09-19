"""Chapter 25, Section 3: chains, forks and colliders in linear structural models.

Simulates each three-variable structure with unit-variance disturbances and compares the
sample correlation and partial correlation of the two end variables with the values
computed from the structural coefficients.
"""
import numpy as np

from regbook import Generated

# <<structures>>
import numpy as np

rng = np.random.default_rng(2503)
n = 100_000
a, b = 1.0, 1.0                                   # the two edge coefficients

def corr_and_partial(x, y, w):
    """Correlation of x and y, and their partial correlation given w (residuals on 1, w)."""
    W = np.column_stack([np.ones(len(w)), w])
    rx = x - W @ np.linalg.lstsq(W, x, rcond=None)[0]
    ry = y - W @ np.linalg.lstsq(W, y, rcond=None)[0]
    return np.corrcoef(x, y)[0, 1], np.corrcoef(rx, ry)[0, 1]

U = rng.normal(size=(3, n))
x = U[0]; m = a * x + U[1]; y = b * m + U[2]      # chain   X -> M -> Y
chain = corr_and_partial(x, y, m)
z = U[0]; x = a * z + U[1]; y = b * z + U[2]      # fork    X <- Z -> Y
fork = corr_and_partial(x, y, z)
x = U[0]; y = U[1]; c = a * x + b * y + U[2]      # collider X -> C <- Y
collider = corr_and_partial(x, y, c)
for name, (r, rp) in [("chain", chain), ("fork", fork), ("collider", collider)]:
    print(f"{name:9s} corr {r:6.3f}   partial corr given middle {rp:6.3f}")
# <</structures>>

# values implied by the coefficients (unit disturbance variances)
chain_r = a * b / np.sqrt(b ** 2 * (a ** 2 + 1) + 1)
fork_r = a * b / np.sqrt((a ** 2 + 1) * (b ** 2 + 1))
coll_p = -a * b / np.sqrt((a ** 2 + 1) * (b ** 2 + 1))
se = 4 / np.sqrt(n)
assert abs(chain[0] - chain_r) < se and abs(chain[1]) < se
assert abs(fork[0] - fork_r) < se and abs(fork[1]) < se
assert abs(collider[0]) < se and abs(collider[1] - coll_p) < se

gen = Generated("ch25", "structures")
gen.num("chainr", chain_r, 3)
gen.num("forkr", fork_r, 3)
gen.num("collp", coll_p, 3)
gen.num("chainrsim", chain[0], 3)
gen.num("forkrsim", fork[0], 3)
gen.num("collpsim", collider[1], 3)
gen.write()
