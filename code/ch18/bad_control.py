"""Chapter 18, Section 1: adjusting for a covariate that the treatment affects.

Synthetic randomized experiment with two arms. The 'covariate' is measured after
treatment and is partly caused by it, so adjusting for it removes part of the effect.
"""
import numpy as np

from regbook import Generated

# <<mediator>>
rng = np.random.default_rng(2025)
n, reps = 200, 20_000
theta, delta, gamma = 1.0, 2.0, 0.75                  # direct effect; T -> w; w -> y
total = theta + gamma * delta                          # effect of assigning treatment
est_unadj, est_adj = np.empty(reps), np.empty(reps)
for r in range(reps):
    T = rng.permutation(np.repeat([0.0, 1.0], n // 2))    # randomized assignment
    w = rng.normal(size=n) + delta * T                    # measured after treatment
    y = theta * T + gamma * w + rng.normal(size=n)
    X = np.column_stack([np.ones(n), T, w])
    est_unadj[r] = y[T == 1].mean() - y[T == 0].mean()
    est_adj[r] = np.linalg.lstsq(X, y, rcond=None)[0][1]
print(f"total effect {total:.2f}: unadjusted mean {est_unadj.mean():.3f},"
      f" adjusted for w mean {est_adj.mean():.3f}")
# <</mediator>>

assert abs(est_unadj.mean() - total) < 0.01
assert abs(est_adj.mean() - theta) < 0.01

gen = Generated("ch18", "bad_control")
gen.num("total", total, 2)
gen.num("unadj", est_unadj.mean(), 3)
gen.num("adj", est_adj.mean(), 3)
gen.write()
