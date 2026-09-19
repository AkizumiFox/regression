"""Chapter 15, Section 5: what balance buys in power.

For a fixed total of n = 40 observations in g = 4 groups, the noncentrality of the F test is
guaranteed to be at least (Delta / sigma)^2 * min_{i != j} n_i n_j / (n_i + n_j) whenever the
group means span a range Delta. Balance maximizes this guarantee.
"""
import itertools

import matplotlib.pyplot as plt
import numpy as np
from scipy import stats

from regbook import COLORS, Generated, figure_path, use_book_style

# <<guarantee>>
def guaranteed_gamma(sizes, delta_over_sigma):
    """Smallest noncentrality over all mean configurations with range at least Delta."""
    h = min(a * b / (a + b) for a, b in itertools.combinations(sizes, 2))
    return delta_over_sigma ** 2 * h

def power_f(gamma, q, nu, alpha=0.05):
    return stats.ncf.sf(stats.f.ppf(1 - alpha, q, nu), q, nu, gamma)

allocations = {"balanced": [10, 10, 10, 10], "mild": [7, 9, 11, 13], "strong": [4, 6, 10, 20]}
for name, sizes in allocations.items():
    gam = guaranteed_gamma(sizes, 1.0)
    print(f"{name:9s} {sizes}: guaranteed gamma = {gam:.3f},"
          f" guaranteed power at Delta = sigma: {power_f(gam, 3, 36):.3f}")
# <</guarantee>>

# check the guarantee by brute force: minimize gamma over configurations with range 1
rng = np.random.default_rng(15)
for sizes in allocations.values():
    w = np.array(sizes, float)
    best = np.inf
    for i, j in itertools.permutations(range(4), 2):
        for _ in range(3000):
            mu = rng.uniform(0, 1, 4)
            mu[i], mu[j] = 1.0, 0.0
            best = min(best, np.sum(w * (mu - w @ mu / w.sum()) ** 2))
    assert best >= guaranteed_gamma(sizes, 1.0) - 1e-12
    assert best < guaranteed_gamma(sizes, 1.0) * 1.05          # and it is nearly attained
# balance is optimal among all allocations of 40 into 4 groups of at least 2
alls = [c for c in itertools.combinations_with_replacement(range(2, 35), 4) if sum(c) == 40]
assert max(alls, key=lambda c: guaranteed_gamma(c, 1.0)) == (10, 10, 10, 10)

gen = Generated("ch15", "allocation_power")
for name, sizes in allocations.items():
    gam = guaranteed_gamma(sizes, 1.0)
    gen.num(f"gamma_{name}", gam, 3)
    gen.num(f"power_{name}", power_f(gam, 3, 36), 3)
# Delta / sigma needed for guaranteed power 0.8
delta80 = {}
for name, sizes in allocations.items():
    d = 0.01
    while power_f(guaranteed_gamma(sizes, d), 3, 36) < 0.8:
        d += 0.001
    delta80[name] = d
    gen.num(f"delta80_{name}", d, 2)
assert max(delta80.values()) < 2.6                          # every curve crosses 0.8 on the plot
gen.write()

use_book_style()
fig, ax = plt.subplots(figsize=(4.4, 2.6))
ds = np.linspace(0, 2.6, 200)
for (name, sizes), col, ls in zip(allocations.items(),
                                  [COLORS["accent"], COLORS["third"], COLORS["second"]],
                                  ["-", "--", "-."]):
    ax.plot(ds, [power_f(guaranteed_gamma(sizes, d), 3, 36) if d > 0 else 0.05 for d in ds],
            color=col, ls=ls, label=f"{name}: " + ", ".join(map(str, sizes)))
    ax.plot([delta80[name]], [0.8], "o", ms=3, color=col)   # where power reaches 0.8
ax.axhline(0.8, color=COLORS["muted"], lw=0.6, ls=":")
ax.set_xlabel(r"range of the group means, $\Delta/\sigma$")
ax.set_ylabel("guaranteed power")
ax.legend(frameon=False, loc="lower right", fontsize=7)
fig.tight_layout()
fig.savefig(figure_path("ch15", "allocation_power"))
