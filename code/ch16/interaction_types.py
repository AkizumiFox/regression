"""Chapter 16, Section 2: three kinds of two-way cell-mean patterns, as interaction plots.

(a) additive: parallel profiles; (b) removable: the profiles fan out, but the cell means are
exp(additive), so the logarithm makes them parallel; (c) qualitative: the profiles cross,
and no increasing transformation can remove the interaction. The script checks each claim,
including the order criterion of the proposition on removable interaction.
"""
import itertools

import matplotlib.pyplot as plt
import numpy as np

from regbook import COLORS, figure_path, use_book_style

# <<patterns>>
additive = np.array([[4.0, 6.0, 7.0], [2.0, 4.0, 5.0]])            # rows: levels of A
removable = np.exp(np.array([[1.0, 1.6, 2.2], [0.4, 1.0, 1.6]]))  # exp of an additive table
crossing = np.array([[3.0, 5.0, 7.0], [6.0, 5.0, 4.0]])


def interaction(mu):
    """mu_ij - mu_i. - mu_.j + mu_.. (equal weights)."""
    return mu - mu.mean(1, keepdims=True) - mu.mean(0, keepdims=True) + mu.mean()


def same_orders(mu):
    """Necessary condition for removability: every column orders the rows alike, and every row the columns."""
    a, b = mu.shape
    rows_ok = all(len(set(np.sign(mu[i] - mu[k]))) == 1 for i, k in itertools.combinations(range(a), 2))
    cols_ok = all(len(set(np.sign(mu[:, j] - mu[:, l]))) == 1 for j, l in itertools.combinations(range(b), 2))
    return rows_ok and cols_ok


for name, mu in [("additive", additive), ("removable", removable), ("crossing", crossing)]:
    print(f"{name:9s} max |interaction| raw {np.abs(interaction(mu)).max():.3f}"
          f"  log {np.abs(interaction(np.log(mu))).max():.3f}   same orders: {same_orders(mu)}")
# <</patterns>>

assert np.abs(interaction(additive)).max() < 1e-12
assert np.abs(interaction(removable)).max() > 0.1 and np.abs(interaction(np.log(removable))).max() < 1e-12
assert same_orders(additive) and same_orders(removable) and not same_orders(crossing)
# the log of the additive table is not additive: interaction depends on the scale in both directions
assert np.abs(interaction(np.log(additive))).max() > 0.05
# weighted row differences in the crossing table range over [min delta, max delta] and change sign
delta = crossing[0] - crossing[1]
assert delta.min() < 0 < delta.max() and np.isclose(delta.mean(), 0.0)

use_book_style()
fig, axes = plt.subplots(1, 3, figsize=(5.8, 2.1), sharey=False)
x = np.arange(1, 4)
for ax, mu, title in [(axes[0], additive, "(a) additive"),
                      (axes[1], removable, "(b) removable"),
                      (axes[2], crossing, "(c) crossing")]:
    for i, colour, mark in [(0, COLORS["accent"], "o"), (1, COLORS["second"], "s")]:
        ax.plot(x, mu[i], color=colour, marker=mark, markersize=3.5, label=f"$A_{i + 1}$")
    ax.set_xticks(x, [r"$B_1$", r"$B_2$", r"$B_3$"])
    ax.set_title(title)
    ax.set_ylabel(r"cell mean $\mu_{ij}$" if ax is axes[0] else "")
axes[0].legend(frameon=False, loc="upper left")
fig.tight_layout()
fig.savefig(figure_path("ch16", "interaction_types"))
