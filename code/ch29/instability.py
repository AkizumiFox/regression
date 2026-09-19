"""Chapter 29, Section 4: how stable is the selected model?

Case-resampling bootstrap of the state data (log violent crime, five candidate regressors,
50 states, public domain): repeat the all-subsets AIC and BIC choices on each resample and
count what comes out.
"""
import itertools
from collections import Counter

import matplotlib.pyplot as plt
import numpy as np

from regbook import COLORS, Generated, figure_path, use_book_style

gen = Generated("ch29", "instability", prefix="ins")

# <<bootstrap>>
import itertools
from collections import Counter
import numpy as np
import statsmodels.api as sm
data = sm.datasets.statecrime.load_pandas().data.drop(index="District of Columbia")
names = ["hs_grad", "poverty", "single", "white", "urban"]
y_all = np.log(data["violent"].to_numpy())
Z_all = data[names].to_numpy()
n = len(y_all)
subsets = [S for k in range(6) for S in itertools.combinations(range(5), k)]


def select(Z, y, penalty):
    def crit(S):
        X = np.column_stack([np.ones(len(y)), Z[:, list(S)]])
        b, *_ = np.linalg.lstsq(X, y, rcond=None)
        return n * np.log(np.sum((y - X @ b) ** 2) / n) + penalty * (len(S) + 2)
    return min(subsets, key=crit)


rng = np.random.default_rng(2911)
B = 200                                                       # the book's figure uses 1000
chosen = {"AIC": Counter(), "BIC": Counter()}
for _ in range(B):
    idx = rng.integers(0, n, size=n)                          # resample whole states
    chosen["AIC"][select(Z_all[idx], y_all[idx], 2.0)] += 1
    chosen["BIC"][select(Z_all[idx], y_all[idx], np.log(n))] += 1
for crit, counts in chosen.items():
    print(crit, "distinct models:", len(counts))
    for S, c in counts.most_common(3):
        print("   ", [names[j] for j in S], c / B)
# <</bootstrap>>

B = 1000
chosen = {"AIC": Counter(), "BIC": Counter()}
for _ in range(B):
    idx = rng.integers(0, n, size=n)
    chosen["AIC"][select(Z_all[idx], y_all[idx], 2.0)] += 1
    chosen["BIC"][select(Z_all[idx], y_all[idx], np.log(n))] += 1
orig = {"AIC": (2, 4), "BIC": (2,)}
assert select(Z_all, y_all, 2.0) == orig["AIC"] and select(Z_all, y_all, np.log(n)) == orig["BIC"]
incl = {}
for crit, counts in chosen.items():
    incl[crit] = np.array([sum(c for S, c in counts.items() if j in S) for j in range(5)]) / B
    gen.int(f"{crit}_distinct", len(counts))
    gen.num(f"{crit}_orig", counts[orig[crit]] / B, 2)
    top = counts.most_common(1)[0]
    gen.num(f"{crit}_top", top[1] / B, 2)
    gen.text(f"{crit}_topmodel", " + ".join(names[j] for j in top[0]))
    for j, nm in enumerate(names):
        gen.num(f"{crit}_{nm}", incl[crit][j], 2)
    print(crit, len(counts), counts[orig[crit]] / B, incl[crit])
assert np.argmax(incl["AIC"]) == 2 and np.argmax(incl["BIC"]) == 2   # single parenthood most often
assert chosen["AIC"][orig["AIC"]] / B < 0.3                 # the AIC model is far from stable
assert len(chosen["AIC"]) > 20 and len(chosen["BIC"]) > 20
gen.int("B", B)
gen.num("r_single_white", np.corrcoef(Z_all[:, 2], Z_all[:, 3])[0, 1], 2)
gen.write()

use_book_style()
fig, ax = plt.subplots(figsize=(4.4, 2.6))
xs = np.arange(5)
ax.bar(xs - 0.18, incl["AIC"], width=0.36, color=COLORS["accent"], label="AIC")
ax.bar(xs + 0.18, incl["BIC"], width=0.36, color=COLORS["second"], label="BIC")
ax.set_xticks(xs, ["high school", "poverty", "single parent", "white", "urban"])
ax.set_ylabel("fraction of resamples\nincluding the regressor")
ax.set_ylim(0, 1.05)
ax.legend(frameon=False, loc="upper left")
fig.tight_layout()
fig.savefig(figure_path("ch29", "instability"))
