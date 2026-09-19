"""Chapter 29, Section 2: all 32 subsets for log violent crime in the 50 states (2009).

Response: log of the violent crime rate per 100,000. Candidate regressors: high-school
graduation, poverty, single-parent households, percentage white, percentage urban; the
intercept is always in. Public-domain data shipped with statsmodels (statecrime); the
District of Columbia is left out, as in Chapter 6.
"""
import itertools

import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import f as fdist

from regbook import COLORS, Generated, figure_path, use_book_style

gen = Generated("ch29", "state_subsets", prefix="ss")

# <<subsets>>
import itertools
import numpy as np
import statsmodels.api as sm
data = sm.datasets.statecrime.load_pandas().data.drop(index="District of Columbia")
names = ["hs_grad", "poverty", "single", "white", "urban"]
y = np.log(data["violent"].to_numpy())
n = len(y)


def sse(cols):
    X = np.column_stack([np.ones(n)] + [data[c].to_numpy() for c in cols])
    b, *_ = np.linalg.lstsq(X, y, rcond=None)
    return np.sum((y - X @ b) ** 2), X.shape[1]


sse_full, P = sse(names)
s2 = sse_full / (n - P)                                  # sigma^2 estimate from the largest model
sst = np.sum((y - y.mean()) ** 2)
rows = []
for k in range(len(names) + 1):
    for cols in itertools.combinations(names, k):
        e, p = sse(cols)
        rows.append(dict(model=cols, p=p, sse=e,
                         cp=e / s2 - n + 2 * p,
                         aic=n * np.log(e / n) + 2 * (p + 1),
                         bic=n * np.log(e / n) + np.log(n) * (p + 1),
                         adjr2=1 - (e / (n - p)) / (sst / (n - 1))))
for crit in ("cp", "aic", "bic"):
    best = min(rows, key=lambda r: r[crit])
    print(f"{crit:4s} chooses {best['model']}")
print("adj R2 chooses", max(rows, key=lambda r: r["adjr2"])["model"])
# <</subsets>>

assert len(rows) == 32
by = {c: min(rows, key=lambda r: r[c])["model"] for c in ("cp", "aic", "bic")}
by["adjr2"] = max(rows, key=lambda r: r["adjr2"])["model"]
assert by["cp"] == ("single", "urban") and by["aic"] == ("single", "urban")
assert by["bic"] == ("single",)
full = [r for r in rows if r["p"] == P][0]
assert np.isclose(full["cp"], P)                           # C_p of the largest model is P

# identity C_p = (P - p)(F - 1) + p with F the test of the submodel against the full model
for r in rows:
    if r["p"] < P:
        F = ((r["sse"] - sse_full) / (P - r["p"])) / s2
        assert np.isclose(r["cp"], (P - r["p"]) * (F - 1) + r["p"])

gen.int("n", n)
gen.int("P", P)
gen.num("s2", s2, 4)
aic_min = min(r["aic"] for r in rows)
bic_min = min(r["bic"] for r in rows)
order = sorted(rows, key=lambda r: r["cp"])
for r in order[:8]:
    tag = "".join(c[0] for c in r["model"])       # e.g. "su" for single + urban
    gen.num(f"{tag}_cp", r["cp"], 2)
    gen.num(f"{tag}_daic", r["aic"] - aic_min, 2)
    gen.num(f"{tag}_dbic", r["bic"] - bic_min, 2)
    gen.num(f"{tag}_adj", r["adjr2"], 3)
    gen.num(f"{tag}_sse", r["sse"], 3)
print([(r["model"], round(r["cp"], 2)) for r in order[:8]])

# the F thresholds for adding one column, at n = 50 and a model with p columns before
p_small = 2
thr_aic = (np.exp(2 / n) - 1) * (n - p_small - 1)
thr_bic = (np.exp(np.log(n) / n) - 1) * (n - p_small - 1)
gen.num("thr_aic", thr_aic, 2)
gen.num("thr_bic", thr_bic, 2)
gen.num("logn", np.log(n), 2)
gen.num("thr_f05", fdist.ppf(0.95, 1, n - p_small - 1), 2)
gen.num("level_aic", fdist.sf(thr_aic, 1, n - p_small - 1), 2)     # the implied test levels
gen.num("level_bic", fdist.sf(thr_bic, 1, n - p_small - 1), 3)
gen.write()

# ---- the C_p plot ------------------------------------------------------------------
use_book_style()
fig, ax = plt.subplots(figsize=(4.4, 3.0))
ps = np.array([r["p"] for r in rows])
cps = np.array([r["cp"] for r in rows])
jit = np.array([(i % 5 - 2) * 0.06 for i in range(len(rows))])
ax.scatter(ps + jit, cps, s=12, color=COLORS["accent"], linewidths=0, zorder=3)
ax.plot([1, P], [1, P], color=COLORS["muted"], linewidth=0.8, linestyle="--")
ax.set_yscale("log")
ax.set_ylim(1, 300)
ax.set_xlabel("p (columns, including the intercept)")
ax.set_ylabel(r"$C_p$ (log scale)")
for r in rows:
    if r["model"] in (("single",), ("single", "urban"), ("poverty", "single", "urban")):
        label = " + ".join(r["model"])
        ax.annotate(label, (r["p"], r["cp"]), xytext=(6, -2), textcoords="offset points", fontsize=7)
fig.tight_layout()
fig.savefig(figure_path("ch29", "cp_plot"))
