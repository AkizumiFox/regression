"""Chapter 41, Section 1: a real record with genuine gaps.

Weekly flask measurements of atmospheric CO2 at Mauna Loa, March 1958 to December
2001 (Keeling and Whorf; public domain, statsmodels.datasets.co2). Of the 2284
weeks in the record, 59 carry no value. The gaps are not spread evenly over the
record: almost all of them fall in the first twelve years. That is missingness
that depends on an observed covariate, and the section uses it to separate the
mechanism from its consequences for an analysis.
"""
import matplotlib.pyplot as plt
import numpy as np
import statsmodels.api as sm

from regbook import COLORS, Generated, figure_path, use_book_style

# <<pattern>>
co2 = sm.datasets.co2.load_pandas().data["co2"]       # one value per week, some absent
t = co2.index.year + (co2.index.dayofyear - 0.5) / 365.25
r = (~co2.isna()).to_numpy().astype(float)            # 1 = recorded, 0 = missing
n, n_mis = len(r), int((1 - r).sum())
print(f"{n} weeks, {n_mis} missing ({100 * n_mis / n:.1f}%)")

era = (t >= 1970).astype(float)                       # one observed covariate: is it 1970 or later?
rate_early = 1 - r[era == 0].mean()
rate_late = 1 - r[era == 1].mean()
print(f"missing before 1970: {100 * rate_early:.1f}%   from 1970 on: {100 * rate_late:.1f}%")
# <</pattern>>

# longest run of consecutive gaps
runs, run = [], 0
for v in r:
    if v == 0:
        run += 1
    elif run:
        runs.append(run)
        run = 0
if run:
    runs.append(run)
assert n == 2284 and n_mis == 59
assert max(runs) == 18 and len(runs) == 22

# <<fit>>
# the model of Section 5.3: quadratic trend plus an annual cycle and its first harmonic
tc = t - 1980.0
X = np.column_stack([np.ones(n), tc, tc**2,
                     np.cos(2 * np.pi * t), np.sin(2 * np.pi * t),
                     np.cos(4 * np.pi * t), np.sin(4 * np.pi * t)])
obs = r == 1
y_obs, X_obs = co2.to_numpy()[obs], X[obs]
beta_cc, *_ = np.linalg.lstsq(X_obs, y_obs, rcond=None)   # complete-case fit
sse = float(np.sum((y_obs - X_obs @ beta_cc) ** 2))
p = X.shape[1]
print(f"slope at 1980: {beta_cc[1]:.4f} ppm per year;  s = {np.sqrt(sse / (obs.sum() - p)):.4f} ppm")
# <</fit>>

# a logistic regression of the missingness indicator on time: the mechanism depends on t
logit = sm.Logit(1 - r, np.column_stack([np.ones(n), tc])).fit(disp=0)
odds_decade = float(np.exp(10 * logit.params[1]))
assert odds_decade < 0.5                                   # gaps become much rarer over time
seasonal_amp = float(np.hypot(beta_cc[3], beta_cc[4]))

gen = Generated("ch41", "mechanisms")
gen.int("n", n)
gen.int("n_mis", n_mis)
gen.int("n_obs", int(obs.sum()))
gen.int("runs", len(runs))
gen.int("longest", max(runs))
gen.num("pct_mis", 100 * n_mis / n, 1)
gen.num("pct_early", 100 * rate_early, 1)
gen.num("pct_late", 100 * rate_late, 1)
gen.num("odds_decade", odds_decade, 3)
gen.num("slope", beta_cc[1], 4)
gen.num("s", np.sqrt(sse / (obs.sum() - p)), 4)
gen.num("amp", seasonal_amp, 3)
gen.write()

# ---- the record and its gaps ------------------------------------------------
use_book_style()
fig, axes = plt.subplots(1, 2, figsize=(5.6, 2.3))
ax = axes[0]
ax.plot(t[obs], co2.to_numpy()[obs], color=COLORS["accent"], linewidth=0.5)
for tm in t[~obs]:
    ax.axvline(tm, color=COLORS["second"], linewidth=0.4, alpha=0.7, zorder=0)
ax.set_xlabel("year")
ax.set_ylabel("CO$_2$ (ppm)")
ax.set_title("(a) weekly record, gaps marked")

ax = axes[1]
years = np.arange(1958, 2002)
frac = np.array([1 - r[(t >= yr) & (t < yr + 1)].mean() for yr in years])
ax.bar(years, 100 * frac, color=COLORS["third"], width=0.8, linewidth=0)
ax.set_xlabel("year")
ax.set_ylabel("weeks missing (%)")
ax.set_title("(b) missing fraction by year")
fig.tight_layout()
fig.savefig(figure_path("ch41", "co2_gaps"))
