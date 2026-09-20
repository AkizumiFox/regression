"""Chapter 35, Section 2: logit, probit and complementary log-log links on the same binary data.

Same public-domain survey as anes_logistic.py (statsmodels.datasets.anes96).
"""
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy import stats

from regbook import COLORS, Generated, figure_path, use_book_style

# two constants used in the text
sd_logistic = np.pi / np.sqrt(3)                 # standard deviation of the standard logistic law
scale_probit = stats.norm.pdf(0) / 0.25          # beta_logit / beta_probit if the slopes match at 1/2
scale_cloglog = 0.25 / (0.5 * np.log(2))         # the same matching for the complementary log-log
assert np.isclose(sd_logistic, 1.8137993642342178)
assert np.isclose(scale_probit, 1.5957691216057308)

# <<fits>>
import numpy as np
import pandas as pd
import statsmodels.api as sm

anes = sm.datasets.anes96.load_pandas().data
X = pd.DataFrame({"const": 1.0, "PID": anes["PID"], "age": anes["age"] / 10,
                  "educ": anes["educ"], "income": anes["income"]})
y = anes["vote"]

links = {"logit": sm.families.links.Logit(), "probit": sm.families.links.Probit(),
         "cloglog": sm.families.links.CLogLog()}
fits = {k: sm.GLM(y, X, family=sm.families.Binomial(link=v)).fit() for k, v in links.items()}

table = pd.DataFrame({k: f.params for k, f in fits.items()})
table["probit x 1.814"] = table["probit"] * np.pi / np.sqrt(3)
print(table.round(3).to_string())
print("maximized log-likelihoods:",
      {k: round(f.llf, 2) for k, f in fits.items()})
# <</fits>>

p = {k: f.fittedvalues.to_numpy() for k, f in fits.items()}
max_probit = np.abs(p["probit"] - p["logit"]).max()
max_cloglog = np.abs(p["cloglog"] - p["logit"]).max()
assert max_probit < max_cloglog                # probit tracks the logit far more closely

# where the two symmetric links disagree: only in the tails
tail = p["logit"] < 0.05
assert np.abs(p["probit"][tail] - p["logit"][tail]).max() < 0.02

# the rescaled probit slopes are close to the logit slopes
ratio = (table["probit x 1.814"] / table["logit"])[["PID", "age", "income"]]
assert ratio.between(0.9, 1.15).all()

gen = Generated("ch35", "links")
gen.num("sd_logistic", sd_logistic, 4)
gen.num("scale_probit", scale_probit, 4)
gen.num("scale_cloglog", scale_cloglog, 4)
gen.num("sd_gumbel", np.pi / np.sqrt(6), 4)
gen.num("euler", 0.5772156649015329, 4)
gen.num("b_pid_logit", table.loc["PID", "logit"], 3)
gen.num("b_pid_probit", table.loc["PID", "probit"], 3)
gen.num("b_pid_probit_scaled", table.loc["PID", "probit x 1.814"], 3)
gen.num("b_pid_cloglog", table.loc["PID", "cloglog"], 3)
gen.num("ll_logit", fits["logit"].llf, 2)
gen.num("ll_probit", fits["probit"].llf, 2)
gen.num("ll_cloglog", fits["cloglog"].llf, 2)
gen.num("aic_gap_probit", 2 * (fits["logit"].llf - fits["probit"].llf), 2)
gen.num("aic_gap_cloglog", 2 * (fits["logit"].llf - fits["cloglog"].llf), 2)
gen.num("max_probit", max_probit, 4)
gen.num("max_cloglog", max_cloglog, 4)
gen.write()

# ---- response curves and fitted probabilities ---------------------------------------------
use_book_style()
fig, axes = plt.subplots(1, 2, figsize=(5.6, 2.5))

eta = np.linspace(-6, 6, 400)
ax = axes[0]
ax.plot(eta, 1 / (1 + np.exp(-eta)), color=COLORS["accent"], label="logit")
ax.plot(eta, stats.norm.cdf(eta / scale_probit), color=COLORS["second"], linestyle="--",
        label="probit, rescaled")
# match the cloglog at the median and give it the same slope there
ax.plot(eta, 1 - np.exp(-np.exp(np.log(np.log(2)) + eta * scale_cloglog)), color=COLORS["third"],
        linestyle=":", label="cloglog, rescaled")
ax.axhline(0.5, color=COLORS["grid"], linewidth=0.6, zorder=0)
ax.set_xlabel(r"$\eta$")
ax.set_ylabel(r"$\pi$")
ax.legend(frameon=False, loc="upper left")
ax.set_title("(a) rescaled response curves")

ax = axes[1]
ax.scatter(p["logit"], p["probit"] - p["logit"], s=3, color=COLORS["second"],
           linewidths=0, alpha=0.7, label="probit")
ax.scatter(p["logit"], p["cloglog"] - p["logit"], s=3, color=COLORS["third"],
           linewidths=0, alpha=0.7, label="cloglog")
ax.axhline(0, color=COLORS["grid"], linewidth=0.6, zorder=0)
ax.set_xlabel(r"fitted $\hat\pi_i$ under the logit link")
ax.set_ylabel("difference in fitted probability")
ax.legend(frameon=False, markerscale=3, loc="lower left")
ax.set_title("(b) the same data, three links")
fig.tight_layout()
fig.savefig(figure_path("ch35", "link_comparison"))
