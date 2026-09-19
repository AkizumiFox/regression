"""Chapter 25, Section 1: a good predictor that is a bad causal model.

Synthetic structural causal model. Illness severity Z (unrecorded) raises the number of
clinic visits X and lowers the health score Y a year later; each visit raises Y by
theta = 2. Observational data give a negative regression slope. A second data set in
which X is set at random by the analyst (an intervention) shows the causal slope.
"""
import matplotlib.pyplot as plt
import numpy as np

from regbook import COLORS, Generated, figure_path, use_book_style

# <<simulate>>
import numpy as np

rng = np.random.default_rng(2501)
n = 2000
alpha, theta, gamma = 1.5, 2.0, -6.0          # Z -> X, X -> Y, Z -> Y

def draw(n, do_x=None):
    """Draw from the model; with do_x, X is set by the analyst instead of by Z."""
    Z = rng.normal(size=n)                     # illness severity, never recorded
    X = 4 + alpha * Z + rng.normal(size=n) if do_x is None else do_x
    Y = 50 + theta * X + gamma * Z + rng.normal(scale=2.0, size=n)
    return X, Y

x_obs, y_obs = draw(n)                                  # observational data
x_do = rng.uniform(x_obs.min(), x_obs.max(), size=n)    # X assigned at random
x_int, y_int = draw(n, do_x=x_do)                       # interventional data

slope_obs, icpt_obs = np.polyfit(x_obs, y_obs, 1)
slope_int, icpt_int = np.polyfit(x_int, y_int, 1)
print(f"slope from observational data:  {slope_obs:.3f}")
print(f"slope from interventional data: {slope_int:.3f}   (theta = {theta})")
# <</simulate>>

# <<mse>>
def mse(x, y, slope, icpt):
    return np.mean((y - icpt - slope * x) ** 2)

causal_icpt = np.mean(y_obs) - theta * np.mean(x_obs)   # causal slope, fitted level
for label, x, y in [("observational", x_obs, y_obs), ("interventional", x_int, y_int)]:
    print(f"{label:15s} MSE: regression line {mse(x, y, slope_obs, icpt_obs):6.2f},"
          f" causal line {mse(x, y, theta, causal_icpt):6.2f}")
# <</mse>>

# population values
var_x = alpha ** 2 + 1
beta_star = theta + gamma * alpha / var_x
var_y = theta ** 2 * var_x + gamma ** 2 + 2 * theta * gamma * alpha + 4.0
r2_pop = beta_star ** 2 * var_x / var_y
r2_obs = np.corrcoef(x_obs, y_obs)[0, 1] ** 2
assert abs(slope_obs - beta_star) < 0.1
assert abs(slope_int - theta) < 0.1
assert abs(r2_obs - r2_pop) < 0.03
# the regression line predicts better where it was fitted, worse under intervention
m_oo, m_oc = mse(x_obs, y_obs, slope_obs, icpt_obs), mse(x_obs, y_obs, theta, causal_icpt)
m_io, m_ic = mse(x_int, y_int, slope_obs, icpt_obs), mse(x_int, y_int, theta, causal_icpt)
assert m_oo < m_oc and m_ic < m_io
# the observational excess is exactly (slope - theta)^2 times the variance of X
assert abs((m_oc - m_oo) - (slope_obs - theta) ** 2 * np.var(x_obs)) < 1e-8
gen = Generated("ch25", "prediction")
gen.int("n", n)
gen.num("betastar", beta_star, 3)
gen.num("slopeobs", slope_obs, 3)
gen.num("slopeint", slope_int, 3)
gen.num("vary", np.var(y_obs), 1)
assert m_oc > np.var(y_obs)                            # worse than predicting the mean
gen.num("mseoo", m_oo, 1)
gen.num("mseoc", m_oc, 1)
gen.num("mseio", m_io, 1)
gen.num("mseic", m_ic, 1)
gen.write()

use_book_style()
fig, axes = plt.subplots(1, 2, figsize=(5.8, 2.6), sharey=True)
xs = np.linspace(x_obs.min(), x_obs.max(), 2)
for ax, x, y, title in [(axes[0], x_obs, y_obs, "(a) observed: X follows illness"),
                        (axes[1], x_int, y_int, "(b) intervened: X set at random")]:
    ax.scatter(x, y, s=3, color=COLORS["accent"], alpha=0.35, linewidths=0)
    ax.plot(xs, icpt_obs + slope_obs * xs, color=COLORS["second"],
            label=f"observational fit, slope {slope_obs:.2f}".replace("-", "\u2212"))
    ax.plot(xs, causal_icpt + theta * xs, color=COLORS["third"], linestyle="--",
            label=f"causal line, slope {theta:.0f}")
    ax.set_xlabel("clinic visits X")
    ax.set_title(title)
axes[0].set_ylabel("health score Y")
axes[1].legend(loc="upper left", frameon=False)
fig.tight_layout()
fig.savefig(figure_path("ch25", "prediction_intervention"))
