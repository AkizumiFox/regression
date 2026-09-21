# Comparing distributional models

Three fits to the same data now sit side by side: a homoscedastic linear model, a
family of quantile regressions, and a location-scale GAMLSS. They estimate different
things, so none can be judged by how well it hits its own target. What they share is
that each delivers, at every covariate value, a predictive distribution for the
response. That is the object to score.

## Scoring rules

::: {#def-qnt-score}
[Scoring rule, propriety]

A **scoring rule** is a function \( S(F,y) \) assigning a penalty when the predictive
distribution \( F \) is issued and the value \( y \) is observed; smaller is better.
Write \( S(F,G)=\E_{Y\sim G}S(F,Y) \) for its expectation when \( Y \) really follows
\( G \). The rule is **proper** relative to a class of distributions if
\[
S(G,G)\le S(F,G)\qquad\text{for all }F,G\text{ in the class},
\]{#eq-qnt-proper}

and **strictly proper** if equality forces \( F=G \). A **consistent scoring
function** for a functional \( T \) is a function \( s(x,y) \) with
\( \E_{G}s\{T(G),Y\}\le\E_{G}s(x,Y) \) for every \( x \); a functional admitting one
is called **elicitable**.
:::

Propriety is the minimal honesty requirement: under a rule that is not proper a
forecaster does better by misreporting, and a comparison of methods rewards whoever
misreports most skilfully. Three consistent scoring functions have already appeared:
squared error for the mean, the check loss for the \( \tau \)-quantile
by @thm-qnt-check(a), and @eq-qnt-expectile-loss for the \( \tau \)-expectile
by @prp-qnt-expectile(a). Each scores a single number; the next result scores a whole
distribution and ties it back to the check loss.

::: {#prp-qnt-comparison}
[The continuous ranked probability score and the probability integral transform]

For a distribution function \( F \) with finite mean and a real \( y \), define
\[
\operatorname{CRPS}(F,y)=\int_{-\infty}^{\infty}\bigl\{F(s)-\mathbf{1}\{y\le s\}\bigr\}^2\,ds .
\]{#eq-qnt-crps}

::: {.enumerate options="label=(\alph*)"}
1. \( \operatorname{CRPS}(F,y) \) is finite, and for any \( G \) with finite mean
   \[
   \E_{Y\sim G}\operatorname{CRPS}(F,Y)
   =\int\bigl\{F(s)-G(s)\bigr\}^2ds
    +\int G(s)\bigl\{1-G(s)\bigr\}\,ds .
   \]{#eq-qnt-crps-decomposition}

   Hence the CRPS is a strictly proper scoring rule on the distributions with finite
   mean.

2. If \( F \) is continuous and strictly increasing with finite mean, then
   \[
   \operatorname{CRPS}(F,y)=2\int_{0}^{1}\rho_{\tau}\bigl\{y-F^{-1}(\tau)\bigr\}\,d\tau :
   \]{#eq-qnt-crps-check}

   the CRPS is the check loss averaged over all levels.

3. If \( Y\sim F \) with \( F \) continuous, then \( U=F(Y) \) is uniform on
   \( (0,1) \). If the predictive distributions \( F_1,\dots,F_n \) are the true
   conditional distributions of independent \( Y_1,\dots,Y_n \) and each is
   continuous, then the **probability integral transforms** \( U_i=F_i(Y_i) \) are
   independent and uniform.
:::

:::

::: {.proof}
(a) The integrand of @eq-qnt-crps is bounded by one and vanishes outside the region
where \( F(s) \) and \( \mathbf{1}\{y\le s\} \) disagree, which lies in the set where
\( \min\{F(s),1-F(s)\}>0 \) together with a bounded interval; a finite mean
makes \( \int\min\{F,1-F\}\,ds<\infty \), so the score is finite. For fixed
\( s \), \( \mathbf{1}\{Y\le s\} \) is Bernoulli with mean \( G(s) \), so
\[
\E_{G}\bigl\{F(s)-\mathbf{1}\{Y\le s\}\bigr\}^2
=\bigl\{F(s)-G(s)\bigr\}^2+G(s)\bigl\{1-G(s)\bigr\} .
\]
Integrating over \( s \) (Tonelli, the integrand being nonnegative)
gives @eq-qnt-crps-decomposition. The second term does not involve \( F \), so the
expectation is minimized exactly when \( \int(F-G)^2=0 \), that is when \( F=G \)
almost everywhere, and two right-continuous distribution functions agreeing almost
everywhere are equal.

(b) Write \( q_{\tau}=F^{-1}(\tau) \) and \( u=F(y) \), and split both sides at
\( y \). For the left side, \( F(s)-\mathbf{1}\{y\le s\} \) equals \( F(s) \) for
\( s<y \) and \( F(s)-1 \) for \( s\ge y \), so
\( \operatorname{CRPS}(F,y)=\int_{-\infty}^{y}F^2+\int_{y}^{\infty}(1-F)^2 \).
Substituting \( s=q_{\tau} \) and integrating by parts,
\[
\begin{gathered}
\int_{-\infty}^{y}F(s)^2ds=\int_0^{u}\tau^2\,dq_{\tau}
   =u^2y-2\int_0^{u}\tau q_{\tau}\,d\tau,\\
\int_{y}^{\infty}\{1-F\}^2ds=2\int_{u}^{1}(1-\tau)q_{\tau}\,d\tau-(1-u)^2y,
\end{gathered}
\]
the boundary terms at \( 0 \) and \( 1 \) vanishing because the mean is finite. For
the right side, \( q_{\tau}\le y \) exactly when \( \tau\le u \), so
\( \rho_{\tau}(y-q_{\tau})=\tau(y-q_{\tau}) \) for \( \tau<u \) and
\( (1-\tau)(q_{\tau}-y) \) for \( \tau>u \). Hence
\[
\begin{gathered}
2\int_0^1\rho_{\tau}(y-q_{\tau})\,d\tau
 =2y\int_0^{u}\tau\,d\tau-2\int_0^{u}\tau q_{\tau}\,d\tau\\
 +\;2\int_{u}^{1}(1-\tau)q_{\tau}\,d\tau-2y\int_{u}^{1}(1-\tau)\,d\tau,
\end{gathered}
\]
and \( 2\int_0^u\tau\,d\tau=u^2 \), \( 2\int_u^1(1-\tau)d\tau=(1-u)^2 \). The two
expressions agree term by term.

(c) For continuous \( F \) and \( u\in(0,1) \), the set \( \{y:F(y)\le u\} \) is the
interval \( (-\infty,q_u] \) up to endpoints of flat stretches, which carry no
probability, so \( P\{F(Y)\le u\}=F(q_u)=u \) by continuity. Independence of the
\( U_i \) follows from independence of the \( Y_i \), each \( F_i \) being a fixed
function once the covariates are given.
:::

Part (b) is the bridge between the two halves of this chapter. A quantile regression
fitted at a dense grid of levels gives a predictive quantile function, evaluated
by @eq-qnt-crps-check with exactly the loss it was fitted with; a GAMLSS gives a
predictive distribution function, evaluated by @eq-qnt-crps. The same number, in the
units of the response, is available for both.

::: {.remark}
[Other proper rules]

The **logarithmic score** \( S(F,y)=-\log f(y) \) is strictly proper and is what
likelihood fitting optimizes, but it needs a density and is unbounded: one
observation where the predictive density is nearly zero dominates everything,
as @exm-qnt-unbounded showed. The CRPS needs no density, is finite for every
observation, and grows only linearly in how far the observation falls from the
forecast where the log score diverges — the safer default across methods. For a
single level, the **pinball loss** \( \rho_{\tau}(y-\hat q_{\tau}) \) reports how
well one quantile is predicted.
:::

## Calibration and sharpness

Part (c) gives the standard graphical check. Compute \( u_i=\hat F_i(y_i) \) for
held-out observations and plot their histogram: under a correct model it is flat, and
departures have a vocabulary. A **hump** in the middle means predictive distributions
that are too wide; a **U shape** too narrow, the common and dangerous failure; a
**slope** a biased location; **spikes at 0 and 1** alone, wrong tails with a sound
body.

For a discrete response \( F_i(Y_i) \) is not uniform, and the histogram must use a
randomized version, drawing \( U_i \) uniformly on \( [F_i(y_i^{-}),F_i(y_i)] \):
the randomized quantile residuals of Dunn and Smyth (1996), exactly normal under a
correct model after applying \( \Phi^{-1} \).

::: {.warning}
[In-sample calibration proves nothing]

The transforms must be computed out of sample. A fitted model is calibrated on its
own training data almost by construction: @prp-bin-calibration showed that in a
logistic regression with an intercept the recalibration slope is exactly one and the
intercept exactly zero, whatever the model's merits. A location-scale fit that has
just matched \( \hat\sigma \) to its own residuals gives a flat in-sample histogram
whether or not its scale curve is right.
:::

Calibration alone is not enough either: the unconditional distribution of the
response is perfectly calibrated and useless. The principle of Gneiting, Balabdaoui
and Raftery (2007) is to maximize **sharpness subject to calibration** — among
predictive distributions that pass the check, prefer the most concentrated. A proper
scoring rule does both at once, which is why the scores come first.

## Three models, one comparison

::: {#exm-qnt-compare}
[Mean model, quantile model and GAMLSS on Engel's budgets]

The \( 235 \) households are split into ten folds. On each fold, three predictive
distributions are built from the other nine:

::: {.enumerate options="label=(\roman*)"}
1. the homoscedastic linear model of Part II, with normal errors and one
   \( \hat\sigma \);

2. quantile regressions at \( 199 \) levels \( \tau=1/200,\dots,199/200 \), sorted at
   each covariate value by @prp-qnt-crossing(b) so that the fitted quantiles are
   ordered;

3. the location-scale GAMLSS of @exm-qnt-engel-gamlss.
:::

All three are scored the same way, by the discretized form of @eq-qnt-crps-check over
that grid of levels, so no model is helped by a formula the others do not get.

| model | CRPS | pinball \( 0.1 \) | pinball \( 0.5 \) | pinball \( 0.9 \) | below \( \hat q_{0.1} \) | below \( \hat q_{0.9} \) |
|---|---|---|---|---|---|---|
| mean model | 59.34 | 20.41 | 38.98 | 21.76 | 0.030 | 0.923 |
| quantile model | 53.00 | 17.69 | 38.05 | 14.93 | 0.089 | 0.902 |
| GAMLSS | 53.12 | 16.26 | 37.53 | 16.49 | 0.140 | 0.885 |

Scores are in francs, and each model gets its own best predictive distribution: for
the mean model the exact \( t \) predictive of
[Section 12.3](../ch12-intervals-and-bands/03-prediction.html). Both distributional
models beat it by about \( 6 \) francs of CRPS per household: the differences are
\( -6.35 \) with a paired standard error of \( 1.75 \), and \( -6.23 \) with
standard error \( 2.12 \). Against each other they differ by \( -0.12 \) with a
paired standard error of \( 1.99 \), which is no difference at all. The mean model
fails where @prp-qnt-beyond-mean predicted, putting only \( 0.030 \) of the
households below its own tenth percentile.

The Kolmogorov–Smirnov distance of the held-out transforms from uniformity is
\( 0.098 \) for the mean model, \( 0.024 \) for the quantile model and \( 0.055 \)
for the GAMLSS ([Figure 45.6.1](06-comparing.html#fig-qnt-pit)). The mean model's hump
is the signature of predictive distributions that are too wide; the GAMLSS histogram
is raised at both ends, the signature of a normal shape where the data want heavier
tails — an actionable criticism, whose answer is a three-parameter family.
:::

::: {when-format="html"}
![**Figure 45.6.1.** Held-out probability integral transforms from ten-fold
cross-validation; a correct model gives a flat histogram (dashed). The mean model is
humped, the quantile model close to flat, the GAMLSS raised at both
ends.](pit_histograms.svg){#fig-qnt-pit width=100%}
:::

::: {when-format="pdf"}
![Held-out probability integral transforms from ten-fold cross-validation; a correct
model gives a flat histogram (dashed). The mean model is humped, the quantile model
close to flat, the GAMLSS raised at both ends.](pit_histograms.pdf){width=100%}
:::

```{.python .run #cell-scores-scores}
import numpy as np
from scipy import stats

LEVELS = np.arange(1, 200) / 200.0          # the tau grid the CRPS integral uses

def check_loss(u, tau):
    """The check function of Section 45.2."""
    return u * (tau - (u < 0))

def crps_from_quantiles(q, y):
    """CRPS = 2 * integral of the check loss over tau, on the grid LEVELS.

    q has one row per observation and one column per level of LEVELS.
    """
    return 2 * np.mean(check_loss(y[:, None] - q, LEVELS[None, :]), axis=1)

def pit_from_quantiles(q, y):
    """The probability integral transform read off a set of fitted quantiles."""
    return np.mean(q <= y[:, None], axis=1)

# one household: a normal predictive with mean 600 and scale 120 francs, observed at 700
q_one = 600.0 + 120.0 * stats.norm.ppf(LEVELS)[None, :]
y_one = np.array([700.0])
print(f"CRPS {crps_from_quantiles(q_one, y_one)[0]:.2f} francs,"
      f" PIT {pit_from_quantiles(q_one, y_one)[0]:.3f}")
```

## Reading a score difference honestly

Three cautions close the comparison, the same three that have followed every model
comparison in this book. A difference needs a standard error: \( 59.34 \) against
\( 53.00 \) is three and a half paired standard errors, \( 53.00 \) against
\( 53.12 \) less than a tenth of one, and the formal test, corrected for dependence
in time series, is Diebold and Mariano (1995). Selecting on the score inflates it:
the GAMLSS smoothing parameters were chosen in @exm-qnt-engel-gamlss on these same
ten folds, so its score is optimistic in the way @prp-sel-selection-bias describes.
And a score is not a purpose: the CRPS weights all levels alike, whereas an insurer
who cares about the upper tail should compare pinball losses at \( \tau=0.95 \) and a
nutritionist concerned with the poorest households at \( \tau=0.05 \). Decide what
the model is for, then choose a scoring function consistent for it.

## Exercises

### A. Check your understanding

::: {#exr-qnt-pit-shapes}
[A1]

A forecaster issues \( \Normal(0,4) \) every day while the truth is
\( \Normal(0,1) \). Sketch the probability integral transform histogram, and say
which of the four shapes listed above appears.
:::

::: {.solution}
\( U=\Phi(Y/2) \) with \( Y\sim\Normal(0,1) \) concentrates near \( 1/2 \), since
\( Y/2 \) rarely leaves \( (-1,1) \) and \( \Phi \) maps that to
\( (0.16,0.84) \): a hump in the middle, the signature of predictive distributions
that are too wide.
:::

### B. Practice

::: {#exr-qnt-crps-normal}
[B1]

Show that for \( F=\Normal(\mu,\sigma^2) \) and \( z=(y-\mu)/\sigma \),
\[
\operatorname{CRPS}(F,y)=\sigma\Bigl\{z\bigl(2\Phi(z)-1\bigr)+2\varphi(z)
   -\tfrac1{\sqrt{\pi}}\Bigr\},
\]
using the identity \( \operatorname{CRPS}(F,y)=\E\lvert X-y\rvert
-\tfrac12\E\lvert X-X'\rvert \) for independent \( X,X'\sim F \).
:::

::: {.solution}
For \( Z\sim\Normal(0,1) \) and a constant \( a \), integrating by parts gives
\( \E\lvert Z-a\rvert=2\varphi(a)+a\{2\Phi(a)-1\} \). Writing \( X=\mu+\sigma Z \),
\( \E\lvert X-y\rvert=\sigma\E\lvert Z-z\rvert
=\sigma\{2\varphi(z)+z(2\Phi(z)-1)\} \). Next \( X-X'\sim\Normal(0,2\sigma^2) \), and
\( \E\lvert\Normal(0,v^2)\rvert=v\sqrt{2/\pi} \), so
\( \E\lvert X-X'\rvert=\sigma\sqrt2\cdot\sqrt{2/\pi}=2\sigma/\sqrt{\pi} \).
Substituting both into the identity gives the stated formula.
By @prp-qnt-comparison(a) the expected score is minimized in \( (\mu,\sigma) \) at the
true values.
:::

::: {#exr-qnt-crps-energy}
[B2]

Prove the identity used in @exr-qnt-crps-normal:
\( \operatorname{CRPS}(F,y)=\E\lvert X-y\rvert-\tfrac12\E\lvert X-X'\rvert \) for
independent \( X,X'\sim F \) with finite mean.
:::

::: {.solution}
The key is \( \lvert a-b\rvert=\int\{\mathbf{1}\{a\le s\}-\mathbf{1}\{b\le s\}\}^2ds \),
since the integrand is one exactly on the interval between \( a \) and \( b \).
Taking \( a=X \), \( b=y \) and expectations, and using
\( \E\{A-c\}^2=\Var(A)+(\E A-c)^2 \) with \( A=\mathbf{1}\{X\le s\} \),
\[
\E\lvert X-y\rvert=\int\Bigl[\bigl\{F(s)-\mathbf{1}\{y\le s\}\bigr\}^2
   +F(s)\{1-F(s)\}\Bigr]ds=\operatorname{CRPS}(F,y)+\int F(1-F)\,ds .
\]
Applying the same identity to the independent pair, \( \E\{A-A'\}^2=2\Var(A) \)
gives \( \E\lvert X-X'\rvert=2\int F(1-F)\,ds \). Subtracting half of the second
display from the first gives the identity; finiteness of the mean makes
\( \int F(1-F)\,ds \) finite.
:::

### C. Going deeper

::: {#exr-qnt-crps-decomposition}
[C1]

Interpret the two terms of @eq-qnt-crps-decomposition. Show that the second depends
only on the truth \( G \), so it is an irreducible floor, and that the first is an
integrated squared error between the forecast and the truth. Relate this to the
bias-variance decomposition of prediction error in
[Section 29.1](../ch29-model-selection/01-prediction-error.html).
:::

