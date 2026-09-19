# Confidence Intervals, Regions and Prediction Bands

A test asks whether one value of a parameter is compatible with the data; an interval reports all such
values, and how precisely they are pinned down. Under normal errors the answers are exact for every sample
size, and they come from three facts: the fitted vector is normal, the residual sum of squares is an
independent chi-squared variable, and their ratio gives \( t \) and \( F \) pivots. One estimable function gives a
\( t \) interval; several give an ellipsoid, the dual of the \( F \) test. A future observation gives a prediction
interval, whose width never falls below the noise and whose validity rests on the normality of a single error.
The whole regression surface gives the Working–Hotelling band, from the Cauchy–Schwarz inequality. Inverting a
prediction band gives a calibration set, sometimes an interval, sometimes two rays, sometimes the whole line,
and no honest procedure can avoid the last two. Finally, in the conjugate model of
[Chapter 7](../ch07-optimality/index.html), credible and predictive intervals have the same shapes, and under the
flat prior they are the frequentist ones.

**What you need.** The sampling distributions of the normal linear model
([Chapter 7](../ch07-optimality/index.html): @thm-opt-sampling, @cor-opt-t and @cor-opt-quadratic) and
the conjugate posterior (@thm-opt-bayes-conjugate, @lem-opt-nig-marginal, @cor-opt-flat-prior).
Projections, leverage and the Frisch–Waugh–Lovell theorem ([Chapter 6](../ch06-projections/index.html):
@thm-proj-M-formula, @prp-proj-leverage, @thm-proj-fwl). Estimability
([Chapter 8](../ch08-estimability/index.html): @def-est-estimable, @thm-est-characterization) and the sum of
squares of a linear hypothesis ([Chapter 9](../ch09-sums-of-squares/index.html): @thm-ss-hypothesis).
Chi-squared and \( F \) distributions of quadratic forms ([Chapter 4](../ch04-quadratic-forms/index.html):
@thm-qf-chisq, @def-qf-noncentral-f, @def-qf-noncentral-t). The tests of
[Chapter 11](../ch11-general-linear-hypothesis/index.html) are the duals of the intervals here, but the
proofs below do not depend on them.

## Roadmap

- [Intervals for estimable functions](01-estimable-intervals.html)
- [Confidence ellipsoids](02-ellipsoids.html)
- [Prediction intervals](03-prediction.html)
- [Simultaneous bands for the regression line and surface](04-bands.html)
- [Inverse prediction and calibration](05-calibration.html)
- [Bayesian thread II: credible intervals and predictive distributions](06-bayes.html)
- [Summary and notes](07-summary-and-notes.html)
