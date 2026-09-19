# Model Selection and Prediction

A fitted model is judged by how well it predicts responses it has not seen, and its error on the responses it was
fitted to flatters it. For any fitting procedure the flattery is twice the covariance between the fit and the data,
\( 2\sigma^2p/n \) for least squares whether the model is right or wrong. Correcting for it gives Mallows' \( C_p \) and AIC.
BIC answers a different question, which model is most probable: it finds a true finite model with probability tending to
one, while AIC predicts almost as well as the best candidate when no finite model is true. Cross-validation estimates
prediction error with fewer assumptions, and for least squares and ridge it needs no refitting. The search over subsets
is unstable, and every estimate and interval computed after it is biased by it, unless the inference is designed to
survive the search.

**What you need.** [Chapter 6](../ch06-projections/index.html) (projections, @thm-proj-nested, and the
Frisch–Waugh–Lovell theorem, @thm-proj-fwl), [Chapter 7](../ch07-optimality/index.html) (the
likelihood of the normal model, @thm-opt-mle), [Chapter 9](../ch09-sums-of-squares/index.html)
(adjusted \( R^2 \), @prp-ss-adjusted-r2), [Chapter 11](../ch11-general-linear-hypothesis/index.html)
(the \( F \) test, @thm-glh-f-test), [Chapter 13](../ch13-multiplicity/index.html) (Scheffé and
Bonferroni, @thm-mc-scheffe and @thm-mc-bonferroni), [Chapter 14](../ch14-correlation-lack-of-fit-prediction/index.html)
(prediction error and its optimism, @prp-cor-optimism and @thm-cor-prediction-error),
[Chapter 19](../ch19-theory-of-departures/index.html) (omitted and irrelevant regressors, @thm-dep-mse) and [Chapter 20](../ch20-residuals-leverage-influence/index.html) (the deletion
formulas and PRESS, @thm-res-deletion). From this part:
[Chapter 26](../ch26-collinearity/index.html) (collinearity) and
[Chapter 27](../ch27-shrinkage/index.html) (ridge regression, thresholding and Stein's lemma, @lem-shr-stein).

## Roadmap

- [Prediction error and optimism](01-prediction-error.html)
- [Mallows' Cp, AIC and BIC](02-cp-aic-bic.html)
- [Cross-validation](03-cross-validation.html)
- [All subsets and stepwise search](04-subset-search.html)
- [Selection bias and honest inference](05-selection-bias.html)
- [Summary and notes](06-summary-and-notes.html)
