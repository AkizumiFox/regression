# Summary and Notes

## Summary

::: {.idea}

1. A model of the conditional mean constrains the whole conditional distribution: in
   a location model every quantile curve is the same curve shifted, and in a
   location-scale model the shape is frozen though the spread
   moves (@prp-qnt-beyond-mean, using @lem-qnt-monotone).

2. The \( \tau \)-quantile is exactly the set of minimizers of the expected check
   loss, for any distribution with a finite mean; the sample problem is a linear
   program, equivariant under scaling, shifts, reparameterization and monotone
   transformation of the response (@thm-qnt-check).

3. A regression quantile interpolates \( p \) observations, and the numbers of
   residual signs satisfy \( n_{-}\le\tau n\le n_{-}+n_0 \) (@prp-qnt-counts).

4. Knight's identity splits the increment of the check loss into a linear term and a
   controllable remainder (@lem-qnt-knight); with the convexity lemma it gives
   asymptotic normality with a sandwich covariance built from the error densities at
   zero, and in the iid case the sparsity function (@thm-qnt-asymptotics).

5. That density can be estimated by a difference quotient of the empirical quantile
   function, or avoided by rank-score tests or the pairs bootstrap; on heteroscedastic
   data only the last two keep their nominal coverage (@exm-qnt-coverage).

6. Fitted quantile curves are ordered at the average covariate vector but can cross
   elsewhere; sorting them can only bring them closer to the truth, in every
   \( L^p \) sense (@prp-qnt-crossing).

7. Squared deviation in place of absolute gives expectiles: unique, increasing in
   \( \tau \), equal to the mean at \( \tau=1/2 \), affinely equivariant, and at the
   quantile level \( \E(\mu-Y)_{+}/\E\lvert Y-\mu\rvert \), which depends on the
   distribution (@prp-qnt-expectile).

8. Expectile regression is weighted least squares with weights determined by the
   signs of the residuals; a step that leaves the weights unchanged has found the
   exact global minimizer (@prp-qnt-iwls).

9. A GAMLSS gives every parameter of a chosen family its own structured additive
   predictor (@def-qnt-gamlss) and maximizes a penalized log-likelihood by block-wise
   Fisher scoring, each step a penalized weighted least squares fit; for the normal
   location-scale model the information is block diagonal (@thm-qnt-gamlss).

10. The likelihood of a flexible location-scale model is unbounded, so the smoothing
    parameters are what make the problem well posed and cannot be chosen by
    maximizing it (@exm-qnt-unbounded).

11. Predictive distributions are compared by proper scoring rules. The continuous
    ranked probability score is strictly proper, splits into a discrepancy and an
    irreducible floor, and is the check loss averaged over all levels; probability
    integral transforms are uniform under a correct model (@prp-qnt-comparison).

12. On Engel's budgets the quantile model and the location-scale GAMLSS beat the
    homoscedastic mean model by about three paired standard errors of the score, and
    are indistinguishable from each other (@exm-qnt-compare).

:::

## Notes and sources

**Quantile regression.** The subject begins with Koenker and Bassett (1978), which
introduced the check loss, the linear programming formulation, the counting
properties and the asymptotics under identically distributed errors, and with Koenker
and Bassett (1982), which brought Engel's budget data into the literature as a test
bed for heteroscedasticity. Koenker (2005) is the standard monograph and the source
for everything in [Section 45.3](03-inference.html) stated here without proof;
Koenker and Hallock (2001) is a short introduction. The median regression ancestor is
far older, as
[Section 5.9](../ch05-model-and-least-squares/09-summary-and-notes.html) records.

The proof in @thm-qnt-asymptotics follows the convexity route of Pollard (1991),
using the exact identity of Knight (1998); the convexity lemma is the one ingredient
taken on trust. The sparsity estimate and its bandwidth are Hall and Sheather (1988);
rank-score inference Gutenbrunner and Jurečková (1992); bootstrap validity Hahn
(1995); computation Koenker and d'Orey (1987) and Portnoy and Koenker (1997). Crossing
and its repairs: He (1997), Bondell, Reich and Wang (2010), and the rearrangement of
Chernozhukov, Fernández-Val and Galichon (2010), whose \( L^p \)-contraction property
is @prp-qnt-crossing(b). The statistic \( R_1(\tau) \) is Koenker and Machado (1999),
nonparametric extensions start from Koenker, Ng and Portnoy (1994), and the Bayesian
route through the asymmetric Laplace is Yu and Moyeed (2001), made practical by
Kozumi and Kobayashi (2011).

**Expectiles.** Newey and Powell (1987) introduced asymmetric least squares, proved
the properties collected in @prp-qnt-expectile, and gave the iterated weighted least
squares algorithm; Efron (1991) developed the regression percentile view and Jones
(1994) the sense in which expectiles are quantiles of a related distribution. The
revival came with smoothing: Schnabel and Eilers (2009), Sobotka and Kneib (2012).
The risk-measurement interest follows Gneiting (2011) and Ziegel (2016).

**GAMLSS.** The framework is Rigby and Stasinopoulos (2005), with software in
Stasinopoulos and Rigby (2007); the ancestor is the LMS method of Cole and Green
(1992). The Bayesian version is Klein, Kneib, Lang and Sohn (2015) and the boosting
version Mayr, Fenske, Hofner, Kneib and Schmid (2012); randomized quantile residuals
are Dunn and Smyth (1996). The survey chapter of Fahrmeir, Kneib, Lang and Marx
(2021) covers all three approaches and supplied this chapter's coverage checklist.

**Scoring rules.** The continuous ranked probability score dates to Matheson and
Winkler (1976); Gneiting and Raftery (2007) is the definitive treatment of propriety,
and Gneiting, Balabdaoui and Raftery (2007) formulated sharpness subject to
calibration. The probability integral transform goes back to Rosenblatt (1952) and
entered forecast assessment through the prequential principle of Dawid (1984).
Comparing predictions by a paired test of their score difference is Diebold and
Mariano (1995).

**Data.** Engel's 1857 budget survey is public domain and distributed with
statsmodels. It appeared in chapters 20 to 23 as an example of a difficulty; here it
is the subject, and nothing in this chapter is a claim about nineteenth-century
Belgium.

## What the ladder reached, and what lies beyond it

This book began with a projection onto a subspace and ends with an estimated
conditional distribution. Each part removed one assumption: Parts I to VI kept the
linear predictor, the normal error and the constant variance and asked what could be
proved; Part VII dropped independence, Part VIII normality, Part IX linearity, and
Part X the mean.

What holds the sequence together is that almost nothing was thrown away at each step.
The normal equations became estimating equations; the projection became a smoother
matrix and then a weighted smoother matrix; the residual sum of squares became a
deviance and then a penalized log-likelihood; the trace of the hat matrix became
effective degrees of freedom, and it is still a trace. A reader who has followed the
argument to here has not learned forty-five methods but one method, restated.

Three directions are left untouched, and it is more useful to name them than to
gesture at them.

**Prediction without a model.** Random forests, gradient boosting and neural networks
predict better than any model in this book on many problems, and they can be made
distributional: quantile regression forests (Meinshausen 2006) return a conditional
distribution rather than a conditional mean. What they do not return is a coefficient
with an interpretation. Conformal prediction (Vovk, Gammerman and Shafer 2005) is the
counterweight, wrapping any predictor in a prediction set with finite-sample coverage
under exchangeability alone — a stronger guarantee than anything proved here, and a
weaker statement, because it says nothing about the mechanism.

**Causal inference.** [Chapter 25](../ch25-causal-interpretation/index.html) went as
far as a regression book can. The modern subject — instrumental variables
beyond @def-eiv-instrument, difference-in-differences, regression discontinuity,
synthetic controls, the semiparametric theory behind doubly robust estimation — is a
separate discipline with its own books, linked to this one by the
Frisch–Waugh–Lovell theorem (@thm-proj-fwl), the engine inside double machine
learning.

**Data that are not vectors, and data that are not complete.** Functional data
analysis (Ramsay and Silverman 2005) rebuilds the linear model when the response or
the covariate is a curve, with the penalized splines of
[Chapter 43](../ch43-smoothing/index.html) as its tool and the same trace as its
degrees of freedom. Conditional transformation models (Hothorn, Kneib and Bühlmann
2014) learn the map to a reference law, a fourth route alongside the three of
[Section 45.1](01-beyond-the-mean.html). Regression for event times has its own
century of theory beginning with Cox (1972); and censored quantile regression (Powell
1986) belongs here too, its point being that quantiles are equivariant under the
censoring map \( y\mapsto\max(y,c) \) by @thm-qnt-check(d) while means are not.

None of these is promised in a later chapter, because there is no later chapter. What
this book offers is the apparatus with which to read them.

## References

- Bondell, Howard D., Reich, Brian J. and Wang, Huixia (2010). Noncrossing Quantile Regression Curve Estimation. *Biometrika* 97(4), 825–838.
- Chernozhukov, Victor, Fernández-Val, Iván and Galichon, Alfred (2010). Quantile and Probability Curves without Crossing. *Econometrica* 78(3), 1093–1125.
- Cole, Tim J. and Green, Peter J. (1992). Smoothing Reference Centile Curves: The LMS Method and Penalized Likelihood. *Statistics in Medicine* 11(10), 1305–1319.
- Cox, David R. (1972). Regression Models and Life-Tables. *Journal of the Royal Statistical Society, Series B* 34(2), 187–220.
- Dawid, A. Philip (1984). Statistical Theory: The Prequential Approach. *Journal of the Royal Statistical Society, Series A* 147(2), 278–292.
- Diebold, Francis X. and Mariano, Roberto S. (1995). Comparing Predictive Accuracy. *Journal of Business and Economic Statistics* 13(3), 253–263.
- Dunn, Peter K. and Smyth, Gordon K. (1996). Randomized Quantile Residuals. *Journal of Computational and Graphical Statistics* 5(3), 236–244.
- Efron, Bradley (1991). Regression Percentiles Using Asymmetric Squared Error Loss. *Statistica Sinica* 1(1), 93–125.
- Fahrmeir, Ludwig, Kneib, Thomas, Lang, Stefan and Marx, Brian D. (2021). *Regression: Models, Methods and Applications*. 2nd edition. Berlin: Springer.
- Gneiting, Tilmann (2011). Making and Evaluating Point Forecasts. *Journal of the American Statistical Association* 106(494), 746–762.
- Gneiting, Tilmann, Balabdaoui, Fadoua and Raftery, Adrian E. (2007). Probabilistic Forecasts, Calibration and Sharpness. *Journal of the Royal Statistical Society, Series B* 69(2), 243–268.
- Gneiting, Tilmann and Raftery, Adrian E. (2007). Strictly Proper Scoring Rules, Prediction, and Estimation. *Journal of the American Statistical Association* 102(477), 359–378.
- Gutenbrunner, Cornelius and Jurečková, Jana (1992). Regression Rank Scores and Regression Quantiles. *The Annals of Statistics* 20(1), 305–330.
- Hahn, Jinyong (1995). Bootstrapping Quantile Regression Estimators. *Econometric Theory* 11(1), 105–121.
- Hall, Peter and Sheather, Simon J. (1988). On the Distribution of a Studentized Quantile. *Journal of the Royal Statistical Society, Series B* 50(3), 381–391.
- He, Xuming (1997). Quantile Curves without Crossing. *The American Statistician* 51(2), 186–192.
- Hothorn, Torsten, Kneib, Thomas and Bühlmann, Peter (2014). Conditional Transformation Models. *Journal of the Royal Statistical Society, Series B* 76(1), 3–27.
- Jones, M. Chris (1994). Expectiles and M-quantiles are Quantiles. *Statistics and Probability Letters* 20(2), 149–153.
- Klein, Nadja, Kneib, Thomas, Lang, Stefan and Sohn, Alexander (2015). Bayesian Structured Additive Distributional Regression with an Application to Regional Income Inequality in Germany. *The Annals of Applied Statistics* 9(2), 1024–1052.
- Knight, Keith (1998). Limiting Distributions for L1 Regression Estimators under General Conditions. *The Annals of Statistics* 26(2), 755–770.
- Koenker, Roger (2005). *Quantile Regression*. Cambridge: Cambridge University Press.
- Koenker, Roger and Bassett, Gilbert (1978). Regression Quantiles. *Econometrica* 46(1), 33–50.
- Koenker, Roger and Bassett, Gilbert (1982). Robust Tests for Heteroscedasticity Based on Regression Quantiles. *Econometrica* 50(1), 43–61.
- Koenker, Roger and d'Orey, Vasco (1987). Computing Regression Quantiles. *Applied Statistics* 36(3), 383–393.
- Koenker, Roger and Hallock, Kevin F. (2001). Quantile Regression. *Journal of Economic Perspectives* 15(4), 143–156.
- Koenker, Roger and Machado, José A. F. (1999). Goodness of Fit and Related Inference Processes for Quantile Regression. *Journal of the American Statistical Association* 94(448), 1296–1310.
- Koenker, Roger, Ng, Pin and Portnoy, Stephen (1994). Quantile Smoothing Splines. *Biometrika* 81(4), 673–680.
- Kozumi, Hideo and Kobayashi, Genya (2011). Gibbs Sampling Methods for Bayesian Quantile Regression. *Journal of Statistical Computation and Simulation* 81(11), 1565–1578.
- Matheson, James E. and Winkler, Robert L. (1976). Scoring Rules for Continuous Probability Distributions. *Management Science* 22(10), 1087–1096.
- Mayr, Andreas, Fenske, Nora, Hofner, Benjamin, Kneib, Thomas and Schmid, Matthias (2012). Generalized Additive Models for Location, Scale and Shape for High Dimensional Data: A Flexible Approach Based on Boosting. *Journal of the Royal Statistical Society, Series C* 61(3), 403–427.
- Meinshausen, Nicolai (2006). Quantile Regression Forests. *Journal of Machine Learning Research* 7, 983–999.
- Newey, Whitney K. and Powell, James L. (1987). Asymmetric Least Squares Estimation and Testing. *Econometrica* 55(4), 819–847.
- Pollard, David (1991). Asymptotics for Least Absolute Deviation Regression Estimators. *Econometric Theory* 7(2), 186–199.
- Portnoy, Stephen and Koenker, Roger (1997). The Gaussian Hare and the Laplacian Tortoise: Computability of Squared-Error versus Absolute-Error Estimators. *Statistical Science* 12(4), 279–300.
- Powell, James L. (1986). Censored Regression Quantiles. *Journal of Econometrics* 32(1), 143–155.
- Ramsay, James O. and Silverman, Bernard W. (2005). *Functional Data Analysis*. 2nd edition. New York: Springer.
- Rigby, Robert A. and Stasinopoulos, D. Mikis (2005). Generalized Additive Models for Location, Scale and Shape. *Journal of the Royal Statistical Society, Series C* 54(3), 507–554.
- Rosenblatt, Murray (1952). Remarks on a Multivariate Transformation. *The Annals of Mathematical Statistics* 23(3), 470–472.
- Schnabel, Sabine K. and Eilers, Paul H. C. (2009). Optimal Expectile Smoothing. *Computational Statistics and Data Analysis* 53(12), 4168–4177.
- Sobotka, Fabian and Kneib, Thomas (2012). Geoadditive Expectile Regression. *Computational Statistics and Data Analysis* 56(4), 755–767.
- Stasinopoulos, D. Mikis and Rigby, Robert A. (2007). Generalized Additive Models for Location Scale and Shape (GAMLSS) in R. *Journal of Statistical Software* 23(7), 1–46.
- Vovk, Vladimir, Gammerman, Alexander and Shafer, Glenn (2005). *Algorithmic Learning in a Random World*. New York: Springer.
- Yu, Keming and Moyeed, Rana A. (2001). Bayesian Quantile Regression. *Statistics and Probability Letters* 54(4), 437–447.
- Ziegel, Johanna F. (2016). Coherence and Elicitability. *Mathematical Finance* 26(4), 901–918.
