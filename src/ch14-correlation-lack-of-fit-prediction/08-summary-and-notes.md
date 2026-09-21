# Summary and notes

## Summary

::: {.idea}

1. If the errors are normal, homoscedastic and independent of the regressors, every \( t \) and \( F \)
   statistic has the same law given every design, so tests and intervals stay exact with random
   regressors of any distribution. Power depends on the design drawn and is averaged over it
   (@lem-cor-conditioning, @thm-cor-conditional).

2. The estimates themselves change law. With normal regressors the slope of a simple regression is a
   scaled \( t(n-1) \) variable with variance \( \sigma^2/(\sigma_X^2(n-3)) \), and with \( k \) regressors the slopes have
   covariance \( \sigma^2\bSigma_{XX}^{-1}/(n-k-2) \) (@prp-cor-slope-unconditional, @lem-cor-inverse-wishart,
   @cor-cor-slope-covariance).

3. The conditional mean is the only predictor whose errors are uncorrelated with every function of the
   regressors, and the correlation ratio bounds the correlation any predictor can reach
   (@prp-cor-best-predictor).

4. A fitted linear predictor pays for estimation. With normal regressors its expected error on a new
   case is \( \sigma^2(1+1/n)(n-2)/(n-k-2) \), while its in-sample error has mean \( \sigma^2(n-k-1)/n \); with a fixed
   design the gap is \( 2r\sigma^2/n \) (@prp-cor-prediction-decomposition, @thm-cor-prediction-error,
   @prp-cor-optimism).

5. When the new error is correlated with the data, the best linear unbiased predictor adds to the
   estimated mean the best linear prediction of that error from the generalized least squares
   residuals (@thm-cor-blup).

6. \( R^2 \) is the maximum likelihood estimate of \( \rho^2_{Y\cdot X} \). Given the design it is a transformed noncentral
   \( F \) variable; with normal regressors its noncentrality is \( \lambda\chi^2(n-1) \), \( \lambda=\rho^2/(1-\rho^2) \), which gives its
   exact law, its bias, the power of the overall \( F \) test and an exact interval for \( \rho^2 \)
   (@prp-cor-mvn-mle, @thm-cor-r2-distribution).

7. The sample correlation has an exact law expressed through noncentral \( t \) probabilities, and it is
   stochastically increasing in \( \rho \). Fisher's \( z=\tanh^{-1}r \) is asymptotically \( \Normal(\zeta,1/n) \) whatever \( \rho \),
   and in practice \( \Normal(\zeta,1/(n-3)) \) (@prp-cor-r-exact, @thm-cor-fisher-z).

8. A sample partial correlation given \( q \) variables behaves exactly like an ordinary correlation from \( n-q \)
   cases, and its test is the \( t \) test of the corresponding regression coefficient
   (@thm-cor-partial-law, @prp-cor-partial-coefficient, @thm-cor-partial-test).

9. With replicated rows, the residual sum of squares splits into pure error and lack of fit, and their
   ratio is an exact \( F \) test of the model against every mean function of the regressors
   (@lem-cor-row-structure, @thm-cor-lack-of-fit).

10. Without replicates, any larger model chosen from \( \X \) alone, any regressors built from the fitted values,
    and any smoother of the residuals give exact tests. Each has power only against departures its
    alternative can represent, and some departures reduce its power below the level
    (@prp-cor-augmented-test, @thm-cor-fitted-regressors, @prp-cor-smoother-null).

:::

## Notes and sources

**Source texts.** The chapter covers Christensen (2020, sections 6.3 to 6.7) on prediction theory,
multiple and partial correlation, best linear unbiased prediction and lack of fit; Rencher and Schaalje
(2008, chapter 10) on regression with random regressors; and Agresti (2015, sections 3.3 and 3.4), whose
worked example includes partial correlations and model checking. Christensen's near-replicate tests and his
view of prediction as the goal of regression shaped [Section 14.2](02-best-linear-prediction.html) and
[Section 14.7](07-near-replicates.html). Rencher and Schaalje give the normal-theory formulas for \( R^2 \), \( z \)
and partial correlations. Agresti's intervals of
section 3.3 belong to [Chapter 12](../ch12-intervals-and-bands/index.html).

**Random regressors and prediction.** Reporting precision conditionally on an ancillary design is Fisher's
conditionality principle; Cox (1958) gave the standard discussion. The Wishart distribution is due to
Wishart (1928), and Anderson (2003) treats it and the mean of its inverse through densities; the proof of
@lem-cor-inverse-wishart here uses only conditioning and symmetry. Breiman and Freedman (1983) use the
prediction error of regression with normal regressors to study how many variables to enter, a question
[Chapter 29](../ch29-model-selection/index.html) takes up. Best linear unbiased prediction is due to Goldberger (1962); [Chapter 31](../ch31-general-gauss-markov/index.html)
and [Chapter 32](../ch32-linear-mixed-models/index.html) develop it for general covariance structures and mixed
models (@thm-mix-blup).

**Correlation coefficients.** Fisher found the exact distribution of \( r \) (1915), introduced \( z \) with its
approximate mean and variance (1921), showed that partial correlations follow the same law with the sample
size reduced by the number of conditioning variables (1924), and derived the distribution of the multiple
correlation (1928). Hotelling (1953) refined the moments of \( r \) and \( z \). @lem-cor-bivariate-moments is a case of
Isserlis's (1918) formula, and the mean of \( R^2 \) in @exm-cor-r2-bias follows from Fisher's (1928) distribution; see Olkin and Pratt (1958), who also give the unbiased estimator of \( \rho^2 \). For the
limit theorems of [Section 14.4](04-correlation-coefficient.html), see van der Vaart (1998).

**Lack of fit.** The pure-error test is Fisher's (1922). Near-replicate tests were proposed by
Shillington (1979) and analysed by Christensen (1989, 1991); the rainbow test is due to Utts (1982).
Ramsey (1969) introduced RESET, and Milliken and Graybill (1970) proved @thm-cor-fitted-regressors, which also
covers Tukey's (1949) one degree of freedom for nonadditivity. Smoothing-based tests were proposed by
Azzalini and Bowman (1993) and are the subject of Hart (1997).

**Data.** The stack loss data are from Brownlee (1965). The RAND Health Insurance Experiment is described by
Manning et al. (1987); the extract and the 2009 state data are public-domain data distributed with statsmodels.

**What is new here.** The conditional derivations of the laws of \( R^2 \), \( r \) and partial correlations (@thm-cor-r2-distribution, @prp-cor-r-exact, @thm-cor-partial-law), the proof of @lem-cor-inverse-wishart,
the organization of lack-of-fit tests around @prp-cor-augmented-test, and the examples and simulations are
this book's.

## References

- Agresti, Alan (2015). *Foundations of Linear and Generalized Linear Models*. Hoboken, NJ: Wiley.
- Anderson, Theodore W. (2003). *An Introduction to Multivariate Statistical Analysis*. 3rd edition. Hoboken, NJ: Wiley.
- Azzalini, Adelchi and Bowman, Adrian W. (1993). On the Use of Nonparametric Regression for Checking Linear Relationships. *Journal of the Royal Statistical Society, Series B* 55(2), 549–557.
- Breiman, Leo and Freedman, David (1983). How Many Variables Should Be Entered in a Regression Equation? *Journal of the American Statistical Association* 78(381), 131–136.
- Brownlee, Kenneth A. (1965). *Statistical Theory and Methodology in Science and Engineering*. 2nd edition. New York: Wiley.
- Christensen, Ronald (1989). Lack-of-Fit Tests Based on Near or Exact Replicates. *The Annals of Statistics* 17(2), 673–683.
- Christensen, Ronald (1991). Small-Sample Characterizations of Near Replicate Lack-of-Fit Tests. *Journal of the American Statistical Association* 86(415), 752–756.
- Christensen, Ronald (2020). *Plane Answers to Complex Questions: The Theory of Linear Models*. 5th edition. Cham: Springer.
- Cox, David R. (1958). Some Problems Connected with Statistical Inference. *The Annals of Mathematical Statistics* 29(2), 357–372.
- Fisher, Ronald A. (1915). Frequency Distribution of the Values of the Correlation Coefficient in Samples from an Indefinitely Large Population. *Biometrika* 10(4), 507–521.
- Fisher, Ronald A. (1921). On the “Probable Error” of a Coefficient of Correlation Deduced from a Small Sample. *Metron* 1(4), 3–32.
- Fisher, Ronald A. (1922). The Goodness of Fit of Regression Formulae, and the Distribution of Regression Coefficients. *Journal of the Royal Statistical Society* 85(4), 597–612.
- Fisher, Ronald A. (1924). The Distribution of the Partial Correlation Coefficient. *Metron* 3, 329–332.
- Fisher, Ronald A. (1928). The General Sampling Distribution of the Multiple Correlation Coefficient. *Proceedings of the Royal Society of London. Series A* 121(788), 654–673.
- Goldberger, Arthur S. (1962). Best Linear Unbiased Prediction in the Generalized Linear Regression Model. *Journal of the American Statistical Association* 57(298), 369–375.
- Hart, Jeffrey D. (1997). *Nonparametric Smoothing and Lack-of-Fit Tests*. New York: Springer.
- Hotelling, Harold (1953). New Light on the Correlation Coefficient and its Transforms. *Journal of the Royal Statistical Society, Series B* 15(2), 193–232.
- Isserlis, Leon (1918). On a Formula for the Product-Moment Coefficient of Any Order of a Normal Frequency Distribution in Any Number of Variables. *Biometrika* 12(1/2), 134–139.
- Manning, Willard G., Newhouse, Joseph P., Duan, Naihua, Keeler, Emmett B., Leibowitz, Arleen and Marquis, M. Susan (1987). Health Insurance and the Demand for Medical Care: Evidence from a Randomized Experiment. *American Economic Review* 77(3), 251–277.
- Milliken, George A. and Graybill, Franklin A. (1970). Extensions of the General Linear Hypothesis Model. *Journal of the American Statistical Association* 65(330), 797–807.
- Olkin, Ingram and Pratt, John W. (1958). Unbiased Estimation of Certain Correlation Coefficients. *The Annals of Mathematical Statistics* 29(1), 201–211.
- Ramsey, James B. (1969). Tests for Specification Errors in Classical Linear Least-Squares Regression Analysis. *Journal of the Royal Statistical Society, Series B* 31(2), 350–371.
- Rencher, Alvin C. and Schaalje, G. Bruce (2008). *Linear Models in Statistics*. 2nd edition. Hoboken, NJ: Wiley.
- Shillington, E. Richard (1979). Testing Lack of Fit in Regression without Replication. *The Canadian Journal of Statistics* 7(2), 137–146.
- Tukey, John W. (1949). One Degree of Freedom for Non-Additivity. *Biometrics* 5(3), 232–242.
- Utts, Jessica M. (1982). The Rainbow Test for Lack of Fit in Regression. *Communications in Statistics. Theory and Methods* 11(24), 2801–2815.
- van der Vaart, Aad W. (1998). *Asymptotic Statistics*. Cambridge: Cambridge University Press.
- Wishart, John (1928). The Generalised Product Moment Distribution in Samples from a Normal Multivariate Population. *Biometrika* 20A(1/2), 32–52.
