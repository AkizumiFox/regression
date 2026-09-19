# Summary and Notes

## Summary

::: {.idea}

1. A strictly increasing transformation of the response acts on linearity, constant variance and normality
   at once. Multiplicative errors make the logarithm exact for the first two, and log-scale coefficients are
   then ratios of quantiles and means (@prp-tr-multiplicative). The mean comes first.

2. The Box–Cox family is continuous in \( \lambda \), increasing, ordered by convexity, and bounded on one side
   except at the logarithm (@def-tr-box-cox and @prp-tr-family).

3. With the geometric-mean normalization the profile log-likelihood of \( \lambda \) is
   \( -\frac n2\log\text{SSE}_z(\lambda) \) plus a constant, which gives the likelihood ratio interval (@thm-tr-box-cox),
   free of the units of \( y \) (@prp-tr-units). The constructed variable locates its stationary points (@prp-tr-constructed).

4. Component-plus-residual plots show the shape of a regressor's effect, distorted by the part of the
   nonlinearity the other regressors reproduce (@prp-tr-cpr). The Box–Tidwell constructed variable gives an
   exact test of a regressor's power, and its iteration is Gauss–Newton (@thm-tr-box-tidwell).

5. \( \int v^{-1/2} \) stabilizes a variance function \( v \), uniquely up to sign and a constant (@thm-tr-variance-stabilizing):
   square roots for counts, arcsine square roots for proportions, logarithms for
   a constant coefficient of variation (@cor-tr-classical).

6. Quantiles pass through a monotone transformation, means do not. The back-transformed fit estimates
   medians; the mean needs a factor, which smearing supplies without normality (@prp-tr-retransformation).
   Prediction intervals and intervals for medians transform exactly (@prp-tr-back-intervals), and
   elasticities follow by differentiating through both scales (@prp-tr-elasticity).

:::

## Notes and sources

**Source books.** Coverage was checked against Sen and Srivastava (1990, chapter 9), Seber and Lee (2003,
section 10.3 and parts of sections 10.4–10.5) and Christensen (2020, section 12.7). Proofs, examples and exercises here are this book's own.

**Power transformations.** Tukey (1957) studied the powers systematically, and Tukey (1977) popularized the
ladder. Box and Cox (1964) introduced the likelihood approach and the geometric-mean normalization, and
Hernandez and Johnson (1980) gave its large-sample theory when no power is exactly right. On whether
inference should allow for estimating \( \lambda \), see Bickel and Doksum (1981), Box and Cox (1982) and Hinkley and
Runger (1984). Carroll and Ruppert (1988) is the standard monograph on transformation and weighting,
including the *transform-both-sides* model, which transforms the response and a known nonlinear mean
function alike. Yeo and Johnson (2000) give a family defined on the whole line.

**Diagnostics.** Constructed variables and their added-variable plots are developed by Atkinson (1985) and
Cook and Weisberg (1982); Atkinson and Cook and Wang (1983) study the influence of single cases on
\( \hat{\lambda} \). Larsen and McCleary (1972) analysed partial residual plots, Mallows (1986) augmented them, and Cook (1993) explained when
they work. Box and Tidwell (1962) proposed the regressor transformation, and Royston and Altman (1994)
fractional polynomials.

**Variance stabilization and retransformation.** Bartlett (1947) reviewed variance-stabilizing
transformations and the spread-versus-level idea, Anscombe (1948) derived the constants \( 3/8 \), and
Grizzle, Starmer and Koch (1969) developed weighted least squares on transformed proportions. The limit
theorems are in van der Vaart (1998). The smearing estimate is due to Duan (1983). Manning (1998) showed
how heteroscedasticity on the log scale biases retransformed means, one reason why generalized linear
models with a log link (Chapter 34) are now the usual alternative for skewed positive responses.

## References

- Anscombe, Francis J. (1948). The Transformation of Poisson, Binomial and Negative-Binomial Data. *Biometrika* 35(3/4), 246–254.
- Atkinson, Anthony C. (1985). *Plots, Transformations, and Regression: An Introduction to Graphical Methods of Diagnostic Regression Analysis*. Oxford: Clarendon Press.
- Bartlett, Maurice S. (1947). The Use of Transformations. *Biometrics* 3(1), 39–52.
- Bickel, Peter J. and Doksum, Kjell A. (1981). An Analysis of Transformations Revisited. *Journal of the American Statistical Association* 76(374), 296–311.
- Box, George E. P. and Cox, David R. (1964). An Analysis of Transformations. *Journal of the Royal Statistical Society, Series B* 26(2), 211–252.
- Box, George E. P. and Cox, David R. (1982). An Analysis of Transformations Revisited, Rebutted. *Journal of the American Statistical Association* 77(377), 209–210.
- Box, George E. P. and Tidwell, Paul W. (1962). Transformation of the Independent Variables. *Technometrics* 4(4), 531–550.
- Carroll, Raymond J. and Ruppert, David (1988). *Transformation and Weighting in Regression*. New York: Chapman and Hall.
- Christensen, Ronald (2020). *Plane Answers to Complex Questions: The Theory of Linear Models*. 5th edition. Cham: Springer.
- Cook, R. Dennis (1993). Exploring Partial Residual Plots. *Technometrics* 35(4), 351–362.
- Cook, R. Dennis and Wang, Paul C. (1983). Transformations and Influential Cases in Regression. *Technometrics* 25(4), 337–343.
- Cook, R. Dennis and Weisberg, Sanford (1982). *Residuals and Influence in Regression*. New York: Chapman and Hall.
- Duan, Naihua (1983). Smearing Estimate: A Nonparametric Retransformation Method. *Journal of the American Statistical Association* 78(383), 605–610.
- Grizzle, James E., Starmer, C. Frank and Koch, Gary G. (1969). Analysis of Categorical Data by Linear Models. *Biometrics* 25(3), 489–504.
- Hernandez, Francisco and Johnson, Richard A. (1980). The Large-Sample Behavior of Transformations to Normality. *Journal of the American Statistical Association* 75(372), 855–861.
- Hinkley, David V. and Runger, George (1984). The Analysis of Transformed Data. *Journal of the American Statistical Association* 79(386), 302–309.
- Larsen, Wayne A. and McCleary, Susan J. (1972). The Use of Partial Residual Plots in Regression Analysis. *Technometrics* 14(3), 781–790.
- Mallows, Colin L. (1986). Augmented Partial Residuals. *Technometrics* 28(4), 313–319.
- Manning, Willard G. (1998). The Logged Dependent Variable, Heteroscedasticity, and the Retransformation Problem. *Journal of Health Economics* 17(3), 283–295.
- Royston, Patrick and Altman, Douglas G. (1994). Regression Using Fractional Polynomials of Continuous Covariates: Parsimonious Parametric Modelling. *Applied Statistics* 43(3), 429–467.
- Seber, George A. F. and Lee, Alan J. (2003). *Linear Regression Analysis*. 2nd edition. Hoboken, NJ: Wiley.
- Sen, Ashish and Srivastava, Muni (1990). *Regression Analysis: Theory, Methods, and Applications*. New York: Springer.
- Tukey, John W. (1957). On the Comparative Anatomy of Transformations. *The Annals of Mathematical Statistics* 28(3), 602–632.
- Tukey, John W. (1977). *Exploratory Data Analysis*. Reading, MA: Addison-Wesley.
- van der Vaart, Aad W. (1998). *Asymptotic Statistics*. Cambridge: Cambridge University Press.
- Yeo, In-Kwon and Johnson, Richard A. (2000). A New Family of Power Transformations to Improve Normality or Symmetry. *Biometrika* 87(4), 954–959.
