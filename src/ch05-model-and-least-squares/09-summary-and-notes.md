# Summary and Notes

## Summary

::: {.idea}

1. The regression function \( m(\bm x)=\E(Y\mid\bm X=\bm x) \) is the best predictor in mean
   squared error, and its error is automatically uncorrelated with every function of the
   regressors (@prp-lm-error-decomposition). A regression model is a claim about \( m \) and
   about the distribution around it, made in layers: linearity, constant variance,
   uncorrelated errors, normality ([Section 5.1](01-what-a-model-claims.html)).

2. The least squares line has slope \( S_{xy}/S_{xx} \) and passes through the point of means
   (@thm-lm-simple-ls, @prp-lm-simple-fit). In standard units its slope is the correlation,
   which is regression toward the mean (@eq-lm-toward-mean).

3. The linear model \( \Y=\X\bbeta+\be \), \( \E(\be)=\bzero \), \( \Cov(\be)=\sigma^2\I \), is linear
   in \( \bbeta \), not in the regressors. Polynomials, harmonics, transformations, indicator
   variables and interactions are all linear models (@def-lm-linear-model).

4. For a full-rank model matrix, least squares has the unique solution
   \( \hbeta=(\X\T\X)^{-1}\X\T\y \) of the normal equations \( \X\T\X\bb=\X\T\y \), and the criterion
   exceeds its minimum by \( (\bb-\hbeta)\T\X\T\X(\bb-\hbeta) \) (@thm-lm-ls-full-rank). The residuals
   are orthogonal to every column (@prp-lm-fit-algebra). Slopes can be computed from centred
   data (@prp-lm-centred), and indicator coefficients are differences of group means
   (@prp-lm-indicator-means).

5. Under the second-moment assumptions, \( \hbeta \) is unbiased with covariance
   \( \sigma^2(\X\T\X)^{-1} \), and so is every linear function \( \bm a\T\hbeta \)
   (@thm-lm-moments, @cor-lm-linear-functions). Unbiasedness needs only the mean assumption;
   a wrong mean produces bias and a wrong covariance produces the sandwich formula
   (@prp-lm-misspecified).

6. \( s^2=\text{SSE}/(n-p) \) is unbiased for \( \sigma^2 \), because the residuals are the errors
   projected into a space with \( n-p \) degrees of freedom (@thm-lm-sigma2). Its variance
   depends on the kurtosis of the errors (@prp-lm-var-s2).

7. A coefficient is the difference in mean response per unit of its regressor between profiles
   that agree on the *other regressors in the model* (@def-lm-coefficient-interpretation). It
   changes with the units (@prp-lm-reparameterization) and with the set of other regressors,
   by exactly the amount given by the omitted-variable formula (@prp-lm-omitted). It is not by
   itself a causal effect.

8. Everything later in the book changes one of three components of this model: the
   distribution of the response, the form of the predictor, or the link between predictor and
   mean ([Section 5.8](08-panorama.html)).

:::

## Notes and sources

**History.** Least squares was published by Legendre (1805) as a method for fitting orbits,
and claimed by Gauss (1809), who gave it a probabilistic justification through the normal
distribution. Half a century earlier Boscovich had fitted a line to measurements of meridian arcs by making
the absolute deviations as small as possible, an idea Laplace took up and that returns in this
book as quantile regression. The priority dispute and the early history are told by Stigler (1981; 1986).
The word *regression* comes from Galton (1886), whose study of the heights of parents and
children found the phenomenon now called regression toward the mean. For Galton the slope of a
fitted line measured an inherited tendency. Only later did the least squares machinery of the
astronomers and the regression of the biometricians become one subject.

**The source books.** The chapter covers the material of Rencher and Schaalje (2008,
chapter 1, chapter 6 and sections 7.1–7.7) and of Sen and Srivastava (1990, chapters 1 and 2),
whose presentations of simple and multiple regression by calculus and matrix algebra are the
closest to the one here. Sen and Srivastava's emphasis on the Gauss–Markov conditions as a
benchmark, and their discussion of centring and scaling, are reflected in
[Section 5.1](01-what-a-model-claims.html) and [Section 5.7](07-interpreting-coefficients.html).
Rencher and Schaalje also treat the geometry, the normal model and \( R^2 \) in their chapter 7;
in this book those topics are in [Chapter 6](../ch06-projections/index.html),
[Chapter 7](../ch07-optimality/index.html) and [Chapter 9](../ch09-sums-of-squares/index.html).
The view of a model as a distribution, a linear predictor and a link, and the discussion of
indicator variables and model matrices, follow Agresti (2015, chapter 1). The panorama of
[Section 5.8](08-panorama.html) is organized like the survey in Fahrmeir, Kneib, Lang and Marx
(2021, chapter 2), whose chapter 3 also discusses what the standard assumptions mean in applications and how
covariates are coded. Seber and Lee (2003, chapter 3) and Christensen (2020,
chapters 1–2) give terser, more general treatments that start from projections.

**Fixed and random regressors.** Conditioning on the regressors is standard in linear-model
theory. Freedman (2009) is a careful discussion of what regression models assume and of the gap
between association and causation that [Section 5.1](01-what-a-model-claims.html) and
[Section 5.7](07-interpreting-coefficients.html) point to.

**Standard errors under departures.** The sandwich formula @eq-lm-sandwich, with the unknown
covariance matrix estimated from squared residuals, is due to White (1980) in econometrics; it is
developed in Chapter 21.

**Beyond the linear model.** Generalized linear models were introduced by Nelder and
Wedderburn (1972) and quantile regression by Koenker and Bassett (1978). Additive models were
developed and popularized by Hastie and Tibshirani (1986; 1990), and generalized additive models
for location, scale and shape were proposed by Rigby and Stasinopoulos (2005). Each is the
subject of a later part of the book.

**Data.** The Niño 1+2 sea surface temperatures are from the US National Oceanic and
Atmospheric Administration. The carbon dioxide record is that of Keeling and Whorf (2004). The
Danish money demand data are those analysed by Johansen and Juselius (1990). The RAND Health
Insurance Experiment is described by Manning et al. (1987); the extract used here is the one
analysed by Cameron and Trivedi (2005). All four data sets are distributed with the Python
package statsmodels.

## References

- Agresti, Alan (2015). *Foundations of Linear and Generalized Linear Models*. Hoboken, NJ: Wiley.
- Cameron, A. Colin and Trivedi, Pravin K. (2005). *Microeconometrics: Methods and Applications*. Cambridge: Cambridge University Press.
- Christensen, Ronald (2020). *Plane Answers to Complex Questions: The Theory of Linear Models*. 5th edition. Cham: Springer.
- Fahrmeir, Ludwig, Kneib, Thomas, Lang, Stefan and Marx, Brian D. (2021). *Regression: Models, Methods and Applications*. 2nd edition. Berlin: Springer.
- Freedman, David A. (2009). *Statistical Models: Theory and Practice*. Revised edition. Cambridge: Cambridge University Press.
- Galton, Francis (1886). Regression towards Mediocrity in Hereditary Stature. *Journal of the Anthropological Institute of Great Britain and Ireland* 15, 246–263.
- Gauss, Carl Friedrich (1809). *Theoria Motus Corporum Coelestium in Sectionibus Conicis Solem Ambientium*. Hamburg: Perthes und Besser.
- Hastie, Trevor J. and Tibshirani, Robert J. (1986). Generalized Additive Models. *Statistical Science* 1(3), 297–310.
- Hastie, Trevor J. and Tibshirani, Robert J. (1990). *Generalized Additive Models*. London: Chapman and Hall.
- Johansen, Søren and Juselius, Katarina (1990). Maximum Likelihood Estimation and Inference on Cointegration, with Applications to the Demand for Money. *Oxford Bulletin of Economics and Statistics* 52(2), 169–210.
- Keeling, Charles D. and Whorf, Timothy P. (2004). Atmospheric CO2 Concentrations Derived from Flask Air Samples at Sites in the SIO Network. In *Trends: A Compendium of Data on Global Change*. Oak Ridge, TN: Carbon Dioxide Information Analysis Center, Oak Ridge National Laboratory.
- Koenker, Roger and Bassett, Gilbert (1978). Regression Quantiles. *Econometrica* 46(1), 33–50.
- Legendre, Adrien-Marie (1805). *Nouvelles méthodes pour la détermination des orbites des comètes*. Paris: Firmin Didot.
- Manning, Willard G., Newhouse, Joseph P., Duan, Naihua, Keeler, Emmett B., Leibowitz, Arleen and Marquis, M. Susan (1987). Health Insurance and the Demand for Medical Care: Evidence from a Randomized Experiment. *American Economic Review* 77(3), 251–277.
- Nelder, John A. and Wedderburn, Robert W. M. (1972). Generalized Linear Models. *Journal of the Royal Statistical Society, Series A* 135(3), 370–384.
- Rencher, Alvin C. and Schaalje, G. Bruce (2008). *Linear Models in Statistics*. 2nd edition. Hoboken, NJ: Wiley.
- Rigby, Robert A. and Stasinopoulos, D. Mikis (2005). Generalized Additive Models for Location, Scale and Shape. *Journal of the Royal Statistical Society, Series C (Applied Statistics)* 54(3), 507–554.
- Seber, George A. F. and Lee, Alan J. (2003). *Linear Regression Analysis*. 2nd edition. Hoboken, NJ: Wiley.
- Sen, Ashish and Srivastava, Muni (1990). *Regression Analysis: Theory, Methods, and Applications*. New York: Springer.
- Stigler, Stephen M. (1981). Gauss and the Invention of Least Squares. *The Annals of Statistics* 9(3), 465–474.
- Stigler, Stephen M. (1986). *The History of Statistics: The Measurement of Uncertainty before 1900*. Cambridge, MA: Belknap Press of Harvard University Press.
- White, Halbert (1980). A Heteroskedasticity-Consistent Covariance Matrix Estimator and a Direct Test for Heteroskedasticity. *Econometrica* 48(4), 817–838.
