# Summary and Notes

## Summary

::: {.idea}

1. Additivity is the structural assumption that makes nonparametric regression possible in
   more than one or two dimensions. The model is @def-add-model; its components are fixed
   only up to constants until centred (@prp-add-identifiability), a constraint best absorbed
   into the basis.

2. In a basis, an additive model is one penalized regression with a block-diagonal penalty:
   @eq-add-criterion and the penalized normal equations @eq-add-normal. The fit is a linear
   smoother, its effective degrees of freedom are \( \tr\bS \), and they split over the
   terms (@prp-add-df).

3. Backfitting — smooth each term against the partial residual, and cycle — is block
   Gauss–Seidel on those equations and block coordinate descent on the criterion, converging
   geometrically to the unique fit when the penalized cross-product is positive definite and
   to *a* minimizer otherwise (@thm-add-backfitting); unpenalized, it is alternating
   projection (@cor-add-projection).

4. Concurvity is collinearity between function spaces (@def-add-concurvity-index). It
   leaves the fitted mean alone, inflates the variance of the components by
   \( 1/(1-\kappa_j) \) in the worst direction (@prp-add-concurvity), and is invisible in a
   correlation matrix.

5. A varying-coefficient term \( f(z)u \) is an ordinary penalized term whose design is the
   basis scaled row by row (@def-add-varying, @prp-add-varying-fit): smooth interactions,
   factor-by-curve models and time-varying coefficients at no extra cost.

6. A spatial effect is a term too (@def-add-spatial): for areal data the penalty
   @eq-add-mrf, with one null direction per connected component and a conditional reading as
   shrinkage towards the neighbour mean (@prp-add-mrf); for coordinates a kriging basis
   penalized by the correlation matrix at the knots (@prp-add-kriging), the mixed model of
   Part VII in disguise.

7. Structured additive regression (@def-add-star) is the whole catalogue in one predictor.
   It is identified whenever the *unpenalized parts* of its terms
   are (@prp-add-star-identified), it is a linear mixed model after one orthogonal
   transformation (@prp-add-mixed), and its smoothing parameters are variance ratios
   estimated by @eq-add-reml-update (@prp-add-reml-update).

8. With a link and an exponential dispersion family the fit is penalized
   IRLS (@thm-add-gam): Chapter 34's algorithm with the penalty added. Degrees of freedom,
   dispersion and Chapter 39's diagnostics carry over, and the dispersion belongs in the
   smoothing-parameter update or every curve is undersmoothed.

9. Inference for a smooth term is approximate, and doubly so because the smoothing
   parameters were chosen from the same data. The \( p \)-values rank terms; they do not
   test them.

:::

## Notes and sources

**Additive models and backfitting.**  The additive model and backfitting come from Friedman
and Stuetzle (1981), where the cycling served projection pursuit regression, and from
Breiman and Friedman (1985). The definitive analysis is Buja, Hastie and Tibshirani (1989):
the Gauss–Seidel reading, the convergence theory for symmetric shrinking smoothers,
concurvity, and the distinction between convergence of the fit and of the components are all
theirs. @thm-add-backfitting is that theory for penalized smoothers, where the argument
reduces to the Ostrowski–Reich theorem (Golub and Van Loan, 2013, chapter 11); Hastie and
Tibshirani (1990) made the models widely used. The order here reverses the historical one:
because [Chapter 43](../ch43-smoothing/index.html) built its smoothers from bases and
penalties, the direct solve is available from the start and backfitting becomes an algorithm
for large problems rather than the definition of the estimator. Wood (2017) takes the same
view and is the standard modern reference for the chapter. A fourth route, met in
[Chapter 30](../ch30-regularization-boosting/index.html) rather than here, is componentwise
boosting: fit a block to the current residual as backfitting does, but damp the step and pick
the block greedily instead of cycling (@thm-reg-boosting).

**Concurvity.**  The word and the phenomenon are from Buja, Hastie and Tibshirani (1989);
measuring it by canonical correlations between the component spaces, as in
@def-add-concurvity-index, and reporting it routinely follow Wood (2017). Ramsay, Burnett and
Krewski (2003) document how badly it distorts component estimates in air-pollution
time-series studies, where pollution, temperature and season are nearly functions of one
another.

**Varying coefficients.**  Hastie and Tibshirani (1993) introduced these models in the
generality used here; the tensor-product construction is developed in Wood (2017).

**Spatial terms.**  The Markov random field penalty and its conditional reading are due to
Besag (1974); the intrinsic form used here, with the convolution model that adds an
unstructured region effect, is Besag, York and Mollié (1991), and Rue and Held (2005) is the
systematic treatment. The kriging side goes back to Matheron (1963); the low-rank penalized
version of @prp-add-kriging, and the name *geoadditive model*, are from Kammann and Wand
(2003), and Wendland (1995) constructed the compactly supported radial functions used in the
example. Reich, Hodges and Zadnik (2006) sharpen the spatial-confounding warning: a spatial
random effect can move a fixed-effect estimate substantially, in a direction that is hard to
predict.

**Structured additive regression.**  The name and the framework are from Fahrmeir, Kneib and
Lang (2004) and Brezger and Lang (2006); the presentation in Fahrmeir, Kneib, Lang and Marx
(2021, chapter 9) is the one this chapter is answerable to. The mixed-model representation is
older: Verbyla et al. (1999) developed it and Ruppert, Wand and Carroll (2003) built a
book on it. The updates @eq-add-reml-update are
Schall (1991), analysed in their modern form by Wood and Fasiolo (2017) as the
Fellner–Schall method. The fully Bayesian version is Lang and Brezger (2004); the warning
about improper posteriors from vague priors on variance components is Hobert and Casella
(1996), with Gelman (2006) on which proper priors behave well.

**Generalized additive models.**  Local scoring is Hastie and Tibshirani (1986, 1990). The
direct penalized-IRLS form of @thm-add-gam, which is what software now uses, is developed in
Wood (2017, chapter 6); the performance iteration is Gu (1992), and the more stable outer
iteration, which differentiates a properly defined restricted likelihood, is Wood (2011). The
approximate reference distribution quoted without proof in @thm-add-gam(d) is derived, with
simulation evidence for its accuracy and its failure modes, in Wood (2013). The bands come
from the Bayesian argument of Wahba (1983) and Nychka (1988) through @thm-smo-mixed; Marra
and Wood (2012) find their coverage survives averaged over the curve but not pointwise. Term
selection by an extra penalty on the null space is Marra and Wood (2011).

**Two threads the book does not develop.**  Smoothing-spline analysis of variance, Gu (2013),
builds additive and interaction models from reproducing-kernel Hilbert spaces rather than
from bases and difference penalties, reaching the same models by a more elegant and more
demanding route. And the theory of *rates*: the optimal rate quoted in
[Section 44.1](01-additive-models.html) is Stone (1980), and the fact that each additive
component can still be estimated at the one-dimensional rate whatever \( q \) is — so that
additivity defeats the curse of dimensionality rather than merely dodging it — is Stone
(1985), with the sharp version for backfitting in Opsomer and Ruppert (1997). Those proofs
need empirical-process tools this book does not assume.

## References

- Besag, Julian (1974). Spatial Interaction and the Statistical Analysis of Lattice Systems. *Journal of the Royal Statistical Society, Series B* 36(2), 192–236.
- Besag, Julian, York, Jeremy and Mollié, Annie (1991). Bayesian Image Restoration, with Two Applications in Spatial Statistics. *Annals of the Institute of Statistical Mathematics* 43(1), 1–20.
- Breiman, Leo and Friedman, Jerome H. (1985). Estimating Optimal Transformations for Multiple Regression and Correlation. *Journal of the American Statistical Association* 80(391), 580–598.
- Brezger, Andreas and Lang, Stefan (2006). Generalized Structured Additive Regression Based on Bayesian P-Splines. *Computational Statistics and Data Analysis* 50(4), 967–991.
- Buja, Andreas, Hastie, Trevor and Tibshirani, Robert (1989). Linear Smoothers and Additive Models. *The Annals of Statistics* 17(2), 453–510.
- Fahrmeir, Ludwig, Kneib, Thomas and Lang, Stefan (2004). Penalized Structured Additive Regression for Space–Time Data: A Bayesian Perspective. *Statistica Sinica* 14(3), 731–761.
- Fahrmeir, Ludwig, Kneib, Thomas, Lang, Stefan and Marx, Brian D. (2021). *Regression: Models, Methods and Applications*. 2nd edition. Berlin: Springer.
- Friedman, Jerome H. and Stuetzle, Werner (1981). Projection Pursuit Regression. *Journal of the American Statistical Association* 76(376), 817–823.
- Gelman, Andrew (2006). Prior Distributions for Variance Parameters in Hierarchical Models. *Bayesian Analysis* 1(3), 515–534.
- Golub, Gene H. and Van Loan, Charles F. (2013). *Matrix Computations*. 4th edition. Baltimore: Johns Hopkins University Press.
- Gu, Chong (1992). Cross-Validating Non-Gaussian Data. *Journal of Computational and Graphical Statistics* 1(2), 169–179.
- Gu, Chong (2013). *Smoothing Spline ANOVA Models*. 2nd edition. New York: Springer.
- Hastie, Trevor and Tibshirani, Robert (1986). Generalized Additive Models. *Statistical Science* 1(3), 297–310.
- Hastie, Trevor and Tibshirani, Robert (1990). *Generalized Additive Models*. London: Chapman and Hall.
- Hastie, Trevor and Tibshirani, Robert (1993). Varying-Coefficient Models. *Journal of the Royal Statistical Society, Series B* 55(4), 757–796.
- Hobert, James P. and Casella, George (1996). The Effect of Improper Priors on Gibbs Sampling in Hierarchical Linear Mixed Models. *Journal of the American Statistical Association* 91(436), 1461–1473.
- Kammann, E. E. and Wand, M. P. (2003). Geoadditive Models. *Journal of the Royal Statistical Society, Series C* 52(1), 1–18.
- Lang, Stefan and Brezger, Andreas (2004). Bayesian P-Splines. *Journal of Computational and Graphical Statistics* 13(1), 183–212.
- Marra, Giampiero and Wood, Simon N. (2011). Practical Variable Selection for Generalized Additive Models. *Computational Statistics and Data Analysis* 55(7), 2372–2387.
- Marra, Giampiero and Wood, Simon N. (2012). Coverage Properties of Confidence Intervals for Generalized Additive Model Components. *Scandinavian Journal of Statistics* 39(1), 53–74.
- Matheron, Georges (1963). Principles of Geostatistics. *Economic Geology* 58(8), 1246–1266.
- Nychka, Douglas (1988). Bayesian Confidence Intervals for Smoothing Splines. *Journal of the American Statistical Association* 83(404), 1134–1143.
- Opsomer, Jean D. and Ruppert, David (1997). Fitting a Bivariate Additive Model by Local Polynomial Regression. *The Annals of Statistics* 25(1), 186–211.
- Ramsay, Timothy O., Burnett, Richard T. and Krewski, Daniel (2003). The Effect of Concurvity in Generalized Additive Models Linking Mortality to Ambient Particulate Matter. *Epidemiology* 14(1), 18–23.
- Reich, Brian J., Hodges, James S. and Zadnik, Vesna (2006). Effects of Residual Smoothing on the Posterior of the Fixed Effects in Disease-Mapping Models. *Biometrics* 62(4), 1197–1206.
- Rue, Håvard and Held, Leonhard (2005). *Gaussian Markov Random Fields: Theory and Applications*. Boca Raton: Chapman and Hall/CRC.
- Ruppert, David, Wand, M. P. and Carroll, R. J. (2003). *Semiparametric Regression*. Cambridge: Cambridge University Press.
- Schall, Robert (1991). Estimation in Generalized Linear Models with Random Effects. *Biometrika* 78(4), 719–727.
- Stone, Charles J. (1980). Optimal Rates of Convergence for Nonparametric Estimators. *The Annals of Statistics* 8(6), 1348–1360.
- Stone, Charles J. (1985). Additive Regression and Other Nonparametric Models. *The Annals of Statistics* 13(2), 689–705.
- Verbyla, Arūnas P., Cullis, Brian R., Kenward, Michael G. and Welham, Sue J. (1999). The Analysis of Designed Experiments and Longitudinal Data by Using Smoothing Splines. *Journal of the Royal Statistical Society, Series C* 48(3), 269–311.
- Wahba, Grace (1983). Bayesian "Confidence Intervals" for the Cross-Validated Smoothing Spline. *Journal of the Royal Statistical Society, Series B* 45(1), 133–150.
- Wendland, Holger (1995). Piecewise Polynomial, Positive Definite and Compactly Supported Radial Functions of Minimal Degree. *Advances in Computational Mathematics* 4(1), 389–396.
- Wood, Simon N. (2011). Fast Stable Restricted Maximum Likelihood and Marginal Likelihood Estimation of Semiparametric Generalized Linear Models. *Journal of the Royal Statistical Society, Series B* 73(1), 3–36.
- Wood, Simon N. (2013). On p-Values for Smooth Components of an Extended Generalized Additive Model. *Biometrika* 100(1), 221–228.
- Wood, Simon N. (2017). *Generalized Additive Models: An Introduction with R*. 2nd edition. Boca Raton: Chapman and Hall/CRC.
- Wood, Simon N. and Fasiolo, Matteo (2017). A Generalized Fellner–Schall Method for Smoothing Parameter Optimization with Application to Tweedie Location, Scale and Shape Models. *Biometrics* 73(4), 1071–1081.
