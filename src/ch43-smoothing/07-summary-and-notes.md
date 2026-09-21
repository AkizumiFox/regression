# Summary and Notes

## Summary

::: {.idea}

1. A local polynomial fit is weighted least squares solved afresh at every point;
   its weights are linear in \( \y \) and reproduce polynomials of the fitted
   degree exactly, whatever the design and the bandwidth
   (@prp-smo-equivalent-kernel).

2. Under a regular design, a bounded-variation kernel and a twice-differentiable
   \( f \), the local linear estimate has bias \( \tfrac12\mu_2f''(x)h^2+o(h^2) \)
   and variance \( \sigma^2R(K)/\{nhg(x)\} \); the local constant estimate adds
   \( \mu_2h^2f'g'/g \) in the interior and a term of order \( h \) at the
   boundary. The optimal bandwidth is of order \( n^{-1/5} \) and the pointwise
   mean squared error of order \( n^{-4/5} \)
   (@lem-smo-moment-sums, @thm-smo-local-poly).

3. B-splines, defined by the de Boor recurrence, are nonnegative, supported on
   \( d+2 \) consecutive knots, sum to one, are \( C^{d-1} \) at simple knots, have
   a derivative that is again a combination of B-splines, and form a basis of the
   spline space of @thm-ply-spline-space (@def-smo-bspline, @prp-smo-bspline);
   coefficients at the Greville abscissae reproduce the identity
   (@prp-smo-greville). The basis matrix is banded and well conditioned where the
   truncated power basis is neither.

4. A penalized spline takes a rich B-spline basis and a difference penalty
   (@def-smo-pspline). The fit is a generalized ridge regression: it shrinks the
   Demmler–Reinsch coordinates by \( (1+\lambda e_j)^{-1} \), spends
   \( \tr(\bS_\lambda)=\sum_j(1+\lambda e_j)^{-1} \) degrees of freedom, and as
   \( \lambda\to\infty \) tends to the least squares fit in the null space of the
   penalty — a constant for \( k=1 \), a line for \( k=2 \)
   (@lem-smo-diagonalize, @prp-smo-pspline-ridge, @prp-smo-penalty-limit). Once the
   penalty is there, the number of knots stops mattering.

5. Among all functions with square-integrable second derivative, the minimizer of
   \( \sum_i\{y_i-f(x_i)\}^2+\lambda\int f''^2 \) is the natural cubic spline with
   knots at the design points (@lem-smo-natural-interpolant,
   @thm-smo-smoothing-spline). Its fit solves a pentadiagonal system in \( O(n) \)
   operations (@prp-smo-reinsch).

6. Every smoother here is a linear operator \( \bS_\lambda \). Its trace measures
   what was spent, \( \sigma^2\bS_\lambda\bS_\lambda\T \) is the covariance of the
   fit, and the obvious pointwise interval is an interval for
   \( \E\{\hat f(x_i)\} \), not for \( f(x_i) \) (@prp-smo-linear-smoother).

7. The smoothing parameter can be chosen by leave-one-out cross-validation, which
   needs no refitting, by generalized cross-validation, which needs only the
   trace, or by a corrected AIC: they agree in large samples with independent
   errors, differ in the upper tail, and all overfit badly when the errors are
   correlated (@thm-smo-lambda).

8. A quadratic penalty is a normal prior. Splitting the coefficients into the
   penalty's null space and its complement turns a penalized spline into a linear
   mixed model with \( \lambda=\sigma^2/\sigma_u^2 \), so that REML estimates the
   smoothing parameter and the posterior supplies a band
   (@thm-smo-mixed). That band is wider than the naive one by roughly the amount
   needed to absorb the bias on average — stated, not proved,
   in @prp-smo-bayes-bands(c) — which buys across-the-function coverage and not
   pointwise coverage.

:::

## Notes and sources

**Kernels and local polynomials.**  The local constant estimator is due
independently to Nadaraya (1964) and Watson (1964), and the optimality of the
kernel named after him to Epanechnikov (1969), though as @exr-smo-efficiency
shows the gain is a few per cent. The modern theory of local *polynomial*
fitting is Fan and Gijbels (1996), whose chapter 3 has the random-design version
of @thm-smo-local-poly with the conditional expansions handled properly; Wand and
Jones (1995) give the background. The design-bias and boundary-bias comparisons
in @thm-smo-local-poly(b) and (c) are why local linear replaced
Nadaraya–Watson, an argument made forcefully by Hastie and Loader (1993). The
fixed-design proof given here, through @lem-smo-moment-sums, makes every error
term explicit; the price is assumption (A), and @exr-smo-random-design says what
changes without it. That the rate \( n^{-4/5} \) cannot be improved for
twice-differentiable \( f \) is Stone (1980, 1982), whose curse of dimensionality
is what makes [Chapter 44](../ch44-additive-models/index.html) necessary. Loess
is Cleveland (1979), extended by Cleveland and Devlin (1988), and plug-in
bandwidths are Ruppert, Sheather and Wand (1995).

**B-splines.**  Splines as a subject begin with Schoenberg (1946). The stable
recurrence @eq-smo-deboor is due to de Boor (1972) and Cox (1972); de Boor
(2001) is the standard reference and contains Marsden's identity, the
Curry–Schoenberg theorem that the B-splines span the spline space for a general
knot sequence, and the general form of @prp-smo-greville. The proofs here need
only the recurrence, the partition of unity and the derivative formula, and are
stated for the clamped sequence; the induction in @prp-smo-bspline(f) is the
shortest route I know to the basis property for it.

**Penalized splines.**  Penalizing differences of a fitted sequence is
Whittaker's (1923) graduation, rediscovered many times; its marriage with the
B-spline basis is Eilers and Marx (1996), which introduced both the name and the
practice of taking many knots and letting the penalty do the work. O'Sullivan
(1986) had earlier used the integral penalty on a B-spline basis, and Wand and
Ormerod (2008) relate the two. Ruppert (2002) is the
source of the rule of thumb for the knot count and Ruppert, Wand and Carroll
(2003) the book-length treatment. The simultaneous
diagonalization of @lem-smo-diagonalize goes back to Demmler and Reinsch (1975).

**Smoothing splines.**  Reinsch (1967) gave the algorithm; the variational
characterization @thm-smo-smoothing-spline and the integration-by-parts argument
behind it are due to Schoenberg and to Reinsch, and the treatment here follows
Green and Silverman (1994, ch. 2), whose band matrices are those of
@prp-smo-reinsch. Wahba (1990) is the definitive account from the
reproducing-kernel side, with the representer theorem of Kimeldorf and Wahba
(1971) at its centre; Eubank (1999) is gentler. Silverman (1984) computed the
equivalent kernel quoted in [Section 43.4](04-smoothing-splines.html), and
Silverman (1985) surveys the approach; the \( O(n) \) computation of the diagonal
of \( \bS_\lambda \) is Hutchinson and de Hoog (1985).

**Choosing the smoothing parameter.**  Generalized cross-validation is Craven
and Wahba (1979), and the leave-one-out shortcut for quadratic penalties
is @thm-sel-loocv. The corrected AIC @eq-smo-aicc is Hurvich, Simonoff and Tsai
(1998), who make the case that its pole is a feature. Opsomer, Wang and Yang
(2001) review what correlated errors do, which is what @exm-smo-correlated
illustrates: no automatic rule separates a smooth mean from smooth noise without
an assumption. Reiss and Ogden (2009) compare REML and GCV. Nothing here
addresses inference *after* the selection except by saying so;
@prp-sel-selection-bias is the warning, and the remedies are those of
[Chapter 23](../ch23-resampling-inference/index.html) and
[Chapter 29](../ch29-model-selection/index.html).

**Mixed models and Bayes.**  The equivalence between a smoothing penalty and a
random effect is old — implicit in Whittaker's graduation and explicit in
Wahba (1978) — and was made into a working method by Ruppert, Wand and Carroll
(2003), whose chapter 4 contains the transformation of @thm-smo-mixed(a) in
essentially the form given here. Wahba (1983) introduced the Bayesian band and Nychka (1988)
analysed its coverage, coining the across-the-function reading that
@prp-smo-bayes-bands(c) states; that part is quoted, not proved, and the
simulation in @exm-smo-coverage is the evidence offered for it. Sun and Loader
(1994) give genuinely simultaneous bands. The adaptive methods named at the end
of [Section 43.6](06-mixed-model-and-bayes.html) are Denison, Mallick and Smith
(1998), Donoho and Johnstone (1994) and Ramsay (1988).

**Coverage.**  This chapter follows the ground of Fahrmeir, Kneib, Lang and Marx
(2021, §8.1 and §8.3); their §8.2, on bivariate and spatial smoothing, belongs to
[Chapter 44](../ch44-additive-models/index.html), as does Hastie and Tibshirani
(1990), the book that made smoothers part of the regression curriculum.

## References

- Cleveland, William S. (1979). Robust Locally Weighted Regression and Smoothing Scatterplots. *Journal of the American Statistical Association* 74(368), 829–836.
- Cleveland, William S. and Devlin, Susan J. (1988). Locally Weighted Regression: An Approach to Regression Analysis by Local Fitting. *Journal of the American Statistical Association* 83(403), 596–610.
- Cox, Maurice G. (1972). The Numerical Evaluation of B-Splines. *Journal of the Institute of Mathematics and Its Applications* 10(2), 134–149.
- Craven, Peter and Wahba, Grace (1979). Smoothing Noisy Data with Spline Functions. *Numerische Mathematik* 31(4), 377–403.
- de Boor, Carl (1972). On Calculating with B-Splines. *Journal of Approximation Theory* 6(1), 50–62.
- de Boor, Carl (2001). *A Practical Guide to Splines*. Revised edition. New York: Springer.
- Demmler, A. and Reinsch, C. (1975). Oscillation Matrices with Spline Smoothing. *Numerische Mathematik* 24(5), 375–382.
- Denison, David G. T., Mallick, Bani K. and Smith, Adrian F. M. (1998). Automatic Bayesian Curve Fitting. *Journal of the Royal Statistical Society, Series B* 60(2), 333–350.
- Donoho, David L. and Johnstone, Iain M. (1994). Ideal Spatial Adaptation by Wavelet Shrinkage. *Biometrika* 81(3), 425–455.
- Eilers, Paul H. C. and Marx, Brian D. (1996). Flexible Smoothing with B-splines and Penalties. *Statistical Science* 11(2), 89–121.
- Epanechnikov, V. A. (1969). Non-parametric Estimation of a Multivariate Probability Density. *Theory of Probability and Its Applications* 14(1), 153–158.
- Eubank, Randall L. (1999). *Nonparametric Regression and Spline Smoothing*. 2nd edition. New York: Marcel Dekker.
- Fahrmeir, Ludwig, Kneib, Thomas, Lang, Stefan and Marx, Brian D. (2021). *Regression: Models, Methods and Applications*. 2nd edition. Berlin: Springer.
- Fan, Jianqing and Gijbels, Irène (1996). *Local Polynomial Modelling and Its Applications*. London: Chapman and Hall.
- Green, Peter J. and Silverman, Bernard W. (1994). *Nonparametric Regression and Generalized Linear Models: A Roughness Penalty Approach*. London: Chapman and Hall.
- Hastie, Trevor and Loader, Clive (1993). Local Regression: Automatic Kernel Carpentry. *Statistical Science* 8(2), 120–129.
- Hastie, Trevor J. and Tibshirani, Robert J. (1990). *Generalized Additive Models*. London: Chapman and Hall.
- Hurvich, Clifford M., Simonoff, Jeffrey S. and Tsai, Chih-Ling (1998). Smoothing Parameter Selection in Nonparametric Regression Using an Improved Akaike Information Criterion. *Journal of the Royal Statistical Society, Series B* 60(2), 271–293.
- Hutchinson, M. F. and de Hoog, F. R. (1985). Smoothing Noisy Data with Spline Functions. *Numerische Mathematik* 47(1), 99–106.
- Kimeldorf, George and Wahba, Grace (1971). Some Results on Tchebycheffian Spline Functions. *Journal of Mathematical Analysis and Applications* 33(1), 82–95.
- Nadaraya, E. A. (1964). On Estimating Regression. *Theory of Probability and Its Applications* 9(1), 141–142.
- Nychka, Douglas (1988). Bayesian Confidence Intervals for Smoothing Splines. *Journal of the American Statistical Association* 83(404), 1134–1143.
- Opsomer, Jean, Wang, Yuedong and Yang, Yuhong (2001). Nonparametric Regression with Correlated Errors. *Statistical Science* 16(2), 134–153.
- O'Sullivan, Finbarr (1986). A Statistical Perspective on Ill-posed Inverse Problems. *Statistical Science* 1(4), 502–518.
- Ramsay, James O. (1988). Monotone Regression Splines in Action. *Statistical Science* 3(4), 425–441.
- Reinsch, Christian H. (1967). Smoothing by Spline Functions. *Numerische Mathematik* 10(3), 177–183.
- Reiss, Philip T. and Ogden, R. Todd (2009). Smoothing Parameter Selection for a Class of Semiparametric Linear Models. *Journal of the Royal Statistical Society, Series B* 71(2), 505–523.
- Ruppert, David (2002). Selecting the Number of Knots for Penalized Splines. *Journal of Computational and Graphical Statistics* 11(4), 735–757.
- Ruppert, David, Sheather, Simon J. and Wand, M. P. (1995). An Effective Bandwidth Selector for Local Least Squares Regression. *Journal of the American Statistical Association* 90(432), 1257–1270.
- Ruppert, David, Wand, M. P. and Carroll, Raymond J. (2003). *Semiparametric Regression*. Cambridge: Cambridge University Press.
- Schoenberg, I. J. (1946). Contributions to the Problem of Approximation of Equidistant Data by Analytic Functions. *Quarterly of Applied Mathematics* 4, 45–99 and 112–141.
- Silverman, Bernard W. (1984). Spline Smoothing: The Equivalent Variable Kernel Method. *The Annals of Statistics* 12(3), 898–916.
- Silverman, Bernard W. (1985). Some Aspects of the Spline Smoothing Approach to Non-parametric Regression Curve Fitting. *Journal of the Royal Statistical Society, Series B* 47(1), 1–52.
- Stone, Charles J. (1980). Optimal Rates of Convergence for Nonparametric Estimators. *The Annals of Statistics* 8(6), 1348–1360.
- Stone, Charles J. (1982). Optimal Global Rates of Convergence for Nonparametric Regression. *The Annals of Statistics* 10(4), 1040–1053.
- Sun, Jiayang and Loader, Clive R. (1994). Simultaneous Confidence Bands for Linear Regression and Smoothing. *The Annals of Statistics* 22(3), 1328–1345.
- Wahba, Grace (1978). Improper Priors, Spline Smoothing and the Problem of Guarding Against Model Errors in Regression. *Journal of the Royal Statistical Society, Series B* 40(3), 364–372.
- Wahba, Grace (1983). Bayesian "Confidence Intervals" for the Cross-validated Smoothing Spline. *Journal of the Royal Statistical Society, Series B* 45(1), 133–150.
- Wahba, Grace (1990). *Spline Models for Observational Data*. Philadelphia: SIAM.
- Wand, M. P. and Jones, M. C. (1995). *Kernel Smoothing*. London: Chapman and Hall.
- Wand, M. P. and Ormerod, J. T. (2008). On Semiparametric Regression with O'Sullivan Penalized Splines. *Australian and New Zealand Journal of Statistics* 50(2), 179–198.
- Watson, Geoffrey S. (1964). Smooth Regression Analysis. *Sankhyā: The Indian Journal of Statistics, Series A* 26(4), 359–372.
- Whittaker, Edmund T. (1923). On a New Method of Graduation. *Proceedings of the Edinburgh Mathematical Society* 41, 63–75.
- Wood, Simon N. (2017). *Generalized Additive Models: An Introduction with R*. 2nd edition. Boca Raton: Chapman and Hall/CRC.
