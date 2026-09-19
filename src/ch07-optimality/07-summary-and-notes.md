# Summary and Notes

## Summary

::: {.idea}

1. A linear estimator \( \mathbf{a}\T\Y \) is unbiased for \( \blambda\T\bbeta \) iff \( \X\T\mathbf{a}=\blambda \). Such estimators exist
   iff \( \blambda\in\C(\X\T) \), the same condition under which least squares estimates \( \blambda\T\bbeta \) uniquely. Their
   coefficient vectors form the flat \( \M\boldsymbol{\uprho}+\C(\X)\perpc \) (@prp-opt-lue and @prp-opt-lue-set).

2. **Gauss–Markov.** Under \( \E(\be)=\bzero \) and \( \Cov(\be)=\sigma^2\I \), in any rank, \( \blambda\T\hbeta \) is the unique BLUE of
   every estimable \( \blambda\T\bbeta \). For any LUE, \( \Var(\mathbf{a}\T\Y)=\Var(\blambda\T\hbeta)+\sigma^2\norm{(\I-\M)\mathbf{a}}^2 \). In
   matrix form, \( \Cov(\A\Y)-\Cov(\bLambda\T\hbeta)=\sigma^2\A(\I-\M)\A\T\succeq\mathbf{0} \) (@thm-opt-gauss-markov and
   @cor-opt-gm-consequences). With \( \Cov(\be)=\sigma^2\V \), generalized least squares is the BLUE (@cor-opt-aitken).

3. The theorem does not cover nonlinear estimators (the midrange beats the mean under
   uniform errors, @exm-opt-midrange) or biased ones (ridge, @exr-opt-ridge). Orthogonal columns minimize coefficient variances
   (@prp-opt-orthogonal-design).

4. Under normality, the MLE of \( \bbeta \) is least squares and the MLE of \( \sigma^2 \) is \( \text{SSE}/n \), biased by the
   factor \( (n-r)/n \) and inconsistent when \( r \) grows with \( n \) (@thm-opt-mle and @exm-opt-neyman-scott). The profile
   likelihood of an estimable function is a decreasing function of \( t^2 \) (@prp-opt-profile-t).

5. \( (\M\Y,\text{SSE}) \) is complete and sufficient. By Lehmann–Scheffé, \( \blambda\T\hbeta \) and \( s^2 \) are the UMVUEs,
   with no linearity restriction (@prp-opt-complete-sufficient and @thm-opt-umvue). Least squares attains the
   Cramér–Rao bound for \( \bbeta \), and \( s^2 \) has efficiency \( (n-p)/n \) (@thm-opt-cramer-rao).

6. Under normality, \( \bLambda\T\hbeta\sim\Normal(\bLambda\T\bbeta,\sigma^2\bLambda\T(\X\T\X)\ginv\bLambda) \),
   \( \text{SSE}/\sigma^2\sim\chi^2(n-r) \), the two are independent, and standardized estimable functions are \( t(n-r) \)
   (@thm-opt-sampling, @cor-opt-t and @cor-opt-quadratic). Without normality the variance of SSE depends on the
   kurtosis (@prp-opt-var-sse).

7. The normal-inverse-gamma prior is conjugate. The posterior mean is a precision-weighted combination of the prior mean and
   \( \hbeta \), and the flat prior \( 1/\sigma^2 \) reproduces least squares and the \( t \) intervals
   (@thm-opt-bayes-conjugate, @prp-opt-shrinkage and @cor-opt-flat-prior).

:::

## Notes and sources

**Gauss and Markov.** Gauss (1823) proved that least squares gives the linear unbiased combination of
observations with the smallest variance, without assuming normal errors. It was his second justification of the
method, after the normal-theory argument of 1809. Markov's lectures on probability (1900 and later editions) restated
the result, and Neyman (1934) attached Markov's name to it. The history is traced by Plackett (1949, 1972) and
Stigler (1986). The projection proof of [Section 7.2](02-gauss-markov.html) splits the coefficient vector of a
competing estimator by Pythagoras, so that its excess variance is \( \sigma^2\norm{(\I-\M)\mathbf{a}}^2 \). This is the
proof of Seber and Lee (2003, section 3.2, Theorem 3.2), stated there for linear functions of the mean vector; we
apply it to estimable functions in the rank-deficient case. Christensen (2020, section 2.3) argues instead through
the covariance with \( \boldsymbol{\uprho}\T\M\Y \). Rencher and Schaalje (2008,
section 7.3) give the matrix form for full-rank models, and Agresti (2015, section 2.7) discusses its scope and its
extension to generalized least squares. The generalized least squares version (@cor-opt-aitken) is due to
Aitken (1935). The characterization of models in which ordinary least squares remains best (@exr-opt-zyskind) is
due to Zyskind (1967) and, in coordinate-free form, to Kruskal (1968).

**What Gauss–Markov does not cover.** That shrinkage can improve on unbiased estimation in three or more
dimensions was discovered by Stein (1956), and James and Stein (1961) gave the explicit estimator. Hoerl and
Kennard (1970) introduced ridge regression. The existence of a better ridge estimator (@exr-opt-ridge) is
their theorem. Wald (1940) proposed the group-means slope estimator of @exm-opt-three-slopes for regressions
with errors in the regressor. The optimality of orthogonal weighing designs goes back to Hotelling (1944).

**Likelihood.** The inconsistency of the maximum likelihood estimate of a variance in the presence of many
nuisance means is the example of Neyman and Scott (1948). Restricted maximum likelihood (@exr-opt-reml) was
developed by Patterson and Thompson (1971) for variance components. Christensen (2020, section 2.4), Seber and Lee
(2003, section 3.5) and Rencher and Schaalje (2008, section 7.6) derive the normal-theory MLEs. The maximization
through \( \log u\le u-1 \) used in @thm-opt-mle follows Seber and Lee (2003, section 3.5).

**Minimum variance.** The Rao–Blackwell theorem is due to Rao (1945) and Blackwell (1947). Completeness and
the uniqueness theorem are due to Lehmann and Scheffé (1950). The two results we quoted, the factorization
criterion and the completeness of full-rank exponential families, are proved in Lehmann and Casella (1998)
and Lehmann and Romano (2005). The same variance bound was found independently by Rao (1945) and
Cramér (1946). Christensen (2020, section 2.5) outlines the minimum variance unbiased argument for the linear model
by reparameterizing to full rank. Seber and Lee (2003, section 3.5) compare the Cramér–Rao bound with the variance of
\( s^2 \). The optimality of \( s^2 \) among nonnegative quadratic unbiased estimators (@exr-opt-quadratic-sigma2) is
due to Atiqullah (1962). That the divisor \( k+2 \) minimizes the mean squared error among multiples of SSE (@prp-opt-divisor)
is due to Theil and Schweitzer (1961). The reflection argument of @exm-opt-reflection is our own illustration.

**Sampling distributions.** The distribution theory of [Section 7.5](05-sampling-distributions.html) follows
Christensen (2020, section 2.6), Seber and Lee (2003, section 3.4) and Rencher and Schaalje (2008, section 7.6.3),
with the rank-deficient forms throughout. The effect of kurtosis on inference about variances, and the relative
robustness of inference about means, were emphasized by Box (1953).

**Bayesian thread.** The conjugate normal-inverse-gamma analysis goes back to Raiffa and Schlaifer (1961) and
is developed for regression by Zellner (1971) and Box and Tiao (1973). O'Hagan and Forster (2004) give a modern
account. Zellner (1986) introduced the \( g \)-prior. Seber and Lee (2003, section 3.12), Christensen (2020, section 2.10) and
Rencher and Schaalje (2008, chapter 11) treat Bayesian estimation in the linear model, and Fahrmeir, Kneib,
Lang and Marx (2021, section 4.4) present the conjugate analysis as the starting point for richer priors. For an
applied account of the Gauss–Markov conditions and what follows from them, see Sen and Srivastava (1990, section 2.5).

**Data.** The stack loss data are from Brownlee (1965) and are distributed with statsmodels.

## References

- Agresti, Alan (2015). *Foundations of Linear and Generalized Linear Models*. Hoboken, NJ: Wiley.
- Aitken, Alexander C. (1935). On Least Squares and Linear Combination of Observations. *Proceedings of the Royal Society of Edinburgh* 55, 42–48.
- Atiqullah, M. (1962). The Estimation of Residual Variance in Quadratically Balanced Least-Squares Problems and the Robustness of the F-Test. *Biometrika* 49(1/2), 83–91.
- Blackwell, David (1947). Conditional Expectation and Unbiased Sequential Estimation. *The Annals of Mathematical Statistics* 18(1), 105–110.
- Box, George E. P. (1953). Non-Normality and Tests on Variances. *Biometrika* 40(3/4), 318–335.
- Box, George E. P. and Tiao, George C. (1973). *Bayesian Inference in Statistical Analysis*. Reading, MA: Addison-Wesley.
- Brownlee, K. A. (1965). *Statistical Theory and Methodology in Science and Engineering*. 2nd edition. New York: Wiley.
- Christensen, Ronald (2020). *Plane Answers to Complex Questions: The Theory of Linear Models*. 5th edition. Cham: Springer.
- Cramér, Harald (1946). *Mathematical Methods of Statistics*. Princeton, NJ: Princeton University Press.
- Fahrmeir, Ludwig, Kneib, Thomas, Lang, Stefan and Marx, Brian D. (2021). *Regression: Models, Methods and Applications*. 2nd edition. Berlin: Springer.
- Gauss, Carl Friedrich (1823). *Theoria Combinationis Observationum Erroribus Minimis Obnoxiae*. Göttingen: Dieterich.
- Hoerl, Arthur E. and Kennard, Robert W. (1970). Ridge Regression: Biased Estimation for Nonorthogonal Problems. *Technometrics* 12(1), 55–67.
- Hotelling, Harold (1944). Some Improvements in Weighing and Other Experimental Techniques. *The Annals of Mathematical Statistics* 15(3), 297–306.
- James, W. and Stein, Charles (1961). Estimation with Quadratic Loss. *Proceedings of the Fourth Berkeley Symposium on Mathematical Statistics and Probability* 1, 361–379.
- Kruskal, William (1968). When Are Gauss–Markov and Least Squares Estimators Identical? A Coordinate-Free Approach. *The Annals of Mathematical Statistics* 39(1), 70–75.
- Lehmann, E. L. and Casella, George (1998). *Theory of Point Estimation*. 2nd edition. New York: Springer.
- Lehmann, E. L. and Romano, Joseph P. (2005). *Testing Statistical Hypotheses*. 3rd edition. New York: Springer.
- Lehmann, E. L. and Scheffé, Henry (1950). Completeness, Similar Regions, and Unbiased Estimation. Part I. *Sankhyā* 10(4), 305–340.
- Neyman, Jerzy (1934). On the Two Different Aspects of the Representative Method. *Journal of the Royal Statistical Society* 97(4), 558–625.
- Neyman, J. and Scott, Elizabeth L. (1948). Consistent Estimates Based on Partially Consistent Observations. *Econometrica* 16(1), 1–32.
- O'Hagan, Anthony and Forster, Jonathan (2004). *Kendall's Advanced Theory of Statistics, Volume 2B: Bayesian Inference*. 2nd edition. London: Arnold.
- Patterson, H. D. and Thompson, Robin (1971). Recovery of Inter-Block Information When Block Sizes Are Unequal. *Biometrika* 58(3), 545–554.
- Plackett, R. L. (1949). A Historical Note on the Method of Least Squares. *Biometrika* 36(3/4), 458–460.
- Plackett, R. L. (1972). Studies in the History of Probability and Statistics. XXIX: The Discovery of the Method of Least Squares. *Biometrika* 59(2), 239–251.
- Raiffa, Howard and Schlaifer, Robert (1961). *Applied Statistical Decision Theory*. Boston: Division of Research, Graduate School of Business Administration, Harvard University.
- Rao, C. Radhakrishna (1945). Information and the Accuracy Attainable in the Estimation of Statistical Parameters. *Bulletin of the Calcutta Mathematical Society* 37, 81–91.
- Rencher, Alvin C. and Schaalje, G. Bruce (2008). *Linear Models in Statistics*. 2nd edition. Hoboken, NJ: Wiley.
- Seber, George A. F. and Lee, Alan J. (2003). *Linear Regression Analysis*. 2nd edition. Hoboken, NJ: Wiley.
- Sen, Ashish and Srivastava, Muni (1990). *Regression Analysis: Theory, Methods, and Applications*. New York: Springer.
- Stein, Charles (1956). Inadmissibility of the Usual Estimator for the Mean of a Multivariate Normal Distribution. *Proceedings of the Third Berkeley Symposium on Mathematical Statistics and Probability* 1, 197–206.
- Stigler, Stephen M. (1986). *The History of Statistics: The Measurement of Uncertainty before 1900*. Cambridge, MA: Harvard University Press.
- Theil, H. and Schweitzer, A. (1961). The Best Quadratic Estimator of the Residual Variance in Regression Analysis. *Statistica Neerlandica* 15(1), 19–23.
- Wald, Abraham (1940). The Fitting of Straight Lines if Both Variables Are Subject to Error. *The Annals of Mathematical Statistics* 11(3), 284–300.
- Zellner, Arnold (1971). *An Introduction to Bayesian Inference in Econometrics*. New York: Wiley.
- Zellner, Arnold (1986). On Assessing Prior Distributions and Bayesian Regression Analysis with g-Prior Distributions. In P. K. Goel and A. Zellner (eds.), *Bayesian Inference and Decision Techniques: Essays in Honor of Bruno de Finetti*, 233–243. Amsterdam: North-Holland.
- Zyskind, George (1967). On Canonical Forms, Non-Negative Covariance Matrices and Best and Simple Least Squares Linear Estimators in Linear Models. *The Annals of Mathematical Statistics* 38(4), 1092–1109.
