# Summary and Notes

## Summary

::: {.idea}

1. Expectation of a random matrix is taken entry by entry, and it commutes with constant
           linear maps, transposes and traces (@prp-rv-expectation-linear).

2. \( \Cov(\Y)=\E[(\Y-\bmu)(\Y-\bmu)\T]=\E(\Y\Y\T)-\bmu\bmu\T \), and
           \( \Cov(\bU,\V) \) is the rectangular cross-covariance, with
           \( \Cov(\V,\bU)=\Cov(\bU,\V)\T \) (@def-rv-mean-cov and @prp-rv-cov-basic).
           Uncorrelated does not imply independent (@exm-rv-circle).

3. Affine maps: \( \E(\A\Y+\bb)=\A\bmu+\bb \), \( \Cov(\A\Y+\bb)=\A\bSigma\A\T \),
           \( \Cov(\A\bU,\B\V)=\A\Cov(\bU,\V)\B\T \) (@thm-rv-linear). In particular
           \( \Var(\mathbf{a}\T\Y)=\mathbf{a}\T\bSigma\mathbf{a} \).

4. Covariance matrices are exactly the nonnegative definite matrices. \( \bSigma\mathbf{a}=\bzero \)
           iff \( \mathbf{a}\T\Y \) is constant with probability one, and \( \Y \) lives in the flat
           \( \bmu+\C(\bSigma) \) (@thm-rv-cov-nnd).

5. Whitening by \( \W \) with \( \W\bSigma\W\T=\I \) is unique up to an orthogonal factor.
           Every choice gives the same Mahalanobis distance \( (\Y-\bmu)\T\bSigma^{-1}(\Y-\bmu) \)
           (@prp-rv-whitening).

6. \( \E(\Y\T\A\Y)=\tr(\A\bSigma)+\bmu\T\A\bmu \) under second moments only
           (@thm-rv-quadform-mean). Hence \( \E(\text{RSS})=\sigma^2(n-r)+\norm{(\I-\M)\boldsymbol{\uptheta}}^2 \)
           (@exm-rv-rss-bias). The variance of a quadratic form needs third and fourth moments
           (@thm-rv-quadform-variance).

7. Correlation matrices are the nonnegative definite matrices with unit diagonal. They are
           invariant to rescaling, and their entries must be jointly consistent
           (@prp-rv-correlation and @prp-rv-three-correlations).

8. For partitioned vectors, \( \C(\bSigma_{XY})\subseteq\C(\bSigma_{XX}) \). The best linear
           predictor \( \bmu_Y+\bSigma_{YX}\bSigma_{XX}\ginv(\X-\bmu_X) \) leaves an error uncorrelated with
           \( \X \) with covariance matrix given by the Schur complement \( \bSigma_{YY\cdot X} \)
           (@thm-rv-partitioned and @thm-rv-blp).

9. Moment generating functions determine distributions and factor exactly under
           independence (@thm-rv-mgf). Pairwise independence is weaker than mutual
           independence (@exm-rv-pairwise).

10. The sample covariance matrix is unbiased under uncorrelated rows and has rank at most
            \( n-1 \). A tiny eigenvalue of the sample correlation matrix flags a nearly constant
            combination (@prp-rv-sample-moments and @exm-rv-longley).

:::

## Notes and sources

**The operator calculus.**  Treating expectation and covariance as operators on random
vectors and matrices is standard in linear-model texts. Rencher and Schaalje (2008, chapter 3)
develops mean vectors, covariance and correlation matrices, partitioned vectors and linear
functions in detail and at a gentle pace. Seber and Lee (2003, §§1.4–1.6) give a
compact operator treatment, first and second moments of quadratic forms, and moment generating
functions with their use in proving independence. This chapter covers the same ground.
Rencher and Schaalje already introduce the standardized distance and the transformation by
\( \bSigma^{-1/2} \). The chapter adds the characterization of singular covariance matrices through
their column space (@thm-rv-cov-nnd), the non-uniqueness of whitening and the comparison
of symmetric and Cholesky whitening, the best linear predictor with generalized inverses, and
the sample covariance matrix. Christensen (2020, chapter 1) states the
covariance rules briefly in the notation of linear models. For the multivariate setting,
Anderson (2003) remains the standard reference.

**Mahalanobis distance and whitening.**  The distance
@eq-rv-mahalanobis is due to Mahalanobis (1936), who introduced it to
compare populations described by several correlated measurements. Kessy et al. (2018) compare
whitening transformations. They show that the symmetric choice \( \bSigma^{-1/2} \) is the one closest to the
original data in mean squared error (@exr-rv-whitening-closest), and set it against
Cholesky and principal-component whitening.

**Quadratic forms.**  The trace argument for @thm-rv-quadform-mean is folklore.
@thm-rv-quadform-variance is stated without proof by Atiqullah (1962),
who used it to study how nonnormal errors affect the residual variance estimate and the \( F \)-test.
The kurtosis term depends on \( \A \) only through its diagonal. For the residual sum of squares of a
correct model it equals \( (\mu_4-3\sigma^4)\sum_i(1-h_{ii})^2 \), which takes its simplest form,
\( (\mu_4-3\sigma^4)(n-r)^2/n \), when all leverages \( h_{ii} \) are equal. These are the quadratically balanced designs of
Atiqullah's title.
Seber and Lee (2003, theorem 1.6) give a proof by direct expansion. The proof here organizes the
same expansion by splitting the form into quadratic, linear and constant parts.
Mathai and Provost (1992) is a book-length treatment of moments and distributions of
quadratic forms, including general covariance structures.

**Correlation, partial covariance and prediction.**  The feasible range in
@prp-rv-three-correlations is a special case of the fact that a symmetric matrix is
nonnegative definite iff a Schur complement is, when the leading block is positive definite.
Horn and Johnson (2013) give the general Schur-complement and Loewner-order results used in @thm-rv-blp.
Stating the best linear predictor with a generalized inverse, so that nothing
needs \( \bSigma_{XX} \) to be nonsingular, is in the spirit of Rao (1973), who uses
generalized inverses throughout.

**Probability background.**  Uniqueness of moment generating functions and
characteristic functions, and the device of @prp-rv-cramer-wold, are proved in
Billingsley (1995). The Longley data come from
Longley (1967), who used them to test the numerical accuracy of regression
programs. Their near-collinearity, visible in @exm-rv-longley, is the reason.

## References

- Anderson, Theodore W. (2003). *An Introduction to Multivariate Statistical Analysis*. 3rd edition. Hoboken, NJ: Wiley.
- Atiqullah, M. (1962). The Estimation of Residual Variance in Quadratically Balanced Least-Squares Problems and the Robustness of the F-Test. *Biometrika* 49(1/2), 83–91.
- Billingsley, Patrick (1995). *Probability and Measure*. 3rd edition. New York: Wiley.
- Christensen, Ronald (2020). *Plane Answers to Complex Questions: The Theory of Linear Models*. 5th edition. Cham: Springer.
- Horn, Roger A. and Johnson, Charles R. (2013). *Matrix Analysis*. 2nd edition. Cambridge: Cambridge University Press.
- Kessy, Agnan, Lewin, Alex and Strimmer, Korbinian (2018). Optimal Whitening and Decorrelation. *The American Statistician* 72(4), 309–314.
- Longley, James W. (1967). An Appraisal of Least Squares Programs for the Electronic Computer from the Point of View of the User. *Journal of the American Statistical Association* 62(319), 819–841.
- Mahalanobis, Prasanta Chandra (1936). On the Generalized Distance in Statistics. *Proceedings of the National Institute of Sciences of India* 2(1), 49–55.
- Mathai, Arak M. and Provost, Serge B. (1992). *Quadratic Forms in Random Variables: Theory and Applications*. New York: Marcel Dekker.
- Rao, C. Radhakrishna (1973). *Linear Statistical Inference and Its Applications*. 2nd edition. New York: Wiley.
- Rencher, Alvin C. and Schaalje, G. Bruce (2008). *Linear Models in Statistics*. 2nd edition. Hoboken, NJ: Wiley.
- Seber, George A. F. and Lee, Alan J. (2003). *Linear Regression Analysis*. 2nd edition. Hoboken, NJ: Wiley.
