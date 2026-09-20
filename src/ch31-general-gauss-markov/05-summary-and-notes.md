# Summary and Notes

## Summary

::: {.idea}

1. The general Gauss–Markov model @eq-ggm-model keeps \( \E(\Y)=\X\bbeta \) and
   replaces \( \sigma^2\I \) by \( \sigma^2\V \). Only the product is identified, so
   \( \V \) carries the *shape* of the covariance and \( \sigma^2 \) its scale (@exr-ggm-scale).

2. With \( \V \) known and positive definite, generalized least squares (@def-ggm-gls)
   is ordinary least squares on the whitened data @eq-ggm-whitened, and everything
   in Parts II and III transfers: the \( \chi^2 \) law of \( \text{SSE}_{\V} \), the
   unbiased \( s_{\V}^2 \), \( t \) intervals and the \( F \) test of a reduced model,
   all on \( n-r \) degrees of freedom (@thm-ggm-inference). What does not transfer
   is the Euclidean reading of the residual, which is \( \V^{-1} \)-orthogonal to
   \( \C(\X) \); \( R^2 \) and leverage lose their meaning.

3. For structured \( \V \) the whitening has a closed form. Equicorrelated clusters
   give the quasi-demeaning @eq-ggm-quasi, which interpolates between doing
   nothing and fully centring each cluster.

4. A linear unbiased estimator \( \mathbf{c}\T\Y \) of \( \blambda\T\bbeta \) is best exactly
   when \( \V\mathbf{c}\in\C(\X) \), and any two best estimators agree on \( \C(\X:\V) \) (@lem-ggm-blue-criterion).
   This one criterion drives the rest of the chapter.

5. Ordinary least squares is best for *everything* iff \( \C(\V\X)\subseteq\C(\X) \),
   iff \( \M\V=\V\M \), iff \( \V=\X\A\X\T+\W\B\W\T \), iff \( \C(\X) \) is spanned by
   eigenvectors of \( \V \) (@thm-ggm-ols-blue). Balanced layouts with a shared
   shock satisfy it; unbalanced ones do not, and then some contrasts are still
   estimated best and others are not.

6. The efficiency lost by ignoring the covariance is bounded below by the
   Kantorovich ratio and is usually small (@prp-ggm-ols-loss). The error in the
   reported standard error is not bounded at all and is usually large. Fix the
   standard error first.

7. When \( \V \) is singular the errors live in \( \C(\V) \), so the data live in
   \( \C(\X:\V) \) and \( r-d \) linear functions of \( \bbeta \) are known exactly (@thm-ggm-singular).
   The BLUE uses those functions to pin down a particular solution and then does
   generalized least squares inside \( \C(\X)\cap\C(\V) \), or comes in closed form
   from Rao's formula @eq-ggm-rao. The residual degrees of freedom are
   \( \rank(\X:\V)-\rank(\X) \), not \( n-r \).

8. When \( \V=\V(\boldsymbol{\uptheta}) \) must be estimated, the two-step estimator
   @eq-ggm-fgls is asymptotically equivalent to generalized least squares under
   stated regularity and a stochastic equicontinuity condition (@thm-ggm-feasible, @lem-ggm-drift).
   In small samples the efficiency comes back almost at once while the reported
   standard error stays badly optimistic (@exm-ggm-fgls-simulation).
:::

## Notes and sources

**Generalized least squares.** The estimator and its optimality are due to
Aitken (1935), who derived it by the same change of variables used here.
Christensen (2020, §2.7) develops it coordinate-freely for arbitrary rank, the
route [Section 6.9](../ch06-projections/09-inner-products.html) follows; Seber and
Lee (2003, §3.10) and Rencher and Schaalje (2008, §7.8) give the full-rank
version. That only \( \sigma^2\V \) is a parameter is standard but often left
implicit, so "normalizing \( \V \)" is sometimes mistaken for an assumption.

**When least squares is best.** The condition \( \C(\V\X)\subseteq\C(\X) \) has been
found many times: Kruskal (1968) gave the coordinate-free statement that the two
*fits* coincide, proved here as @thm-proj-kruskal; Zyskind (1967) characterized
through eigenvectors the covariance structures for which least squares is best,
which is condition (v); and Rao (1967) gave the structural form
\( \V=\X\A\X\T+\W\B\W\T \) now usually called simple covariance structure.
Watson (1955) first asked how much efficiency is lost when the condition fails;
the sharp scalar bound follows from Kantorovich's inequality, proved as
@thm-dep-efficiency in
[Section 19.2](../ch19-theory-of-departures/02-wrong-covariance.html), and the
determinant bound quoted in [Section 31.2](02-ols-blue.html) was proved
independently by Bloomfield and Watson (1975) and by Knott (1975). Puntanen and
Styan (1989) survey the whole equivalence. Recovery of interblock information in
incomplete block designs, the applied face of the failure of the condition, goes
back to Yates (1940).

**Singular covariance matrices.** The systematic theory is due to
Rao (1971, 1973), whose unified theory of linear estimation produced the formula
@eq-ggm-rao with \( \bT=\V+\X\bS\X\T \) for any nonnegative definite \( \bS \)
with \( \C(\X)\subseteq\C(\bT) \), and to Rao and Mitra (1971) for the
generalized-inverse machinery. Christensen (2020, ch. 10) gives the
vector-space treatment, including consistent estimation and the
degrees-of-freedom count \( \rank(\X:\V)-\rank(\X) \). The proof given here is
different: it splits the observation space orthogonally into \( \C(\V) \) and its
complement, turning the model into an ordinary generalized least squares problem
subject to exact linear constraints, and then verifies the two conditions
of @lem-ggm-blue-criterion. Our @thm-ggm-singular(c) uses the Moore–Penrose
inverse throughout, which keeps it short; with an arbitrary generalized inverse of
\( \bT \) the same statements hold, but \( \bT\bT\ginv\X=\X \) and the invariance of
@eq-ggm-rao then have to be established separately.

**Feasible generalized least squares.** Two-step estimation with an estimated
covariance is old; the econometric literature traces it through Cochrane and
Orcutt (1949) for autoregressive errors and Zellner (1962) for seemingly
unrelated regressions. Amemiya (1985, ch. 6) states and proves the asymptotic
equivalence for several concrete families, and Carroll and Ruppert (1988, ch. 3)
for smooth variance functions; both need a uniformity condition of the kind
labelled (G5) in @thm-ggm-feasible, and our @lem-ggm-drift isolates the
cancellation that makes such a condition plausible. The finite-sample deficiency
of the plug-in covariance is quantified by Kackar and Harville (1984), whose
correction underlies the Kenward–Roger adjustment of
[Chapter 32](../ch32-linear-mixed-models/index.html). The complementary strategy of
leaving the estimator alone and correcting the standard error is
White (1980) and, for clustered data, Liang and Zeger (1986); it appears in this
book as @thm-het-sandwich.

**What comes next.** [Chapter 32](../ch32-linear-mixed-models/index.html) makes
\( \V \) the marginal covariance of a model with random effects, which supplies
both a reason for a particular \( \V(\boldsymbol{\uptheta}) \) and better estimators of
\( \boldsymbol{\uptheta} \) than the moment estimators used here.
[Chapter 33](../ch33-clustered-longitudinal-splitplot/index.html) applies the
whole apparatus to split-plot, repeated-measures and longitudinal designs, where
@thm-ggm-ols-blue explains which classical analyses survive and which do not.

## References

- Aitken, Alexander C. (1935). On Least Squares and Linear Combination of Observations. *Proceedings of the Royal Society of Edinburgh* 55, 42–48.
- Amemiya, Takeshi (1985). *Advanced Econometrics*. Cambridge, MA: Harvard University Press.
- Bloomfield, Peter and Watson, Geoffrey S. (1975). The Inefficiency of Least Squares. *Biometrika* 62(1), 121–128.
- Carroll, Raymond J. and Ruppert, David (1988). *Transformation and Weighting in Regression*. New York: Chapman and Hall.
- Christensen, Ronald (2020). *Plane Answers to Complex Questions: The Theory of Linear Models*. 5th edition. Cham: Springer.
- Cochrane, Donald and Orcutt, Guy H. (1949). Application of Least Squares Regression to Relationships Containing Auto-Correlated Error Terms. *Journal of the American Statistical Association* 44(245), 32–61.
- Kackar, Raghu N. and Harville, David A. (1984). Approximations for Standard Errors of Estimators of Fixed and Random Effects in Mixed Linear Models. *Journal of the American Statistical Association* 79(388), 853–862.
- Knott, Martin (1975). On the Minimum Efficiency of Least Squares. *Biometrika* 62(1), 129–132.
- Kruskal, William (1968). When Are Gauss–Markov and Least Squares Estimators Identical? A Coordinate-Free Approach. *The Annals of Mathematical Statistics* 39(1), 70–75.
- Liang, Kung-Yee and Zeger, Scott L. (1986). Longitudinal Data Analysis Using Generalized Linear Models. *Biometrika* 73(1), 13–22.
- Puntanen, Simo and Styan, George P. H. (1989). The Equality of the Ordinary Least Squares Estimator and the Best Linear Unbiased Estimator. *The American Statistician* 43(3), 153–161.
- Rao, C. Radhakrishna (1967). Least Squares Theory Using an Estimated Dispersion Matrix and Its Application to Measurement of Signals. *Proceedings of the Fifth Berkeley Symposium on Mathematical Statistics and Probability* 1, 355–372.
- Rao, C. Radhakrishna (1971). Unified Theory of Linear Estimation. *Sankhyā, Series A* 33(4), 371–394.
- Rao, C. Radhakrishna (1973). *Linear Statistical Inference and Its Applications*. 2nd edition. New York: Wiley.
- Rao, C. Radhakrishna and Mitra, Sujit Kumar (1971). *Generalized Inverse of Matrices and Its Applications*. New York: Wiley.
- Rencher, Alvin C. and Schaalje, G. Bruce (2008). *Linear Models in Statistics*. 2nd edition. Hoboken, NJ: Wiley.
- Seber, George A. F. and Lee, Alan J. (2003). *Linear Regression Analysis*. 2nd edition. Hoboken, NJ: Wiley.
- Watson, Geoffrey S. (1955). Serial Correlation in Regression Analysis I. *Biometrika* 42(3–4), 327–341.
- White, Halbert (1980). A Heteroskedasticity-Consistent Covariance Matrix Estimator and a Direct Test for Heteroskedasticity. *Econometrica* 48(4), 817–838.
- Yates, Frank (1940). The Recovery of Inter-Block Information in Balanced Incomplete Block Designs. *Annals of Eugenics* 10(4), 317–325.
- Zellner, Arnold (1962). An Efficient Method of Estimating Seemingly Unrelated Regressions and Tests for Aggregation Bias. *Journal of the American Statistical Association* 57(298), 348–368.
- Zyskind, George (1967). On Canonical Forms, Non-Negative Covariance Matrices and Best and Simple Least Squares Linear Estimators in Linear Models. *The Annals of Mathematical Statistics* 38(4), 1092–1109.
