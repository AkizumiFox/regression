# Summary and Notes

## Summary

::: {.idea}

1. The bootstrap estimates the distribution of \( \hat{\theta}-\theta \) by that of \( \hat{\theta}^*-\hat{\theta} \) under an estimated
   mechanism: resampled residuals with the design fixed, or resampled cases (@def-bs-bootstrap).

2. The residual bootstrap has covariance \( \hat{\sigma}_*^2(\X\T\X)^{-1} \), biased by \( (n-p)/n \) (@prp-bs-residual-moments), and is
   consistent for contrasts under exchangeable errors when no case dominates (@thm-bs-residual); it fails under
   heteroscedasticity.

3. The case bootstrap, linearized, has the HC0 sandwich covariance (@prp-bs-case-linear). The wild bootstrap reproduces the
   sandwich exactly (@prp-bs-wild-moments) and is consistent under heteroscedasticity (@thm-bs-wild-consistency).

4. Percentile, basic and studentized intervals are first-order correct (@cor-bs-interval-consistency); studentizing reproduces
   exact pivots exactly (@prp-bs-parametric) and, by Hall's Edgeworth analysis (not proved here), gives higher-order accuracy.

5. Permutation tests are exact for any statistic under invariance (@thm-bs-permutation): in simple regression with exchangeable
   errors, and with nuisance regressors only for permutations that fix them (@prp-bs-invariant-nuisance). Heteroscedasticity turns
   them into tests of independence unless the statistic is studentized.

6. The Freedman–Lane scheme is exact when an exact test exists, ignores the nuisance coefficients, and its permuted
   statistics are asymptotically those of an exact oracle (@prp-bs-freedman-lane(d)); its asymptotic size is due to Anderson and
   Robinson (2001). Kennedy's scheme is more liberal (@prp-bs-kennedy).

7. Resampling fails for extremes (@prp-bs-maximum), heavy tails, dominant cases, large \( p/n \), dependence and, unless studentized,
   small samples ([Section 23.6](06-failures.html)).

:::

## Notes and sources

**The bootstrap.** Efron (1979) introduced the bootstrap. Freedman (1981) set out residual and case resampling for regression
and proved both consistent; @thm-bs-residual is a version for linear contrasts. Bickel and Freedman (1981) developed the general
theory, and Wu (1986) showed that the residual bootstrap fails under heteroscedasticity. Efron and Tibshirani (1993) and Davison
and Hinkley (1997) are the standard books; the latter uses the leverage-adjusted residuals adopted here. Sen and Srivastava
(1990, section 5.4) describe the residual bootstrap test of a linear hypothesis.

**The wild bootstrap.** The idea is in Wu (1986); Liu (1988) and Mammen (1993) developed its theory, and Mammen introduced the
two-point law. Davidson and Flachaire (2008) recommend Rademacher weights with restricted residuals for tests.

**Intervals.** Hall (1988) compared bootstrap intervals by Edgeworth expansions, and Hall (1992) is the full account, including
the iterated bootstrap of @exr-bs-iterated. Efron (1987) introduced BCa; DiCiccio and Efron (1996) review the methods. Prediction intervals for a new response, which add a resampled
error to the resampled estimate in place of the normal-theory interval of @thm-ci-prediction-interval, are treated in Davison and
Hinkley (1997, chapter 6).

**Permutation tests.** Fisher (1935) and Pitman (1937) introduced permutation tests. The proof of @thm-bs-permutation follows
Lehmann and Romano (2005, chapter 15); on the Monte Carlo count see Phipson and Smyth (2010). DiCiccio and Romano (2017) proved
the validity of studentized permutation tests under dependence without correlation.

**Partial tests.** The schemes are due to Freedman and Lane (1983), ter Braak (1992), Kennedy (1995), Draper and Stoneman (1966)
and Manly (1997). Anderson and Legendre (1999) compared them by simulation, Anderson and Robinson (2001) analytically, and
Winkler and coauthors (2014) review them for the general linear model.

**Failures.** Bickel, Götze and van Zwet (1997) analyse the \( m \)-out-of-\( n \) bootstrap, Athreya (1987) the infinite-variance
mean, and El Karoui and Purdom (2018) large \( p/n \). Künsch (1989) introduced the moving-block bootstrap; Lahiri (2003) treats
dependent data at length.

## References

- Anderson, Marti J. and Legendre, Pierre (1999). An Empirical Comparison of Permutation Methods for Tests of Partial Regression Coefficients in a Linear Model. *Journal of Statistical Computation and Simulation* 62(3), 271–303.
- Anderson, Marti J. and Robinson, John (2001). Permutation Tests for Linear Models. *Australian and New Zealand Journal of Statistics* 43(1), 75–88.
- Athreya, Krishna B. (1987). Bootstrap of the Mean in the Infinite Variance Case. *The Annals of Statistics* 15(2), 724–731.
- Bickel, Peter J. and Freedman, David A. (1981). Some Asymptotic Theory for the Bootstrap. *The Annals of Statistics* 9(6), 1196–1217.
- Bickel, Peter J., Götze, Friedrich and van Zwet, Willem R. (1997). Resampling Fewer than n Observations: Gains, Losses, and Remedies for Losses. *Statistica Sinica* 7(1), 1–31.
- Davidson, Russell and Flachaire, Emmanuel (2008). The Wild Bootstrap, Tamed at Last. *Journal of Econometrics* 146(1), 162–169.
- Davison, Anthony C. and Hinkley, David V. (1997). *Bootstrap Methods and Their Application*. Cambridge: Cambridge University Press.
- DiCiccio, Cyrus J. and Romano, Joseph P. (2017). Robust Permutation Tests for Correlation and Regression Coefficients. *Journal of the American Statistical Association* 112(519), 1211–1220.
- DiCiccio, Thomas J. and Efron, Bradley (1996). Bootstrap Confidence Intervals. *Statistical Science* 11(3), 189–228.
- Draper, Norman R. and Stoneman, David M. (1966). Testing for the Inclusion of Variables in Linear Regression by a Randomisation Technique. *Technometrics* 8(4), 695–699.
- Efron, Bradley (1979). Bootstrap Methods: Another Look at the Jackknife. *The Annals of Statistics* 7(1), 1–26.
- Efron, Bradley (1987). Better Bootstrap Confidence Intervals. *Journal of the American Statistical Association* 82(397), 171–185.
- Efron, Bradley and Tibshirani, Robert (1993). *An Introduction to the Bootstrap*. New York: Chapman and Hall.
- El Karoui, Noureddine and Purdom, Elizabeth (2018). Can We Trust the Bootstrap in High-Dimensions? The Case of Linear Models. *Journal of Machine Learning Research* 19, 1–66.
- Fisher, Ronald A. (1935). *The Design of Experiments*. Edinburgh: Oliver and Boyd.
- Freedman, David A. (1981). Bootstrapping Regression Models. *The Annals of Statistics* 9(6), 1218–1228.
- Freedman, David and Lane, David (1983). A Nonstochastic Interpretation of Reported Significance Levels. *Journal of Business and Economic Statistics* 1(4), 292–298.
- Hall, Peter (1988). Theoretical Comparison of Bootstrap Confidence Intervals. *The Annals of Statistics* 16(3), 927–953.
- Hall, Peter (1992). *The Bootstrap and Edgeworth Expansion*. New York: Springer.
- Kennedy, Peter E. (1995). Randomization Tests in Econometrics. *Journal of Business and Economic Statistics* 13(1), 85–94.
- Künsch, Hans R. (1989). The Jackknife and the Bootstrap for General Stationary Observations. *The Annals of Statistics* 17(3), 1217–1241.
- Lahiri, Soumendra N. (2003). *Resampling Methods for Dependent Data*. New York: Springer.
- Lehmann, Erich L. and Romano, Joseph P. (2005). *Testing Statistical Hypotheses*. 3rd edition. New York: Springer.
- Liu, Regina Y. (1988). Bootstrap Procedures under Some Non-I.I.D. Models. *The Annals of Statistics* 16(4), 1696–1708.
- Mammen, Enno (1993). Bootstrap and Wild Bootstrap for High Dimensional Linear Models. *The Annals of Statistics* 21(1), 255–285.
- Manly, Bryan F. J. (1997). *Randomization, Bootstrap and Monte Carlo Methods in Biology*. 2nd edition. London: Chapman and Hall.
- Phipson, Belinda and Smyth, Gordon K. (2010). Permutation P-values Should Never Be Zero: Calculating Exact P-values When Permutations Are Randomly Drawn. *Statistical Applications in Genetics and Molecular Biology* 9(1), Article 39.
- Pitman, Edwin J. G. (1937). Significance Tests Which May Be Applied to Samples from Any Populations. *Supplement to the Journal of the Royal Statistical Society* 4(1), 119–130.
- Sen, Ashish and Srivastava, Muni (1990). *Regression Analysis: Theory, Methods, and Applications*. New York: Springer.
- ter Braak, Cajo J. F. (1992). Permutation Versus Bootstrap Significance Tests in Multiple Regression and ANOVA. In K.-H. Jöckel, G. Rothe and W. Sendler (eds.), *Bootstrapping and Related Techniques*. Berlin: Springer.
- van der Vaart, Aad W. (1998). *Asymptotic Statistics*. Cambridge: Cambridge University Press.
- Winkler, Anderson M., Ridgway, Gerard R., Webster, Matthew A., Smith, Stephen M. and Nichols, Thomas E. (2014). Permutation Inference for the General Linear Model. *NeuroImage* 92, 381–397.
- Wu, C. F. Jeff (1986). Jackknife, Bootstrap and Other Resampling Methods in Regression Analysis. *The Annals of Statistics* 14(4), 1261–1295.
