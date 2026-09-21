# Summary and Notes

## Summary

::: {.idea}

1. The pattern of missingness is data. Writing the joint law as a model for the
   complete array times a mechanism @eq-mis-selection separates what the analyst
   wants from how the holes arose, and sorts mechanisms into MCAR, MAR and MNAR
   (@def-mis-mechanisms).

2. Nothing in the observed data distinguishes MAR from MNAR: every choice of the
   conditional distribution of the unrecorded values given the recorded ones gives
   the same observed-data likelihood (@prp-mis-untestable). MCAR against MAR is
   testable; MAR against MNAR is not.

3. A complete-case analysis estimates functionals of the law given completeness. It
   is valid under MCAR, and — the useful case — valid for a correctly specified
   regression whenever the missingness depends on the covariates but not on the
   response (@prp-mis-complete-case). It is not generally valid under MAR.

4. Mean imputation deflates the variance by the fraction recorded
   (@prp-mis-mean-imputation); regression imputation without noise inflates
   correlations; any single imputation understates every standard error;
   available-case correlation matrices need not be nonnegative definite.

5. Under MAR with distinct parameters the mechanism is ignorable: likelihood,
   likelihood ratio tests and Bayesian posteriors may be computed from the
   observed-data likelihood alone (@thm-mis-ignorable). The licence costs a model
   for any variable that can go missing, including covariates. For a monotone
   pattern that likelihood factors into conditional models on nested subsamples
   (@exr-mis-monotone of [Section 41.1](01-mechanisms.html) and
   @exr-mis-normal-monotone of [Section 41.3](03-likelihood.html)).

6. EM maximizes the observed-data likelihood by alternating a conditional
   expectation of the complete-data log-likelihood with a complete-data fit
   (@def-mis-em). Every sweep increases the observed-data likelihood
   (@thm-mis-em), and fixed points are stationary points, not necessarily maxima.

7. For missing responses in a linear model, EM returns the observed-data least
   squares estimate and \( \mathrm{SSE}_o/(n-m) \) (@prp-mis-em-responses),
   deriving the degrees-of-freedom correction that @thm-dsn-missing(e) asserted. For
   a missing covariate the E step must condition on the response
   (@prp-mis-em-covariate).

8. Observed information equals expected complete-data information minus the
   variance of the complete-data score (@prp-mis-louis). That difference sets both
   EM's rate of convergence and the amount by which the complete-data information
   understates the variance; standard errors come from Louis's formula, numerical
   differentiation or the bootstrap, never from the matrix the M step inverts.

9. Multiple imputation draws \( M \) completions from the posterior predictive
   distribution (@def-mis-multiple-imputation), analyses each, and combines by
   \( \bar Q \) and \( T=\bar U+(1+M^{-1})B \) on \( (M-1)/\gamma^2 \) degrees of
   freedom (@thm-mis-rubin). A finite \( M \) inflates the standard error by
   \( (1+\gamma_\infty/M)^{1/2} \). The imputation model must contain everything
   the analysis will use, above all the response.

10. Under MNAR the selection and pattern-mixture factorizations describe the same
    unidentified object (@prp-mis-sensitivity). A delta-adjustment names the
    unidentified part, every value of it fits the observed data equally well, and
    the tipping point reports how far the assumption must move before the
    conclusion does.

:::

## Notes and sources

**The framework.** The classification into MCAR, MAR and MNAR, the selection
factorization and the ignorability theorem are due to Rubin (1976); the book-length
development is Little and Rubin (2019), whose organization this chapter follows, the
proofs, examples and numbers here being the book's own. Schafer (1997) treats the
multivariate normal and categorical cases, including the monotone factorizations of
[Section 41.3](03-likelihood.html), whose closed form is Anderson's (1957). The
realized-value and everywhere versions of MAR are separated by Seaman, Galati,
Jackson and Carlin (2013); the test of MCAR is Little (1988), and the general form of
@prp-mis-untestable is Molenberghs, Beunckens, Sotto and Kenward (2008).

**Complete cases.** Part (b) of @prp-mis-complete-case has been rediscovered many
times; the clearest statements are Little (1992) and White and Carlin (2010), who
also note that complete-case analysis can beat multiple imputation for efficiency
when the mechanism depends only on covariates. Inverse-probability weighting descends
from Horvitz and Thompson (1952), the doubly robust augmentation from Robins,
Rotnitzky and Zhao (1994).

**EM.** The algorithm was named and given its general form by Dempster, Laird and
Rubin (1977), who also identified the rate of convergence as the largest fraction of
missing information; special cases are older, the iteration of @thm-dsn-missing(c)
being the missing-plot scheme of Healy and Westmacott (1956), built on Yates's (1933)
formula. The missing-information principle is Orchard and Woodbury (1972), and
@eq-mis-louis is Louis's (1982). Wu (1983) corrected the convergence claims of the
original paper; Meng and Rubin (1991) extract the observed information from the EM
iterations themselves, and McLachlan and Krishnan (2008) is the book-length
treatment.

**Multiple imputation.** The idea is Rubin's (1978) and the definitive statement is
Rubin (1987), where the combining rules, the degrees of freedom and the efficiency
calculation of @thm-mis-rubin all appear; the small-sample degrees of freedom are
Barnard and Rubin (1999). Congeniality is Meng's (1994), whose discussion in the same
issue remains the best entry point to when Rubin's variance estimator is consistent.
Chained equations were developed by Van Buuren, Boshuizen and Knook (1999) and by
Raghunathan, Lepkowski, Van Hoewyk and Solenberger (2001); van Buuren (2018) is the
practical handbook.

**MNAR and sensitivity.** Selection models in the econometric tradition begin with
Heckman (1979), pattern-mixture models with Little (1993), who introduced the
complete-case missing-value restriction of @exr-mis-ccmv. The non-identifiability in
@prp-mis-sensitivity(c) is why Manski's (1990) assumption-free bounds, the subject of
@exr-mis-worst-case, remain worth computing. Delta-adjustment and tipping-point
analyses became standard in clinical trials after the National Research
Council (2010) report, which is short and unusually clear about the difference
between preventing missing data and repairing it. The reference-based schemes named
in [Section 41.6](06-sensitivity.html) are due to Carpenter, Roger and
Kenward (2013).

**Where this goes next.** Nothing here is special to a normal response: ignorability,
EM and Rubin's rules apply verbatim to the generalized linear models of Chapters 34
to 37, with harder integrals and an E step usually done by simulation. The
estimating-equation methods of
[Chapter 40](../ch40-glmm-gee/index.html) stay consistent when the missingness is
completely at random or depends only on modelled, fully observed covariates, and
become inconsistent when it depends on previously observed responses — the MAR case
in which likelihood remains valid — because the estimating function
\( \E\{R_i(Y_i-\mu_i)\} \) is then no longer zero. The repair is inverse-probability
weighting.

## References

- Anderson, Theodore W. (1957). Maximum Likelihood Estimates for a Multivariate Normal Distribution When Some Observations Are Missing. *Journal of the American Statistical Association* 52(278), 200–203.
- Barnard, John and Rubin, Donald B. (1999). Small-Sample Degrees of Freedom with Multiple Imputation. *Biometrika* 86(4), 948–955.
- Carpenter, James R., Roger, James H. and Kenward, Michael G. (2013). Analysis of Longitudinal Trials with Protocol Deviation: A Framework for Relevant, Accessible Assumptions, and Inference via Multiple Imputation. *Journal of Biopharmaceutical Statistics* 23(6), 1352–1371.
- Dempster, Arthur P., Laird, Nan M. and Rubin, Donald B. (1977). Maximum Likelihood from Incomplete Data via the EM Algorithm. *Journal of the Royal Statistical Society, Series B* 39(1), 1–38.
- Healy, Michael and Westmacott, Michael (1956). Missing Values in Experiments Analysed on Automatic Computers. *Journal of the Royal Statistical Society, Series C* 5(3), 203–206.
- Heckman, James J. (1979). Sample Selection Bias as a Specification Error. *Econometrica* 47(1), 153–161.
- Horvitz, Daniel G. and Thompson, D. J. (1952). A Generalization of Sampling Without Replacement from a Finite Universe. *Journal of the American Statistical Association* 47(260), 663–685.
- Little, Roderick J. A. (1988). A Test of Missing Completely at Random for Multivariate Data with Missing Values. *Journal of the American Statistical Association* 83(404), 1198–1202.
- Little, Roderick J. A. (1992). Regression with Missing X's: A Review. *Journal of the American Statistical Association* 87(420), 1227–1237.
- Little, Roderick J. A. (1993). Pattern-Mixture Models for Multivariate Incomplete Data. *Journal of the American Statistical Association* 88(421), 125–134.
- Little, Roderick J. A. and Rubin, Donald B. (2019). *Statistical Analysis with Missing Data*. 3rd edition. Hoboken, NJ: Wiley.
- Louis, Thomas A. (1982). Finding the Observed Information Matrix When Using the EM Algorithm. *Journal of the Royal Statistical Society, Series B* 44(2), 226–233.
- Manski, Charles F. (1990). Nonparametric Bounds on Treatment Effects. *American Economic Review* 80(2), 319–323.
- McLachlan, Geoffrey J. and Krishnan, Thriyambakam (2008). *The EM Algorithm and Extensions*. 2nd edition. Hoboken, NJ: Wiley.
- Meng, Xiao-Li (1994). Multiple-Imputation Inferences with Uncongenial Sources of Input. *Statistical Science* 9(4), 538–558.
- Meng, Xiao-Li and Rubin, Donald B. (1991). Using EM to Obtain Asymptotic Variance–Covariance Matrices: The SEM Algorithm. *Journal of the American Statistical Association* 86(416), 899–909.
- Molenberghs, Geert, Beunckens, Caroline, Sotto, Cristina and Kenward, Michael G. (2008). Every Missingness Not at Random Model Has a Missingness at Random Counterpart with Equal Fit. *Journal of the Royal Statistical Society, Series B* 70(2), 371–388.
- National Research Council (2010). *The Prevention and Treatment of Missing Data in Clinical Trials*. Washington, DC: National Academies Press.
- Orchard, Terence and Woodbury, Max A. (1972). A Missing Information Principle: Theory and Applications. *Proceedings of the Sixth Berkeley Symposium on Mathematical Statistics and Probability* 1, 697–715.
- Raghunathan, Trivellore E., Lepkowski, James M., Van Hoewyk, John and Solenberger, Peter (2001). A Multivariate Technique for Multiply Imputing Missing Values Using a Sequence of Regression Models. *Survey Methodology* 27(1), 85–95.
- Robins, James M., Rotnitzky, Andrea and Zhao, Lue Ping (1994). Estimation of Regression Coefficients When Some Regressors Are Not Always Observed. *Journal of the American Statistical Association* 89(427), 846–866.
- Rubin, Donald B. (1976). Inference and Missing Data. *Biometrika* 63(3), 581–592.
- Rubin, Donald B. (1978). Multiple Imputations in Sample Surveys — A Phenomenological Bayesian Approach to Nonresponse. *Proceedings of the Survey Research Methods Section, American Statistical Association*, 20–34.
- Rubin, Donald B. (1987). *Multiple Imputation for Nonresponse in Surveys*. New York: Wiley.
- Schafer, Joseph L. (1997). *Analysis of Incomplete Multivariate Data*. London: Chapman and Hall.
- Seaman, Shaun R., Galati, John, Jackson, Dan and Carlin, John (2013). What Is Meant by "Missing at Random"? *Statistical Science* 28(2), 257–268.
- Van Buuren, Stef (2018). *Flexible Imputation of Missing Data*. 2nd edition. Boca Raton, FL: Chapman and Hall/CRC.
- Van Buuren, Stef, Boshuizen, Hendriek C. and Knook, Dick L. (1999). Multiple Imputation of Missing Blood Pressure Covariates in Survival Analysis. *Statistics in Medicine* 18(6), 681–694.
- White, Ian R. and Carlin, John B. (2010). Bias and Efficiency of Multiple Imputation Compared with Complete-Case Analysis for Missing Covariate Values. *Statistics in Medicine* 29(28), 2920–2931.
- Wu, C. F. Jeff (1983). On the Convergence Properties of the EM Algorithm. *The Annals of Statistics* 11(1), 95–103.
- Yates, Frank (1933). The Analysis of Replicated Experiments When the Field Results Are Incomplete. *Empire Journal of Experimental Agriculture* 1(2), 129–142.
