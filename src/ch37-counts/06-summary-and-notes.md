# Summary and Notes

## Summary

::: {.idea}

1. The Poisson is the exponential dispersion family with \( V(\mu)=\mu \) and no
   free dispersion; its canonical link is the logarithm, so coefficients are log
   rate ratios (@def-cnt-poisson).

2. The likelihood equations are \( \X\T(\y-\hat{\bmu})=\bzero \), the information
   is \( \X\T\W\X \) with \( \W=\diag(\hat\mu_i) \), and the log-likelihood is
   strictly concave at full column rank (@thm-cnt-score). Fitted counts reproduce
   the observed margins spanned by the model matrix (@thm-cnt-marginals).

3. The maximum likelihood estimate exists exactly when no direction
   \( \mathbf{d} \) satisfies @eq-cnt-divergent (@prp-cnt-existence); a ridge
   penalty always restores one.

4. Exposure enters as an offset, a column with coefficient fixed at one. Grouping
   records with identical regressors and carrying group size as exposure leaves the
   fit unchanged; estimating the exposure coefficient costs precision, fixing a
   wrong one costs bias (@prp-cnt-offsets).

5. Independent Poisson counts conditioned on their total are multinomial, and the
   two fits agree whenever the model has a free parameter for every total the
   design fixed (@thm-cnt-poisson-multinomial). Independence is the no-interaction
   model (@cor-cnt-independence), and interaction parameters are log odds ratios
   (@prp-cnt-interaction).

6. A gamma mixture of Poissons is negative binomial, with variance
   \( \mu+\mu^2/\kappa \) (@lem-cnt-gamma-mixture, @def-cnt-negbin); any mixing
   distribution inflates the variance. At fixed \( \kappa \) the NB2 model is a
   generalized linear model, \( \bbeta \) and \( \kappa \) are orthogonal, and the
   test of \( \alpha=1/\kappa=0 \) is on the boundary, so its null distribution is
   \( \tfrac12\chi^2(0)+\tfrac12\chi^2(1) \); part (f) gives the NB1 score and its
   relation to quasi-Poisson (@prp-cnt-negbin). The score test @eq-cnt-score-test
   needs only the Poisson fit (@thm-cnt-dispersion-score).

7. Zero-inflated models add a mass at zero; hurdle models cut the distribution at
   zero and fit the two pieces separately, which factorizes the likelihood
   (@def-cnt-zero-models, @prp-cnt-zero-inflation). The Vuong statistic
   @eq-cnt-vuong is not valid for the nested comparison it is usually applied to.

:::

## Notes and sources

**Poisson regression.** The log-linear model for counts predates the generalized
linear model, but Nelder and Wedderburn (1972) put it in its modern setting, and
McCullagh and Nelder (1989) is the standard reference for the deviance, residuals
and fitting algorithm. Agresti (2015, chapter 7) and Fahrmeir et al. (2021,
section 5.2) are the sources against which this chapter's coverage was checked;
Cameron and Trivedi (2013) and Winkelmann (2008) are book-length econometric
treatments. The existence condition @prp-cnt-existence is the Poisson analogue of
the separation conditions of
[Chapter 35](../ch35-binary-responses/index.html), proved here from scratch. The
Anscombe residual of @exr-cnt-anscombe is from Anscombe (1948).

**Rates.** Treating survival data with piecewise-constant hazards as Poisson
counts with an offset is due to Holford (1980) and Laird and Olivier (1981); it is
the reason that much of the early software for survival analysis was software for
log-linear models.

**Contingency tables.** The Poisson–multinomial correspondence and the maximum
likelihood theory of hierarchical log-linear models are in Birch (1963) and,
comprehensively, Bishop, Fienberg and Holland (1975). Graphical log-linear models
are due to Darroch, Lauritzen and Speed (1980). The chi-squared test for
independence is Pearson (1900); @cor-cnt-independence is its likelihood-based
restatement.

**The negative binomial.** The gamma mixture derivation goes back to Greenwood
and Yule (1920), who used it for accident proneness. The NB1 and NB2
parameterizations and their different estimating equations are discussed by
Cameron and Trivedi (2013, section 3.3) and Hilbe (2011). The score test
@eq-cnt-score-test is due to Dean and Lawless (1989); the derivation here, by
expanding the log-likelihood in \( \alpha \), is the usual one. The boundary
problem is the setting of Self and Liang (1987), invoked by @prp-cnt-negbin(e);
@thm-mix-boundary computes a comparable case exactly.

**Zeros.** The zero-inflated Poisson model and its EM algorithm are Lambert
(1992), written for defect counts in manufacturing; hurdle models are Mullahy
(1986). Vuong (1989) introduced the model-selection statistic, for strictly
non-nested models; Wilson (2015) documents its misuse for zero-inflation testing,
and Schennach and Wilhelm (2017) analyse its size and propose a repair. The
definitions and likelihoods of the two model families are @def-cnt-zero-models and
the preamble to @prp-cnt-zero-inflation. [Chapter 38](../ch38-quasi-likelihood/index.html) returns to
overdispersion without a likelihood, and [Chapter 40](../ch40-glmm-gee/index.html) to counts that are
not independent.

## References

- Agresti, Alan (2015). *Foundations of Linear and Generalized Linear Models*. Hoboken, NJ: Wiley.
- Anscombe, F. J. (1948). The Transformation of Poisson, Binomial and Negative-Binomial Data. *Biometrika* 35(3–4), 246–254.
- Birch, M. W. (1963). Maximum Likelihood in Three-Way Contingency Tables. *Journal of the Royal Statistical Society, Series B* 25(1), 220–233.
- Bishop, Yvonne M. M., Fienberg, Stephen E. and Holland, Paul W. (1975). *Discrete Multivariate Analysis: Theory and Practice*. Cambridge, MA: MIT Press.
- Cameron, A. Colin and Trivedi, Pravin K. (2013). *Regression Analysis of Count Data*. 2nd edition. Cambridge: Cambridge University Press.
- Cox, David R. and Reid, Nancy (1987). Parameter Orthogonality and Approximate Conditional Inference. *Journal of the Royal Statistical Society, Series B* 49(1), 1–39.
- Darroch, J. N., Lauritzen, S. L. and Speed, T. P. (1980). Markov Fields and Log-Linear Interaction Models for Contingency Tables. *The Annals of Statistics* 8(3), 522–539.
- Dean, C. and Lawless, J. F. (1989). Tests for Detecting Overdispersion in Poisson Regression Models. *Journal of the American Statistical Association* 84(406), 467–472.
- Fahrmeir, Ludwig, Kneib, Thomas, Lang, Stefan and Marx, Brian D. (2021). *Regression: Models, Methods and Applications*. 2nd edition. Berlin: Springer.
- Greenwood, Major and Yule, G. Udny (1920). An Inquiry into the Nature of Frequency Distributions Representative of Multiple Happenings. *Journal of the Royal Statistical Society* 83(2), 255–279.
- Hilbe, Joseph M. (2011). *Negative Binomial Regression*. 2nd edition. Cambridge: Cambridge University Press.
- Holford, Theodore R. (1980). The Analysis of Rates and of Survivorship Using Log-Linear Models. *Biometrics* 36(2), 299–305.
- Laird, Nan and Olivier, Donald (1981). Covariance Analysis of Censored Survival Data Using Log-Linear Analysis Techniques. *Journal of the American Statistical Association* 76(374), 231–240.
- Lambert, Diane (1992). Zero-Inflated Poisson Regression, with an Application to Defects in Manufacturing. *Technometrics* 34(1), 1–14.
- McCullagh, Peter and Nelder, John A. (1989). *Generalized Linear Models*. 2nd edition. London: Chapman and Hall.
- Mullahy, John (1986). Specification and Testing of Some Modified Count Data Models. *Journal of Econometrics* 33(3), 341–365.
- Nelder, John A. and Wedderburn, Robert W. M. (1972). Generalized Linear Models. *Journal of the Royal Statistical Society, Series A* 135(3), 370–384.
- Pearson, Karl (1900). On the Criterion that a Given System of Deviations from the Probable in the Case of a Correlated System of Variables is Such that It Can Be Reasonably Supposed to Have Arisen from Random Sampling. *Philosophical Magazine, Series 5* 50(302), 157–175.
- Schennach, Susanne M. and Wilhelm, Daniel (2017). A Simple Parametric Model Selection Test. *Journal of the American Statistical Association* 112(520), 1663–1674.
- Self, Steven G. and Liang, Kung-Yee (1987). Asymptotic Properties of Maximum Likelihood Estimators and Likelihood Ratio Tests under Nonstandard Conditions. *Journal of the American Statistical Association* 82(398), 605–610.
- Vuong, Quang H. (1989). Likelihood Ratio Tests for Model Selection and Non-Nested Hypotheses. *Econometrica* 57(2), 307–333.
- Wilson, Paul (2015). The Misuse of the Vuong Test for Non-Nested Models to Test for Zero-Inflation. *Economics Letters* 127, 51–53.
- Winkelmann, Rainer (2008). *Econometric Analysis of Count Data*. 5th edition. Berlin: Springer.
