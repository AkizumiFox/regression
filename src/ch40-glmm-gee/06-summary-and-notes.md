# Summary and Notes

## Summary

::: {.idea}

1. A nonlinear link and a random effect do not commute, so a clustered model has two sets
   of coefficients. The identity and log links are the exceptions; for the probit the
   marginal coefficient is \( \beta^{C}/\sqrt{1+\tau^{2}} \) exactly, and for the logit it
   is strictly attenuated with local factor
   \( 1-\Var(P_\theta)/\{\E P_\theta\E(1-P_\theta)\} \) (@prp-gmm-marginal).

2. The generalized linear mixed model puts random effects in the linear predictor of a
   generalized linear model (@def-gmm-glmm). Its marginal moments follow by conditioning:
   a binomial response acquires the beta-binomial overdispersion factor, and a
   Poisson–lognormal response acquires *exactly* the NB2 variance
   function (@prp-gmm-moments).

3. A cluster-level random intercept can only induce nonnegative within-cluster correlation,
   because two nondecreasing functions of one random variable have nonnegative
   covariance (@lem-gmm-association). Data with competition inside a cluster need a
   different tool.

4. Predicting a random effect in a GLMM is empirical Bayes rather than BLUP, but the
   shrinkage survives: the posterior mode lies strictly between zero and the cluster's own
   estimate, with the familiar factor \( \tau^{2}/(\tau^{2}+I_i^{-1}) \) to first
   order (@prp-gmm-prediction).

5. The integrated likelihood factorizes over independent clusters. The Laplace
   approximation has relative error \( O(n_i^{-1}) \) — governed by the observations *per
   cluster*, not by the sample size — and adaptive Gauss–Hermite quadrature with \( K=1 \)
   node is exactly the Laplace approximation (@thm-gmm-likelihood).

6. Penalized quasi-likelihood is Henderson's mixed-model equations applied to the working
   response, and is inconsistent: it attenuates both \( \hbeta \) and
   \( \hat{\boldsymbol{\uptheta}} \), badly for binary data in small clusters (@prp-gmm-pql, @exm-gmm-compare).

7. Generalized estimating equations solve \( \sum_i\bD_i\T\V_i^{-1}(\Y_i-\bmu_i)=\bzero \)
   with a working correlation the analyst chooses and does not have to believe (@def-gmm-gee).
   The estimating function is unbiased for every working correlation, so
   \( \hbeta \) is consistent for the marginal coefficient and asymptotically normal with
   the sandwich covariance \( \A^{-1}\B\A^{-1} \) (@thm-gmm-gee).

8. GEE gives no likelihood, no random-effect predictions and no distribution; a
   non-diagonal working correlation also needs the mean model to hold conditionally on the
   *whole* covariate history, which feedback destroys; and it fails under missingness that
   depends on previously observed responses — the one case in which a likelihood-based
   analysis still works ([Section 40.4](04-gee.html)).

9. A wrong working correlation costs efficiency for covariates that vary within clusters
   and nothing at all for covariates fixed within a cluster; the sandwich needs a
   small-sample correction below about forty clusters; and QIC, the quasi-likelihood
   analogue of AIC, is sound for mean models and weak for correlation structures (@prp-gmm-working).

10. Because every working correlation targets the same \( \bbeta^{0} \), a large gap between
    two GEE fits of the same mean model is evidence against the *mean model* (@exm-gmm-qic).

:::

## Notes and sources

**Coverage.** The topics follow Agresti (2015, sections 9.4–9.6) and Fahrmeir, Kneib, Lang
and Marx (2021, sections 7.5–7.6), but put the marginal–conditional distinction first,
before either estimation method, because that is where the choice between the two traditions
is actually made. Book-length treatments are Diggle, Heagerty, Liang and Zeger (2002) for
longitudinal data, McCulloch, Searle and Neuhaus (2008) for the mixed-model theory, and
Molenberghs and Verbeke (2005) for discrete longitudinal responses.

**Marginal versus conditional.** The distinction was made precise for binary data by
Neuhaus, Kalbfleisch and Hauck (1991), still the clearest statement of what each coefficient
means and when the two agree. The approximation factor \( (1+c^{2}\tau^{2})^{-1/2} \) with
\( c=16\sqrt3/(15\pi) \) is from Zeger, Liang and Albert (1988); the constant comes from
matching the logistic and normal distributions and predates its use here by a long way. The
proof of @prp-gmm-marginal(d) given above, through the derivative identity
@eq-gmm-attenuation-derivative, seems to be the shortest route to the strict inequality, and
needs no assumption on the random-effect distribution beyond non-degeneracy. The median odds
ratio of @exr-gmm-mor was popularized in epidemiology by Larsen and Merlo (2005).

**Fitting.** Breslow and Clayton (1993) introduced penalized quasi-likelihood under that
name and set the pattern for GLMM software for a decade; Breslow and Lin (1995) and Lin and
Breslow (1996) quantified its bias and proposed corrections. Adaptive Gauss–Hermite
quadrature for this problem is due to Liu and Pierce (1994) and, for nonlinear mixed models,
Pinheiro and Bates (1995), who also observed that one adaptive node reproduces the Laplace
approximation. The Laplace expansion itself is classical; the form used in
@thm-gmm-likelihood(b) is the one in Tierney and Kadane (1986) and Barndorff-Nielsen and Cox
(1989, chapter 3). The asymptotics of the integrated likelihood quoted in
[Section 40.2](02-glmms.html) are in McCulloch, Searle and Neuhaus (2008, chapter 14). Monte
Carlo EM for this problem is Wei and Tanner (1990), with automated sample-size rules in
Booth and Hobert (1999); the MCMC route with iteratively reweighted proposals is Gamerman
(1997).

**Estimating equations.** Liang and Zeger (1986) and Zeger and Liang (1986) introduced GEE,
and the pair is worth reading together: the first gives the theory, the second the practice.
The conditions and the expansion in @thm-gmm-gee follow theirs, with the cluster sandwich
argument taken from @thm-het-sandwich as
[Section 33.2](../ch33-clustered-longitudinal-splitplot/02-clustered-data.html) already
arranged it. The warning about time-varying covariates and the full-covariate conditional
mean is due to Pepe and Anderson (1994), and is treated at length in Diggle, Heagerty, Liang
and Zeger (2002, chapter 12). Small-sample corrections to the sandwich are due to Mancl and
DeRouen (2001), Kauermann and Carroll (2001) and Fay and Graubard (2001); the first is used
in @exm-gmm-coverage. QIC is Pan (2001); the criticism of it as a selector of working
correlations, and the CIC variant, are in Hin and Wang (2009). Alternating logistic
regressions are Carey, Zeger and Diggle (1993). The failure of GEE under missingness that
depends on observed responses, and the inverse-probability-weighted repair, are from Robins,
Rotnitzky and Zhao (1995).

**What is new here.** @lem-gmm-association and its consequence — that a cluster-level random
intercept cannot produce negative within-cluster correlation, so that the travel-mode data
of @exm-gmm-modechoice are outside the reach of any random-intercept GLMM — is elementary,
and the association inequality behind it is standard, but we have not seen the two put
together as a reason to prefer estimating equations. The efficiency comparison of
@exm-gmm-efficiency separates within-cluster from cluster-level covariates, making the
often-repeated advice "the working correlation hardly matters" precise: it hardly matters
for cluster-level covariates and matters by a factor of nearly two for within-cluster ones.

**Data and code.** Grunfeld's investment panel is the public-domain eleven-firm version
distributed with statsmodels, also used in [Chapter 32](../ch32-linear-mixed-models/index.html)
and [Chapter 33](../ch33-clustered-longitudinal-splitplot/index.html); its provenance is
Grunfeld (1958), with the version history discussed by Kleiber and Zeileis (2010). The
travel-mode choice data are the public-domain sample of \( 210 \) non-business trips between
Sydney, Canberra and Melbourne from a 1987 intercity mode-choice study, also distributed with
statsmodels; the sample over-represents the rarer modes, as @exm-gmm-modechoice notes. The
code in `code/ch40/` is written from the formulas in this chapter and checked against
statsmodels where an equivalent fit exists; the simulations use fixed seeds.

## References

- Agresti, Alan (2015). *Foundations of Linear and Generalized Linear Models*. Hoboken, NJ: Wiley.
- Barndorff-Nielsen, Ole E. and Cox, David R. (1989). *Asymptotic Techniques for Use in Statistics*. London: Chapman and Hall.
- Booth, James G. and Hobert, James P. (1999). Maximizing Generalized Linear Mixed Model Likelihoods with an Automated Monte Carlo EM Algorithm. *Journal of the Royal Statistical Society, Series B* 61(1), 265–285.
- Breslow, Norman E. and Clayton, David G. (1993). Approximate Inference in Generalized Linear Mixed Models. *Journal of the American Statistical Association* 88(421), 9–25.
- Breslow, Norman E. and Lin, Xihong (1995). Bias Correction in Generalised Linear Mixed Models with a Single Component of Dispersion. *Biometrika* 82(1), 81–91.
- Carey, Vincent, Zeger, Scott L. and Diggle, Peter (1993). Modelling Multivariate Binary Data with Alternating Logistic Regressions. *Biometrika* 80(3), 517–526.
- Diggle, Peter J., Heagerty, Patrick, Liang, Kung-Yee and Zeger, Scott L. (2002). *Analysis of Longitudinal Data*. 2nd edition. Oxford: Oxford University Press.
- Fahrmeir, Ludwig, Kneib, Thomas, Lang, Stefan and Marx, Brian D. (2021). *Regression: Models, Methods and Applications*. 2nd edition. Berlin: Springer.
- Fay, Michael P. and Graubard, Barry I. (2001). Small-Sample Adjustments for Wald-Type Tests Using Sandwich Estimators. *Biometrics* 57(4), 1198–1206.
- Gamerman, Dani (1997). Sampling from the Posterior Distribution in Generalized Linear Mixed Models. *Statistics and Computing* 7(1), 57–68.
- Grunfeld, Yehuda (1958). *The Determinants of Corporate Investment*. PhD thesis, University of Chicago.
- Hin, Lin-Yee and Wang, You-Gan (2009). Working-Correlation-Structure Identification in Generalized Estimating Equations. *Statistics in Medicine* 28(4), 642–658.
- Kauermann, Göran and Carroll, Raymond J. (2001). A Note on the Efficiency of Sandwich Covariance Matrix Estimation. *Journal of the American Statistical Association* 96(456), 1387–1396.
- Kleiber, Christian and Zeileis, Achim (2010). The Grunfeld Data at 50. *German Economic Review* 11(4), 404–417.
- Larsen, Klaus and Merlo, Juan (2005). Appropriate Assessment of Neighborhood Effects on Individual Health: Integrating Random and Fixed Effects in Multilevel Logistic Regression. *American Journal of Epidemiology* 161(1), 81–88.
- Liang, Kung-Yee and Zeger, Scott L. (1986). Longitudinal Data Analysis Using Generalized Linear Models. *Biometrika* 73(1), 13–22.
- Lin, Xihong and Breslow, Norman E. (1996). Bias Correction in Generalized Linear Mixed Models with Multiple Components of Dispersion. *Journal of the American Statistical Association* 91(435), 1007–1016.
- Liu, Qing and Pierce, Donald A. (1994). A Note on Gauss–Hermite Quadrature. *Biometrika* 81(3), 624–629.
- Mancl, Lloyd A. and DeRouen, Timothy A. (2001). A Covariance Estimator for GEE with Improved Small-Sample Properties. *Biometrics* 57(1), 126–134.
- McCulloch, Charles E., Searle, Shayle R. and Neuhaus, John M. (2008). *Generalized, Linear, and Mixed Models*. 2nd edition. Hoboken, NJ: Wiley.
- Molenberghs, Geert and Verbeke, Geert (2005). *Models for Discrete Longitudinal Data*. New York: Springer.
- Neuhaus, John M., Kalbfleisch, John D. and Hauck, Walter W. (1991). A Comparison of Cluster-Specific and Population-Averaged Approaches for Analyzing Correlated Binary Data. *International Statistical Review* 59(1), 25–35.
- Pan, Wei (2001). Akaike's Information Criterion in Generalized Estimating Equations. *Biometrics* 57(1), 120–125.
- Pepe, Margaret Sullivan and Anderson, Garnet L. (1994). A Cautionary Note on Inference for Marginal Regression Models with Longitudinal Data and General Correlated Response Data. *Communications in Statistics — Simulation and Computation* 23(4), 939–951.
- Pinheiro, José C. and Bates, Douglas M. (1995). Approximations to the Log-Likelihood Function in the Nonlinear Mixed-Effects Model. *Journal of Computational and Graphical Statistics* 4(1), 12–35.
- Robins, James M., Rotnitzky, Andrea and Zhao, Lue Ping (1995). Analysis of Semiparametric Regression Models for Repeated Outcomes in the Presence of Missing Data. *Journal of the American Statistical Association* 90(429), 106–121.
- Tierney, Luke and Kadane, Joseph B. (1986). Accurate Approximations for Posterior Moments and Marginal Densities. *Journal of the American Statistical Association* 81(393), 82–86.
- Wei, Greg C. G. and Tanner, Martin A. (1990). A Monte Carlo Implementation of the EM Algorithm and the Poor Man's Data Augmentation Algorithms. *Journal of the American Statistical Association* 85(411), 699–704.
- Zeger, Scott L. and Liang, Kung-Yee (1986). Longitudinal Data Analysis for Discrete and Continuous Outcomes. *Biometrics* 42(1), 121–130.
- Zeger, Scott L., Liang, Kung-Yee and Albert, Paul S. (1988). Models for Longitudinal Data: A Generalized Estimating Equation Approach. *Biometrics* 44(4), 1049–1060.
