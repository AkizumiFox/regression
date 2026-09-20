# Summary and Notes

## Summary

::: {.idea}

1. A factor is random when its levels are a sample and the conclusion is meant for
   the population of levels. The one-way random effects model has covariance
   \( \sigma^2\I+\sigma_a^2\Z\Z\T \), intraclass correlation
   \( \sigma_a^2/(\sigma_a^2+\sigma^2) \), and unbiased moment estimates from the
   analysis of variance (@def-mix-oneway and @prp-mix-oneway).

2. In a balanced normal layout the between and within sums of squares are
   independent scaled chi-squares, so \( \text{MSB}/\text{MSE} \) has an exact
   \( F \) distribution for every value of the variance ratio. Inverting it gives an
   exact interval for \( \sigma_a^2/\sigma^2 \) and for the intraclass
   correlation (@cor-mix-icc-interval).

3. The moment estimate of a variance can be negative, with probability
   \( \Pr\{F(g-1,n-g)<1/(1+m\gamma)\} \), which is just over one half when
   \( \sigma_a^2=0 \). Reporting, truncating and re-thinking the model are all
   defensible responses (@eq-mix-negative).

4. The linear mixed model \( \Y=\X\bbeta+\Z\bu+\be \) has marginal mean
   \( \X\bbeta \) and covariance \( \V=\Z\G\Z\T+\R \), a structured member of the
   general Gauss–Markov family. Woodbury's identity inverts \( \V \) at the cost of a
   \( q\times q \) inverse (@def-mix-model and @prp-mix-vinverse).

5. Ordinary least squares is efficient exactly when \( \C(\V\X)\subseteq\C(\X) \),
   which holds for cluster-constant regressors in a balanced design and fails as soon
   as a regressor varies within clusters; the mixed estimate then interpolates
   between the within and between estimators (@exm-mix-ols-blue
   and @exm-mix-ols-not-blue).

6. The best linear unbiased predictor of \( \blambda\T\bbeta+\mathbf{c}\T\bu \) is
   \( \blambda\T\hbeta+\mathbf{c}\T\G\Z\T\V^{-1}(\y-\X\hbeta) \): the best linear
   predictor with \( \bbeta \) replaced by its generalized least squares estimate. Its
   prediction error variance has two terms, one for predicting \( \bu \) and one for
   estimating \( \bbeta \) (@thm-mix-blup, which is @thm-cor-blup of
   [Section 14.2](../ch14-correlation-lack-of-fit-prediction/02-best-linear-prediction.html)
   written for mixed functions).

7. In the balanced one-way model the BLUP is \( B(\bar y_k-\bar y) \) with
   \( B=m\gamma/(1+m\gamma) \): shrinkage towards the overall mean, of the empirical
   Bayes kind, with a genuine gain in mean squared error and a conditional bias that
   matters for ranking (@eq-mix-shrinkage).

8. Henderson's mixed model equations solve for \( (\hbeta,\hat{\bu}) \) in one sparse
   system of size \( p+q \); they are the normal equations of a penalized least
   squares criterion with penalty \( \bu\T\G^{-1}\bu \), and the inverse of their
   coefficient matrix holds both \( \Cov(\hbeta) \) and the prediction error
   variance (@thm-mix-henderson).

9. Restricted maximum likelihood is the likelihood of a maximal set of error
   contrasts. It does not depend on which contrasts are used, it equals the profile
   likelihood plus \( -\tfrac12\log\det(\X\T\V^{-1}\X) \) up to constants, and it
   reduces to \( s^2=\text{SSE}/(n-p) \) when \( \V=\sigma^2\I \) (@def-mix-reml
   and @thm-mix-reml).

10. In the balanced one-way model, the interior maximum likelihood estimate
    underestimates \( \sigma_a^2 \) by exactly \( \lambda_1/(gm) \), while REML
    reproduces the
    unbiased moment estimates. The two score equations differ only in
    \( \tr(\V^{-1}\V_j) \) against \( \tr(\bP\V_j) \) (@prp-mix-oneway-reml
    and @prp-mix-scores).

11. Feasible generalized least squares is still unbiased when the variance estimate
    is translation invariant and even. Its Wald statistic is not normal: in a
    balanced cluster-randomized trial it is exactly \( t(g-2) \), and treating it as
    normal costs real coverage. Satterthwaite and Kenward–Roger approximate the right
    reference distribution in general (@prp-mix-inference).

12. Testing \( \sigma_a^2=0 \) puts the null on the boundary. In the balanced one-way
    model the restricted likelihood ratio is an explicit increasing function of the
    \( F \) ratio, zero with probability \( \Pr(F\le1) \), and its exact null
    distribution follows; asymptotically it is the chi-bar-squared mixture
    \( \tfrac12\chi^2(0)+\tfrac12\chi^2(1) \) (@thm-mix-boundary, the
    variance-component companion of @prp-mix-inference).

13. Restricted likelihoods may be compared only across covariance structures with the
    same fixed effects; comparisons of fixed effects need the ordinary likelihood
    ([Section 32.5](05-inference.html)).

14. In the hierarchical reading, the BLUP is a posterior mean and Henderson's
    coefficient matrix is a posterior precision; the restricted likelihood is the
    marginal posterior of the variance components under a flat prior on \( \bbeta \).
    A plug-in analysis omits the variance of the conditional mean, so its intervals
    are too short (@prp-mix-bayes).

:::

## Notes and sources

**Fixed and random effects.**  The separation of variability into components goes
back to Fisher (1918), who decomposed the correlation between relatives. Eisenhart
(1947) fixed the terminology "model I" and "model II" and made explicit that the
choice changes the expected mean squares and therefore which ratio is a test.
[Section 32.1](01-random-effects.html) follows the expected-mean-squares route of
[Chapter 9](../ch09-sums-of-squares/index.html), and @eq-mix-m0 is the first of
Henderson's (1953) three methods for unbalanced data. Searle, Casella and McCulloch
(1992) is the standard reference for variance component estimation, including methods
this chapter omits such as Rao's (1971) MINQUE; Burdick and Graybill (1992) treats
intervals for variance components at book length.

**Prediction.**  @thm-mix-blup is
[Section 14.2](../ch14-correlation-lack-of-fit-prediction/02-best-linear-prediction.html)'s
theorem with a random effect in place of the new observation's error, which is why
its proof is a substitution. Best linear unbiased prediction
is due to Henderson, in a sequence
of papers beginning with the abstract Henderson (1950) and developed in Henderson,
Kempthorne, Searle and von Krosigk (1959) and Henderson (1975); the motivating problem
was ranking dairy sires from unbalanced records, where shrinkage is the difference
between a usable and an unusable ranking. Harville (1976, 1977) gave the modern
treatment, including the relation to Bayesian prediction. Robinson (1991) is the best
single essay on what BLUP is and why several derivations converge on it, and the
warning in [Section 32.3](03-blup.html) about ranking shrunk predictions is developed
there. The equivalence with penalized least squares, and hence with ridge
regression (@thm-shr-ridge), is the bridge to smoothing: Ruppert, Wand and Carroll (2003) build
semiparametric regression on it, and Chapter 43 uses it to estimate a smoothing
parameter by REML.

**REML.**  Restricted maximum likelihood was introduced by Patterson and Thompson
(1971) for recovering inter-block information in incomplete block designs, the problem
Yates (1940) had posed; [Section 18.5](../ch18-covariance-and-design/05-latin-squares-incomplete-blocks.html)
met it in its fixed-effects form. Harville (1974) proved the marginal-posterior
characterization that @prp-mix-bayes(b) restates, and Harville (1977) surveyed the
area; Hartley and Rao (1967) had treated maximum likelihood for the mixed model a few
years earlier, and Corbeil and Searle (1976) established REML's advantage for small
designs by simulation. The determinant identity @eq-mix-reml-determinant and the identity
\( \mathbf{K}(\mathbf{K}\T\V\mathbf{K})^{-1}\mathbf{K}\T=\bP \) are standard; the
proof of the latter given here, through a constrained minimization, seems the
shortest route and makes @lem-mix-P(c) immediate. On algorithms, Dempster, Laird and
Rubin (1977) supplied EM, Laird and Ware (1982) brought it to longitudinal data, and
Gilmour, Thompson and Cullis (1995) introduced average information REML. Bates,
Mächler, Bolker and Walker (2015) describe a sparse-matrix implementation in which the
profiled deviance is evaluated by a single sparse Cholesky update.

**Inference.**  The unbiasedness result @prp-mix-inference(a) is the covariance analogue
of the weighted least squares argument of
[Section 21.3](../ch21-nonnormality-heteroscedasticity-serial/03-weighted-least-squares.html);
Kackar and Harville (1981) proved it in this generality, and Kackar and Harville
(1984) gave the approximation to the extra variance of
\( \hbeta(\hat{\boldsymbol{\uptheta}}) \) that Kenward and Roger (1997) turned into a
usable small-sample correction. The degrees-of-freedom
approximation of @eq-mix-satterthwaite is Satterthwaite (1946); Giesbrecht and Burns
(1985) extended it to mixed models. For the boundary problem, Self and Liang (1987)
gave the general asymptotic theory of likelihood ratios when the true parameter is on
the boundary, and Stram and Lee (1994) worked out the mixed-model cases most often
needed. Crainiceanu and Ruppert (2004) showed that for a single variance component the
exact finite-sample null distribution of the restricted likelihood ratio can be
computed, and documented how poor the chi-bar-squared approximation can be; the exact
balanced calculation of @thm-mix-boundary is a special case that can be done by hand,
and the proof sketch of the chi-bar-squared limit given after it is our own.
One question left to the econometric literature is whether the random effects are
uncorrelated with the regressors at all, as @def-mix-model assumes; if they are not,
the mixed fit is inconsistent and the within fit of @exm-mix-ols-not-blue is not.
Hausman (1978) tests this by comparing the two, using the fact that when one estimator
is the BLUE the covariance of the difference is the difference of the covariances.
Verbeke and Molenberghs (2000) and Pinheiro and Bates (2000) are the standard applied
references, the first for the theory of inference in linear mixed models and the second
for practice; McCulloch, Searle and Neuhaus (2008) covers both linear and generalized
mixed models.

**Bayesian mixed models.**  Gelman (2006) changed practice on priors for variance
parameters, showing that the inverse gamma with small parameters is not a safe default
and recommending half-\( t \) priors on standard deviations; the scale-mixture
representation used in our sampler is from that paper.

**Sources for the chapter.**  The coverage checklist was built from
Rencher and Schaalje (2008, chapter 17), Agresti (2015, sections 9.2–9.3),
Fahrmeir, Kneib, Lang and Marx (2021, sections 7.1–7.4) and
Sen and Srivastava (1990, sections 7.4–7.5). All the text, proofs, examples, data and
exercises here are our own. The proficiency study is simulated from stated parameter
values with a fixed seed; the investment panel is Grunfeld's, in the public-domain
reconstruction distributed with statsmodels, and continues the analysis begun in
[Section 6.6](../ch06-projections/06-fwl.html).

## References

- Agresti, Alan (2015). *Foundations of Linear and Generalized Linear Models*. Hoboken, NJ: Wiley.
- Bates, Douglas, Mächler, Martin, Bolker, Ben and Walker, Steve (2015). Fitting Linear Mixed-Effects Models Using lme4. *Journal of Statistical Software* 67(1), 1–48.
- Burdick, Richard K. and Graybill, Franklin A. (1992). *Confidence Intervals on Variance Components*. New York: Marcel Dekker.
- Corbeil, Robert R. and Searle, Shayle R. (1976). Restricted Maximum Likelihood (REML) Estimation of Variance Components in the Mixed Model. *Technometrics* 18(1), 31–38.
- Crainiceanu, Ciprian M. and Ruppert, David (2004). Likelihood Ratio Tests in Linear Mixed Models with One Variance Component. *Journal of the Royal Statistical Society, Series B* 66(1), 165–185.
- Dempster, Arthur P., Laird, Nan M. and Rubin, Donald B. (1977). Maximum Likelihood from Incomplete Data via the EM Algorithm. *Journal of the Royal Statistical Society, Series B* 39(1), 1–38.
- Eisenhart, Churchill (1947). The Assumptions Underlying the Analysis of Variance. *Biometrics* 3(1), 1–21.
- Fahrmeir, Ludwig, Kneib, Thomas, Lang, Stefan and Marx, Brian D. (2021). *Regression: Models, Methods and Applications*. 2nd edition. Berlin: Springer.
- Fisher, Ronald A. (1918). The Correlation between Relatives on the Supposition of Mendelian Inheritance. *Transactions of the Royal Society of Edinburgh* 52(2), 399–433.
- Gelman, Andrew (2006). Prior Distributions for Variance Parameters in Hierarchical Models. *Bayesian Analysis* 1(3), 515–534.
- Giesbrecht, Francis G. and Burns, Joseph C. (1985). Two-Stage Analysis Based on a Mixed Model: Large-Sample Asymptotic Theory and Small-Sample Simulation Results. *Biometrics* 41(2), 477–486.
- Gilmour, Arthur R., Thompson, Robin and Cullis, Brian R. (1995). Average Information REML: An Efficient Algorithm for Variance Parameter Estimation in Linear Mixed Models. *Biometrics* 51(4), 1440–1450.
- Hartley, Herman O. and Rao, Jon N. K. (1967). Maximum-Likelihood Estimation for the Mixed Analysis of Variance Model. *Biometrika* 54(1/2), 93–108.
- Harville, David A. (1974). Bayesian Inference for Variance Components Using Only Error Contrasts. *Biometrika* 61(2), 383–385.
- Harville, David A. (1976). Extension of the Gauss–Markov Theorem to Include the Estimation of Random Effects. *The Annals of Statistics* 4(2), 384–395.
- Harville, David A. (1977). Maximum Likelihood Approaches to Variance Component Estimation and to Related Problems. *Journal of the American Statistical Association* 72(358), 320–338.
- Hausman, Jerry A. (1978). Specification Tests in Econometrics. *Econometrica* 46(6), 1251–1271.
- Henderson, Charles R. (1950). Estimation of Genetic Parameters (abstract). *The Annals of Mathematical Statistics* 21(2), 309–310.
- Henderson, Charles R. (1953). Estimation of Variance and Covariance Components. *Biometrics* 9(2), 226–252.
- Henderson, Charles R. (1975). Best Linear Unbiased Estimation and Prediction under a Selection Model. *Biometrics* 31(2), 423–447.
- Henderson, Charles R., Kempthorne, Oscar, Searle, Shayle R. and von Krosigk, C. M. (1959). The Estimation of Environmental and Genetic Trends from Records Subject to Culling. *Biometrics* 15(2), 192–218.
- Kackar, Raghu N. and Harville, David A. (1981). Unbiasedness of Two-Stage Estimation and Prediction Procedures for Mixed Linear Models. *Communications in Statistics — Theory and Methods* 10(13), 1249–1261.
- Kackar, Raghu N. and Harville, David A. (1984). Approximations for Standard Errors of Estimators of Fixed and Random Effects in Mixed Linear Models. *Journal of the American Statistical Association* 79(388), 853–862.
- Kenward, Michael G. and Roger, James H. (1997). Small Sample Inference for Fixed Effects from Restricted Maximum Likelihood. *Biometrics* 53(3), 983–997.
- Laird, Nan M. and Ware, James H. (1982). Random-Effects Models for Longitudinal Data. *Biometrics* 38(4), 963–974.
- McCulloch, Charles E., Searle, Shayle R. and Neuhaus, John M. (2008). *Generalized, Linear, and Mixed Models*. 2nd edition. Hoboken, NJ: Wiley.
- Patterson, H. Desmond and Thompson, Robin (1971). Recovery of Inter-Block Information when Block Sizes are Unequal. *Biometrika* 58(3), 545–554.
- Pinheiro, José C. and Bates, Douglas M. (2000). *Mixed-Effects Models in S and S-PLUS*. New York: Springer.
- Rao, C. Radhakrishna (1971). Estimation of Variance and Covariance Components — MINQUE Theory. *Journal of Multivariate Analysis* 1(3), 257–275.
- Rencher, Alvin C. and Schaalje, G. Bruce (2008). *Linear Models in Statistics*. 2nd edition. Hoboken, NJ: Wiley.
- Robinson, George K. (1991). That BLUP is a Good Thing: The Estimation of Random Effects. *Statistical Science* 6(1), 15–32.
- Ruppert, David, Wand, M. P. and Carroll, Raymond J. (2003). *Semiparametric Regression*. Cambridge: Cambridge University Press.
- Satterthwaite, Franklin E. (1946). An Approximate Distribution of Estimates of Variance Components. *Biometrics Bulletin* 2(6), 110–114.
- Searle, Shayle R., Casella, George and McCulloch, Charles E. (1992). *Variance Components*. New York: Wiley.
- Self, Steven G. and Liang, Kung-Yee (1987). Asymptotic Properties of Maximum Likelihood Estimators and Likelihood Ratio Tests under Nonstandard Conditions. *Journal of the American Statistical Association* 82(398), 605–610.
- Sen, Ashish and Srivastava, Muni (1990). *Regression Analysis: Theory, Methods, and Applications*. New York: Springer.
- Stram, Daniel O. and Lee, Jae Won (1994). Variance Components Testing in the Longitudinal Mixed Effects Model. *Biometrics* 50(4), 1171–1177.
- Verbeke, Geert and Molenberghs, Geert (2000). *Linear Mixed Models for Longitudinal Data*. New York: Springer.
- Yates, Frank (1940). The Recovery of Inter-Block Information in Balanced Incomplete Block Designs. *Annals of Eugenics* 10(4), 317–325.
