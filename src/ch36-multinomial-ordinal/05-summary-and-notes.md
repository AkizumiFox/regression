# Summary and Notes

## Summary

::: {.idea}

1. The multinomial distribution is an exponential family with natural parameters the
           baseline-category logits, cumulant \( b(\boldsymbol{\upeta})=\log(1+\sum_s e^{\eta_s}) \),
           mean \( m\nabla b \) and covariance
           \( m(\diag(\boldsymbol{\uppi})-\boldsymbol{\uppi}\boldsymbol{\uppi}\T) \)
           (@lem-mlt-multinomial). The chapter then chooses which \( c-1 \) logits to make
           linear.

2. The baseline-category logit model (@def-mlt-baseline) is the canonical-link
           generalized linear model for that family: score \( \X\T(\y_r-\bmu_r) \) per category,
           observed information equal to expected, and a log-likelihood strictly concave when
           \( \rank(\X)=p \), so the maximum is unique when it exists (@thm-mlt-score).

3. Changing the baseline is a linear reparameterization, leaving estimates, fitted
           values, the maximized likelihood and the standard errors of contrasts unchanged
           (@prp-mlt-invariance). A coefficient is a log odds ratio against the baseline,
           not a derivative of a probability (@prp-mlt-derivative).

4. For an ordered response, modelling \( P(Y\le r) \) uses the ordering
           (@def-mlt-cumulative). Fitted probabilities are positive only where the linear
           predictors stay in order: automatic with common slopes, a real restriction
           without them (@prp-mlt-ordering). Score and expected information are
           \( \sum_i\bD_i\T\V_i^{-1}(\tilde{\y}_i-m_i\tilde{\boldsymbol{\uppi}}_i) \) and
           \( \sum_i m_i\bD_i\T\V_i^{-1}\bD_i \) (@prp-mlt-cumulative-score), and the
           log-likelihood is concave when the link density is log-concave.

5. The proportional odds model has a latent-variable derivation, slopes that survive
           any collapsing of adjacent categories, cutpoints that never need constraining, and
           one odds ratio applying at every cut (@thm-mlt-proportional-odds). The collapsing
           invariance is also how it is checked.

6. That restriction is tested by a score statistic computed from the restricted fit
           alone, on \( (c-2)|A| \) degrees of freedom, and diagnosed by plotting the slopes of
           the \( c-1 \) dichotomized fits (@prp-mlt-checks). When it fails: partial
           proportional odds, a different link, a scale effect, a different family, or an
           honest table of slopes.

7. A sequential model (@def-mlt-sequential) models the probability of stopping at
           each stage. Its likelihood factors into \( c-1 \) binomial likelihoods, so with
           stage-specific slopes it is \( c-1 \) separate binary regressions with asymptotically
           independent estimates (@thm-mlt-factorization), and with a common slope one binary
           regression on expanded data (@cor-mlt-expanded). No ordering constraint is needed.

8. With the complementary log-log link the two ordinal families coincide, cutpoints
           related by a cumulative sum on the exponential scale (@prp-mlt-cloglog): the
           grouped-data proportional hazards model. With the logit link they differ, and the
           data can tell them apart.

:::

## Notes and sources

**Coverage.** The topics follow Agresti (2015, chapter 6) and Fahrmeir, Kneib, Lang and
Marx (2021, §6.2–6.3); the treatment, proofs, examples and exercises are our own.

**Baseline-category logits.** The model grew up independently in biostatistics and
econometrics. McFadden (1974) derived it from random utility maximization, which is
the source of the name *multinomial logit* and of the discrete-choice literature; the
underlying axiom is Luce's (1959), and the independence of irrelevant alternatives in
@exr-mlt-iia is its most famous consequence and problem. McFadden's
conditional logit model, in which regressors vary by category, is the standard
extension, not treated here. Existence of the maximum likelihood estimate is settled by
Albert and Anderson (1984); @thm-bin-separation is the binary case. The individualized
fits of the remark in [Section 36.1](01-baseline-category-logits.html) are analysed by
Begg and Gray (1984).

**Cumulative models.** Modelling cumulative probabilities with a common slope goes back
to Walker and Duncan (1967). McCullagh (1980) named the proportional odds model and gave
the latent-variable derivation, the general cumulative link family and the invariance
properties of @thm-mlt-proportional-odds; McCullagh and Nelder (1989, chapter 5) is a
compressed version. Anderson (1984)
introduced a different ordinal family, the *stereotype* model, which estimates category
scores rather than assuming them. That the log-likelihood of a cumulative link model is
concave whenever the link density is log-concave was proved by Pratt (1981) and, for
grouped data specifically, by Burridge (1981); it is the reason Fisher scoring is as
docile here as it is. The asymptotics of [Section 36.2](02-cumulative-logit-models.html)
are imported rather than proved, this model not being a generalized linear model; the
standard accounts are van der Vaart (1998) and Lehmann and Casella (1998).
Agresti (2010) is the book-length treatment, and covers the scale-effect models of
[Section 36.3](03-proportional-odds.html) and the connection with rank statistics
in @exr-mlt-po-rank.

**Checking proportional odds.** The score test of @prp-mlt-checks is the one most
software reports; Brant (1990) gave the Wald version based on the separate dichotomized
fits, which is where the habit of plotting them comes from; it tests the same
hypothesis as the score statistic but is a different statistic. Partial proportional odds
models are due to Peterson and Harrell (1990), who also noted that such fits can
produce negative probabilities — the failure of @prp-mlt-ordering.

**Sequential models.** Continuation-ratio logits go back to Fienberg (1980); the
sequential formulation used here, with a link and a linear predictor at each stage, is
Tutz's (1991). Fahrmeir and Tutz (2001) and Tutz (2012) give the textbook treatments, the
latter with a taxonomy of ordinal models built out of binary ones. The factorization
of @thm-mlt-factorization is the algebra that underlies discrete-time survival
analysis, where the stages are time intervals and the continuation ratios are hazards;
the continuous-time ancestor is Cox (1972), and
@prp-mlt-cloglog, due to Läärä and Matthews (1985), connects the two.

**Data.** Both examples use public-domain data shipped with statsmodels: the travel-mode
data of @exm-mlt-mode come from the online archive of Greene's econometrics text, and
@exm-mlt-selfplacement uses an extract of the 1996 American National Election Study,
whose education and income variables are ordered scales treated here as numerical
scores. Treating them as factors changes the coefficients but no conclusion about
proportional odds. Neither example is a causal statement, for the reasons in
[Chapter 25](../ch25-causal-interpretation/index.html).

**What is missing.** Responses clustered, repeated or overdispersed relative to
@lem-mlt-multinomial need more than this chapter offers: Chapter 38 takes up
quasi-likelihood and overdispersion, Chapter 40 random effects and estimating equations,
and Chapter 39 Bayesian fitting, where the latent variable of
@thm-mlt-proportional-odds(a) becomes a data augmentation scheme.

## References

- Agresti, Alan (2010). *Analysis of Ordinal Categorical Data*. 2nd edition. Hoboken, NJ: Wiley.
- Agresti, Alan (2015). *Foundations of Linear and Generalized Linear Models*. Hoboken, NJ: Wiley.
- Albert, Adelin and Anderson, John A. (1984). On the Existence of Maximum Likelihood Estimates in Logistic Regression Models. *Biometrika* 71(1), 1–10.
- Anderson, John A. (1984). Regression and Ordered Categorical Variables. *Journal of the Royal Statistical Society, Series B* 46(1), 1–30.
- Begg, Colin B. and Gray, Robert (1984). Calculation of Polychotomous Logistic Regression Parameters Using Individualized Regressions. *Biometrika* 71(1), 11–18.
- Brant, Rollin (1990). Assessing Proportionality in the Proportional Odds Model for Ordinal Logistic Regression. *Biometrics* 46(4), 1171–1178.
- Burridge, J. (1981). A Note on Maximum Likelihood Estimation for Regression Models Using Grouped Data. *Journal of the Royal Statistical Society, Series B* 43(1), 41–45.
- Cox, David R. (1972). Regression Models and Life-Tables. *Journal of the Royal Statistical Society, Series B* 34(2), 187–220.
- Fahrmeir, Ludwig and Tutz, Gerhard (2001). *Multivariate Statistical Modelling Based on Generalized Linear Models*. 2nd edition. New York: Springer.
- Fahrmeir, Ludwig, Kneib, Thomas, Lang, Stefan and Marx, Brian D. (2021). *Regression: Models, Methods and Applications*. 2nd edition. Berlin: Springer.
- Fienberg, Stephen E. (1980). *The Analysis of Cross-Classified Categorical Data*. 2nd edition. Cambridge, MA: MIT Press.
- Läärä, Esa and Matthews, John N. S. (1985). The Equivalence of Two Models for Ordinal Data. *Biometrika* 72(1), 206–207.
- Lehmann, E. L. and Casella, George (1998). *Theory of Point Estimation*. 2nd edition. New York: Springer.
- Luce, R. Duncan (1959). *Individual Choice Behavior: A Theoretical Analysis*. New York: Wiley.
- McCullagh, Peter (1980). Regression Models for Ordinal Data. *Journal of the Royal Statistical Society, Series B* 42(2), 109–142.
- McCullagh, Peter and Nelder, John A. (1989). *Generalized Linear Models*. 2nd edition. London: Chapman and Hall.
- McFadden, Daniel (1974). Conditional Logit Analysis of Qualitative Choice Behavior. In Paul Zarembka (editor), *Frontiers in Econometrics*, 105–142. New York: Academic Press.
- Peterson, Bercedis and Harrell, Frank E. (1990). Partial Proportional Odds Models for Ordinal Response Variables. *Journal of the Royal Statistical Society, Series C* 39(2), 205–217.
- Pratt, John W. (1981). Concavity of the Log Likelihood. *Journal of the American Statistical Association* 76(373), 103–106.
- Tutz, Gerhard (1991). Sequential Models in Categorical Regression. *Computational Statistics and Data Analysis* 11(3), 275–295.
- Tutz, Gerhard (2012). *Regression for Categorical Data*. Cambridge: Cambridge University Press.
- Walker, Strother H. and Duncan, David B. (1967). Estimation of the Probability of an Event as a Function of Several Independent Variables. *Biometrika* 54(1–2), 167–179.
- van der Vaart, Aad W. (1998). *Asymptotic Statistics*. Cambridge: Cambridge University Press.
