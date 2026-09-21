# Summary and Notes

## Summary

::: {.idea}

1. A **variance specification** (@def-ql-variance-function) gives the mean and the variance
   of each response and nothing else: \( \E Y_i=h(\x_{(i)}\T\bbeta) \) and
   \( \Var(Y_i)=\phi V(\mu_i)/w_i \), with \( \phi \) free even for a variance function
   borrowed from a family that fixed it. No distribution need exist with these moments, and
   for ungrouped binary data none does with \( \phi\ne1 \). Overdispersion of a group
   proportion has two shapes, a constant multiple of the binomial variance and the factor
   \( 1+\rho(m-1) \) of exchangeable trials (@lem-ql-exchangeable), agreeing at one group
   size only.

2. A variance function fixes the estimating equations, not the fit: the NB1 model shares
   the quasi-Poisson variance function and gives a different
   \( \hbeta \) (@exm-ql-strikes-choice).

3. The **quasi-likelihood** \( Q(\mu;y)=w\phi^{-1}\int_y^{\mu}(y-t)V(t)^{-1}dt \) of
   @def-ql-quasi has the score, the two Bartlett identities and the saturated maximum of a
   log-likelihood, and equals the log-likelihood minus its saturated value for an
   exponential dispersion family (@prp-ql-properties); its \( -2\phi \) multiple is the
   **quasi-deviance**.

4. The **quasi-score** is the score @eq-glm-score, so the estimate is the usual iteratively
   reweighted least squares fit (@thm-ql-score(a)), unbiased as soon as the mean model is
   right whatever the variances; the identity between its variance and its sensitivity
   needs the variance model too (@thm-ql-score(b)–(e)).

5. Among **linear estimating functions** \( \mathbf{C}(\bbeta)(\y-\bmu) \) the
   quasi-score maximizes the Godambe information, and nothing else does except its
   nonsingular linear transformations (@prp-ql-optimality). The proof is Gauss–Markov with
   an orthogonal projection in place of the completion of the square.

6. Under (Q1)–(Q5), \( \hbeta \) is asymptotically normal with the **sandwich** covariance
   \( \A_n^{-1}\B_n\A_n^{-1} \), whose middle is estimated from the
   residuals (@thm-ql-asymptotics); a correct variance specification makes
   \( \B_n=\A_n \) and collapses it to the model-based \( \phi(\X\T\W\X)^{-1} \), then
   valid and more efficient.

7. The dispersion is estimated by \( X^2/(n-p) \), of expectation exactly \( \phi \) at
   the true mean, while the deviance version is not in general
   consistent (@prp-ql-dispersion(a),(b)). Nested models are compared by the
   quasi-\( F \) statistic @eq-ql-quasi-f, which needs a correct variance function; the
   Akaike criterion has no analogue here beyond a convention (@prp-ql-dispersion(c),(d)).

8. Overdispersion comes from unobserved heterogeneity, clustering and contagion, each
   computable exactly (@prp-ql-overdispersion(a)–(c)). Ignoring a dispersion \( \phi_0 \)
   turns a nominal \( 5\% \) test of one coefficient into one of size
   \( \Pr\{\chi^2(1)>3.841/\phi_0\} \), and of six coefficients into something far
   worse (@eq-ql-size). The repairs assume different things — a mean and a variance shape,
   a full distribution, a mixing law — and the last changes the estimand.

:::

## Notes and sources

**Wedderburn's paper.** Quasi-likelihood is due to Wedderburn (1974), written two years
after the generalized linear model itself (Nelder and Wedderburn 1972) and in the same
spirit: notice which part of the structure the algorithm uses, and keep only that. The
definition @eq-ql-quasi, the two Bartlett identities and the table of quasi-likelihoods are
all in that paper, as is the name. McCullagh and Nelder (1989, chapter 9) is the standard
extended treatment; their chapter 10 develops the joint model for the mean and the
dispersion, which this chapter does not cover. Multiplying Poisson or binomial standard
errors by \( \sqrt{X^2/(n-p)} \) predates the theory — Finney (1947) used it for probit
analysis — and Wedderburn said what it estimates.

**Estimating functions and optimality.** The general theory of unbiased estimating
functions is Godambe (1960), who showed the likelihood score optimal among them in the
sense used in @prp-ql-optimality; Godambe and Heyde (1987) survey the subject. The
projection proof given here makes the analogy exact: "optimal among linear estimating
functions" stands to "optimal among all" as "best linear unbiased" stands to "minimum
variance unbiased" in @thm-opt-gauss-markov. McCullagh (1983) proved the asymptotics under
conditions close to (Q1)–(Q5).

**The sandwich.** Its history is traced in
[Section 21.4](../ch21-nonnormality-heteroscedasticity-serial/04-sandwich-estimators.html):
Eicker (1963), Huber (1967) and White (1980) in three literatures at once. The
generalized-linear-model version, with the working weights in the bread, is in
Royall (1986), who argued for it as the default; Liang and Zeger (1986) made it the
foundation of generalized estimating equations. Its small-sample bias is analysed by
Kauermann and Carroll (2001), and the leverage corrections of @exm-ql-coverage are the HC2
and HC3 of MacKinnon and White (1985).

**Binomial overdispersion.** The exchangeable-trials calculation @eq-ql-exchangeable is the
design effect of cluster sampling (Kish 1965), which
[Chapter 33](../ch33-clustered-longitudinal-splitplot/index.html) met as
\( 1+(m-1)\rho \). The beta-binomial distribution is due to Skellam (1948), and
Williams (1982) proposed fitting the variance function
\( \{1+\rho(m_i-1)\}\mu(1-\mu)/m_i \) by alternating a quasi-likelihood step with a
moment equation setting \( X^2 \) equal to its degrees of freedom.

**Which device to use.** Ver Hoef and Boveng (2007) is a clear applied statement of the
quasi-Poisson against negative binomial question, showing that the choice is between the
NB1 and NB2 variance functions and matters most where the fitted means are most spread
out — what @exm-ql-strikes-choice shows in miniature. Cameron and Trivedi (2013) covers the
quasi-likelihood, negative binomial and sandwich routes together. The "QAIC" convention
of @prp-ql-dispersion(d) is from Lebreton, Burnham, Clobert and Anderson (1992), widely
used in ecology and resting on no theorem.

**The data.** The strike durations of @exm-ql-strikes-fit are the subset of Kennan's (1985)
contract-strike series distributed with statsmodels, from Bureau of Labor Statistics
records. Kennan fitted duration models of a different kind; the Poisson mean model here
illustrates the methods of this chapter and is not an account of strike durations.

**What comes next.** [Chapter 39](../ch39-glms-in-practice-bayes/index.html) returns to models with likelihoods, taking up model
building, diagnostics and the Bayesian treatment. [Chapter 40](../ch40-glmm-gee/index.html) extends the estimating-equation
idea to correlated responses, where the working variance becomes a working *correlation
matrix*. [Chapter 41](../ch41-missing-data/index.html) asks what happens when some responses are not observed at all.

## References

- Cameron, A. Colin and Trivedi, Pravin K. (2013). *Regression Analysis of Count Data*. 2nd edition. Cambridge: Cambridge University Press.
- Eicker, Friedhelm (1963). Asymptotic Normality and Consistency of the Least Squares Estimators for Families of Linear Regressions. *The Annals of Mathematical Statistics* 34(2), 447–456.
- Fahrmeir, Ludwig and Kaufmann, Heinz (1985). Consistency and Asymptotic Normality of the Maximum Likelihood Estimator in Generalized Linear Models. *The Annals of Statistics* 13(1), 342–368.
- Finney, David J. (1947). *Probit Analysis: A Statistical Treatment of the Sigmoid Response Curve*. Cambridge: Cambridge University Press.
- Godambe, V. P. (1960). An Optimum Property of Regular Maximum Likelihood Estimation. *The Annals of Mathematical Statistics* 31(4), 1208–1211.
- Godambe, V. P. and Heyde, C. C. (1987). Quasi-Likelihood and Optimal Estimation. *International Statistical Review* 55(3), 231–244.
- Huber, Peter J. (1967). The Behavior of Maximum Likelihood Estimates under Nonstandard Conditions. In *Proceedings of the Fifth Berkeley Symposium on Mathematical Statistics and Probability*, volume 1, 221–233. Berkeley: University of California Press.
- Kauermann, Göran and Carroll, Raymond J. (2001). A Note on the Efficiency of Sandwich Covariance Matrix Estimation. *Journal of the American Statistical Association* 96(456), 1387–1396.
- Kennan, John (1985). The Duration of Contract Strikes in US Manufacturing. *Journal of Econometrics* 28(1), 5–28.
- Kish, Leslie (1965). *Survey Sampling*. New York: Wiley.
- Lebreton, Jean-Dominique, Burnham, Kenneth P., Clobert, Jean and Anderson, David R. (1992). Modeling Survival and Testing Biological Hypotheses Using Marked Animals: A Unified Approach with Case Studies. *Ecological Monographs* 62(1), 67–118.
- Liang, Kung-Yee and Zeger, Scott L. (1986). Longitudinal Data Analysis Using Generalized Linear Models. *Biometrika* 73(1), 13–22.
- MacKinnon, James G. and White, Halbert (1985). Some Heteroskedasticity-Consistent Covariance Matrix Estimators with Improved Finite Sample Properties. *Journal of Econometrics* 29(3), 305–325.
- McCullagh, Peter (1983). Quasi-Likelihood Functions. *The Annals of Statistics* 11(1), 59–67.
- McCullagh, Peter and Nelder, John A. (1989). *Generalized Linear Models*. 2nd edition. London: Chapman and Hall.
- Nelder, John A. and Wedderburn, Robert W. M. (1972). Generalized Linear Models. *Journal of the Royal Statistical Society, Series A* 135(3), 370–384.
- Royall, Richard M. (1986). Model Robust Confidence Intervals Using Maximum Likelihood Estimators. *International Statistical Review* 54(2), 221–226.
- Skellam, J. G. (1948). A Probability Distribution Derived from the Binomial Distribution by Regarding the Probability of Success as Variable between the Sets of Trials. *Journal of the Royal Statistical Society, Series B* 10(2), 257–261.
- Ver Hoef, Jay M. and Boveng, Peter L. (2007). Quasi-Poisson vs. Negative Binomial Regression: How Should We Model Overdispersed Count Data? *Ecology* 88(11), 2766–2772.
- Wedderburn, Robert W. M. (1974). Quasi-Likelihood Functions, Generalized Linear Models, and the Gauss–Newton Method. *Biometrika* 61(3), 439–447.
- White, Halbert (1980). A Heteroskedasticity-Consistent Covariance Matrix Estimator and a Direct Test for Heteroskedasticity. *Econometrica* 48(4), 817–838.
- Williams, D. A. (1982). Extra-Binomial Variation in Logistic Linear Models. *Journal of the Royal Statistical Society, Series C (Applied Statistics)* 31(2), 144–148.
