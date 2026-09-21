# Summary and Notes

## Summary

::: {.idea}

1. Omitting regressors biases \( \hbeta_1^{\mathrm{S}} \) by \( \boldsymbol{\Pi}\bbeta_2 \) and inflates \( s^2 \), while the covariances look
   exactly as the short model claims (@thm-dep-omitted). Irrelevant regressors cost only variance (@thm-dep-overfit).

2. The short regression has the smaller mean squared error matrix iff \( \gamma\le1 \); at the design points the
   threshold is \( p_2 \) (@thm-dep-mse).

3. Under \( \Cov(\Y)=\sigma^2\V \) the usual variance estimate is off by a factor bounded by eigenvalue ratios of \( \V \) (@thm-dep-covariance), while the efficiency loss is bounded by Kantorovich's constant (@thm-dep-efficiency).

4. With non-normal errors and \( \max_ih_{ii}\to0 \), linear estimates are asymptotically normal and \( t \) and \( F \) tests
   asymptotically valid (@thm-dep-nonnormal); leverage bounds the finite-sample skewness and kurtosis (@prp-dep-cumulants).

5. Intervals for \( \sigma^2 \), variance-ratio tests and prediction intervals stay wrong in large samples (@prp-dep-sigma-interval, @prp-ci-pi-limit); quadratic balance protects \( F \) tests (@prp-dep-balance).

6. One contaminated case moves \( \hbeta \) by a leverage factor times its residual, hides itself in the fit (@thm-dep-outlier), and can move the estimate arbitrarily far (@prp-dep-breakdown).

7. With random regressors and a curved mean, least squares estimates the population projection with a sandwich
   covariance that the usual formula misses (@prp-dep-design-target, @thm-dep-random-x).

8. Collinearity breaks no assumption but inflates the variance of coefficients and of predictions against the
   pattern of the data (@prp-dep-collinear-directions).

:::

## Notes and sources

**Coverage.** The chapter covers the topics of Seber and Lee (2003, chapter 9) with its own arguments, and adds
@thm-dep-mse, @prp-dep-cumulants, the limiting coverages and the breakdown point. Measurement error is in
[Chapter 24](../ch24-errors-in-variables/index.html) (rounding in [Section 19.5](05-random-regressors.html)), perturbation theory in [Chapter 10](../ch10-computation/index.html).

**Underfitting and overfitting.** The omitted-variable formula is classical; it is @prp-lm-omitted in [Chapter 5](../ch05-model-and-least-squares/index.html). The
criterion \( \gamma\le1 \) for the restricted estimator to have the smaller mean squared error matrix is due to
Toro-Vizcarrondo and Wallace (1968) (in their notation, with the factor \( 1/2 \) in the
noncentrality, the condition reads \( \lambda\le1/2 \)), who also proposed testing it with the noncentral \( F \) distribution. Pretest estimators, which choose between the two models from the data, are taken up with model selection in [Section 29.2](../ch29-model-selection/02-cp-aic-bic.html) and [Section 29.5](../ch29-model-selection/05-selection-bias.html).

**Wrong covariance.** The bounds on the bias of the usual variance estimate are Swindel's (1968). The efficiency
bound of @thm-dep-efficiency rests on Kantorovich's inequality. The sharp bound for the determinant of the whole
covariance matrix is due to Bloomfield and Watson (1975) and Knott (1975). Grenander (1954) showed that least squares is
asymptotically efficient for polynomial trends under stationary errors, the phenomenon of @exm-dep-ar1. Kruskal's
(1968) condition for exact equality is @thm-proj-kruskal.

**Non-normality.** The Lindeberg–Feller theorem is in Billingsley (1995), and Cramér–Wold, Slutsky and continuous
mapping in van der Vaart (1998). The leverage condition \( h^*_n\to0 \) is the one used by Huber (1973) in his asymptotic
theory of robust regression. That tests on variances are sensitive to non-normality while tests on means are not was
emphasized by Box (1953). Box and Watson (1962) and Atiqullah (1962) studied the regression \( F \) test through its
first two moments, and Atiqullah introduced quadratic balance. @prp-dep-balance is a compact version of their
calculation.

**Outliers and breakdown.** The breakdown point in its finite-sample form is from Donoho and Huber (1983), and the
influence function from Hampel (1974). Regression estimators with breakdown point near one half exist, least median of
squares (Rousseeuw 1984) being the first. The stack loss data are from Brownlee (1965), and their last observation
is a much-analysed outlier. [Chapter 20](../ch20-residuals-leverage-influence/index.html) develops diagnostics for such cases.

**Random regressors.** The view of least squares as estimating a best linear approximation, with a sandwich covariance
that the usual formula gets wrong, is due to White (1980). Berkson (1950) distinguished controlled regressors from
regressors measured with error. Longley's (1967) data are the standard example of severe collinearity.

## References

- Atiqullah, M. (1962). The Estimation of Residual Variance in Quadratically Balanced Least-Squares Problems and the Robustness of the F-Test. *Biometrika* 49(1/2), 83–91.
- Berkson, Joseph (1950). Are There Two Regressions? *Journal of the American Statistical Association* 45(250), 164–180.
- Billingsley, Patrick (1995). *Probability and Measure*. 3rd edition. New York: Wiley.
- Bloomfield, Peter and Watson, Geoffrey S. (1975). The Inefficiency of Least Squares. *Biometrika* 62(1), 121–128.
- Box, George E. P. (1953). Non-Normality and Tests on Variances. *Biometrika* 40(3/4), 318–335.
- Box, George E. P. and Watson, Geoffrey S. (1962). Robustness to Non-Normality of Regression Tests. *Biometrika* 49(1/2), 93–106.
- Brownlee, K. A. (1965). *Statistical Theory and Methodology in Science and Engineering*. 2nd edition. New York: Wiley.
- Donoho, David L. and Huber, Peter J. (1983). The Notion of Breakdown Point. In P. J. Bickel, K. A. Doksum and J. L. Hodges (eds), *A Festschrift for Erich L. Lehmann*, 157–184. Belmont, CA: Wadsworth.
- Grenander, Ulf (1954). On the Estimation of Regression Coefficients in the Case of an Autocorrelated Disturbance. *The Annals of Mathematical Statistics* 25(2), 252–272.
- Hampel, Frank R. (1974). The Influence Curve and Its Role in Robust Estimation. *Journal of the American Statistical Association* 69(346), 383–393.
- Huber, Peter J. (1973). Robust Regression: Asymptotics, Conjectures and Monte Carlo. *The Annals of Statistics* 1(5), 799–821.
- Knott, M. (1975). On the Minimum Efficiency of Least Squares. *Biometrika* 62(1), 129–132.
- Kruskal, William (1968). When Are Gauss–Markov and Least Squares Estimators Identical? A Coordinate-Free Approach. *The Annals of Mathematical Statistics* 39(1), 70–75.
- Longley, James W. (1967). An Appraisal of Least Squares Programs for the Electronic Computer from the Point of View of the User. *Journal of the American Statistical Association* 62(319), 819–841.
- Rousseeuw, Peter J. (1984). Least Median of Squares Regression. *Journal of the American Statistical Association* 79(388), 871–880.
- Seber, George A. F. and Lee, Alan J. (2003). *Linear Regression Analysis*. 2nd edition. Hoboken, NJ: Wiley.
- Swindel, Benee F. (1968). On the Bias of Some Least-Squares Estimators of Variance in a General Linear Model. *Biometrika* 55(2), 313–316.
- Toro-Vizcarrondo, Carlos and Wallace, T. Dudley (1968). A Test of the Mean Square Error Criterion for Restrictions in Linear Regression. *Journal of the American Statistical Association* 63(322), 558–572.
- van der Vaart, Aad W. (1998). *Asymptotic Statistics*. Cambridge: Cambridge University Press.
- White, Halbert (1980). Using Least Squares to Approximate Unknown Regression Functions. *International Economic Review* 21(1), 149–170.
