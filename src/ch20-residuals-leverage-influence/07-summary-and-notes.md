# Summary and Notes

## Summary

::: {.idea}

1. Residuals are the projection of the errors onto \( \C(\X)\perpc \), with variances \( \sigma^2(1-h_{ii}) \). The internally
   studentized residual \( r_i \) has a scaled beta law and is bounded by \( \sqrt{n-p} \); the externally studentized residual
   \( t_i \) uses a scale estimated without case \( i \) (@def-res-residuals and @prp-res-internal).

2. Residual plots show only the part of a misspecified mean that is orthogonal to \( \C(\X) \): they can reveal wrong fitted
   values, not biased coefficients. The added-variable and component-plus-residual plots share their line and their residuals
   and differ in the horizontal axis (@prp-res-plots and @prp-res-partial-residual).

3. Leverage is potential influence. Under a Gaussian design \( h_{ii}-1/n \) is a scaled beta variable, so the \( 2p/n \) rule fires at
   a rate that depends on \( p \). With finite second moments the largest leverage vanishes as \( n \) grows; with heavy-tailed
   regressors it need not (@prp-res-leverage-beta and @prp-res-max-leverage).

4. Partial leverage isolates the pull of a case on one coefficient, and \( h_{ii} \) is the leverage without \( \x_j \) plus the partial
   leverage for \( \x_j \) (@def-res-partial-leverage and @prp-res-partial-leverage).

5. Deleting one case changes the coefficients, fitted values and residual variance by explicit formulas in \( \hat{\varepsilon}_i \)
   and \( h_{ii} \), and PRESS overestimates fixed-design prediction error slightly (@thm-res-deletion and @prp-res-press).

6. \( t_i \) is the \( t \) statistic of the mean-shift outlier model and has the \( t(n-p-1,\Delta\sqrt{1-h_{ii}}/\sigma) \) law; the maximum
   over cases, with a Bonferroni correction, is an outlier test of level at most \( \alpha \) (@thm-res-external-t).

7. Cook's distance, DFFITS and DFBETAS are products of a residual and a leverage factor. DFFITS is the largest standardized change
   in any linear combination of the coefficients. Their conventional cut-offs are scaled conventions, not tests
   (@def-res-cooks, @thm-res-influence and @prp-res-cook-null).

8. A set of cases can be deleted exactly, and the mean-shift \( F \) test checks a prespecified set (@thm-res-group-deletion).

9. A cluster of \( m \) similar bad cases masks itself: each has leverage below \( 1/m \), its residual shrinks by \( 1+mh_\star \), single
   deletion recovers a fraction \( 1/[m(1+(m-1)h_\star)] \) of the group effect, and clean cases are swamped (@prp-res-masking).

10. Huber and bisquare M-estimates are computed by IRLS, which never increases the objective with a fixed scale. They resist
    outlying responses at cases of moderate leverage, but not bad leverage points. High-breakdown fits such as LTS, computed by concentration steps, expose
    clusters that least squares and M-estimates absorb (@prp-res-irls, @prp-dep-breakdown and @prp-res-cstep).

:::

## Notes and sources

**Sources.** Seber and Lee (2003, sections 10.2 and 10.6) derive the studentized residuals and their laws, the
leave-one-out diagnostics including COVRATIO and the Andrews–Pregibon statistic, the mean-shift test for a set of cases,
and brief accounts of masking, swamping and robust residuals. Christensen (2020, sections 12.1 and 12.5–12.6) develops
leverage from the projection \( \M \) and derives the updating formulae and predicted residuals. Sen and Srivastava (1990,
chapter 8) and Rencher and Schaalje (2008, chapter 9) are briefer. The proofs here
take the book's own route: the mean-shift model and the Frisch–Waugh–Lovell theorem for the law of \( t_i \) and for group
deletion, orthogonal projections for the law of \( r_i \), and partial leverage for DFBETAS.

**Residuals and influence.** Hoaglin and Welsch (1978) proposed the \( 2p/n \) rule. The \( t \) law of the externally
studentized residual is due to Beckman and Trussell (1974). Cook (1977) introduced his distance and its ellipsoid reading,
Belsley, Kuh and Welsch (1980) DFFITS, DFBETAS and COVRATIO with their size-adjusted cut-offs, and Andrews and Pregibon (1978)
a determinant-based measure for sets of cases. Book-length treatments are Cook and Weisberg (1982), Atkinson (1985) and
Chatterjee and Hadi (1988); the review of Beckman and Cook (1983) surveys outlier testing, including the Bonferroni approach.
PRESS is from Allen (1974), and partial residual plots were popularized by Larsen and McCleary (1972). @prp-res-leverage-beta
is the regression form of the classical beta law for the Mahalanobis distance of a normal observation from the sample mean.
Kiefer and Wolfowitz (1960) founded optimal design, which this book does not develop.

**Masking and robust fitting.** Beckman and Cook (1983) discuss masking and swamping at length; @prp-res-masking is a formal
version for the simplest configuration. Hadi and Simonoff (1993), Atkinson and Riani (2000) and Rousseeuw and van Zomeren
(1990) give the procedures mentioned in [Section 20.5](05-masking.html). M-estimation began with Huber (1964); Huber (1973)
treated regression, under the leverage condition of @prp-res-max-leverage. The bisquare is from Beaton and Tukey (1974), IRLS
for robust regression was studied by Holland and Welsch (1977), and the finite-sample breakdown point is from Donoho and Huber
(1983). LMS and LTS are from Rousseeuw (1984), with Rousseeuw and Leroy (1987) the standard account and Rousseeuw and Van
Driessen (2006) the concentration step. The MM-estimator is Yohai's (1987). Maronna, Martin and Yohai (2006) give the theory; median and quantile regression are the subject of [Chapter 45](../ch45-quantile-gamlss/index.html) (@def-qnt-quantile).

**Data.** Engel's budgets were made widely available by Koenker and Bassett (1982) and are distributed with statsmodels; the
state data are those of [Chapter 6](../ch06-projections/index.html).

## References

- Allen, David M. (1974). The Relationship Between Variable Selection and Data Augmentation and a Method for Prediction. *Technometrics* 16(1), 125–127.
- Andrews, David F. and Pregibon, Daryl (1978). Finding the Outliers that Matter. *Journal of the Royal Statistical Society, Series B* 40(1), 85–93.
- Atkinson, Anthony C. (1985). *Plots, Transformations, and Regression: An Introduction to Graphical Methods of Diagnostic Regression Analysis*. Oxford: Clarendon Press.
- Atkinson, Anthony C. and Riani, Marco (2000). *Robust Diagnostic Regression Analysis*. New York: Springer.
- Beaton, Albert E. and Tukey, John W. (1974). The Fitting of Power Series, Meaning Polynomials, Illustrated on Band-Spectroscopic Data. *Technometrics* 16(2), 147–185.
- Beckman, Richard J. and Cook, R. Dennis (1983). Outlier..........s. *Technometrics* 25(2), 119–149.
- Beckman, Richard J. and Trussell, H. Joel (1974). The Distribution of an Arbitrary Studentized Residual and the Effects of Updating in Multiple Regression. *Journal of the American Statistical Association* 69(345), 199–201.
- Belsley, David A., Kuh, Edwin and Welsch, Roy E. (1980). *Regression Diagnostics: Identifying Influential Data and Sources of Collinearity*. New York: Wiley.
- Chatterjee, Samprit and Hadi, Ali S. (1988). *Sensitivity Analysis in Linear Regression*. New York: Wiley.
- Christensen, Ronald (2020). *Plane Answers to Complex Questions: The Theory of Linear Models*. 5th edition. Cham: Springer.
- Cook, R. Dennis (1977). Detection of Influential Observation in Linear Regression. *Technometrics* 19(1), 15–18.
- Cook, R. Dennis and Weisberg, Sanford (1982). *Residuals and Influence in Regression*. New York: Chapman and Hall.
- Donoho, David L. and Huber, Peter J. (1983). The Notion of Breakdown Point. In P. J. Bickel, K. A. Doksum and J. L. Hodges (eds.), *A Festschrift for Erich L. Lehmann*, 157–184. Belmont, CA: Wadsworth.
- Hadi, Ali S. and Simonoff, Jeffrey S. (1993). Procedures for the Identification of Multiple Outliers in Linear Models. *Journal of the American Statistical Association* 88(424), 1264–1272.
- Hoaglin, David C. and Welsch, Roy E. (1978). The Hat Matrix in Regression and ANOVA. *The American Statistician* 32(1), 17–22.
- Holland, Paul W. and Welsch, Roy E. (1977). Robust Regression Using Iteratively Reweighted Least-Squares. *Communications in Statistics – Theory and Methods* 6(9), 813–827.
- Huber, Peter J. (1964). Robust Estimation of a Location Parameter. *The Annals of Mathematical Statistics* 35(1), 73–101.
- Huber, Peter J. (1973). Robust Regression: Asymptotics, Conjectures and Monte Carlo. *The Annals of Statistics* 1(5), 799–821.
- Kiefer, Jack and Wolfowitz, Jacob (1960). The Equivalence of Two Extremum Problems. *Canadian Journal of Mathematics* 12, 363–366.
- Koenker, Roger and Bassett, Gilbert (1982). Robust Tests for Heteroscedasticity Based on Regression Quantiles. *Econometrica* 50(1), 43–61.
- Larsen, Wayne A. and McCleary, Susan J. (1972). The Use of Partial Residual Plots in Regression Analysis. *Technometrics* 14(3), 781–790.
- Maronna, Ricardo A., Martin, R. Douglas and Yohai, Víctor J. (2006). *Robust Statistics: Theory and Methods*. Chichester: Wiley.
- Rencher, Alvin C. and Schaalje, G. Bruce (2008). *Linear Models in Statistics*. 2nd edition. Hoboken, NJ: Wiley.
- Rousseeuw, Peter J. (1984). Least Median of Squares Regression. *Journal of the American Statistical Association* 79(388), 871–880.
- Rousseeuw, Peter J. and Leroy, Annick M. (1987). *Robust Regression and Outlier Detection*. New York: Wiley.
- Rousseeuw, Peter J. and Van Driessen, Katrien (2006). Computing LTS Regression for Large Data Sets. *Data Mining and Knowledge Discovery* 12(1), 29–45.
- Rousseeuw, Peter J. and van Zomeren, Bert C. (1990). Unmasking Multivariate Outliers and Leverage Points. *Journal of the American Statistical Association* 85(411), 633–639.
- Seber, George A. F. and Lee, Alan J. (2003). *Linear Regression Analysis*. 2nd edition. Hoboken, NJ: Wiley.
- Sen, Ashish and Srivastava, Muni (1990). *Regression Analysis: Theory, Methods, and Applications*. New York: Springer.
- Yohai, Víctor J. (1987). High Breakdown-Point and High Efficiency Robust Estimates for Regression. *The Annals of Statistics* 15(2), 642–656.
