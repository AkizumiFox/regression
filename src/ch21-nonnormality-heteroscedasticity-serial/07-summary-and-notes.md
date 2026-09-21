# Summary and Notes

## Summary

::: {.idea}

1. Residuals mix all the errors. Each differs from its error with variance \( \sigma^2h_{ii} \), and residuals keep only
   part of the errors' skewness and kurtosis, so they look more normal than the errors when \( p/n \) is not
   small (@prp-het-normal-check).

2. Normal plots of studentized residuals with simulated envelopes are the main check of normality. Coefficient
   inference is protected by the central limit theorem; prediction intervals are not, and their coverage tends to
   \( \Pr(\lvert\varepsilon_0\rvert\le z\sigma) \) (@prp-ci-pi-limit), which
   [Table 21.1.1](01-checking-normality.html#tab-het-prediction) evaluates for five error laws.

3. The Breusch–Pagan statistic is half the explained sum of squares of \( \hat{\varepsilon}_i^2/\tilde{\sigma}^2 \) on the
   variance variables, whatever the variance function. Nonnormal errors inflate its null law by \( (\kappa-1)/2 \), with \( \kappa \) the raw kurtosis, and
   Koenker's \( nR^2 \) corrects this (@thm-het-breusch-pagan, using @lem-het-clt). The Goldfeld–Quandt test is exact
   under normality (@prp-het-goldfeld-quandt).

4. Weighted least squares is optimal with known weights, unbiased with wrong ones, unbiased under symmetric errors
   with weights from residuals, and asymptotically optimal with weights from well-replicated
   groups (@thm-het-wls). Weights from a few replicates can be worse than none.

5. HC0–HC3 estimate the covariance of least squares without a variance model (@def-het-hc). HC2 is unbiased under
   homoscedasticity and HC3 is the jackknife (@prp-het-hc-facts). All are consistent (@thm-het-sandwich), but high
   leverage makes them unreliable in small samples.

6. Given \( \X \) and normal errors, the Durbin–Watson \( d\approx2(1-r_1) \) is a ratio of quadratic forms with a computable
   law, bracketed by design-free bounds when there is an intercept (@thm-het-durbin-watson,
   @lem-het-dw-eigen, @lem-het-interlace).

7. AR(1) errors are whitened by the Prais–Winsten transformation, and Cochrane–Orcutt iteration decreases the
   conditional sum of squares (@def-het-ar1, @thm-het-ar1). Smooth regressors make least squares efficient but its
   standard errors wrong. Newey–West estimates are nonnegative definite (@prp-het-newey-west). Residual
   autocorrelation often signals missing dynamics.

:::

## Notes and sources

**Sources.** The chapter covers Sen and Srivastava (1990, chapters 5–6 and sections 7.1–7.3 and 7.6), Seber and
Lee (2003, sections 10.4–10.5) and Christensen (2020, sections 12.2–12.4). It adds proofs of the facts about residuals
behind the diagnostics, the null limits of both forms of the Breusch–Pagan test, the properties of feasible
weighting in the grouped case, the consistency of the sandwich estimators, and the Durbin–Watson bounds through
interlacing. The limit theorems are from van der Vaart (1998).

**Normality.** Blom (1958) gave the plotting positions, Shapiro and Wilk (1965) and Shapiro and Francia (1972) the
two correlation-type tests, and Jarque and Bera (1980) the skewness–kurtosis score test. Lilliefors (1967)
tabulated the Kolmogorov–Smirnov test with estimated parameters, and Anderson and Darling (1952) introduced
their tail-weighted distance. Simulated envelopes are
due to Atkinson (1981). That residuals look more normal than errors is an old observation;
@prp-het-normal-check(c) quantifies it.

**Heteroscedasticity.** The score test is due to Breusch and Pagan (1979) and to Cook and Weisberg (1983); the score
principle is Rao's (1948). Koenker (1981) exposed the sensitivity to kurtosis and gave the studentized form. White
(1980) proposed the direct test and the covariance estimator, and Goldfeld and Quandt (1965) the exact \( F \) test.
Levene (1960) and Brown and Forsythe (1974) gave robust tests for groups, and Box (1953) showed the non-robustness of
normal-theory variance tests. Carroll and Ruppert (1988) is the standard monograph on weighting; Carroll (1982)
showed that weights estimated by smoothing can be fully efficient asymptotically. The symmetry argument of
@thm-het-wls(c) goes back to Kakwani (1967).

**Sandwich estimators.** Eicker (1963) gave conditions for asymptotic normality of least squares with unequal
variances and Eicker (1967) the consistency of the covariance estimator; Huber (1967) introduced the sandwich form for
general estimators, and White (1980) made it standard in econometrics. MacKinnon and White (1985) proposed HC2 and
HC3. Long and Ervin (2000) recommended HC3, and Imbens and Kolesár (2016) recommend HC2 with the degrees-of-freedom
adjustment of Bell and McCaffrey (2002).

**Serial correlation.** Durbin and Watson (1950, 1951, 1971) introduced the statistic and the bounds and studied its
exact distribution; Imhof (1961) gave the inversion formula. Durbin (1970) treated lagged responses, Breusch (1978)
and Godfrey (1978) the score test for higher-order errors, and Ljung and Box (1978) the portmanteau statistic. The
transformations are those of Cochrane and Orcutt (1949) and Prais and Winsten (1954). Grenander (1954) proved the
asymptotic efficiency of least squares for polynomial trends. Newey and West (1987) proposed the Bartlett-weighted
estimator and Andrews (1991) studied its bandwidth. Hendry and Mizon (1978) argued for testing the common factor
restriction, and Granger and Newbold (1974) exposed spurious regressions. Fuller (1996) gives the asymptotic theory
of regression with autocorrelated errors in full.

## References

- Anderson, Theodore W. and Darling, Donald A. (1952). Asymptotic Theory of Certain "Goodness of Fit" Criteria Based on Stochastic Processes. *The Annals of Mathematical Statistics* 23(2), 193–212.
- Andrews, Donald W. K. (1991). Heteroskedasticity and Autocorrelation Consistent Covariance Matrix Estimation. *Econometrica* 59(3), 817–858.
- Atkinson, Anthony C. (1981). Two Graphical Displays for Outlying and Influential Observations in Regression. *Biometrika* 68(1), 13–20.
- Bell, Robert M. and McCaffrey, Daniel F. (2002). Bias Reduction in Standard Errors for Linear Regression with Multi-Stage Samples. *Survey Methodology* 28(2), 169–181.
- Blom, Gunnar (1958). *Statistical Estimates and Transformed Beta-Variables*. New York: Wiley.
- Box, George E. P. (1953). Non-Normality and Tests on Variances. *Biometrika* 40(3/4), 318–335.
- Breusch, Trevor S. (1978). Testing for Autocorrelation in Dynamic Linear Models. *Australian Economic Papers* 17, 334–355.
- Breusch, Trevor S. and Pagan, Adrian R. (1979). A Simple Test for Heteroscedasticity and Random Coefficient Variation. *Econometrica* 47(5), 1287–1294.
- Brown, Morton B. and Forsythe, Alan B. (1974). Robust Tests for the Equality of Variances. *Journal of the American Statistical Association* 69(346), 364–367.
- Carroll, Raymond J. (1982). Adapting for Heteroscedasticity in Linear Models. *The Annals of Statistics* 10(4), 1224–1233.
- Carroll, Raymond J. and Ruppert, David (1988). *Transformation and Weighting in Regression*. New York: Chapman and Hall.
- Christensen, Ronald (2020). *Plane Answers to Complex Questions: The Theory of Linear Models*. 5th edition. Cham: Springer.
- Cochrane, Donald and Orcutt, Guy H. (1949). Application of Least Squares Regression to Relationships Containing Auto-Correlated Error Terms. *Journal of the American Statistical Association* 44(245), 32–61.
- Cook, R. Dennis and Weisberg, Sanford (1983). Diagnostics for Heteroscedasticity in Regression. *Biometrika* 70(1), 1–10.
- Durbin, James (1970). Testing for Serial Correlation in Least-Squares Regression when Some of the Regressors are Lagged Dependent Variables. *Econometrica* 38(3), 410–421.
- Durbin, James and Watson, Geoffrey S. (1950). Testing for Serial Correlation in Least Squares Regression. I. *Biometrika* 37(3/4), 409–428.
- Durbin, James and Watson, Geoffrey S. (1951). Testing for Serial Correlation in Least Squares Regression. II. *Biometrika* 38(1/2), 159–178.
- Durbin, James and Watson, Geoffrey S. (1971). Testing for Serial Correlation in Least Squares Regression. III. *Biometrika* 58(1), 1–19.
- Eicker, Friedhelm (1963). Asymptotic Normality and Consistency of the Least Squares Estimators for Families of Linear Regressions. *The Annals of Mathematical Statistics* 34(2), 447–456.
- Eicker, Friedhelm (1967). Limit Theorems for Regressions with Unequal and Dependent Errors. In *Proceedings of the Fifth Berkeley Symposium on Mathematical Statistics and Probability*, Vol. 1, 59–82. Berkeley: University of California Press.
- Fuller, Wayne A. (1996). *Introduction to Statistical Time Series*. 2nd edition. New York: Wiley.
- Godfrey, Leslie G. (1978). Testing Against General Autoregressive and Moving Average Error Models when the Regressors Include Lagged Dependent Variables. *Econometrica* 46(6), 1293–1301.
- Goldfeld, Stephen M. and Quandt, Richard E. (1965). Some Tests for Homoscedasticity. *Journal of the American Statistical Association* 60(310), 539–547.
- Granger, Clive W. J. and Newbold, Paul (1974). Spurious Regressions in Econometrics. *Journal of Econometrics* 2(2), 111–120.
- Grenander, Ulf (1954). On the Estimation of Regression Coefficients in the Case of an Autocorrelated Disturbance. *The Annals of Mathematical Statistics* 25(2), 252–272.
- Hendry, David F. and Mizon, Grayham E. (1978). Serial Correlation as a Convenient Simplification, Not a Nuisance: A Comment on a Study of the Demand for Money by the Bank of England. *The Economic Journal* 88(351), 549–563.
- Huber, Peter J. (1967). The Behavior of Maximum Likelihood Estimates under Nonstandard Conditions. In *Proceedings of the Fifth Berkeley Symposium on Mathematical Statistics and Probability*, Vol. 1, 221–233. Berkeley: University of California Press.
- Imbens, Guido W. and Kolesár, Michal (2016). Robust Standard Errors in Small Samples: Some Practical Advice. *The Review of Economics and Statistics* 98(4), 701–712.
- Imhof, J. P. (1961). Computing the Distribution of Quadratic Forms in Normal Variables. *Biometrika* 48(3/4), 419–426.
- Jarque, Carlos M. and Bera, Anil K. (1980). Efficient Tests for Normality, Homoscedasticity and Serial Independence of Regression Residuals. *Economics Letters* 6(3), 255–259.
- Kakwani, Nanak C. (1967). The Unbiasedness of Zellner's Seemingly Unrelated Regression Equations Estimators. *Journal of the American Statistical Association* 62(317), 141–142.
- Koenker, Roger (1981). A Note on Studentizing a Test for Heteroscedasticity. *Journal of Econometrics* 17(1), 107–112.
- Levene, Howard (1960). Robust Tests for Equality of Variances. In *Contributions to Probability and Statistics: Essays in Honor of Harold Hotelling*, ed. I. Olkin et al., 278–292. Stanford: Stanford University Press.
- Lilliefors, Hubert W. (1967). On the Kolmogorov–Smirnov Test for Normality with Mean and Variance Unknown. *Journal of the American Statistical Association* 62(318), 399–402.
- Ljung, Greta M. and Box, George E. P. (1978). On a Measure of Lack of Fit in Time Series Models. *Biometrika* 65(2), 297–303.
- Long, J. Scott and Ervin, Laurie H. (2000). Using Heteroscedasticity Consistent Standard Errors in the Linear Regression Model. *The American Statistician* 54(3), 217–224.
- MacKinnon, James G. and White, Halbert (1985). Some Heteroskedasticity-Consistent Covariance Matrix Estimators with Improved Finite Sample Properties. *Journal of Econometrics* 29(3), 305–325.
- Newey, Whitney K. and West, Kenneth D. (1987). A Simple, Positive Semi-Definite, Heteroskedasticity and Autocorrelation Consistent Covariance Matrix. *Econometrica* 55(3), 703–708.
- Prais, Sigbert J. and Winsten, Christopher B. (1954). Trend Estimators and Serial Correlation. Cowles Commission Discussion Paper No. 383, Chicago.
- Rao, C. Radhakrishna (1948). Large Sample Tests of Statistical Hypotheses Concerning Several Parameters with Applications to Problems of Estimation. *Proceedings of the Cambridge Philosophical Society* 44(1), 50–57.
- Seber, George A. F. and Lee, Alan J. (2003). *Linear Regression Analysis*. 2nd edition. Hoboken, NJ: Wiley.
- Sen, Ashish and Srivastava, Muni (1990). *Regression Analysis: Theory, Methods, and Applications*. New York: Springer.
- Shapiro, Samuel S. and Francia, R. S. (1972). An Approximate Analysis of Variance Test for Normality. *Journal of the American Statistical Association* 67(337), 215–216.
- Shapiro, Samuel S. and Wilk, Martin B. (1965). An Analysis of Variance Test for Normality (Complete Samples). *Biometrika* 52(3/4), 591–611.
- van der Vaart, Aad W. (1998). *Asymptotic Statistics*. Cambridge: Cambridge University Press.
- White, Halbert (1980). A Heteroskedasticity-Consistent Covariance Matrix Estimator and a Direct Test for Heteroskedasticity. *Econometrica* 48(4), 817–838.
