# Summary and Notes

## Summary

::: {.idea}

1. The training error of any fitting procedure understates its in-sample prediction error by
   \( 2n^{-1}\sum_i\Cov(\hat{\mu}_i,Y_i) \) on average. For least squares on \( p \) columns this is \( 2\sigma^2p/n \), whether
   or not the model is correct, and for a linear smoother \( 2\sigma^2\tr(\bS)/n \) (@thm-sel-optimism).

2. The degrees of freedom of a procedure are \( \sigma^{-2}\sum_i\Cov(\hat{\mu}_i,Y_i) \) (@def-sel-df). Under normality
   they equal the expected divergence of the fit, by Stein's lemma (@lem-shr-stein, @prp-sel-divergence). Soft thresholding costs its
   expected number of nonzero coefficients, hard thresholding and subset search cost more (@prp-sel-hard-df and @exm-sel-subset-df).

3. Mallows' \( C_p \) is an unbiased estimate of the scaled risk of each candidate, and it prefers a larger model iff the
   added columns have \( F>2 \) (@def-sel-cp and @prp-sel-cp). AIC estimates a Kullback–Leibler discrepancy, and AICc
   estimates it without bias in the normal linear model (@def-sel-aic-bic and @thm-sel-aicc). BIC approximates a log marginal
   likelihood, with penalty \( \log n \) per parameter.

4. Criteria with a growing penalty such as BIC find a true finite model with probability tending to one. AIC keeps each
   irrelevant column with probability about \( \Pr\{\chi^2(1)>2\} \) however large \( n \) is, but it is efficient for
   prediction when no candidate is true (@prp-sel-aic-bic).

5. In an orthonormal design selection is hard thresholding, and no threshold dominates another (@prp-sel-orthonormal and @eq-sel-hard-risk).

6. \( K \)-fold cross-validation estimates the expected error of the procedure fitted to \( n-n/K \) cases, not the error of the
   model at hand (@prp-sel-cv-target). For quadratic penalties, leave-one-out needs no refitting, and GCV replaces the
   leverages by their average (@thm-sel-loocv and @def-sel-gcv). Every data-driven step must be repeated inside the folds (@exm-sel-nested).

7. Criteria that depend on \( \text{SSE} \) agree within a size (@prp-sel-same-size). Branch and bound finds the best subset
   exactly with far fewer fits (@prp-sel-branch-bound). Stepwise steps maximize partial correlations (@prp-sel-forward-step) and can miss the best subset, and the chosen subset is unstable (@exm-sel-instability).

8. After selection, estimates are biased away from zero, minimized error estimates are optimistic, and naive intervals
   under-cover (@prp-sel-selection-bias). Sample splitting (@prp-sel-splitting), simultaneous intervals over all submodels (@prp-sel-simultaneous) and conditional intervals (@exm-sel-truncated) restore validity at different prices.

:::

## Notes and sources

**Coverage.** The chapter covers Seber and Lee (2003, chapter 12), Christensen (2020, chapter 14), Sen and Srivastava
(1990, chapter 11) and Fahrmeir, Kneib, Lang and Marx (2021, section 3.4) with its own arguments and examples. The
comparison of short and long regressions behind selection is @thm-dep-mse. Shrinkage, the high-dimensional case, and
penalized and boosted fits are in [Chapter 27](../ch27-shrinkage/index.html),
[Chapter 28](../ch28-high-dimensional/index.html) and [Chapter 30](../ch30-regularization-boosting/index.html); bagging and
random forests are outside the book's scope.

**Optimism and degrees of freedom.** The covariance form of the optimism is due to Efron (1986, 2004), and the
covariance definition of degrees of freedom, with the observation that selection adds to them, to Ye (1998). Stein
(1981) proved the identity behind @prp-sel-divergence and the unbiased risk estimate. Zou, Hastie and Tibshirani (2007) showed that the number
of nonzero lasso coefficients is unbiased for its degrees of freedom, and Tibshirani (2015) studied the search degrees
of freedom of best subset selection, of which @prp-sel-hard-df is the orthonormal case. Breiman (1992) estimated the
prediction error of a selected model by resampling residuals, and Breiman and Spector (1992) compared fixed and random
designs and leave-one-out with five- and ten-fold cross-validation.

**Criteria.** \( C_p \) and its plot are due to Mallows (1973), AIC to Akaike (1973, 1974), and AICc to Sugiura (1978) and
Hurvich and Tsai (1989), whose derivation @thm-sel-aicc follows. BIC is from Schwarz (1978); Kass and Raftery (1995) and
Kass and Wasserman (1995) discuss its accuracy. Consistency for penalties with \( c_n\to\infty \), \( c_n/n\to0 \) is due to
Nishii (1984), efficiency to Shibata (1981) and Li (1987), and their incompatibility to Yang (2005). Konishi and
Kitagawa (2008) treat information criteria at book length. The pretest estimator goes back to Bancroft (1944).

**Cross-validation.** Cross-validation is due to Stone (1974) and Geisser (1975), PRESS to Allen (1974), the equivalence
with AIC to Stone (1977), and the need for large validation sets for consistency to Shao (1993); Shao (1997) unifies these
results. GCV is from Craven and Wahba (1979) and Golub, Heath and Wahba (1979), and the one-standard-error rule from
Breiman, Friedman, Olshen and Stone (1984). Bengio and Grandvalet (2004) and Bates, Hastie and Tibshirani (2024)
analyse its variance and its target. Varma and Simon (2006) and Hastie, Tibshirani and Friedman (2009, section 7.10)
discuss cross-validation after screening.

**Search.** Enumeration by sweeps is due to Schatzoff, Tsao and Fienberg (1968), branch and bound to Furnival and Wilson
(1974), and stepwise regression to Efroymson (1960). Miller (2002) is the standard monograph on subset selection.
Natarajan (1995) proved NP-hardness, Bertsimas, King and Mazumder (2016) used mixed-integer optimization, and Hastie,
Tibshirani and Tibshirani (2020) compare best subset, forward stepwise and the lasso. Breiman (1996b) analysed
instability, Breiman (1996a) proposed bagging, Hoeting, Madigan, Raftery and Volinsky (1999) review model averaging, and
Meinshausen and Bühlmann (2010) introduced stability selection.

**Inference after selection.** The screening paradox is Freedman's (1983). Leeb and Pötscher (2005, 2006) showed that
post-selection distributions cannot be estimated uniformly. Sample splitting goes back to Cox (1975), with modern
versions by Wasserman and Roeder (2009) and Rinaldo, Wasserman and G'Sell (2019). @prp-sel-simultaneous is the PoSI
approach of Berk, Brown, Buja, Zhang and Zhao (2013). Exact conditional inference is due to Lee, Sun, Sun and Taylor
(2016) and Tibshirani, Taylor, Lockhart and Tibshirani (2016), and data carving to Fithian, Sun and Taylor (2014).

## References

- Akaike, Hirotugu (1973). Information Theory and an Extension of the Maximum Likelihood Principle. In Petrov, B. N. and Csáki, F. (eds.), *Second International Symposium on Information Theory*, 267–281. Budapest: Akadémiai Kiadó.
- Akaike, Hirotugu (1974). A New Look at the Statistical Model Identification. *IEEE Transactions on Automatic Control* 19(6), 716–723.
- Allen, David M. (1974). The Relationship Between Variable Selection and Data Augmentation and a Method for Prediction. *Technometrics* 16(1), 125–127.
- Bancroft, T. A. (1944). On Biases in Estimation Due to the Use of Preliminary Tests of Significance. *The Annals of Mathematical Statistics* 15(2), 190–204.
- Bates, Stephen, Hastie, Trevor and Tibshirani, Robert (2024). Cross-Validation: What Does It Estimate and How Well Does It Do It? *Journal of the American Statistical Association* 119(546), 1434–1445.
- Bengio, Yoshua and Grandvalet, Yves (2004). No Unbiased Estimator of the Variance of K-Fold Cross-Validation. *Journal of Machine Learning Research* 5, 1089–1105.
- Berk, Richard, Brown, Lawrence, Buja, Andreas, Zhang, Kai and Zhao, Linda (2013). Valid Post-Selection Inference. *The Annals of Statistics* 41(2), 802–837.
- Bertsimas, Dimitris, King, Angela and Mazumder, Rahul (2016). Best Subset Selection via a Modern Optimization Lens. *The Annals of Statistics* 44(2), 813–852.
- Breiman, Leo (1992). The Little Bootstrap and Other Methods for Dimensionality Selection in Regression: X-Fixed Prediction Error. *Journal of the American Statistical Association* 87(419), 738–754.
- Breiman, Leo (1996a). Bagging Predictors. *Machine Learning* 24(2), 123–140.
- Breiman, Leo (1996b). Heuristics of Instability and Stabilization in Model Selection. *The Annals of Statistics* 24(6), 2350–2383.
- Breiman, Leo, Friedman, Jerome H., Olshen, Richard A. and Stone, Charles J. (1984). *Classification and Regression Trees*. Belmont, CA: Wadsworth.
- Breiman, Leo and Spector, Philip (1992). Submodel Selection and Evaluation in Regression. The X-Random Case. *International Statistical Review* 60(3), 291–319.
- Christensen, Ronald (2020). *Plane Answers to Complex Questions: The Theory of Linear Models*. 5th edition. Cham: Springer.
- Cox, D. R. (1975). A Note on Data-Splitting for the Evaluation of Significance Levels. *Biometrika* 62(2), 441–444.
- Craven, Peter and Wahba, Grace (1979). Smoothing Noisy Data with Spline Functions. *Numerische Mathematik* 31(4), 377–403.
- Donoho, David L. and Johnstone, Iain M. (1995). Adapting to Unknown Smoothness via Wavelet Shrinkage. *Journal of the American Statistical Association* 90(432), 1200–1224.
- Efron, Bradley (1986). How Biased Is the Apparent Error Rate of a Prediction Rule? *Journal of the American Statistical Association* 81(394), 461–470.
- Efron, Bradley (2004). The Estimation of Prediction Error: Covariance Penalties and Cross-Validation. *Journal of the American Statistical Association* 99(467), 619–632.
- Efroymson, M. A. (1960). Multiple Regression Analysis. In Ralston, Anthony and Wilf, Herbert S. (eds.), *Mathematical Methods for Digital Computers*, 191–203. New York: Wiley.
- Fahrmeir, Ludwig, Kneib, Thomas, Lang, Stefan and Marx, Brian D. (2021). *Regression: Models, Methods and Applications*. 2nd edition. Berlin: Springer.
- Fithian, William, Sun, Dennis and Taylor, Jonathan (2014). Optimal Inference After Model Selection. arXiv:1410.2597.
- Freedman, David A. (1983). A Note on Screening Regression Equations. *The American Statistician* 37(2), 152–155.
- Furnival, George M. and Wilson, Robert W. (1974). Regressions by Leaps and Bounds. *Technometrics* 16(4), 499–511.
- Geisser, Seymour (1975). The Predictive Sample Reuse Method with Applications. *Journal of the American Statistical Association* 70(350), 320–328.
- Golub, Gene H., Heath, Michael and Wahba, Grace (1979). Generalized Cross-Validation as a Method for Choosing a Good Ridge Parameter. *Technometrics* 21(2), 215–223.
- Hastie, Trevor, Tibshirani, Robert and Friedman, Jerome (2009). *The Elements of Statistical Learning*. 2nd edition. New York: Springer.
- Hastie, Trevor, Tibshirani, Robert and Tibshirani, Ryan (2020). Best Subset, Forward Stepwise or Lasso? Analysis and Recommendations Based on Extensive Comparisons. *Statistical Science* 35(4), 579–592.
- Hoeting, Jennifer A., Madigan, David, Raftery, Adrian E. and Volinsky, Chris T. (1999). Bayesian Model Averaging: A Tutorial. *Statistical Science* 14(4), 382–417.
- Hurvich, Clifford M. and Tsai, Chih-Ling (1989). Regression and Time Series Model Selection in Small Samples. *Biometrika* 76(2), 297–307.
- Kass, Robert E. and Raftery, Adrian E. (1995). Bayes Factors. *Journal of the American Statistical Association* 90(430), 773–795.
- Kass, Robert E. and Wasserman, Larry (1995). A Reference Bayesian Test for Nested Hypotheses and Its Relationship to the Schwarz Criterion. *Journal of the American Statistical Association* 90(431), 928–934.
- Konishi, Sadanori and Kitagawa, Genshiro (2008). *Information Criteria and Statistical Modeling*. New York: Springer.
- Kullback, S. and Leibler, R. A. (1951). On Information and Sufficiency. *The Annals of Mathematical Statistics* 22(1), 79–86.
- Lee, Jason D., Sun, Dennis L., Sun, Yuekai and Taylor, Jonathan E. (2016). Exact Post-Selection Inference, with Application to the Lasso. *The Annals of Statistics* 44(3), 907–927.
- Leeb, Hannes and Pötscher, Benedikt M. (2005). Model Selection and Inference: Facts and Fiction. *Econometric Theory* 21(1), 21–59.
- Leeb, Hannes and Pötscher, Benedikt M. (2006). Can One Estimate the Conditional Distribution of Post-Model-Selection Estimators? *The Annals of Statistics* 34(5), 2554–2591.
- Li, Ker-Chau (1987). Asymptotic Optimality for \( C_p \), \( C_L \), Cross-Validation and Generalized Cross-Validation: Discrete Index Set. *The Annals of Statistics* 15(3), 958–975.
- Mallows, C. L. (1973). Some Comments on \( C_p \). *Technometrics* 15(4), 661–675.
- Meinshausen, Nicolai and Bühlmann, Peter (2010). Stability Selection. *Journal of the Royal Statistical Society, Series B* 72(4), 417–473.
- Miller, Alan (2002). *Subset Selection in Regression*. 2nd edition. Boca Raton, FL: Chapman and Hall/CRC.
- Natarajan, B. K. (1995). Sparse Approximate Solutions to Linear Systems. *SIAM Journal on Computing* 24(2), 227–234.
- Nishii, Ryuei (1984). Asymptotic Properties of Criteria for Selection of Variables in Multiple Regression. *The Annals of Statistics* 12(2), 758–765.
- Rinaldo, Alessandro, Wasserman, Larry and G'Sell, Max (2019). Bootstrapping and Sample Splitting for High-Dimensional, Assumption-Lean Inference. *The Annals of Statistics* 47(6), 3438–3469.
- Schatzoff, Martin, Tsao, R. and Fienberg, Stephen (1968). Efficient Calculation of All Possible Regressions. *Technometrics* 10(4), 769–779.
- Schwarz, Gideon (1978). Estimating the Dimension of a Model. *The Annals of Statistics* 6(2), 461–464.
- Seber, George A. F. and Lee, Alan J. (2003). *Linear Regression Analysis*. 2nd edition. Hoboken, NJ: Wiley.
- Sen, Ashish and Srivastava, Muni (1990). *Regression Analysis: Theory, Methods, and Applications*. New York: Springer.
- Shao, Jun (1993). Linear Model Selection by Cross-Validation. *Journal of the American Statistical Association* 88(422), 486–494.
- Shao, Jun (1997). An Asymptotic Theory for Linear Model Selection. *Statistica Sinica* 7(2), 221–264.
- Shibata, Ritei (1981). An Optimal Selection of Regression Variables. *Biometrika* 68(1), 45–54.
- Stein, Charles M. (1981). Estimation of the Mean of a Multivariate Normal Distribution. *The Annals of Statistics* 9(6), 1135–1151.
- Stone, M. (1974). Cross-Validatory Choice and Assessment of Statistical Predictions. *Journal of the Royal Statistical Society, Series B* 36(2), 111–147.
- Stone, M. (1977). An Asymptotic Equivalence of Choice of Model by Cross-Validation and Akaike's Criterion. *Journal of the Royal Statistical Society, Series B* 39(1), 44–47.
- Sugiura, Nariaki (1978). Further Analysts of the Data by Akaike's Information Criterion and the Finite Corrections. *Communications in Statistics: Theory and Methods* 7(1), 13–26.
- Tibshirani, Ryan J. (2015). Degrees of Freedom and Model Search. *Statistica Sinica* 25(3), 1265–1296.
- Tibshirani, Ryan J., Taylor, Jonathan, Lockhart, Richard and Tibshirani, Robert (2016). Exact Post-Selection Inference for Sequential Regression Procedures. *Journal of the American Statistical Association* 111(514), 600–620.
- Tierney, Luke and Kadane, Joseph B. (1986). Accurate Approximations for Posterior Moments and Marginal Densities. *Journal of the American Statistical Association* 81(393), 82–86.
- Varma, Sudhir and Simon, Richard (2006). Bias in Error Estimation When Using Cross-Validation for Model Selection. *BMC Bioinformatics* 7, 91.
- Wasserman, Larry and Roeder, Kathryn (2009). High-Dimensional Variable Selection. *The Annals of Statistics* 37(5A), 2178–2201.
- Yang, Yuhong (2005). Can the Strengths of AIC and BIC Be Shared? A Conflict Between Model Identification and Regression Estimation. *Biometrika* 92(4), 937–950.
- Ye, Jianming (1998). On Measuring and Correcting the Effects of Data Mining and Model Selection. *Journal of the American Statistical Association* 93(441), 120–131.
- Zou, Hui, Hastie, Trevor and Tibshirani, Robert (2007). On the "Degrees of Freedom" of the Lasso. *The Annals of Statistics* 35(5), 2173–2192.
