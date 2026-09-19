# Summary and Notes

## Summary

::: {.idea}

1. When \( p>n \), least squares interpolates and no single coefficient is estimable. The minimum-norm
   solution \( \X^+\y \) is the ridgeless limit of ridge and the limit of gradient descent from zero; its error
   splits into the invisible null-space part of \( \bbeta \) and the variance \( \sigma^2\sum_id_i^{-2} \) (@thm-hd-min-norm).

2. For isotropic Gaussian regressors its risk is exact (@prp-hd-isotropic-risk). Fitting \( p \) of \( D \)
   features gives the double descent curve (@cor-hd-double-descent), whose peak at \( p=n \) belongs to the
   estimator and disappears with a tuned ridge penalty.

3. Sparsity restores identifiability (@prp-hd-sparse-identifiable). The compatibility and restricted
   eigenvalue conditions ask for curvature only on a cone of approximately sparse vectors (@def-hd-re); they
   transfer from populations to samples with a loss of \( 16s\delta \) (@prp-hd-re-basic), and can be small in
   samples even when they are large in the population.

4. On the event \( \norm{\X\T\be/n}_\infty\le\lambda_0=\sigma\sqrt{2\log(2p/\delta)/n} \), which has probability at
   least \( 1-\delta \) (@lem-hd-gaussian-max), every lasso estimate with \( \lambda\ge2\lambda_0 \) has its error in the
   cone and prediction error at most \( 3\lambda\norm{\bbeta}_1 \) without conditions and \( 9\lambda^2s/\phi^2(S) \)
   under compatibility (@thm-hd-lasso-prediction). Its \( \ell_1 \) and \( \ell_2 \) errors are of order
   \( \lambda s \) and \( \lambda\sqrt s \) (@thm-hd-lasso-estimation).

5. Exact sign recovery holds under the irrepresentable and beta-min conditions (@thm-hd-support). Without
   the irrepresentable condition it has probability at most \( \frac12 \) for every sample
   size (@prp-hd-irrepresentable-necessary); screening needs only the restricted eigenvalue and beta-min
   conditions (@prp-hd-screening).

6. Intervals computed after selection under-cover. The debiased lasso adds a one-step correction; its error
   is exactly normal plus a remainder bounded by \( \sqrt n\mu_j\norm{\hat{\bbeta}_\lambda-\bbeta}_1 \), which vanishes when
   \( s\log p/\sqrt n\to0 \) (@thm-hd-debiased). With node-wise lasso rows it is the Frisch–Waugh–Lovell estimator
   with a lasso residual (@lem-hd-nodewise).

:::

## Notes and sources

**Scope.** None of the planning sources treats \( p>n \) in depth; Agresti (2015, section 11.2) notes it briefly,
through a Bayesian prior. Bühlmann and van de Geer (2011), Hastie, Tibshirani and Wainwright (2015) and
Wainwright (2019) develop all the topics of the chapter at length.

**Interpolation.** The term *double descent* is from Belkin, Hsu, Ma and Mandal (2019). Belkin, Hsu and Xu (2020)
give exact risk formulas for Gaussian features, including omitted features as in @cor-hd-double-descent.
Hastie, Montanari, Rosset and Tibshirani (2022) treat general covariances and ridge asymptotically, and
Bartlett, Long, Lugosi and Tsigler (2020) characterize benign overfitting.

**Design conditions and the lasso.** The lasso is due to Tibshirani (1996); Tibshirani (2013) studies uniqueness.
Greenshtein and Ritov (2004) introduced persistence, behind the slow rate. The restricted eigenvalue condition
and the fast rates are from Bickel, Ritov and Tsybakov (2009), with the Dantzig selector of Candès and Tao (2007);
the basic-inequality presentation follows Bühlmann and van de Geer (2011, chapter 6), and van de Geer and
Bühlmann (2009) compare the conditions. Donoho and Huo (2001) used mutual coherence, and Cohen, Dahmen and
DeVore (2009) the null space property. Raskutti, Wainwright and Yu (2010, 2011) prove the random-design bound
and the minimax lower bound. Belloni, Chernozhukov and Wang (2011) and Sun and Zhang (2012) free \( \lambda \) from \( \sigma \).

**Selection and inference.** The irrepresentable condition is due to Zhao and Yu (2006) and Meinshausen and
Bühlmann (2006), the primal–dual witness proof to Wainwright (2009); Zou (2006) and Meinshausen and Bühlmann (2010)
propose the adaptive lasso and stability selection. Leeb and Pötscher (2005, 2006) analyse post-selection
distributions and prove that they cannot be estimated uniformly. The debiased lasso is due to Zhang and Zhang (2014), van de Geer, Bühlmann, Ritov and Dezeure
(2014) and Javanmard and Montanari (2014). The alternatives are due to Wasserman and Roeder (2009), Meinshausen,
Meier and Bühlmann (2009), Lee, Sun, Sun and Taylor (2016), Berk, Brown, Buja, Zhang and Zhao (2013) and Belloni,
Chernozhukov and Hansen (2014). [Chapter 29](../ch29-model-selection/index.html) treats selection bias for classical
rules, and [Chapter 30](../ch30-regularization-boosting/index.html) extends the lasso to other penalties.

## References

- Agresti, Alan (2015). *Foundations of Linear and Generalized Linear Models*. Hoboken, NJ: Wiley.
- Bartlett, Peter L., Long, Philip M., Lugosi, Gábor and Tsigler, Alexander (2020). Benign Overfitting in Linear Regression. *Proceedings of the National Academy of Sciences* 117(48), 30063–30070.
- Belkin, Mikhail, Hsu, Daniel, Ma, Siyuan and Mandal, Soumik (2019). Reconciling Modern Machine-Learning Practice and the Classical Bias–Variance Trade-off. *Proceedings of the National Academy of Sciences* 116(32), 15849–15854.
- Belkin, Mikhail, Hsu, Daniel and Xu, Ji (2020). Two Models of Double Descent for Weak Features. *SIAM Journal on Mathematics of Data Science* 2(4), 1167–1180.
- Belloni, Alexandre, Chernozhukov, Victor and Hansen, Christian (2014). Inference on Treatment Effects after Selection among High-Dimensional Controls. *The Review of Economic Studies* 81(2), 608–650.
- Belloni, Alexandre, Chernozhukov, Victor and Wang, Lie (2011). Square-Root Lasso: Pivotal Recovery of Sparse Signals via Conic Programming. *Biometrika* 98(4), 791–806.
- Berk, Richard, Brown, Lawrence, Buja, Andreas, Zhang, Kai and Zhao, Linda (2013). Valid Post-Selection Inference. *The Annals of Statistics* 41(2), 802–837.
- Bickel, Peter J., Ritov, Ya'acov and Tsybakov, Alexandre B. (2009). Simultaneous Analysis of Lasso and Dantzig Selector. *The Annals of Statistics* 37(4), 1705–1732.
- Bühlmann, Peter and van de Geer, Sara (2011). *Statistics for High-Dimensional Data: Methods, Theory and Applications*. Berlin: Springer.
- Candès, Emmanuel and Tao, Terence (2007). The Dantzig Selector: Statistical Estimation When p Is Much Larger than n. *The Annals of Statistics* 35(6), 2313–2351.
- Cohen, Albert, Dahmen, Wolfgang and DeVore, Ronald (2009). Compressed Sensing and Best k-Term Approximation. *Journal of the American Mathematical Society* 22(1), 211–231.
- Donoho, David L. and Huo, Xiaoming (2001). Uncertainty Principles and Ideal Atomic Decomposition. *IEEE Transactions on Information Theory* 47(7), 2845–2862.
- Greenshtein, Eitan and Ritov, Ya'acov (2004). Persistence in High-Dimensional Linear Predictor Selection and the Virtue of Overparametrization. *Bernoulli* 10(6), 971–988.
- Hastie, Trevor, Montanari, Andrea, Rosset, Saharon and Tibshirani, Ryan J. (2022). Surprises in High-Dimensional Ridgeless Least Squares Interpolation. *The Annals of Statistics* 50(2), 949–986.
- Hastie, Trevor, Tibshirani, Robert and Wainwright, Martin (2015). *Statistical Learning with Sparsity: The Lasso and Generalizations*. Boca Raton, FL: CRC Press.
- Javanmard, Adel and Montanari, Andrea (2014). Confidence Intervals and Hypothesis Testing for High-Dimensional Regression. *Journal of Machine Learning Research* 15, 2869–2909.
- Lee, Jason D., Sun, Dennis L., Sun, Yuekai and Taylor, Jonathan E. (2016). Exact Post-Selection Inference, with Application to the Lasso. *The Annals of Statistics* 44(3), 907–927.
- Leeb, Hannes and Pötscher, Benedikt M. (2005). Model Selection and Inference: Facts and Fiction. *Econometric Theory* 21(1), 21–59.
- Leeb, Hannes and Pötscher, Benedikt M. (2006). Can One Estimate the Conditional Distribution of Post-Model-Selection Estimators? *The Annals of Statistics* 34(5), 2554–2591.
- Meinshausen, Nicolai and Bühlmann, Peter (2006). High-Dimensional Graphs and Variable Selection with the Lasso. *The Annals of Statistics* 34(3), 1436–1462.
- Meinshausen, Nicolai and Bühlmann, Peter (2010). Stability Selection. *Journal of the Royal Statistical Society, Series B* 72(4), 417–473.
- Meinshausen, Nicolai, Meier, Lukas and Bühlmann, Peter (2009). p-Values for High-Dimensional Regression. *Journal of the American Statistical Association* 104(488), 1671–1681.
- Raskutti, Garvesh, Wainwright, Martin J. and Yu, Bin (2010). Restricted Eigenvalue Properties for Correlated Gaussian Designs. *Journal of Machine Learning Research* 11, 2241–2259.
- Raskutti, Garvesh, Wainwright, Martin J. and Yu, Bin (2011). Minimax Rates of Estimation for High-Dimensional Linear Regression over \( \ell_q \)-Balls. *IEEE Transactions on Information Theory* 57(10), 6976–6994.
- Sun, Tingni and Zhang, Cun-Hui (2012). Scaled Sparse Linear Regression. *Biometrika* 99(4), 879–898.
- Tibshirani, Robert (1996). Regression Shrinkage and Selection via the Lasso. *Journal of the Royal Statistical Society, Series B* 58(1), 267–288.
- Tibshirani, Ryan J. (2013). The Lasso Problem and Uniqueness. *Electronic Journal of Statistics* 7, 1456–1490.
- van de Geer, Sara and Bühlmann, Peter (2009). On the Conditions Used to Prove Oracle Results for the Lasso. *Electronic Journal of Statistics* 3, 1360–1392.
- van de Geer, Sara, Bühlmann, Peter, Ritov, Ya'acov and Dezeure, Ruben (2014). On Asymptotically Optimal Confidence Regions and Tests for High-Dimensional Models. *The Annals of Statistics* 42(3), 1166–1202.
- Wainwright, Martin J. (2009). Sharp Thresholds for High-Dimensional and Noisy Sparsity Recovery Using \( \ell_1 \)-Constrained Quadratic Programming (Lasso). *IEEE Transactions on Information Theory* 55(5), 2183–2202.
- Wainwright, Martin J. (2019). *High-Dimensional Statistics: A Non-Asymptotic Viewpoint*. Cambridge: Cambridge University Press.
- Wasserman, Larry and Roeder, Kathryn (2009). High-Dimensional Variable Selection. *The Annals of Statistics* 37(5A), 2178–2201.
- Zhang, Cun-Hui and Zhang, Stephanie S. (2014). Confidence Intervals for Low Dimensional Parameters in High Dimensional Linear Models. *Journal of the Royal Statistical Society, Series B* 76(1), 217–242.
- Zhao, Peng and Yu, Bin (2006). On Model Selection Consistency of Lasso. *Journal of Machine Learning Research* 7, 2541–2563.
- Zou, Hui (2006). The Adaptive Lasso and Its Oracle Properties. *Journal of the American Statistical Association* 101(476), 1418–1429.
