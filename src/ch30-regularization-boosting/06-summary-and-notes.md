# Summary and Notes

## Summary

::: {.idea}

1. A penalized criterion (@def-reg-penalized) has a minimizer when it is coercive; with a convex penalty all minimizers share
   their fit and penalty, and strict convexity or full rank gives uniqueness (@thm-reg-existence).

2. For convex penalties, \( \X\T(\y-\X\hbeta)\in\lambda\,\partial P(\hbeta) \) characterizes minimizers (@lem-reg-optimality,
   @lem-reg-norm-subgradient), and penalized estimates are posterior modes (@prp-reg-map).

3. SCAD and MCP are sparse, continuous and nearly unbiased; MCP is firm thresholding (@prp-reg-mcp).

4. The elastic net is unique, a lasso on augmented data, and groups correlated regressors (@prp-reg-elastic-net,
   @prp-reg-rescaled).

5. The group lasso selects whole groups by block soft thresholding (@prp-reg-group-lasso) and, orthonormalized, ignores the
   coding of factors (@prp-reg-group-invariance).

6. Componentwise boosting converges geometrically to least squares (@thm-reg-boosting); early stopping is a spectral filter
   like ridge (@prp-reg-boosting-operator), and the operator trace understates its degrees of freedom.

7. Forward stagewise tracks the lasso path (@prp-reg-stagewise-orthonormal); its limit is the lasso when paths are monotone
   and the monotone lasso otherwise (@prp-reg-stagewise).

:::

## Notes and sources

**Coverage.** Fahrmeir, Kneib, Lang and Marx (2021, §4.2–4.3) treat penalties and componentwise boosting in applied terms;
the theory and proofs are added here. Hastie, Tibshirani and Friedman (2009, chapters 3, 10 and 16) cover the same ground from
the machine-learning side.

**Penalties.** Ridge regression is due to Hoerl and Kennard (1970), the lasso to Tibshirani (1996) and bridge penalties to Frank
and Friedman (1993); uniqueness of the lasso in general position is from Tibshirani (2013). Boyd and Vandenberghe (2004) cover
subgradients, and Park and Casella (2008) develop the Bayesian lasso. SCAD and the oracle property are from Fan and Li (2001), MCP
from Zhang (2010).

**Elastic net and groups.** The elastic net, its grouping effect and rescaling are from Zou and Hastie (2005); our grouping
bound adds the residual-norm term. Pathwise coordinate descent is from Friedman, Hastie, Höfling and Tibshirani (2007) and
Friedman, Hastie and Tibshirani (2010). The group lasso is due to Yuan and Lin (2006), its orthonormalized form is discussed by
Simon and Tibshirani (2012), and the sparse group lasso is from Simon, Friedman, Hastie and Tibshirani (2013).

**Boosting.** AdaBoost is due to Freund and Schapire (1997); Friedman, Hastie and Tibshirani (2000) and Friedman (2001) gave the
statistical and gradient view, and matching pursuit is from Mallat and Zhang (1993). \( L_2 \) boosting and the boosting operator
are from Bühlmann and Yu (2003), componentwise boosting in high dimensions from Bühlmann (2006), and a survey is Bühlmann and
Hothorn (2007). Our proof of @thm-reg-boosting is the standard argument for greedy coordinate descent on a quadratic. Efron,
Hastie, Johnstone and Tibshirani (2004) proved that FS\(_0\) is the lasso path under their positive cone condition, and Hastie,
Taylor, Tibshirani and Walther (2007) described FS\(_0\) in general as the monotone lasso; @prp-reg-stagewise-orthonormal is a
special case proved directly.

## References

- Boyd, Stephen and Vandenberghe, Lieven (2004). *Convex Optimization*. Cambridge: Cambridge University Press.
- Bühlmann, Peter (2006). Boosting for High-Dimensional Linear Models. *The Annals of Statistics* 34(2), 559–583.
- Bühlmann, Peter and Hothorn, Torsten (2007). Boosting Algorithms: Regularization, Prediction and Model Fitting. *Statistical Science* 22(4), 477–505.
- Bühlmann, Peter and Yu, Bin (2003). Boosting with the L2 Loss: Regression and Classification. *Journal of the American Statistical Association* 98(462), 324–339.
- Efron, Bradley, Hastie, Trevor, Johnstone, Iain and Tibshirani, Robert (2004). Least Angle Regression. *The Annals of Statistics* 32(2), 407–499.
- Fahrmeir, Ludwig, Kneib, Thomas, Lang, Stefan and Marx, Brian D. (2021). *Regression: Models, Methods and Applications*. 2nd edition. Berlin: Springer.
- Fan, Jianqing and Li, Runze (2001). Variable Selection via Nonconcave Penalized Likelihood and Its Oracle Properties. *Journal of the American Statistical Association* 96(456), 1348–1360.
- Frank, Ildiko E. and Friedman, Jerome H. (1993). A Statistical View of Some Chemometrics Regression Tools. *Technometrics* 35(2), 109–135.
- Freund, Yoav and Schapire, Robert E. (1997). A Decision-Theoretic Generalization of On-Line Learning and an Application to Boosting. *Journal of Computer and System Sciences* 55(1), 119–139.
- Friedman, Jerome H. (2001). Greedy Function Approximation: A Gradient Boosting Machine. *The Annals of Statistics* 29(5), 1189–1232.
- Friedman, Jerome, Hastie, Trevor, Höfling, Holger and Tibshirani, Robert (2007). Pathwise Coordinate Optimization. *The Annals of Applied Statistics* 1(2), 302–332.
- Friedman, Jerome, Hastie, Trevor and Tibshirani, Robert (2000). Additive Logistic Regression: A Statistical View of Boosting. *The Annals of Statistics* 28(2), 337–407.
- Friedman, Jerome, Hastie, Trevor and Tibshirani, Robert (2010). Regularization Paths for Generalized Linear Models via Coordinate Descent. *Journal of Statistical Software* 33(1), 1–22.
- George, Edward I. and McCulloch, Robert E. (1993). Variable Selection via Gibbs Sampling. *Journal of the American Statistical Association* 88(423), 881–889.
- Gertheiss, Jan and Tutz, Gerhard (2010). Sparse Modeling of Categorial Explanatory Variables. *The Annals of Applied Statistics* 4(4), 2150–2180.
- Hastie, Trevor, Taylor, Jonathan, Tibshirani, Robert and Walther, Guenther (2007). Forward Stagewise Regression and the Monotone Lasso. *Electronic Journal of Statistics* 1, 1–29.
- Hastie, Trevor, Tibshirani, Robert and Friedman, Jerome (2009). *The Elements of Statistical Learning*. 2nd edition. New York: Springer.
- Hoerl, Arthur E. and Kennard, Robert W. (1970). Ridge Regression: Biased Estimation for Nonorthogonal Problems. *Technometrics* 12(1), 55–67.
- Leeb, Hannes and Pötscher, Benedikt M. (2005). Model Selection and Inference: Facts and Fiction. *Econometric Theory* 21(1), 21–59.
- Mallat, Stéphane G. and Zhang, Zhifeng (1993). Matching Pursuits with Time-Frequency Dictionaries. *IEEE Transactions on Signal Processing* 41(12), 3397–3415.
- Park, Trevor and Casella, George (2008). The Bayesian Lasso. *Journal of the American Statistical Association* 103(482), 681–686.
- Simon, Noah, Friedman, Jerome, Hastie, Trevor and Tibshirani, Robert (2013). A Sparse-Group Lasso. *Journal of Computational and Graphical Statistics* 22(2), 231–245.
- Simon, Noah and Tibshirani, Robert (2012). Standardization and the Group Lasso Penalty. *Statistica Sinica* 22(3), 983–1001.
- Tibshirani, Robert (1996). Regression Shrinkage and Selection via the Lasso. *Journal of the Royal Statistical Society, Series B* 58(1), 267–288.
- Tibshirani, Robert, Saunders, Michael, Rosset, Saharon, Zhu, Ji and Knight, Keith (2005). Sparsity and Smoothness via the Fused Lasso. *Journal of the Royal Statistical Society, Series B* 67(1), 91–108.
- Tibshirani, Ryan J. (2013). The Lasso Problem and Uniqueness. *Electronic Journal of Statistics* 7, 1456–1490.
- Tseng, Paul (2001). Convergence of a Block Coordinate Descent Method for Nondifferentiable Minimization. *Journal of Optimization Theory and Applications* 109(3), 475–494.
- Yuan, Ming and Lin, Yi (2006). Model Selection and Estimation in Regression with Grouped Variables. *Journal of the Royal Statistical Society, Series B* 68(1), 49–67.
- Zhang, Cun-Hui (2010). Nearly Unbiased Variable Selection under Minimax Concave Penalty. *The Annals of Statistics* 38(2), 894–942.
- Zhang, Tong and Yu, Bin (2005). Boosting with Early Stopping: Convergence and Consistency. *The Annals of Statistics* 33(4), 1538–1579.
- Zou, Hui (2006). The Adaptive Lasso and Its Oracle Properties. *Journal of the American Statistical Association* 101(476), 1418–1429.
- Zou, Hui and Hastie, Trevor (2005). Regularization and Variable Selection via the Elastic Net. *Journal of the Royal Statistical Society, Series B* 67(2), 301–320.
- Zou, Hui, Hastie, Trevor and Tibshirani, Robert (2007). On the "Degrees of Freedom" of the Lasso. *The Annals of Statistics* 35(5), 2173–2192.
- Zou, Hui and Li, Runze (2008). One-Step Sparse Estimates in Nonconcave Penalized Likelihood Models. *The Annals of Statistics* 36(4), 1509–1533.
