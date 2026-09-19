# Summary and Notes

## Summary

::: {.idea}

1. Mean squared error is variance plus squared bias. A linear estimator beats least squares in
   the matrix sense iff its squared bias, measured in the metric of the covariance it saves, is at
   most one (@thm-shr-bias-variance); none wins for every \( \bbeta \) (@exr-shr-no-linear-dominance).

2. Shrinking the \( j \)th canonical coordinate is governed by its signal-to-noise ratio
   \( \tau_j^2=d_j^2\alpha_j^2/\sigma^2 \): the oracle factor is \( \tau_j^2/(1+\tau_j^2) \), and dropping the
   direction helps iff \( \tau_j^2<1 \) (@prp-shr-canonical).

3. Ridge regression shrinks the \( j \)th direction by \( d_j^2/(d_j^2+\lambda) \), has effective degrees of
   freedom \( \sum_jd_j^2/(d_j^2+\lambda) \), and is the posterior mean under an exchangeable normal prior
   (@thm-shr-ridge). Every \( \lambda<2\sigma^2/\max_j\alpha_j^2 \) beats least squares, and every
   \( \lambda\le2\sigma^2/\norm{\bbeta}^2 \) does so in the matrix sense (@thm-shr-ridge-dominates).

4. Principal component regression keeps the leading \( k \) directions, and beats least squares in
   the matrix sense iff the dropped directions carry total signal \( \sum_{j>k}\tau_j^2\le1 \). The
   directions worth dropping are those with small \( \tau_j^2 \), not small \( d_j^2 \) (@prp-shr-pcr,
   @exm-shr-longley-pcr). Tuning or selecting from the data makes the minimized leave-one-out error
   optimistic; nested cross-validation is honest.

5. The lasso minimizes \( \tfrac12\norm{\y-\X\bb}^2+\lambda\norm{\bb}_1 \). Its solutions are the solutions of
   the KKT conditions and share one fitted vector (@def-shr-lasso, @thm-shr-lasso-kkt). Under an orthonormal design it
   is soft thresholding, against proportional shrinkage for ridge and hard thresholding for subset
   selection (@thm-shr-lasso-orthonormal). Coordinate descent computes it by repeated soft
   thresholding (@lem-shr-coordinate).

6. By simulation, best subset selection wins when a few coefficients are large, ridge wins when
   many are moderate, and the lasso sits between (@exm-shr-lasso-simulation).

7. With Stein's lemma (@lem-shr-stein), the James–Stein estimator \( (1-(p-2)\sigma^2/\norm{\Z}^2)\Z \) has
   risk below \( p\sigma^2 \) for every \( \boldsymbol{\uptheta} \) when \( p\ge3 \) (@thm-shr-james-stein), its
   positive part does better still (@prp-shr-positive-part), and in regression, with an estimated
   variance, the factor is about \( 1-1/F \) (@cor-shr-js-regression).

8. James–Stein is an empirical Bayes rule: it estimates the shrinkage factor
   \( \sigma^2/(\sigma^2+\omega^2) \) of a normal prior without bias, and pays at most \( 2\sigma^2 \) for not
   knowing \( \omega^2 \) (@prp-shr-empirical-bayes). Ridge with a fixed \( \lambda \) is the same rule
   with a guessed prior variance.

:::

## Notes and sources

**Coverage.** Christensen (2020, sections 13.3–13.5) treats canonical form, principal component
and ridge regression; Sen and Srivastava (1990, chapter 12) add a Stein-type shrinkage of least
squares; Fahrmeir, Kneib, Lang and Marx (2021, section 4.2) present ridge and the lasso as penalized
least squares. This chapter organizes
them around the shrinkage factors of @prp-shr-canonical. Examples, simulations and exercises are
the book's own.

**Mean squared error.** The matrix comparison of estimators and its equivalence with every
quadratic loss are due to Theobald (1974), who used them for the matrix-sense part of @thm-shr-ridge-dominates.

**Ridge regression.** Hoerl and Kennard (1970) named the method, introduced the ridge trace and
proved the existence theorem. The same regularization of ill-posed problems is due to Tikhonov
(Tikhonov and Arsenin 1977). Generalized cross-validation for ridge is due to Golub, Heath and
Wahba (1979). Hoerl, Kennard and Baldwin (1975) proposed the plug-in choice of \( \lambda \), and Draper and Van
Nostrand (1979) review the debate over data-chosen \( \lambda \).

**Principal components.** Hotelling (1933) introduced principal components and Massy (1965)
regression on them. The warnings about low-variance components are from Jolliffe (1982) and Hadi
and Ling (1998); Jolliffe (2002) is the standard monograph. Frank and Friedman (1993) compare
principal components, partial least squares and ridge regression.

**The lasso.** Tibshirani (1996) introduced the lasso, Chen, Donoho and Saunders (1998) the same
criterion as basis pursuit, and Park and Casella (2008) a Bayesian version. Tibshirani (2013) treats uniqueness and sparsity of solutions (@exr-shr-sparsity-bound). Path algorithms are due to Osborne, Presnell and Turlach (2000) and
Efron, Hastie, Johnstone and Tibshirani (2004); coordinate descent to Fu (1998), Friedman, Hastie,
Höfling and Tibshirani (2007) and Friedman, Hastie and Tibshirani (2010), with convergence theory
by Tseng (2001). Donoho and Johnstone (1994) analysed soft thresholding, and Hastie, Tibshirani
and Tibshirani (2020) compare subset selection and the lasso at scale. Hastie, Tibshirani and
Wainwright (2015) and Hastie, Tibshirani and Friedman (2009, chapter 3) cover these methods for
prediction.

**Stein's paradox.** The inadmissibility of \( \Z \) for \( p\ge3 \) is due to Stein (1956), and the
explicit estimator, including an estimated variance, to James and Stein (1961). The proof through
@lem-shr-stein follows Stein (1981); the empirical Bayes reading is from Efron and Morris (1973). Lehmann and Casella (1998, chapter 5) discuss admissibility and the positive-part
estimator; Efron and Hastie (2016, chapter 7) link James–Stein and ridge regression.

## References

- Boyd, Stephen and Vandenberghe, Lieven (2004). *Convex Optimization*. Cambridge: Cambridge University Press.
- Chen, Scott Shaobing, Donoho, David L. and Saunders, Michael A. (1998). Atomic Decomposition by Basis Pursuit. *SIAM Journal on Scientific Computing* 20(1), 33–61.
- Christensen, Ronald (2020). *Plane Answers to Complex Questions: The Theory of Linear Models*. 5th edition. Cham: Springer.
- Donoho, David L. and Johnstone, Iain M. (1994). Ideal Spatial Adaptation by Wavelet Shrinkage. *Biometrika* 81(3), 425–455.
- Draper, Norman R. and Van Nostrand, R. Craig (1979). Ridge Regression and James–Stein Estimation: Review and Comments. *Technometrics* 21(4), 451–466.
- Efron, Bradley and Hastie, Trevor (2016). *Computer Age Statistical Inference: Algorithms, Evidence, and Data Science*. Cambridge: Cambridge University Press.
- Efron, Bradley, Hastie, Trevor, Johnstone, Iain and Tibshirani, Robert (2004). Least Angle Regression. *The Annals of Statistics* 32(2), 407–499.
- Efron, Bradley and Morris, Carl (1973). Stein's Estimation Rule and Its Competitors—An Empirical Bayes Approach. *Journal of the American Statistical Association* 68(341), 117–130.
- Fahrmeir, Ludwig, Kneib, Thomas, Lang, Stefan and Marx, Brian D. (2021). *Regression: Models, Methods and Applications*. 2nd edition. Berlin: Springer.
- Frank, Ildiko E. and Friedman, Jerome H. (1993). A Statistical View of Some Chemometrics Regression Tools. *Technometrics* 35(2), 109–135.
- Friedman, Jerome, Hastie, Trevor, Höfling, Holger and Tibshirani, Robert (2007). Pathwise Coordinate Optimization. *The Annals of Applied Statistics* 1(2), 302–332.
- Friedman, Jerome, Hastie, Trevor and Tibshirani, Robert (2010). Regularization Paths for Generalized Linear Models via Coordinate Descent. *Journal of Statistical Software* 33(1), 1–22.
- Fu, Wenjiang J. (1998). Penalized Regressions: The Bridge versus the Lasso. *Journal of Computational and Graphical Statistics* 7(3), 397–416.
- Golub, Gene H., Heath, Michael and Wahba, Grace (1979). Generalized Cross-Validation as a Method for Choosing a Good Ridge Parameter. *Technometrics* 21(2), 215–223.
- Hadi, Ali S. and Ling, Robert F. (1998). Some Cautionary Notes on the Use of Principal Components Regression. *The American Statistician* 52(1), 15–19.
- Hastie, Trevor, Tibshirani, Robert and Friedman, Jerome (2009). *The Elements of Statistical Learning*. 2nd edition. New York: Springer.
- Hastie, Trevor, Tibshirani, Robert and Tibshirani, Ryan (2020). Best Subset, Forward Stepwise or Lasso? Analysis and Recommendations Based on Extensive Comparisons. *Statistical Science* 35(4), 579–592.
- Hastie, Trevor, Tibshirani, Robert and Wainwright, Martin (2015). *Statistical Learning with Sparsity: The Lasso and Generalizations*. Boca Raton, FL: CRC Press.
- Hoerl, Arthur E. and Kennard, Robert W. (1970). Ridge Regression: Biased Estimation for Nonorthogonal Problems. *Technometrics* 12(1), 55–67.
- Hoerl, Arthur E., Kennard, Robert W. and Baldwin, Kent F. (1975). Ridge Regression: Some Simulations. *Communications in Statistics* 4(2), 105–123.
- Hotelling, Harold (1933). Analysis of a Complex of Statistical Variables into Principal Components. *Journal of Educational Psychology* 24, 417–441 and 498–520.
- James, W. and Stein, Charles (1961). Estimation with Quadratic Loss. *Proceedings of the Fourth Berkeley Symposium on Mathematical Statistics and Probability* 1, 361–379.
- Jolliffe, Ian T. (1982). A Note on the Use of Principal Components in Regression. *Journal of the Royal Statistical Society. Series C (Applied Statistics)* 31(3), 300–303.
- Jolliffe, Ian T. (2002). *Principal Component Analysis*. 2nd edition. New York: Springer.
- Lehmann, E. L. and Casella, George (1998). *Theory of Point Estimation*. 2nd edition. New York: Springer.
- Longley, James W. (1967). An Appraisal of Least Squares Programs for the Electronic Computer from the Point of View of the User. *Journal of the American Statistical Association* 62(319), 819–841.
- Massy, William F. (1965). Principal Components Regression in Exploratory Statistical Research. *Journal of the American Statistical Association* 60(309), 234–256.
- Osborne, Michael R., Presnell, Brett and Turlach, Berwin A. (2000). A New Approach to Variable Selection in Least Squares Problems. *IMA Journal of Numerical Analysis* 20(3), 389–403.
- Park, Trevor and Casella, George (2008). The Bayesian Lasso. *Journal of the American Statistical Association* 103(482), 681–686.
- Sen, Ashish and Srivastava, Muni (1990). *Regression Analysis: Theory, Methods, and Applications*. New York: Springer.
- Stein, Charles (1956). Inadmissibility of the Usual Estimator for the Mean of a Multivariate Normal Distribution. *Proceedings of the Third Berkeley Symposium on Mathematical Statistics and Probability* 1, 197–206.
- Stein, Charles M. (1981). Estimation of the Mean of a Multivariate Normal Distribution. *The Annals of Statistics* 9(6), 1135–1151.
- Theobald, C. M. (1974). Generalizations of Mean Square Error Applied to Ridge Regression. *Journal of the Royal Statistical Society. Series B* 36(1), 103–106.
- Tibshirani, Robert (1996). Regression Shrinkage and Selection via the Lasso. *Journal of the Royal Statistical Society. Series B* 58(1), 267–288.
- Tibshirani, Ryan J. (2013). The Lasso Problem and Uniqueness. *Electronic Journal of Statistics* 7, 1456–1490.
- Tikhonov, Andrey N. and Arsenin, Vasiliy Y. (1977). *Solutions of Ill-Posed Problems*. Washington, DC: Winston.
- Tseng, Paul (2001). Convergence of a Block Coordinate Descent Method for Nondifferentiable Minimization. *Journal of Optimization Theory and Applications* 109(3), 475–494.
