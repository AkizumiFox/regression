# Summary and Notes

## Summary

::: {.idea}

1. Collinearity breaks no assumption. It inflates the variance of \( \mathbf{a}\T\hbeta \) when \( \mathbf{a} \) points along an
   eigenvector of \( \X\T\X \) with a small eigenvalue; the worst-to-best ratio is \( \kappa(\X)^2 \), and fitted values keep
   total variance \( \sigma^2p \) (@thm-col-variance).

2. Near a rank-deficient design, the functions with bounded variance are those estimable in the limit
   (@prp-col-near-estimable).

3. The variance inflation factor is the variance ratio to an orthogonal design of the same spread: scale-free,
   but a ratio, not a verdict (@def-col-vif, @prp-col-vif). Large factors that vanish under centring belong to
   questions about the origin (@exm-col-quadratic).

4. For a block of columns, the generalized factor \( \prod_i1/(1-\rho_i^2) \) measures the inflation of the joint
   confidence ellipsoid and does not depend on coding (@def-col-gvif, @prp-col-gvif).

5. Condition indices of the unit-length scaled design and the variance-decomposition proportions locate
   near-dependencies (@def-col-condition-index, @thm-col-decomposition). Centring can only hide them, namely
   those involving the intercept (@prp-col-scaling), and only high-leverage cases can mask them
   (@prp-col-deletion).

6. Predictions inside the convex hull, or along the pattern of the regressors, are unharmed
   (@prp-col-prediction).

7. Remedies: data along the weak direction (@prp-col-new-observation), a question the data can answer,
   outside information, or biased estimation ([Chapter 27](../ch27-shrinkage/index.html)).

:::

## Notes and sources

**The problem.** Frisch introduced the word "multicollinearity" in the 1930s. The statistical treatment rests on
the observation that the eigenvectors of \( \X\T\X \) sort linear functions of \( \bbeta \) into well and badly estimated
ones; Silvey (1969) set this out, including the advice to take new observations along the weak eigenvectors
behind @prp-col-new-observation. Goldberger (1991) insisted that collinearity is a shortage of data, not a
disease of the model. Stewart (1987) gave a numerical analyst's account.

**Diagnostics.** The name "variance inflation factor" is usually traced to Marquardt (1970). The generalized
factor and its \( 1/(2q) \) power are due to Fox and Monette (1992); the canonical-correlation form of
@prp-col-gvif(b) is implicit there. O'Brien (2007) examined the rules of thumb. Condition indices and
variance-decomposition proportions, with the working rule of indices above about \( 30 \) and two or more
proportions above one half, are due to Belsley, Kuh and Welsch (1980, chapter 3); Belsley (1991) is a
book-length development, and Belsley (1984), with its discussion, argues for scaling without centring. The
bounds @thm-col-decomposition(c) and @prp-col-deletion are stated in the form this book needs and follow from
standard eigenvalue inequalities. Among the textbooks used for coverage, Seber and Lee (2003, sections 9.7 and
10.7) treat variance inflation factors, eigenvalues, centring and influential cases; Christensen (2020,
sections 13.1–13.2) develops collinearity through ill-defined directions and tolerances; and Sen and
Srivastava (1990, chapter 10) give detailed variance-decomposition examples.

**Data and remedies.** The US macroeconomic series and Grunfeld's panel are public-domain data distributed
with statsmodels; Longley (1967) built his data to test regression programs. Mixed estimation is due to Theil
and Goldberger (1961). Shrinkage is the subject of [Chapter 27](../ch27-shrinkage/index.html), and selection
among collinear regressors of [Chapter 29](../ch29-model-selection/index.html).

## References

- Belsley, David A. (1984). Demeaning Conditioning Diagnostics through Centering. *The American Statistician* 38(2), 73–77.
- Belsley, David A. (1991). *Conditioning Diagnostics: Collinearity and Weak Data in Regression*. New York: Wiley.
- Belsley, David A., Kuh, Edwin and Welsch, Roy E. (1980). *Regression Diagnostics: Identifying Influential Data and Sources of Collinearity*. New York: Wiley.
- Christensen, Ronald (2020). *Plane Answers to Complex Questions: The Theory of Linear Models*. 5th edition. Cham: Springer.
- Fox, John and Monette, Georges (1992). Generalized Collinearity Diagnostics. *Journal of the American Statistical Association* 87(417), 178–183.
- Goldberger, Arthur S. (1991). *A Course in Econometrics*. Cambridge, MA: Harvard University Press.
- Longley, James W. (1967). An Appraisal of Least Squares Programs for the Electronic Computer from the Point of View of the User. *Journal of the American Statistical Association* 62(319), 819–841.
- Marquardt, Donald W. (1970). Generalized Inverses, Ridge Regression, Biased Linear Estimation, and Nonlinear Estimation. *Technometrics* 12(3), 591–612.
- O'Brien, Robert M. (2007). A Caution Regarding Rules of Thumb for Variance Inflation Factors. *Quality & Quantity* 41(5), 673–690.
- Seber, George A. F. and Lee, Alan J. (2003). *Linear Regression Analysis*. 2nd edition. Hoboken, NJ: Wiley.
- Sen, Ashish and Srivastava, Muni (1990). *Regression Analysis: Theory, Methods, and Applications*. New York: Springer.
- Silvey, S. D. (1969). Multicollinearity and Imprecise Estimation. *Journal of the Royal Statistical Society, Series B* 31(3), 539–552.
- Stewart, G. W. (1987). Collinearity and Least Squares Regression. *Statistical Science* 2(1), 68–84.
- Theil, Henri and Goldberger, Arthur S. (1961). On Pure and Mixed Statistical Estimation in Economics. *International Economic Review* 2(1), 65–78.
