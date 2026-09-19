# Summary and Notes

## Summary

::: {.idea}

1. \( \Y \) is multivariate normal when every linear combination \( \bm a\T\Y \) is normal,
           constants included (@def-mvn). The definition covers singular covariance
           matrices, which regression needs for fitted values and residuals.

2. \( \bmu+\A\Z \), with \( \Z \) standard normal, is \( \Normal(\bmu,\A\A\T) \). Every mean and
           every nonnegative definite covariance occurs (@thm-mvn-existence).

3. The moment generating function is \( \exp(\bm t\T\bmu+\tfrac12\bm t\T\bSigma\bm t) \), so
           a normal distribution is determined by \( \bmu \) and \( \bSigma \) (@thm-mvn-mgf).
           When \( \bSigma \) is positive definite there is a density whose contours are the
           ellipsoids of constant Mahalanobis distance, and \( \Delta^2(\Y)\sim\chi^2(n) \)
           (@thm-mvn-density and @prp-mvn-mahalanobis).

4. Affine maps of normal vectors are normal, with no rank condition, and so are all
           subvectors (@thm-mvn-linear and @cor-mvn-marginals). A normal vector of rank \( r \) is an
           exact affine image of \( r \) independent standard normals
           (@prp-mvn-rank-representation).

5. Normal margins do not imply joint normality, and zero correlation between normal
           variables does not imply independence without it
           (@exm-mvn-sign-flip and @exm-mvn-mixture).

6. For jointly normal blocks, zero cross-covariance is equivalent to independence.
           \( \A\Y \) and \( \B\Y \) are independent iff \( \A\bSigma\B\T=\bzero \)
           (@thm-mvn-independence and @cor-mvn-AY-BY).

7. \( \Y_1\mid\Y_2=\y_2\sim\Normal\bigl(\bmu_1+\bSigma_{12}\bSigma_{22}\ginv(\y_2-\bmu_2),\
        \bSigma_{11\cdot2}\bigr) \): linear mean, constant covariance, reduced variance
           (@thm-mvn-conditional). The conditional mean is the best predictor of all, and
           coincides with the best linear predictor (@prp-mvn-best-predictor).

8. The partial correlation is the correlation of the residuals from linear prediction
           on the conditioning variables. Under normality it is the conditional correlation,
           and it is zero iff the variables are conditionally independent
           (@def-mvn-partial-correlation and @prp-mvn-partial-meaning). It can be computed
           by a recursion, from the inverse covariance matrix, or by FWL residuals
           (@prp-mvn-partial-recursion and @prp-mvn-precision).

:::

## Notes and sources

**Three definitions.**  Linear-model texts reach the multivariate normal by
different routes, and comparing them explains the choices made here.
Rencher and Schaalje (2008, chapter 4) start from the positive definite case, obtain the
density by transforming independent standard normals with \( \bSigma^{1/2} \), and prove
the properties of linear functions and independence with moment generating functions,
for full-rank transformations.
Seber and Lee (2003, chapter 2) also begin with the density, then extend the
definition to singular covariance matrices through the representation \( \A\Z+\bmu \) and
prove that normality of all linear combinations characterizes the family. Their
moment-generating-function proof of the independence theorem is the model for
[Section 3.4](04-independence.html).
Christensen (2020, §1.2) takes \( \A\Z+\bmu \) as the definition from the outset,
shows with characteristic functions that the law depends only on \( \bmu \) and \( \bSigma \), and
proves the independence theorem by building an explicit independent version. Defining the
distribution through linear combinations, as this chapter does, goes back at least to
Rao (1973, chapter 8). It has the advantage that the theorem on affine
transformations needs no rank condition and a one-line proof.

**Cramér–Wold.**  The device of reducing a vector to its linear combinations is
from Cramér and Wold (1936), where it was used to reduce multivariate limit theorems to
univariate ones. [Chapter 2](../ch02-random-vectors/index.html) proves it, together with the uniqueness and
factorization properties of moment generating functions used in [Section 3.2](02-mgf-density.html).
Billingsley (1995) is the reference for these facts and for the
conditioning property of independent variables used in the proof of
@thm-mvn-conditional.

**Conditional distributions and partial correlation.**  The proof of
@thm-mvn-conditional through the independent error \( \W \) borrows one idea from
Rencher and Schaalje (2008) and Seber and Lee (2003): both pass through the transformed
vector \( \Y_1-\B\Y_2 \), which is independent of \( \Y_2 \), although both still obtain the
conditional density as a ratio of densities in the positive definite case. The proof here
uses only independence and a generalized inverse, so it covers singular \( \bSigma_{22} \),
an extension that is standard in multivariate analysis (Rao 1973). Anderson (2003)
is the classical reference for the multivariate normal in statistics, including the
distribution theory of sample partial and multiple correlation coefficients that
Chapter 14 uses. Rencher and Schaalje (2008, §4.5) define
partial correlation as a conditional correlation, as in
@prp-mvn-partial-meaning(b). The residual interpretation (a) does not require
normality and is the one that connects to least squares. Regression toward the mean was
described by Galton (1886). The link between zeros of the inverse covariance
matrix and conditional independence was made the basis of a modelling strategy by
Dempster (1972), and Lauritzen (1996) develops the theory of
graphical models that grew from it.

**Data.**  @exm-mvn-longley uses the macroeconomic series published by
Longley (1967) as a test of numerical accuracy, available in
`statsmodels`. The same series appear in [Chapter 2](../ch02-random-vectors/index.html), as an example of
sample moments, and in [Chapter 6](../ch06-projections/index.html), to illustrate ill-conditioning.

## References

- Anderson, T. W. (2003). *An Introduction to Multivariate Statistical Analysis*. 3rd edition. Hoboken, NJ: Wiley.
- Billingsley, Patrick (1995). *Probability and Measure*. 3rd edition. New York: Wiley.
- Christensen, Ronald (2020). *Plane Answers to Complex Questions: The Theory of Linear Models*. 5th edition. Cham: Springer.
- Cramér, Harald and Wold, Herman (1936). Some Theorems on Distribution Functions. *Journal of the London Mathematical Society* s1-11(4), 290–294.
- Dempster, A. P. (1972). Covariance Selection. *Biometrics* 28(1), 157–175.
- Frisch, Ragnar and Waugh, Frederick V. (1933). Partial Time Regressions as Compared with Individual Trends. *Econometrica* 1(4), 387–401.
- Galton, Francis (1886). Regression Towards Mediocrity in Hereditary Stature. *The Journal of the Anthropological Institute of Great Britain and Ireland* 15, 246–263.
- Lauritzen, Steffen L. (1996). *Graphical Models*. Oxford: Oxford University Press.
- Longley, James W. (1967). An Appraisal of Least Squares Programs for the Electronic Computer from the Point of View of the User. *Journal of the American Statistical Association* 62(319), 819–841.
- Rao, C. Radhakrishna (1973). *Linear Statistical Inference and Its Applications*. 2nd edition. New York: Wiley.
- Rencher, Alvin C. and Schaalje, G. Bruce (2008). *Linear Models in Statistics*. 2nd edition. Hoboken, NJ: Wiley.
- Seber, George A. F. and Lee, Alan J. (2003). *Linear Regression Analysis*. 2nd edition. Hoboken, NJ: Wiley.
