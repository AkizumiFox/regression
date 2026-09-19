# Summary and Notes

## Summary

::: {.idea}

1. The data determine the mean vector \( \X\bbeta \) and nothing more. A function of \( \bbeta \) is identifiable iff it is a
           function of \( \X\bbeta \), and \( \bbeta \) itself is identifiable iff \( \rank(\X)=p \) (@def-est-identifiable,
           @thm-est-identifiable-mean and @prp-est-identifiable-rank).

2. A single coefficient is identifiable iff its column is not a combination of the others (@prp-est-coordinate). Rank
           deficiency comes from construction, from identities among variables (age, period and cohort), or from the design
           (@exm-est-apc).

3. A linear function \( \blambda\T\bbeta \) has a linear unbiased estimator iff it is identifiable iff
           \( \blambda\in\C(\X\T) \) (@def-est-estimable and @thm-est-estimable-identifiable). The estimable functions form a
           space of dimension \( \rank(\X) \) (@cor-est-estimable-space).

4. Eight equivalent tests of estimability, including \( \blambda\perp\Null(\X) \),
           \( \blambda\T\G\X\T\X=\blambda\T \), and invariance over least squares solutions (@thm-est-characterization). The
           estimate is \( \blambda\T\hbeta \) for any solution, with variance \( \sigma^2\blambda\T\G\blambda \) for any generalized
           inverse (@prp-est-ls-estimator).

5. For a nonestimable \( \blambda \), a generalized-inverse solution estimates a different, estimable function
           \( \blambda\T\G\X\T\X\bbeta \), and the reported variance can be any number (@prp-est-H and
           @exr-est-any-variance). Decide estimability with the SVD of \( \X \), and treat near-nonestimable functions as
           unestimated (@exm-est-near-estimable).

6. Full-rank reparameterizations are choices of basis for the estimable space (@thm-est-reparameterization). Side conditions
           \( \bT\bb=\bzero \) select a unique solution iff \( \rank(\bT)=p-r \) and no combination of them is estimable
           (@thm-est-side-conditions). The constrained coefficients are then particular estimable functions
           (@prp-est-side-meaning). The minimum-norm solution is itself a side condition (@prp-est-min-norm-side).

7. A coding of a factor is a nonsingular matrix \( \mathbf{K} \), and its coefficients estimate \( \mathbf{K}^{-1}\bmu \), the rows of
           \( \mathbf{K}^{-1} \) applied to the level means (@def-est-coding and @prp-est-coding).

8. In the additive two-way model the rank is \( a+b \) minus the number of connected components of the design, and a
           difference of levels is estimable iff the levels are connected (@thm-est-connected). In the interaction model,
           exactly the combinations of occupied cell means are estimable (@prp-est-interaction).

9. With a covariate, adjusted differences between levels subtract the slope times the difference in covariate means,
           the slope being estimated within levels (@prp-est-adjusted).

10. Contrasts are the estimable functions of a factor's effects (@def-est-contrast and @thm-est-oneway-contrasts).
            Orthogonal contrasts, in the design inner product \( \sum_kc_kd_k/n_k \), decompose the between-level sum of
            squares (@prp-est-contrast-ss). Coding coefficients equal scaled contrasts only for orthogonal codings
            (@prp-est-coding-dual).

:::

## Notes and sources

**Identifiability and estimability.** Christensen (2020, §2.1) develops estimability from identifiability, as this
chapter does, and stresses that the estimable functions are exactly the identifiable linear ones. Agresti (2015,
§1.4) gives a short treatment aimed at generalized linear models, where the same definitions apply to the linear predictor.
Rencher and Schaalje (2008, ch. 12) is a detailed textbook account of the non-full-rank model: estimable functions and
their tests, the estimators \( \mathbf{r}\T\X\T\y \) and \( \blambda\T\hbeta \), reparameterization and side conditions, with small
worked layouts. The notion of an estimable function goes back to R. C. Bose's work on the design of experiments in the
1940s. Scheffé (1959) and Searle (1971) made it standard, and Rao (1973) treats it through generalized inverses. The
characterization @thm-est-characterization collects conditions that are scattered through these sources. That a
generalized-inverse solution estimates the estimable function \( \blambda\T\G\X\T\X\bbeta \) (@prp-est-H) follows at once
from \( \E(\G\X\T\Y)=\G\X\T\X\bbeta \), but its consequence for software output deserves more emphasis than it usually gets.

**Side conditions and software.** The exact conditions of @thm-est-side-conditions, and the formula
\( (\X\T\X+\bT\T\bT)^{-1}\X\T\y \), are standard (see Rencher and Schaalje 2008, §12.6, and Searle 1971). The coding
matrices, called “contrasts”, that S and R use for factors are described by Chambers and Hastie (1992), and Venables and
Ripley (2002) discuss the treatment, sum, Helmert and polynomial codings of R. The minimum-norm solution rests on the Moore–Penrose inverse (Penrose 1955), and column-pivoted QR is described in
Golub and Van Loan (2013).

**Indicator variables.** Sen and Srivastava (1990, ch. 4) is a practical account of indicator variables: the two-sample
comparison as a regression, polychotomous factors, hierarchical codings, factors with continuous covariates, broken-line
regression and aggregated indicators as responses. Where this chapter meets those themes it takes its own route: the
two-sample comparison through the coding proposition (@exr-est-two-sample), the indicator of a single observation through
the Frisch–Waugh–Lovell theorem (@exr-est-single-dummy), and count regressors without intercept in a new setting
(@exr-est-counts).

**Connectedness, empty cells and cohorts.** Connectedness of two-way designs, and its role in estimability, is treated at
length by Searle (1971, 1987). Searle (1987) is also the standard reference for unbalanced data and empty cells, which
[Chapter 17](../ch17-unbalanced-data/index.html) takes up. The age–period–cohort problem of @exm-est-apc is a well-known instance of a dependency among
variables. Holford (1983) showed which functions of the age, period and cohort effects are estimable, among them the
second differences of @exr-est-apc-factors.

**Data.** The survey examples use the extract of the 1996 American National Election Study distributed with statsmodels
(Seabold and Perktold 2010) as a public-domain data set. The other examples use simulated data, generated by the scripts
in `code/ch08/`.

## References

- Agresti, Alan (2015). *Foundations of Linear and Generalized Linear Models*. Hoboken, NJ: Wiley.
- Chambers, John M. and Hastie, Trevor J., eds. (1992). *Statistical Models in S*. Pacific Grove, CA: Wadsworth & Brooks/Cole.
- Christensen, Ronald (2020). *Plane Answers to Complex Questions: The Theory of Linear Models*. 5th edition. Cham: Springer.
- Golub, Gene H. and Van Loan, Charles F. (2013). *Matrix Computations*. 4th edition. Baltimore: Johns Hopkins University Press.
- Holford, Theodore R. (1983). The Estimation of Age, Period and Cohort Effects for Vital Rates. *Biometrics* 39(2), 311–324.
- Penrose, Roger (1955). A Generalized Inverse for Matrices. *Mathematical Proceedings of the Cambridge Philosophical Society* 51(3), 406–413.
- Rao, C. Radhakrishna (1973). *Linear Statistical Inference and Its Applications*. 2nd edition. New York: Wiley.
- Rencher, Alvin C. and Schaalje, G. Bruce (2008). *Linear Models in Statistics*. 2nd edition. Hoboken, NJ: Wiley.
- Scheffé, Henry (1959). *The Analysis of Variance*. New York: Wiley.
- Seabold, Skipper and Perktold, Josef (2010). Statsmodels: Econometric and Statistical Modeling with Python. *Proceedings of the 9th Python in Science Conference*, 92–96.
- Searle, Shayle R. (1971). *Linear Models*. New York: Wiley.
- Searle, Shayle R. (1987). *Linear Models for Unbalanced Data*. New York: Wiley.
- Sen, Ashish and Srivastava, Muni (1990). *Regression Analysis: Theory, Methods, and Applications*. New York: Springer.
- Venables, William N. and Ripley, Brian D. (2002). *Modern Applied Statistics with S*. 4th edition. New York: Springer.
