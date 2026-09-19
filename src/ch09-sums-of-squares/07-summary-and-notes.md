# Summary and notes

## Summary

::: {.idea}

1. An analysis of variance is a choice of mutually orthogonal subspaces filling the model
   space. Sums of squares are squared lengths of the projections of \( \y \), degrees of freedom are
   dimensions, and both add up by Pythagoras. Decompositions and chains of nested models are the
   same thing (@def-ss-anova-decomposition, @thm-ss-decomposition).

2. The regression table comes from the chain \( \spn(\bone)\subseteq\C(\X) \). The corrected identity
   \( \text{SST}=\text{SSR}+\text{SSE} \) needs \( \bone\in\C(\X) \); without it an extra term
   \( -2\bar{y}\sum_i\hat{\varepsilon}_i \) appears (@prp-ss-computing-formulas, @prp-ss-no-intercept).

3. The extra sum of squares \( \text{SS}(\X_2\mid\X_1)=\y\T(\M-\M_1)\y \) depends only on column
   spaces. Adding regressors to a fitted model needs only a \( t\times t \) inverse, and never lowers
   the variances of the old coefficients. For one column it equals \( t_j^2s^2 \)
   (@def-ss-extra, @thm-ss-adding, @cor-ss-single-column).

4. The sum of squares of an estimable hypothesis \( \bLambda\T\bbeta=\bm d \) is
   \( (\bLambda\T\hbeta-\bm d)\T(\bLambda\T\G\bLambda)^{-1}(\bLambda\T\hbeta-\bm d) \), whatever the rank.
   A nonestimable constraint restricts only through its estimable part
   (@thm-ss-hypothesis, @thm-ss-restricted, @prp-ss-effective-df).

5. Sequential (Type I) sums of squares decompose the regression sum of squares but depend on the
   order. Partial, Type II and Type III sums of squares do not decompose anything, and in models
   with interactions they address different hypotheses: weighted row means, additive main effects,
   unweighted row means, or, with reference coding, means at a single level of the other factor
   (@def-ss-sequential, @def-ss-partial, @prp-ss-tested-regression, @thm-ss-two-way-hypotheses).

6. Order stops mattering exactly when the centred subspaces of the terms are orthogonal. For two
   factors this is the proportional-frequency condition. Orthogonal polynomials and orthogonal
   contrasts manufacture orthogonality inside a single term
   (@thm-ss-orthogonal-design, @thm-ss-proportional, @thm-ss-single-df).

7. \( R \) is the largest correlation of \( \y \) with a linear combination of the regressors. Under
   the null \( R^2\sim\mathrm{Beta}\bigl(\tfrac{r-1}2,\tfrac{n-r}2\bigr) \), with mean \( (r-1)/(n-1) \). The
   adjusted \( R^2 \) has mean zero under the null and rises exactly when an added term has
   \( F>1 \). Partial \( R^2 \) equals \( t^2/(t^2+n-r) \) (@thm-ss-r2-max-correlation, @thm-ss-r2-null,
   @prp-ss-adjusted-r2, @thm-ss-partial-r2).

8. \( R^2 \) depends on the design, and it is not a test of the model, a comparison across responses,
   or a measure of predictive accuracy or of cause ([Section 9.5](05-r-squared.html)).

9. Every mean square has expectation \( \sigma^2+\norm{\bP\bmu}^2/\rank(\bP) \), using only second
   moments. Under normality the pieces of a decomposition are independent noncentral chi-squared
   variables, and an \( F \) ratio compares two estimates of \( \sigma^2 \) that agree exactly when its
   hypothesis \( \bP\bmu=\bzero \) holds (@thm-ss-expected-mean-squares, @cor-ss-ms-distributions).

:::

## Notes and sources

**Analysis of variance.**  The decomposition of a total sum of squares into components with
their own degrees of freedom is Fisher's. An analysis of variance already appears in Fisher and Mackenzie (1923), and Fisher set out the
analysis of variance table for general use in *Statistical Methods for Research Workers* (Fisher 1925), after making the notion of degrees of
freedom precise in Fisher (1922). The distribution theory that makes the components independent
chi-squared variables is Cochran's theorem (Cochran 1934), proved in
[Chapter 4](../ch04-quadratic-forms/index.html). The coordinate-free view, in which a table is a list of orthogonal subspaces,
follows Christensen (2020). His section 3.6, which splits a test space into orthogonal lines, is the
source of @thm-ss-single-df, and his section 7.4 treats proportional numbers. Scheffé (1959)
remains a clear account of the classical theory, including expected mean squares for random
effects.

**Extra sums of squares, adding regressors and restrictions.**  Seber and Lee (2003, sections
3.7–3.9) give the updating formulas of @thm-ss-adding, restricted least squares, and their
rank-deficient versions. Our proofs go through the Frisch–Waugh–Lovell theorem and Pythagoras
instead of Lagrange multipliers. @prp-ss-effective-df makes precise the remark, common to the
sources, that only the estimable part of a hypothesis can be tested.

**Types of sums of squares.**  Agresti (2015, sections 2.3–2.4) presents sequential and partial
sums of squares through projection matrices, with the decomposition of the regression sum of
squares when the regressors are uncorrelated, and notes the Type I, II and III terminology. The
problem of main effects in unbalanced layouts with interaction goes back to Yates (1934), whose
“weighted squares of means” is closely related to what software now calls Type III. The names of
the types come from statistical software, and programs have used them inconsistently, which is
why this chapter defines each type by its projection. Nelder (1977) argued that tests
of main effects in the presence of interactions violate the marginality of the model terms; his
position underlies the preference for Type II sums of squares when interactions are absent.
Speed, Hocking and Hackney (1978) identify the cell-mean hypotheses tested by the various sums of
squares, and Searle (1987) treats unbalanced data exhaustively in the same cell-mean terms.

**The coefficient of determination.**  The distribution of the sample multiple correlation
coefficient, under the null and in general for random normal regressors, is due to Fisher (1928);
@thm-ss-r2-null is the simpler fixed-regressor version. Kvålseth (1985) surveys the many
different formulas that software has used for \( R^2 \), especially for models without an
intercept, and Anscombe (1973) gives four data sets with the same \( R^2 \) and very different
relationships, a standing warning that \( R^2 \) is not a test of the model. 

**Data.**  The state crime, election study, Grunfeld and Nile data are public-domain data sets
distributed with statsmodels. The Nile series is the one analysed for a change point by
Cobb (1978).

**What is new here.**  The organization of the chapter around subspaces, the proofs of
@thm-ss-orthogonal-design and @thm-ss-proportional, the reference-coding case of
@thm-ss-two-way-hypotheses, and the numerical examples and simulations are this book's.

## References

- Agresti, Alan (2015). *Foundations of Linear and Generalized Linear Models*. Hoboken, NJ: Wiley.
- Anscombe, Francis J. (1973). Graphs in Statistical Analysis. *The American Statistician* 27(1), 17–21.
- Christensen, Ronald (2020). *Plane Answers to Complex Questions: The Theory of Linear Models*. 5th edition. Cham: Springer.
- Cobb, George W. (1978). The Problem of the Nile: Conditional Solution to a Changepoint Problem. *Biometrika* 65(2), 243–251.
- Cochran, William G. (1934). The Distribution of Quadratic Forms in a Normal System, with Applications to the Analysis of Covariance. *Proceedings of the Cambridge Philosophical Society* 30(2), 178–191.
- Fisher, Ronald A. (1922). On the Interpretation of \( \chi^2 \) from Contingency Tables, and the Calculation of P. *Journal of the Royal Statistical Society* 85(1), 87–94.
- Fisher, Ronald A. (1925). *Statistical Methods for Research Workers*. Edinburgh: Oliver and Boyd.
- Fisher, Ronald A. (1928). The General Sampling Distribution of the Multiple Correlation Coefficient. *Proceedings of the Royal Society of London. Series A* 121(788), 654–673.
- Fisher, Ronald A. and Mackenzie, Winifred A. (1923). Studies in Crop Variation. II. The Manurial Response of Different Potato Varieties. *Journal of Agricultural Science* 13(3), 311–320.
- Kvålseth, Tarald O. (1985). Cautionary Note about \( R^2 \). *The American Statistician* 39(4), 279–285.
- Nelder, John A. (1977). A Reformulation of Linear Models. *Journal of the Royal Statistical Society. Series A* 140(1), 48–77.
- Scheffé, Henry (1959). *The Analysis of Variance*. New York: Wiley.
- Searle, Shayle R. (1987). *Linear Models for Unbalanced Data*. New York: Wiley.
- Seber, George A. F. and Lee, Alan J. (2003). *Linear Regression Analysis*. 2nd edition. Hoboken, NJ: Wiley.
- Speed, F. M., Hocking, Ronald R. and Hackney, O. P. (1978). Methods of Analysis of Linear Models with Unbalanced Data. *Journal of the American Statistical Association* 73(361), 105–112.
- Yates, Frank (1934). The Analysis of Multiple Classifications with Unequal Numbers in the Different Classes. *Journal of the American Statistical Association* 29(185), 51–66.
