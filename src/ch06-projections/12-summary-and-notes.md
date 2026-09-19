# Summary and Notes

## Summary

::: {.idea}

1. Least squares is the nearest-point problem in observation space. The fitted vector
           \( \hY=\M\y \) is the orthogonal projection of \( \y \) onto \( \C(\X) \), and the residual is
           orthogonal to \( \C(\X) \) (@thm-proj-ls-projection).

2. Orthogonal projections are exactly the symmetric idempotent matrices. Their rank
           equals their trace (@thm-proj-sym-idem and @prp-proj-trace-rank).

3. \( \M=\X(\X\T\X)\ginv\X\T \) for *every* generalized inverse
           (@thm-proj-M-formula). Fitted values, residuals and every function of
           them are unique, whatever the rank. Coefficients are unique only when
           \( \rank(\X)=p \), and \( \blambda\T\hbeta \) is invariant iff \( \blambda\in\C(\X\T) \)
           (@thm-proj-normal-equations and @thm-proj-invariant-functions).

4. Nested model spaces give orthogonal decompositions
           \( \I=\Mo+(\M-\Mo)+(\I-\M) \). The drop in residual sum of squares is
           \( \norm{(\M-\Mo)\y}^2 \), with \( r-r_0 \) degrees of freedom
           (@thm-proj-nested).

5. A multiple-regression coefficient is the slope on the part of its regressor
           orthogonal to the others. Fixed-effects estimators and added-variable plots are
           applications (@thm-proj-fwl).

6. Anything that depends only on \( \C(\X) \) is invariant to reparameterization. \( R^2 \)
           is the squared cosine between centred data and centred fit
           (@thm-proj-reparam and @thm-proj-r2-cosine).

7. Leverage \( h_{ii}=m_{ii} \) is the squared length of the projection of a coordinate
           vector. It lies in \( [1/n,1] \) with an intercept and averages \( r/n \)
           (@prp-proj-leverage).

8. Replacing the inner product by \( \bu\T\V^{-1}\bv \) turns least squares into
           generalized least squares. The projection is then idempotent and \( \V^{-1} \)-self-adjoint,
           but not symmetric (@thm-proj-A-projection and @thm-proj-kruskal).

9. Compute projections with orthonormal bases (QR, SVD), never by inverting
           \( \X\T\X \), which squares the condition number ([Section 6.10](10-computation.html)).

10. In the population, least squares estimates the \( L^2 \) projection of \( Y \) onto
            linear functions of \( \mathbf{X} \), whether or not the regression function is linear
            (@thm-proj-blp and @prp-proj-consistency).

:::

## Notes and sources

**The geometric tradition.**  Seeing least squares as projection is as old as
least squares itself, but it became the organizing principle of linear-model theory
only in the twentieth century. Herr (1980) traces the history, from
Fisher's geometric derivations of sampling distributions through the coordinate-free
treatments of the 1960s. Book-length geometric developments at different levels include
Saville and Wood (1991), written for a first course, and
Wickens (1995), for multivariate methods.
Christensen (2020) is the most sustained vector-space treatment of linear
model theory. Its chapter 2 and appendices A–B cover the material of
[Section 6.2](02-subspaces.html), [Section 6.3](03-orthogonal-projection.html) and [Section 6.4](04-least-squares.html) for
the rank-deficient case throughout, and influenced this chapter's choice to state
results for general rank from the start. Seber and Lee (2003, appendix B) gives a
compact account of projections, including the constraint-space result that appears here
as @thm-proj-constraint-space. Agresti (2015, chapter 2) develops
projection matrices through worked null-model, one-way and two-way examples, and is a
gentle companion to [Section 6.5](05-nested.html). For the underlying linear algebra,
Halmos (1958) remains a model of coordinate-free exposition.

**Generalized inverses.**  The four-condition inverse is due to
Penrose (1955). Its statistical uses, and the theory of general
\( \A\ginv \) including the characterization in @exr-proj-all-ginverses, are developed
by Rao (1973). The invariance of \( \X(\X\T\X)\ginv\X\T \) (@lem-proj-ginverse-invariance)
is the fact that makes generalized inverses
harmless in linear models. Every quantity of statistical interest is sandwiched in
this way.

**Frisch–Waugh–Lovell.**  The theorem is due to Frisch and Waugh (1933)
for time trends and Lovell (1963) for seasonal dummies and general
partitions. It is also known under Lovell's name alone, and in econometrics as the
“partialling-out” result. It underlies added-variable plots, double/debiased
machine-learning estimators, and the fixed-effects (“within”) estimator for
panel data.

**Leverage.**  Hoaglin and Welsch (1978) introduced the name “hat matrix”
into common use and established the basic properties in @prp-proj-leverage.
Diagnostic use of leverage is taken up in
Chapter 20.

**Other inner products.**  Generalized least squares is due to
Aitken (1935). The coordinate-free condition for equality of ordinary and
generalized least squares (@thm-proj-kruskal) is from Kruskal (1968).
Christensen (2020, chapter 10) develops the consequences for the general
Gauss–Markov model, which Chapter 31 follows.

**Computation.**  Golub and Van Loan (2013) is the standard reference for QR,
the SVD and their perturbation theory. Higham (2002) gives the precise
statement of the conditioning of least squares problems, including the
\( \kappa^2\norm{\he} \) term mentioned in [Section 6.10](10-computation.html).
Seber and Lee (2003, chapter 11) surveys algorithms from a statistician's
perspective, and [Chapter 10](../ch10-computation/index.html) builds on it.

## References

- Agresti, Alan (2015). *Foundations of Linear and Generalized Linear Models*. Hoboken, NJ: Wiley.
- Aitken, Alexander C. (1935). On Least Squares and Linear Combination of Observations. *Proceedings of the Royal Society of Edinburgh* 55, 42–48.
- Christensen, Ronald (2020). *Plane Answers to Complex Questions: The Theory of Linear Models*. 5th edition. Cham: Springer.
- Frisch, Ragnar and Waugh, Frederick V. (1933). Partial Time Regressions as Compared with Individual Trends. *Econometrica* 1(4), 387–401.
- Golub, Gene H. and Van Loan, Charles F. (2013). *Matrix Computations*. 4th edition. Baltimore: Johns Hopkins University Press.
- Halmos, Paul R. (1958). *Finite-Dimensional Vector Spaces*. 2nd edition. Princeton, NJ: Van Nostrand.
- Herr, David G. (1980). On the History of the Use of Geometry in the General Linear Model. *The American Statistician* 34(1), 43–47.
- Higham, Nicholas J. (2002). *Accuracy and Stability of Numerical Algorithms*. 2nd edition. Philadelphia: SIAM.
- Hoaglin, David C. and Welsch, Roy E. (1978). The Hat Matrix in Regression and ANOVA. *The American Statistician* 32(1), 17–22.
- Kruskal, William (1968). When Are Gauss–Markov and Least Squares Estimators Identical? A Coordinate-Free Approach. *The Annals of Mathematical Statistics* 39(1), 70–75.
- Lovell, Michael C. (1963). Seasonal Adjustment of Economic Time Series and Multiple Regression Analysis. *Journal of the American Statistical Association* 58(304), 993–1010.
- Penrose, Roger (1955). A Generalized Inverse for Matrices. *Mathematical Proceedings of the Cambridge Philosophical Society* 51(3), 406–413.
- Rao, C. Radhakrishna (1973). *Linear Statistical Inference and Its Applications*. 2nd edition. New York: Wiley.
- Saville, David J. and Wood, Graham R. (1991). *Statistical Methods: The Geometric Approach*. New York: Springer.
- Seber, George A. F. and Lee, Alan J. (2003). *Linear Regression Analysis*. 2nd edition. Hoboken, NJ: Wiley.
- Wickens, Thomas D. (1995). *The Geometry of Multivariate Statistics*. Hillsdale, NJ: Lawrence Erlbaum.
