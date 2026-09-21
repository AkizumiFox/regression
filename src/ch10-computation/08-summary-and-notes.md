# Summary and Notes

## Summary

::: {.idea}

1. Floating-point arithmetic commits a relative error of at most \( u=2^{-53} \) per operation
           (@def-cmp-standard-model). Inner products are backward stable (@prp-cmp-inner-product), but
           forming \( \X\T\X \) destroys the information in singular values below about \( u^{1/2}\sigma_1 \)
           (@prp-cmp-gram-rounding, @exm-cmp-lauchli).

2. The Cholesky factorization of the augmented cross-product matrix gives \( \hbeta \), the SSE, the
           sequential sums of squares, the leverages and the standard errors in \( np^2+p^3/3 \) flops
           (@thm-cmp-cholesky, @prp-cmp-augmented), accurately enough when the column-scaled condition number
           is modest. Sweeping fits nested models (@prp-cmp-sweep); centring needs two passes or updating
           (@prp-cmp-welford).

3. Householder QR reduces \( \X \) to a triangle by reflections in \( 2np^2-2p^3/3 \) flops and is backward
           stable with columnwise errors (@thm-cmp-householder, @prp-cmp-qr-quantities). Givens rotations work
           one row at a time, for streaming and sparse data (@thm-cmp-givens).

4. The SVD gives every least squares solution, the minimum-norm one, and the precision \( \sigma/\sigma_i \) of each
           combination \( \bv_i\T\bbeta \) (@thm-cmp-svd-ls). Singular values are perturbed by at most \( \norm{\mathbf{E}}_2 \)
           (@prp-cmp-weyl), so rank is meaningful only relative to a tolerance (@def-cmp-numerical-rank). Truncated
           SVD trades bias for variance (@prp-cmp-tsvd).

5. The least squares solution has sensitivity \( \kappa+\kappa^2\eta \), where \( \eta \) measures the residual
           relative to the fit; residuals have sensitivity \( \kappa \) only (@thm-cmp-perturbation). The normal
           equations add \( \kappa^2 \) regardless. Equilibrating the columns is nearly optimal among scalings
           (@prp-cmp-van-der-sluis, @exm-cmp-longley).

6. A triangular factor can be updated or downdated by one row in \( O(p^2) \) (@thm-cmp-updating), and by one
           column in \( O(np) \) (@prp-cmp-add-column). All leave-one-out residuals, coefficients and sums of squares
           follow from a single fit (@prp-cmp-loo).

7. QR with column pivoting orders the columns by how much each adds to the others, and a small trailing pivot
           certifies a small singular value (@thm-cmp-rank-revealing). Basic and minimum-norm solutions agree on
           everything estimable (@exm-cmp-grunfeld); the tolerance should reflect the accuracy of the data
           (@exm-cmp-tolerance).

:::

## Notes and sources

**Source coverage.** The chapter covers Seber and Lee (2003, chapter 11) on solving the normal equations,
QR by Gram–Schmidt, Householder and Givens, the SVD, centring and leverages (sections 11.2–11.4, 11.7
and 11.10), with parts of sections 11.6, 11.8 and 11.9 on updating, comparing methods and rank
deficiency.

**Floating point and stability.** The standard model and the \( \gamma_n \) notation follow Higham (2002), whose
chapters 3, 10, 19 and 20 are the definitive treatment of inner products, Cholesky, QR and least squares. Goldberg
(1991) is a readable introduction to IEEE arithmetic, now specified by IEEE Std 754-2019. The distinction
between a problem's condition and an algorithm's backward stability goes back to Wilkinson (1963). Trefethen and Bau (1997) give a short modern account.
Läuchli (1961) introduced the matrix of
@exm-cmp-lauchli to show the loss in forming normal equations.

**Algorithms.** Householder (1958) introduced the reflections and Golub (1965) applied them to least
squares; Givens (1958) used plane rotations. Gentleman (1973) developed square-root-free rotations and the
row-wise algorithm, turned into widely used routines by Miller (1992). Björck (1967) analysed modified
Gram–Schmidt, and Björck and Paige (1992) explained its good behaviour for least squares. The SVD algorithm
is due to Golub and Kahan (1965) and Golub and Reinsch (1970), column pivoting to Businger and Golub (1965),
and the counterexample of @exm-cmp-kahan to Kahan (1966). Chan (1987) gave the first rank-revealing QR algorithm, Hong and Pan (1992)
proved that factorizations with polynomial bounds exist, and Gu and Eisenstat (1996) showed how to
compute them efficiently. Goodnight (1979) describes the sweep operator for
statisticians. Welford (1962) gave the updating recursion, and Chan, Golub and LeVeque (1983) compare
variance algorithms. Downdating appears in LINPACK (Dongarra et al. 1979) and was analysed by Stewart
(1979). The standard monographs are Björck (1996), Lawson and Hanson (1974) and Golub and Van Loan (2013);
Thisted (1988) and Kennedy and Gentle (1980) take the statistician's side.

**Perturbation theory.** Golub and Wilkinson (1966) first observed that the least squares solution can be
sensitive to \( \kappa^2 \) when the residual is large. Wedin (1973) gave sharp perturbation bounds for pseudoinverses
and least squares. The equilibration result of @prp-cmp-van-der-sluis is from van der Sluis (1969). The inequality
of @prp-cmp-weyl is the singular-value form of Weyl's (1912) inequality for eigenvalues. Beaton, Rubin and Barone (1976) argued that errors from rounding the data
should be the yardstick for computational accuracy, the point made at the end of
[Section 10.5](05-conditioning.html).

**Statistics.** Recursive least squares via the Sherman–Morrison formula goes back to Plackett (1950). The PRESS
statistic is due to Allen (1974). The recursive residuals of @exr-cmp-recursive-residuals are those of Brown, Durbin
and Evans (1975). Deletion diagnostics are developed in [Chapter 20](../ch20-residuals-leverage-influence/index.html), collinearity in [Chapter 26](../ch26-collinearity/index.html), and ridge and
principal components regression in [Chapter 27](../ch27-shrinkage/index.html).

## References

- Allen, David M. (1974). The Relationship Between Variable Selection and Data Augmentation and a Method for Prediction. *Technometrics* 16(1), 125–127.
- Beaton, Albert E., Rubin, Donald B. and Barone, John L. (1976). The Acceptability of Regression Solutions: Another Look at Computational Accuracy. *Journal of the American Statistical Association* 71(353), 158–168.
- Björck, Åke (1967). Solving Linear Least Squares Problems by Gram–Schmidt Orthogonalization. *BIT* 7(1), 1–21.
- Björck, Åke (1996). *Numerical Methods for Least Squares Problems*. Philadelphia: SIAM.
- Björck, Åke and Paige, Christopher C. (1992). Loss and Recapture of Orthogonality in the Modified Gram–Schmidt Algorithm. *SIAM Journal on Matrix Analysis and Applications* 13(1), 176–190.
- Brown, R. L., Durbin, James and Evans, J. M. (1975). Techniques for Testing the Constancy of Regression Relationships over Time. *Journal of the Royal Statistical Society, Series B* 37(2), 149–192.
- Businger, Peter and Golub, Gene H. (1965). Linear Least Squares Solutions by Householder Transformations. *Numerische Mathematik* 7(3), 269–276.
- Chan, Tony F. (1987). Rank Revealing QR Factorizations. *Linear Algebra and Its Applications* 88/89, 67–82.
- Chan, Tony F., Golub, Gene H. and LeVeque, Randall J. (1983). Algorithms for Computing the Sample Variance: Analysis and Recommendations. *The American Statistician* 37(3), 242–247.
- Dongarra, Jack J., Bunch, James R., Moler, Cleve B. and Stewart, Gilbert W. (1979). *LINPACK Users' Guide*. Philadelphia: SIAM.
- Gentleman, W. Morven (1973). Least Squares Computations by Givens Transformations Without Square Roots. *Journal of the Institute of Mathematics and Its Applications* 12(3), 329–336.
- Givens, Wallace (1958). Computation of Plane Unitary Rotations Transforming a General Matrix to Triangular Form. *Journal of the Society for Industrial and Applied Mathematics* 6(1), 26–50.
- Goldberg, David (1991). What Every Computer Scientist Should Know About Floating-Point Arithmetic. *ACM Computing Surveys* 23(1), 5–48.
- Golub, Gene H. (1965). Numerical Methods for Solving Linear Least Squares Problems. *Numerische Mathematik* 7(3), 206–216.
- Golub, Gene H. and Kahan, William (1965). Calculating the Singular Values and Pseudo-Inverse of a Matrix. *Journal of the Society for Industrial and Applied Mathematics, Series B: Numerical Analysis* 2(2), 205–224.
- Golub, Gene H. and Reinsch, Christian (1970). Singular Value Decomposition and Least Squares Solutions. *Numerische Mathematik* 14(5), 403–420.
- Golub, Gene H. and Van Loan, Charles F. (2013). *Matrix Computations*. 4th edition. Baltimore: Johns Hopkins University Press.
- Golub, Gene H. and Wilkinson, James H. (1966). Note on the Iterative Refinement of Least Squares Solution. *Numerische Mathematik* 9(2), 139–148.
- Goodnight, James H. (1979). A Tutorial on the SWEEP Operator. *The American Statistician* 33(3), 149–158.
- Gu, Ming and Eisenstat, Stanley C. (1996). Efficient Algorithms for Computing a Strong Rank-Revealing QR Factorization. *SIAM Journal on Scientific Computing* 17(4), 848–869.
- Higham, Nicholas J. (2002). *Accuracy and Stability of Numerical Algorithms*. 2nd edition. Philadelphia: SIAM.
- Hong, Yoo Pyo and Pan, Ching-Tsuan (1992). Rank-Revealing QR Factorizations and the Singular Value Decomposition. *Mathematics of Computation* 58(197), 213–232.
- Householder, Alston S. (1958). Unitary Triangularization of a Nonsymmetric Matrix. *Journal of the ACM* 5(4), 339–342.
- IEEE (2019). *IEEE Standard for Floating-Point Arithmetic* (IEEE Std 754-2019). New York: IEEE.
- Kahan, William (1966). Numerical Linear Algebra. *Canadian Mathematical Bulletin* 9(5), 757–801.
- Kennedy, William J. and Gentle, James E. (1980). *Statistical Computing*. New York: Marcel Dekker.
- Läuchli, Peter (1961). Jordan-Elimination und Ausgleichung nach kleinsten Quadraten. *Numerische Mathematik* 3(1), 226–240.
- Lawson, Charles L. and Hanson, Richard J. (1974). *Solving Least Squares Problems*. Englewood Cliffs, NJ: Prentice-Hall. Reprinted Philadelphia: SIAM, 1995.
- Longley, James W. (1967). An Appraisal of Least Squares Programs for the Electronic Computer from the Point of View of the User. *Journal of the American Statistical Association* 62(319), 819–841.
- Miller, Alan J. (1992). Algorithm AS 274: Least Squares Routines to Supplement Those of Gentleman. *Journal of the Royal Statistical Society, Series C* 41(2), 458–478.
- Plackett, Robin L. (1950). Some Theorems in Least Squares. *Biometrika* 37(1/2), 149–157.
- Seber, George A. F. and Lee, Alan J. (2003). *Linear Regression Analysis*. 2nd edition. Hoboken, NJ: Wiley.
- Stewart, Gilbert W. (1979). The Effects of Rounding Error on an Algorithm for Downdating a Cholesky Factorization. *Journal of the Institute of Mathematics and Its Applications* 23(2), 203–213.
- Thisted, Ronald A. (1988). *Elements of Statistical Computing: Numerical Computation*. New York: Chapman and Hall.
- Trefethen, Lloyd N. and Bau, David (1997). *Numerical Linear Algebra*. Philadelphia: SIAM.
- van der Sluis, A. (1969). Condition Numbers and Equilibration of Matrices. *Numerische Mathematik* 14(1), 14–23.
- Wedin, Per-Åke (1973). Perturbation Theory for Pseudo-Inverses. *BIT* 13(2), 217–232.
- Welford, B. P. (1962). Note on a Method for Calculating Corrected Sums of Squares and Products. *Technometrics* 4(3), 419–420.
- Weyl, Hermann (1912). Das asymptotische Verteilungsgesetz der Eigenwerte linearer partieller Differentialgleichungen. *Mathematische Annalen* 71(4), 441–479.
- Wilkinson, James H. (1963). *Rounding Errors in Algebraic Processes*. London: Her Majesty's Stationery Office.
