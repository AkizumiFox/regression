# Summary and Notes

## Summary

::: {.idea}

1. Row rank equals column rank, and \( \rank(\A)+\dim\Null(\A) \) is the number of
           columns. \( \C(\A\T\A)=\C(\A\T) \), so the normal equations are always consistent
           (@thm-mat-row-col-rank, @thm-mat-rank-nullity and @prp-mat-rank-product).

2. A generalized inverse satisfies \( \A\A\ginv\A=\A \). Every matrix has one, built from any
           nonsingular \( r\times r \) submatrix. The Moore–Penrose inverse is the unique one that is
           also reflexive and makes \( \A\A^+ \) and \( \A^+\A \) symmetric
           (@def-mat-ginverse, @thm-mat-ginverse-exists, @prp-mat-nonsingular-submatrix and @def-mat-moore-penrose).

3. \( \A\x=\bb \) is consistent iff \( \A\A\ginv\bb=\bb \), and then its solutions are
           \( \A\ginv\bb+(\I-\A\ginv\A)\bz \). A linear function \( \bm q\T\x \) is the same for all
           solutions iff \( \bm q\in\C(\A\T) \) (@thm-mat-consistency and @cor-mat-invariant).

4. Orthogonal matrices preserve lengths and angles. Gram–Schmidt gives the QR
           factorization, and \( \Real^n=\C(\A\T)\dirsum\Null(\A) \) (@prp-mat-orthogonal, @thm-mat-qr and @prp-mat-complement).

5. \( \tr(\A\B)=\tr(\B\A) \). Block elimination gives the determinant and inverse of a
           partitioned matrix through Schur complements, and the Sherman–Morrison–Woodbury
           identity as a special case
           (@thm-mat-trace-cyclic, @thm-mat-block-determinant, @thm-mat-partitioned-inverse and @thm-mat-woodbury).

6. A symmetric matrix is \( \Q\bLambda\Q\T \) with \( \Q \) orthogonal. Positive definiteness is
           equivalent to positive eigenvalues, to a Cholesky factorization, and to positive leading
           minors. Nonnegative definite matrices have a unique nonnegative definite square root
           (@thm-mat-spectral, @thm-mat-pd-characterizations and @thm-mat-square-root).

7. The Rayleigh quotient \( \x\T\A\x/\x\T\x \) ranges over \( [\lambda_{\min},\lambda_{\max}] \)
           (@thm-mat-extremal-rayleigh and @cor-mat-generalized-rayleigh).

8. Every matrix is \( \bU\bD\V\T \). The condition number \( \kappa=\sigma_{\max}/\sigma_{\min} \)
           bounds error growth, is squared by forming \( \X\T\X \), and its smallest singular vector
           exhibits near-dependences among the columns
           (@thm-mat-svd, @def-mat-condition-number and @prp-mat-svd-norms).

9. Idempotent matrices have eigenvalues \( 0 \) and \( 1 \) and rank equal to trace. Symmetric
           ones are \( \Q_1\Q_1\T \). If symmetric matrices summing to \( \I \) have ranks summing to \( n \), they are
           idempotent and mutually orthogonal
           (@prp-mat-idempotent-basic, @thm-mat-idempotent and @thm-mat-idempotent-sum).

10. Kronecker products factor balanced designs. \( \vecop(\A\B\bm C)=(\bm C\T\otimes\A)\vecop(\B) \)
            (@prp-mat-kronecker).

11. \( \partial(\x\T\A\x)/\partial\x=(\A+\A\T)\x \). Quadratic functions with nonnegative definite
            matrices are minimized by solving a linear system, with or without linear constraints
            (@thm-mat-quadform-derivative, @prp-mat-quadratic-min and @prp-mat-lagrange).

:::

## Notes and sources

**Coverage and framing.**  The selection of results follows what the rest of the book
uses. Rencher and Schaalje (2008, chapter 2) was the checklist for topics: its sections on
positive definite matrices, generalized inverses, orthogonal matrices, idempotent matrices and
matrix calculus correspond to [Section 1.7](07-eigen.html), [Section 1.3](03-inverses.html), [Section 1.5](05-orthogonality.html), [Section 1.9](09-idempotent.html) and [Section 1.11](11-calculus.html).
Several results that Rencher and Schaalje state without proof or with a reference, such as the
general solution of a consistent system, the characterizations of positive definiteness and
the sums of idempotent matrices, are proved here. The vector-space point of view, in which
column spaces, null spaces and orthogonal complements come first and matrices second, is that
of Christensen (2020, appendices A and B). Seber and Lee (2003, appendix A)
is a compact list of further matrix facts. The Kronecker product results in
[Section 1.10](10-kronecker.html) are also in Christensen's appendix B.

**Linear algebra for statisticians.**  Harville (1997) is the most
complete treatment written for statisticians, with full proofs of the determinant facts
assumed in @prp-mat-det and a long development of generalized inverses.
Searle (1982) is more elementary and has many worked examples from linear models.
Horn and Johnson (2013) is the standard reference for eigenvalue inequalities, positive
definite matrices and norms. Halmos (1958) gives the coordinate-free view of
the vector-space facts used in [Section 1.2](02-rank.html).

**Generalized inverses.**  The four conditions of @def-mat-moore-penrose are due to
Penrose (1955), who proved existence and uniqueness. E. H. Moore had
described an equivalent inverse decades earlier. Rao (1973) made the weaker
inverse \( \A\A\ginv\A=\A \) the basic tool of linear-model theory, and it is his approach that
[Section 1.4](04-systems.html) and [Chapter 8](../ch08-estimability/index.html) follow. The construction from a nonsingular
submatrix in @thm-mat-ginverse-exists is a matrix-free way of stating an algorithm
that is often given as a sequence of steps (Searle 1982).

**Schur complements and updating formulas.**  The name “Schur complement” and the
notation for it vary. Zhang (2005) collects the history and many applications.
The rank-one formula @eq-mat-sherman-morrison is from Sherman and Morrison (1950),
and the general identity @eq-mat-woodbury from Woodbury (1950).
Henderson and Searle (1981) trace earlier appearances and compare several derivations.
Obtaining the identity by comparing two expressions for one block of a partitioned inverse,
as in @thm-mat-woodbury, is a classical route.

**Singular values and conditioning.**  Stewart (1993) tells the early history of the singular value decomposition,
beginning with Beltrami and Jordan in the 1870s. The low-rank approximation theorem,
@prp-mat-svd-norms(c) and its Frobenius form, is usually credited to
Eckart and Young (1936), although Stewart points out that Schmidt had proved it
earlier in the setting of integral equations. Golub and Van Loan (2013) is the reference for computing the decomposition and for perturbation bounds such as
@eq-mat-perturbation, and Higham (2002) for their refinements. The use of
scaled condition numbers and singular vectors to diagnose collinearity is due to
Belsley et al. (1980). Chapter 26 returns to it.

**Kronecker products and calculus.**  Henderson and Searle (1979) survey the vec operator
and its uses in multivariate statistics. Magnus and Neudecker (2019) develop matrix calculus
systematically through differentials, which avoid the bookkeeping of the entrywise
definitions used in [Section 1.11](11-calculus.html) and handle symmetric matrix arguments cleanly.

## References

- Belsley, David A., Kuh, Edwin and Welsch, Roy E. (1980). *Regression Diagnostics: Identifying Influential Data and Sources of Collinearity*. New York: Wiley.
- Christensen, Ronald (2020). *Plane Answers to Complex Questions: The Theory of Linear Models*. 5th edition. Cham: Springer.
- Eckart, Carl and Young, Gale (1936). The Approximation of One Matrix by Another of Lower Rank. *Psychometrika* 1(3), 211–218.
- Golub, Gene H. and Van Loan, Charles F. (2013). *Matrix Computations*. 4th edition. Baltimore: Johns Hopkins University Press.
- Halmos, Paul R. (1958). *Finite-Dimensional Vector Spaces*. 2nd edition. Princeton, NJ: Van Nostrand.
- Harville, David A. (1997). *Matrix Algebra From a Statistician's Perspective*. New York: Springer.
- Henderson, Harold V. and Searle, Shayle R. (1979). Vec and Vech Operators for Matrices, with Some Uses in Jacobians and Multivariate Statistics. *Canadian Journal of Statistics* 7(1), 65–81.
- Henderson, Harold V. and Searle, Shayle R. (1981). On Deriving the Inverse of a Sum of Matrices. *SIAM Review* 23(1), 53–60.
- Higham, Nicholas J. (2002). *Accuracy and Stability of Numerical Algorithms*. 2nd edition. Philadelphia: SIAM.
- Horn, Roger A. and Johnson, Charles R. (2013). *Matrix Analysis*. 2nd edition. Cambridge: Cambridge University Press.
- Magnus, Jan R. and Neudecker, Heinz (2019). *Matrix Differential Calculus with Applications in Statistics and Econometrics*. 3rd edition. Hoboken, NJ: Wiley.
- Penrose, Roger (1955). A Generalized Inverse for Matrices. *Mathematical Proceedings of the Cambridge Philosophical Society* 51(3), 406–413.
- Rao, C. Radhakrishna (1973). *Linear Statistical Inference and Its Applications*. 2nd edition. New York: Wiley.
- Rencher, Alvin C. and Schaalje, G. Bruce (2008). *Linear Models in Statistics*. 2nd edition. Hoboken, NJ: Wiley.
- Searle, Shayle R. (1982). *Matrix Algebra Useful for Statistics*. New York: Wiley.
- Seber, George A. F. and Lee, Alan J. (2003). *Linear Regression Analysis*. 2nd edition. Hoboken, NJ: Wiley.
- Sherman, Jack and Morrison, Winifred J. (1950). Adjustment of an Inverse Matrix Corresponding to a Change in One Element of a Given Matrix. *Annals of Mathematical Statistics* 21(1), 124–127.
- Stewart, G. W. (1993). On the Early History of the Singular Value Decomposition. *SIAM Review* 35(4), 551–566.
- Woodbury, Max A. (1950). *Inverting Modified Matrices*. Princeton, NJ: .
- Zhang, Fuzhen (2005). *The Schur Complement and Its Applications*. New York: Springer.
