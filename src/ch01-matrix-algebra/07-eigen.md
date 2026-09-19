# Eigenvalues, the spectral theorem and positive definiteness

## Eigenvalues

A scalar \( \lambda \) is an **eigenvalue** of the \( n\times n \) matrix \( \A \), with
**eigenvector** \( \x\neq\bzero \), if \( \A\x=\lambda\x \). Equivalently, \( \A-\lambda\I \) is
singular, that is, \( \lambda \) is a root of the **characteristic polynomial**
\( p(\lambda)=\det(\lambda\I-\A) \). This polynomial has degree \( n \) and leading coefficient
one, so over the complex numbers it factors as \( \prod_{i=1}^n(\lambda-\lambda_i) \). The
eigenvalues \( \lambda_1,\dots,\lambda_n \) are listed with multiplicity and may be complex
when \( \A \) is not symmetric.

::: {#prp-mat-eigen-basic}
Let \( \A \) be \( n\times n \) with eigenvalues \( \lambda_1,\dots,\lambda_n \).

::: {.enumerate options="label=(\alph*)"}
1. \( \det\A=\prod_i\lambda_i \) and \( \tr\A=\sum_i\lambda_i \).

2. If \( \A\x=\lambda\x \), then \( f(\A)\x=f(\lambda)\x \) for every polynomial \( f \). If
           \( \A \) is nonsingular, \( \A^{-1}\x=\lambda^{-1}\x \).

3. \( \A \) and \( \bm P^{-1}\A\bm P \) have the same characteristic polynomial for nonsingular \( \bm P \).

4. If \( \bm C \) is \( m\times n \) and \( \bm D \) is \( n\times m \), then
           \( \lambda^n\det(\lambda\I_m-\bm C\bm D)=\lambda^m\det(\lambda\I_n-\bm D\bm C) \). So \( \bm C\bm D \) and
           \( \bm D\bm C \) have the same nonzero eigenvalues, with the same multiplicities.
:::

:::

::: {.proof}
(a) Setting \( \lambda=0 \) in \( \det(\lambda\I-\A)=\prod(\lambda-\lambda_i) \) gives
\( (-1)^n\det\A=(-1)^n\prod\lambda_i \). For the trace, compare coefficients of
\( \lambda^{n-1} \). On the right it is \( -\sum\lambda_i \). On the left, in the permutation
expansion of \( \det(\lambda\I-\A) \), only the product of diagonal entries
\( \prod_i(\lambda-a_{ii}) \) contains \( \lambda^{n-1} \), because every other term omits at least
two diagonal entries. Its coefficient is \( -\sum a_{ii} \).
(b) \( \A^k\x=\lambda^k\x \) by induction, and \( \x=\A^{-1}\A\x=\lambda\A^{-1}\x \) with
\( \lambda\ne0 \). (c) \( \det(\lambda\I-\bm P^{-1}\A\bm P)=\det(\bm P^{-1}(\lambda\I-\A)\bm P)=\det(\lambda\I-\A) \).
(d) For \( \lambda\ne0 \), @thm-mat-block-determinant(c) gives
\( \det(\I_m-\lambda^{-1}\bm C\bm D)=\det(\I_n-\lambda^{-1}\bm D\bm C) \). Multiplying by \( \lambda^{m+n} \)
gives the identity for \( \lambda\ne0 \). Both sides are polynomials, so it holds for all
\( \lambda \).
:::

## The spectral theorem

Symmetric matrices have the best possible eigenstructure, and almost every matrix whose
eigenvalues matter in statistics is symmetric.

::: {#thm-mat-spectral}
[Spectral theorem]

Let \( \A \) be an \( n\times n \) symmetric matrix. Then all eigenvalues of \( \A \) are real, and
there is an orthogonal matrix \( \Q=[\bm q_1,\dots,\bm q_n] \) with
\[
\A=\Q\bLambda\Q\T=\sum_{i=1}^n\lambda_i\bm q_i\bm q_i\T,\qquad
\bLambda=\diag(\lambda_1,\dots,\lambda_n).
\]{#eq-mat-spectral}

The columns of \( \Q \) are an orthonormal basis of eigenvectors, \( \A\bm q_i=\lambda_i\bm q_i \).
Eigenvectors belonging to distinct eigenvalues are orthogonal.
:::

::: {.proof}
*Real eigenvalues.* Let \( \A\bz=\lambda\bz \) with \( \bz\in\mathbb C^n \) nonzero, and let
\( \bar{\bz} \) be the complex conjugate. The scalar \( s=\bar{\bz}\T\A\bz \) is real, since its
conjugate is \( \bz\T\A\bar{\bz}=(\bz\T\A\bar{\bz})\T=\bar{\bz}\T\A\bz=s \) by symmetry. Also
\( s=\lambda\bar{\bz}\T\bz=\lambda\sum_i\lvert z_i\rvert^2 \), with the sum real and positive. So
\( \lambda \) is real. The real singular matrix \( \A-\lambda\I \) then has a real nonzero null
vector, so there is a real eigenvector.

*Diagonalization.* By induction on \( n \), the case \( n=1 \) being trivial. Let
\( \bm q_1 \) be a real unit eigenvector for \( \lambda_1 \) and extend it to an orthogonal matrix
\( [\bm q_1,\Q_2] \) ([Section 1.5](05-orthogonality.html)). Since
\( \Q_2\T\A\bm q_1=\lambda_1\Q_2\T\bm q_1=\bzero \) and \( \A \) is symmetric,
\[
[\bm q_1,\Q_2]\T\A[\bm q_1,\Q_2]=\begin{pmatrix}\lambda_1&\bzero\T\\\bzero&\Q_2\T\A\Q_2\end{pmatrix}.
\]
The \( (n-1)\times(n-1) \) block is symmetric, so by induction
\( \Q_2\T\A\Q_2=\bm P\bLambda_2\bm P\T \) with \( \bm P \) orthogonal. Then
\( \Q=[\bm q_1,\Q_2\bm P] \) is orthogonal (@prp-mat-orthogonal(b)) and
\( \Q\T\A\Q=\diag(\lambda_1,\bLambda_2) \).

*Orthogonality.* If \( \A\bu=\lambda\bu \) and \( \A\bv=\mu\bv \), then
\( \lambda\bu\T\bv=\bu\T\A\bv=\mu\bu\T\bv \), so \( \bu\T\bv=0 \) when \( \lambda\ne\mu \).
:::

Unless stated otherwise, the eigenvalues of symmetric \( \A \) are ordered as
\( \lambda_1\ge\lambda_2\ge\dots\ge\lambda_n \), and we write
\( \lambda_{\max}=\lambda_1 \) and \( \lambda_{\min}=\lambda_n \). The following consequences are
used without comment later on.

::: {#cor-mat-spectral-consequences}
Let \( \A=\Q\bLambda\Q\T \) be symmetric.

::: {.enumerate options="label=(\alph*)"}
1. \( \rank(\A) \) counts the nonzero eigenvalues. \( \C(\A) \) is spanned by the
           \( \bm q_i \) with \( \lambda_i\ne0 \), and \( \Null(\A) \) by those with \( \lambda_i=0 \).

2. \( \A^k=\Q\bLambda^k\Q\T \) for \( k=1,2,\dots \). If \( \A \) is nonsingular,
           \( \A^{-1}=\Q\bLambda^{-1}\Q\T \), and \( \A^+=\Q\bLambda^+\Q\T \) in general, where
           \( \bLambda^+ \) inverts the nonzero diagonal entries.

3. The eigenspace \( \{\x:\A\x=\lambda\x\} \) has dimension equal to the multiplicity of
           \( \lambda \) as a root of the characteristic polynomial.
:::

:::

::: {.proof}
Multiplying by the nonsingular \( \Q \) and \( \Q\T \) changes neither rank
(@prp-mat-rank-product(b)) nor, after the change of variable \( \x=\Q\bm c \), the
structure of null spaces and eigenspaces. So each statement reduces to the diagonal
matrix \( \bLambda \), where it is immediate. For (b), \( \Q\T\Q=\I \) collapses the middle factors,
and the four Penrose conditions for \( \Q\bLambda^+\Q\T \) reduce to those for \( \bLambda^+ \).
:::

::: {#exm-mat-equicorrelation-eigen}
[Eigenvalues of the equicorrelation matrix]

For \( \bm E=(1-\rho)\I+\rho\bm J_n \) of @exm-mat-equicorrelation,
\( \bm E\bone=(1+(n-1)\rho)\bone \), and \( \bm E\bv=(1-\rho)\bv \) for every \( \bv\perp\bone \). So the
eigenvalues are \( 1+(n-1)\rho \) once and \( 1-\rho \) with multiplicity \( n-1 \), with
eigenvectors \( \bone/\sqrt n \) and the last \( n-1 \) rows of the Helmert matrix
(@exm-mat-helmert). The determinant and inverse found earlier follow at once from
@prp-mat-eigen-basic and @cor-mat-spectral-consequences(b).
:::

## Positive definite and nonnegative definite matrices

::: {#def-mat-nnd}
A symmetric matrix \( \A \) is **positive definite** if \( \x\T\A\x>0 \) for every
\( \x\ne\bzero \), and **nonnegative definite** if \( \x\T\A\x\ge0 \) for every \( \x \). For
symmetric \( \A,\B \) of the same size we write \( \A\succeq\B \) if \( \A-\B \) is nonnegative
definite, and \( \A\succ\B \) if it is positive definite (the **Loewner order**).
:::

In this book, “nonnegative definite” includes positive definite. Some texts,
Rencher and Schaalje among them, reserve “positive semidefinite” for matrices that are
nonnegative definite but singular. We use the terms only for symmetric matrices.
Covariance matrices are nonnegative definite ([Chapter 2](../ch02-random-vectors/index.html)), and
*Gram matrices* \( \X\T\X \) are too, since \( \x\T\X\T\X\x=\norm{\X\x}^2 \).

::: {#thm-mat-pd-characterizations}
[Characterizations]

For a symmetric \( n\times n \) matrix \( \A \) the following are equivalent:

::: {.enumerate options="label=(\alph*)"}
1. \( \A \) is positive definite;

2. every eigenvalue of \( \A \) is positive;

3. \( \A=\B\T\B \) for some matrix \( \B \) with \( n \) columns and full column rank;

4. \( \A=\bL\bL\T \) with \( \bL \) lower triangular with positive diagonal entries
           (the **Cholesky factorization**; \( \bL \) is unique);

5. every leading principal minor \( \det\A_{[k]} \), \( k=1,\dots,n \), is positive, where
           \( \A_{[k]} \) is the upper left \( k\times k \) block.
:::

Likewise, the following are equivalent:
(a\( ' \)) \( \A \) is nonnegative definite; (b\( ' \)) every eigenvalue is nonnegative;
(c\( ' \)) \( \A=\B\T\B \) for some \( \B \), which can be taken of size \( \rank(\A)\times n \);
(e\( ' \)) every principal minor (determinant of a submatrix with the same row and column
indices) is nonnegative.
:::

::: {.proof}
Write \( \A=\Q\bLambda\Q\T \).
(a)\( \Rightarrow \)(b): \( \lambda_i=\bm q_i\T\A\bm q_i>0 \).
(b)\( \Rightarrow \)(c): \( \B=\bLambda^{1/2}\Q\T \) is nonsingular.
(c)\( \Rightarrow \)(a): \( \x\T\B\T\B\x=\norm{\B\x}^2 \), which is positive for \( \x\neq\bzero \) by full
column rank. (d)\( \Rightarrow \)(c) is immediate, with \( \B=\bL\T \).

(a)\( \Rightarrow \)(d), by induction on \( n \). Write
\( \A=\begin{psmallmatrix}a&\bb\T\\\bb&\A_2\end{psmallmatrix} \). Then \( a=\bm e_1\T\A\bm e_1>0 \).
The Schur complement \( \A_2-\bb\bb\T/a \) is positive definite: for \( \bv\ne\bzero \), put
\( \x=(-\bb\T\bv/a,\ \bv\T)\T \) and expand to get
\( \x\T\A\x=\bv\T(\A_2-\bb\bb\T/a)\bv>0 \). By induction
\( \A_2-\bb\bb\T/a=\bL_2\bL_2\T \), and
\( \bL=\begin{psmallmatrix}\sqrt a&\bzero\T\\\bb/\sqrt a&\bL_2\end{psmallmatrix} \) satisfies
\( \bL\bL\T=\A \). For uniqueness, if \( \bL\bL\T=\bm K\bm K\T \), then
\( \bm K^{-1}\bL=\bm K\T(\bL\T)^{-1} \) is both lower and upper triangular, hence diagonal, say
\( \bD \). From \( \bm K^{-1}\bL=\bD \) and \( \bL\bL\T=\bm K\bm K\T \) we get \( \bD\bD\T=\I \), so \( \bD=\I \)
because its diagonal is positive.

(a)\( \Rightarrow \)(e): \( \A_{[k]}=\bm E\T\A\bm E \) with \( \bm E=[\bm e_1,\dots,\bm e_k] \) is positive definite,
since \( \x\T\A_{[k]}\x=(\bm E\x)\T\A(\bm E\x) \). By (b) and @prp-mat-eigen-basic(a) its
determinant is positive.
(e)\( \Rightarrow \)(a), by induction on \( n \). The block \( \A_{[n-1]} \) satisfies (e) and so is positive
definite, so by (a)\( \Rightarrow \)(d), which is already proved, \( \A_{[n-1]}=\bL\bL\T \). Write \( \A=\begin{psmallmatrix}\A_{[n-1]}&\bm c\\\bm c\T&d\end{psmallmatrix} \) and
\( s=d-\bm c\T\A_{[n-1]}^{-1}\bm c \). By @thm-mat-block-determinant(b),
\( s=\det\A/\det\A_{[n-1]}>0 \). With
\( \bT=\begin{psmallmatrix}\I&\A_{[n-1]}^{-1}\bm c\\\bzero\T&1\end{psmallmatrix} \), a direct check gives
\( \A=\bT\T\diag(\A_{[n-1]},s)\bT=\B\T\B \) with \( \B=\diag(\bL\T,\sqrt s)\,\bT \) nonsingular,
so \( \A \) is positive definite by (c)\( \Rightarrow \)(a).

The primed statements: (a\( ' \))\( \Leftrightarrow \)(b\( ' \))\( \Leftrightarrow \)(c\( ' \)) exactly as
above, taking for \( \B \) the rows of \( \bLambda^{1/2}\Q\T \) that belong to nonzero eigenvalues.
(a\( ' \))\( \Rightarrow \)(e\( ' \)) as above, with an arbitrary index set in place of \( \{1,\dots,k\} \).
For (e\( ' \))\( \Rightarrow \)(a\( ' \)), fix \( t>0 \). Each column of \( \A_{[k]}+t\I \) is a column of
\( \A_{[k]} \) plus \( t \) times a coordinate vector. Expanding the determinant by linearity in
each column gives
\( \det(\A_{[k]}+t\I)=\sum_{S}t^{\,k-\lvert S\rvert}\det\A_{SS} \), the sum over all subsets
\( S\subseteq\{1,\dots,k\} \), with \( \det\A_{\emptyset\emptyset}=1 \). Every term is nonnegative
and the term for \( S=\emptyset \) is \( t^k>0 \). So \( \A+t\I \) is positive definite by
(e)\( \Rightarrow \)(a), and \( \x\T\A\x=\lim_{t\to0}\x\T(\A+t\I)\x\ge0 \).
:::

Two warnings. Leading principal minors do not suffice for nonnegative definiteness:
\( \diag(0,-1) \) has leading minors \( 0 \) and \( 0 \). And (c\( ' \)) holds with many matrices \( \B \). The
Cholesky factor and the symmetric square root below are two convenient choices.

::: {#prp-mat-pd-properties}
Let \( \A \) be \( n\times n \) nonnegative definite.

::: {.enumerate options="label=(\alph*)"}
1. \( a_{ii}\ge0 \), and if \( a_{ii}=0 \) then row \( i \) and column \( i \) of \( \A \) are zero.

2. \( \x\T\A\x=0 \) iff \( \A\x=\bzero \).

3. \( \B\T\A\B \) is nonnegative definite for any \( n\times q \) matrix \( \B \). If \( \A \) is positive
           definite, \( \B\T\A\B \) is positive definite iff \( \B \) has full column rank.

4. Every principal submatrix of \( \A \) is nonnegative definite, and positive definite if
           \( \A \) is.

5. If \( \A \) is positive definite, it is nonsingular and \( \A^{-1} \) is positive definite.
           If also \( \A \) is partitioned as in @eq-mat-partition, both Schur complements
           \( \A_{11\cdot2} \) and \( \A_{22\cdot1} \) are positive definite.

6. If \( \A\succeq\B\succ\bzero \), then \( \B^{-1}\succeq\A^{-1} \).
:::

:::

::: {.proof}
(b) Write \( \A=\B\T\B \) by @thm-mat-pd-characterizations(c\( ' \)). Then \( \x\T\A\x=\norm{\B\x}^2 \),
which vanishes iff \( \B\x=\bzero \), which implies \( \A\x=\bzero \).
(a) \( a_{ii}=\bm e_i\T\A\bm e_i\ge0 \). If \( a_{ii}=0 \) then \( \A\bm e_i=\bzero \) by (b), so column \( i \)
is zero, and row \( i \) by symmetry.
(c) \( \x\T\B\T\A\B\x=(\B\x)\T\A(\B\x) \) is nonnegative, and when \( \A \) is positive definite it is
zero only if \( \B\x=\bzero \).
(d) is (c) with \( \B \) a set of columns of \( \I \).
(e) The eigenvalues of \( \A^{-1} \) are \( 1/\lambda_i>0 \). The inverse in
@eq-mat-partitioned-inverse has \( \A_{22\cdot1}^{-1} \) as a principal submatrix, which is
positive definite by (d). The other Schur complement follows by symmetry of the argument.
(f) By @thm-mat-square-root below, \( \B=\B^{1/2}\B^{1/2} \). Then
\( \bm C=\B^{-1/2}\A\B^{-1/2}\succeq\I \) by (c), so every eigenvalue of \( \bm C \) is at least \( 1 \)
and every eigenvalue of \( \bm C^{-1}=\B^{1/2}\A^{-1}\B^{1/2} \) is at most \( 1 \). So
\( \I\succeq\B^{1/2}\A^{-1}\B^{1/2} \), and multiplying on both sides by \( \B^{-1/2} \) gives
\( \B^{-1}\succeq\A^{-1} \) by (c).
:::

Part (f) is used in [Chapter 7](../ch07-optimality/index.html): an estimator whose covariance matrix is smaller in
the Loewner order has smaller variance for every linear combination of its components.

::: {#thm-mat-square-root}
[Square root]

A nonnegative definite matrix \( \A=\Q\bLambda\Q\T \) has exactly one nonnegative definite
matrix \( \B \) with \( \B^2=\A \), namely
\( \A^{1/2}=\Q\bLambda^{1/2}\Q\T \). It satisfies \( \C(\A^{1/2})=\C(\A) \). If \( \A \) is positive definite,
\( \A^{1/2} \) is positive definite and its inverse is \( \A^{-1/2}=\Q\bLambda^{-1/2}\Q\T \).
:::

::: {.proof}
\( \A^{1/2} \) is symmetric with nonnegative eigenvalues and \( (\Q\bLambda^{1/2}\Q\T)^2=\Q\bLambda\Q\T \).
Its column space is spanned by the \( \bm q_i \) with \( \lambda_i>0 \), as is \( \C(\A) \)
(@cor-mat-spectral-consequences(a)). For uniqueness, let \( \B \) be nonnegative definite with
\( \B^2=\A \), and write \( \B=\bm P\bm\Delta\bm P\T \) with \( \bm\Delta=\diag(\delta_i) \), \( \delta_i\ge0 \).
Then \( \A=\bm P\bm\Delta^2\bm P\T \). Let \( \A\bv=\lambda\bv \) and expand \( \bv=\sum_ic_i\bm p_i \).
Comparing \( \A\bv=\sum_ic_i\delta_i^2\bm p_i \) with \( \lambda\bv \) gives \( c_i(\delta_i^2-\lambda)=0 \), so
\( c_i\ne0 \) only if \( \delta_i=\sqrt\lambda \). Hence \( \B\bv=\sum_ic_i\delta_i\bm p_i=\sqrt\lambda\,\bv \).
Applying this to each \( \bm q_j \) gives \( \B\Q=\Q\bLambda^{1/2} \), that is, \( \B=\A^{1/2} \).
The last statement is @cor-mat-spectral-consequences(b).
:::

## Extremal properties of eigenvalues

::: {#thm-mat-extremal-rayleigh}
[Rayleigh quotient]

Let \( \A \) be symmetric with eigenvalues \( \lambda_1\ge\dots\ge\lambda_n \) and orthonormal
eigenvectors \( \bm q_1,\dots,\bm q_n \). For every \( \x\ne\bzero \),
\[
\lambda_n\le\frac{\x\T\A\x}{\x\T\x}\le\lambda_1 ,
\]{#eq-mat-rayleigh}

with equality on the right at \( \x=\bm q_1 \) and on the left at \( \x=\bm q_n \). More generally, the
maximum of \( \x\T\A\x/\x\T\x \) over nonzero \( \x \) orthogonal to \( \bm q_1,\dots,\bm q_{k-1} \) is
\( \lambda_k \), attained at \( \bm q_k \).
:::

::: {.proof}
Put \( \bm c=\Q\T\x \), so that \( \x\T\x=\bm c\T\bm c \) and
\[
\frac{\x\T\A\x}{\x\T\x}=\frac{\sum_i\lambda_ic_i^2}{\sum_ic_i^2},
\]
a weighted average of the eigenvalues with weights \( c_i^2/\sum_jc_j^2 \). It lies between the
smallest and largest eigenvalue, and equals \( \lambda_1 \) when all weight is on \( c_1 \), that is,
at \( \x=\bm q_1 \). The condition \( \x\perp\bm q_1,\dots,\bm q_{k-1} \) means \( c_1=\dots=c_{k-1}=0 \),
so the average is then over \( \lambda_k,\dots,\lambda_n \), whose largest member is \( \lambda_k \).
:::

[Figure 1.7.1](07-eigen.html#fig-mat-rayleigh) shows the theorem for
\( \A=\begin{psmallmatrix}2.5&1\\1&1.5\end{psmallmatrix} \). The eigenvalues are
\( 3.118 \) and \( 0.882 \), with eigenvectors at angles
\( 31.7^\circ \) and \( 121.7^\circ \). As \( \x \) turns around the unit circle,
the quotient swings between the two eigenvalues. The ellipse \( \x\T\A\x=1 \) has its
principal axes along the eigenvectors, with semi-axis lengths \( \lambda_i^{-1/2} \).
Consequently the Rayleigh quotient is largest in the direction in which the ellipse is
narrowest.

::: {when-format="html"}
![**Figure 1.7.1.** The Rayleigh quotient of a \( 2\times2 \) symmetric matrix. (a) The level set
\( \x\T\A\x=1 \), an ellipse with axes along the orthonormal eigenvectors \( \bm q_1,\bm q_2 \),
and the unit circle. (b) The quotient \( \bu\T\A\bu \) along the unit circle. Its maximum
\( \lambda_1 \) and minimum \( \lambda_2 \) are attained at the eigenvector directions
(@thm-mat-extremal-rayleigh).](rayleigh.svg){#fig-mat-rayleigh width=100%}
:::

::: {when-format="pdf"}
![The Rayleigh quotient of a \( 2\times2 \) symmetric matrix. (a) The level set
\( \x\T\A\x=1 \), an ellipse with axes along the orthonormal eigenvectors \( \bm q_1,\bm q_2 \),
and the unit circle. (b) The quotient \( \bu\T\A\bu \) along the unit circle. Its maximum
\( \lambda_1 \) and minimum \( \lambda_2 \) are attained at the eigenvector directions
(@thm-mat-extremal-rayleigh).](rayleigh.pdf){width=100%}
:::

::: {#cor-mat-generalized-rayleigh}
[Relative to a positive definite matrix]

Let \( \A \) be symmetric and \( \B \) positive definite, both \( n\times n \).

::: {.enumerate options="label=(\alph*)"}
1. \( \max_{\x\ne\bzero}\x\T\A\x/\x\T\B\x \) is the largest eigenvalue of \( \B^{-1}\A \),
           which has the same eigenvalues as the symmetric \( \B^{-1/2}\A\B^{-1/2} \).

2. There is a nonsingular \( \bm P \) with \( \bm P\T\B\bm P=\I \) and \( \bm P\T\A\bm P \) diagonal.

3. For \( \bm a\in\Real^n \),
           \( \max_{\x\ne\bzero}(\bm a\T\x)^2/\x\T\B\x=\bm a\T\B^{-1}\bm a \), attained at \( \x=\B^{-1}\bm a \).
:::

:::

::: {.proof}
Substitute \( \x=\B^{-1/2}\bz \). The quotient in (a) becomes
\( \bz\T\B^{-1/2}\A\B^{-1/2}\bz/\bz\T\bz \), to which @thm-mat-extremal-rayleigh applies, and
\( \B^{-1}\A=\B^{-1/2}(\B^{-1/2}\A\B^{-1/2})\B^{1/2} \) is similar to it
(@prp-mat-eigen-basic(c)). For (b), if \( \B^{-1/2}\A\B^{-1/2}=\Q\bLambda\Q\T \), take
\( \bm P=\B^{-1/2}\Q \). For (c), the matrix \( \A=\bm a\bm a\T \) has
\( \B^{-1/2}\bm a\bm a\T\B^{-1/2} \) of rank at most one, whose only nonzero eigenvalue is its trace
\( \bm a\T\B^{-1}\bm a \). Substituting \( \x=\B^{-1}\bm a \) attains it.
:::

Part (c) is the Cauchy–Schwarz inequality in the inner product \( \x\T\B\bv \). It is the
algebra behind simultaneous confidence intervals (Chapter 12), and
part (b) is used for linear discriminants and for comparing two covariance matrices.

## Exercises

### A. Check your understanding

::: {#exr-mat-2x2-definite}
[A1]

Show that \( \begin{psmallmatrix}a&b\\b&c\end{psmallmatrix} \) is positive definite iff \( a>0 \) and
\( ac>b^2 \), and nonnegative definite iff \( a\ge0 \), \( c\ge0 \) and \( ac\ge b^2 \). Why is \( a\ge0 \) and
\( ac\ge b^2 \) not enough?
:::

### B. Practice

::: {#exr-mat-loewner-det}
[B1]

Let \( \A\succeq\B\succ\bzero \). Show that \( \det\A\ge\det\B \) and \( \tr\A\ge\tr\B \).
:::

::: {#exr-mat-adding-regressors}
[B2]

Let \( \X=[\X_1,\X_2] \) have full column rank. Show that the leading block of \( (\X\T\X)^{-1} \)
satisfies \( \bigl[(\X\T\X)^{-1}\bigr]_{11}\succeq(\X_1\T\X_1)^{-1} \), with equality iff
\( \X_1\T\X_2=\bzero \). Interpret this for the variances of least squares coefficients when regressors
are added to a model.
:::

::: {.solution}
Write \( \A=\X\T\X \) in blocks. By
@thm-mat-partitioned-inverse, \( [\A^{-1}]_{11}=\A_{11\cdot2}^{-1} \). Now
\( \A_{11}-\A_{11\cdot2}=\A_{12}\A_{22}^{-1}\A_{21} \) is nonnegative definite by
@prp-mat-pd-properties(c,e), and \( \A_{11\cdot2} \) is positive definite. So
\( \A_{11}\succeq\A_{11\cdot2}\succ\bzero \), and @prp-mat-pd-properties(f) gives
\( \A_{11\cdot2}^{-1}\succeq\A_{11}^{-1} \). Equality forces \( \A_{12}\A_{22}^{-1}\A_{21}=\bzero \), hence
\( \A_{21}\x=\bzero \) for all \( \x \) by @prp-mat-pd-properties(b), that is, \( \X_2\T\X_1=\bzero \).
When \( \Cov(\Y)=\sigma^2\I \), the covariance matrix of the coefficients of \( \X_1 \) is
\( \sigma^2[\A^{-1}]_{11} \) in the larger model and \( \sigma^2\A_{11}^{-1} \) in the smaller one. Adding
regressors never decreases the variance of any linear combination of the original coefficients,
and leaves it unchanged only when the new regressors are orthogonal to the old ones.
:::

### C. Going deeper

::: {#exr-mat-courant-fischer}
[C1]

Let \( \A \) be symmetric with eigenvalues \( \lambda_1\ge\dots\ge\lambda_n \). Prove the Courant–Fischer
theorem
\[
\lambda_k=\max_{\dim\mathcal S=k}\ \min_{\bzero\ne\x\in\mathcal S}\frac{\x\T\A\x}{\x\T\x}.
\]
Deduce that if \( \B \) is obtained from \( \A \) by deleting one row and the corresponding column, its
eigenvalues \( \mu_1\ge\dots\ge\mu_{n-1} \) interlace: \( \lambda_k\ge\mu_k\ge\lambda_{k+1} \).
:::

::: {#exr-mat-normal-mle}
[C2]

Let \( \bS \) be positive definite of size \( p\times p \). Show that
\( g(\bSigma)=-\log\det\bSigma-\tr(\bSigma^{-1}\bS) \), over positive definite \( \bSigma \), is maximized
uniquely at \( \bSigma=\bS \). *Hint:* write \( g \) in terms of the eigenvalues of
\( \bS^{-1/2}\bSigma\bS^{-1/2} \).
:::

::: {.solution}
Let \( \bm M=\bS^{-1/2}\bSigma\bS^{-1/2} \), which is positive
definite with eigenvalues \( \mu_1,\dots,\mu_p>0 \). Then \( \det\bSigma=\det\bS\det\bm M \) and
\( \tr(\bSigma^{-1}\bS)=\tr(\bS^{1/2}\bSigma^{-1}\bS^{1/2})=\tr(\bm M^{-1}) \). So
\[
g(\bSigma)=-\log\det\bS-\sum_{i=1}^p\Bigl(\log\mu_i+\frac1{\mu_i}\Bigr).
\]
The function \( \mu\mapsto\log\mu+1/\mu \) has derivative \( (\mu-1)/\mu^2 \) and a unique minimum at
\( \mu=1 \). So \( g \) is maximized iff every \( \mu_i=1 \), that is, \( \bm M=\I \) and \( \bSigma=\bS \). The maximum
is \( -\log\det\bS-p \).
:::
