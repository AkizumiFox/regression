# Idempotent matrices

A square matrix \( \A \) is **idempotent** if \( \A^2=\A \). Applying it twice has the same effect
as applying it once, which is how a projection behaves. Examples are \( \I \), \( \bzero \), the
centring matrix of @exm-mat-centring, the products \( \A\G \) and \( \G\A \) for a generalized
inverse \( \G \) (@prp-mat-ginverse-props(a)), and the non-symmetric
\( \begin{psmallmatrix}1&1\\0&0\end{psmallmatrix} \). In linear models, sums of squares are
quadratic forms in *symmetric* idempotent matrices, and their degrees of freedom are ranks.

::: {#prp-mat-idempotent-basic}
[General idempotent matrices]

Let \( \A \) be \( n\times n \) and idempotent of rank \( r \).

::: {.enumerate options="label=(\alph*)"}
1. Every eigenvalue of \( \A \) is \( 0 \) or \( 1 \).

2. \( \rank(\A)=\tr(\A) \).

3. \( \I-\A \) is idempotent, \( \A(\I-\A)=(\I-\A)\A=\bzero \), \( \C(\I-\A)=\Null(\A) \), and
           \( \rank(\A)+\rank(\I-\A)=n \).

4. An idempotent matrix that is nonsingular equals \( \I \).

5. \( \mathbf{P}^{-1}\A\mathbf{P} \) is idempotent for nonsingular \( \mathbf{P} \), and \( \A\T \) is idempotent.
:::

:::

::: {.proof}
(a) If \( \A\x=\lambda\x \) with \( \x\ne\bzero \), then \( \lambda\x=\A\x=\A^2\x=\lambda^2\x \), so
\( \lambda^2=\lambda \).
(b) If \( r=0 \) both sides vanish. Otherwise take a rank factorization \( \A=\B\mathbf{F} \)
(@thm-mat-row-col-rank). Then \( \B\mathbf{F}\B\mathbf{F}=\B\mathbf{F} \). Since \( \B \) has full column rank
and \( \mathbf{F} \) has full row rank, \( \B \) has a left inverse and \( \mathbf{F} \) a right inverse
([Section 1.3](03-inverses.html)), and cancelling them gives \( \mathbf{F}\B=\I_r \). So
\( \tr(\A)=\tr(\B\mathbf{F})=\tr(\mathbf{F}\B)=r \) by @thm-mat-trace-cyclic.
(c) The first two claims are one-line expansions. \( \C(\I-\A)\subseteq\Null(\A) \) because
\( \A(\I-\A)=\bzero \), and \( \bv\in\Null(\A) \) gives \( \bv=(\I-\A)\bv \). The rank statement is then
rank–nullity.
(d) Multiply \( \A^2=\A \) by \( \A^{-1} \).
(e) \( (\mathbf{P}^{-1}\A\mathbf{P})^2=\mathbf{P}^{-1}\A^2\mathbf{P} \) and \( (\A\T)^2=(\A^2)\T \).
:::

By (c), every \( \x \) splits as \( \A\x+(\I-\A)\x \) with the first part in \( \C(\A) \) and the second in
\( \Null(\A) \). So \( \Real^n=\C(\A)\dirsum\Null(\A) \), and \( \A \) projects onto \( \C(\A) \) along
\( \Null(\A) \). The projection is orthogonal exactly when the two subspaces are orthogonal, which
happens iff \( \A \) is symmetric (@exr-mat-oblique).

::: {#thm-mat-idempotent}
[Symmetric idempotent matrices]

Let \( \A \) be an \( n\times n \) symmetric matrix of rank \( r \). Then \( \A \) is idempotent iff all its
eigenvalues are \( 0 \) or \( 1 \). In that case

::: {.enumerate options="label=(\alph*)"}
1. \( \A=\Q_1\Q_1\T \), where the \( n\times r \) matrix \( \Q_1 \) has orthonormal columns spanning
           \( \C(\A) \);

2. \( \rank(\A)=\tr(\A)=r \), the multiplicity of the eigenvalue \( 1 \);

3. \( \A \) and \( \I-\A \) are nonnegative definite, and \( \x\T\A\x=\norm{\A\x}^2\le\norm{\x}^2 \).
:::

:::

::: {.proof}
If \( \A \) is idempotent its eigenvalues are \( 0 \) or \( 1 \) by @prp-mat-idempotent-basic(a).
Conversely, if \( \A=\Q\bLambda\Q\T \) with each \( \lambda_i\in\{0,1\} \), then
\( \A^2=\Q\bLambda^2\Q\T=\Q\bLambda\Q\T=\A \). In that case, order the eigenvalues so that the
\( r \) ones come first and let \( \Q_1 \) hold the corresponding eigenvectors. Then
\( \A=\Q\bLambda\Q\T=\Q_1\Q_1\T \), and \( \C(\A)=\C(\Q_1) \) by
@cor-mat-spectral-consequences(a). This is (a). The trace is the sum of the
eigenvalues, which is (b). For (c), the eigenvalues of \( \A \) and of \( \I-\A \) are nonnegative, and
\( \x\T\A\x=\x\T\A\T\A\x \). Also \( \norm{\x}^2-\norm{\A\x}^2=\x\T(\I-\A)\x\ge0 \).
:::

[Chapter 6](../ch06-projections/index.html) shows that symmetric idempotent matrices are exactly the orthogonal
projections, and \( \Q_1\Q_1\T \) in (a) is the projection onto \( \C(\A) \). In
[Chapter 4](../ch04-quadratic-forms/index.html), part (b) becomes the degrees of freedom of a chi-squared
distribution.

::: {#thm-mat-idempotent-sum}
[Sums of symmetric idempotent matrices]

Let \( \A_1,\dots,\A_k \) be symmetric \( n\times n \) matrices and \( \A=\sum_i\A_i \). Consider

::: {.enumerate options="label=(\roman*)"}
1. \( \A \) is idempotent;

2. every \( \A_i \) is idempotent;

3. \( \A_i\A_j=\bzero \) for all \( i\ne j \).
:::

Any two of these imply the third. Moreover, if \( \A=\I_n \) and
\( \sum_i\rank(\A_i)=n \), then (ii) and (iii) hold.
:::

::: {.proof}
(ii) and (iii) imply (i): \( \A^2=\sum_i\A_i^2+\sum_{i\ne j}\A_i\A_j=\sum_i\A_i=\A \).

(i) and (ii) imply (iii): By @thm-mat-idempotent(c), \( \norm{\x}^2\ge\x\T\A\x \) for all
\( \x \). Fix \( j \) and take \( \x=\A_j\mathbf{y} \), so that \( \A_j\x=\x \) and \( \x\T\A_j\x=\norm{\x}^2 \). Then
\[
\norm{\x}^2\ge\x\T\A\x=\norm{\x}^2+\sum_{i\ne j}\x\T\A_i\x=\norm{\x}^2+\sum_{i\ne j}\norm{\A_i\x}^2 .
\]
So \( \A_i\A_j\mathbf{y}=\A_i\x=\bzero \) for every \( i\ne j \) and every \( \mathbf{y} \).

(i) and (iii) imply (ii): By (iii), \( \A\A_i=\A_i^2 \). Using \( \A^2=\A \),
\( \A_i^2=\A\A_i=\A^2\A_i=\A\A_i^2=\A_i^3 \). The eigenvalues of the symmetric \( \A_i \) therefore
satisfy \( \lambda^2=\lambda^3 \), so they lie in \( \{0,1\} \), and \( \A_i \) is idempotent by @thm-mat-idempotent.

For the last statement, let \( r_i=\rank(\A_i) \). For each \( i \),
\( \I-\A_i=\sum_{j\ne i}\A_j \) has rank at most \( \sum_{j\ne i}r_j=n-r_i \)
(@prp-mat-rank-product(d)). So the eigenvalue \( 1 \) of \( \A_i \), whose eigenspace is
\( \Null(\I-\A_i) \), has multiplicity at least \( r_i \). But \( \A_i \) has exactly \( r_i \) nonzero
eigenvalues, so all of them equal \( 1 \) and \( \A_i \) is idempotent. Now (i) holds trivially, and
(iii) follows from the second part of the proof.
:::

The last statement is the algebraic core of Cochran's theorem ([Chapter 4](../ch04-quadratic-forms/index.html)): if
\( \y\T\y=\sum_i\y\T\A_i\y \) and the ranks add up to \( n \), the quadratic forms are built from
mutually orthogonal projections.

## Exercises

### A. Check your understanding

::: {#exr-mat-idempotent-products}
[A1]

Let \( \A,\B \) be idempotent of the same size. Show that \( \A\B \) is idempotent if \( \A\B=\B\A \), and
that \( \A+\B \) is idempotent iff \( \A\B+\B\A=\bzero \). Show that the latter
forces \( \A\B=\bzero \) (symmetry is not needed).
:::

### B. Practice

::: {#exr-mat-oblique}
[B1]

Let \( \A \) be idempotent. Show that \( \A \) is symmetric iff \( \C(\A)\perp\Null(\A) \).
:::

::: {.solution}
Suppose \( \A \) is symmetric. For \( \bu=\A\x\in\C(\A) \) and
\( \bv\in\Null(\A) \), \( \bu\T\bv=\x\T\A\bv=0 \). Conversely, suppose \( \C(\A)\perp\Null(\A) \). For all
\( \x,\mathbf{y} \), \( (\I-\A)\mathbf{y}\in\Null(\A) \), so \( \x\T\A\T(\I-\A)\mathbf{y}=0 \). Hence \( \A\T(\I-\A)=\bzero \), that is,
\( \A\T=\A\T\A \). The right side is symmetric, so \( \A\T \) is symmetric, and so is \( \A \).
:::
