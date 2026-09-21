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

### C. Going deeper

::: {#exr-mat-projection-order}
[C1]

Let \( \A \) and \( \B \) be symmetric idempotent of the same order. Show that the following are
equivalent: (i) \( \A-\B \) is nonnegative definite; (ii) \( \A-\B \) is idempotent; (iii)
\( \A\B=\B \); (iv) \( \B\A=\B \); (v) \( \C(\B)\subseteq\C(\A) \). Deduce that a difference of two
sums of squares is again a sum of squares exactly when one model is nested in the other, the
situation of [Chapter 11](../ch11-general-linear-hypothesis/index.html).
:::

::: {.solution}
First note that \( \A\x=\x \) for every \( \x\in\C(\A) \): if \( \x=\A\mathbf{y} \), then
\( \A\x=\A^2\mathbf{y}=\A\mathbf{y}=\x \).

(v)\( \Rightarrow \)(iii): each column of \( \B \) lies in \( \C(\B)\subseteq\C(\A) \) and is therefore
fixed by \( \A \), so \( \A\B=\B \). (iii)\( \Rightarrow \)(v): \( \C(\B)=\C(\A\B)\subseteq\C(\A) \).
(iii)\( \Leftrightarrow \)(iv): transpose, using the symmetry of \( \A \) and \( \B \).
(iii)\( \Rightarrow \)(ii): with \( \A\B=\B\A=\B \),
\( (\A-\B)^2=\A-\A\B-\B\A+\B=\A-\B \). (ii)\( \Rightarrow \)(i): \( \A-\B \) is symmetric and
idempotent, hence nonnegative definite by @thm-mat-idempotent(c).

(i)\( \Rightarrow \)(iii): let \( \x\in\C(\B) \), so \( \B\x=\x \) and \( \x\T\B\x=\norm{\x}^2 \). Then
\( \x\T\A\x\ge\x\T\B\x=\norm{\x}^2 \), while \( \x\T\A\x\le\norm{\x}^2 \) by @thm-mat-idempotent(c).
So \( \x\T(\I-\A)\x=0 \), and since \( \I-\A \) is nonnegative definite,
@prp-mat-pd-properties(b) gives \( (\I-\A)\x=\bzero \), that is, \( \A\x=\x \). Applying this to
the columns of \( \B \) gives \( \A\B=\B \).

For the last statement, if \( \A \) and \( \B \) are the projections attached to two models, then
\( \y\T\A\y-\y\T\B\y \) is again a quadratic form in a symmetric idempotent matrix — a sum of
squares with its own degrees of freedom \( \rank(\A)-\rank(\B) \) — exactly when
\( \C(\B)\subseteq\C(\A) \).
:::

::: {#exr-mat-commuting-projections}
[C2]

Let \( \A \) and \( \B \) be symmetric idempotent with \( \A\B=\B\A \). Show that \( \A\B \) is
symmetric idempotent with \( \C(\A\B)=\C(\A)\cap\C(\B) \), and that \( \A+\B-\A\B \) is symmetric
idempotent with column space \( \C(\A)+\C(\B) \); deduce
\( \rank(\A+\B-\A\B)=\rank(\A)+\rank(\B)-\rank(\A\B) \). Give two symmetric idempotent matrices of
order \( 2 \) that do not commute, and show that for them \( \A\B \) is not a projection onto the
intersection.
:::

::: {.solution}
Symmetry: \( (\A\B)\T=\B\T\A\T=\B\A=\A\B \). Idempotence:
\( (\A\B)^2=\A(\B\A)\B=\A^2\B^2=\A\B \). If \( \x\in\C(\A\B) \) then \( \x=\A(\B\mathbf{y})\in\C(\A) \) and,
using \( \A\B=\B\A \), also \( \x=\B(\A\mathbf{y})\in\C(\B) \). Conversely, if \( \x \) lies in both column
spaces, then \( \A\x=\x \) and \( \B\x=\x \), so \( \A\B\x=\x \) and \( \x\in\C(\A\B) \).

Put \( \mathbf{N}=\A+\B-\A\B \). It is symmetric, and expanding with \( \A\B=\B\A \) and
\( \A^2=\A \), \( \B^2=\B \) gives
\( \mathbf{N}^2=\A+\A\B-\A\B+\A\B+\B-\A\B-\A\B-\A\B+\A\B=\A+\B-\A\B \). Its column space lies in
\( \C(\A)+\C(\B) \); conversely \( \mathbf{N}\x=\x \) for \( \x\in\C(\A) \) (then
\( \mathbf{N}\x=\x+\B\x-\B\x \)) and for \( \x\in\C(\B) \), so \( \C(\A)+\C(\B)\subseteq\C(\mathbf{N}) \). The
rank statement is @prp-mat-idempotent-basic(b) applied to the three idempotent matrices,
since the trace is linear.

For a counterexample take \( \A=\mathbf{e}_1\mathbf{e}_1\T \) and \( \B=\frac12\mathbf{J}_2 \), the projections onto two
different lines. Then \( \A\B=\frac12\begin{psmallmatrix}1&1\\0&0\end{psmallmatrix} \) and
\( \B\A=\frac12\begin{psmallmatrix}1&0\\1&0\end{psmallmatrix} \) differ, \( \A\B \) is neither
symmetric nor idempotent, and \( \C(\A)\cap\C(\B)=\{\bzero\} \) although \( \A\B\ne\bzero \).
:::

::: {#exr-mat-idempotent-similar}
[C3]

Let \( \A \) be idempotent of order \( n \) and rank \( r \), not necessarily symmetric. Show that
\( \mathbf{P}^{-1}\A\mathbf{P}=\begin{psmallmatrix}\I_r&\bzero\\\bzero&\bzero\end{psmallmatrix} \) for some
nonsingular \( \mathbf{P} \), and conversely. Read off from this a second proof that \( \tr(\A)=r \), and
find the characteristic polynomial of \( \A \) and the dimensions of its eigenspaces.
:::

::: {.solution}
By @prp-mat-idempotent-basic(c), \( \Real^n=\C(\A)\dirsum\Null(\A) \), with the two
dimensions \( r \) and \( n-r \). Take a basis \( \mathbf{p}_1,\dots,\mathbf{p}_r \) of \( \C(\A) \) and a basis
\( \mathbf{p}_{r+1},\dots,\mathbf{p}_n \) of \( \Null(\A) \); together they are a basis of \( \Real^n \), so
\( \mathbf{P}=[\mathbf{p}_1,\dots,\mathbf{p}_n] \) is nonsingular. Every \( \x\in\C(\A) \) satisfies \( \A\x=\x \), as
in @exr-mat-projection-order, so \( \A\mathbf{p}_i=\mathbf{p}_i \) for \( i\le r \) and \( \A\mathbf{p}_j=\bzero \)
for \( j>r \). Column by column this says
\( \A\mathbf{P}=\mathbf{P}\begin{psmallmatrix}\I_r&\bzero\\\bzero&\bzero\end{psmallmatrix} \). Conversely, a
matrix similar to that block matrix is idempotent by @prp-mat-idempotent-basic(e), and has
the same rank.

Since the trace is unchanged by similarity (@eq-mat-trace-tricks), \( \tr(\A)=r \), which proves
@prp-mat-idempotent-basic(b) again without a rank factorization. Similar matrices have the
same characteristic polynomial (@prp-mat-eigen-basic(c)), so
\( \det(\lambda\I-\A)=\lambda^{n-r}(\lambda-1)^r \). The eigenspaces are \( \C(\A) \) for the
eigenvalue \( 1 \) and \( \Null(\A) \) for the eigenvalue \( 0 \), of dimensions \( r \) and \( n-r \): an
idempotent matrix is diagonalizable even when it is not symmetric.
:::
