# Inverses and generalized inverses

## The inverse of a square matrix

An \( n\times n \) matrix \( \A \) is **nonsingular** if there is a matrix \( \A^{-1} \) with
\( \A\A^{-1}=\A^{-1}\A=\I_n \). The inverse is unique: if \( \B \) and \( \bm C \) both qualify, then
\( \B=\B(\A\bm C)=(\B\A)\bm C=\bm C \).

::: {#prp-mat-inverse}
For an \( n\times n \) matrix \( \A \) the following are equivalent: (a) \( \A \) is nonsingular;
(b) \( \rank(\A)=n \); (c) \( \Null(\A)=\{\bzero\} \); (d) \( \B\A=\I \) for some \( \B \);
(e) \( \A\bm C=\I \) for some \( \bm C \). When they hold, \( \B=\bm C=\A^{-1} \), and
\[
(\A\T)^{-1}=(\A^{-1})\T,\qquad (\A\B)^{-1}=\B^{-1}\A^{-1}
\]
for nonsingular \( \A,\B \) of the same size. A nonsingular symmetric matrix has a symmetric inverse.
:::

::: {.proof}
(b)\( \Leftrightarrow \)(c) is rank–nullity. (a)\( \Rightarrow \)(d) and (a)\( \Rightarrow \)(e) are
trivial. (d)\( \Rightarrow \)(c): \( \A\x=\bzero \) gives \( \x=\B\A\x=\bzero \). (e)\( \Rightarrow \)(b):
\( \I=\A\bm C \) has rank \( n \), and \( \rank(\A\bm C)\le\rank(\A) \). (b)\( \Rightarrow \)(a): since
\( \C(\A)=\Real^n \), we can solve \( \A\bm c_j=\bm e_j \) for each \( j \), giving \( \A\bm C=\I \). Then
\( \A(\bm C\A-\I)=\bzero \), and because \( \Null(\A)=\{\bzero\} \) each column of \( \bm C\A-\I \) is
zero. So \( \bm C\A=\I \) as well. If \( \B\A=\I \) then \( \B=\B\A\A^{-1}=\A^{-1} \), and similarly
for \( \bm C \). The formulas are verified by multiplication:
\( \A\T(\A^{-1})\T=(\A^{-1}\A)\T=\I \) and \( (\A\B)(\B^{-1}\A^{-1})=\I \). If \( \A=\A\T \), then
\( \A^{-1}=(\A\T)^{-1}=(\A^{-1})\T \).
:::

A rectangular matrix has no two-sided inverse. If \( \X \) is \( n\times p \) with full column
rank, then \( \X\T\X \) is nonsingular by @prp-mat-rank-product(c), and
\( (\X\T\X)^{-1}\X\T \) is a *left inverse*: it multiplies \( \X \) on the left to give \( \I_p \).
This is the matrix that computes least squares coefficients.

## Generalized inverses

When \( \A \) is singular or rectangular we still want a matrix that inverts \( \A \) on the
vectors where inversion makes sense, namely those in \( \C(\A) \). The following weak
requirement turns out to be exactly what solving equations needs.

::: {#def-mat-ginverse}
[Generalized inverse]

A **generalized inverse** of an \( m\times n \) matrix \( \A \) is any \( n\times m \) matrix
\( \A\ginv \) with
\[
\A\A\ginv\A=\A .
\]
:::

Read columnwise, the definition says that \( \A\A\ginv \) acts as the identity on every column
of \( \A \), hence on all of \( \C(\A) \): if \( \bu=\A\bb \) then \( \A(\A\ginv\bu)=\A\A\ginv\A\bb=\bu \). So
\( \A\ginv\bu \) is a solution of \( \A\x=\bu \) whenever one exists. If \( \A \) is nonsingular,
multiplying the definition by \( \A^{-1} \) on both sides gives \( \A\ginv=\A^{-1} \). In
general there are many generalized inverses. [Chapter 6](../ch06-projections/index.html) uses them to write
projections, and [Chapter 8](../ch08-estimability/index.html) uses them to describe estimable functions in
models of less than full rank.

::: {#thm-mat-ginverse-exists}
[Existence]

Every matrix has a generalized inverse. Explicitly, if \( \rank(\A)=r\ge1 \) and \( \A_{IJ} \) is
the nonsingular \( r\times r \) submatrix of @prp-mat-nonsingular-submatrix, the
\( n\times m \) matrix \( \G \) that holds \( \A_{IJ}^{-1} \) in the rows indexed by \( J \) and the
columns indexed by \( I \), and zeros elsewhere, is a generalized inverse of \( \A \).
:::

::: {.proof}
For \( \A=\bzero \) any matrix of the right size works. Otherwise let \( \bm E_I \) and \( \bm E_J \)
be the columns of \( \I_m \) and \( \I_n \) indexed by \( I \) and \( J \). Then \( \G=\bm E_J\A_{IJ}^{-1}\bm
E_I\T \), \( \A\bm E_J=\A_{\cdot J} \) and \( \bm E_I\T\A=\A_{I\cdot} \). So
\( \A\G\A=\A_{\cdot J}\A_{IJ}^{-1}\A_{I\cdot}=\A \) by @eq-mat-skeleton.
:::

::: {#prp-mat-ginverse-props}
[Properties of generalized inverses]

Let \( \A \) be \( m\times n \) of rank \( r \) and let \( \G \) be a generalized inverse of \( \A \).

::: {.enumerate options="label=(\alph*)"}
1. \( \A\G \) and \( \G\A \) are idempotent. Moreover \( \C(\A\G)=\C(\A) \) and
           \( \Null(\G\A)=\Null(\A) \), and \( \rank(\A\G)=\rank(\G\A)=r\le\rank(\G) \).

2. \( \G\T \) is a generalized inverse of \( \A\T \).

3. If \( \G_1,\G_2 \) are generalized inverses of \( \A \), so is \( \G_1\A\G_2 \). In particular
           \( \G\A\G \) is a generalized inverse which also satisfies \( (\G\A\G)\A(\G\A\G)=\G\A\G \).

4. If \( \A \) is symmetric, then \( \G\T \) and \( \frac12(\G+\G\T) \) are generalized inverses
           of \( \A \), so a symmetric generalized inverse always exists.

5. \( \G+(\I-\G\A)\bU+\V(\I-\A\G) \) is a generalized inverse of \( \A \) for every \( \bU \)
           and \( \V \) of size \( n\times m \).
:::

:::

::: {.proof}
(a) \( (\A\G)^2=(\A\G\A)\G=\A\G \), and likewise for \( \G\A \). Next,
\( \C(\A)=\C(\A\G\A)\subseteq\C(\A\G)\subseteq\C(\A) \), and
\( \Null(\A)\subseteq\Null(\G\A)\subseteq\Null(\A\G\A)=\Null(\A) \). The ranks follow from
these and rank–nullity, and \( r=\rank(\A\G)\le\rank(\G) \).
(b) Transpose \( \A\G\A=\A \). (c) \( \A\G_1\A\G_2\A=\A\G_2\A=\A \). The second statement is
\( \G\A\G\A\G\A\G=\G(\A\G\A)\G\A\G=\G\A\G\A\G=\G\A\G \).
(d) By (b), \( \G\T \) is a generalized inverse of \( \A\T=\A \). The condition \( \A\G\A=\A \) is linear
in \( \G \), so it is preserved by averaging.
(e) Multiply out, using \( \A(\I-\G\A)=\bzero \) and \( (\I-\A\G)\A=\bzero \).
:::

A generalized inverse with the extra property in (c), \( \G\A\G=\G \), is called
**reflexive**. @exr-mat-reflexive shows that the reflexive ones are exactly those of rank
\( r \). Part (e) in fact produces *every* generalized inverse from a single one as \( \bU \)
and \( \V \) vary. Idempotence in (a) means that \( \A\G \) and \( \G\A \) are projections, generally
not orthogonal ones. [Chapter 6](../ch06-projections/index.html) shows that the particular combination
\( \X(\X\T\X)\ginv\X\T \) is symmetric, does not depend on the choice of generalized inverse,
and is the orthogonal projection onto \( \C(\X) \).

## The Moore–Penrose inverse

One generalized inverse is singled out by symmetry requirements.

::: {#def-mat-moore-penrose}
[Moore–Penrose inverse]

The **Moore–Penrose inverse** of \( \A \) is the matrix \( \A^+ \) satisfying
\[
\A\A^+\A=\A,\qquad \A^+\A\A^+=\A^+,\qquad (\A\A^+)\T=\A\A^+,\qquad (\A^+\A)\T=\A^+\A .
\]
:::

There is at most one such matrix. If \( \G \) and \( \bm H \) both satisfy the four conditions, then
\[
\G=\G\A\G=\G(\A\G)\T=\G\G\T(\A\bm H\A)\T=\G(\A\G)\T(\A\bm H)\T=\G\A\G\A\bm H=\G\A\bm H ,
\]
where the steps use the second and third conditions for \( \G \), the first for \( \bm H \), and the
third for both. The same argument with the fourth condition in place of the third gives
\[
\bm H=(\bm H\A)\T\bm H=\A\T\bm H\T\bm H=(\A\G\A)\T\bm H\T\bm H=(\G\A)\T(\bm H\A)\T\bm H
=\G\A\bm H\A\bm H=\G\A\bm H,
\]
so \( \G=\bm H \). Existence follows from the
singular value decomposition, which gives \( \A^+ \) explicitly (@thm-mat-svd).
