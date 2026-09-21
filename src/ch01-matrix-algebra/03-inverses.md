# Inverses and generalized inverses

## The inverse of a square matrix

An \( n\times n \) matrix \( \A \) is **nonsingular** if there is a matrix \( \A^{-1} \) with
\( \A\A^{-1}=\A^{-1}\A=\I_n \). The inverse is unique: if \( \B \) and \( \mathbf{C} \) both qualify, then
\( \B=\B(\A\mathbf{C})=(\B\A)\mathbf{C}=\mathbf{C} \).

::: {#prp-mat-inverse}
For an \( n\times n \) matrix \( \A \) the following are equivalent: (a) \( \A \) is nonsingular;
(b) \( \rank(\A)=n \); (c) \( \Null(\A)=\{\bzero\} \); (d) \( \B\A=\I \) for some \( \B \);
(e) \( \A\mathbf{C}=\I \) for some \( \mathbf{C} \). When they hold, \( \B=\mathbf{C}=\A^{-1} \), and
\[
(\A\T)^{-1}=(\A^{-1})\T,\qquad (\A\B)^{-1}=\B^{-1}\A^{-1}
\]
for nonsingular \( \A,\B \) of the same size. A nonsingular symmetric matrix has a symmetric inverse.
:::

::: {.proof}
(b)\( \Leftrightarrow \)(c) is rank–nullity. (a)\( \Rightarrow \)(d) and (a)\( \Rightarrow \)(e) are
trivial. (d)\( \Rightarrow \)(c): \( \A\x=\bzero \) gives \( \x=\B\A\x=\bzero \). (e)\( \Rightarrow \)(b):
\( \I=\A\mathbf{C} \) has rank \( n \), and \( \rank(\A\mathbf{C})\le\rank(\A) \). (b)\( \Rightarrow \)(a): since
\( \C(\A)=\Real^n \), we can solve \( \A\mathbf{c}_j=\mathbf{e}_j \) for each \( j \), giving \( \A\mathbf{C}=\I \). Then
\( \A(\mathbf{C}\A-\I)=\bzero \), and because \( \Null(\A)=\{\bzero\} \) each column of \( \mathbf{C}\A-\I \) is
zero. So \( \mathbf{C}\A=\I \) as well. If \( \B\A=\I \) then \( \B=\B\A\A^{-1}=\A^{-1} \), and similarly
for \( \mathbf{C} \). The formulas are verified by multiplication:
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
For \( \A=\bzero \) any matrix of the right size works. Otherwise let \( \mathbf{E}_I \) and \( \mathbf{E}_J \)
be the columns of \( \I_m \) and \( \I_n \) indexed by \( I \) and \( J \). Then \( \G=\mathbf{E}_J\A_{IJ}^{-1}\mathbf{E}_I\T \), \( \A\mathbf{E}_J=\A_{\cdot J} \) and \( \mathbf{E}_I\T\A=\A_{I\cdot} \). So
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

There is at most one such matrix. If \( \G \) and \( \mathbf{H} \) both satisfy the four conditions, then
\[
\G=\G\A\G=\G(\A\G)\T=\G\G\T(\A\mathbf{H}\A)\T=\G(\A\G)\T(\A\mathbf{H})\T=\G\A\G\A\mathbf{H}=\G\A\mathbf{H} ,
\]
where the steps use the second and third conditions for \( \G \), the first for \( \mathbf{H} \), and the
third for both. The same argument with the fourth condition in place of the third gives
\[
\mathbf{H}=(\mathbf{H}\A)\T\mathbf{H}=\A\T\mathbf{H}\T\mathbf{H}=(\A\G\A)\T\mathbf{H}\T\mathbf{H}=(\G\A)\T(\mathbf{H}\A)\T\mathbf{H}
=\G\A\mathbf{H}\A\mathbf{H}=\G\A\mathbf{H},
\]
so \( \G=\mathbf{H} \). Existence follows from the
singular value decomposition, which gives \( \A^+ \) explicitly (@thm-mat-svd).

## Exercises

### A. Check your understanding

::: {#exr-mat-ginverse-compute}
[A1]

Let
\( \A=\begin{psmallmatrix}1&2&3\\2&4&6\end{psmallmatrix} \).
(a) Check that \( \rank(\A)=1 \) and use the recipe of @thm-mat-ginverse-exists with
\( I=\{1\} \) and \( J=\{1\} \) to write down a generalized inverse \( \G \); verify \( \A\G\A=\A \).
(b) Show that \( \A\G \) is not symmetric, so \( \G\ne\A^+ \), although \( \G\A\G=\G \).
(c) Verify that
\( \A^+=\frac{1}{70}\begin{psmallmatrix}1&2\\2&4\\3&6\end{psmallmatrix} \)
satisfies the four conditions of @def-mat-moore-penrose.
:::

::: {.solution}
(a) The second row is twice the first, so the rank is one, and the nonsingular \( 1\times1 \)
submatrix at \( I=J=\{1\} \) is \( (1) \). The recipe puts \( 1 \) in the entry of the \( 3\times2 \)
matrix \( \G \) indexed by row \( 1 \) and column \( 1 \), and zeros elsewhere. Then
\( \A\G=\begin{psmallmatrix}1&0\\2&0\end{psmallmatrix} \) and \( \A\G\A=\A \).
(b) The matrix \( \A\G \) is visibly not symmetric. And \( \G(\A\G)=\G \), because \( \G \) has only
its \( (1,1) \) entry nonzero, so \( \G(\A\G) \) has the first row \( (1,0) \) of \( \A\G \) as its own
first row and zeros elsewhere.
(c) Write \( \A=\bu\bv\T \) with \( \bu=(1,2)\T \) and \( \bv=(1,2,3)\T \), so
\( \bu\T\bu=5 \), \( \bv\T\bv=14 \) and \( \A^+=\bv\bu\T/70 \). Then \( \A\A^+=\bu\bu\T/5 \) and
\( \A^+\A=\bv\bv\T/14 \), both symmetric; \( \A\A^+\A=\bu(\bu\T\bu)(\bv\T\bv)\bv\T/70=\bu\bv\T=\A \);
and \( \A^+\A\A^+=\bv(\bv\T\bv)(\bu\T\bu)\bu\T/70^2=\A^+ \).
:::

::: {#exr-mat-canonical-ginverses}
[A2]

Let \( \A=\begin{psmallmatrix}\I_r&\bzero\\\bzero&\bzero\end{psmallmatrix} \) be \( m\times n \).
Show that the generalized inverses of \( \A \) are exactly the \( n\times m \) matrices
\( \G=\begin{psmallmatrix}\I_r&\bU\\\V&\W\end{psmallmatrix} \) with \( \bU,\V,\W \) arbitrary; that
\( \G \) is reflexive iff \( \W=\V\bU \); and that \( \A^+=\A\T \). How many generalized inverses
does this \( \A \) have when \( r<\min(m,n) \), and how many when \( r=m=n \)?
:::

::: {.solution}
Partition \( \G \) conformably. Block multiplication gives
\( \A\G\A=\begin{psmallmatrix}\G_{11}&\bzero\\\bzero&\bzero\end{psmallmatrix} \), which equals
\( \A \) iff \( \G_{11}=\I_r \); the other blocks are unconstrained. With \( \G_{11}=\I_r \),
\( \G\A=\begin{psmallmatrix}\I_r&\bzero\\\V&\bzero\end{psmallmatrix} \) and
\( \G\A\G=\begin{psmallmatrix}\I_r&\bU\\\V&\V\bU\end{psmallmatrix} \), so \( \G\A\G=\G \) iff
\( \W=\V\bU \). For the Moore–Penrose inverse, \( \A\G \) and \( \G\A \) are symmetric only if
\( \bU=\bzero \) and \( \V=\bzero \), and then reflexivity forces \( \W=\bzero \). The remaining matrix
is \( \A\T \). When \( r<\min(m,n) \) at least one of \( \bU,\V,\W \) is a nonempty block, so there
are infinitely many generalized inverses; when \( r=m=n \) the three blocks are empty and
\( \G=\I_n=\A^{-1} \) is the only one, as it must be.
:::

::: {#exr-mat-jn-plus}
[A3]

Show that a symmetric idempotent matrix is its own Moore–Penrose inverse, and deduce
\( (\I_n-n^{-1}\mathbf{J}_n)^+=\I_n-n^{-1}\mathbf{J}_n \) for the centring matrix of @exm-mat-centring.
Show also that \( \mathbf{J}_n^+=n^{-2}\mathbf{J}_n \).
:::

### B. Practice

::: {#exr-mat-left-inverses}
[B1]

Let \( \X \) be \( n\times p \) with full column rank and write
\( \M=\X(\X\T\X)^{-1}\X\T \). Show that the left inverses of \( \X \), that is, the \( p\times n \)
matrices \( \bL \) with \( \bL\X=\I_p \), are exactly
\[
\bL=(\X\T\X)^{-1}\X\T+\mathbf{K}(\I_n-\M),\qquad \mathbf{K}\ \text{arbitrary},
\]
that every left inverse is a generalized inverse of \( \X \), and that \( (\X\T\X)^{-1}\X\T \) is the
only left inverse whose rows lie in \( \C(\X) \).
:::

::: {.solution}
Since \( \M\X=\X \), every matrix of the stated form satisfies
\( \bL\X=\I_p+\mathbf{K}(\X-\X)=\I_p \). Conversely, if \( \bL\X=\I_p \), then
\( \bL=\bL\M+\bL(\I-\M)=(\bL\X)(\X\T\X)^{-1}\X\T+\bL(\I-\M) \), which is of the stated form with
\( \mathbf{K}=\bL \). A left inverse satisfies \( \X\bL\X=\X \), so it is a generalized inverse. The rows
of \( (\X\T\X)^{-1}\X\T \) are the vectors \( \X(\X\T\X)^{-1}\mathbf{e}_j \), which lie in \( \C(\X) \); the rows
of \( \mathbf{K}(\I-\M) \) are of the form \( (\I-\M)\mathbf{c} \), and \( \X\T(\I-\M)=\bzero \) shows that such a
vector lies in \( \C(\X) \) only if it is zero, since \( \X\bb=(\I-\M)\mathbf{c} \) gives
\( \X\T\X\bb=\bzero \) and hence \( \X\bb=\bzero \) by @exr-mat-cancellation.
:::

::: {#exr-mat-xtx-ginverse}
[B2]

Let \( \X \) be \( n\times p \) of any rank and let \( \G \) be any generalized inverse of \( \X\T\X \).
Show that (a) \( \X\G\X\T\X=\X \), so \( \G\X\T \) is a generalized inverse of \( \X \);
(b) \( \X\T\X\G\X\T=\X\T \); and (c) \( \M=\X\G\X\T \) is symmetric and idempotent and is the same
matrix for every choice of \( \G \). *Hints:* for (a) use @exr-mat-cancellation; for (b) apply (a)
to \( \G\T \), which is also a generalized inverse of \( \X\T\X \) by @prp-mat-ginverse-props(b).
The matrix \( \M \) is the orthogonal projection onto \( \C(\X) \) of
[Chapter 6](../ch06-projections/index.html), and this exercise is what makes it well defined when
\( \X \) has less than full rank.
:::

::: {.solution}
(a) Because \( \G \) is a generalized inverse of \( \X\T\X \), we have
\( \X\T\X(\G\X\T\X)=\X\T\X\I_p \). Cancelling one factor \( \X\T \) on the left, which
@exr-mat-cancellation permits, gives \( \X\G\X\T\X=\X \). Hence
\( \X(\G\X\T)\X=\X \).
(b) Since \( \X\T\X \) is symmetric, \( \G\T \) is a generalized inverse of it, so (a) gives
\( \X\G\T\X\T\X=\X \); transposing yields \( \X\T\X\G\X\T=\X\T \).
(c) Write \( \M=\X\G\X\T \), so \( \M\T=\X\G\T\X\T \). Using (b),
\( \M\T\M=\X\G\T(\X\T\X\G\X\T)=\X\G\T\X\T=\M\T \). The left-hand side is symmetric, so \( \M\T \) is
symmetric, that is, \( \M\T=\M \); and then \( \M^2=\M\T\M=\M \). If \( \G_1,\G_2 \) are two
generalized inverses of \( \X\T\X \), then by (a) for \( \G_1 \) and (b) for \( \G_2 \),
\[
\X\G_2\X\T=(\X\G_1\X\T\X)\G_2\X\T=\X\G_1(\X\T\X\G_2\X\T)=\X\G_1\X\T .
\]
:::

::: {#exr-mat-all-ginverses}
[B3]

The text says that @prp-mat-ginverse-props(e) produces *every* generalized inverse from a
single one. Prove it: if \( \G_0 \) and \( \G \) are generalized inverses of \( \A \), find \( \bU \) and
\( \V \) with \( \G=\G_0+(\I-\G_0\A)\bU+\V(\I-\A\G_0) \).
:::

::: {.solution}
Put \( \mathbf{Z}=\G-\G_0 \), so that \( \A\mathbf{Z}\A=\A-\A=\bzero \). Take \( \bU=\mathbf{Z} \) and
\( \V=\G_0\A\mathbf{Z} \). Then
\[
(\I-\G_0\A)\mathbf{Z}+\G_0\A\mathbf{Z}(\I-\A\G_0)
=\mathbf{Z}-\G_0\A\mathbf{Z}+\G_0\A\mathbf{Z}-\G_0(\A\mathbf{Z}\A)\G_0=\mathbf{Z},
\]
so \( \G_0+(\I-\G_0\A)\bU+\V(\I-\A\G_0)=\G_0+\mathbf{Z}=\G \).
:::

### C. Going deeper

::: {#exr-mat-least-squares-ginverse}
[C1]

Call \( \G \) a **least-squares generalized inverse** of \( \A \) if it satisfies the first and third
conditions of @def-mat-moore-penrose, namely \( \A\G\A=\A \) and \( (\A\G)\T=\A\G \). Show that
(a) \( \G_0=(\A\T\A)\ginv\A\T \) is one, for every choice of generalized inverse of \( \A\T\A \);
(b) if \( \G \) is a least-squares generalized inverse, then \( \x=\G\bb \) minimizes
\( \norm{\A\x-\bb} \), for every \( \bb \); and (c) conversely, if \( \G\bb \) minimizes
\( \norm{\A\x-\bb} \) for every \( \bb \), then \( \G \) is a least-squares generalized inverse. So these
are exactly the matrices that solve least squares problems, whatever the rank of \( \A \).
:::

::: {.solution}
(a) is @exr-mat-xtx-ginverse with \( \X=\A \): part (a) there gives \( \A\G_0\A=\A \) and part (c)
gives the symmetry of \( \A\G_0 \).

(b) First, \( \A\T(\A\G-\I)=(\A\G)\T\A\T-\A\T=(\A\G\A)\T-\A\T=\bzero \), using the symmetry of
\( \A\G \) and then \( \A\G\A=\A \). Hence for every \( \x \) and \( \bb \) the cross term in
\[
\norm{\A\x-\bb}^2=\norm{\A(\x-\G\bb)}^2+2(\x-\G\bb)\T\A\T(\A\G-\I)\bb+\norm{\A\G\bb-\bb}^2
\]
vanishes, and \( \norm{\A\x-\bb}^2=\norm{\A(\x-\G\bb)}^2+\norm{\A\G\bb-\bb}^2\ge\norm{\A\G\bb-\bb}^2 \),
with equality at \( \x=\G\bb \).

(c) Let \( \G_0 \) be the matrix of (a). By (b) the minimum of \( \norm{\A\x-\bb} \) equals
\( \norm{\A\G_0\bb-\bb} \), so the hypothesis says \( \norm{\A\G\bb-\bb}=\norm{\A\G_0\bb-\bb} \) for
every \( \bb \). Applying the displayed decomposition with \( \G_0 \) in place of \( \G \) and
\( \x=\G\bb \),
\[
\norm{\A\G\bb-\bb}^2=\norm{\A(\G-\G_0)\bb}^2+\norm{\A\G_0\bb-\bb}^2 ,
\]
so \( \A(\G-\G_0)\bb=\bzero \) for every \( \bb \), that is, \( \A\G=\A\G_0 \). Then
\( \A\G\A=\A\G_0\A=\A \) and \( (\A\G)\T=(\A\G_0)\T=\A\G_0=\A\G \).
:::

::: {#exr-mat-reverse-order}
[C2]

Let \( \A \) be \( m\times k \) with full column rank and \( \B \) be \( k\times n \) with full row rank.
Check the four conditions of @def-mat-moore-penrose to show that
\( \A^+=(\A\T\A)^{-1}\A\T \) and \( \B^+=\B\T(\B\B\T)^{-1} \), and then that the *reverse-order law*
\( (\A\B)^+=\B^+\A^+ \) holds. (Without the rank conditions it can fail:
@exr-mat-mp-properties asks for an example.)
:::

::: {.solution}
Both \( \A\T\A \) and \( \B\B\T \) are nonsingular by @prp-mat-rank-product(c). Writing
\( \G=(\A\T\A)^{-1}\A\T \), we have \( \G\A=\I_k \), so \( \A\G\A=\A \) and \( \G\A\G=\G \), while
\( \A\G=\A(\A\T\A)^{-1}\A\T \) is symmetric and \( \G\A=\I_k \) is symmetric. Hence \( \G=\A^+ \). Transposing the four
conditions for \( \A^+ \) gives the four conditions for \( (\A^+)\T \) as a Moore–Penrose inverse of
\( \A\T \), so \( (\A\T)^+=(\A^+)\T \); applying the first formula to \( \B\T \), which has full
column rank, gives \( \B^+=\bigl((\B\T)^+\bigr)\T=\B\T(\B\B\T)^{-1} \).

Now put \( \mathbf{P}=\B^+\A^+=\B\T(\B\B\T)^{-1}(\A\T\A)^{-1}\A\T \). Then
\( \A\B\mathbf{P}=\A(\A\T\A)^{-1}\A\T \) and \( \mathbf{P}\A\B=\B\T(\B\B\T)^{-1}\B \), both symmetric, because
\( \B\B\T(\B\B\T)^{-1}=\I \) and \( (\A\T\A)^{-1}\A\T\A=\I \) collapse the middle factors. The same
collapsing gives \( \A\B\mathbf{P}\A\B=\A(\A\T\A)^{-1}\A\T\A\B=\A\B \) and
\( \mathbf{P}\A\B\mathbf{P}=\B\T(\B\B\T)^{-1}\B\B\T(\B\B\T)^{-1}(\A\T\A)^{-1}\A\T=\mathbf{P} \). All four
conditions hold, so \( \mathbf{P}=(\A\B)^+ \).
:::
