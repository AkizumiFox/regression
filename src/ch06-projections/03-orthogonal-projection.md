# Orthogonal projection

## Nearest points

::: {#thm-proj-projection-theorem}
[Projection theorem]

Let \( \mathcal S \) be a subspace of \( \Real^n \) and \( \y\in\Real^n \). There is exactly
one \( \hat{\y}\in\mathcal S \) such that
\[
\norm{\y-\hat{\y}}\le\norm{\y-\bu}\qquad\text{for all }\bu\in\mathcal S .
\]
It is the unique \( \hat{\y}\in\mathcal S \) with \( \y-\hat{\y}\perp\mathcal S \).
Moreover, for every \( \bu\in\mathcal S \),
\[
\norm{\y-\bu}^2=\norm{\y-\hat{\y}}^2+\norm{\hat{\y}-\bu}^2 .
\]{#eq-proj-distance-split}

:::

::: {.proof}
By @thm-proj-direct-sum there is a unique decomposition \( \y=\hat{\y}+\bv \)
with \( \hat{\y}\in\mathcal S \) and \( \bv\in\mathcal S\perpc \). For \( \bu\in\mathcal S \),
\( \y-\bu=\bv+(\hat{\y}-\bu) \), where \( \bv\perp\mathcal S \) and
\( \hat{\y}-\bu\in\mathcal S \). The Pythagorean identity @eq-proj-pythagoras
then gives @eq-proj-distance-split. The right-hand side is at least
\( \norm{\y-\hat{\y}}^2 \), with equality iff \( \bu=\hat{\y} \). So \( \hat{\y} \) is the unique
nearest point. If some other \( \bw\in\mathcal S \) also had \( \y-\bw\perp\mathcal S \),
then \( \y=\bw+(\y-\bw) \) would be a second decomposition, contradicting uniqueness.
:::

The theorem has two halves, and it helps to keep them apart. The
*variational* half says \( \hat{\y} \) minimizes a distance. The
*orthogonality* half says the error \( \y-\hat{\y} \) is perpendicular to the
subspace. The first half explains why least squares is a sensible estimation
principle. The second is what we compute with, because orthogonality is a set of
linear equations.

## The projection as a matrix

Write \( \bP_{\mathcal S}\y \) for the nearest point \( \hat{\y} \). The map
\( \y\mapsto\bP_{\mathcal S}\y \) is linear: if \( \y_k=\hat{\y}_k+\bv_k \) are the
decompositions of \( \y_1 \) and \( \y_2 \), then
\( a\y_1+b\y_2=(a\hat{\y}_1+b\hat{\y}_2)+(a\bv_1+b\bv_2) \) is a decomposition of the
same kind, and by uniqueness it is *the* decomposition. A linear map on
\( \Real^n \) is given by a unique \( n\times n \) matrix, since two matrices that agree
on every vector agree on the standard basis vectors and are therefore equal.

::: {#def-proj-projection}
[Orthogonal projection matrix]

The **orthogonal projection** onto a subspace \( \mathcal S \) is the unique
\( n\times n \) matrix \( \bP_{\mathcal S} \) such that, for every \( \y \),
\( \bP_{\mathcal S}\y\in\mathcal S \) and \( \y-\bP_{\mathcal S}\y\perp\mathcal S \).
When \( \mathcal S=\C(\X) \) for the model matrix of a linear model, we write
\( \M=\bP_{\C(\X)} \).
:::

The most useful formula for \( \bP_{\mathcal S} \) uses an orthonormal basis.

::: {#prp-proj-orthonormal-formula}
Let the columns of the \( n\times r \) matrix \( \Q \) be an orthonormal basis for
\( \mathcal S \), so that \( \Q\T\Q=\I_r \). Then
\[
\bP_{\mathcal S}=\Q\Q\T=\sum_{i=1}^r\bu_i\bu_i\T,
\]
where \( \bu_1,\dots,\bu_r \) are the columns of \( \Q \).
:::

::: {.proof}
For any \( \y \), \( \Q\Q\T\y=\Q(\Q\T\y)\in\C(\Q)=\mathcal S \), and
\( \Q\T(\y-\Q\Q\T\y)=\Q\T\y-\Q\T\y=\bzero \), so the difference is orthogonal to
every column of \( \Q \) and hence to \( \mathcal S \). By uniqueness,
\( \Q\Q\T=\bP_{\mathcal S} \).
:::

::: {#exm-proj-line}
[Projections onto a line]

For \( \mathcal S=\spn(\bv) \) with \( \bv\ne\bzero \), the single vector
\( \bv/\norm{\bv} \) is an orthonormal basis, so
\[
\bP_{\mathcal S}=\frac{\bv\bv\T}{\bv\T\bv},
\qquad
\bP_{\mathcal S}\y=\frac{\bv\T\y}{\bv\T\bv}\,\bv .
\]
With \( \bv=\bone \) this gives \( \bP=n^{-1}\bone\bone\T \), the matrix with every entry
\( 1/n \), and \( \bP\y=\bar y\,\bone \). *The sample mean is a projection*: fitting
the intercept-only model replaces every observation by \( \bar y \), and the
residuals \( y_i-\bar y \) form a vector in \( \bone\perpc \), the set of vectors whose
entries sum to zero.
:::

## Characterizing projection matrices

Which matrices are orthogonal projections? @prp-proj-orthonormal-formula
shows that every projection has the form \( \Q\Q\T \), so it is *symmetric*, and
\( \Q\Q\T\Q\Q\T=\Q\I\Q\T=\Q\Q\T \), so it is *idempotent*. The converse is the
basic recognition theorem.

::: {#thm-proj-sym-idem}
An \( n\times n \) matrix \( \bP \) is the orthogonal projection onto some subspace iff
\( \bP\T=\bP \) and \( \bP^2=\bP \). In that case the subspace is \( \C(\bP) \), and
\( \I-\bP \) is the orthogonal projection onto \( \C(\bP)\perpc \).
:::

::: {.proof}
Necessity was shown above. For sufficiency, let \( \bP \) be symmetric and
idempotent and let \( \y\in\Real^n \). Certainly \( \bP\y\in\C(\bP) \). For any
\( \bu=\bP\bb\in\C(\bP) \),
\[
\bu\T(\y-\bP\y)=\bb\T\bP\T(\I-\bP)\y=\bb\T(\bP-\bP^2)\y=0,
\]
so \( \y-\bP\y\perp\C(\bP) \) and \( \bP=\bP_{\C(\bP)} \). For the last claim, \( \I-\bP \)
is symmetric and idempotent, since \( (\I-\bP)^2=\I-2\bP+\bP^2=\I-\bP \). It sends
\( \y \) to \( \y-\bP\y\in\C(\bP)\perpc \), and \( \y-(\I-\bP)\y=\bP\y\perp\C(\bP)\perpc \).
:::

Before this theorem, checking that a matrix computes a projection meant
thinking about subspaces. Afterwards it is two matrix identities. Also, a
projection is determined by its range: if \( \bP_1,\bP_2 \) are symmetric idempotent
with \( \C(\bP_1)=\C(\bP_2) \), then \( \bP_1=\bP_2 \), because both are the projection
onto the same subspace.

::: {#prp-proj-trace-rank}
[Spectrum, rank and trace]

Let \( \bP \) be the orthogonal projection onto an \( r \)-dimensional subspace
\( \mathcal S \). Then:

::: {.enumerate options="label=(\alph*)"}
1. every eigenvalue of \( \bP \) is \( 0 \) or \( 1 \), and \( 1 \) has multiplicity \( r \);

2. \( \rank(\bP)=\tr(\bP)=r \);

3. \( \bP\bv=\bv \) iff \( \bv\in\mathcal S \), and \( \bP\bv=\bzero \) iff \( \bv\in\mathcal S\perpc \);

4. \( \bv\T\bP\bv=\norm{\bP\bv}^2\le\norm{\bv}^2 \), so \( \bP \) and \( \I-\bP \) are nonnegative definite.
:::

:::

::: {.proof}
Write \( \bP=\Q\Q\T \) as in @prp-proj-orthonormal-formula, and extend \( \Q \) to
an orthogonal matrix \( \vect{U}=[\Q,\Q_\perp] \). Then
\( \bP=\vect{U}\diag(\I_r,\bzero)\vect{U}\T \) is a spectral decomposition, which gives (a).
Rank and trace of a symmetric matrix are the number of nonzero eigenvalues and
their sum, which gives (b). Directly, \( \tr(\Q\Q\T)=\tr(\Q\T\Q)=\tr(\I_r)=r \).
Part (c) restates the projection theorem: the nearest point of \( \mathcal S \) to a
vector already in \( \mathcal S \) is the vector itself, and \( \bv\perp\mathcal S \) iff
\( \bzero \) is its nearest point. For (d), \( \bv\T\bP\bv=\bv\T\bP\T\bP\bv \), and
\( \norm{\bv}^2=\norm{\bP\bv}^2+\norm{(\I-\bP)\bv}^2 \) by Pythagoras.
:::

The identity \( \rank=\tr \) gets a lot of use. The dimension of a subspace, which is
the “degrees of freedom” of every sum of squares in this book, can be read off as a
sum of diagonal entries without finding a basis.

## Generalized inverses

@prp-proj-orthonormal-formula needs an orthonormal basis. In statistics the
subspace usually comes as the column space of a model matrix \( \X \) whose columns
are not orthonormal and may not even be independent. To write \( \M \) directly in
terms of \( \X \) we need a substitute for \( (\X\T\X)^{-1} \) that exists even when
\( \X\T\X \) is singular.

A **generalized inverse** of \( \A \) is any matrix \( \A\ginv \) with
\( \A\A\ginv\A=\A \) (@def-mat-ginverse), and every matrix has one: if \( \A \) has
rank \( r \), invert a nonsingular \( r\times r \) submatrix and pad with zeros
(@thm-mat-ginverse-exists). If \( \A \) is square and nonsingular the condition forces
\( \A\ginv=\A^{-1} \); when \( \A \) is singular there are infinitely many, and the
whole family is described in @exr-proj-all-ginverses.

Among the many generalized inverses, one is singled out by extra symmetry
conditions. The **Moore–Penrose inverse** \( \A^{+} \) is the unique matrix
satisfying \( \A\A^+\A=\A \), \( \A^+\A\A^+=\A^+ \), and both \( \A\A^+ \) and \( \A^+\A \)
symmetric (Penrose 1955). [Section 6.10](10-computation.html) computes
it from the singular value decomposition. For projections we will need
none of these extra conditions, and that is the point of the next result.

::: {#lem-proj-ginverse-invariance}
Let \( \X \) be \( n\times p \) and let \( \G \) be any generalized inverse of \( \X\T\X \).
Then:

::: {.enumerate options="label=(\alph*)"}
1. \( \X\G\X\T\X=\X \) and \( \X\T\X\G\X\T=\X\T \);

2. \( \X\G\X\T \) does not depend on the choice of generalized inverse;

3. \( \X\G\X\T \) is symmetric, even when \( \G \) is not.
:::

:::

::: {.proof}
By @cor-proj-gram-colspace(b), \( \C(\X\T)=\C(\X\T\X) \), so there is a
\( p\times n \) matrix \( \B \) with \( \X\T=\X\T\X\B \). Transposing, \( \X=\B\T\X\T\X \). Then
\[
\X\G\X\T\X=\B\T(\X\T\X)\G(\X\T\X)=\B\T\X\T\X=\X,
\]
which is the first identity in (a). The second follows the same way, using
\( \X\T=\X\T\X\B \) on the right. For (b) and (c),
\[
\X\G\X\T=\B\T(\X\T\X)\G(\X\T\X)\B=\B\T\X\T\X\B ,
\]
which involves no \( \G \) at all, and \( (\B\T\X\T\X\B)\T=\B\T\X\T\X\B \).
:::

::: {#thm-proj-M-formula}
For any \( n\times p \) matrix \( \X \) and any generalized inverse \( (\X\T\X)\ginv \),
\[
\M=\X(\X\T\X)\ginv\X\T
\]
is the orthogonal projection onto \( \C(\X) \). In particular \( \rank(\M)=\tr(\M)=\rank(\X) \).
If \( \X \) has full column rank, \( \M=\X(\X\T\X)^{-1}\X\T \).
:::

::: {.proof}
By @lem-proj-ginverse-invariance, \( \M \) is symmetric and
\( \M^2=\X(\X\T\X)\ginv\bigl[\X\T\X(\X\T\X)\ginv\X\T\bigr]=\X(\X\T\X)\ginv\X\T=\M \).
By @thm-proj-sym-idem, \( \M \) is the projection onto \( \C(\M) \). Clearly
\( \C(\M)\subseteq\C(\X) \), and \( \M\X=\X \) (again by the lemma) shows
\( \C(\X)\subseteq\C(\M) \). The rank statement is
@prp-proj-trace-rank(b).
:::

It is worth pausing on how strong @thm-proj-M-formula is. The matrix
\( (\X\T\X)\ginv \) can be almost anything: symmetric or not, reflexive or not, with
entries of any size. The sandwich \( \X(\cdot)\X\T \) removes all of that freedom. We
will see in [Section 6.4](04-least-squares.html) that this is exactly the difference
between coefficients and fitted values, made concrete.

::: {#exm-proj-oneway-M}
[The one-way layout]

In the setting of @exm-proj-oneway-rank, sort the observations by group and
let group \( k \) have \( n_k \) members. The model space is spanned by the orthogonal
vectors \( \bz_1,\dots,\bz_g \), so by @prp-proj-orthonormal-formula with
\( \bu_k=\bz_k/\sqrt{n_k} \),
\[
\M=\sum_{k=1}^g\frac{\bz_k\bz_k\T}{n_k}
=\begin{pmatrix}
    n_1^{-1}\bone_{n_1}\bone_{n_1}\T & & \\
    & \ddots & \\
    & & n_g^{-1}\bone_{n_g}\bone_{n_g}\T
\end{pmatrix},
\]
a block-diagonal matrix of within-group averaging operators. \( \M\y \) replaces each
observation by its group mean. The intercept column adds nothing to \( \C(\X) \), so
it does not appear, and \( \tr(\M)=\sum_k n_k/n_k=g=\rank(\X) \). We got here without
a generalized inverse. @thm-proj-M-formula guarantees that any generalized
inverse of the singular \( (g+1)\times(g+1) \) matrix \( \X\T\X \) produces this same
matrix. [Section 6.4](04-least-squares.html) checks this numerically.
:::

## Sums and products of projections

Complicated subspaces are often built from simpler ones, and their projections
are built the same way.

::: {#thm-proj-sum}
Let \( \bP_1,\bP_2 \) be orthogonal projections. Then \( \bP_1+\bP_2 \) is an
orthogonal projection iff \( \bP_1\bP_2=\bzero \), that is, iff
\( \C(\bP_1)\perp\C(\bP_2) \). In that case it projects onto
\( \C(\bP_1)\dirsum\C(\bP_2) \).
:::

::: {.proof}
\( \bP_1+\bP_2 \) is always symmetric, and
\( (\bP_1+\bP_2)^2=\bP_1+\bP_2+\bP_1\bP_2+\bP_2\bP_1 \). So it is idempotent iff
\( \bP_1\bP_2+\bP_2\bP_1=\bzero \). If \( \bP_1\bP_2=\bzero \), then transposing gives
\( \bP_2\bP_1=\bzero \) and the condition holds. Conversely, suppose
\( \bP_1\bP_2+\bP_2\bP_1=\bzero \). Multiplying on the left by \( \bP_1 \) gives
\( \bP_1\bP_2+\bP_1\bP_2\bP_1=\bzero \), and on the right by \( \bP_1 \) gives
\( \bP_1\bP_2\bP_1+\bP_2\bP_1=\bzero \). Subtracting, \( \bP_1\bP_2=\bP_2\bP_1 \), and
then \( 2\bP_1\bP_2=\bzero \). The equivalence with orthogonality of ranges holds
because \( \bP_1\bP_2=\bzero \) says exactly that every vector of \( \C(\bP_2) \) is sent
to \( \bzero \) by \( \bP_1 \), that is (@prp-proj-trace-rank(c)), lies in
\( \C(\bP_1)\perpc \). Finally, \( \C(\bP_1+\bP_2)\subseteq\C(\bP_1)+\C(\bP_2) \), and for
\( \bu=\bP_1\bu+\bP_2\bu \) in that sum we have \( (\bP_1+\bP_2)\bu=\bu \), so the ranges
agree.
:::

By induction, if \( \bP_1,\dots,\bP_k \) are projections with \( \bP_i\bP_j=\bzero \) for
\( i\ne j \), their sum projects onto the orthogonal direct sum of their ranges. When
the sum is \( \I \), we get an orthogonal decomposition
\[
\y=\bP_1\y+\dots+\bP_k\y,\qquad
\norm{\y}^2=\norm{\bP_1\y}^2+\dots+\norm{\bP_k\y}^2 .
\]
Every analysis of variance table in this book is a table of such a decomposition.

Products behave differently. \( \bP_1\bP_2 \) need not be symmetric or idempotent, and it is
generally not a projection. @exr-proj-product shows that it is one
iff \( \bP_1 \) and \( \bP_2 \) commute, in which case it projects onto
\( \C(\bP_1)\cap\C(\bP_2) \). The most important instance is the nested case
\( \C(\bP_2)\subseteq\C(\bP_1) \), treated in [Section 6.5](05-nested.html).

## Exercises

### A. Check your understanding

::: {#exr-proj-small}
[A1]

Let \( \mathcal S=\spn\{(1,1,0)\T,(0,1,1)\T\} \). Compute \( \bP_{\mathcal S} \) in two ways: by
Gram–Schmidt and @prp-proj-orthonormal-formula, and by @thm-proj-M-formula.
Verify symmetry, idempotence and \( \tr\bP_{\mathcal S}=2 \).
Find \( \mathcal S\perpc \) and check that \( \I-\bP_{\mathcal S} \) projects onto it.
:::

::: {#exr-proj-spectral-converse}
[A2]

Show that a symmetric matrix whose eigenvalues all lie in \( \{0,1\} \) is an orthogonal
projection. Give a non-symmetric \( 2\times2 \) matrix with eigenvalues \( 0 \) and \( 1 \) that is
not idempotent, and one that is idempotent.
:::

::: {#exr-proj-affine-projection}
[A3]

Let \( \bP \) be an orthogonal projection and \( a,b \) real numbers. Show that the eigenvalues of
\( a\bP+b\I \) are \( b \) and \( a+b \). Show that when \( b\ne0 \) and \( a+b\ne0 \),
\[
(a\bP+b\I)^{-1}=\frac1b\Bigl(\I-\frac{a}{a+b}\bP\Bigr).
\]
Use this to invert the equicorrelation matrix \( (1-\rho)\I+\rho\bone\bone\T \).
:::

::: {#exr-proj-ginverse-transpose}
[A4]

Show that if \( \G \) is a generalized inverse of the symmetric matrix \( \X\T\X \), then so is
\( \G\T \), and so is \( \G\X\T\X\G\T \). Which of the four Penrose conditions does the last one
satisfy?
:::

### B. Practice

::: {#exr-proj-product}
[B1]

Let \( \bP_1,\bP_2 \) be orthogonal projections. Show that \( \bP_1\bP_2 \) is an orthogonal
projection iff \( \bP_1\bP_2=\bP_2\bP_1 \), and that in this case it projects onto
\( \C(\bP_1)\cap\C(\bP_2) \). *Hint:* for the range, show
\( \bP_1\bP_2\bu=\bu \) iff \( \bu \) lies in both subspaces.
:::

::: {.solution}
If \( \bP_1\bP_2=\bP_2\bP_1 \), then
\( (\bP_1\bP_2)\T=\bP_2\bP_1=\bP_1\bP_2 \), and
\( (\bP_1\bP_2)^2=\bP_1\bP_2\bP_2\bP_1=\bP_1\bP_2\bP_1=\bP_1\bP_1\bP_2=\bP_1\bP_2 \). So the
product is symmetric and idempotent. Conversely, if \( \bP_1\bP_2 \) is an orthogonal
projection, it is symmetric, so \( \bP_1\bP_2=(\bP_1\bP_2)\T=\bP_2\bP_1 \). For the range: if
\( \bu \) lies in both subspaces, then \( \bP_1\bP_2\bu=\bP_1\bu=\bu \). If \( \bP_1\bP_2\bu=\bu \),
then \( \bu\in\C(\bP_1) \), and also \( \bu=\bP_2\bP_1\bu\in\C(\bP_2) \). A projection's range
is the set of vectors it fixes, so \( \C(\bP_1\bP_2)=\C(\bP_1)\cap\C(\bP_2) \).
:::

::: {#exr-proj-all-ginverses}
[B2]

Let \( \G_0 \) be a generalized inverse of \( \A \).

::: {.enumerate options="label=(\alph*)"}
1. Show that \( \G_0+(\I-\G_0\A)\bU+\V(\I-\A\G_0) \) is a generalized inverse of \( \A \) for all
           conformable \( \bU,\V \).

2. Show that every generalized inverse \( \G \) of \( \A \) has this form. *Hint:*
           take \( \bU=\G-\G_0 \) and \( \V=\G_0\A(\G-\G_0) \), and use \( \A(\G-\G_0)\A=\bzero \).

3. Conclude that a matrix has a unique generalized inverse iff it is square and
           nonsingular.
:::
:::

::: {.solution}
*Parts a,b.* For (a), multiply out:
\( \A\bigl[\G_0+(\I-\G_0\A)\bU+\V(\I-\A\G_0)\bigr]\A
=\A\G_0\A+(\A-\A\G_0\A)\bU\A+\A\V(\A-\A\G_0\A)=\A \). For (b), put \( \boldsymbol{\Delta}=\G-\G_0 \),
so that \( \A\boldsymbol{\Delta}\A=\A-\A=\bzero \). With \( \bU=\boldsymbol{\Delta} \) and \( \V=\G_0\A\boldsymbol{\Delta} \),
\[
(\I-\G_0\A)\boldsymbol{\Delta}+\G_0\A\boldsymbol{\Delta}(\I-\A\G_0)
=\boldsymbol{\Delta}-\G_0\A\boldsymbol{\Delta}+\G_0\A\boldsymbol{\Delta}-\G_0(\A\boldsymbol{\Delta}\A)\G_0=\boldsymbol{\Delta} ,
\]
so the formula reproduces \( \G_0+\boldsymbol{\Delta}=\G \).
:::

::: {#exr-proj-penrose}
[B3]

Using the singular value decomposition (@thm-mat-svd), show that \( (\X\T\X)^+\X\T=\X^+ \) and
\( \X\X^+=\M \). Show also that \( \X^+=(\X\T\X)^{-1}\X\T \) when \( \X \) has full column rank.
:::

### C. Going deeper

::: {#exr-proj-M-unique-characterization}
[C1]

Show that \( \M \) is the only \( n\times n \) matrix \( \bP \) satisfying \( \bP\X=\X \),
\( \bP\T=\bP \) and \( \rank(\bP)=\rank(\X) \). Is the rank condition necessary? Is symmetry?
:::
