# Vectors, matrices and the operations we use

This chapter works in \( \Real^n \) and with real matrices only. It fixes the conventions,
states the results later chapters cite, and proves them. A reader who already knows
linear algebra can go straight to the results that are specific to statistics:
generalized inverses ([Section 1.3](03-inverses.html) and [Section 1.4](04-systems.html)), Schur complements and the
Woodbury identity ([Section 1.6](06-partitioned.html)), nonnegative definite matrices
([Section 1.7](07-eigen.html)), conditioning ([Section 1.8](08-svd.html)) and idempotent matrices
([Section 1.9](09-idempotent.html)). Orthogonal projections, the geometric heart of least
squares, get only the minimum here. [Chapter 6](../ch06-projections/index.html) develops them.

## Conventions

An \( m\times n \) matrix \( \A=(a_{ij}) \) has \( m \) rows and \( n \) columns. We write
\( \A=[\mathbf{a}_1,\dots,\mathbf{a}_n] \) for its columns and \( \mathbf{a}_{(i)}\T \) for its \( i \)th row. A
**vector** is always a column, so \( \x\in\Real^n \) is an \( n\times1 \) matrix and \( \x\T \) is
a row. Bold upright letters denote vectors and matrices, and italic letters denote
scalars. The special matrices we use constantly are the identity \( \I_n \), the zero matrix
\( \bzero \), the vector of ones \( \bone_n \), the matrix of ones \( \mathbf{J}_n=\bone_n\bone_n\T \), and
the coordinate vectors \( \mathbf{e}_1,\dots,\mathbf{e}_n \), the columns of \( \I_n \). Subscripts giving
sizes are dropped when the size is clear. \( \diag(d_1,\dots,d_n) \) is the diagonal matrix
with the given diagonal entries. A square matrix is **symmetric** if \( \A\T=\A \),
*upper triangular* if \( a_{ij}=0 \) for \( i>j \), and *lower triangular* if
\( \A\T \) is upper triangular.

The **transpose** satisfies \( (\A\T)\T=\A \), \( (\A+\B)\T=\A\T+\B\T \) and, most
usefully,
\[
(\A\B)\T=\B\T\A\T .
\]{#eq-mat-transpose-product}

## Products as linear combinations

If \( \A \) is \( m\times n \) and \( \B \) is \( n\times q \), the \( (i,j) \) entry of \( \A\B \) is
\( \sum_k a_{ik}b_{kj} \). Three equivalent readings of the product carry most of the
arguments in this book:
\[
\A\bb=\sum_{j=1}^n b_j\mathbf{a}_j,\qquad
\A\B=[\A\bb_1,\dots,\A\bb_q],\qquad
\A\B=\sum_{k=1}^n\mathbf{a}_k\mathbf{b}_{(k)}\T .
\]{#eq-mat-product-readings}

The first says that \( \A\bb \) combines the columns of \( \A \), with weights taken from
\( \bb \). The second says that each column of \( \A\B \) is such a combination.
The third writes a product as a sum of \( n \) rank-one matrices (outer products), and is the
form in which spectral and singular value decompositions are usually used. Transposing
the first reading shows that each row of \( \A\B \) is a combination of the rows of \( \B \).

Matrices partitioned into conformable blocks multiply as if the blocks were scalars, as
long as the order of factors is kept:
\[
\begin{pmatrix}\A_{11}&\A_{12}\\\A_{21}&\A_{22}\end{pmatrix}
\begin{pmatrix}\B_{11}&\B_{12}\\\B_{21}&\B_{22}\end{pmatrix}
=\begin{pmatrix}\A_{11}\B_{11}+\A_{12}\B_{21}&\A_{11}\B_{12}+\A_{12}\B_{22}\\
\A_{21}\B_{11}+\A_{22}\B_{21}&\A_{21}\B_{12}+\A_{22}\B_{22}\end{pmatrix}.
\]
Multiplication is associative and distributive but not commutative, and \( \A\B=\bzero \)
does not force \( \A=\bzero \) or \( \B=\bzero \). The **Hadamard product**
\( \A\circ\B=(a_{ij}b_{ij}) \) multiplies entrywise. It appears only occasionally.

## Inner products, lengths and the Cauchy–Schwarz inequality

The **inner product** of \( \bu,\bv\in\Real^n \) is \( \inner{\bu}{\bv}=\bu\T\bv=\sum_iu_iv_i \),
and the **Euclidean norm** is \( \norm{\bu}=(\bu\T\bu)^{1/2} \). The entries of \( \A\T\A \) are
the inner products of the columns of \( \A \), and those of \( \A\A\T \) are the inner products of
its rows. In particular
\[
\sum_{i,j}a_{ij}^2=\sum_j\norm{\mathbf{a}_j}^2 ,
\]{#eq-mat-frobenius}

so \( \A\T\A=\bzero \) forces every column to be zero, that is, \( \A=\bzero \). The square root of
@eq-mat-frobenius is the **Frobenius norm** \( \norm{\A}_F \).

::: {#prp-mat-cauchy-schwarz}
[Cauchy–Schwarz]

For \( \bu,\bv\in\Real^n \), \( (\bu\T\bv)^2\le(\bu\T\bu)(\bv\T\bv) \), with equality iff one
vector is a multiple of the other. Consequently
\( \norm{\bu+\bv}\le\norm{\bu}+\norm{\bv} \).
:::

::: {.proof}
If \( \bv=\bzero \) there is nothing to prove. Otherwise, the quadratic
\( q(t)=\norm{\bu-t\bv}^2=\bu\T\bu-2t\,\bu\T\bv+t^2\bv\T\bv \) is nonnegative for all \( t \), so
its discriminant \( 4(\bu\T\bv)^2-4(\bu\T\bu)(\bv\T\bv) \) is at most zero. Equality holds iff
\( q \) has a real root, that is, iff \( \bu=t\bv \) for some \( t \). The triangle inequality
follows by expanding \( \norm{\bu+\bv}^2 \) and bounding the cross term.
:::

## Quadratic forms

A **quadratic form** in \( \x\in\Real^n \) is a function
\( \x\T\A\x=\sum_{i,j}a_{ij}x_ix_j \). Since \( \x\T\A\x \) is a scalar it equals its own
transpose \( \x\T\A\T\x \), so
\[
\x\T\A\x=\x\T\Bigl(\tfrac12(\A+\A\T)\Bigr)\x .
\]{#eq-mat-symmetrize}

Every quadratic form therefore has a symmetric matrix, and that matrix is unique: if
\( \x\T\bS\x=0 \) for all \( \x \) with \( \bS \) symmetric, then taking \( \x=\mathbf{e}_i \) gives
\( s_{ii}=0 \), and \( \x=\mathbf{e}_i+\mathbf{e}_j \) gives \( 2s_{ij}=0 \). From now on we always take the matrix
in a quadratic form to be symmetric.

::: {#exm-mat-centring}
[Sums of squares as quadratic forms]

The sum of squared deviations of \( \y\in\Real^n \) from its mean is
\[
\sum_{i=1}^n(y_i-\bar{y})^2=\y\T\y-n\bar{y}^2
=\y\T\Bigl(\I_n-\tfrac1n\mathbf{J}_n\Bigr)\y ,
\]
because \( n\bar{y}^2=(\bone\T\y)^2/n=\y\T\mathbf{J}_n\y/n \). The **centring matrix**
\( \I_n-n^{-1}\mathbf{J}_n \) is symmetric and satisfies
\( (\I_n-n^{-1}\mathbf{J}_n)^2=\I_n-n^{-1}\mathbf{J}_n \), because \( \mathbf{J}_n^2=n\mathbf{J}_n \). Every sum of squares
in an analysis of variance is a quadratic form whose matrix has these two properties
([Section 1.9](09-idempotent.html)).
:::

## Exercises

### A. Check your understanding

::: {#exr-mat-readings}
[A1]

Let \( \A \) be \( m\times n \) and let \( \bD=\diag(d_1,\dots,d_n) \). Using
@eq-mat-product-readings, identify \( \A\mathbf{e}_j \), \( \mathbf{e}_i\T\A \) and \( \mathbf{e}_i\T\A\mathbf{e}_j \), and show
that \( \A\bD \) multiplies the \( j \)th column of \( \A \) by \( d_j \), while a diagonal matrix
multiplying \( \A \) on the left scales its rows. What are \( \A\bone_n \), \( \bone_m\T\A \) and
\( \bone_m\T\A\bone_n \)?
:::

::: {#exr-mat-polarization}
[A2]

Show the parallelogram law \( \norm{\bu+\bv}^2+\norm{\bu-\bv}^2=2\norm{\bu}^2+2\norm{\bv}^2 \) and the
polarization identity \( \bu\T\bv=\tfrac14\bigl(\norm{\bu+\bv}^2-\norm{\bu-\bv}^2\bigr) \). Deduce that a
matrix \( \bT \) with \( \norm{\bT\x}=\norm{\x} \) for every \( \x \) satisfies
\( (\bT\bu)\T(\bT\bv)=\bu\T\bv \) for every \( \bu,\bv \): preserving lengths preserves angles.
:::

::: {.solution}
Expanding \( \norm{\bu\pm\bv}^2=\norm{\bu}^2\pm2\bu\T\bv+\norm{\bv}^2 \) and adding gives the first
identity; subtracting gives the second. If \( \bT \) preserves lengths, then both norms on the
right of the polarization identity are unchanged when \( \bu,\bv \) are replaced by
\( \bT\bu,\bT\bv \), since \( \bT\bu\pm\bT\bv=\bT(\bu\pm\bv) \). So \( (\bT\bu)\T(\bT\bv)=\bu\T\bv \).
:::

::: {#exr-mat-mean-bound}
[A3]

Apply @prp-mat-cauchy-schwarz with \( \bu=\bone_n \) and \( \bv=\y \) to show
\( \bigl(\sum_iy_i\bigr)^2\le n\sum_iy_i^2 \), with equality iff all the \( y_i \) are equal. Which
property of the centring matrix of @exm-mat-centring is this?
:::

::: {.solution}
Cauchy–Schwarz gives \( (\bone\T\y)^2\le(\bone\T\bone)(\y\T\y)=n\y\T\y \), with equality iff \( \y \)
is a multiple of \( \bone \). Rearranged, this says \( \y\T\y-n\bar{y}^2\ge0 \), that is, the quadratic
form of the centring matrix \( \I_n-n^{-1}\mathbf{J}_n \) is nonnegative, and vanishes only for constant
\( \y \). The sum of squared deviations is never negative, and is zero only when there is no
variation at all.
:::

### B. Practice

::: {#exr-mat-skew}
[B1]

Let \( \A \) be \( n\times n \). (a) Show that \( \x\T\A\bv=0 \) for all \( \x,\bv \) iff \( \A=\bzero \).
(b) Show that \( \x\T\A\x=0 \) for all \( \x \) iff \( \A\T=-\A \) (such an \( \A \) is called
*skew-symmetric*), and give a nonzero \( 2\times2 \) example. (c) Deduce that two matrices define
the same quadratic form iff they differ by a skew-symmetric matrix, so that the symmetric matrix
in @eq-mat-symmetrize is the only symmetric one.
:::

::: {.solution}
(a) Taking \( \x=\mathbf{e}_i \) and \( \bv=\mathbf{e}_j \) gives \( a_{ij}=0 \); the converse is clear.
(b) By @eq-mat-symmetrize, \( \x\T\A\x=\x\T\bS\x \) with \( \bS=\frac12(\A+\A\T) \) symmetric, and the
uniqueness argument of this section shows that \( \x\T\bS\x\equiv0 \) with \( \bS \) symmetric forces
\( \bS=\bzero \), that is, \( \A\T=-\A \). Conversely \( \x\T\A\x=\frac12\x\T(\A+\A\T)\x=0 \) for
skew-symmetric \( \A \). The matrix \( \begin{psmallmatrix}0&1\\-1&0\end{psmallmatrix} \) is an
example: it rotates \( \x \) by a quarter turn, and a vector is orthogonal to its own rotation.
(c) \( \x\T\A\x=\x\T\B\x \) for all \( \x \) iff \( \A-\B \) is skew-symmetric by (b). The
matrices with a given quadratic form are therefore \( \bS+\mathbf{K} \) with \( \mathbf{K} \)
skew-symmetric, and only \( \mathbf{K}=\bzero \) gives a symmetric one, since a matrix that is both
symmetric and skew-symmetric is zero.
:::

::: {#exr-mat-frobenius-submult}
[B2]

Show that \( \norm{\A\B}_F\le\norm{\A}_F\norm{\B}_F \) whenever the product is defined, and in
particular \( \norm{\A\x}\le\norm{\A}_F\norm{\x} \). For fixed \( \A\ne\bzero \), when does the second
inequality become an equality?
:::

::: {.solution}
The \( (i,j) \) entry of \( \A\B \) is \( \mathbf{a}_{(i)}\T\mathbf{b}_j \), so @prp-mat-cauchy-schwarz gives
\( (\mathbf{a}_{(i)}\T\mathbf{b}_j)^2\le\norm{\mathbf{a}_{(i)}}^2\norm{\mathbf{b}_j}^2 \). Summing over \( i \) and \( j \) and
using @eq-mat-frobenius for the rows of \( \A \) and the columns of \( \B \),
\[
\norm{\A\B}_F^2\le\Bigl(\sum_i\norm{\mathbf{a}_{(i)}}^2\Bigr)\Bigl(\sum_j\norm{\mathbf{b}_j}^2\Bigr)
=\norm{\A}_F^2\norm{\B}_F^2 .
\]
Taking \( \B=\x \) gives the second inequality. Equality there needs equality in Cauchy–Schwarz
for every row, that is, every row of \( \A \) is a multiple of \( \x\T \). So equality holds iff
\( \A=\mathbf{c}\x\T \) for some vector \( \mathbf{c} \), and in particular only for \( \A \) of rank one.
:::

::: {#exr-mat-outer-rows}
[B3]

Let \( \X \) be \( n\times p \) with rows \( \x_{(i)}\T \), and let \( \y\in\Real^n \). Use the third
reading in @eq-mat-product-readings to show
\[
\X\T\X=\sum_{i=1}^n\x_{(i)}\x_{(i)}\T,\qquad \X\T\y=\sum_{i=1}^ny_i\x_{(i)} .
\]
Deduce that deleting the \( i \)th observation subtracts \( \x_{(i)}\x_{(i)}\T \) and
\( y_i\x_{(i)} \), and that splitting the rows into groups adds the groups' cross-products. This is
the form in which @exm-mat-deletion updates a fit.
:::

::: {.solution}
The columns of \( \X\T \) are the vectors \( \x_{(i)} \) and the rows of \( \X \) are the \( \x_{(i)}\T \),
so the third reading of the product writes \( \X\T\X \) as the sum of the outer products
\( \x_{(i)}\x_{(i)}\T \). The same reading applied to \( \X\T\y \), with \( \y \) a single column,
gives \( \sum_iy_i\x_{(i)} \). Both are sums over observations of terms that involve one
observation each, so removing an observation removes its term and pooling two groups adds their
sums.
:::

### C. Going deeper

::: {#exr-mat-lagrange-identity}
[C1]

Prove Lagrange's identity
\[
\Bigl(\sum_iu_i^2\Bigr)\Bigl(\sum_iv_i^2\Bigr)-\Bigl(\sum_iu_iv_i\Bigr)^2
=\sum_{i<j}(u_iv_j-u_jv_i)^2 ,
\]
and read off from it a second proof of @prp-mat-cauchy-schwarz together with its equality case.
Deduce that \( \norm{\bu\bv\T-\bv\bu\T}_F^2=2\bigl(\norm{\bu}^2\norm{\bv}^2-(\bu\T\bv)^2\bigr) \), and
show that equality in the triangle inequality, \( \norm{\bu+\bv}=\norm{\bu}+\norm{\bv} \), holds iff
one of the two vectors is a *nonnegative* multiple of the other.
:::

::: {.solution}
Write the left-hand side as a double sum and symmetrize:
\[
\sum_{i,j}u_i^2v_j^2-\sum_{i,j}u_iv_iu_jv_j
=\tfrac12\sum_{i,j}\bigl(u_i^2v_j^2+u_j^2v_i^2-2u_iv_ju_jv_i\bigr)
=\tfrac12\sum_{i,j}(u_iv_j-u_jv_i)^2 ,
\]
and the terms with \( i=j \) vanish while each pair \( i\ne j \) is counted twice, which gives the
sum over \( i<j \). The right-hand side is nonnegative, so
\( (\bu\T\bv)^2\le\norm{\bu}^2\norm{\bv}^2 \), with equality iff \( u_iv_j=u_jv_i \) for all \( i,j \),
which for \( \bv\ne\bzero \) means \( \bu=(\bu\T\bv/\norm{\bv}^2)\bv \), a multiple of \( \bv \). The entries of
\( \bu\bv\T-\bv\bu\T \) are \( u_iv_j-u_jv_i \), so its squared Frobenius norm is the full double sum,
namely twice the sum over \( i<j \).

For the triangle inequality, squaring \( \norm{\bu+\bv}=\norm{\bu}+\norm{\bv} \) gives
\( \bu\T\bv=\norm{\bu}\norm{\bv} \), so Cauchy–Schwarz holds with equality and \( \bu=t\bv \) for some
\( t \) (or \( \bv=\bzero \)). Then \( \bu\T\bv=t\norm{\bv}^2 \) and
\( \norm{\bu}\norm{\bv}=\lvert t\rvert\norm{\bv}^2 \), so \( t\ge0 \). Conversely \( \bu=t\bv \) with
\( t\ge0 \) gives \( \norm{\bu+\bv}=(1+t)\norm{\bv}=\norm{\bu}+\norm{\bv} \).
:::

::: {#exr-mat-nearest-symmetric}
[C2]

For \( m\times n \) matrices put \( \inner{\A}{\B}_F=\sum_{i,j}a_{ij}b_{ij} \), so that
\( \inner{\A}{\A}_F=\norm{\A}_F^2 \) by @eq-mat-frobenius. Show that this is the ordinary inner
product of the two matrices read as vectors of length \( mn \), so that @prp-mat-cauchy-schwarz
gives \( \lvert\inner{\A}{\B}_F\rvert\le\norm{\A}_F\norm{\B}_F \). Now let \( \A \) be square and write
\( \A=\bS+\mathbf{K} \) with \( \bS=\frac12(\A+\A\T) \) and \( \mathbf{K}=\frac12(\A-\A\T) \). Show that
\( \inner{\bS_0}{\mathbf{K}_0}_F=0 \) whenever \( \bS_0 \) is symmetric and \( \mathbf{K}_0 \) is skew-symmetric, and
deduce that \( \bS \) is the unique symmetric matrix closest to \( \A \) in Frobenius norm. So the
symmetrization in @eq-mat-symmetrize is not merely convenient: it drops exactly the part of
\( \A \) that is orthogonal to every symmetric matrix.
:::

::: {.solution}
Listing the entries of a matrix in any fixed order turns \( \inner{\A}{\B}_F \) into
\( \sum_ka_kb_k \) for the two resulting vectors, so all the properties of the Euclidean inner
product, Cauchy–Schwarz among them, carry over.

If \( \bS_0\T=\bS_0 \) and \( \mathbf{K}_0\T=-\mathbf{K}_0 \), then relabelling the summation indices gives
\[
\inner{\bS_0}{\mathbf{K}_0}_F=\sum_{i,j}s_{ij}k_{ij}=\sum_{i,j}s_{ji}k_{ji}
=\sum_{i,j}s_{ij}(-k_{ij})=-\inner{\bS_0}{\mathbf{K}_0}_F ,
\]
so the inner product is zero. Let \( \bS_0 \) be any symmetric matrix. Then
\( \A-\bS_0=(\bS-\bS_0)+\mathbf{K} \) with the first term symmetric and the second skew-symmetric, so by
the orthogonality just proved
\[
\norm{\A-\bS_0}_F^2=\norm{\bS-\bS_0}_F^2+\norm{\mathbf{K}}_F^2\ge\norm{\mathbf{K}}_F^2=\norm{\A-\bS}_F^2 ,
\]
with equality iff \( \bS_0=\bS \). The decomposition \( \A=\bS+\mathbf{K} \) is an orthogonal one, and
\( \bS \) is the orthogonal projection of \( \A \) onto the subspace of symmetric matrices, in the
sense that [Chapter 6](../ch06-projections/index.html) develops for subspaces of \( \Real^n \).
:::
