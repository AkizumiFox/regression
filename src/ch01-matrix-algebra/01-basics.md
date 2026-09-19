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
\( \A=[\bm a_1,\dots,\bm a_n] \) for its columns and \( \bm a_{(i)}\T \) for its \( i \)th row. A
**vector** is always a column, so \( \x\in\Real^n \) is an \( n\times1 \) matrix and \( \x\T \) is
a row. Bold upright letters denote vectors and matrices, and italic letters denote
scalars. The special matrices we use constantly are the identity \( \I_n \), the zero matrix
\( \bzero \), the vector of ones \( \bone_n \), the matrix of ones \( \bm J_n=\bone_n\bone_n\T \), and
the coordinate vectors \( \bm e_1,\dots,\bm e_n \), the columns of \( \I_n \). Subscripts giving
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
\A\bb=\sum_{j=1}^n b_j\bm a_j,\qquad
\A\B=[\A\bb_1,\dots,\A\bb_q],\qquad
\A\B=\sum_{k=1}^n\bm a_k\bm b_{(k)}\T .
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
\sum_{i,j}a_{ij}^2=\sum_j\norm{\bm a_j}^2 ,
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
\( \x\T\bS\x=0 \) for all \( \x \) with \( \bS \) symmetric, then taking \( \x=\bm e_i \) gives
\( s_{ii}=0 \), and \( \x=\bm e_i+\bm e_j \) gives \( 2s_{ij}=0 \). From now on we always take the matrix
in a quadratic form to be symmetric.

::: {#exm-mat-centring}
[Sums of squares as quadratic forms]

The sum of squared deviations of \( \y\in\Real^n \) from its mean is
\[
\sum_{i=1}^n(y_i-\bar{y})^2=\y\T\y-n\bar{y}^2
=\y\T\Bigl(\I_n-\tfrac1n\bm J_n\Bigr)\y ,
\]
because \( n\bar{y}^2=(\bone\T\y)^2/n=\y\T\bm J_n\y/n \). The **centring matrix**
\( \I_n-n^{-1}\bm J_n \) is symmetric and satisfies
\( (\I_n-n^{-1}\bm J_n)^2=\I_n-n^{-1}\bm J_n \), because \( \bm J_n^2=n\bm J_n \). Every sum of squares
in an analysis of variance is a quadratic form whose matrix has these two properties
([Section 1.9](09-idempotent.html)).
:::
