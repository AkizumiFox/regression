# Kronecker products and the vec operator

Balanced designs have model matrices with a product structure, and multivariate
observations have covariance matrices with one. The Kronecker product expresses both.

::: {#def-mat-kronecker}
The **Kronecker product** of an \( m\times n \) matrix \( \A \) and a \( p\times q \) matrix \( \B \) is the
\( mp\times nq \) block matrix \( \A\otimes\B=(a_{ij}\B) \). The **vec operator** places the \( n \) columns of
an \( m\times n \) matrix one below another, giving a vector of length \( mn \): \( \vecop(\A)=(\mathbf{a}_1\T,\dots,\mathbf{a}_n\T)\T \).
:::

::: {#prp-mat-kronecker}
For matrices of compatible sizes:

::: {.enumerate options="label=(\alph*)"}
1. \( (\A\otimes\B)(\mathbf{C}\otimes\mathbf{D})=\A\mathbf{C}\otimes\B\mathbf{D} \) (mixed-product rule);

2. \( (\A\otimes\B)\T=\A\T\otimes\B\T \), and \( (\A\otimes\B)^{-1}=\A^{-1}\otimes\B^{-1} \) for
           nonsingular \( \A,\B \). More generally \( \A\ginv\otimes\B\ginv \) is a generalized inverse
           of \( \A\otimes\B \);

3. \( \tr(\A\otimes\B)=\tr(\A)\tr(\B) \) and \( \rank(\A\otimes\B)=\rank(\A)\rank(\B) \);

4. if \( \A \) (\( m\times m \)) and \( \B \) (\( p\times p \)) are symmetric with eigenvalues \( \lambda_i \) and
           \( \mu_j \), the eigenvalues of \( \A\otimes\B \) are the products \( \lambda_i\mu_j \), and
           \( \det(\A\otimes\B)=(\det\A)^p(\det\B)^m \);

5. Kronecker products of orthogonal, symmetric, idempotent, or nonnegative definite
           matrices are of the same kind;

6. \( \vecop(\A\B\mathbf{C})=(\mathbf{C}\T\otimes\A)\vecop(\B) \) and \( \tr(\A\T\B)=\vecop(\A)\T\vecop(\B) \).
:::

:::

::: {.proof}
(a) The \( (i,j) \) block of the left side is \( \sum_ka_{ik}c_{kj}\B\mathbf{D} \), which is the \( (i,j) \)
block of the right side. (b) Transposition is direct, and the rest follows from (a), for example
\( (\A\otimes\B)(\A\ginv\otimes\B\ginv)(\A\otimes\B)=\A\A\ginv\A\otimes\B\B\ginv\B \).
(d) With \( \A=\mathbf{P}\bLambda\mathbf{P}\T \) and \( \B=\Q\mathbf{M}\Q\T \), (a) and (b) give
\( \A\otimes\B=(\mathbf{P}\otimes\Q)(\bLambda\otimes\mathbf{M})(\mathbf{P}\otimes\Q)\T \), a spectral decomposition
because \( \mathbf{P}\otimes\Q \) is orthogonal and \( \bLambda\otimes\mathbf{M} \) is diagonal with entries \( \lambda_i\mu_j \).
The determinant is the product of these. (c) The trace statement is direct. For the rank,
by (a) \( (\A\otimes\B)(\A\ginv\otimes\B\ginv)=\A\A\ginv\otimes\B\B\ginv \) is idempotent with
the same rank as \( \A\otimes\B \) (@prp-mat-ginverse-props(a)). Its trace is
\( \tr(\A\A\ginv)\tr(\B\B\ginv)=\rank(\A)\rank(\B) \) by @prp-mat-idempotent-basic(b).
(e) follows from (a), (b) and (d). (f) The \( j \)th column of \( \A\B\mathbf{C} \) is
\( \A\B\mathbf{c}_j=\sum_kc_{kj}\A\mathbf{b}_k \), which is the \( j \)th block of
\( (\mathbf{C}\T\otimes\A)\vecop(\B) \). The trace identity is \( \sum_{ij}a_{ij}b_{ij} \) written two ways.
:::

The determinant formula in (d) holds for all square \( \A \) and \( \B \), not only symmetric ones. We
need only the symmetric case.

::: {#exm-mat-two-way-kronecker}
[The balanced two-way layout]

Let observations \( y_{ij} \), \( i=1,\dots,a \), \( j=1,\dots,b \), be stacked with \( j \) running
inside \( i \). The additive model \( \E(y_{ij})=\mu+\alpha_i+\beta_j \) has model matrix
\[
\X=\bigl[\bone_a\otimes\bone_b,\ \I_a\otimes\bone_b,\ \bone_a\otimes\I_b\bigr],
\]
because \( \I_a\otimes\bone_b \) has a \( 1 \) in column \( i \) exactly for the \( b \) observations in row
\( i \). With \( \bar{\mathbf{J}}_k=k^{-1}\mathbf{J}_k \), a symmetric idempotent averaging matrix, the mixed-product
rule shows that
\( \I_a\otimes\bar{\mathbf{J}}_b \), \( \bar{\mathbf{J}}_a\otimes\I_b \) and \( \bar{\mathbf{J}}_a\otimes\bar{\mathbf{J}}_b \)
replace each observation by its row mean, its column mean and the grand mean. The matrix
\( \M=\I_a\otimes\bar{\mathbf{J}}_b+\bar{\mathbf{J}}_a\otimes\I_b-\bar{\mathbf{J}}_a\otimes\bar{\mathbf{J}}_b \) is symmetric
and idempotent, with trace \( a+b-1=\rank(\X) \). It is the projection onto
\( \C(\X) \) (@thm-tw-additive in [Chapter 16](../ch16-multiway-layouts/index.html)). The script checks this for \( a=3 \) and
\( b=4 \), together with each part of @prp-mat-kronecker.
:::

```{.python .run #cell-kronecker-kron}
import numpy as np
rng = np.random.default_rng(11)
def vec(M):
    return M.reshape(-1, order="F")                     # stack the columns

A, B, C = rng.normal(size=(3, 4)), rng.normal(size=(4, 2)), rng.normal(size=(2, 5))
assert np.allclose(vec(A @ B @ C), np.kron(C.T, A) @ vec(B))          # vec(ABC)

# balanced two-way layout: a = 3 rows, b = 4 columns, one observation per cell
a, b = 3, 4
one = lambda k: np.ones((k, 1))
X = np.hstack([np.kron(one(a), one(b)),        # intercept
               np.kron(np.eye(a), one(b)),     # row indicators
               np.kron(one(a), np.eye(b))])    # column indicators
print("rank(X) =", np.linalg.matrix_rank(X), "= a + b - 1 =", a + b - 1)
```

## Exercises

### A. Check your understanding

::: {#exr-mat-kron-indicators}
[A1]

With observations stacked as in @exm-mat-two-way-kronecker, the row and column indicator
matrices are \( \mathbf{R}=\I_a\otimes\bone_b \) and \( \mathbf{C}=\bone_a\otimes\I_b \). Use the mixed-product
rule to compute \( \mathbf{R}\T\mathbf{R} \), \( \mathbf{C}\T\mathbf{C} \) and \( \mathbf{R}\T\mathbf{C} \), and identify the
cross-product matrix \( \X\T\X \) of the balanced additive model.
:::

::: {.solution}
By @prp-mat-kronecker(a) and (b),
\( \mathbf{R}\T\mathbf{R}=(\I_a\T\I_a)\otimes(\bone_b\T\bone_b)=b\,\I_a \) and
\( \mathbf{C}\T\mathbf{C}=(\bone_a\T\bone_a)\otimes\I_b=a\,\I_b \), while
\( \mathbf{R}\T\mathbf{C}=(\I_a\T\bone_a)\otimes(\bone_b\T\I_b)=\bone_a\bone_b\T \), the \( a\times b \) matrix
of ones. With \( \bone_{ab}=\bone_a\otimes\bone_b=\mathbf{R}\bone_a=\mathbf{C}\bone_b \), the cross-product
matrix of \( \X=[\bone_{ab},\mathbf{R},\mathbf{C}] \) is
\[
\X\T\X=\begin{pmatrix}ab&b\bone_a\T&a\bone_b\T\\
b\bone_a&b\I_a&\bone_a\bone_b\T\\
a\bone_b&\bone_b\bone_a\T&a\I_b\end{pmatrix}.
\]
Every cell of the layout is filled once, so each pair of a row level and a column level occurs
exactly once, which is what the block of ones records.
:::

::: {#exr-mat-kron-orthonormal}
[A2]

Show that \( (\A\otimes\B)\T(\A\otimes\B)=(\A\T\A)\otimes(\B\T\B) \), and deduce that the Kronecker
product of two matrices with orthonormal columns again has orthonormal columns. Show also that
\( \norm{\A}_F=\norm{\vecop(\A)} \), so that \( \vecop \) turns the Frobenius norm into the Euclidean
norm.
:::

::: {#exr-mat-kron-separable}
[A3]

Let \( \bSigma=\bSigma_1\otimes\bSigma_2 \) with \( \bSigma_1 \) (\( m\times m \)) and \( \bSigma_2 \)
(\( p\times p \)) positive definite — a *separable* covariance matrix, one variance structure
across rows and another across columns. Show that \( \bSigma \) is positive definite, that
\( \bSigma^{-1}=\bSigma_1^{-1}\otimes\bSigma_2^{-1} \), that
\( \bSigma^{1/2}=\bSigma_1^{1/2}\otimes\bSigma_2^{1/2} \) in the sense of @thm-mat-square-root, and
that \( \log\det\bSigma=p\log\det\bSigma_1+m\log\det\bSigma_2 \).
:::

### B. Practice

::: {#exr-mat-vec-multivariate}
[B1]

Let \( \mathbf{Y} \) be \( n\times q \), let \( \X \) be \( n\times p \) with full column rank and let
\( \mathbf{B} \) be \( p\times q \), so that \( q \) responses are regressed on the same \( \X \). Using
@prp-mat-kronecker(f), show that \( \vecop(\X\mathbf{B})=(\I_q\otimes\X)\vecop(\mathbf{B}) \) and that
minimizing \( \norm{\mathbf{Y}-\X\mathbf{B}}_F^2 \) leads to
\[
(\I_q\otimes\X\T\X)\vecop(\mathbf{B})=\vecop(\X\T\mathbf{Y}).
\]
Deduce that a multivariate regression on a common model matrix is \( q \) separate univariate
regressions, one per column of \( \mathbf{Y} \).
:::

::: {.solution}
Part (f) with \( \X\mathbf{B}\I_q \) gives \( \vecop(\X\mathbf{B})=(\I_q\T\otimes\X)\vecop(\mathbf{B}) \). By
@exr-mat-kron-orthonormal, \( \norm{\mathbf{Y}-\X\mathbf{B}}_F^2=\norm{\vecop(\mathbf{Y})-(\I_q\otimes\X)\vecop(\mathbf{B})}^2 \),
an ordinary least squares criterion with model matrix \( \I_q\otimes\X \). Its normal equations
are
\[
(\I_q\otimes\X)\T(\I_q\otimes\X)\vecop(\mathbf{B})=(\I_q\otimes\X)\T\vecop(\mathbf{Y}),
\]
and the mixed-product rule turns the left side into \( (\I_q\otimes\X\T\X)\vecop(\mathbf{B}) \) and the
right side into \( \vecop(\X\T\mathbf{Y}) \). The matrix \( \I_q\otimes\X\T\X \) is block diagonal with
identical blocks, so the system separates into \( \X\T\X\mathbf{b}_j=\X\T\y_j \) for each column
\( \y_j \) of \( \mathbf{Y} \). Nothing is gained by treating the responses jointly, as long as they
share a model matrix and the criterion is the unweighted sum of squares.
:::

::: {#exr-mat-kron-identifiability}
[B2]

Show that \( \A\otimes\B=\bzero \) iff \( \A=\bzero \) or \( \B=\bzero \), and that if
\( \A\otimes\B=\mathbf{C}\otimes\mathbf{D} \) with \( \A,\mathbf{C} \) of one size and \( \B,\mathbf{D} \) of another, all four
nonzero, then \( \mathbf{C}=c^{-1}\A \) and \( \mathbf{D}=c\B \) for some scalar \( c\ne0 \). A separable
covariance matrix therefore determines its two factors only up to reciprocal scalars.
:::

::: {.solution}
The \( (i,j) \) block of \( \A\otimes\B \) is \( a_{ij}\B \), so the product is zero iff every
\( a_{ij}\B=\bzero \). If \( \A\otimes\B=\mathbf{C}\otimes\mathbf{D} \), then \( a_{ij}\B=c_{ij}\mathbf{D} \) for all
\( i,j \). Choose \( k,l \) with \( a_{kl}\ne0 \). If \( c_{kl}=0 \) then \( \B=\bzero \), which is
excluded, so \( \mathbf{D}=c\B \) with \( c=a_{kl}/c_{kl}\ne0 \). Substituting back,
\( a_{ij}\B=c\,c_{ij}\B \) for all \( i,j \), and \( \B\ne\bzero \) forces \( a_{ij}=c\,c_{ij} \), that
is, \( \mathbf{C}=c^{-1}\A \). Conversely \( (c^{-1}\A)\otimes(c\B)=\A\otimes\B \).
:::

::: {#exr-mat-kron-anova}
[B3]

Write \( \bar{\mathbf{J}}_k=k^{-1}\mathbf{J}_k \) and \( \Q_k=\I_k-\bar{\mathbf{J}}_k \), both symmetric idempotent, of
ranks \( 1 \) and \( k-1 \). Expanding \( \I_{ab}=(\bar{\mathbf{J}}_a+\Q_a)\otimes(\bar{\mathbf{J}}_b+\Q_b) \), show
that the four matrices
\[
\bar{\mathbf{J}}_a\otimes\bar{\mathbf{J}}_b,\quad \Q_a\otimes\bar{\mathbf{J}}_b,\quad
\bar{\mathbf{J}}_a\otimes\Q_b,\quad \Q_a\otimes\Q_b
\]
are symmetric idempotent, sum to \( \I_{ab} \) and multiply to zero in pairs, with ranks
\( 1 \), \( a-1 \), \( b-1 \) and \( (a-1)(b-1) \). Show that the matrix \( \M \) of
@exm-mat-two-way-kronecker is the sum of the first three, so that \( \I-\M=\Q_a\otimes\Q_b \).
These are the grand mean, row, column and interaction parts of a two-way analysis of variance,
and their ranks are its degrees of freedom.
:::

::: {.solution}
Each of the four is a Kronecker product of symmetric idempotent matrices, hence symmetric
idempotent by @prp-mat-kronecker(e), with rank the product of the ranks by
@prp-mat-kronecker(c). Expanding the product of the two decompositions of the identity gives
\( \I_a\otimes\I_b \) as their sum. Any two of them have a factor \( \bar{\mathbf{J}}\Q=\bzero \) or
\( \Q\bar{\mathbf{J}}=\bzero \) in one of the two slots, so their product vanishes by the mixed-product
rule. The four ranks add to \( 1+(a-1)+(b-1)+(a-1)(b-1)=ab \), as @thm-mat-idempotent-sum
requires. Finally
\[
\M=\I_a\otimes\bar{\mathbf{J}}_b+\bar{\mathbf{J}}_a\otimes\I_b-\bar{\mathbf{J}}_a\otimes\bar{\mathbf{J}}_b
=(\bar{\mathbf{J}}_a+\Q_a)\otimes\bar{\mathbf{J}}_b+\bar{\mathbf{J}}_a\otimes(\bar{\mathbf{J}}_b+\Q_b)
-\bar{\mathbf{J}}_a\otimes\bar{\mathbf{J}}_b,
\]
which collects into \( \bar{\mathbf{J}}_a\otimes\bar{\mathbf{J}}_b+\Q_a\otimes\bar{\mathbf{J}}_b+\bar{\mathbf{J}}_a\otimes\Q_b \),
of rank \( 1+(a-1)+(b-1)=a+b-1 \). Subtracting from \( \I_{ab} \) leaves \( \Q_a\otimes\Q_b \).
:::

### C. Going deeper

::: {#exr-mat-commutation}
[C1]

Show that for each \( m,n \) there is exactly one \( mn\times mn \) matrix \( \mathbf{K}_{m,n} \) with
\( \mathbf{K}_{m,n}\vecop(\A)=\vecop(\A\T) \) for every \( m\times n \) matrix \( \A \), that it is a
permutation matrix, and that \( \mathbf{K}_{m,n}\T=\mathbf{K}_{m,n}^{-1}=\mathbf{K}_{n,m} \). Then show that for
\( \A \) of size \( m\times n \) and \( \B \) of size \( p\times q \),
\[
\mathbf{K}_{p,m}(\A\otimes\B)=(\B\otimes\A)\mathbf{K}_{q,n}.
\]
This **commutation matrix** is how the two orders of a Kronecker product are related, and how
derivatives with respect to a matrix are rearranged.
:::

::: {.solution}
The map \( \vecop(\A)\mapsto\vecop(\A\T) \) permutes the \( mn \) coordinates: the entry
\( a_{ij} \), which sits in position \( i+(j-1)m \) of \( \vecop(\A) \), sits in position
\( j+(i-1)n \) of \( \vecop(\A\T) \). A linear map on \( \Real^{mn} \) is determined by its values on
the coordinate vectors, and every coordinate vector is \( \vecop(\mathbf{e}_i\mathbf{e}_j\T) \) for some
\( i,j \), so \( \mathbf{K}_{m,n} \) exists and is unique, and it is the permutation matrix just
described. Applying the rule twice returns \( \vecop(\A) \), so
\( \mathbf{K}_{n,m}\mathbf{K}_{m,n}=\I \); a permutation matrix is orthogonal, so its inverse is its
transpose.

For the last identity, let \( \mathbf{C} \) be \( q\times n \), so that both sides may be applied to
\( \vecop(\mathbf{C}) \). By @prp-mat-kronecker(f), \( (\A\otimes\B)\vecop(\mathbf{C})=\vecop(\B\mathbf{C}\A\T) \),
a \( p\times m \) matrix, so the left side sends \( \vecop(\mathbf{C}) \) to
\( \vecop\bigl((\B\mathbf{C}\A\T)\T\bigr)=\vecop(\A\mathbf{C}\T\B\T) \). The right side sends it to
\( (\B\otimes\A)\vecop(\mathbf{C}\T)=\vecop(\A\mathbf{C}\T\B\T) \) by the same rule. The two matrices agree
on every \( \vecop(\mathbf{C}) \), that is, on all of \( \Real^{nq} \), so they are equal.
:::

::: {#exr-mat-two-way-random}
[C2]

In a balanced two-way layout with random row and column effects, the response vector has
covariance
\[
\V=\sigma_\alpha^2(\I_a\otimes\mathbf{J}_b)+\sigma_\beta^2(\mathbf{J}_a\otimes\I_b)+\sigma^2\I_{ab},
\qquad \sigma^2>0 .
\]
Using the four idempotent matrices of @exr-mat-kron-anova, show that
\[
\V=(b\sigma_\alpha^2+a\sigma_\beta^2+\sigma^2)\,\bar{\mathbf{J}}_a\otimes\bar{\mathbf{J}}_b
+(b\sigma_\alpha^2+\sigma^2)\,\Q_a\otimes\bar{\mathbf{J}}_b
+(a\sigma_\beta^2+\sigma^2)\,\bar{\mathbf{J}}_a\otimes\Q_b
+\sigma^2\,\Q_a\otimes\Q_b ,
\]
read off the eigenvalues of \( \V \) with their multiplicities, and write down \( \det\V \) and
\( \V^{-1} \). On which subspace does \( \V \) act as \( \sigma^2\I \)?
:::

::: {.solution}
Since \( \mathbf{J}_k=k\bar{\mathbf{J}}_k \) and \( \I_k=\bar{\mathbf{J}}_k+\Q_k \),
\[
\I_a\otimes\mathbf{J}_b=b(\bar{\mathbf{J}}_a\otimes\bar{\mathbf{J}}_b+\Q_a\otimes\bar{\mathbf{J}}_b),\qquad
\mathbf{J}_a\otimes\I_b=a(\bar{\mathbf{J}}_a\otimes\bar{\mathbf{J}}_b+\bar{\mathbf{J}}_a\otimes\Q_b),
\]
and \( \I_{ab} \) is the sum of all four. Collecting coefficients gives the displayed form. The
four matrices are symmetric idempotent with pairwise products zero and sum \( \I_{ab} \), so by
@thm-mat-idempotent(a) the columns of an orthonormal basis of each column space are
eigenvectors of \( \V \), with eigenvalue the coefficient in front. The eigenvalues are therefore
\( b\sigma_\alpha^2+a\sigma_\beta^2+\sigma^2 \) once, \( b\sigma_\alpha^2+\sigma^2 \) with
multiplicity \( a-1 \), \( a\sigma_\beta^2+\sigma^2 \) with multiplicity \( b-1 \) and \( \sigma^2 \)
with multiplicity \( (a-1)(b-1) \). The determinant is their product,
\[
\det\V=(b\sigma_\alpha^2+a\sigma_\beta^2+\sigma^2)(b\sigma_\alpha^2+\sigma^2)^{a-1}
(a\sigma_\beta^2+\sigma^2)^{b-1}(\sigma^2)^{(a-1)(b-1)},
\]
and the inverse replaces each coefficient in the displayed sum by its reciprocal, since the
four terms multiply to zero in pairs and each is idempotent. Finally
\( \V\bu=\sigma^2\bu \) for every \( \bu\in\C(\Q_a\otimes\Q_b) \): on the interaction space the
row and column variance components have no effect at all, which is why that space supplies the
error term of the balanced two-way analysis of variance.
:::
