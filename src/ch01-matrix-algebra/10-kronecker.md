# Kronecker products and the vec operator

Balanced designs have model matrices with a product structure, and multivariate
observations have covariance matrices with one. The Kronecker product expresses both.

::: {#def-mat-kronecker}
The **Kronecker product** of an \( m\times n \) matrix \( \A \) and a \( p\times q \) matrix \( \B \) is the
\( mp\times nq \) block matrix \( \A\otimes\B=(a_{ij}\B) \). The **vec operator** places the \( n \) columns of
an \( m\times n \) matrix one below another, giving a vector of length \( mn \): \( \vecop(\A)=(\bm a_1\T,\dots,\bm a_n\T)\T \).
:::

::: {#prp-mat-kronecker}
For matrices of compatible sizes:

::: {.enumerate options="label=(\alph*)"}
1. \( (\A\otimes\B)(\bm C\otimes\bm D)=\A\bm C\otimes\B\bm D \) (mixed-product rule);

2. \( (\A\otimes\B)\T=\A\T\otimes\B\T \), and \( (\A\otimes\B)^{-1}=\A^{-1}\otimes\B^{-1} \) for
           nonsingular \( \A,\B \). More generally \( \A\ginv\otimes\B\ginv \) is a generalized inverse
           of \( \A\otimes\B \);

3. \( \tr(\A\otimes\B)=\tr(\A)\tr(\B) \) and \( \rank(\A\otimes\B)=\rank(\A)\rank(\B) \);

4. if \( \A \) (\( m\times m \)) and \( \B \) (\( p\times p \)) are symmetric with eigenvalues \( \lambda_i \) and
           \( \mu_j \), the eigenvalues of \( \A\otimes\B \) are the products \( \lambda_i\mu_j \), and
           \( \det(\A\otimes\B)=(\det\A)^p(\det\B)^m \);

5. Kronecker products of orthogonal, symmetric, idempotent, or nonnegative definite
           matrices are of the same kind;

6. \( \vecop(\A\B\bm C)=(\bm C\T\otimes\A)\vecop(\B) \) and \( \tr(\A\T\B)=\vecop(\A)\T\vecop(\B) \).
:::

:::

::: {.proof}
(a) The \( (i,j) \) block of the left side is \( \sum_ka_{ik}c_{kj}\B\bm D \), which is the \( (i,j) \)
block of the right side. (b) Transposition is direct, and the rest follows from (a), for example
\( (\A\otimes\B)(\A\ginv\otimes\B\ginv)(\A\otimes\B)=\A\A\ginv\A\otimes\B\B\ginv\B \).
(d) With \( \A=\bm P\bLambda\bm P\T \) and \( \B=\Q\bm M\Q\T \), (a) and (b) give
\( \A\otimes\B=(\bm P\otimes\Q)(\bLambda\otimes\bm M)(\bm P\otimes\Q)\T \), a spectral decomposition
because \( \bm P\otimes\Q \) is orthogonal and \( \bLambda\otimes\bm M \) is diagonal with entries \( \lambda_i\mu_j \).
The determinant is the product of these. (c) The trace statement is direct. For the rank,
by (a) \( (\A\otimes\B)(\A\ginv\otimes\B\ginv)=\A\A\ginv\otimes\B\B\ginv \) is idempotent with
the same rank as \( \A\otimes\B \) (@prp-mat-ginverse-props(a)). Its trace is
\( \tr(\A\A\ginv)\tr(\B\B\ginv)=\rank(\A)\rank(\B) \) by @prp-mat-idempotent-basic(b).
(e) follows from (a), (b) and (d). (f) The \( j \)th column of \( \A\B\bm C \) is
\( \A\B\bm c_j=\sum_kc_{kj}\A\bm b_k \), which is the \( j \)th block of
\( (\bm C\T\otimes\A)\vecop(\B) \). The trace identity is \( \sum_{ij}a_{ij}b_{ij} \) written two ways.
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
\( i \). With \( \bar{\bm J}_k=k^{-1}\bm J_k \), a symmetric idempotent averaging matrix, the mixed-product
rule shows that
\( \I_a\otimes\bar{\bm J}_b \), \( \bar{\bm J}_a\otimes\I_b \) and \( \bar{\bm J}_a\otimes\bar{\bm J}_b \)
replace each observation by its row mean, its column mean and the grand mean. The matrix
\( \M=\I_a\otimes\bar{\bm J}_b+\bar{\bm J}_a\otimes\I_b-\bar{\bm J}_a\otimes\bar{\bm J}_b \) is symmetric
and idempotent, with trace \( a+b-1=\rank(\X) \). It is the projection onto \( \C(\X) \)
(Chapter 15). The script checks this for \( a=3 \) and
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
