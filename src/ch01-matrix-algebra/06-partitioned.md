# Trace, determinants and partitioned matrices

## Trace

The **trace** of a square matrix is the sum of its diagonal entries,
\( \tr(\A)=\sum_ia_{ii} \). It is linear, \( \tr(a\A+b\B)=a\tr(\A)+b\tr(\B) \), and
\( \tr(\A\T)=\tr(\A) \).

::: {#thm-mat-trace-cyclic}
[Cyclic property]

If \( \A \) is \( m\times n \) and \( \B \) is \( n\times m \), then \( \tr(\A\B)=\tr(\B\A) \).
:::

::: {.proof}
\( \tr(\A\B)=\sum_{i=1}^m\sum_{k=1}^na_{ik}b_{ki}=\sum_{k=1}^n\sum_{i=1}^mb_{ki}a_{ik}=\tr(\B\A) \).
:::

The two products need not have the same size. Applied repeatedly, the theorem allows
cyclic permutations such as \( \tr(\A\B\mathbf{C})=\tr(\mathbf{C}\A\B)=\tr(\B\mathbf{C}\A) \), but not arbitrary
reorderings. Three consequences occur again and again:
\[
\tr(\A\T\A)=\sum_{i,j}a_{ij}^2=\norm{\A}_F^2,\qquad
\x\T\A\x=\tr(\A\x\x\T),\qquad
\tr(\mathbf{P}^{-1}\A\mathbf{P})=\tr(\A).
\]{#eq-mat-trace-tricks}

The middle one turns a scalar quadratic form into the trace of a matrix that is linear in
\( \x\x\T \). Taking expectations then gives the mean of a random quadratic form
([Chapter 2](../ch02-random-vectors/index.html)).

## Determinants

The **determinant** of an \( n\times n \) matrix is
\[
\det\A=\sum_{\pi}\operatorname{sgn}(\pi)\,a_{1\pi(1)}a_{2\pi(2)}\cdots a_{n\pi(n)},
\]
the sum over all permutations \( \pi \) of \( \{1,\dots,n\} \), where \( \operatorname{sgn}(\pi)=\pm1 \)
according to whether \( \pi \) is a product of an even or odd number of transpositions. We
collect the standard properties. Parts (a)–(d) are proved in any linear algebra text
(Harville 1997, chapter 13), and the others follow from them.

::: {#prp-mat-det}
[Determinants]

Let \( \A,\B \) be \( n\times n \).

::: {.enumerate options="label=(\alph*)"}
1. \( \det\A\T=\det\A \).

2. \( \det\A \) is a linear function of each column separately. It changes sign when two
           columns are interchanged, and is unchanged when a multiple of one column is
           added to another. The same holds for rows.

3. A triangular matrix has determinant equal to the product of its diagonal entries.

4. \( \det(\A\B)=\det\A\det\B \).

5. \( \A \) is nonsingular iff \( \det\A\neq0 \), and then \( \det\A^{-1}=1/\det\A \).

6. \( \det(c\A)=c^n\det\A \), and \( \det\Q=\pm1 \) for orthogonal \( \Q \).
:::

:::

::: {.proof}
[Proof of (e) and (f)]

If \( \A \) is nonsingular, \( 1=\det\I=\det\A\det\A^{-1} \) by (d). If \( \A \) is singular, some
column is a combination of the others. Subtracting that combination leaves a zero column
without changing the determinant, by (b), and the determinant of a matrix with a zero
column is zero, since every term of the sum contains a factor from that column.
For (f), apply (b) to each of the \( n \) columns, and use (a) and (d) on \( \Q\T\Q=\I \).
:::

## Schur complements and block determinants

Partitioned matrices arise whenever the regressors split into two groups, the model
matrix \( \X=[\X_1,\X_2] \) giving
\[
\X\T\X=\begin{pmatrix}\X_1\T\X_1&\X_1\T\X_2\\\X_2\T\X_1&\X_2\T\X_2\end{pmatrix},
\]
and whenever a random vector splits into two blocks ([Chapter 3](../ch03-multivariate-normal/index.html)).
In what follows
\[
\A=\begin{pmatrix}\A_{11}&\A_{12}\\\A_{21}&\A_{22}\end{pmatrix},
\]{#eq-mat-partition}

with \( \A_{11} \) of size \( k\times k \) and \( \A_{22} \) of size \( (n-k)\times(n-k) \).

::: {#def-mat-schur}
[Schur complement]

If \( \A_{11} \) is nonsingular, the **Schur complement** of \( \A_{11} \) in \( \A \) is
\( \A_{22\cdot1}=\A_{22}-\A_{21}\A_{11}^{-1}\A_{12} \). If \( \A_{22} \) is nonsingular, the Schur
complement of \( \A_{22} \) is \( \A_{11\cdot2}=\A_{11}-\A_{12}\A_{22}^{-1}\A_{21} \).
:::

Everything in this subsection comes from one identity, block Gaussian elimination. If
\( \A_{11} \) is nonsingular, multiplying out shows
\[
\A=\begin{pmatrix}\I&\bzero\\\A_{21}\A_{11}^{-1}&\I\end{pmatrix}
\begin{pmatrix}\A_{11}&\bzero\\\bzero&\A_{22\cdot1}\end{pmatrix}
\begin{pmatrix}\I&\A_{11}^{-1}\A_{12}\\\bzero&\I\end{pmatrix}.
\]{#eq-mat-block-ldu}

The outer factors are block triangular with identity blocks on the diagonal, and are
nonsingular with inverses obtained by changing the sign of the off-diagonal block.

::: {#thm-mat-block-determinant}
[Block determinants]

Let \( \A \) be partitioned as in @eq-mat-partition.

::: {.enumerate options="label=(\alph*)"}
1. If \( \A_{12}=\bzero \) or \( \A_{21}=\bzero \), then \( \det\A=\det\A_{11}\det\A_{22} \).

2. If \( \A_{11} \) is nonsingular, \( \det\A=\det\A_{11}\det\A_{22\cdot1} \). If \( \A_{22} \) is
           nonsingular, \( \det\A=\det\A_{22}\det\A_{11\cdot2} \).

3. For \( m\times n \) matrices \( \bU \) and \( n\times m \) matrices \( \V \),
           \( \det(\I_m+\bU\V)=\det(\I_n+\V\bU) \). If \( \B \) is \( m\times m \) and nonsingular,
           \( \det(\B+\bU\V)=\det\B\,\det(\I_n+\V\B^{-1}\bU) \).
:::

:::

::: {.proof}
(a) Suppose \( \A_{12}=\bzero \) (otherwise transpose). In the permutation expansion, a term
that is nonzero cannot use an entry of the zero block \( \A_{12} \). So each of the first \( k \)
rows takes its entry from the first \( k \) columns, and \( \pi \) maps \( \{1,\dots,k\} \) onto
itself. Such a permutation splits into a permutation of the first \( k \) indices and one of
the rest. Its sign is the product of their signs, so the sum factors as
\( \det\A_{11}\det\A_{22} \).
(b) Take determinants in @eq-mat-block-ldu, using (a) and @prp-mat-det(d).
The outer factors have determinant \( 1 \). The second statement follows in the same way
from the elimination that starts with \( \A_{22} \).
(c) Apply (b) in both ways to
\( \begin{psmallmatrix}\I_m&-\bU\\\V&\I_n\end{psmallmatrix} \), whose Schur complements are
\( \I_n+\V\bU \) and \( \I_m+\bU\V \). For the second statement, write
\( \B+\bU\V=\B(\I_m+\B^{-1}\bU\V) \) and apply the first.
:::

Part (c) with \( n=1 \) is the **matrix determinant lemma**
\( \det(\B+\bu\bv\T)=\det\B\,(1+\bv\T\B^{-1}\bu) \).

## Inverses of partitioned matrices

::: {#thm-mat-partitioned-inverse}
[Partitioned inverse]

Let \( \A \) be partitioned as in @eq-mat-partition with \( \A_{11} \) nonsingular. Then
\( \A \) is nonsingular iff \( \A_{22\cdot1} \) is, and in that case
\[
\A^{-1}=\begin{pmatrix}
    \A_{11}^{-1}+\A_{11}^{-1}\A_{12}\A_{22\cdot1}^{-1}\A_{21}\A_{11}^{-1} &
    -\A_{11}^{-1}\A_{12}\A_{22\cdot1}^{-1}\\
    -\A_{22\cdot1}^{-1}\A_{21}\A_{11}^{-1} & \A_{22\cdot1}^{-1}
\end{pmatrix}.
\]{#eq-mat-partitioned-inverse}

If \( \A_{22} \) is also nonsingular, the leading block of \( \A^{-1} \) equals \( \A_{11\cdot2}^{-1} \):
\[
\bigl(\A_{11}-\A_{12}\A_{22}^{-1}\A_{21}\bigr)^{-1}
=\A_{11}^{-1}+\A_{11}^{-1}\A_{12}\A_{22\cdot1}^{-1}\A_{21}\A_{11}^{-1}.
\]{#eq-mat-leading-block}

:::

::: {.proof}
The nonsingularity statement is @thm-mat-block-determinant(b). Inverting the three
factors of @eq-mat-block-ldu and multiplying them in reverse order gives
\[
\A^{-1}=\begin{pmatrix}\I&-\A_{11}^{-1}\A_{12}\\\bzero&\I\end{pmatrix}
\begin{pmatrix}\A_{11}^{-1}&\bzero\\\bzero&\A_{22\cdot1}^{-1}\end{pmatrix}
\begin{pmatrix}\I&\bzero\\-\A_{21}\A_{11}^{-1}&\I\end{pmatrix},
\]
which multiplies out to @eq-mat-partitioned-inverse. Eliminating with \( \A_{22} \)
first gives a second expression for \( \A^{-1} \), whose leading block is \( \A_{11\cdot2}^{-1} \).
The inverse is unique, so the two leading blocks agree.
:::

For a symmetric positive definite \( \X\T\X \) ([Section 1.7](07-eigen.html)), the theorem says that
the block of \( (\X\T\X)^{-1} \) belonging to \( \X_1 \) is
\( \bigl(\X_1\T\X_1-\X_1\T\X_2(\X_2\T\X_2)^{-1}\X_2\T\X_1\bigr)^{-1} \). [Chapter 6](../ch06-projections/index.html)
derives the same formula geometrically, from the Frisch–Waugh–Lovell theorem.

Identity @eq-mat-leading-block is a disguised form of a result that deserves its own
name. Rename \( \A_{11}=\B \), \( \A_{12}=\bU \), \( \A_{21}=\V\T \) and \( \A_{22}=-\mathbf{C}^{-1} \).

::: {#thm-mat-woodbury}
[Sherman–Morrison–Woodbury]

Let \( \B \) (\( n\times n \)) and \( \mathbf{C} \) (\( k\times k \)) be nonsingular, and let \( \bU,\V \) be
\( n\times k \). If \( \mathbf{C}^{-1}+\V\T\B^{-1}\bU \) is nonsingular, then so is \( \B+\bU\mathbf{C}\V\T \), and
\[
(\B+\bU\mathbf{C}\V\T)^{-1}=\B^{-1}-\B^{-1}\bU\bigl(\mathbf{C}^{-1}+\V\T\B^{-1}\bU\bigr)^{-1}\V\T\B^{-1}.
\]{#eq-mat-woodbury}

In particular, for vectors \( \bu,\bv \) with \( 1+\bv\T\B^{-1}\bu\ne0 \),
\[
(\B+\bu\bv\T)^{-1}=\B^{-1}-\frac{\B^{-1}\bu\bv\T\B^{-1}}{1+\bv\T\B^{-1}\bu}.
\]{#eq-mat-sherman-morrison}

:::

::: {.proof}
Apply @thm-mat-partitioned-inverse to
\( \mathbf{N}=\begin{psmallmatrix}\B&\bU\\\V\T&-\mathbf{C}^{-1}\end{psmallmatrix} \). Both diagonal blocks
are nonsingular. The Schur complement of \( \B \) is \( -(\mathbf{C}^{-1}+\V\T\B^{-1}\bU) \), which is
nonsingular by hypothesis, so \( \mathbf{N} \) is nonsingular. Then the Schur complement of the
lower block, \( \B+\bU\mathbf{C}\V\T \), is nonsingular by the same theorem with the roles of the
blocks exchanged, and @eq-mat-leading-block becomes @eq-mat-woodbury. For
@eq-mat-sherman-morrison take \( k=1 \) and \( \mathbf{C}=1 \).
:::

::: {.idea}
Every result in this section comes from one block elimination,
@eq-mat-block-ldu: a partitioned matrix is a unit triangular matrix times
\( \diag(\A_{11},\A_{22\cdot1}) \) times another. Taking determinants gives
@thm-mat-block-determinant, inverting the factors gives
@thm-mat-partitioned-inverse, and eliminating in the other order and comparing gives
the Woodbury identity. When a formula about a block is forgotten, redo the elimination.
:::

The identity replaces the inversion of an \( n\times n \) matrix by that of a \( k\times k \)
matrix. That matters when \( \B^{-1} \) is already known or cheap, as in the following
two examples.

::: {#exm-mat-deletion}
[Deleting an observation]

Removing the \( i \)th row \( \x_i\T \) from a model matrix changes \( \X\T\X \) into
\( \X\T\X-\x_i\x_i\T \). With \( \bu=-\x_i \), \( \bv=\x_i \) and
\( h_{i}=\x_i\T(\X\T\X)^{-1}\x_i \), @eq-mat-sherman-morrison gives
\[
\bigl(\X\T\X-\x_i\x_i\T\bigr)^{-1}=(\X\T\X)^{-1}+\frac{(\X\T\X)^{-1}\x_i\x_i\T(\X\T\X)^{-1}}{1-h_{i}},
\]
valid whenever \( h_i<1 \), and the determinant lemma gives
\( \det(\X\T\X-\x_i\x_i\T)=(1-h_i)\det(\X\T\X) \). The quantity \( h_i \) is the leverage of observation \( i \)
([Chapter 6](../ch06-projections/index.html)), and these formulas are what make the deletion diagnostics of
[Chapter 20](../ch20-residuals-leverage-influence/index.html) cheap (@thm-res-deletion). For Brownlee's stack-loss data
(\( n=21 \), \( p=4 \)), the listing updates the inverse for each of
the \( 21 \) deletions. The largest relative discrepancy from direct inversion is
\( 1\times 10^{-13} \). Observation \( 17 \) has the largest leverage,
\( 0.412 \), and removing it multiplies \( \det(\X\T\X) \) by
\( 0.588 \).
:::

```{.python .run #cell-woodbury-deletion}
import numpy as np
import statsmodels.api as sm
data = sm.datasets.stackloss.load_pandas().data
X = np.column_stack([np.ones(len(data)), data[["AIRFLOW", "WATERTEMP", "ACIDCONC"]]])
n, p = X.shape

A_inv = np.linalg.inv(X.T @ X)
h = np.einsum("ij,jk,ik->i", X, A_inv, X)             # leverages x_i' (X'X)^{-1} x_i

worst = 0.0
for i in range(n):
    x = X[i]
    # Sherman-Morrison with u = -x, v = x:  (X'X - x x')^{-1}
    B_update = A_inv + np.outer(A_inv @ x, x @ A_inv) / (1 - h[i])
    X_del = np.delete(X, i, axis=0)
    B_direct = np.linalg.inv(X_del.T @ X_del)
    worst = max(worst, np.abs(B_update - B_direct).max() / np.abs(B_direct).max())
print(f"largest relative discrepancy over all {n} deletions: {worst:.1e}")
```

::: {#exm-mat-factor}
[Diagonal plus low rank]

In factor models and in the linear mixed models of Chapter 32, covariance
matrices have the form \( \bSigma=\bD+\bL\bL\T \) with \( \bD \) diagonal and positive and \( \bL \) of size
\( N\times k \), \( k\ll N \). By @eq-mat-woodbury and @thm-mat-block-determinant(c),
\[
\bSigma^{-1}=\bD^{-1}-\bD^{-1}\bL(\I_k+\bL\T\bD^{-1}\bL)^{-1}\bL\T\bD^{-1},\qquad
\det\bSigma=\det\bD\,\det(\I_k+\bL\T\bD^{-1}\bL).
\]
Each needs only a \( k\times k \) factorization. The script checks this with \( N=1500 \)
and \( k=4 \). The solution of \( \bSigma\bw=\bv \) agrees with a direct solve to a
relative error of \( 8\times 10^{-14} \).
:::

::: {#exm-mat-equicorrelation}
[Equicorrelation]

The matrix \( \mathbf{E}=(1-\rho)\I_n+\rho\mathbf{J}_n \) has unit diagonal and every off-diagonal entry
equal to \( \rho \). It is the correlation matrix of exchangeable observations. Taking
\( \B=(1-\rho)\I \) and \( \bu=\rho\bone \), \( \bv=\bone \), for \( \rho\ne1 \) and \( 1+(n-1)\rho\ne0 \),
\[
\mathbf{E}^{-1}=\frac{1}{1-\rho}\Bigl(\I_n-\frac{\rho}{1+(n-1)\rho}\mathbf{J}_n\Bigr),\qquad
\det\mathbf{E}=(1-\rho)^{n-1}\bigl(1+(n-1)\rho\bigr).
\]
:::

## Exercises

### A. Check your understanding

::: {#exr-mat-trace-order}
[A1]

Show that \( \tr(\A\B\mathbf{C})=\tr(\A\mathbf{C}\B) \) when \( \A,\B,\mathbf{C} \) are symmetric, and find \( 2\times2 \) matrices
for which the two traces differ.
:::

::: {.solution}
For symmetric matrices,
\( \tr(\A\B\mathbf{C})=\tr\bigl((\A\B\mathbf{C})\T\bigr)=\tr(\mathbf{C}\B\A)=\tr(\A\mathbf{C}\B) \), using transposition and
then the cyclic property. For a counterexample take
\( \A=\begin{psmallmatrix}0&1\\0&0\end{psmallmatrix} \),
\( \B=\begin{psmallmatrix}0&0\\1&0\end{psmallmatrix} \) and
\( \mathbf{C}=\begin{psmallmatrix}1&0\\0&0\end{psmallmatrix} \). Then \( \A\B\mathbf{C}=\mathbf{C} \) has trace \( 1 \),
while \( \A\mathbf{C}=\bzero \) gives \( \tr(\A\mathbf{C}\B)=0 \).
:::

::: {#exr-mat-bordered}
[A2]

Let \( \A \) be nonsingular, \( \mathbf{a} \) a vector and \( c \) a scalar with \( s=c-\mathbf{a}\T\A^{-1}\mathbf{a}\ne0 \). Write
down the inverse and the determinant of
\( \begin{psmallmatrix}\A&\mathbf{a}\\\mathbf{a}\T&c\end{psmallmatrix} \). Use it to update \( (\X\T\X)^{-1} \) when a
column is added to \( \X \).
:::

### B. Practice

::: {#exr-mat-leave-one-out}
[B1]

Let \( \hbeta=(\X\T\X)^{-1}\X\T\y \), \( e_i=y_i-\x_i\T\hbeta \) and \( h_i=\x_i\T(\X\T\X)^{-1}\x_i<1 \). Show that
the least squares coefficients computed without observation \( i \) are
\[
\hbeta_{(i)}=\hbeta-\frac{(\X\T\X)^{-1}\x_ie_i}{1-h_i}.
\]
:::

::: {.solution}
Write \( \mathbf{K}=(\X\T\X)^{-1} \) and \( \x=\x_i \). Deleting
row \( i \) turns \( \X\T\y \) into \( \X\T\y-\x y_i \), so by @exm-mat-deletion
\[
\hbeta_{(i)}=\Bigl(\mathbf{K}+\frac{\mathbf{K}\x\x\T\mathbf{K}}{1-h_i}\Bigr)(\X\T\y-\x y_i)
=\hbeta-\mathbf{K}\x y_i+\frac{\mathbf{K}\x\,(\x\T\hbeta-h_iy_i)}{1-h_i}
=\hbeta-\mathbf{K}\x\,\frac{y_i-\x\T\hbeta}{1-h_i},
\]
since \( (1-h_i)y_i-\x\T\hbeta+h_iy_i=y_i-\x\T\hbeta \).
:::

::: {#exr-mat-push-through}
[B2]

For an \( n\times p \) matrix \( \X \) and \( \lambda>0 \), show that
\( (\X\T\X+\lambda\I_p)^{-1}\X\T=\X\T(\X\X\T+\lambda\I_n)^{-1} \). Explain why the right-hand side is
preferable when \( p \) is much larger than \( n \), and derive the identity a second way from @thm-mat-woodbury.
:::

::: {.solution}
Both \( \X\T\X+\lambda\I \) and \( \X\X\T+\lambda\I \) are
positive definite. The identity
\( (\X\T\X+\lambda\I_p)\X\T=\X\T\X\X\T+\lambda\X\T=\X\T(\X\X\T+\lambda\I_n) \) gives the result after
multiplying by the two inverses. When \( p\gg n \), the right-hand side inverts an \( n\times n \) matrix
instead of a \( p\times p \) one.
For the second derivation, apply @thm-mat-woodbury with \( \B=\lambda\I_p \),
\( \bU=\X\T \), \( \mathbf{C}=\I_n \) and \( \V=\X\T \):
\( (\lambda\I_p+\X\T\X)^{-1}=\lambda^{-1}\I_p-\lambda^{-2}\X\T\mathbf{W}^{-1}\X \) with
\( \mathbf{W}=\I_n+\X\X\T/\lambda \). Multiplying on the right by \( \X\T \) gives
\( \lambda^{-1}\X\T\bigl[\I_n-\mathbf{W}^{-1}(\X\X\T/\lambda)\bigr]
=\lambda^{-1}\X\T\mathbf{W}^{-1}\bigl[\mathbf{W}-\X\X\T/\lambda\bigr]=\lambda^{-1}\X\T\mathbf{W}^{-1}
=\X\T(\X\X\T+\lambda\I_n)^{-1} \).
:::
