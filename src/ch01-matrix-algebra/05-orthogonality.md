# Orthogonality and orthogonal matrices

## Orthogonal vectors and complements

Vectors \( \bu,\bv \) are **orthogonal**, written \( \bu\perp\bv \), if \( \bu\T\bv=0 \). For
orthogonal vectors the Pythagorean identity \( \norm{\bu+\bv}^2=\norm{\bu}^2+\norm{\bv}^2 \) holds.
A set of nonzero, mutually orthogonal vectors is independent: if \( \sum_jc_j\bv_j=\bzero \),
taking the inner product with \( \bv_k \) leaves \( c_k\norm{\bv_k}^2=0 \). A set is
**orthonormal** if, in addition, each vector has unit length. The columns of an
\( n\times r \) matrix \( \Q \) are orthonormal iff \( \Q\T\Q=\I_r \).

The **orthogonal complement** of a subspace \( \mathcal S\subseteq\Real^n \) is
\( \mathcal S\perpc=\{\bv:\bv\T\bu=0\text{ for all }\bu\in\mathcal S\} \). The next
proposition is all that this chapter needs about complements. [Chapter 6](../ch06-projections/index.html)
builds the theory of orthogonal projections on it.

::: {#prp-mat-complement}
Let \( \A \) be \( m\times n \) and let \( \mathcal S\subseteq\Real^n \) be a subspace.

::: {.enumerate options="label=(\alph*)"}
1. \( \Null(\A)=\C(\A\T)\perpc \) and \( \Null(\A\T)=\C(\A)\perpc \).

2. \( \dim\mathcal S\perpc=n-\dim\mathcal S \), and \( \Real^n=\mathcal S\dirsum\mathcal S\perpc \).
           In particular \( \Real^n=\C(\A\T)\dirsum\Null(\A) \) and \( \Real^m=\C(\A)\dirsum\Null(\A\T) \).

3. \( (\mathcal S\perpc)\perpc=\mathcal S \).
:::

:::

::: {.proof}
(a) \( \A\bv=\bzero \) says that \( \bv \) is orthogonal to every row of \( \A \), hence to their span.
Apply this to \( \A\T \) for the second identity.
(b) Put a basis of \( \mathcal S \) (\( k \) vectors) into the columns of \( \B \). Then
\( \mathcal S\perpc=\Null(\B\T) \), which has dimension \( n-k \) by rank–nullity. If
\( \bu\in\mathcal S\cap\mathcal S\perpc \) then \( \bu\T\bu=0 \), so the sum is direct and has
dimension \( n \), so it is \( \Real^n \).
(c) \( \mathcal S\subseteq(\mathcal S\perpc)\perpc \), and both have dimension \( \dim\mathcal S \) by (b).
:::

## Gram–Schmidt and the QR factorization

::: {#thm-mat-qr}
[Gram–Schmidt; QR factorization]

Let \( \X=[\x_1,\dots,\x_p] \) be \( n\times p \) with full column rank. Define recursively
\[
\bw_j=\x_j-\sum_{k<j}(\bu_k\T\x_j)\,\bu_k,\qquad \bu_j=\bw_j/\norm{\bw_j}.
\]
Then \( \bu_1,\dots,\bu_p \) are orthonormal and, for each \( j \), span the same subspace as
\( \x_1,\dots,\x_j \). Equivalently, \( \X=\Q\R \) with \( \Q=[\bu_1,\dots,\bu_p] \) satisfying
\( \Q\T\Q=\I_p \) and \( \R \) upper triangular with positive diagonal entries.
:::

::: {.proof}
By induction on \( j \). Suppose \( \bu_1,\dots,\bu_{j-1} \) are orthonormal and span
\( \spn(\x_1,\dots,\x_{j-1}) \). Then \( \bw_j\ne\bzero \), since otherwise \( \x_j \) would lie in that
span, contradicting independence. For \( i<j \),
\( \bu_i\T\bw_j=\bu_i\T\x_j-\bu_i\T\x_j=0 \). So \( \bu_1,\dots,\bu_j \) are orthonormal, and they
span \( \spn(\x_1,\dots,\x_j) \) because \( \x_j=\norm{\bw_j}\bu_j+\sum_{k<j}(\bu_k\T\x_j)\bu_k \). That
last identity, read column by column, is \( \X=\Q\R \) with \( r_{jj}=\norm{\bw_j}>0 \) and
\( r_{kj}=\bu_k\T\x_j \) for \( k<j \).
:::

The factorization is unique (@exr-mat-qr-unique). Numerical software computes it
by a more stable route, but the object is the same, and [Chapter 6](../ch06-projections/index.html) and Chapter 10
rely on it. Two consequences are used repeatedly. First, any orthonormal set
\( \bu_1,\dots,\bu_r \) in \( \Real^n \) extends to an orthonormal basis: apply Gram–Schmidt to
\( \bu_1,\dots,\bu_r,\bm e_1,\dots,\bm e_n \), discarding each vector that depends on its
predecessors. Second, every subspace has an orthonormal basis.

## Orthogonal matrices

A square matrix \( \Q \) is **orthogonal** if \( \Q\T\Q=\I \). By @prp-mat-inverse,
\( \Q\T \) is then the inverse of \( \Q \), so also \( \Q\Q\T=\I \): the rows of an orthogonal matrix are
orthonormal as well as its columns.

::: {#prp-mat-orthogonal}
Let \( \Q \) and \( \bm P \) be \( n\times n \) orthogonal matrices.

::: {.enumerate options="label=(\alph*)"}
1. \( (\Q\x)\T(\Q\bv)=\x\T\bv \) and \( \norm{\Q\x}=\norm{\x} \) for all \( \x,\bv \).

2. \( \Q\T \) and \( \Q\bm P \) are orthogonal.

3. Every entry satisfies \( \lvert q_{ij}\rvert\le1 \), and \( \det\Q=\pm1 \).

4. If \( [\Q_1,\Q_2] \) is orthogonal with \( \Q_1 \) of size \( n\times r \), then
           \( \Q_1\Q_1\T+\Q_2\Q_2\T=\I_n \), and \( \C(\Q_2)=\C(\Q_1)\perpc \).
:::

:::

::: {.proof}
(a) \( (\Q\x)\T\Q\bv=\x\T\Q\T\Q\bv=\x\T\bv \). (b) \( (\Q\bm P)\T\Q\bm P=\bm P\T\bm P=\I \).
(c) Each column has unit length. For the determinant, see
@prp-mat-det: \( 1=\det(\Q\T\Q)=(\det\Q)^2 \). (d) The first statement is
\( \Q\Q\T=\I \) written in blocks. For the second, \( \Q_1\T\Q_2=\bzero \) gives
\( \C(\Q_2)\subseteq\C(\Q_1)\perpc \), and both have dimension \( n-r \).
:::

Part (a) says an orthogonal matrix is a rigid motion that fixes the origin: a rotation, a
reflection, or a composition of the two. Changing coordinates by an orthogonal
matrix changes neither lengths nor angles. This is why the distribution theory of
[Chapter 4](../ch04-quadratic-forms/index.html) can rotate a normal vector with independent components into a
convenient coordinate system without changing its distribution.

::: {#exm-mat-helmert}
[Helmert matrices]

The \( n\times n \) **Helmert matrix** \( \bm H \) has first row \( n^{-1/2}\bone\T \) and, for
\( k=2,\dots,n \), \( k \)th row
\[
\bm h_k\T=\frac{1}{\sqrt{k(k-1)}}\bigl(\underbrace{1,\dots,1}_{k-1},\,-(k-1),\,0,\dots,0\bigr).
\]
Each row has unit length, since \( (k-1)+(k-1)^2=k(k-1) \). Row \( k \) is orthogonal to row \( 1 \)
because its entries sum to zero. For \( 2\le j<k \), the nonzero entries of \( \bm h_j \) sit in
positions \( 1,\dots,j \), where \( \bm h_k \) is constant, so \( \bm h_j\T\bm h_k \) is proportional to
the sum of the entries of \( \bm h_j \), which is zero. So \( \bm H \) is orthogonal. For
\( \y\in\Real^n \), \( (\bm H\y)_1=\sqrt n\,\bar{y} \) and, by @prp-mat-orthogonal(a),
\[
\sum_{k=2}^n(\bm h_k\T\y)^2=\norm{\bm H\y}^2-n\bar{y}^2=\sum_{i=1}^n(y_i-\bar{y})^2 .
\]
The sum of squared deviations is a sum of \( n-1 \) squares of orthonormal contrasts. With
independent normal observations this gives a short proof that the sample mean and variance
are independent (@exr-mvn-rotation).
:::

Two other families appear later. **Permutation matrices**, the identity with its rows
reordered, are orthogonal. For a unit vector \( \bu \), the **Householder reflection**
\( \I-2\bu\bu\T \) is symmetric and orthogonal and reflects \( \Real^n \) in the hyperplane
\( \bu\perpc \). Householder reflections are how QR factorizations are computed in practice
(@exr-mat-householder).

## Exercises

### A. Check your understanding

::: {#exr-mat-rotations}
[A1]

Show that every \( 2\times2 \) orthogonal matrix is either a rotation
\( \begin{psmallmatrix}\cos\theta&-\sin\theta\\\sin\theta&\cos\theta\end{psmallmatrix} \)
(determinant \( 1 \)) or a reflection
\( \begin{psmallmatrix}\cos\theta&\sin\theta\\\sin\theta&-\cos\theta\end{psmallmatrix} \)
(determinant \( -1 \)). Find the line fixed by the reflection.
:::

### B. Practice

::: {#exr-mat-householder}
[B1]

Let \( \bu \) be a unit vector and \( \bm H=\I-2\bu\bu\T \). Show that \( \bm H \) is symmetric and orthogonal,
that \( \bm H\bu=-\bu \), and that \( \bm H\bv=\bv \) for \( \bv\perp\bu \). Given \( \x \) not a nonnegative multiple
of \( \bm e_1 \), find \( \bu \) with \( \bm H\x=\norm{\x}\bm e_1 \). Explain how repeating this column by
column produces a QR factorization.
:::

::: {.solution}
\( \bm H\T=\bm H \), and
\( \bm H\T\bm H=\I-4\bu\bu\T+4\bu(\bu\T\bu)\bu\T=\I \). Also \( \bm H\bu=\bu-2\bu=-\bu \), and \( \bm H\bv=\bv \)
when \( \bu\T\bv=0 \). Let \( \bw=\x-\norm{\x}\bm e_1\ne\bzero \) and \( \bu=\bw/\norm{\bw} \). Then
\( \bw\T\bw=2\norm{\x}^2-2\norm{\x}x_1=2\bw\T\x \), so
\( \bm H\x=\x-2\bw(\bw\T\x)/(\bw\T\bw)=\x-\bw=\norm{\x}\bm e_1 \). Applying such a reflection \( \bm H_1 \) to
the first column of \( \X \), then a reflection acting on coordinates \( 2,\dots,n \) to the second
column of \( \bm H_1\X \), and so on, gives \( \bm H_p\cdots\bm H_1\X=\begin{psmallmatrix}\R\\\bzero\end{psmallmatrix} \)
with \( \R \) upper triangular. Hence \( \X=\Q\R \), where \( \Q \) holds the first \( p \) columns of the orthogonal
matrix \( \bm H_1\cdots\bm H_p \). In floating point one uses \( \bw=\x+\operatorname{sgn}(x_1)\norm{\x}\bm e_1 \)
instead, which avoids cancellation. It gives
\( \bm H\x=-\operatorname{sgn}(x_1)\norm{\x}\bm e_1 \), so diagonal entries of \( \R \) may come out
negative; changing the sign of such a row of \( \R \) and of the matching column of \( \Q \) at the end
restores the positive diagonal of @thm-mat-qr.
:::

::: {#exr-mat-qr-unique}
[B2]

Show that the QR factorization of a matrix with full column rank, with \( \R \) having positive
diagonal, is unique, and that \( \R\T \) is the Cholesky factor of \( \X\T\X \).
:::

::: {.solution}
If \( \X=\Q_1\R_1=\Q_2\R_2 \), then
\( \X\T\X=\R_1\T\R_1=\R_2\T\R_2 \). Each \( \R_k\T \) is lower triangular with positive diagonal, so both
are the Cholesky factor of \( \X\T\X \), which is unique by
@thm-mat-pd-characterizations(d). So \( \R_1=\R_2 \) and \( \Q_1=\X\R_1^{-1}=\Q_2 \).
:::
