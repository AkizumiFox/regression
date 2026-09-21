# Rank, column space and null space

## Subspaces and dimension

A nonempty set \( \mathcal S\subseteq\Real^n \) is a **subspace** if \( \bu,\bv\in\mathcal S \)
and \( a,b\in\Real \) imply \( a\bu+b\bv\in\mathcal S \). Vectors \( \bv_1,\dots,\bv_k \) are
**linearly independent** if \( \sum_jc_j\bv_j=\bzero \) only for \( c_1=\dots=c_k=0 \), that is,
if the matrix \( [\bv_1,\dots,\bv_k] \) sends only the zero vector to zero. Their
**span** is the set of all their linear combinations. A **basis** of \( \mathcal S \) is
an independent set that spans it.

We use the following facts of elementary linear algebra without proof
(Halmos 1958, chapter 1). Every subspace of \( \Real^n \) has a basis. All bases of
a subspace have the same number of elements, its **dimension** \( \dim\mathcal S \). In an
\( r \)-dimensional subspace, any \( r \) independent vectors form a basis, and any independent set
can be extended to a basis. If \( \mathcal S\subseteq\mathcal T \) and
\( \dim\mathcal S=\dim\mathcal T \), then \( \mathcal S=\mathcal T \). Finally, if
\( \mathcal S\cap\mathcal T=\{\bzero\} \), the union of a basis of \( \mathcal S \) and a basis of
\( \mathcal T \) is independent, so the sum \( \mathcal S+\mathcal T=\{\bu+\bv:\bu\in\mathcal
S,\bv\in\mathcal T\} \) has dimension \( \dim\mathcal S+\dim\mathcal T \). Such a sum is
called **direct** and written \( \mathcal S\dirsum\mathcal T \).

## The four subspaces of a matrix

An \( m\times n \) matrix \( \A \) carries two subspaces of \( \Real^m \) and two of \( \Real^n \):
\[
\begin{aligned}  \C(\A)&=\{\A\bb:\bb\in\Real^n\}\subseteq\Real^m, &
\Null(\A)&=\{\bb\in\Real^n:\A\bb=\bzero\}\subseteq\Real^n,\\
\C(\A\T)&\subseteq\Real^n, & \Null(\A\T)&\subseteq\Real^m .\end{aligned}
\]
By @eq-mat-product-readings, \( \C(\A) \) is the span of the columns, the
**column space**, and \( \C(\A\T) \) is the span of the rows, the **row space**.
\( \Null(\A) \) is the **null space**. The **rank** of \( \A \) is \( \rank(\A)=\dim\C(\A) \),
the largest number of independent columns. The number of independent rows,
\( \dim\C(\A\T) \), looks like a different quantity. It is not.

::: {#thm-mat-row-col-rank}
[Row rank equals column rank]

For every matrix \( \A \), \( \dim\C(\A\T)=\dim\C(\A) \). That is, \( \rank(\A\T)=\rank(\A) \).
Moreover, if \( \rank(\A)=r\ge1 \) there are an \( m\times r \) matrix \( \B \) and an \( r\times n \)
matrix \( \mathbf{F} \), both of rank \( r \), with \( \A=\B\mathbf{F} \).
:::

::: {.proof}
Let \( r=\rank(\A) \) and let the columns of the \( m\times r \) matrix \( \B \) be a basis of \( \C(\A) \).
Each column \( \mathbf{a}_j \) is a combination \( \B\mathbf{f}_j \) of these, so \( \A=\B\mathbf{F} \) with
\( \mathbf{F}=[\mathbf{f}_1,\dots,\mathbf{f}_n] \). By the transposed reading of
@eq-mat-product-readings, every row of \( \A=\B\mathbf{F} \) is a combination of the \( r \) rows
of \( \mathbf{F} \), so \( \dim\C(\A\T)\le r=\dim\C(\A) \). Applying this to \( \A\T \) gives the reverse
inequality. Finally \( \rank(\mathbf{F})\le r \) because \( \mathbf{F} \) has \( r \) rows, and
\( r=\dim\C(\A\T)\le\dim\C(\mathbf{F}\T)=\rank(\mathbf{F}) \), so \( \rank(\mathbf{F})=r \).
:::

The factorization \( \A=\B\mathbf{F} \) is a **rank factorization**: \( \B \) has
**full column rank** (its columns are independent) and \( \mathbf{F} \) has **full row
rank**. It is a convenient tool for proofs, as @prp-mat-idempotent-basic shows.

::: {#thm-mat-rank-nullity}
[Rank–nullity]

For an \( m\times n \) matrix \( \A \),  \( \rank(\A)+\dim\Null(\A)=n \).
:::

::: {.proof}
Let \( \bv_1,\dots,\bv_k \) be a basis of \( \Null(\A) \) and extend it to a basis
\( \bv_1,\dots,\bv_k,\bw_1,\dots,\bw_{n-k} \) of \( \Real^n \). Every \( \A\bb \) is a combination of
\( \A\bw_1,\dots,\A\bw_{n-k} \), since \( \A\bv_i=\bzero \). These vectors are independent: if
\( \sum_jc_j\A\bw_j=\bzero \), then \( \sum_jc_j\bw_j\in\Null(\A) \) is a combination of the
\( \bv_i \), and independence of the full basis forces every \( c_j=0 \). So
\( \rank(\A)=n-k \).
:::

In a linear model with \( n\times p \) model matrix \( \X \), \( \dim\Null(\X)=p-\rank(\X) \) counts the
directions in which the coefficient vector can move without changing the mean
\( \X\bbeta \). \( \X \) has full column rank iff \( \Null(\X)=\{\bzero\} \).

## Rank of products and of \( \A\T\A \)

::: {#prp-mat-rank-product}
Let \( \A \) be \( m\times n \) and \( \B \) be \( n\times q \).

::: {.enumerate options="label=(\alph*)"}
1. \( \C(\A\B)\subseteq\C(\A) \) and \( \Null(\B)\subseteq\Null(\A\B) \), so
           \( \rank(\A\B)\le\min\{\rank(\A),\rank(\B)\} \).

2. If \( \B \) has full row rank, then \( \C(\A\B)=\C(\A) \). If \( \A \) has full column rank,
           then \( \Null(\A\B)=\Null(\B) \). So \( \rank(\A\B)=\rank(\A) \) in the first case and
           \( \rank(\A\B)=\rank(\B) \) in the second. In particular, multiplying by a nonsingular
           matrix does not change rank.

3. \( \Null(\A\T\A)=\Null(\A) \) and \( \C(\A\T\A)=\C(\A\T) \), so
           \( \rank(\A\T\A)=\rank(\A\A\T)=\rank(\A) \).

4. \( \rank(\A+\mathbf{C})\le\rank(\A)+\rank(\mathbf{C}) \) for \( \mathbf{C} \) of the same size as \( \A \).
:::

:::

::: {.proof}
(a) is @eq-mat-product-readings, together with \( \rank(\A\B)=\rank(\B\T\A\T)\le\rank(\B\T) \)
by the first inclusion and @thm-mat-row-col-rank. (b) If \( \B \) has full row rank,
\( \C(\B)=\Real^n \) because it is an \( n \)-dimensional subspace of \( \Real^n \), so every \( \A\bu \)
equals \( \A\B\bv \) for some \( \bv \). If \( \A \) has full column rank, \( \A\B\bv=\bzero \) forces
\( \B\bv=\bzero \). The rank statements follow from rank–nullity. (c) If \( \A\T\A\bv=\bzero \)
then \( \norm{\A\bv}^2=\bv\T\A\T\A\bv=0 \), so \( \A\bv=\bzero \). Thus the null spaces agree, and
by @thm-mat-rank-nullity \( \rank(\A\T\A)=\rank(\A) \). Since
\( \C(\A\T\A)\subseteq\C(\A\T) \) and both have dimension \( \rank(\A) \), they are equal.
Applying this to \( \A\T \) gives \( \rank(\A\A\T)=\rank(\A) \). (d) \( \C(\A+\mathbf{C})\subseteq\C(\A)+\C(\mathbf{C}) \),
and the dimension of a sum is at most the sum of the dimensions.
:::

Part (c) is the reason the normal equations \( \X\T\X\bb=\X\T\y \) always have a solution: the
right-hand side lies in \( \C(\X\T)=\C(\X\T\X) \) ([Section 1.4](04-systems.html)).

## Nonsingular submatrices

::: {#prp-mat-nonsingular-submatrix}
Let \( \A \) have rank \( r\ge1 \). Choose any \( r \) linearly independent rows, indexed by \( I \), and any
\( r \) linearly independent columns, indexed by \( J \). Then the \( r\times r \) submatrix
\( \A_{IJ} \) at their intersection is nonsingular, and
\[
\A=\A_{\cdot J}\,\A_{IJ}^{-1}\,\A_{I\cdot},
\]{#eq-mat-skeleton}

where \( \A_{\cdot J} \) holds the columns in \( J \) and \( \A_{I\cdot} \) the rows in \( I \).
:::

::: {.proof}
The columns in \( J \) are a basis of \( \C(\A) \), so \( \A=\A_{\cdot J}\mathbf{F} \) for some \( r\times n \)
matrix \( \mathbf{F} \). Keeping only the rows in \( I \) gives \( \A_{I\cdot}=\A_{IJ}\mathbf{F} \). The left side
has rank \( r \), so by @prp-mat-rank-product(a) the \( r\times r \) matrix \( \A_{IJ} \) has rank
\( r \) and is nonsingular. Then \( \mathbf{F}=\A_{IJ}^{-1}\A_{I\cdot} \), which gives @eq-mat-skeleton.
:::

Conversely, a nonsingular \( k\times k \) submatrix makes the \( k \) columns of \( \A \) that contain it
linearly independent, so \( \rank(\A)\ge k \). In particular, a matrix has rank \( r \) iff it has a
nonsingular \( r\times r \) submatrix and no nonsingular larger one. This proposition gives an explicit generalized inverse in @thm-mat-ginverse-exists.

## Exercises

### A. Check your understanding

::: {#exr-mat-pairwise}
[A1]

Show that \( \sum_{i<j}(y_i-y_j)^2=n\sum_i(y_i-\bar{y})^2 \), and write both sides as
\( \y\T\A\y \) with \( \A \) symmetric. Find the symmetric matrix of the quadratic form
\( 2x_1^2-x_2^2+6x_1x_2-4x_2x_3 \).
:::

::: {#exr-mat-cancellation}
[A2]

Show that \( \X\T\X\B=\X\T\X\mathbf{C} \) implies \( \X\B=\X\mathbf{C} \). Deduce that \( \X\T\X\bb=\bzero \) iff
\( \X\bb=\bzero \).
:::

::: {.solution}
Put \( \mathbf{D}=\B-\mathbf{C} \), so that \( \X\T\X\mathbf{D}=\bzero \). Then
\( (\X\mathbf{D})\T(\X\mathbf{D})=\mathbf{D}\T\X\T\X\mathbf{D}=\bzero \), and by @eq-mat-frobenius a matrix
\( \mathbf{Z} \) with \( \mathbf{Z}\T\mathbf{Z}=\bzero \) is zero. So \( \X\mathbf{D}=\bzero \). The second statement is the case
\( \mathbf{C}=\bzero \), with the converse obvious.
:::

::: {#exr-mat-two-way-rank}
[A3]

For the additive two-way layout with two rows, three columns and one observation per cell, write
the \( 6\times6 \) model matrix \( \X=[\bone,\mathbf{R},\mathbf{C}] \), where \( \mathbf{R} \) and \( \mathbf{C} \) hold the row and
column indicators. Find \( \rank(\X) \) and a basis of \( \Null(\X) \). Compare \( [\mathbf{R},\mathbf{C}]\T \) with the
matrix \( \A \) of @exm-mat-margins.
:::

### B. Practice

::: {#exr-mat-sylvester}
[B1]

Let \( \A \) be \( m\times n \) and \( \B \) be \( n\times q \). Prove
\( \rank(\A\B)\ge\rank(\A)+\rank(\B)-n \). *Hint:* apply rank–nullity to \( \A\B_0 \), where the
columns of \( \B_0 \) are a basis of \( \C(\B) \).
:::

::: {.solution}
Let the \( n\times s \) matrix \( \B_0 \), \( s=\rank(\B) \), have columns
forming a basis of \( \C(\B) \). Then \( \C(\A\B)=\C(\A\B_0) \), so by rank–nullity
\( \rank(\A\B)=s-\dim\Null(\A\B_0) \). The map \( \mathbf{c}\mapsto\B_0\mathbf{c} \) is one-to-one and sends
\( \Null(\A\B_0) \) into \( \Null(\A) \), so \( \dim\Null(\A\B_0)\le\dim\Null(\A)=n-\rank(\A) \). Hence
\( \rank(\A\B)\ge\rank(\B)-n+\rank(\A) \).
:::

::: {#exr-mat-rank-factorization}
[B2]

Show that if \( \A=\B_1\mathbf{F}_1=\B_2\mathbf{F}_2 \) are two rank factorizations, then \( \B_2=\B_1\bT \) and
\( \mathbf{F}_2=\bT^{-1}\mathbf{F}_1 \) for a nonsingular \( \bT \).
:::

### C. Going deeper

::: {#exr-mat-frobenius-rank}
[C1]

First prove a lemma: if \( \mathcal S\subseteq\Real^n \) is a subspace and \( \A \) is \( m\times n \),
then the image \( \A\mathcal S=\{\A\x:\x\in\mathcal S\} \) has
\( \dim\A\mathcal S=\dim\mathcal S-\dim(\mathcal S\cap\Null(\A)) \). Use it to prove Frobenius's
inequality
\[
\rank(\A\B)+\rank(\B\mathbf{C})\le\rank(\B)+\rank(\A\B\mathbf{C})
\]
for \( \A \) of size \( m\times n \), \( \B \) of size \( n\times q \) and \( \mathbf{C} \) of size \( q\times s \),
and recover @exr-mat-sylvester as the case \( \B=\I_n \).
:::

::: {.solution}
For the lemma, let the columns of \( \bS_0 \) be a basis of \( \mathcal S \), so that
\( \A\mathcal S=\C(\A\bS_0) \) and \( \A\bS_0 \) has \( \dim\mathcal S \) columns. By
@thm-mat-rank-nullity, \( \dim\C(\A\bS_0)=\dim\mathcal S-\dim\Null(\A\bS_0) \). The map
\( \mathbf{c}\mapsto\bS_0\mathbf{c} \) is one-to-one, because the columns of \( \bS_0 \) are independent, and
it carries \( \Null(\A\bS_0) \) onto \( \mathcal S\cap\Null(\A) \). So the two have the same
dimension, which proves the lemma.

Apply it twice, to \( \mathcal S=\C(\B) \) and to \( \mathcal S=\C(\B\mathbf{C}) \):
\[
\rank(\A\B)=\rank(\B)-\dim\bigl(\C(\B)\cap\Null(\A)\bigr),\qquad
\rank(\A\B\mathbf{C})=\rank(\B\mathbf{C})-\dim\bigl(\C(\B\mathbf{C})\cap\Null(\A)\bigr).
\]
Since \( \C(\B\mathbf{C})\subseteq\C(\B) \), the second intersection is contained in the first, so its
dimension is no larger. Hence
\( \rank(\A\B\mathbf{C})\ge\rank(\B\mathbf{C})-\dim(\C(\B)\cap\Null(\A))=\rank(\B\mathbf{C})-\rank(\B)+\rank(\A\B) \),
which is the inequality. With \( \B=\I_n \) it reads
\( \rank(\A)+\rank(\mathbf{C})\le n+\rank(\A\mathbf{C}) \).
:::

::: {#exr-mat-dim-sum}
[C2]

Prove that any two subspaces \( \mathcal S,\mathcal T\subseteq\Real^n \) satisfy
\( \dim(\mathcal S+\mathcal T)=\dim\mathcal S+\dim\mathcal T-\dim(\mathcal S\cap\mathcal T) \),
which the text states only when the intersection is \( \{\bzero\} \). Deduce that
\[
\rank[\X_1,\X_2]=\rank(\X_1)+\rank(\X_2)-\dim\bigl(\C(\X_1)\cap\C(\X_2)\bigr),
\]
and use it to recompute the rank in @exr-mat-two-way-rank.
:::

::: {.solution}
Let \( \bu_1,\dots,\bu_d \) be a basis of \( \mathcal S\cap\mathcal T \), extended to a basis
\( \bu_1,\dots,\bu_d,\bv_1,\dots,\bv_{k-d} \) of \( \mathcal S \) (\( k=\dim\mathcal S \)) and to a
basis \( \bu_1,\dots,\bu_d,\bw_1,\dots,\bw_{l-d} \) of \( \mathcal T \) (\( l=\dim\mathcal T \)). The
\( k+l-d \) vectors \( \bu_i,\bv_j,\bw_r \) clearly span \( \mathcal S+\mathcal T \). They are
independent: suppose \( \sum_ia_i\bu_i+\sum_jb_j\bv_j+\sum_rc_r\bw_r=\bzero \). Then
\( \sum_rc_r\bw_r=-\sum_ia_i\bu_i-\sum_jb_j\bv_j \) lies in \( \mathcal S \) and in \( \mathcal T \),
hence in \( \mathcal S\cap\mathcal T \), so it is a combination of the \( \bu_i \). As
\( \bu_1,\dots,\bu_d,\bw_1,\dots,\bw_{l-d} \) are independent, every \( c_r=0 \); and then every
\( a_i \) and \( b_j \) vanishes as well. So the dimension of the sum is \( k+l-d \).

Since \( \C[\X_1,\X_2]=\C(\X_1)+\C(\X_2) \), the rank formula follows. In the two-way layout with
two rows and three columns, \( \C(\mathbf{R}) \) consists of the vectors that are constant within a
row and \( \C(\mathbf{C}) \) of those constant within a column, of dimensions \( 2 \) and \( 3 \); a vector
in both is constant, so the intersection is \( \spn(\bone) \), of dimension one. Hence
\( \rank[\mathbf{R},\mathbf{C}]=2+3-1=4 \), and adjoining \( \bone\in\C(\mathbf{R}) \) does not change it, so
\( \rank(\X)=4 \).
:::

::: {#exr-mat-rank-square}
[C3]

Show that \( \rank(\A^2)=\rank(\A) \) for a square \( \A \) iff \( \C(\A)\cap\Null(\A)=\{\bzero\} \), and
that this in turn holds iff \( \Real^n=\C(\A)\dirsum\Null(\A) \). Give a matrix with
\( \rank(\A^2)<\rank(\A) \). Show that every symmetric matrix, and every idempotent matrix,
satisfies the condition.
:::

::: {.solution}
The lemma of @exr-mat-frobenius-rank with \( \mathcal S=\C(\A) \) gives
\( \rank(\A^2)=\rank(\A)-\dim(\C(\A)\cap\Null(\A)) \), which proves the first equivalence. By
@thm-mat-rank-nullity the two subspaces have dimensions adding to \( n \), so their sum is all
of \( \Real^n \) exactly when it is direct, that is, exactly when they meet only at \( \bzero \).
For \( \A=\begin{psmallmatrix}0&1\\0&0\end{psmallmatrix} \), \( \rank(\A)=1 \) but \( \A^2=\bzero \);
here \( \C(\A)=\Null(\A)=\spn(\mathbf{e}_1) \). If \( \A \) is symmetric and \( \x=\A\mathbf{y} \) lies in
\( \Null(\A) \), then \( \x\T\x=\mathbf{y}\T\A\T\x=\mathbf{y}\T\A\x=0 \), so \( \x=\bzero \). If \( \A \) is
idempotent, then \( \A^2=\A \) and the ranks are equal by definition.
:::
