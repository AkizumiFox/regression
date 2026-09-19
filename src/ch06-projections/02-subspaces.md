# The subspaces we need

This section gathers the vector-space facts used in the rest of the
chapter. Readers who have worked through [Chapter 1](../ch01-matrix-algebra/index.html) can skim the
definitions. The results on orthogonal complements, especially
@lem-proj-null-colspace and @thm-proj-direct-sum, are used constantly
afterwards and are proved in full.

## Subspaces, spans and dimension

A nonempty set \( \mathcal S\subseteq\Real^n \) is a **subspace** if it is closed
under addition and scalar multiplication: \( \bu,\bv\in\mathcal S \) and
\( a,b\in\Real \) imply \( a\bu+b\bv\in\mathcal S \). Every subspace contains \( \bzero \).
The **span** of vectors \( \bv_1,\dots,\bv_k \) is the set of all their linear
combinations. It is the smallest subspace containing them. A set of vectors is
*linearly independent* if the only linear combination equal to \( \bzero \) is
the one with all coefficients zero. A **basis** for \( \mathcal S \) is a linearly
independent set that spans \( \mathcal S \). All bases of \( \mathcal S \) have the same
number of elements, the **dimension** \( \dim\mathcal S \).

For an \( n\times p \) matrix \( \A \) there are two subspaces we use constantly:
\[
\begin{aligned}  \C(\A) &= \{\A\bb:\bb\in\Real^p\}\subseteq\Real^n &&\text{(column space)},\\
\Null(\A) &= \{\bb\in\Real^p:\A\bb=\bzero\}\subseteq\Real^p &&\text{(null space)}.\end{aligned}
\]
The **rank** of \( \A \) is \( \rank(\A)=\dim\C(\A) \). Two standard facts from
[Chapter 1](../ch01-matrix-algebra/index.html): row rank equals column rank, so
\( \rank(\A)=\rank(\A\T) \); and the rank–nullity theorem
\[
\rank(\A)+\dim\Null(\A)=p .
\]{#eq-proj-rank-nullity}

In regression, \( \C(\X) \) is the **model space**. \( \Null(\X) \) is the set of
directions in which \( \bbeta \) can move without changing \( \X\bbeta \). It is
\( \{\bzero\} \) exactly when \( \X \) has full column rank.

::: {#exm-proj-oneway-rank}
[A rank-deficient layout]

Suppose \( n \) observations fall into \( g \) groups and
\( \X=[\bone,\,\bz_1,\dots,\bz_g] \), where \( \bz_k \) is the indicator of group \( k \).
Since \( \bz_1+\dots+\bz_g=\bone \), the vector \( (1,-1,\dots,-1)\T \) lies in
\( \Null(\X) \). The \( g \) indicators are linearly independent (their supports are
disjoint), so \( \rank(\X)=g \) while \( p=g+1 \). By
@eq-proj-rank-nullity, \( \Null(\X) \) is exactly the line spanned by
\( (1,-1,\dots,-1)\T \). The model space \( \C(\X)=\C([\bz_1,\dots,\bz_g]) \) is the set
of vectors that are constant within groups.
:::

## Inner products, orthogonality and complements

Unless stated otherwise \( \Real^n \) carries the Euclidean inner product
\( \inner{\bu}{\bv}=\bu\T\bv \) and norm \( \norm{\bu}=(\bu\T\bu)^{1/2} \). Vectors
are **orthogonal**, written \( \bu\perp\bv \), if \( \bu\T\bv=0 \). The
Cauchy–Schwarz inequality \( \lvert\bu\T\bv\rvert\le\norm{\bu}\norm{\bv} \), with
equality iff one vector is a multiple of the other, lets us define the angle
\( \theta\in[0,\pi] \) between nonzero vectors by
\[
\cos\theta=\frac{\bu\T\bv}{\norm{\bu}\,\norm{\bv}} .
\]{#eq-proj-cosine}

The angle comes back in [Section 6.7](07-reparameterization.html), where \( R^2 \) turns out
to be a squared cosine. The Pythagorean identity
\[
\bu\perp\bv \quad\Longrightarrow\quad \norm{\bu+\bv}^2=\norm{\bu}^2+\norm{\bv}^2
\]{#eq-proj-pythagoras}

follows by expanding \( (\bu+\bv)\T(\bu+\bv) \). The cross terms are the ones that vanish.

::: {#def-proj-complement}
[Orthogonal complement]

The **orthogonal complement** of a subspace \( \mathcal S\subseteq\Real^n \) is
\[
\mathcal S\perpc=\{\bv\in\Real^n:\bv\T\bu=0\text{ for all }\bu\in\mathcal S\}.
\]
:::

It is immediate that \( \mathcal S\perpc \) is a subspace and that
\( \mathcal S\cap\mathcal S\perpc=\{\bzero\} \): a vector orthogonal to itself has
zero length. To test whether \( \bv\in\mathcal S\perpc \) it is enough to check
orthogonality to a spanning set of \( \mathcal S \). In particular,
\( \bv\in\C(\A)\perpc \) iff \( \A\T\bv=\bzero \).

::: {#lem-proj-null-colspace}
For any \( n\times p \) matrix \( \A \),
\[
\Null(\A\T)=\C(\A)\perpc
\qquad\text{and}\qquad
\Null(\A)=\C(\A\T)\perpc .
\]
:::

::: {.proof}
\( \bv\in\Null(\A\T) \) means \( \A\T\bv=\bzero \), that is, \( \bv \) is orthogonal to every
column of \( \A \), which is the same as \( \bv\perp\C(\A) \). The second identity is the
first applied to \( \A\T \).
:::

## Orthonormal bases and Gram–Schmidt

Vectors \( \bu_1,\dots,\bu_k \) are **orthonormal** if \( \bu_i\T\bu_j \) is \( 1 \) for
\( i=j \) and \( 0 \) otherwise. Orthonormal vectors are linearly independent: taking the
inner product of \( \sum_j a_j\bu_j=\bzero \) with \( \bu_i \) gives \( a_i=0 \). The
coordinates of a vector \( \bv\in\spn(\bu_1,\dots,\bu_k) \) are simply the inner
products \( a_i=\bu_i\T\bv \), with no linear system to solve. This is why
orthonormal bases are so convenient.

::: {#prp-proj-gram-schmidt}
[Gram–Schmidt]

Let \( \bv_1,\dots,\bv_k \) be linearly independent. Define recursively
\[
\bw_j=\bv_j-\sum_{i<j}(\bu_i\T\bv_j)\,\bu_i,\qquad
\bu_j=\bw_j/\norm{\bw_j},\qquad j=1,\dots,k.
\]
Then each \( \bw_j\neq\bzero \). The vectors \( \bu_1,\dots,\bu_k \) are orthonormal, and
\( \spn(\bu_1,\dots,\bu_j)=\spn(\bv_1,\dots,\bv_j) \) for every \( j \).
:::

::: {.proof}
Induction on \( j \). Suppose the claims hold for \( j-1 \). The vector \( \bw_j \) is \( \bv_j \)
minus an element of \( \spn(\bu_1,\dots,\bu_{j-1})=\spn(\bv_1,\dots,\bv_{j-1}) \). If
\( \bw_j \) were zero, \( \bv_j \) would lie in that span, contradicting independence.
For \( l<j \),
\( \bu_l\T\bw_j=\bu_l\T\bv_j-\sum_{i<j}(\bu_i\T\bv_j)\,\bu_l\T\bu_i
=\bu_l\T\bv_j-\bu_l\T\bv_j=0 \), so \( \bu_j \) is orthogonal to its predecessors.
Finally, \( \bu_j\in\spn(\bu_1,\dots,\bu_{j-1},\bv_j)=\spn(\bv_1,\dots,\bv_j) \) and
\( \bv_j\in\spn(\bu_1,\dots,\bu_j) \), so the two spans agree.
:::

Two consequences are used repeatedly. Every subspace has an orthonormal basis
(apply Gram–Schmidt to any basis). And any orthonormal basis of a subspace
\( \mathcal S \) can be extended to an orthonormal basis of \( \Real^n \): extend it to a
basis of \( \Real^n \) in any way, then apply Gram–Schmidt, which leaves the first
vectors unchanged. [Section 6.10](10-computation.html) shows that Gram–Schmidt, written
in matrix form, is the QR factorization.

## Direct sums and the orthogonal decomposition

For subspaces \( \mathcal S,\mathcal T \), the sum
\( \mathcal S+\mathcal T=\{\bu+\bv:\bu\in\mathcal S,\bv\in\mathcal T\} \) is a
subspace. When \( \mathcal S\cap\mathcal T=\{\bzero\} \) the sum is *direct*,
written \( \mathcal S\dirsum\mathcal T \). Every element then decomposes as \( \bu+\bv \)
in exactly one way: if \( \bu+\bv=\bu'+\bv' \) then \( \bu-\bu'=\bv'-\bv \) lies in both
subspaces and so is zero. When in addition \( \mathcal S\perp\mathcal T \), the sum is an
*orthogonal* direct sum.

::: {#thm-proj-direct-sum}
[Orthogonal decomposition]

For every subspace \( \mathcal S\subseteq\Real^n \),
\[
\Real^n=\mathcal S\dirsum\mathcal S\perpc,
\qquad
\dim\mathcal S+\dim\mathcal S\perpc=n,
\qquad
(\mathcal S\perpc)\perpc=\mathcal S .
\]
Hence every \( \y\in\Real^n \) can be written uniquely as \( \y=\bu+\bv \) with
\( \bu\in\mathcal S \) and \( \bv\perp\mathcal S \).
:::

::: {.proof}
Let \( \bu_1,\dots,\bu_r \) be an orthonormal basis of \( \mathcal S \), extended to an
orthonormal basis \( \bu_1,\dots,\bu_n \) of \( \Real^n \). A vector
\( \y=\sum_{i=1}^n a_i\bu_i \) is orthogonal to \( \mathcal S \) iff
\( a_1=\dots=a_r=0 \), so \( \mathcal S\perpc=\spn(\bu_{r+1},\dots,\bu_n) \) has
dimension \( n-r \). Writing
\( \y=\sum_{i\le r}a_i\bu_i+\sum_{i>r}a_i\bu_i \) exhibits the decomposition, and it
is unique because \( \mathcal S\cap\mathcal S\perpc=\{\bzero\} \). By the same
argument with the roles of the two index sets swapped,
\( (\mathcal S\perpc)\perpc=\spn(\bu_1,\dots,\bu_r)=\mathcal S \).
:::

Combining @thm-proj-direct-sum with @lem-proj-null-colspace gives a
matrix fact used in almost every later chapter.

::: {#cor-proj-gram-colspace}
For any \( n\times p \) matrix \( \A \):

::: {.enumerate options="label=(\alph*)"}
1. \( \Null(\A\T\A)=\Null(\A) \);

2. \( \C(\A\T\A)=\C(\A\T) \) and \( \C(\A\A\T)=\C(\A) \);

3. \( \rank(\A\T\A)=\rank(\A\A\T)=\rank(\A) \);

4. \( \C(\A)\perpc=\Null(\A\T) \) and \( \Real^p=\C(\A\T)\dirsum\Null(\A) \).
:::

:::

::: {.proof}
(a) If \( \A\bb=\bzero \) then \( \A\T\A\bb=\bzero \). Conversely, \( \A\T\A\bb=\bzero \) gives
\( \bb\T\A\T\A\bb=\norm{\A\bb}^2=0 \), so \( \A\bb=\bzero \).
(b) By @lem-proj-null-colspace applied to the symmetric matrix \( \A\T\A \),
and then (a) and @lem-proj-null-colspace again,
\( \C(\A\T\A)=\Null(\A\T\A)\perpc=\Null(\A)\perpc=(\C(\A\T)\perpc)\perpc=\C(\A\T) \),
using \( (\mathcal S\perpc)\perpc=\mathcal S \) in the last step. The second identity
is the first applied to \( \A\T \).
(c) Take dimensions in (b) and use \( \rank(\A\T)=\rank(\A) \).
(d) This is @lem-proj-null-colspace and @thm-proj-direct-sum.
:::

Part (b) is the reason the normal equations \( \X\T\X\bb=\X\T\y \) can always be
solved: the right-hand side lies in \( \C(\X\T)=\C(\X\T\X) \). We come back to this in
[Section 6.4](04-least-squares.html).

::: {#lem-proj-complement-intersection}
For subspaces \( \mathcal S_1,\mathcal S_2\subseteq\Real^n \),
\[
(\mathcal S_1+\mathcal S_2)\perpc=\mathcal S_1\perpc\cap\mathcal S_2\perpc
\qquad\text{and}\qquad
(\mathcal S_1\cap\mathcal S_2)\perpc=\mathcal S_1\perpc+\mathcal S_2\perpc .
\]
:::

::: {.proof}
A vector is orthogonal to \( \mathcal S_1+\mathcal S_2 \) iff it is orthogonal to
every \( \bu_1+\bu_2 \), which (taking \( \bu_2=\bzero \) or \( \bu_1=\bzero \)) is the case
iff it is orthogonal to both subspaces. For the second identity, apply the first
to \( \mathcal S_1\perpc \) and \( \mathcal S_2\perpc \):
\( (\mathcal S_1\perpc+\mathcal S_2\perpc)\perpc=\mathcal S_1\cap\mathcal S_2 \), and
take complements of both sides.
:::

@lem-proj-complement-intersection is what we need to describe reduced
models defined by linear constraints (@thm-proj-constraint-space).

## Exercises

### A. Check your understanding

::: {#exr-proj-complement-basic}
[A1]

Show that \( \mathcal S\perpc \) is a subspace for any set \( \mathcal S \) (not necessarily a
subspace), and that \( \mathcal S\subseteq\mathcal T \) implies
\( \mathcal T\perpc\subseteq\mathcal S\perpc \). Show that \( \mathcal S\perpc=\spn(\mathcal S)\perpc \).
:::

### B. Practice

::: {#exr-proj-oneway-spaces}
[B1]

A one-way layout has three groups of sizes \( 2,2,1 \) and model matrix
\( \X=[\bone,\bz_1,\bz_2,\bz_3] \). Find bases for \( \C(\X) \), \( \Null(\X) \), \( \C(\X\T) \) and
\( \C(\X)\perpc \), and check that the dimensions agree with
@eq-proj-rank-nullity and @thm-proj-direct-sum.
:::

::: {#exr-proj-colspace-product}
[B2]

Show that \( \C(\A\B)\subseteq\C(\A) \), and deduce
\( \rank(\A\B)\le\min\{\rank(\A),\rank(\B)\} \). Show that \( \C(\A\B)=\C(\A) \) whenever \( \B \)
has full row rank.
:::

### C. Going deeper

::: {#exr-proj-dim-formula}
[C1]

For subspaces \( \mathcal S,\mathcal T\subseteq\Real^n \), prove
\( \dim(\mathcal S+\mathcal T)=\dim\mathcal S+\dim\mathcal T-\dim(\mathcal S\cap\mathcal T) \).
Use it to show that two \( r \)-dimensional subspaces of \( \Real^n \) with \( 2r>n \) must share a
nonzero vector.
:::
