# Systems of equations and consistency

A system \( \A\x=\bb \) with \( \A \) of size \( m\times n \) is **consistent** if it has at least
one solution. Least squares leads to such systems in which \( \A=\X\T\X \) may be singular, so
we need to know when solutions exist, what all of them look like, and which features are
common to all of them.

::: {#thm-mat-consistency}
[Consistency and the general solution]

Let \( \G \) be any generalized inverse of \( \A \). The following are equivalent:

::: {.enumerate options="label=(\alph*)"}
1. \( \A\x=\bb \) is consistent;

2. \( \bb\in\C(\A) \);

3. \( \rank[\A,\bb]=\rank(\A) \);

4. \( \A\G\bb=\bb \).
:::

If the system is consistent, \( \x_0=\G\bb \) is a solution, and the set of all solutions is
\[
\{\G\bb+(\I-\G\A)\bz:\bz\in\Real^n\}=\x_0+\Null(\A).
\]{#eq-mat-general-solution}

:::

::: {.proof}
(a)\( \Leftrightarrow \)(b) is the first reading in @eq-mat-product-readings.
(b)\( \Leftrightarrow \)(c): appending \( \bb \) as a column leaves the column space unchanged iff
\( \bb \) is already in it, and a subspace of \( \C([\A,\bb]) \) with the same dimension equals
it. (d)\( \Rightarrow \)(a): \( \G\bb \) is a solution. (b)\( \Rightarrow \)(d): if \( \bb=\A\bu \), then
\( \A\G\bb=\A\G\A\bu=\A\bu=\bb \).

The solutions of a consistent system are \( \x_0+\Null(\A) \), since \( \A\x=\bb \) iff
\( \A(\x-\x_0)=\bzero \). It remains to show \( \Null(\A)=\C(\I-\G\A) \). The inclusion
\( \supseteq \) holds because \( \A(\I-\G\A)=\bzero \). Conversely, if \( \A\bv=\bzero \), then
\( \bv=(\I-\G\A)\bv \).
:::

Criterion (d) is useful in practice: it tests consistency with a single matrix product,
for any generalized inverse. By @prp-mat-ginverse-props(a) and rank–nullity,
\( \rank(\I-\G\A)=n-r \), so the solution set is an affine subspace of dimension \( n-r \). The
solution is unique iff \( \rank(\A)=n \), and then \( \G\bb \) is the same for every \( \G \).

::: {.idea}
Any generalized inverse does three jobs for \( \A\x=\bb \). It tests consistency
(\( \A\G\bb=\bb \)), produces one solution (\( \G\bb \)), and parametrizes all of them
(\( \G\bb+(\I-\G\A)\bz \)). Which \( \G \) is used changes the particular solution but never
the solution set, and a linear function \( \mathbf{q}\T\x \) is the same across solutions exactly
when \( \mathbf{q}\in\C(\A\T) \).
:::

::: {#cor-mat-invariant}
[Invariant linear functions]

Let \( \A\x=\bb \) be consistent and let \( \mathbf{q}\in\Real^n \). The value \( \mathbf{q}\T\x \) is the same
for every solution \( \x \) iff \( \mathbf{q}\in\C(\A\T) \), that is, iff \( \mathbf{q}\T=\mathbf{q}\T\G\A \) for one
(and then every) generalized inverse \( \G \).
:::

::: {.proof}
By @eq-mat-general-solution, \( \mathbf{q}\T\x \) is constant over solutions iff
\( \mathbf{q}\T(\I-\G\A)\bz=0 \) for all \( \bz \), that is, iff \( \mathbf{q}\T=\mathbf{q}\T\G\A \). That equation puts
\( \mathbf{q} \) in \( \C(\A\T) \). Conversely, if \( \mathbf{q}=\A\T\bw \), then
\( \mathbf{q}\T\G\A=\bw\T\A\G\A=\bw\T\A=\mathbf{q}\T \).
:::

Applied to the normal equations \( \X\T\X\bb=\X\T\y \), which are consistent by
@prp-mat-rank-product(c), the corollary says that \( \mathbf{q}\T\hbeta \) is the same for all
least squares solutions iff \( \mathbf{q}\in\C(\X\T\X)=\C(\X\T) \). This is the algebra behind
estimable functions (@thm-est-characterization).

::: {#prp-mat-min-norm}
[Minimum-norm solution]

If \( \A\x=\bb \) is consistent, then \( \A^+\bb \) is its unique solution lying in \( \C(\A\T) \),
and \( \norm{\A^+\bb}<\norm{\x} \) for every other solution \( \x \).
:::

::: {.proof}
\( \A^+\bb \) is a solution because \( \A^+ \) is a generalized inverse. From the second and fourth
Penrose conditions, \( \A^+=(\A^+\A)\A^+=(\A^+\A)\T\A^+=\A\T(\A^+)\T\A^+ \), so
\( \A^+\bb\in\C(\A\T) \). Every vector of \( \C(\A\T) \) is orthogonal to every
\( \bv\in\Null(\A) \), since \( (\A\T\bw)\T\bv=\bw\T\A\bv=0 \). Any other solution is
\( \x=\A^+\bb+\bv \) with \( \bzero\neq\bv\in\Null(\A) \), and then
\( \norm{\x}^2=\norm{\A^+\bb}^2+\norm{\bv}^2 \). The same computation shows that \( \x\notin\C(\A\T) \),
because \( \x\T\bv=\norm{\bv}^2\ne0 \).
:::

::: {#exm-mat-margins}
[Filling a table from its margins]

Suppose the six entries of a \( 2\times3 \) table are unknown, but its row totals \( (9,6) \) and
column totals \( (4,5,6) \) are known. With the cells ordered
\( x_{11},x_{12},x_{13},x_{21},x_{22},x_{23} \), the equations are \( \A\x=\bb \) with
\[
\A=\begin{pmatrix}
    1&1&1&0&0&0\\ 0&0&0&1&1&1\\ 1&0&0&1&0&0\\ 0&1&0&0&1&0\\ 0&0&1&0&0&1
\end{pmatrix},\qquad
\bb=\begin{pmatrix}9\\6\\4\\5\\6\end{pmatrix}.
\]
The rows of \( \A \) satisfy one relation: the two row-total rows and the three
column-total rows both add to \( \bone_6\T \). So \( \rank(\A)=4 \) and there is a
\( 6-4=2 \)-dimensional family of solutions. The same relation must hold for \( \bb \). By
@thm-mat-consistency(c), the system is consistent iff the row totals and the column
totals have the same sum, and here both sums are \( 15 \).

The listing builds three generalized inverses: \( \G_1 \) from the nonsingular submatrix in
rows 1–4 and columns 1–4 (@thm-mat-ginverse-exists), the Moore–Penrose inverse
\( \G_2=\A^+ \), and \( \G_3=\G_1+(\I-\G_1\A)\mathbf{Z} \) for an integer matrix \( \mathbf{Z} \)
(@prp-mat-ginverse-props(e)). They give the solutions
\[
\begin{aligned}  \G_1\bb&=(-2,\ 5,\ 6,\ 6,\ 0,\ 0)\T, & \norm{\G_1\bb}&=10.050,\\
\G_2\bb&=(5/2,\ 3,\ 7/2,\ 3/2,\ 2,\ 5/2)\T, & \norm{\G_2\bb}&=6.325,\\
\G_3\bb&=(21,\ -}39,\ 27,\ \ensuremath{-}17,\ 44,\ \ensuremath{-21)\T, & \norm{\G_3\bb}&=73.192.\end{aligned}
\]
The individual cell \( x_{11} \) takes the values \( -2 \), \( 5/2 \) and
\( 21 \): margins do not determine cells. But
\( x_{12}+x_{13}-x_{21} \) equals \( 5 \) for all three, as
@cor-mat-invariant predicts, because its coefficient vector is the first row of \( \A \)
minus the third. The minimum-norm solution has the additive form
\( x_{ij}=r_i/3+c_j/2-15/6 \) in terms of the row totals \( r_i \) and column totals \( c_j \)
(@exr-mat-margins-min-norm). If the last column total is changed to \( 7 \), then
\( \A\G\bb\neq\bb \) for each of the three matrices, and the script confirms that the test
detects the inconsistency. The script also checks all \( 60 \) submatrices
formed by \( 4 \) independent rows and \( 4 \) independent columns, and finds each one nonsingular, as
@prp-mat-nonsingular-submatrix requires.
:::

```{.python .run #cell-ginverse-system-ginverses}
from fractions import Fraction
from itertools import combinations
import numpy as np
m, n = A.shape
# cells ordered (1,1), (1,2), (1,3), (2,1), (2,2), (2,3)
A = np.array([[1, 1, 1, 0, 0, 0],      # row 1 total
              [0, 0, 0, 1, 1, 1],      # row 2 total
              [1, 0, 0, 1, 0, 0],      # column 1 total
              [0, 1, 0, 0, 1, 0],      # column 2 total
              [0, 0, 1, 0, 0, 1]],     # column 3 total
             dtype=float)
b = np.array([9, 6, 4, 5, 6], dtype=float)       # both sets of totals add to 15
r = np.linalg.matrix_rank(A)
print("rank(A) =", r, " rank([A, b]) =", np.linalg.matrix_rank(np.column_stack([A, b])))

def ginverse_from_submatrix(A, rows, cols):
    """Place the inverse of the nonsingular block A[rows, cols] at [cols, rows]."""
    G = np.zeros(A.T.shape)
    G[np.ix_(cols, rows)] = np.linalg.inv(A[np.ix_(rows, cols)])
    return G

rows, cols = [0, 1, 2, 3], [0, 1, 2, 3]                  # r independent rows and columns
G1 = ginverse_from_submatrix(A, rows, cols)
G2 = np.linalg.pinv(A)                                   # Moore-Penrose inverse
rng = np.random.default_rng(1)
Z = rng.integers(-2, 3, size=A.T.shape).astype(float)
G3 = G1 + (np.eye(n) - G1 @ A) @ Z                       # another, non-reflexive one

for G in (G1, G2, G3):
    assert np.allclose(A @ G @ A, A)                     # each is a g-inverse
    assert np.allclose(A @ G @ b, b)                     # consistency test passes
    x = G @ b
    print(np.round(x, 4), " A x - b =", np.abs(A @ x - b).max())
```

## Exercises

### A. Check your understanding

::: {#exr-mat-vector-ginverse}
[A1]

Find all generalized inverses of a nonzero vector \( \mathbf{a}\in\Real^n \), regarded as an
\( n\times1 \) matrix. Which one is \( \mathbf{a}^+ \)?
:::

::: {#exr-mat-margins-min-norm}
[A2]

Generalize @exm-mat-margins to an \( a\times b \) table with row totals \( r_i \), column totals
\( c_j \) and grand total \( T=\sum_ir_i=\sum_jc_j \). Show that \( x_{ij}=r_i/b+c_j/a-T/(ab) \) is a solution
lying in \( \C(\A\T) \), and hence is the minimum-norm solution.
:::

::: {.solution}
Row \( i \) of the proposed table sums to
\( r_i+\sum_jc_j/a-bT/(ab)=r_i \), and column \( j \) sums to \( \sum_ir_i/b+c_j-T/b=c_j \). The rows of \( \A \)
are the row indicators and the column indicators, so \( \C(\A\T) \) consists of the tables of the
form \( x_{ij}=\alpha_i+\beta_j \). The proposed solution has this form. By
@prp-mat-min-norm, it is \( \A^+\bb \). For \( a=2 \), \( b=3 \) and \( T=15 \) this is the additive
formula quoted in @exm-mat-margins.
:::

::: {#exr-mat-consistency-numeric}
[A3]

For which \( c \) and \( d \) is the system \( x_1+x_2+x_3=2 \), \( x_1-x_2=c \), \( 3x_1+x_2+2x_3=d \) consistent?
When it is, give its general solution in the form @eq-mat-general-solution.
:::

### B. Practice

::: {#exr-mat-reflexive}
[B1]

Let \( \G \) be a generalized inverse of \( \A \). Show that \( \G\A\G=\G \) iff \( \rank(\G)=\rank(\A) \).
:::

::: {.solution}
If \( \G\A\G=\G \), then \( \rank(\G)\le\rank(\A) \) by
@prp-mat-rank-product(a), and \( \rank(\A)\le\rank(\G) \) by
@prp-mat-ginverse-props(a). Conversely, suppose \( \rank(\G)=\rank(\A)=r \). Since
\( \C(\G\A)\subseteq\C(\G) \) and \( \rank(\G\A)=r \), the two column spaces coincide, so
\( \G=\G\A\mathbf{K} \) for some \( \mathbf{K} \). Then \( \G\A\G=\G\A\G\A\mathbf{K}=\G\A\mathbf{K}=\G \), using
\( \A\G\A=\A \).
:::

::: {#exr-mat-mp-factorization}
[B2]

Let \( \A=\B\mathbf{F} \) be a rank factorization. Show that
\( \A^+=\mathbf{F}\T(\mathbf{F}\mathbf{F}\T)^{-1}(\B\T\B)^{-1}\B\T \). Deduce that \( \A^+=(\A\T\A)^{-1}\A\T \) when \( \A \)
has full column rank.
:::

::: {.solution}
Write \( \G \) for the proposed matrix. Then
\( \A\G=\B(\B\T\B)^{-1}\B\T \) and \( \G\A=\mathbf{F}\T(\mathbf{F}\mathbf{F}\T)^{-1}\mathbf{F} \), both symmetric. Also
\( \A\G\A=\B(\B\T\B)^{-1}\B\T\B\mathbf{F}=\A \) and
\( \G\A\G=\mathbf{F}\T(\mathbf{F}\mathbf{F}\T)^{-1}\mathbf{F}\mathbf{F}\T(\mathbf{F}\mathbf{F}\T)^{-1}(\B\T\B)^{-1}\B\T=\G \). By uniqueness
\( \G=\A^+ \). If \( \A \) has full column rank, take \( \B=\A \) and \( \mathbf{F}=\I \).
:::

::: {#exr-mat-independent-solutions}
[B3]

Let \( \A \) be \( m\times n \) of rank \( r \) and let \( \A\x=\bb \) be consistent with \( \bb\ne\bzero \). Show that
the system has \( n-r+1 \) linearly independent solutions and no more, and that every solution is
an affine combination of them.
:::

### C. Going deeper

::: {#exr-mat-mp-properties}
[C1]

Show that \( (\A^+)^+=\A \), \( (\A\T)^+=(\A^+)\T \), \( (\A\T\A)^+=\A^+(\A^+)\T \) and
\( \A^+=(\A\T\A)^+\A\T \). Give an example with \( (\A\B)^+\ne\B^+\A^+ \).
:::
