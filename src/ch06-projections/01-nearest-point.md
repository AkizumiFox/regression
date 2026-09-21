# Least squares as a nearest-point problem

[Chapter 5](../ch05-model-and-least-squares/index.html) introduced least squares as a calculus problem:
write down the sum of squared errors as a function of the coefficients,
differentiate, set the derivative to zero. That route gets the right answer, but
it hides the structure of the answer. It treats the coefficient vector \( \bbeta \)
as the basic unknown. It needs \( \X\T\X \) to be invertible before it can say
anything. And it gives no hint of why some quantities, such as fitted values,
residuals and sums of squares, behave well in every model, while others, such as
individual coefficients, can be fragile or not defined at all.

This chapter starts from a different place. Instead of the \( p \)-dimensional
space where \( \bbeta \) lives, it works in the \( n \)-dimensional space where the
data live.

## Two pictures of the same data

A regression data set with \( n \) observations and a single regressor is usually
drawn as a scatterplot: \( n \) points in the plane, one for each observation, with
the regressor on one axis and the response on the other. Call this the
**variable space** picture. It is the right picture for seeing the shape of a
relationship, but it scales badly. With \( k \) regressors it needs \( k+1 \) dimensions,
and it has no natural place for the model's coefficients.

The other picture turns this around. Stack the \( n \) responses into a single
vector \( \y=(y_1,\dots,y_n)\T \), and regard it as *one* point in \( \Real^n \).
Each column of the model matrix \( \X \) is also a vector in \( \Real^n \): the
intercept column is \( \bone \), and the column for a regressor holds its \( n \)
observed values. This is the **observation space** picture. Its dimension is
the sample size, which is usually far more than three. But it contains only
\( p+1 \) vectors of interest, and those always span a subspace of dimension at most
\( p+1 \). Most of the geometry happens in a subspace small enough to draw.

In observation space, the linear model
\[
\E(\Y) = \X\bbeta, \qquad \bbeta\in\Real^p,
\]
says one thing: the mean vector \( \E(\Y) \) lies in the set
\[
\C(\X) = \{\X\bb : \bb\in\Real^p\},
\]
the **column space** of \( \X \). That set is a subspace of \( \Real^n \). The
coefficients \( \bbeta \) are the coordinates of \( \E(\Y) \) with respect to the
columns of \( \X \), but the claim itself is only about *where* \( \E(\Y) \) is. It
does not depend on which spanning vectors are used to describe the subspace.

## The least squares principle, restated

The observed \( \y \) is almost never in \( \C(\X) \), because noise moves it off the
subspace. A natural estimate of \( \E(\Y) \) is the point of \( \C(\X) \) that is closest
to \( \y \). With Euclidean distance, this is exactly least squares:
\[
\norm{\y-\X\bb}^2 = \sum_{i=1}^n \bigl(y_i-\x_{(i)}\T\bb\bigr)^2,
\]
where \( \x_{(i)}\T \) is the \( i \)th row of \( \X \). Minimizing over \( \bb \) is the same as
searching \( \C(\X) \) for the point nearest to \( \y \).

::: {#def-proj-ls}
[Least squares]

A vector \( \hbeta\in\Real^p \) is a **least squares estimate** of \( \bbeta \) if
\[
\norm{\y-\X\hbeta} \le \norm{\y-\X\bb} \qquad\text{for every } \bb\in\Real^p.
\]
The vector \( \hY=\X\hbeta \) is the **fitted value** vector, and
\( \he=\y-\hY \) is the **residual** vector.
:::

The definition does not assume that the minimizer is unique, and deliberately
so. In [Section 6.4](04-least-squares.html) we will see that the *nearest point*
\( \hY \) is always unique. The coefficient vector \( \hbeta \) that produces it is
unique only when the columns of \( \X \) are linearly independent.

## A picture with three observations

The smallest case that still shows all the geometry has \( n=3 \) and \( p=2 \). Then
\( \C(\X) \) is a plane through the origin in \( \Real^3 \), and \( \y \) is a point
somewhere off the plane.

::: {#exm-proj-three-points}
[Three observations]

Take responses \( \y=(1,4,2)\T \) at regressor values \( \x=(0,1,3)\T \), with an
intercept, so that \( \X=[\bone,\x] \). The least squares fit (computed by the
listing below) is
\[
\hY = (2.143,\ 2.286,\ 2.571)\T,\qquad
\he = (-1.143,\ 1.714,\ -0.571)\T,
\]
with intercept \( 2.143 \) and slope \( 0.143 \). Two facts
stand out. First, the residual is orthogonal to both columns of \( \X \): the
entries of \( \he \) sum to zero, and \( \x\T\he = 0 \). Second, the squared lengths add:
\[
\norm{\y}^2=21.000 = 16.429 + 4.571
= \norm{\hY}^2+\norm{\he}^2 .
\]
[Figure 6.1.1](01-nearest-point.html#fig-proj-three-points) shows the configuration. The fitted vector is the
foot of the perpendicular dropped from \( \y \) to the plane.
:::

::: {when-format="html"}
![**Figure 6.1.1.** Least squares with \( n=3 \) observations. The model's column space
\( \C(\X) \) is the plane spanned by \( \bone \) and \( \x \). The fitted vector
\( \hY=\M\y \) is the point of the plane nearest to \( \y \), and the residual
\( \he=\y-\hY \) is perpendicular to the plane.](ls_geometry_3d.svg){#fig-proj-three-points width=58%}
:::

::: {when-format="pdf"}
![Least squares with \( n=3 \) observations. The model's column space
\( \C(\X) \) is the plane spanned by \( \bone \) and \( \x \). The fitted vector
\( \hY=\M\y \) is the point of the plane nearest to \( \y \), and the residual
\( \he=\y-\hY \) is perpendicular to the plane.](ls_geometry_3d.pdf){width=58%}
:::

```{.python .run #cell-ls-geometry-3d-setup}
import numpy as np

y = np.array([1.0, 4.0, 2.0])
X = np.column_stack([np.ones(3), [0.0, 1.0, 3.0]])   # intercept and one regressor

beta_hat, *_ = np.linalg.lstsq(X, y, rcond=None)
y_hat = X @ beta_hat              # the point of C(X) nearest to y
e_hat = y - y_hat                 # the residual vector

print("fitted  ", y_hat)
print("residual", e_hat)
print("X^T e   ", X.T @ e_hat)    # zero: the residual is orthogonal to C(X)
```

The two facts in @exm-proj-three-points are not coincidences of this data
set. They are the whole theory in miniature. The residual is orthogonal to the
subspace *because* \( \hY \) is the nearest point, and the squared lengths add
*because* the two pieces are orthogonal. The rest of the chapter makes
these two statements precise and follows their consequences.

::: {.idea}
Fitted values are determined by the *subspace* \( \C(\X) \). Coefficients are
determined by the *basis* used to describe it. Anything that depends only
on \( \C(\X) \) is unchanged by reparameterization, by rescaling or recoding
regressors, and by rank deficiency. This covers fitted values, residuals, sums
of squares, \( R^2 \) and leverages.
:::

## Plan of the chapter

[Section 6.2](02-subspaces.html) collects the facts about subspaces that the
geometry needs. [Section 6.3](03-orthogonal-projection.html) proves that nearest points
exist and are unique, and identifies the matrices that compute them: the
symmetric idempotent matrices, including \( \X(\X\T\X)\ginv\X\T \) for every
generalized inverse. [Section 6.4](04-least-squares.html) applies this to least
squares. [Section 6.5](05-nested.html) studies nested subspaces, which is where sums
of squares come from. [Section 6.6](06-fwl.html) proves the Frisch–Waugh–Lovell
theorem, which says what a single coefficient in a multiple regression measures.
[Section 6.7](07-reparameterization.html) covers reparameterization and \( R^2 \),
and [Section 6.8](08-leverage.html) covers leverage. [Section 6.9](09-inner-products.html)
replaces the Euclidean inner product with a general one, the step that leads to
generalized least squares. [Section 6.10](10-computation.html) shows how the projection
should be computed in floating point.

## Exercises

### A. Check your understanding

::: {#exr-proj-three-points-exact}
[A1]

Redo @exm-proj-three-points in exact arithmetic. From the simple regression formulas of
@thm-lm-simple-ls, show that the intercept is \( 15/7 \) and the slope \( 1/7 \), so that
\( \hY=(15,16,18)\T/7 \) and \( \he=(-8,12,-4)\T/7 \). Check exactly that \( \bone\T\he=0 \),
\( \x\T\he=0 \) and \( \norm{\hY}^2+\norm{\he}^2=\norm{\y}^2=21 \).
:::

::: {.solution}
Here \( \sum x_i=4 \), \( \sum x_i^2=10 \), \( \sum y_i=7 \) and \( \sum x_iy_i=10 \), so
\( S_{xx}=10-4^2/3=14/3 \) and \( S_{xy}=10-4\cdot7/3=2/3 \). The slope is
\( S_{xy}/S_{xx}=1/7 \) and the intercept \( \bar y-\bar x/7=7/3-4/21=15/7 \). Evaluating
\( 15/7+x_i/7 \) at \( x=0,1,3 \) gives \( \hY \), and \( \he=\y-\hY \) as stated. Then
\( \bone\T\he=(-8+12-4)/7=0 \) and \( \x\T\he=(0\cdot(-8)+1\cdot12+3\cdot(-4))/7=0 \). Finally
\( \norm{\hY}^2=(225+256+324)/49=805/49 \) and \( \norm{\he}^2=(64+144+16)/49=224/49 \), which
sum to \( 1029/49=21=1+16+4 \).
:::

::: {#exr-proj-duplicate-columns}
[A2]

Let \( \y=(1,3)\T \) and let \( \X \) be the \( 2\times2 \) matrix whose two columns both equal
\( \bone \). Show that \( \norm{\y-\X\bb}^2 \) depends on \( \bb \) only through \( b_1+b_2 \), and that
the least squares estimates of @def-proj-ls are exactly the \( \bb \) with \( b_1+b_2=2 \). Which
of \( \hbeta \), \( \hY \), \( \he \) does the definition determine uniquely here?
:::

::: {.solution}
\( \X\bb=(b_1+b_2)\bone \), so with \( t=b_1+b_2 \) the criterion is
\( (1-t)^2+(3-t)^2=2(t-2)^2+2 \), minimized exactly at \( t=2 \). The estimates therefore form
the line \( \{\bb:b_1+b_2=2\} \), while \( \hY=(2,2)\T \) and \( \he=(-1,1)\T \) are the same for every
one of them. The fitted vector and the residual are determined; the coefficient vector is not.
:::

::: {#exr-proj-equivariance}
[A3]

Let \( \hY \) be a nearest point to \( \y \) in \( \C(\X) \). Show directly from @def-proj-ls, with no
formula for \( \hY \), that \( c\,\hY \) is a nearest point to \( c\,\y \) for every \( c\in\Real \), and
that \( \hY+\bu \) is a nearest point to \( \y+\bu \) for every \( \bu\in\C(\X) \). Deduce that adding
\( \bu \) leaves the residual unchanged.
:::

### B. Practice

::: {#exr-proj-unique-fit}
[B1]

The nearest point is unique, and this can be proved before any of the machinery of
[Section 6.3](03-orthogonal-projection.html). Let \( \bu,\bv\in\C(\X) \) both attain
\( d=\min_{\bb}\norm{\y-\X\bb} \). Expand inner products, or use the parallelogram law of
@exr-mat-polarization, to verify
\[
\Bigl\lVert\y-\tfrac{\bu+\bv}{2}\Bigr\rVert^2
=\tfrac12\norm{\y-\bu}^2+\tfrac12\norm{\y-\bv}^2-\tfrac14\norm{\bu-\bv}^2 ,
\]
and conclude that \( \bu=\bv \). Where does the argument use that \( \C(\X) \) contains the
midpoint of any two of its elements?
:::

::: {.solution}
Write \( \y-(\bu+\bv)/2=\tfrac12(\y-\bu)+\tfrac12(\y-\bv) \) and expand:
\[
\Bigl\lVert\y-\tfrac{\bu+\bv}{2}\Bigr\rVert^2
=\tfrac14\norm{\y-\bu}^2+\tfrac14\norm{\y-\bv}^2+\tfrac12(\y-\bu)\T(\y-\bv).
\]
Also \( \norm{\bu-\bv}^2=\norm{(\y-\bv)-(\y-\bu)}^2
=\norm{\y-\bu}^2+\norm{\y-\bv}^2-2(\y-\bu)\T(\y-\bv) \). Eliminating the cross term between the
two displays gives the identity. Now \( \norm{\y-\bu}=\norm{\y-\bv}=d \), so the left side is
\( d^2-\tfrac14\norm{\bu-\bv}^2 \). The midpoint lies in \( \C(\X) \), so the left side is at least
\( d^2 \) by the definition of \( d \); hence \( \norm{\bu-\bv}=0 \). Without the midpoint the
inequality \( \lVert\y-(\bu+\bv)/2\rVert\ge d \) would not be available.
:::

::: {#exr-proj-no-intercept-shift}
[B2]

Adding a constant to every response is the commonest change of units, and only a model with
an intercept handles it in the obvious way.

::: {.enumerate options="label=(\alph*)"}
1. Suppose \( \bone\in\C(\X) \). Using @exr-proj-equivariance, show that the fit for
   \( \y+c\bone \) is \( \hY+c\bone \) and that the residual is unchanged.

2. Take \( n=2 \), \( \X=\x=(1,2)\T \) and no intercept. Show that the fit for \( \y=(1,1)\T \)
   is \( (3/5,6/5)\T \), and that the fit for \( \y+\bone=(2,2)\T \) is \( (6/5,12/5)\T \), not
   \( (8/5,11/5)\T \). Which residual changes, and why is there no contradiction with (a)?
:::
:::

::: {.solution}
(a) When \( \bone\in\C(\X) \), the map \( \bv\mapsto\bv+c\bone \) is a bijection of \( \C(\X) \)
onto itself, and \( \norm{(\y+c\bone)-(\bv+c\bone)}=\norm{\y-\bv} \), so the distances to be
minimized match up in pairs. The nearest point to \( \y+c\bone \) is therefore \( \hY+c\bone \),
and subtracting leaves \( \he \) unchanged. (This is @exr-proj-equivariance with
\( \bu=c\bone \).)

(b) Minimizing \( (1-b)^2+(1-2b)^2 \) gives \( 10b=6 \), \( b=3/5 \), and
\( b\x=(3/5,6/5)\T \). For \( \y=(2,2)\T \), \( 10b=12 \), \( b=6/5 \), and \( b\x=(6/5,12/5)\T \). The
residual changes from \( (2/5,-1/5)\T \) to \( (4/5,-2/5)\T \). There is no contradiction:
\( \bone\notin\C(\x) \), so \( \y+\bone \) is not \( \y \) plus an element of the model space, and (a)
does not apply. The intercept column is what makes a shift of the response a movement inside
the model space.
:::

::: {#exr-proj-nested-improvement}
[B3]

Let \( \C(\X_1)\subseteq\C(\X_2) \), and let \( d_k=\min\{\norm{\y-\bv}:\bv\in\C(\X_k)\} \).

::: {.enumerate options="label=(\alph*)"}
1. Show \( d_2\le d_1 \), with equality iff some nearest point of \( \y \) in \( \C(\X_1) \) is also a
   nearest point in \( \C(\X_2) \). Conclude that appending a column to \( \X \) can never raise the
   residual sum of squares.

2. Equality does not require the new column to be redundant. With \( n=3 \), \( \X=\bone \) and
   \( \y=(1,2,3)\T \), the fit is \( \hY=(2,2,2)\T \) and \( \he=(-1,0,1)\T \). Show that appending
   \( \bz=(1,-2,1)\T \) enlarges the column space but leaves \( \hY \) and \( d \) unchanged.
:::
:::

### C. Going deeper

::: {#exr-proj-existence}
[C1]

@def-proj-ls takes for granted that a minimizer exists. Prove it, without orthonormal bases.

::: {.enumerate options="label=(\alph*)"}
1. Reduce to the case in which the columns of \( \X \) are linearly independent.

2. Show that \( c=\min\{\norm{\X\bb}:\norm{\bb}=1\} \) is positive, so \( \norm{\X\bb}\ge c\norm{\bb} \)
   for every \( \bb \).

3. Deduce that the infimum of \( \norm{\y-\X\bb} \) over \( \Real^p \) equals its infimum over a
   closed ball, and that it is attained.

4. Show that the conclusion can fail for sets that are not subspaces: give a nonempty
   \( \mathcal S\subseteq\Real^2 \) and a \( \y \) with no nearest point in \( \mathcal S \).
:::
:::

::: {.solution}
(a) Discard columns that are linear combinations of the others until the rest are independent.
This changes neither \( \C(\X) \) nor the set of attainable values \( \X\bb \), so it changes
neither the infimum nor whether it is attained.

(b) The function \( \bb\mapsto\norm{\X\bb} \) is continuous and the unit sphere of \( \Real^p \) is
closed and bounded, so the minimum \( c \) is attained at some \( \bb_0 \) with \( \norm{\bb_0}=1 \).
If \( c=0 \) then \( \X\bb_0=\bzero \) with \( \bb_0\ne\bzero \), contradicting independence. Homogeneity
gives \( \norm{\X\bb}\ge c\norm{\bb} \) for all \( \bb \).

(c) At \( \bb=\bzero \) the criterion is \( \norm{\y} \). If \( \norm{\bb}>2\norm{\y}/c \), then
\( \norm{\y-\X\bb}\ge\norm{\X\bb}-\norm{\y}\ge c\norm{\bb}-\norm{\y}>\norm{\y} \). So the infimum
over \( \Real^p \) equals the infimum over the closed ball \( \norm{\bb}\le2\norm{\y}/c \), which is
closed and bounded; a continuous function attains its minimum there.

(d) Let \( \mathcal S=\{\bv:\norm{\bv}<1\} \) and \( \y=(2,0)\T \). Distances from \( \y \) to
\( \mathcal S \) have infimum \( 1 \), approached along \( \bv=(t,0)\T \) as \( t\uparrow1 \), but no point of
\( \mathcal S \) is at distance \( 1 \). What fails is that \( \mathcal S \) is not closed; it is
both convex and bounded.
:::

::: {#exr-proj-coefficient-set}
[C2]

Write \( \Null(\X)=\{\mathbf{d}\in\Real^p:\X\mathbf{d}=\bzero\} \), and let \( \hbeta_0 \) be any least
squares estimate.

::: {.enumerate options="label=(\alph*)"}
1. Using @exr-proj-unique-fit, show that the set of least squares estimates is exactly
   \( \hbeta_0+\Null(\X) \).

2. Deduce that the estimate is unique iff \( \X \) has full column rank, and that otherwise the
   solution set is a translate of a subspace of dimension \( p-\rank(\X) \) (@thm-mat-rank-nullity).

3. Show that for \( \mathbf{c}\in\Real^p \) the number \( \mathbf{c}\T\hbeta \) is the same for every least
   squares estimate iff \( \mathbf{c}\T\mathbf{d}=0 \) for every \( \mathbf{d}\in\Null(\X) \).

4. Which \( \mathbf{c} \) qualify in @exr-proj-duplicate-columns?
:::

[Chapter 8](../ch08-estimability/index.html) takes the condition in (c) as the definition of an
estimable function.
:::

::: {.solution}
(a) By @exr-proj-unique-fit there is a single vector \( \hY \) attaining the minimum, so \( \bb \)
is a least squares estimate iff \( \X\bb=\hY=\X\hbeta_0 \), that is, iff
\( \X(\bb-\hbeta_0)=\bzero \).

(b) The set is a single point iff \( \Null(\X)=\{\bzero\} \), which is full column rank; its
dimension is \( p-\rank(\X) \) by the rank–nullity theorem.

(c) Every estimate is \( \hbeta_0+\mathbf{d} \), so the values taken by \( \mathbf{c}\T\hbeta \) are
\( \mathbf{c}\T\hbeta_0+\mathbf{c}\T\mathbf{d} \) as \( \mathbf{d} \) ranges over \( \Null(\X) \). These agree for all
\( \mathbf{d} \) iff \( \mathbf{c}\T\mathbf{d}=0 \) throughout \( \Null(\X) \); if \( \mathbf{c}\T\mathbf{d}\ne0 \) for some
\( \mathbf{d} \), the two estimates \( \hbeta_0 \) and \( \hbeta_0+\mathbf{d} \) give different values.

(d) There \( \Null(\X) \) is spanned by \( (1,-1)\T \), so the condition is \( c_1=c_2 \): the sum
\( \beta_1+\beta_2 \) is determined by the data, and no other linear function is, apart from its
multiples.
:::
