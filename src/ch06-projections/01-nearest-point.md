# Least squares as a nearest-point problem

Chapter 5 introduced least squares as a calculus problem:
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
