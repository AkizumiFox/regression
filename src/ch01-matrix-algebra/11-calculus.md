# Vector and matrix calculus

## Gradients of linear and quadratic functions

For a differentiable scalar function \( f \) of \( \x\in\Real^n \), the **gradient**
\( \partial f/\partial\x \) is the column vector with entries \( \partial f/\partial x_i \), and the
**Hessian** \( \partial^2f/\partial\x\,\partial\x\T \) is the \( n\times n \) matrix of second
derivatives. For a vector-valued \( \bm g(\x)\in\Real^m \), the \( m\times n \) *Jacobian*
\( \partial\bm g/\partial\x\T \) has \( (i,j) \) entry \( \partial g_i/\partial x_j \).

::: {#thm-mat-quadform-derivative}
[Derivatives of linear and quadratic forms]

Let \( \bm a\in\Real^n \), let \( \A \) be \( n\times n \) and \( \B \) be \( m\times n \). Then
\[
\frac{\partial(\bm a\T\x)}{\partial\x}=\bm a,\qquad
\frac{\partial(\B\x)}{\partial\x\T}=\B,\qquad
\frac{\partial(\x\T\A\x)}{\partial\x}=(\A+\A\T)\x,\qquad
\frac{\partial^2(\x\T\A\x)}{\partial\x\,\partial\x\T}=\A+\A\T .
\]
In particular \( \partial(\x\T\A\x)/\partial\x=2\A\x \) when \( \A \) is symmetric.
:::

::: {.proof}
The first two are immediate from \( \bm a\T\x=\sum_ia_ix_i \) and \( (\B\x)_i=\sum_jb_{ij}x_j \). For the
quadratic form, expand
\[
(\x+\bm h)\T\A(\x+\bm h)=\x\T\A\x+\bm h\T(\A+\A\T)\x+\bm h\T\A\bm h .
\]
The middle term is linear in \( \bm h \) and the last is \( O(\norm{\bm h}^2) \), so \( (\A+\A\T)\x \) is
the gradient. The gradient is linear in \( \x \) with matrix \( \A+\A\T \), which is the Hessian by
the second formula.
:::

::: {#exm-mat-ls-gradient}
[The least squares criterion]

Let \( S(\bb)=\norm{\y-\X\bb}^2=\y\T\y-2\y\T\X\bb+\bb\T\X\T\X\bb \). By
@thm-mat-quadform-derivative,
\[
\frac{\partial S}{\partial\bb}=-2\X\T\y+2\X\T\X\bb=-2\X\T(\y-\X\bb),\qquad
\frac{\partial^2S}{\partial\bb\,\partial\bb\T}=2\X\T\X .
\]
Setting the gradient to zero gives the normal equations \( \X\T\X\bb=\X\T\y \). The Hessian is
nonnegative definite, so \( S \) is convex and every solution is a global minimizer. The next
result proves this without calculus, and covers the rank-deficient case.
:::

::: {#prp-mat-quadratic-min}
[Minimizing a quadratic function]

Let \( \A \) be \( n\times n \) nonnegative definite, \( \bb\in\Real^n \) and
\( f(\x)=\x\T\A\x-2\bb\T\x \). If \( \bb\in\C(\A) \), the minimizers of \( f \) are exactly the solutions
of \( \A\x=\bb \), and the minimum value is \( -\bb\T\A\ginv\bb \) for any generalized inverse
\( \A\ginv \). If \( \bb\notin\C(\A) \), then \( f \) is unbounded below.
:::

::: {.proof}
Let \( \bb\in\C(\A) \) and let \( \x_0 \) solve \( \A\x_0=\bb \). Completing the square,
\[
f(\x)=(\x-\x_0)\T\A(\x-\x_0)-\x_0\T\A\x_0 .
\]
The first term is nonnegative, and by @prp-mat-pd-properties(b) it is zero iff
\( \A(\x-\x_0)=\bzero \), that is, iff \( \A\x=\bb \). The minimum is
\( -\x_0\T\A\x_0=-\x_0\T\A\A\ginv\A\x_0=-\bb\T\A\ginv\bb \). If \( \bb\notin\C(\A) \), write
\( \bb=\bb_1+\bb_2 \) with \( \bb_1\in\C(\A) \) and \( \bzero\ne\bb_2\in\C(\A)\perpc=\Null(\A) \)
(@prp-mat-complement). Then \( f(t\bb_2)=-2t\norm{\bb_2}^2\to-\infty \) as \( t\to\infty \).
:::

With \( \A=\X\T\X \) and \( \bb=\X\T\y \), the proposition gives least squares for every model matrix:
the minimizers are the solutions of the normal equations, and the minimum of \( S \) is
\( \y\T\y-\y\T\X(\X\T\X)\ginv\X\T\y \).

## Constrained minimization

To minimize a differentiable \( f(\x) \) subject to linear constraints \( \bm H\x=\bm h \), with
\( \bm H \) of size \( q\times n \), the method of **Lagrange multipliers** looks for stationary
points of \( f(\x)+2\bm\theta\T(\bm H\x-\bm h) \) in \( (\x,\bm\theta) \). The factor \( 2 \) is only
for convenience. For quadratic \( f \) the stationary point can be written down and checked
directly.

::: {#prp-mat-lagrange}
Let \( \A \) be positive definite, \( \bm H \) of size \( q\times n \) with full row rank, and
\( f(\x)=\x\T\A\x-2\bb\T\x \). The unique minimizer of \( f \) subject to \( \bm H\x=\bm h \) is
\[
\x_*=\x_0-\A^{-1}\bm H\T\bigl(\bm H\A^{-1}\bm H\T\bigr)^{-1}(\bm H\x_0-\bm h),
\qquad \x_0=\A^{-1}\bb .
\]
:::

::: {.proof}
\( \bm H\A^{-1}\bm H\T \) is positive definite by @prp-mat-pd-properties(c,e), so \( \x_* \) is
well defined, and \( \bm H\x_*=\bm h \) by direct substitution. Stationarity of the Lagrangian gives
\( 2\A\x-2\bb+2\bm H\T\bm\theta=\bzero \). Indeed \( \A\x_*-\bb=-\bm H\T\bm\theta_* \) with
\( \bm\theta_*=(\bm H\A^{-1}\bm H\T)^{-1}(\bm H\x_0-\bm h) \). For any feasible \( \x \), put \( \bm d=\x-\x_* \),
so that \( \bm H\bm d=\bzero \). Then
\[
f(\x)=f(\x_*)+2\bm d\T(\A\x_*-\bb)+\bm d\T\A\bm d=f(\x_*)-2(\bm H\bm d)\T\bm\theta_*+\bm d\T\A\bm d
=f(\x_*)+\bm d\T\A\bm d ,
\]
which exceeds \( f(\x_*) \) unless \( \bm d=\bzero \).
:::

With \( \A=\X\T\X \) and \( \bb=\X\T\y \), this is least squares subject to a linear hypothesis
\( \bm H\bbeta=\bm h \), the starting point of Chapter 11.

## Derivatives with respect to matrices and scalars

Likelihoods for normal models involve \( \log\det\bSigma \) and \( \tr(\bSigma^{-1}\bS) \). For a
scalar function \( f \) of an \( m\times n \) matrix \( \W \) with functionally independent entries,
\( \partial f/\partial\W \) is the \( m\times n \) matrix with entries \( \partial f/\partial w_{ij} \).

::: {#prp-mat-matrix-derivatives}

::: {.enumerate options="label=(\alph*)"}
1. \( \partial\tr(\A\W)/\partial\W=\A\T \).

2. If \( \det\W>0 \), then \( \partial\log\det\W/\partial\W=(\W^{-1})\T \).

3. If \( \A(t) \) is a nonsingular matrix whose entries are differentiable functions of a
           scalar \( t \), with derivative \( \dot{\A}=\partial\A/\partial t \), then
           \[
\frac{\partial\A^{-1}}{\partial t}=-\A^{-1}\dot{\A}\A^{-1},\qquad
          \frac{\partial\log\lvert\det\A\rvert}{\partial t}=\tr\bigl(\A^{-1}\dot{\A}\bigr).
\]
:::

:::

::: {.proof}
(a) \( \tr(\A\W)=\sum_{i,j}a_{ji}w_{ij} \). (c) Differentiate \( \A\A^{-1}=\I \) by the product rule to
get \( \dot{\A}\A^{-1}+\A\,\partial\A^{-1}/\partial t=\bzero \). For the determinant,
\( \A(t+h)=\A(t)\bigl(\I+h\A^{-1}\dot{\A}+o(h)\bigr) \). In the permutation expansion of
\( \det(\I+h\bm C) \) only the diagonal product contributes terms of first order in \( h \), so
\( \det(\I+h\bm C)=1+h\tr\bm C+O(h^2) \). Hence
\( \det\A(t+h)=\det\A(t)\bigl(1+h\tr(\A^{-1}\dot{\A})+o(h)\bigr) \), and taking logarithms gives the
formula. (b) Apply (c) with \( \A(t)=\W+t\bm e_i\bm e_j\T \). Then
\( \tr(\W^{-1}\bm e_i\bm e_j\T)=\bm e_j\T\W^{-1}\bm e_i \), the \( (j,i) \) entry of \( \W^{-1} \).
:::

When \( \W \) is constrained to be symmetric, \( w_{ij} \) and \( w_{ji} \) are the same variable, so
the derivatives with respect to off-diagonal entries double. Some books (Rencher and Schaalje 2008)
report formulas of the type \( 2\W^{-1}-\diag(\W^{-1}) \) for this reason. The
simplest way to avoid the issue is to differentiate along a symmetric direction \( \W+t\bm E \),
as in (c). This is the convenient route for the normal likelihoods of
Chapter 32. The listing checks the gradient formulas against central finite
differences, and the script also checks part (c) and @prp-mat-quadratic-min. The largest
discrepancy in the listing is \( 3\times 10^{-9} \), which is the size one
expects from the finite-difference approximation.

```{.python .run #cell-calculus-check-gradients}
import numpy as np
rng = np.random.default_rng(5)
p = 4
def num_grad(f, x, h=1e-6):
    """Central-difference gradient of a scalar function of a vector or matrix."""
    g = np.zeros_like(x)
    for idx in np.ndindex(x.shape):
        e = np.zeros_like(x)
        e[idx] = h
        g[idx] = (f(x + e) - f(x - e)) / (2 * h)
    return g

a = rng.normal(size=p)
A = rng.normal(size=(p, p))                         # not symmetric
x = rng.normal(size=p)
checks = {
    "a'x":        (num_grad(lambda v: a @ v, x), a),
    "x'Ax":       (num_grad(lambda v: v @ A @ v, x), (A + A.T) @ x),
}
X = rng.normal(size=(20, p))
y = rng.normal(size=20)
b = rng.normal(size=p)
grad_ls = -2 * X.T @ (y - X @ b)
checks["||y - Xb||^2"] = (num_grad(lambda v: np.sum((y - X @ v) ** 2), b), grad_ls)

S = rng.normal(size=(p, p)); S = S @ S.T + p * np.eye(p)   # positive definite
checks["log det S"] = (num_grad(lambda W: np.linalg.slogdet(W)[1], S), np.linalg.inv(S).T)
checks["tr(AS)"] = (num_grad(lambda W: np.trace(A @ W), S), A.T)
for name, (numeric, formula) in checks.items():
    print(f"{name:14s} max error {np.abs(numeric - formula).max():.1e}")
```

## Exercises

### A. Check your understanding

::: {#exr-mat-clustered}
[A1]

Observations fall in \( a \) clusters of size \( b \), and \( \bSigma=\sigma^2\,\I_a\otimes\bm E_b \), where
\( \bm E_b=(1-\rho)\I_b+\rho\bm J_b \). Find the eigenvalues, determinant and inverse of \( \bSigma \) using
@prp-mat-kronecker, @exm-mat-equicorrelation and @exm-mat-equicorrelation-eigen.
:::

::: {#exr-mat-ridge}
[A2]

Use @prp-mat-quadratic-min to show that
\( \norm{\y-\X\bb}^2+\lambda\norm{\bb}^2 \) with \( \lambda>0 \) has the unique minimizer
\( (\X\T\X+\lambda\I)^{-1}\X\T\y \), whatever the rank of \( \X \).
:::

### B. Practice

::: {#exr-mat-constrained-ls}
[B1]

Let \( \X \) have full column rank, \( \hbeta=(\X\T\X)^{-1}\X\T\y \), and let \( \bm H \) be \( q\times p \) with full row
rank. Show that the minimizer of \( \norm{\y-\X\bb}^2 \) subject to \( \bm H\bb=\bm h \) is
\[
\hbeta_H=\hbeta-(\X\T\X)^{-1}\bm H\T\bm W^{-1}(\bm H\hbeta-\bm h),\qquad
\bm W=\bm H(\X\T\X)^{-1}\bm H\T,
\]
and that the increase in the residual sum of squares is
\( (\bm H\hbeta-\bm h)\T\bm W^{-1}(\bm H\hbeta-\bm h) \).
:::

::: {.solution}
Apply @prp-mat-lagrange with \( \A=\X\T\X \) and
\( \bb=\X\T\y \). Since \( \norm{\y-\X\bb}^2=\y\T\y+f(\bb) \) and \( \x_0=\hbeta \), this gives \( \hbeta_H \). By
completing the square as in @prp-mat-quadratic-min,
\( \norm{\y-\X\bb}^2=\norm{\y-\X\hbeta}^2+(\bb-\hbeta)\T\X\T\X(\bb-\hbeta) \) for every \( \bb \). With
\( \bb-\hbeta=-(\X\T\X)^{-1}\bm H\T\bm W^{-1}(\bm H\hbeta-\bm h) \), the increase is
\[
(\bm H\hbeta-\bm h)\T\bm W^{-1}\bm H(\X\T\X)^{-1}\bm H\T\bm W^{-1}(\bm H\hbeta-\bm h)
=(\bm H\hbeta-\bm h)\T\bm W^{-1}(\bm H\hbeta-\bm h).
\]
This quadratic form is the numerator of the \( F \) statistic in Chapter 11.
:::
