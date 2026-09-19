# Why not invert the cross-product matrix

[Chapter 5](../ch05-model-and-least-squares/index.html) found the least squares estimate by setting a gradient to zero, which gave
\( \hbeta=(\X\T\X)^{-1}\X\T\y \) when \( \X \) has full column rank (@thm-lm-ls-full-rank).
Read as a recipe, the formula says: form the \( p\times p \) matrix \( \X\T\X \), invert it,
and multiply. [Section 6.10](../ch06-projections/10-computation.html) already showed that the
recipe can lose twice as many digits as necessary. This section explains why.

## Floating-point numbers

A double-precision number is zero or a number \( \pm m\times2^{e-52} \), where \( m \) is an
integer with \( 2^{52}\le m<2^{53} \) and \( e \) is an integer in a fixed range
(\( -1022\le e\le1023 \)). Between \( 2^e \) and \( 2^{e+1} \) the doubles are therefore equally
spaced, \( 2^{e-52} \) apart. Just above \( 1 \) the spacing is
\( \epsilon_M=2^{-52}\approx2.22\times 10^{-16} \), the **machine epsilon**. Rounding a real
number \( x \) in the range of the doubles to the nearest double gives a number
\( \operatorname{fl}(x) \) with
\[
\operatorname{fl}(x)=x(1+\delta),\qquad \lvert\delta\rvert\le u=2^{-53}\approx1.11\times 10^{-16}.
\]
The **unit roundoff** \( u=\epsilon_M/2 \) is half the gap, because the nearest double is at
most half a gap away. Rounding errors are *relative*: a double represents large and small
numbers with the same number of significant digits, about sixteen.

The IEEE 754 standard, which essentially all current hardware follows, requires the basic
operations to be correctly rounded. The computed result of \( x+y \) is the exact sum,
rounded. That is the only property of the hardware that the analysis below uses.

::: {#def-cmp-standard-model}
[Standard model of floating-point arithmetic]

For floating-point numbers \( x,y \) and \( \circ\in\{+,-,\times,/\} \), the computed result
satisfies
\[
\operatorname{fl}(x\circ y)=(x\circ y)(1+\delta),\qquad\lvert\delta\rvert\le u ,
\]
and \( \operatorname{fl}(\sqrt x)=\sqrt x\,(1+\delta) \) with \( \lvert\delta\rvert\le u \),
provided no result overflows or underflows.
:::

Two consequences are worth seeing at once. The number \( 1+u \) lies exactly halfway
between the doubles \( 1 \) and \( 1+\epsilon_M \), and the tie is broken towards \( 1 \), so
adding \( u \) to \( 1 \) changes nothing. And most decimal fractions are not doubles, so
familiar identities fail at the level of \( u \):

```{.python .run #cell-floating-point-unit}
import numpy as np

eps = np.finfo(float).eps          # spacing of the doubles just above 1: 2^-52
u = eps / 2                        # unit roundoff: 2^-53
print("eps =", eps, "  u =", u)
print("1 + u == 1:", 1.0 + u == 1.0, "   1 + eps == 1:", 1.0 + eps == 1.0)
print("0.1 + 0.2 == 0.3:", 0.1 + 0.2 == 0.3)
```

Trouble comes from accumulation and above all from subtraction. The difference of nearly equal numbers
is computed exactly but exposes the errors they already carried: if \( a \) and \( b \) agree in ten digits, a
rounding error in either becomes a large relative error in \( a-b \). This **cancellation** is behind every
numerical difficulty in this chapter.

## Accumulating rounding errors

Products of many factors \( 1+\delta_i \) appear whenever errors accumulate. The following
lemma packages them.

::: {#lem-cmp-gamma}
Let \( \lvert\delta_i\rvert\le u \) and \( \rho_i\in\{1,-1\} \) for \( i=1,\dots,k \), and suppose
\( ku<1 \). Then
\[
\prod_{i=1}^k(1+\delta_i)^{\rho_i}=1+\theta_k,\qquad
\lvert\theta_k\rvert\le\gamma_k=\frac{ku}{1-ku}.
\]{#eq-cmp-gamma}

:::

::: {.proof}
Each factor lies between \( 1-u \) and \( 1/(1-u) \): for \( \rho_i=1 \) because
\( 1+u\le1/(1-u) \), and for \( \rho_i=-1 \) because \( 1/(1+\delta_i)\ge1/(1+u)\ge1-u \). So the
product lies between \( (1-u)^k \) and \( (1-u)^{-k} \). By Bernoulli's inequality
\( (1-u)^k\ge1-ku \), which gives the lower end \( 1-ku\ge1-\gamma_k \), and also the upper end
\( (1-u)^{-k}\le1/(1-ku)=1+\gamma_k \).
:::

The constant \( \gamma_k \) is essentially \( ku \): for \( k=10^6 \) it is about
\( 1.1\times 10^{-10} \). This is a worst case; in practice errors of different signs partly
cancel and grow more like \( \sqrt k\,u \).

The basic step of every least squares algorithm is an inner product.

::: {#prp-cmp-inner-product}
[Rounding error in an inner product]

Let \( \x,\y\in\Real^n \) have floating-point entries and \( nu<1 \). The inner product computed
by recursive summation, \( s_1=\operatorname{fl}(x_1y_1) \) and
\( s_k=\operatorname{fl}(s_{k-1}+\operatorname{fl}(x_ky_k)) \), satisfies
\[
s_n=\sum_{i=1}^nx_iy_i(1+\theta_i),\qquad\lvert\theta_i\rvert\le\gamma_n .
\]
Consequently \( \lvert s_n-\x\T\y\rvert\le\gamma_n\sum_i\lvert x_iy_i\rvert \), and
\( s_n=(\x+\Delta\x)\T\y \) for a vector \( \Delta\x \) with
\( \lvert\Delta x_i\rvert\le\gamma_n\lvert x_i\rvert \).
:::

::: {.proof}
By the standard model, \( \operatorname{fl}(x_ky_k)=x_ky_k(1+\mu_k) \) and
\( s_k=\bigl(s_{k-1}+x_ky_k(1+\mu_k)\bigr)(1+\sigma_k) \) for \( k\ge2 \), with all
\( \lvert\mu_k\rvert,\lvert\sigma_k\rvert\le u \). Unrolling the recursion,
\[
s_n=\sum_{i=1}^nx_iy_i(1+\mu_i)\prod_{k=\max(i,2)}^n(1+\sigma_k).
\]
The \( i \)th term carries at most \( n \) factors of the form \( 1+\delta \), so by
@lem-cmp-gamma it equals \( x_iy_i(1+\theta_i) \) with \( \lvert\theta_i\rvert\le\gamma_n \)
(since \( \gamma_k\le\gamma_n \) for \( k\le n \)). The two consequences follow by taking
\( \Delta x_i=\theta_ix_i \).
:::

The first conclusion is a *forward* bound, which can be large relative to a result whose terms cancel.
The second is a *backward* statement: the computed number is the exact inner product of slightly
perturbed data. Over 200 random inner
products of length up to 400, checked in rational arithmetic, the largest ratio of error to bound is
\( 0.149 \).

```{.python .run #cell-floating-point-inner}
from fractions import Fraction

def gamma(k, u=2.0 ** -53):
    return k * u / (1 - k * u)

def dot_recursive(x, y):
    """Inner product by recursive summation, one rounding per operation."""
    s = 0.0
    for a, b in zip(x, y):
        s = s + a * b
    return s

rng = np.random.default_rng(101)
worst = 0.0
for trial in range(200):
    n = int(rng.integers(2, 400))
    x = rng.normal(size=n) * 10.0 ** rng.integers(-3, 4, size=n)
    y = rng.normal(size=n)
    exact = sum(Fraction(a) * Fraction(b) for a, b in zip(x, y))       # exact rational value
    error = abs(Fraction(dot_recursive(x, y)) - exact)
    bound = Fraction(gamma(n)) * sum(abs(Fraction(a) * Fraction(b)) for a, b in zip(x, y))
    worst = max(worst, float(error / bound))
print(f"largest error / bound over 200 inner products: {worst:.3f}")
```

## Conditioning and stability

A computation turns data \( d \) into a result \( f(d) \). The computed result can be inaccurate
because the *problem* is sensitive, or because the *algorithm* adds more error than the data
justify. The first is measured by a **condition number**, the largest factor by which a small
relative change in \( d \) can be magnified in \( f(d) \); the second by the following standard.

::: {#def-cmp-backward-stable}
[Backward stability]

An algorithm that computes \( f(d) \) is **backward stable** if, for every input \( d \), its
computed output \( \hat f(d) \) is the exact result for nearby data:
\( \hat f(d)=f(d+\Delta d) \) with \( \norm{\Delta d}\le c\,u\norm{d} \), where \( c \) is a
modest constant depending only on the dimensions of the problem.
:::

@prp-cmp-inner-product says that recursive summation is a backward stable way to
compute inner products, with \( c=\gamma_n/u\approx n \). Backward stability is the most one can ask, because storing the data as doubles already perturbs
them by \( u \). The two notions combine, to first order in \( u \), as
\[
\text{relative forward error}\;\lesssim\;\text{condition number}\times\text{backward error}.
\]{#eq-cmp-rule}

For least squares the condition number, which involves \( \kappa(\X) \) and sometimes \( \kappa(\X)^2 \), is the
subject of [Section 10.5](05-conditioning.html). Sections 10.2 to 10.4 ask whether algorithms are backward
stable *for the least squares problem*, with perturbations of \( \X \) and \( \y \).

## What forming the cross-product matrix loses

The normal-equations recipe is not backward stable for least squares, and the damage is
done in its very first step. Here is what rounding does to \( \X\T\X \).

::: {#prp-cmp-gram-rounding}
[Rounding in the cross-product matrix]

Let \( \X \) be \( n\times p \) with full column rank and let \( \hat{\A} \) be \( \X\T\X \) computed
entry by entry as in @prp-cmp-inner-product. Then
\( \hat{\A}=\X\T\X+\Delta\A \) with \( \lvert\Delta\A\rvert\le\gamma_n\lvert\X\rvert\T\lvert\X\rvert \)
entrywise, and
\[
\norm{\Delta\A}_2\le p\,\gamma_n\norm{\X}_2^2 .
\]
In particular \( \hat{\A} \) is positive definite whenever \( p\,\gamma_n\,\kappa(\X)^2<1 \).
:::

::: {.proof}
The entrywise bound is @prp-cmp-inner-product applied to each pair of columns.
The spectral norm is at most the Frobenius norm, the Frobenius norm is monotone in the
absolute values of the entries and submultiplicative, and
\( \norm{\X}_F^2=\sum_i\sigma_i^2\le p\,\sigma_1^2 \) by @prp-mat-svd-norms(a). Hence
\[
\norm{\Delta\A}_2\le\gamma_n\bigl\lVert\,\lvert\X\rvert\T\lvert\X\rvert\,\bigr\rVert_F
\le\gamma_n\norm{\X}_F^2\le p\,\gamma_n\sigma_1^2 .
\]
For a unit vector \( \bv \),
\( \bv\T\hat{\A}\bv=\norm{\X\bv}^2+\bv\T\Delta\A\bv\ge\sigma_p^2-p\,\gamma_n\sigma_1^2 \),
which is positive when \( p\,\gamma_n\sigma_1^2/\sigma_p^2<1 \).
:::

The information about \( \X \) that sits in its smallest singular value \( \sigma_p \) appears in
\( \X\T\X \) as an eigenvalue \( \sigma_p^2 \). Rounding perturbs \( \X\T\X \) by an amount of order
\( u\sigma_1^2 \). Once \( \sigma_p^2 \) is below that level, the computed cross-product matrix may
not determine the direction \( \bv_p \) at all. Even rounding the *exact* \( \X\T\X \) to doubles
does this much damage. The loss becomes total around \( \kappa(\X)\approx u^{-1/2}=9.5\times 10^{7} \), and it
is already serious well before that: the normal equations are then a problem with condition
number \( \kappa(\X\T\X)=\kappa(\X)^2 \) (@eq-mat-kappa-squared), and
@eq-cmp-rule predicts relative errors of order \( \kappa(\X)^2u \) in \( \hbeta \) whatever
method is then used to solve them.

::: {#exm-cmp-lauchli}
[A cross-product matrix that rounds to a singular one]

Take \( \delta=10^{-9} \) and
\[
\X=\begin{pmatrix}1&1\\\delta&0\\0&\delta\end{pmatrix},\qquad
\X\T\X=\begin{pmatrix}1+\delta^2&1\\1&1+\delta^2\end{pmatrix}.
\]
The columns of \( \X \) are independent, and \( \kappa(\X)=\sqrt{2+\delta^2}/\delta=1.41\times 10^{9} \).
But \( \delta^2=10^{-18} \) is smaller than \( u \), so \( 1+\delta^2 \) rounds to \( 1 \) and the
computed \( \X\T\X \) is the singular matrix of ones. The Cholesky factorization of
[Section 10.2](02-cholesky.html) fails, and no method applied to the computed normal
equations can recover the solution. For \( \y=\X(1,1)\T \), the exact least squares estimate is
\( (1,1)\T \) with zero residual. QR, which never forms \( \X\T\X \), returns it with relative
error \( 3.9\times 10^{-16} \).
:::

```{.python .run #cell-floating-point-lauchli}
delta = 1e-9
X = np.array([[1.0, 1.0],
              [delta, 0.0],
              [0.0, delta]])
y = X @ np.array([1.0, 1.0])            # exact data: the solution is b = (1, 1), residual 0

A = X.T @ X                             # 1 + delta^2 rounds to 1
print("computed X'X =\n", A)
print("cond(X) =", np.linalg.cond(X))
try:
    np.linalg.cholesky(A)
except np.linalg.LinAlgError as err:
    print("Cholesky of the computed X'X fails:", err)

Q, R = np.linalg.qr(X)
b_qr = np.linalg.solve(R, Q.T @ y)
print("QR solution:", b_qr)
```

## Why not the inverse either

Even when \( \X\T\X \) is safely positive definite, the inverse costs \( p^3 \) flops against \( p^3/3 \) for a
Cholesky factorization plus \( 2p^2 \) for two triangular solves, and it is less accurate. Statistics needs
only the diagonal of \( \Cov(\hbeta)=\sigma^2(\X\T\X)^{-1} \) (@thm-lm-moments): if \( \X\T\X=\R\T\R \), its
\( j \)th entry is \( \sigma^2 \) times the squared length of row \( j \) of \( \R^{-1} \).

## Counting the work

Cost is traditionally counted in floating-point operations (**flops**), each addition,
subtraction, multiplication and division counting one. Flop counts ignore memory traffic and
parallelism, so they are only a first guide.

| Task | Leading flop count | Where |
|---|---|---|
| Form \( \X\T\X \) and \( \X\T\y \) (using symmetry) | \( np^2 \) | this section |
| Cholesky factorization of \( \X\T\X \) | \( p^3/3 \) | [Section 10.2](02-cholesky.html) |
| Two triangular solves | \( 2p^2 \) | [Section 10.2](02-cholesky.html) |
| Householder QR of \( \X \), with \( \Q\T\y \) | \( 2np^2-2p^3/3 \) | [Section 10.3](03-householder-givens.html) |
| Givens QR of \( \X \) | \( 3np^2-p^3 \) | [Section 10.3](03-householder-givens.html) |
| Modified Gram–Schmidt | \( 2np^2 \) | [Section 10.3](03-householder-givens.html) |
| Singular values and \( \V \), QR first (R-SVD) | \( 2np^2+11p^3 \) | [Section 10.4](04-svd.html) |

When \( n\gg p \), as usual in regression, the normal equations cost about \( np^2 \) and Householder QR about
\( 2np^2 \): the accurate method is at most twice as expensive. And the one-pass, \( O(p^2) \)-memory
accumulation of cross-products has a QR counterpart in row-wise Givens QR
([Section 10.3](03-householder-givens.html)).

::: {.idea}
Forming \( \X\T\X \) squares the condition number before any equation is solved, and
inverting it adds cost and error. Factor \( \X \) itself, by an orthogonal transformation, and
read everything off the factors.
:::

## Exercises

### A. Check your understanding

::: {#exr-cmp-doubles}
[A1]

Show that the doubles in \( [1,2) \) are exactly the numbers \( 1+k\,2^{-52} \), \( k=0,1,\dots,2^{52}-1 \).
What is the smallest double greater than \( 1 \)? What is the largest integer \( N \) such that
every integer from \( 0 \) to \( N \) is a double?
:::

::: {#exr-cmp-cancellation}
[A2]

For \( x=10^{-10} \), the computer evaluates \( (1+x)-1 \). Use the standard model to show that the
result is \( x+\delta(1+x) \) for some \( \lvert\delta\rvert\le u \), where the subtraction itself
is exact. Estimate the relative error, and check it numerically. Which of the two
operations is responsible?
:::

::: {.solution}
The addition gives \( \operatorname{fl}(1+x)=(1+x)(1+\delta) \). This number is within a factor of
two of \( 1 \), and the difference of two doubles within a factor of two of each other is itself a
double, so the subtraction is exact:
\( \operatorname{fl}(1+x)-1=x+\delta(1+x) \). The relative error is
\( \lvert\delta\rvert(1+x)/x\approx u/x\approx10^{-6} \). All of it comes from rounding \( 1+x \),
which kept only about six significant digits of \( x \). The subtraction then exposed that loss.
:::

### B. Practice

::: {#exr-cmp-summation}
[B1]

Show that recursive summation of \( x_1,\dots,x_n \) returns
\( \sum_ix_i(1+\theta_i) \) with \( \lvert\theta_i\rvert\le\gamma_{n-1} \), and that pairwise
summation returns such a sum with \( \lvert\theta_i\rvert\le\gamma_{\lceil\log_2n\rceil} \).
:::

::: {.solution}
With \( s_1=x_1 \) and \( s_k=(s_{k-1}+x_k)(1+\sigma_k) \),
\( s_n=\sum_ix_i\prod_{k=\max(i,2)}^n(1+\sigma_k) \). Term \( i \) has at most \( n-1 \) factors, so
@lem-cmp-gamma gives \( \lvert\theta_i\rvert\le\gamma_{n-1} \). In pairwise summation each
\( x_i \) passes through one addition per level of the tree, and a balanced tree over \( n \) leaves
has \( \lceil\log_2n\rceil \) levels. The same unrolling gives at most
\( \lceil\log_2n\rceil \) factors per term.
:::

::: {#exr-cmp-matvec}
[B2]

Let \( \A \) be \( m\times n \) and \( \bb\in\Real^n \). Show that the computed product
\( \operatorname{fl}(\A\bb) \), each entry an inner product, equals \( (\A+\Delta\A)\bb \) with
\( \lvert\Delta\A\rvert\le\gamma_n\lvert\A\rvert \). Deduce
\( \norm{\operatorname{fl}(\A\bb)-\A\bb}\le\gamma_n\sqrt n\,\norm{\A}_2\norm{\bb} \).
:::

### C. Going deeper

::: {#exr-cmp-lauchli-family}
[C1]

For \( 0<\delta\le2^{-27} \) let \( \X_\delta \) be the matrix of @exm-cmp-lauchli, and let
\( \y=\X_\delta(3,1)\T \) and \( \y'=\X_\delta(1,3)\T \). Show that the computed cross-product
matrices of \( [\X_\delta,\y] \) and \( [\X_\delta,\y'] \) are identical (and do not depend on
\( \delta \)), although the exact least squares estimates are \( (3,1)\T \) and \( (1,3)\T \).
Conclude that no algorithm that starts from the computed cross-product matrix can solve
these problems. Is anything special about \( 2^{-27} \)?
:::
