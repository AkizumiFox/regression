# Smoothing splines

A P-spline starts from a basis and adds a penalty to it. The smoothing spline
does the opposite: it asks a question about functions, with no basis in sight —
among *all* twice-differentiable functions, which balances fit against roughness
best? — and the answer turns out to be a spline, by the shortest and most
satisfying derivation in the subject.

## The variational problem

Measure the roughness of a function by \( \int_a^b\{f''(t)\}^2dt \): zero exactly
for straight lines, large for a function that turns sharply, and, unlike the
number of knots, defined for every function smooth enough to have it. Let
\[
\mathcal W=\Bigl\{f:[a,b]\to\Real\ \Big|\ f'\text{ is absolutely continuous and}
\int_a^b\{f''(t)\}^2dt<\infty\Bigr\}
\]
be the set of such functions, and consider
\[
S_\lambda(f)=\sum_{i=1}^{n}\{y_i-f(x_i)\}^2+\lambda\int_a^b\{f''(t)\}^2\,dt ,
\qquad \lambda>0 .
\]{#eq-smo-variational}

Take the design points distinct and ordered, \( a\le x_1<\dots<x_n\le b \). (If
some coincide, replace each group by its mean with a weight equal to its size;
@exr-smo-ties checks that this changes nothing.)

The two extremes are clear. As \( \lambda\to0 \) the criterion is minimized by any
function through the data, so the problem is ill-posed without the penalty; as
\( \lambda\to\infty \) the penalty forces \( f''\equiv0 \) and the minimizer is the
least squares line. In between the minimizer is some compromise between an
interpolating curve and a line, and it is not obvious in advance that it is
finite-dimensional.

::: {#def-smo-natural}
[Natural cubic spline]

A **natural cubic spline** with knots \( x_1<\dots<x_n \) is a function
\( g:[a,b]\to\Real \) that is a cubic polynomial on each \( [x_i,x_{i+1}] \), has two
continuous derivatives on \( [a,b] \), and is linear on \( [a,x_1] \) and on
\( [x_n,b] \). Equivalently: on \( [x_1,x_n] \), \( g \) is a cubic spline with
interior knots \( x_2,\dots,x_{n-1} \) satisfying the **natural boundary
conditions** \( g''(x_1)=g''(x_n)=0 \), extended to \( [a,b] \) by its tangent
lines at \( x_1 \) and \( x_n \). (The two conditions are exactly what makes the
extension \( C^2 \).)
:::

## Interpolation with the least roughness

Everything follows from one integration by parts.

::: {#lem-smo-natural-interpolant}
[The natural interpolant, and why it is the smoothest]

Let \( n\ge2 \) and \( x_1<\dots<x_n \).

::: {.enumerate options="label=(\alph*)"}
1. If \( g \) is a natural cubic spline with these knots and \( h\in\mathcal W \)
   satisfies \( h(x_i)=0 \) for every \( i \), then \( \int_a^bg''h''=0 \).

2. For every \( \bz\in\Real^n \) there is exactly one natural cubic spline \( g \)
   with \( g(x_i)=z_i \); the space of natural cubic splines with these knots has
   dimension \( n \).

3. That \( g \) is the unique minimizer of \( \int_a^bf''^2 \) among all
   \( f\in\mathcal W \) with \( f(x_i)=z_i \).
:::

:::

::: {.proof}
(a) Integrating by parts,
\( \int_a^bg''h''=[g''h']_a^b-\int_a^bg'''h' \). The boundary term vanishes because
\( g \) is linear near \( a \) and \( b \), so \( g''(a)=g''(b)=0 \). On \( [a,x_1] \)
and \( [x_n,b] \) we have \( g'''=0 \), and on \( (x_i,x_{i+1}) \) the third
derivative of a cubic is a constant \( c_i \). Hence
\[
\int_a^bg'''h'=\sum_{i=1}^{n-1}c_i\int_{x_i}^{x_{i+1}}h'
=\sum_{i=1}^{n-1}c_i\{h(x_{i+1})-h(x_i)\}=0 .
\]

(b) A cubic spline on \( [x_1,x_n] \) with the \( n-2 \) interior knots
\( x_2,\dots,x_{n-1} \) lives in a space of dimension \( 4+(n-2)=n+2
\) (@thm-ply-spline-space with \( d=3 \)), and extending it linearly outside is
determined once \( g''(x_1)=g''(x_n)=0 \). Those are two linear constraints, so the
natural cubic splines form a space of dimension at least \( n \). The evaluation
map \( g\mapsto\{g(x_1),\dots,g(x_n)\} \) is linear into \( \Real^n \); if \( g \) is
in its kernel, apply (a) with \( h=g \) to get \( \int_a^bg''^2=0 \), so \( g''\equiv0 \)
and \( g \) is linear, and a line vanishing at \( n\ge2 \) distinct points is zero.
The map is therefore injective, so the dimension is at most \( n \) and hence
exactly \( n \), and the map is a bijection.

(c) Let \( f\in\mathcal W \) interpolate \( \bz \) and put \( h=f-g \), so
\( h(x_i)=0 \). By (a) the cross term vanishes, leaving
\[
\int_a^bf''^2=\int_a^bg''^2+2\int_a^bg''h''+\int_a^bh''^2
=\int_a^bg''^2+\int_a^bh''^2 ,
\]
which is at least \( \int_a^bg''^2 \), with equality iff \( h''\equiv0 \), that is
iff \( h \) is linear; and a line vanishing at \( n\ge2 \) points is zero.
:::

::: {#thm-smo-smoothing-spline}
[The smoothing spline]

Let \( \lambda>0 \) and \( n\ge2 \) distinct design points. Then \( S_\lambda \) of
@eq-smo-variational has a unique minimizer over \( \mathcal W \), and it is the
natural cubic spline with knots at \( x_1,\dots,x_n \) whose vector of values is
\[
\hat{\mathbf f}=(\I+\lambda\bP)^{-1}\y ,
\]{#eq-smo-reinsch-fit}

where \( \bP \) is the nonnegative definite matrix representing the roughness,
\( \int_a^bg''^2=\bz\T\bP\bz \), of the natural cubic spline with values \( \bz \).
:::

::: {.proof}
Let \( f\in\mathcal W \) be arbitrary and let \( g \) be the natural cubic spline
with \( g(x_i)=f(x_i) \), which exists and is unique
by @lem-smo-natural-interpolant(b). The two functions have the same residual sum of
squares, and by @lem-smo-natural-interpolant(c) the penalty of \( g \) is no
larger, strictly smaller unless \( f=g \). Hence
\( S_\lambda(g)\le S_\lambda(f) \) with equality only if \( f=g \), and it suffices
to minimize over the \( n \)-dimensional space of natural cubic splines.

Parameterize that space by \( \bz\in\Real^n \) through the bijection
of @lem-smo-natural-interpolant(b). Since \( g \) depends linearly on \( \bz \), so
does \( g'' \), and the penalty \( \int_a^bg''^2 \) is a quadratic form
\( \bz\T\bP\bz \) with \( \bP \) symmetric and nonnegative definite.
@eq-smo-variational becomes
\( \norm{\y-\bz}^2+\lambda\bz\T\bP\bz \), a strictly convex quadratic in \( \bz \)
because \( \I+\lambda\bP \) is positive definite; by @prp-mat-quadratic-min its
unique minimizer is @eq-smo-reinsch-fit.
:::

::: {.idea}
The theorem is a *representer* result: an optimization over an
infinite-dimensional space of functions has a finite-dimensional solution, with
one basis function per observation. The reason is the orthogonality in
@lem-smo-natural-interpolant(a): any direction in which \( f \) differs from a
natural spline without changing the fit costs roughness and buys nothing.
Kimeldorf and Wahba (1971) proved the general version; see the notes.
:::

## The Reinsch algorithm

The proof used \( \bP \) without exhibiting it. Doing so
turns @eq-smo-reinsch-fit, an \( n\times n \) dense solve, into a banded one.

::: {#prp-smo-reinsch}
[Band matrices for the roughness]

Let \( h_i=x_{i+1}-x_i \) for \( i=1,\dots,n-1 \), and let \( n\ge3 \). Define
\( \Q \), of size \( n\times(n-2) \), and \( \R \), symmetric tridiagonal of order
\( n-2 \), by
\[
q_{j-1,j}=\frac1{h_{j-1}},\quad
q_{j,j}=-\frac1{h_{j-1}}-\frac1{h_j},\quad
q_{j+1,j}=\frac1{h_j},
\]
\[
r_{jj}=\frac{h_{j-1}+h_j}{3},\qquad r_{j,j+1}=r_{j+1,j}=\frac{h_j}{6},
\]
for \( j=2,\dots,n-1 \) (columns and rows indexed by the interior knots), all other
entries zero. Then:

::: {.enumerate options="label=(\alph*)"}
1. \( \R \) is positive definite.

2. If \( g \) is the natural cubic spline with values \( \bz \) and
   \( \bgamma=\{g''(x_2),\dots,g''(x_{n-1})\}\T \), then \( \R\bgamma=\Q\T\bz \) and
   \( \int_a^bg''^2=\bgamma\T\R\bgamma \), so that
   \( \bP=\Q\R^{-1}\Q\T \).

3. The fit of @thm-smo-smoothing-spline satisfies
   \[
   (\R+\lambda\Q\T\Q)\,\hat{\bgamma}=\Q\T\y,
   \qquad \hat{\mathbf f}=\y-\lambda\Q\hat{\bgamma} ,
   \]{#eq-smo-reinsch}

   and \( \R+\lambda\Q\T\Q \) has bandwidth two, so the fit costs \( O(n) \)
   operations.
:::

:::

::: {.proof}
(a) \( \R \) is symmetric with positive diagonal, and
\( r_{j,j-1}+r_{j,j+1}=(h_{j-1}+h_j)/6<r_{jj} \), so it is strictly diagonally
dominant and therefore positive definite.

(b) On \( [x_j,x_{j+1}] \), \( g'' \) is linear with endpoint values
\( \gamma_j \) and \( \gamma_{j+1} \) (setting \( \gamma_1=\gamma_n=0 \)). Writing
\( u=t-x_j \) and integrating twice,
\[
g(x_j+u)=z_j+g'(x_j)u+\tfrac12\gamma_ju^2+\frac{\gamma_{j+1}-\gamma_j}{6h_j}u^3 .
\]
Setting \( u=h_j \) and solving for \( g'(x_j) \) gives
\( g'(x_j)=(z_{j+1}-z_j)/h_j-h_j(2\gamma_j+\gamma_{j+1})/6 \), and differentiating
at \( u=h_j \),
\[
g'(x_{j+1}^-)=\frac{z_{j+1}-z_j}{h_j}+\frac{h_j(\gamma_j+2\gamma_{j+1})}{6} .
\]
Equating this with \( g'(x_{j+1}^+) \), read off the interval
\( [x_{j+1},x_{j+2}] \), and writing \( i=j+1 \),
\[
\begin{aligned}
\frac{z_{i+1}-z_i}{h_i}-\frac{z_i-z_{i-1}}{h_{i-1}}
&=\frac{h_{i-1}\gamma_{i-1}+2(h_{i-1}+h_i)\gamma_i+h_i\gamma_{i+1}}{6} ,
\end{aligned}
\]
which is \( (\Q\T\bz)_i=(\R\bgamma)_i \). For the penalty,
\( \int_{x_j}^{x_{j+1}}g''^2=h_j(\gamma_j^2+\gamma_j\gamma_{j+1}+\gamma_{j+1}^2)/3 \)
by direct integration of a linear function squared, and summing over \( j \)
reproduces \( \bgamma\T\R\bgamma \) term by term. Since
\( \bgamma=\R^{-1}\Q\T\bz \), the penalty is \( \bz\T\Q\R^{-1}\Q\T\bz \).

(c) Put \( \bgamma=\R^{-1}\Q\T\hat{\mathbf f} \). From @eq-smo-reinsch-fit,
\( \hat{\mathbf f}+\lambda\Q\bgamma=\y \), and multiplying by \( \Q\T \) and using
\( \R\bgamma=\Q\T\hat{\mathbf f} \) gives
\( \R\bgamma=\Q\T\y-\lambda\Q\T\Q\bgamma \), which is @eq-smo-reinsch. \( \R \) is
tridiagonal and \( \Q \) has three nonzero entries per column in consecutive rows,
so \( \R+\lambda\Q\T\Q \) has bandwidth two. A Cholesky factor inherits the
bandwidth: in @eq-cmp-cholesky, if \( r_{ki}=0 \) whenever \( i-k>2 \) for every
\( k<j \), then for \( i-j>2 \) every product \( r_{kj}r_{ki} \) in the sum
vanishes (a nonzero \( r_{kj} \) needs \( j-k\le2 \), whence \( i-k>2 \)) and so
does \( a_{ji} \), so \( r_{ji}=0 \). The factorization of @thm-cmp-cholesky
therefore visits \( O(n) \) entries instead of \( O(n^3) \), and so do the two
multiplications by \( \Q \).
:::

Reinsch (1967) gave this algorithm before the variational theorem was widely
known; it is why smoothing splines remain the cheapest smoother for a single
covariate. Nothing in it forms \( \bS_\lambda=(\I+\lambda\bP)^{-1} \), which is
dense; the diagonal of \( \bS_\lambda \), needed for cross-validation in
[Section 43.5](05-choosing-lambda.html), can also be got in \( O(n) \) from a
banded inverse (Hutchinson and de Hoog, 1985).

::: {#exm-smo-spline-fits}
[A smoothing spline, checked three ways]

Take \( n=120 \) design points drawn uniformly on \( [0,1] \)
and the mean function of @exm-smo-bandwidth with
\( \sigma=0.25 \). At \( \lambda=1\times 10^{-4} \)
the fit uses \( 12.32 \) effective degrees of freedom.
Three independent computations agree with it. The roughness
\( \bgamma\T\R\bgamma \) of the natural spline interpolating the true curve,
\( 1069.5 \), matches \( \int_a^bg''^2 \) from a
separately constructed natural cubic interpolant to a relative error of
\( 1.6\times 10^{-12} \); the banded solve
@eq-smo-reinsch agrees with the dense form @eq-smo-reinsch-fit to
\( 1.3\times 10^{-10} \); and minimizing
\( \norm{\y-\B\bgamma}^2+\lambda\bgamma\T\boldsymbol{\Omega}\bgamma \) over the
*whole* space of cubic splines with a knot at every design point — dimension
\( 122 \), with
\( \boldsymbol{\Omega}_{jk}=\int B_j''B_k'' \) computed exactly — gives fitted
values differing from the smoothing spline by
\( 4.4\times 10^{-10} \). That is
@thm-smo-smoothing-spline in numerical form: the minimizer over the larger space
is the natural spline.
:::

::: {when-format="html"}
![**Figure 43.4.1.** (a) Smoothing splines at three values of \( \lambda \),
labelled by effective degrees of freedom. (b) Three rows of \( \bS_\lambda \) at
df \( =9 \), the weights producing \( \hat f(x) \) at \( x=0.02 \),
\( 0.29 \) and \( 0.62 \): local, summing to one, slightly negative beside the
peak.](smoothing_spline.svg){#fig-smo-spline width=100%}
:::

::: {when-format="pdf"}
![(a) Smoothing splines at three values of \( \lambda \),
labelled by effective degrees of freedom. (b) Three rows of \( \bS_\lambda \) at
df \( =9 \), the weights producing \( \hat f(x) \) at \( x=0.02 \),
\( 0.29 \) and \( 0.62 \): local, summing to one, slightly negative beside the
peak.](smoothing_spline.pdf){width=100%}
:::

```{.python .run #cell-smoothing-spline-reinsch}
import numpy as np


def f_true(x):
    """The mean function of Section 43.1."""
    return 2 * x + np.exp(-25 * (x - 0.35) ** 2) - 0.75 * np.exp(-50 * (x - 0.80) ** 2)


def reinsch(t):
    """The band matrices Q (N x (N-2)) and R ((N-2) x (N-2)) of the Reinsch form."""
    N = len(t)
    h = np.diff(t)
    Q = np.zeros((N, N - 2))
    R = np.zeros((N - 2, N - 2))
    for j in range(N - 2):                       # column j belongs to interior knot t[j+1]
        Q[j, j], Q[j + 1, j], Q[j + 2, j] = 1 / h[j], -1 / h[j] - 1 / h[j + 1], 1 / h[j + 1]
        R[j, j] = (h[j] + h[j + 1]) / 3
        if j + 1 < N - 2:
            R[j, j + 1] = R[j + 1, j] = h[j + 1] / 6
    return Q, R


def smoothing_spline(t, y, lam):
    """Fitted values of the smoothing spline, by solving a pentadiagonal system."""
    Q, R = reinsch(t)
    gamma = np.linalg.solve(R + lam * Q.T @ Q, Q.T @ y)   # the second derivatives
    return y - lam * Q @ gamma, gamma


def smoother_matrix(t, lam):
    """The matrix S_lambda with fitted values S_lambda y (formed only to look at it)."""
    Q, R = reinsch(t)
    return np.eye(len(t)) - lam * Q @ np.linalg.solve(R + lam * Q.T @ Q, Q.T)


n, sigma = 120, 0.25
rng = np.random.default_rng(430430)
t = np.sort(rng.uniform(0, 1, n))
y = f_true(t) + sigma * rng.normal(size=n)

lam = 1e-4
fit, _ = smoothing_spline(t, y, lam)
S = smoother_matrix(t, lam)
df = np.trace(S)
print("effective degrees of freedom at lambda = 1e-4:", round(df, 2))
print("rows of S sum to one, largest departure:", np.max(np.abs(S.sum(axis=1) - 1)))
```

## Variations, and the relation to penalized splines

**Other orders.** Replacing \( f'' \) by \( f^{(m)} \) in @eq-smo-variational gives a
**natural spline of degree \( 2m-1 \)**: the proof of
@lem-smo-natural-interpolant goes through with \( m \) integrations by parts and
the boundary conditions \( f^{(j)}(x_1)=f^{(j)}(x_n)=0 \) for
\( j=m,\dots,2m-2 \). The case \( m=1 \) gives a piecewise linear fit,
\( m=2 \) the cubic one above, \( m=3 \) a quintic. In practice \( m=2 \) is used
almost always; the extra smoothness of \( m=3 \) rarely shows.

**Reproducing kernels.** The penalty \( \int f''^2 \) is the squared norm of the
part of \( f \) orthogonal to the straight lines, in a Hilbert space of functions
in which evaluation \( f\mapsto f(x) \) is a bounded linear functional and so
representable as an inner product with a function \( R_x \), the *reproducing
kernel*. In that language @thm-smo-smoothing-spline is a special case of the
representer theorem of Kimeldorf and Wahba (1971): the minimizer of
\( \sum_i\{y_i-f(x_i)\}^2+\lambda\norm{Pf}^2 \) lies in the span of
\( R_{x_1},\dots,R_{x_n} \) plus the null space of \( P \). That is the route of
Wahba (1990), and it covers thin-plate splines, splines on a sphere and the
kernel methods of machine learning. This book does not develop it; the proof
above is self-contained for a single covariate.

**Smoothing splines and P-splines.** A smoothing spline is a penalized spline
with a knot at every distinct design point and the exact derivative penalty, as
@exm-smo-spline-fits confirms numerically, and a P-spline with \( K=20 \) knots is
a *low-rank* approximation to it: both fits are determined, through
@eq-smo-smoother, by how the penalty ranks a handful of smooth directions, so the
two agree once \( m \) comfortably exceeds the effective degrees of freedom. The
smoothing spline has no knot choice and costs \( O(n) \); the P-spline carries an
explicit coefficient vector and generalizes to several covariates and to
nonnormal responses, which is why
[Chapter 44](../ch44-additive-models/index.html) uses it.

::: {.remark}
[The equivalent kernel of a smoothing spline]

Panel (b) of [Figure 43.4.1](04-smoothing-splines.html#fig-smo-spline) shows the
rows of \( \bS_\lambda \), the smoothing spline's equivalent kernel in the sense
of [Section 43.1](01-kernels.html). Silverman (1984) proved that at an interior
point, for large \( n \) and small \( \lambda \), the weight given to a design
point at distance \( u \) is approximately
\( \{n b(x)g(x)\}^{-1}\,\tfrac12\exp(-\lvert v\rvert/\sqrt2)
\sin(\lvert v\rvert/\sqrt2+\pi/4) \) with \( v=u/b(x) \) and local bandwidth
\( b(x)=\{\lambda/(ng(x))\}^{1/4} \). Two features of that formula show in the
picture: the weights decay exponentially, and the sine makes them change sign, so
a smoothing spline is a kernel smoother with a small negative side lobe — here
the most negative weight is \( -0.0034 \) against a
peak of \( 0.0778 \). The bandwidth adapts to the
design density by itself, as \( g(x)^{-1/4} \). We state this without proof; the
argument is a Green's-function calculation.
:::

## Exercises

### A. Check your understanding

::: {#exr-smo-limits}
[A1]

Show directly from @eq-smo-reinsch that \( \hat{\mathbf f}\to\y \) as
\( \lambda\to0 \) and that \( \hat{\mathbf f} \) tends to the least squares line as
\( \lambda\to\infty \). (For the second, use @prp-smo-reinsch(b): what does
\( \Q\T\bz=\bzero \) say about \( \bz \)?)
:::

::: {.solution}
As \( \lambda\to0 \), @eq-smo-reinsch gives \( \hat{\bgamma}\to\R^{-1}\Q\T\y \),
bounded, so \( \hat{\mathbf f}=\y-\lambda\Q\hat{\bgamma}\to\y \). As
\( \lambda\to\infty \), the penalty \( \bz\T\bP\bz \) must tend to zero, and
\( \bP=\Q\R^{-1}\Q\T \) with \( \R \) positive definite, so \( \Q\T\bz\to\bzero \).
Now \( (\Q\T\bz)_i \) is a second divided difference of \( \bz \) at
\( x_{i-1},x_i,x_{i+1} \), and it vanishes for all \( i \) exactly when the points
\( (x_i,z_i) \) are collinear. So the limit lies in the two-dimensional space
\( \C([\bone,\x]) \), and since the residual sum of squares is still being
minimized there, it is the least squares line.
:::

### B. Practice

::: {#exr-smo-ties}
[B1]

Suppose the design has ties: the distinct values are \( t_1<\dots<t_N \) and
\( n_k \) observations share \( t_k \), with mean \( \bar y_k \). Show that
\( \sum_i\{y_i-f(x_i)\}^2=\sum_kn_k\{\bar y_k-f(t_k)\}^2+C \) with \( C \) free of
\( f \), so that @thm-smo-smoothing-spline holds with knots at the distinct values
and \( \I \) replaced by \( \diag(n_1,\dots,n_N) \) in @eq-smo-reinsch-fit.
:::

::: {#exr-smo-reinsch-check}
[B2]

Using the cell above, build \( \Q \) and \( \R \) for \( n=8 \) equally spaced
points on \( [0,1] \) and verify by hand that \( \R \) is tridiagonal with
\( r_{jj}=2h/3 \) and \( r_{j,j+1}=h/6 \), and that \( \Q\T\bone=\bzero \) and
\( \Q\T\x=\bzero \). What does the second pair of identities say about
\( \bS_\lambda\bone \) and \( \bS_\lambda\x \)?
:::

### C. Going deeper

::: {#exr-smo-higher-m}
[C1]

State and prove the analogue of @lem-smo-natural-interpolant(a) for the penalty
\( \int_a^b\{f^{(m)}\}^2 \), with \( g \) a natural spline of degree \( 2m-1 \).
Which boundary conditions are needed to kill the boundary terms after \( m \)
integrations by parts, and what replaces "a line vanishing at \( n\ge2 \) points is
zero"?
:::

::: {#exr-smo-interpolation-limit}
[C2]

Show that if \( \lambda>0 \) then \( \hat f \) never interpolates the data unless
the data already lie on a straight line. (Use @eq-smo-reinsch-fit and the fact
that \( \bP\bz=\bzero \) iff \( \bz\in\C([\bone,\x]) \), which
@exr-smo-limits establishes.) Contrast this with a regression spline with
\( K=n-4 \) knots, and say what the difference means for the leave-one-out
shortcut of @thm-sel-loocv.
:::

::: {.solution}
\( \hat{\mathbf f}=\y \) requires \( (\I+\lambda\bP)\y=\y \), that is
\( \bP\y=\bzero \), that is \( \y \) in \( \C([\bone,\x]) \). So for
\( \lambda>0 \) the fit leaves a nonzero residual unless the data are exactly
linear, and every \( s_{ii} \) is strictly less than one, which is what makes the
shortcut of @thm-sel-loocv well defined. An unpenalized regression spline with
\( n \) basis functions interpolates, every \( s_{ii}=1 \), the shortcut divides
by zero, and cross-validation is useless — one more reason the dial should be a
penalty and not a knot count.
:::
