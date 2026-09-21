# Polynomial regression and its instability

The simplest way to let a regression curve bend is to add powers of the
regressor as columns. Nothing in Parts II and III changes, because nothing there
ever asked what the columns *were*.

::: {#def-ply-polynomial}
[Polynomial regression]

Let \( x_1,\dots,x_n \) be observed values of a single regressor. The
**polynomial regression model of degree \( d \)** is
\[
\E(Y_i)=\beta_0+\beta_1x_i+\beta_2x_i^2+\dots+\beta_dx_i^d,
\qquad i=1,\dots,n,
\]
a linear model with the \( n\times(d+1) \) **Vandermonde** matrix
\( \X=[\bone,\x,\x^2,\dots,\x^d] \), where \( \x^j \) has entries \( x_i^j \). More
generally, any functions \( \phi_0,\dots,\phi_d \) give columns \( \phi_j(\x) \) and
a fitted curve \( \hat f(t)=\sum_j\hat{\beta}_j\phi_j(t) \); they are a **basis**,
and the model is the space \( \mathcal P_d=\spn\{\phi_0,\dots,\phi_d\} \), not the
list.
:::

## Everything from Parts II and III applies

Write \( \X \) for the Vandermonde matrix. Then \( \hbeta \) is unbiased as soon
as the \( x_i \) take \( d+1 \) distinct values (@exr-ply-vandermonde-rank); the
fitted curve \( \hat f(t)=\mathbf{x}(t)\T\hbeta \), with
\( \mathbf{x}(t)=(1,t,\dots,t^d)\T \), is estimable at every \( t \), so
@thm-ci-estimable-interval gives it an interval and @thm-ci-working-hotelling a
band; raising the degree is a nested comparison (@thm-glh-f-test); replicated
\( x \) values give the lack-of-fit test of @thm-cor-lack-of-fit against the pure
error of @def-cor-pure-error; and the diagnostics of
[Chapter 20](../ch20-residuals-leverage-influence/index.html) are untouched. By
@thm-proj-reparam none of this depends on the basis, only on
\( \mathcal P_d \).

::: {.warning}
What does *not* carry over is the reading of a coefficient: there is no way to
change \( x \) while holding \( x^2 \) fixed, so nothing beyond the intercept has
the "other regressors held fixed" meaning of @thm-proj-fwl. Only the curve and
its contrasts \( f(t_2)-f(t_1) \) mean anything, and whether those describe an
intervention is the question of
[Chapter 25](../ch25-causal-interpretation/index.html).
:::

## What the degree buys

Raising the degree adds a column, so bias falls and variance rises: the
trade-off of @exm-dep-polynomial-risk. The useful question is how the two move
in practice.

::: {#exm-ply-co2-degree}
[Degree in the carbon dioxide series]

The monthly mean carbon dioxide at Mauna Loa, March 1958 to December 2001,
gives \( n=521 \) observations; @exm-lm-co2 modelled it with a
quadratic trend and two seasonal harmonics. Keep the four harmonic columns and
replace the quadratic by a polynomial of degree \( d \) in time rescaled to
\( [-1,1] \). The residual standard deviation falls from
\( 0.7602 \) ppm at \( d=2 \) (\( p=7 \)) to
\( 0.4874 \) at \( d=8 \) and \( 0.4300 \) at
\( d=16 \), while no harmonic coefficient moves by more than
\( 0.0150 \) ppm.

The improvement is real: panel (a) of [Figure 42.1.1](#fig-ply-degree) shows the
residuals of the quadratic fit wandering by a part per million over decades, and
the higher-degree trends tracking that wandering. The price is in panel (b). In
the middle of the range the standard error of the fitted mean is much the same
at all three degrees — \( 0.0830 \),
\( 0.0676 \), \( 0.0730 \) ppm at 1980 —
but at the left-hand edge it is \( 0.1217 \),
\( 0.2011 \), \( 0.2899 \), so the ratio
of edge to centre rises from \( 1.47 \) to
\( 2.97 \) to \( 3.97 \). Leverages tell
the same story: average \( 0.0134 \),
\( 0.0250 \), \( 0.0403 \) against
maximum \( 0.0256 \), \( 0.1702 \),
\( 0.4546 \). At degree 16 one month, the first of the series,
carries almost half of its own fitted value.
:::

::: {when-format="html"}
![**Figure 42.1.1.** Polynomial trends for the Mauna Loa series, seasonal columns fixed.
(a) Residuals of the quadratic fit, with the degree-8 and degree-16 trends measured from
it. (b) Standard error of the fitted mean.](polynomial_degree.svg){#fig-ply-degree width=100%}
:::

::: {when-format="pdf"}
![Polynomial trends for the Mauna Loa series, seasonal columns fixed.
(a) Residuals of the quadratic fit, with the degree-8 and degree-16 trends measured from
it. (b) Standard error of the fitted mean.](polynomial_degree.pdf){width=100%}
:::

```{.python .run #cell-polynomial-degree-setup}
import numpy as np
import statsmodels.api as sm

co2 = sm.datasets.co2.load_pandas().data["co2"]
monthly = co2.resample("MS").mean().dropna()         # weekly values -> monthly means
y = monthly.to_numpy()
t = (monthly.index.year + (monthly.index.month - 0.5) / 12).to_numpy()
n = len(y)


def harmonics(s):
    """The seasonal columns of Example 5.3, held fixed throughout the chapter."""
    return np.column_stack([np.cos(2 * np.pi * s), np.sin(2 * np.pi * s),
                            np.cos(4 * np.pi * s), np.sin(4 * np.pi * s)])


lo, hi = t.min(), t.max()
u = lambda s: 2 * (s - lo) / (hi - lo) - 1           # time rescaled to [-1, 1]

print(f"n = {n}, time from {lo:.2f} to {hi:.2f}")
```

```{.python .run #cell-polynomial-degree-degree}
def trend_design(s, d):
    """Model matrix: powers 1, u, ..., u^d of rescaled time, then the four harmonics."""
    return np.column_stack([np.vander(u(s), d + 1, increasing=True), harmonics(s)])


for d in (2, 8, 16):
    X = trend_design(t, d)
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    resid = y - X @ beta
    p = X.shape[1]
    lev = np.einsum("ij,ij->i", X @ np.linalg.pinv(X.T @ X), X)   # leverages h_ii
    print(f"degree {d:2d}:  p = {p:2d}   sigma_hat = {np.sqrt(resid @ resid / (n - p)):.4f}"
          f"   mean h = {lev.mean():.4f}   max h = {lev.max():.4f}")
```

## Four ways a polynomial misbehaves

The pattern in @exm-ply-co2-degree is not particular to that series: it is the
visible part of four properties of the monomial basis and of the polynomial
class itself.

::: {#prp-ply-instability}
[Instability of a high-degree polynomial fit]

Let \( \X_d=[\bone,\x,\dots,\x^d] \) be the Vandermonde matrix of a design with
at least \( d+1 \) distinct points, and let \( \hat f_d \) be the least squares
polynomial of degree \( d \).

::: {.enumerate options="label=(\alph*)"}
1. *(Conditioning.)* If the design points are a sample from a distribution on
   \( [0,1] \) with density \( g \), then
   \( n^{-1}\X_d\T\X_d\to\mathbf{G} \) almost surely, where
   \( G_{jk}=\int_0^1x^{j+k}g(x)\,dx \). For the uniform density this limit is the
   \( (d+1)\times(d+1) \) Hilbert matrix \( G_{jk}=1/(j+k+1) \), whose condition number
   grows exponentially in \( d \).

2. *(Oscillation.)* Interpolation by a polynomial of degree \( n-1 \) at \( n \)
   equally spaced points need not converge as \( n \) grows, even for a function
   that is analytic on the interval; the error can diverge. The failure is a
   property of the node placement, not of polynomials: at Chebyshev nodes the
   same interpolants converge geometrically.

3. *(Extrapolation.)* Let \( \tilde{\x}_d=(\I-\M_{d-1})\x^d \) be the part of the
   top power orthogonal to the lower ones. Then \( t\mapsto\Var\{\hat f_d(t)\} \)
   is a polynomial of degree exactly \( 2d \) with leading coefficient
   \( \sigma^2/\norm{\tilde{\x}_d}^2 \), so
   \[
   \frac{\Var\{\hat f_d(t)\}}{\sigma^2t^{2d}}\longrightarrow
   \frac{1}{\norm{\tilde{\x}_d}^2}\qquad\text{as }\lvert t\rvert\to\infty .
   \]
   Each extra degree therefore multiplies the standard error of the
   extrapolated curve by a further factor \( \lvert t\rvert \), up to the
   constant \( \norm{\tilde{\x}_d}/\norm{\tilde{\x}_{d+1}} \).

4. *(Non-locality.)* Changing one response \( y_i \) by \( \Delta \) changes the
   fitted curve at every \( t \) by \( \Delta\,\ell_i(t) \), where
   \( \ell_i(t)=\mathbf{x}(t)\T(\X_d\T\X_d)^{-1}\mathbf{x}_i \) is itself a polynomial of
   degree at most \( d \) and never identically zero, so it vanishes at most
   \( d \) times: no observation can be given a purely local influence.
:::

:::

::: {.proof}
(a) The \( (j,k) \) entry of \( n^{-1}\X_d\T\X_d \) is \( n^{-1}\sum_ix_i^{j+k} \), an
average of bounded independent and identically distributed variables, which by
the strong law converges almost surely to \( \E(X^{j+k})=\int x^{j+k}g \); for
\( g\equiv1 \) that is \( 1/(j+k+1) \). The exponential growth of
\( \kappa(\mathbf{G}) \) is classical and not reproduced here: Todd (1954) shows it
is of order \( (1+\sqrt2)^{4d} \) up to algebraic factors, and Higham (2002)
treats the Hilbert matrix and Vandermonde systems. By @eq-cmp-rule the normal
equations (@thm-proj-normal-equations) then lose about
\( 2\log_{10}\kappa(\X_d) \) of the sixteen digits of double precision.

(b) This is Runge's phenomenon, stated here without proof: Runge (1901) gave
the example and Trefethen (2019) the modern account, Chebyshev nodes included.
Nothing later depends on the divergence itself, only on the weaker statement
that a polynomial forced to follow the data in one region must do something in
the others, and that something is an oscillation largest at the ends.

(c) With \( \mathbf{x}(t)=(1,t,\dots,t^d)\T \),
\[
\begin{aligned}
\Var\{\hat f_d(t)\}
&=\sigma^2\mathbf{x}(t)\T(\X_d\T\X_d)^{-1}\mathbf{x}(t)\\
&=\sigma^2\sum_{j=0}^{d}\sum_{k=0}^{d}\bigl[(\X_d\T\X_d)^{-1}\bigr]_{jk}t^{j+k},
\end{aligned}
\]
a polynomial in \( t \) whose coefficient of \( t^{2d} \) is
\( \sigma^2[(\X_d\T\X_d)^{-1}]_{dd} \). By @eq-proj-partitioned-inverse with
\( \X_1=\x^d \) and \( \X_2=\X_{d-1} \), that entry is
\( 1/\norm{(\I-\M_{d-1})\x^d}^2=1/\norm{\tilde{\x}_d}^2 \), which is positive
because \( \x^d\notin\C(\X_{d-1}) \) when the design has at least \( d+1 \) distinct
points (@exr-ply-vandermonde-rank). A polynomial with positive leading
coefficient \( c \) satisfies \( q(t)/t^{2d}\to c \), which is the claim.

(d) Least squares is linear in \( \y \), so the fitted curve changes by
\( \Delta\mathbf{x}(t)\T(\X_d\T\X_d)^{-1}\mathbf{x}_i \) when \( y_i \) changes by
\( \Delta \). In \( t \) this is a combination of \( 1,t,\dots,t^d \), hence a
polynomial of degree at most \( d \), and never the zero polynomial, since
\( (\X_d\T\X_d)^{-1}\mathbf{x}_i=\bzero \) would force \( \mathbf{x}_i=\bzero \)
while \( \mathbf{x}_i \) has leading entry one. Such a polynomial has at most
\( d \) real roots.
:::

Part (a) explains the numbers in @exm-proj-qr-polynomial, and part (b) is drawn
in panel (a) of [Figure 42.1.2](#fig-ply-runge). Part (d) is why a polynomial fit
is never local: moving the December 2001 carbon dioxide value by \( 10 \) ppm
moves the fitted trend by \( 0.170 \) ppm at the end of the
series and \( 0.059 \) ppm in March 1958 when \( d=2 \), and by
\( 2.749 \) and \( 0.197 \) ppm when
\( d=12 \), the degree-12 weight changing sign \( 12 \) times on the way.

::: {#exm-ply-runge}
[Runge's function]

Interpolate \( f(x)=1/(1+25x^2) \) at \( m \) equally spaced points of \( [-1,1] \) by
the polynomial of degree \( m-1 \). The largest error on the interval is
\( 1.92 \) for \( m=11 \), \( 59.82 \) for
\( m=21 \) and \( 2388.03 \) for \( m=31 \): the function is
smooth and bounded by one, and the interpolant does not converge to it. Move the
same nodes to the Chebyshev positions \( \cos\{(2k+1)\pi/2m\} \), which crowd
towards the ends, and the errors become \( 0.10915 \),
\( 0.01533 \), \( 0.00206 \). Only where
the information sits has changed.
:::

Trouble at the ends is a shortage of design points relative to what the fitted
function is asked to do there. Either add points, rarely an option, or ask the
function to do less — which is what knots and penalties do.

::: {#exm-ply-extrapolation}
[Eleven years out of sample]

Fit the same models to the \( 389 \) months up to the end of 1990 and predict the
\( 132 \) months of 1991–2001. Inside the fitting range the residual standard
deviation falls steadily with the degree — \( 0.439 \),
\( 0.428 \), \( 0.417 \),
\( 0.406 \) ppm at degrees 2, 4, 6, 8 — while outside it the
root mean squared prediction errors are \( 3.99 \),
\( 2.63 \), \( 13.92 \) and
\( 188.23 \) ppm. The degree-8 fit, the best of the four
in-sample, is out by \( -540.7 \) ppm in December 2001 against
an observed value near \( 371 \): it predicts a concentration
several hundred parts per million below zero.
:::

::: {when-format="html"}
![**Figure 42.1.2.** (a) Runge's phenomenon: equally spaced interpolation of \( 1/(1+25x^2) \),
11 and 21 nodes. (b) Polynomial trends fitted up to 1990, continued to 2002.](runge_extrapolation.svg){#fig-ply-runge width=100%}
:::

::: {when-format="pdf"}
![(a) Runge's phenomenon: equally spaced interpolation of \( 1/(1+25x^2) \),
11 and 21 nodes. (b) Polynomial trends fitted up to 1990, continued to 2002.](runge_extrapolation.pdf){width=100%}
:::

```{.python .run #cell-runge-extrapolation-runge}
import numpy as np

runge = lambda x: 1.0 / (1.0 + 25.0 * x ** 2)
grid = np.linspace(-1.0, 1.0, 2001)

for m in (5, 11, 21, 31):
    equal = np.linspace(-1.0, 1.0, m)                          # m equally spaced nodes
    cheb = np.cos((2 * np.arange(m) + 1) * np.pi / (2 * m))     # m Chebyshev nodes
    err = [np.max(np.abs(np.polyval(np.polyfit(z, runge(z), m - 1), grid) - runge(grid)))
           for z in (equal, cheb)]
    print(f"{m:2d} nodes, degree {m - 1:2d}:  equally spaced {err[0]:10.4f}   Chebyshev {err[1]:8.5f}")
```

## Centring and scaling are a partial repair

Part (a) is about the *basis*, and a basis can be changed without touching the
model. The cheapest change is affine: fit in
\( u=\{2x-\max x-\min x\}/(\max x-\min x) \), which runs over \( [-1,1] \), or at
least in \( x-\bar x \). By @thm-proj-reparam the fit is identical and only the
arithmetic differs. For the carbon dioxide design the condition number of the
trend columns is

| Degree | raw years | centred | rescaled to \( [-1,1] \) |
|---|---|---|---|
| 2 | \( 1.08\times 10^{11} \) | \( 3.19\times 10^{2} \) | \( 3.76\times 10^{0} \) |
| 4 | \( 1.36\times 10^{22} \) | \( 1.42\times 10^{5} \) | \( 1.90\times 10^{1} \) |
| 8 | \( 9.87\times 10^{40} \) | \( 3.11\times 10^{10} \) | \( 5.49\times 10^{2} \) |
| 12 | \( 1.48\times 10^{54} \) | \( 7.17\times 10^{15} \) | \( 1.69\times 10^{4} \) |
| 16 | \( 2.63\times 10^{67} \) | \( 3.25\times 10^{21} \) | \( 5.38\times 10^{5} \) |

Centring gains nine orders of magnitude at degree 2 and forty-six at degree 16,
because raising a number near 1980 to the sixteenth power is what made the
columns incommensurable; rescaling gains a further two to sixteen, close to the
best any diagonal rescaling can do (@prp-cmp-van-der-sluis). But the rescaled
design is still at \( 5.38\times 10^{5} \) by degree 16 and
grows geometrically after that, high powers of a number in \( [-1,1] \) all
looking like a spike at the ends. The next section removes that dependence
exactly.

```{.python .run #cell-polynomial-degree-conditioning}
for d in (2, 4, 8, 12, 16):
    raw = np.linalg.cond(np.vander(t, d + 1, increasing=True))
    centred = np.linalg.cond(np.vander(t - 1980.0, d + 1, increasing=True))
    scaled = np.linalg.cond(np.vander(u(t), d + 1, increasing=True))
    print(f"degree {d:2d}:  raw {raw:.2e}   centred {centred:.2e}   scaled {scaled:.2e}")
```

Two complaints have now been separated. *The basis is bad* is a numerical
statement, fixed exactly by a change of basis that no fitted value feels:
[Section 42.2](02-orthogonal-polynomials.html). *The family is bad* is a statement
about \( \mathcal P_d \) itself, which only a different space of functions can
fix: [Section 42.3](03-piecewise-polynomials.html).

## Exercises

### A. Check your understanding

::: {#exr-ply-vandermonde-rank}
[A1]

Show that \( \X_d=[\bone,\x,\dots,\x^d] \) has full column rank if and only if the
\( x_i \) take at least \( d+1 \) distinct values, and deduce that a polynomial of degree
\( n-1 \) through \( n \) points with distinct \( x_i \) exists, is unique and interpolates
exactly.
:::

::: {.solution}
If \( \X_d\mathbf{c}=\bzero \) then \( q(t)=\sum_jc_jt^j \), of degree at most \( d \),
vanishes at every distinct design value; with \( k\ge d+1 \) of them a nonzero \( q \)
would have more roots than its degree, so \( \mathbf{c}=\bzero \). If \( k\le d \), the
nonzero \( q(t)=\prod_{l\le k}(t-x_{(l)}) \) vanishes on the design. Taking \( d=n-1 \)
makes \( \X \) square of full rank, so the residuals vanish.
:::

### B. Practice

::: {#exr-ply-hilbert}
[B1]

For the design \( x_i=(i-1)/(n-1) \), show that the \( (j,k) \) entry of
\( n^{-1}\X_d\T\X_d \) converges to \( 1/(j+k+1) \) for each fixed \( j,k\le d \): part (a)
of @prp-ply-instability for a fixed rather than a random design.
:::

::: {#exr-ply-hat-at-endpoint}
[B2]

For the design \( x_i=-1+2(i-1)/(n-1) \) with \( n=50 \), compute the leverage of the
leftmost point for \( d=1,3,\dots,15 \) and plot it against \( d \), comparing with the
average \( (d+1)/n \) and the bound of @prp-proj-leverage.
:::

::: {#exr-ply-centring-invariance}
[B3]

Let \( \mathbf{u}=a\x+b\bone \) with \( a\neq0 \). Show that
\( \C([\bone,\mathbf{u},\dots,\mathbf{u}^d])=\C([\bone,\x,\dots,\x^d]) \) and that the
change of basis matrix is triangular with diagonal entries \( a^j \). Deduce that
centring and scaling leave \( \hY \), \( \he \) and every leverage unchanged, and say why
they nevertheless change \( \kappa(\X) \).
:::

::: {.solution}
Expanding \( (ax+b)^j \) writes \( \mathbf{u}^j \) in terms of \( \bone,\x,\dots,\x^j \)
with leading coefficient \( a^j \), so \( \mathbf{T} \) is upper triangular with
\( T_{jj}=a^j\neq0 \), hence nonsingular, and the column spaces agree. By
@thm-proj-reparam \( \M \) is unchanged, so fitted values, residuals and leverages are
too. The condition number is a property of \( \X \), not of \( \C(\X) \).
:::

::: {#exr-ply-influence-polynomial}
[B4]

Verify part (d) of @prp-ply-instability numerically: for the degree-12 trend on the
carbon dioxide design compute \( \ell_i(t) \) for the last observation and count its sign
changes over the data, then repeat for \( d=2 \). Why is the count for \( d=2 \) smaller
than \( 2 \) in some designs?
:::

::: {.solution}
The counts are \( 12 \) and \( 2 \). A polynomial of degree \( d \) has at most \( d \) real
roots but may have fewer, and roots outside the range of the data are not counted: if the
parabola's vertex lies outside \( [\min x,\max x] \), only one sign change is visible.
:::

### C. Going deeper

::: {#exr-ply-weierstrass}
[C1]

The Weierstrass approximation theorem says that every continuous \( f \) on \( [a,b] \) is a
uniform limit of polynomials. Using parts (a) and (c) of @prp-ply-instability, explain why
this is not an argument for fitting high-degree polynomials to data. What does the theorem
guarantee about the bias of the degree-\( d \) fit as \( d \) grows, and what about the
variance?
:::

::: {#exr-ply-lebesgue}
[C2]

For interpolation at nodes \( x_1<\dots<x_m \) the interpolant is
\( \sum_kf(x_k)L_k(t) \), with \( L_k \) the Lagrange basis, and the **Lebesgue constant**
is \( \Lambda=\max_t\sum_k\lvert L_k(t)\rvert \). Show that the interpolation error is at
most \( (1+\Lambda) \) times that of the best uniform approximation of degree \( m-1 \),
compute \( \Lambda \) for \( m=5,9,13,17 \) equally spaced and Chebyshev nodes, and relate
the two sequences to @exm-ply-runge.
:::
