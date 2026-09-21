# B-splines

[Chapter 42](../ch42-polynomials-piecewise/index.html) fixed the space we shall work
in. Given a degree \( d\ge0 \), an interval \( [a,b] \) and \( K \) interior knots
\( a<t_1<\dots<t_K<b \), the **spline space** is the set of functions
that are polynomials of degree at most \( d \) on each interval between
consecutive knots and have \( d-1 \) continuous derivatives at each
knot (@def-ply-piecewise). Its dimension is \( d+1+K \), and the truncated power
functions
\( 1,x,\dots,x^d,(x-t_1)_+^d,\dots,(x-t_K)_+^d \) are a
basis (@thm-ply-spline-space). (Chapter 42 writes \( \kappa_i \) for these knots;
here \( \kappa \) is needed for the longer sequence built next.)

That basis is a bad one to compute with, for the reason @prp-ply-knots gives:
its members are global, nearly collinear, and of wildly different sizes once the
knots are close together. This section builds the basis everything else in the
chapter uses. Its members are polynomials of the same space, so nothing about the
*fit* changes; only the numbers do, and they change dramatically.

::: {.warning}
From here to the end of the chapter, \( K \) is the number of interior knots, not
a kernel. The kernel \( K(\cdot) \) of [Section 43.1](01-kernels.html) does not
reappear.
:::

## The knot sequence and the recurrence

B-splines are defined by a recurrence on the degree, and the recurrence needs
knots outside \( [a,b] \) as well as inside it. Collect the interior knots and
\( d+1 \) copies of each endpoint into a nondecreasing **extended knot sequence**
\( \kappa_1\le\kappa_2\le\dots\le\kappa_{m+d+1} \) with \( m=K+d+1 \), namely
\[
\kappa_1=\dots=\kappa_{d+1}=a,\qquad \kappa_{d+1+i}=t_i\ \ (i=1,\dots,K),\qquad
\kappa_{d+K+2}=\dots=\kappa_{m+d+1}=b .
\]{#eq-smo-clamped}

This is the **clamped** sequence; an alternative is discussed below. In either
case there are \( m+d+1=K+2d+2 \) entries and, as we shall see, \( m=K+d+1 \)
B-splines, exactly the dimension of the spline space. From here on \( \kappa_j \)
always means an entry of the extended sequence.

::: {#def-smo-bspline}
[B-splines]

Given a nondecreasing sequence \( \kappa_1\le\dots\le\kappa_{m+d+1} \), define
\[
B_j^0(x)=\begin{cases}1,&\kappa_j\le x<\kappa_{j+1},\\0,&\text{otherwise},\end{cases}
\qquad j=1,\dots,m+d,
\]
and recursively, for \( r=1,\dots,d \) and \( j=1,\dots,m+d+1-r \),
\[
B_j^{r}(x)=\omega_j^{r}(x)\,B_j^{r-1}(x)
+\{1-\omega_{j+1}^{r}(x)\}\,B_{j+1}^{r-1}(x),
\]{#eq-smo-deboor}

where \( \omega_j^{r}(x)=(x-\kappa_j)/(\kappa_{j+r}-\kappa_j) \).

A term whose denominator vanishes is read as zero. The functions
\( B_1^d,\dots,B_m^d \) are the **B-splines of degree \( d \)** on the sequence,
and the matrix \( \B \) with \( \B_{ij}=B_j^d(x_i) \) is the **B-spline basis
matrix** of a design.
:::

The recurrence is due to de Boor (1972) and Cox (1972). It has an immediate
reading: \( B_j^r \) is a weighted average of two B-splines of one degree lower,
with weights that ramp linearly from \( 0 \) to \( 1 \) across the support. The
piece of the definition that does the work is the convention about zero
denominators, which is what lets the endpoints be repeated; a repeated knot makes
the corresponding lower-degree B-spline identically zero, and the offending term
disappears with it.

## What B-splines are like

::: {#prp-smo-bspline}
[Properties of the B-spline basis]

Let \( B_1^d,\dots,B_m^d \) be as in @def-smo-bspline. Then:

::: {.enumerate options="label=(\alph*)"}
1. *(Local support.)* \( B_j^d(x)=0 \) for \( x\notin[\kappa_j,\kappa_{j+d+1}) \),
   and, for the sequence @eq-smo-clamped with distinct interior knots,
   \( B_j^d(x)>0 \) for \( \kappa_j<x<\kappa_{j+d+1} \). In particular at most
   \( d+1 \) of the \( B_j^d \) are nonzero at any point.

2. *(Nonnegativity.)* \( B_j^d\ge0 \) everywhere.

3. *(Partition of unity.)* For the sequence @eq-smo-clamped,
   \( \sum_{j=1}^mB_j^d(x)=1 \) for every \( x\in[a,b) \).

4. *(Piecewise polynomial.)* On each interval between consecutive distinct knots,
   \( B_j^d \) is a polynomial of degree at most \( d \).

5. *(Derivative.)* For \( d\ge1 \),
   \[
   \frac{d}{dx}B_j^d(x)=d\Bigl\{\frac{B_j^{d-1}(x)}{\kappa_{j+d}-\kappa_j}
   -\frac{B_{j+1}^{d-1}(x)}{\kappa_{j+d+1}-\kappa_{j+1}}\Bigr\},
   \]{#eq-smo-bderiv}

   with the same convention on zero denominators. Consequently \( B_j^d \) has
   \( d-1 \) continuous derivatives at every simple knot.

6. *(Basis.)* For the sequence @eq-smo-clamped with distinct interior knots,
   \( B_1^d,\dots,B_m^d \) are linearly independent on \( [a,b] \) and span the
   spline space of @thm-ply-spline-space, of which they are therefore a basis.
   (With a repeated interior knot they span a larger space, one with a weaker
   continuity requirement there; see @exr-smo-multiplicity.)
:::

:::

::: {.proof}
(a) and (b) go together, by induction on \( r \). For \( r=0 \) both are clear.
Suppose they hold for \( r-1 \). By @eq-smo-deboor, \( B_j^r \) is a combination of
\( B_j^{r-1} \), supported in \( [\kappa_j,\kappa_{j+r}) \), and \( B_{j+1}^{r-1} \),
supported in \( [\kappa_{j+1},\kappa_{j+r+1}) \), so \( B_j^r \) vanishes outside
\( [\kappa_j,\kappa_{j+r+1}) \). On that interval
\( \omega_j^r(x)=(x-\kappa_j)/(\kappa_{j+r}-\kappa_j)\ge0 \) and
\( 1-\omega_{j+1}^r(x)=(\kappa_{j+r+1}-x)/(\kappa_{j+r+1}-\kappa_{j+1})\ge0 \),
so \( B_j^r\ge0 \). For strict positivity, let
\( \kappa_j<x<\kappa_{j+r+1} \). If \( x<\kappa_{j+r} \) then \( \kappa_j<\kappa_{j+r} \),
so \( \omega_j^r(x)>0 \) and \( B_j^{r-1}(x)>0 \) by the induction hypothesis, and
the second term is nonnegative. If \( x\ge\kappa_{j+r} \) then
\( \kappa_{j+1}<\kappa_{j+r+1} \) — otherwise \( \kappa_{j+r}=\kappa_{j+r+1} \), which
contradicts \( \kappa_{j+r}\le x<\kappa_{j+r+1} \) — so
\( 1-\omega_{j+1}^r(x)>0 \) and \( B_{j+1}^{r-1}(x)>0 \). Finally, a point lies in
at most \( d+1 \) of the intervals \( [\kappa_j,\kappa_{j+d+1}) \).

(c) Induct on \( r \), proving \( \sum_{j=1}^{m+d+1-r}B_j^{r}(x)=1 \) for
\( x\in[a,b) \). For \( r=0 \) the intervals \( [\kappa_j,\kappa_{j+1}) \) that are
nonempty partition \( [a,b) \). For the step, substitute @eq-smo-deboor and
collect terms with the same lower index:
\[
\begin{aligned}
\sum_{j}B_j^{r}
&=\omega_1^{r}B_1^{r-1}
+\sum_{j\ge2}\bigl\{\omega_j^{r}+1-\omega_j^{r}\bigr\}B_j^{r-1}\\
&\qquad+\{1-\omega_{L+1}^{r}\}B_{L+1}^{r-1},
\end{aligned}
\]
where \( L=m+d+1-r \) is the largest index of a degree-\( r \) B-spline. The first
and last terms vanish on \( [a,b) \): for @eq-smo-clamped,
\( B_1^{r-1} \) is supported in \( [\kappa_1,\kappa_{r})=[a,a) \), which is empty
when \( r\le d \), and \( B_{L+1}^{r-1} \) in
\( [\kappa_{L+1},\kappa_{L+r})=[b,b) \). What is left is
\( \sum_{j=2}^{L}B_j^{r-1} \), which is the partition of unity one degree lower.

(d) By induction: \( B_j^0 \) is constant on each such interval, and
@eq-smo-deboor multiplies polynomials of degree \( r-1 \) by linear functions.

(e) Induct on \( d \). Write \( D_j^{r}=(\kappa_{j+r}-\kappa_j)^{-1} \), so that
\( \omega_j^{r}(x)=(x-\kappa_j)D_j^{r} \) and
\( \{\omega_j^{r}\}'=D_j^{r} \). For \( d=1 \), \( B_j^1 \) is continuous and
piecewise linear, and differentiating @eq-smo-deboor on the interior of each knot
interval gives \( D_j^1B_j^0-D_{j+1}^1B_{j+1}^0 \), which is @eq-smo-bderiv.
For the step, differentiate @eq-smo-deboor:
\[
\begin{aligned}
(B_j^d)'&=D_j^{d}B_j^{d-1}-D_{j+1}^{d}B_{j+1}^{d-1}\\
&\qquad+\omega_j^{d}(B_j^{d-1})'+\{1-\omega_{j+1}^{d}\}(B_{j+1}^{d-1})' ,
\end{aligned}
\]
so it suffices to show that the last two terms equal
\( (d-1)\{D_j^{d}B_j^{d-1}-D_{j+1}^{d}B_{j+1}^{d-1}\} \). Substituting the
induction hypothesis on the left and @eq-smo-deboor on the right, both sides
become combinations of \( B_j^{d-2} \), \( B_{j+1}^{d-2} \) and \( B_{j+2}^{d-2} \),
and it is enough to match the three coefficients. For \( B_j^{d-2} \):
\( \omega_j^{d}D_j^{d-1}=(x-\kappa_j)D_j^{d}D_j^{d-1} \) on the left and
\( D_j^{d}\omega_j^{d-1}=D_j^{d}(x-\kappa_j)D_j^{d-1} \) on the right. For
\( B_{j+2}^{d-2} \): \( -\{1-\omega_{j+1}^{d}\}D_{j+2}^{d-1}
=-(\kappa_{j+d+1}-x)D_{j+1}^{d}D_{j+2}^{d-1} \) on the left, and
\( -D_{j+1}^{d}\{1-\omega_{j+2}^{d-1}\}=-D_{j+1}^{d}(\kappa_{j+d+1}-x)D_{j+2}^{d-1} \)
on the right. For \( B_{j+1}^{d-2} \) the two coefficients are
\( D_{j+1}^{d-1}\{1-\omega_j^{d}-\omega_{j+1}^{d}\} \) and
\( D_j^{d}-\omega_{j+1}^{d-1}(D_j^{d}+D_{j+1}^{d}) \), and their difference
reduces, after cancelling \( (x-\kappa_{j+1})D_{j+1}^{d-1}D_{j+1}^{d} \) and
collecting the terms in \( D_j^{d} \), to
\( D_j^{d}\{D_{j+1}^{d-1}(\kappa_{j+d}-\kappa_{j+1})-1\} \), which is zero because
\( D_{j+1}^{d-1}=(\kappa_{j+d}-\kappa_{j+1})^{-1} \). Smoothness now follows by a
second induction: \( B_j^1 \) is continuous, and if every B-spline of degree
\( d-1 \) is \( C^{d-2} \) at simple knots, then @eq-smo-deboor makes \( B_j^d \) at
least \( C^{d-2} \), while @eq-smo-bderiv makes its derivative a combination of
\( C^{d-2} \) functions, so \( B_j^d\in C^{d-1} \).

(f) By (d) and (e) each \( B_j^d \) belongs to the spline space, which by
@thm-ply-spline-space has dimension \( d+1+K=m \). Since there are \( m \) of them,
it is enough to prove independence, and we do so by induction on \( d \). For
\( d=0 \), \( m=K+1 \) and the \( B_j^0 \) are indicators of disjoint intervals.
Let \( d\ge1 \) and suppose \( s=\sum_{j=1}^mc_jB_j^d\equiv0 \) on \( [a,b] \).
By @eq-smo-bderiv, with \( c_0=c_{m+1}=0 \),
\[
s'=d\sum_{j=2}^{m}\frac{c_j-c_{j-1}}{\kappa_{j+d}-\kappa_j}\,B_j^{d-1},
\]{#eq-smo-sprime}

the terms \( j=1 \) and \( j=m+1 \) having vanished because
\( \kappa_{1+d}=\kappa_1=a \) and \( \kappa_{m+1+d}=\kappa_{m+1}=b \) for the
sequence @eq-smo-clamped. The knots \( \kappa_2,\dots,\kappa_{m+d} \) are the
clamped sequence of degree \( d-1 \) for the same \( K \) interior knots, and
\( B_2^{d-1},\dots,B_m^{d-1} \) are its \( K+d=m-1 \) B-splines. By the induction
hypothesis they are independent, so every \( c_j-c_{j-1}=0 \) for
\( j=2,\dots,m \) (the denominators are positive, since no \( d+1 \) consecutive
knots of @eq-smo-clamped coincide). Hence all \( c_j \) are equal to some \( c \),
and \( 0=s=c\sum_jB_j^d=c \) on \( [a,b) \) by (c). So \( c=0 \).
:::

::: {.idea}
Local support is the whole point. Since only \( d+1 \) B-splines are nonzero at
any \( x \), the matrix \( \B \) has \( d+1 \) nonzero entries per row and
\( \B\T\B \) is banded of bandwidth \( d \), and moving one coefficient
\( \gamma_j \) changes the curve only on \( [\kappa_j,\kappa_{j+d+1}) \). Neither
statement holds for the truncated power basis, where changing the coefficient of
\( (x-t_1)_+^d \) moves everything to the right of \( t_1 \).
:::

## The coefficients live at the Greville abscissae

A spline written in the B-spline basis, \( f=\sum_j\gamma_jB_j^d \), has
coefficients that are not values of \( f \), but they are close to being so, and
they sit at identifiable places.

::: {#prp-smo-greville}
[Greville abscissae]

Let \( \xi_j=(\kappa_{j+1}+\dots+\kappa_{j+d})/d \) for \( j=1,\dots,m \), the
**Greville abscissae** of the sequence. For @eq-smo-clamped,
\[
\sum_{j=1}^m\xi_j\,B_j^d(x)=x\qquad\text{for all }x\in[a,b] .
\]{#eq-smo-greville}

Consequently \( \sum_j(\alpha+\beta\xi_j)B_j^d(x)=\alpha+\beta x \): coefficients
that are an affine function of the Greville abscissae produce an affine function
of \( x \).
:::

::: {.proof}
Let \( s(x)=\sum_j\xi_jB_j^d(x) \). Since
\( \xi_j-\xi_{j-1}=(\kappa_{j+d}-\kappa_j)/d \), @eq-smo-sprime gives
\( s'(x)=\sum_{j=2}^{m}B_j^{d-1}(x)=1 \) on \( [a,b) \), by the partition of unity
of degree \( d-1 \) established in the proof of @prp-smo-bspline(c). So
\( s(x)=x+c \) for a constant \( c \), and \( s \) is continuous for \( d\ge2 \)
(for \( d=1 \) the identity is immediate). At \( x=a \) only \( B_1^d \) is nonzero
by @prp-smo-bspline(a), so \( B_1^d(a)=1 \) by the partition of unity, and
\( s(a)=\xi_1=(\kappa_2+\dots+\kappa_{d+1})/d=a \). Hence \( c=0 \). The second
statement follows by adding \( \alpha \) times the partition of unity.
:::

::: {.remark}
[A uniform knot sequence]

@eq-smo-clamped is not the only choice. Penalized splines
([Section 43.3](03-penalized-splines.html)) use instead a sequence that continues
the equal spacing \( \delta=(b-a)/(K+1) \) of the interior knots \( d \) steps
beyond each end:
\( a-d\delta,\dots,a-\delta,a,t_1,\dots,t_K,b,b+\delta,\dots,b+d\delta \).
There are again \( K+2d+2 \) knots and \( m=K+d+1 \) B-splines, and restricted to
\( [a,b] \) they span the same spline space, so every fitted curve is identical.
The advantage is that the Greville abscissae are then equally spaced,
\( \xi_j=a+\{j-(d+1)/2\}\delta \), consecutive ones \( \delta \) apart, which is
what makes the difference penalty of [Section 43.3](03-penalized-splines.html)
behave exactly as advertised. That this sequence too gives a basis of the spline
space on \( [a,b] \), and that @eq-smo-greville continues to hold for it, are
quoted from de Boor (2001, chs. IX–X) — the proofs of (c) and (f) above used the
repeated endpoints — and confirmed numerically in @exm-smo-conditioning.
:::

## Conditioning

::: {#exm-smo-conditioning}
[Two bases for the same space]

Take \( n=300 \) design points drawn uniformly on \( [0,1] \), cubic splines
(\( d=3 \)) and \( K=6 \) equally spaced interior knots, so the spline space has
dimension \( 10 \). Least squares in the truncated power basis
and least squares in the B-spline basis give fitted values that agree to
\( 1.5\times 10^{-13} \), and the
uniform sequence of the remark above agrees with the clamped one to
\( 1.8\times 10^{-15} \): three bases,
one fit. The condition numbers @def-mat-condition-number of the two design
matrices do not agree at all — here \( 5.3 \) and
\( 2.1\times 10^{4} \). The table shows how they move
with the number of knots.

| \( K \) | \( \kappa(\B) \), B-spline | \( \kappa \), truncated power |
|---|---|---|
| 5 | 5.3 | \( 1.25\times 10^{4} \) |
| 10 | 5.8 | \( 1.08\times 10^{5} \) |
| 20 | 5.5 | \( 1.01\times 10^{6} \) |
| 40 | 7.4 | \( 1.45\times 10^{7} \) |

[Figure 43.2.1](02-b-splines.html#fig-smo-bspline) draws the basis and a fit in
it. The B-spline condition number is essentially constant; the truncated power
one grows by a decade for every doubling of the knot count, and squaring it to
form a cross-product matrix costs twice as many digits again
([Section 6.10](../ch06-projections/10-computation.html)). At \( K=40 \) the
truncated power fit has lost most of double precision, the B-spline fit a
digit.
:::

::: {when-format="html"}
![**Figure 43.2.1.** (a) The ten cubic B-splines on \( [0,1] \) with six equally
spaced interior knots (vertical lines); the dashed line is their
sum (@prp-smo-bspline(c)). (b) A least squares fit in this basis: thin curves
\( \gamma_jB_j^3 \), thick curve their sum, dots the coefficients
\( \gamma_j \) at the Greville abscissae of @prp-smo-greville.](bspline_basis.svg){#fig-smo-bspline width=100%}
:::

::: {when-format="pdf"}
![(a) The ten cubic B-splines on \( [0,1] \) with six equally
spaced interior knots (vertical lines); the dashed line is their
sum (@prp-smo-bspline(c)). (b) A least squares fit in this basis: thin curves
\( \gamma_jB_j^3 \), thick curve their sum, dots the coefficients
\( \gamma_j \) at the Greville abscissae of @prp-smo-greville.](bspline_basis.pdf){width=100%}
:::

The whole recurrence is a dozen lines of numpy.

```{.python .run #cell-bsplines-basis}
import numpy as np


def knot_vector(a, b, K, d):
    """d+1 copies of a, K equally spaced interior knots, d+1 copies of b."""
    interior = a + (b - a) * np.arange(1, K + 1) / (K + 1)
    return np.concatenate([np.full(d + 1, a), interior, np.full(d + 1, b)])


def bspline_basis(x, kn, d):
    """The B-splines of degree d on the knot vector kn, by the de Boor recurrence."""
    x = np.atleast_1d(np.asarray(x, float))
    B = np.array([(x >= kn[j]) & (x < kn[j + 1]) for j in range(len(kn) - 1)], float).T
    last = np.max(np.nonzero(kn[:-1] < kn[1:])[0])       # close the rightmost interval
    B[x == kn[-1], last] = 1.0
    for deg in range(1, d + 1):
        C = np.zeros((len(x), B.shape[1] - 1))
        for j in range(C.shape[1]):
            if kn[j + deg] > kn[j]:                       # a term with a zero denominator
                C[:, j] += (x - kn[j]) / (kn[j + deg] - kn[j]) * B[:, j]
            if kn[j + deg + 1] > kn[j + 1]:               # contributes nothing
                C[:, j] += (kn[j + deg + 1] - x) / (kn[j + deg + 1] - kn[j + 1]) * B[:, j + 1]
        B = C
    return B


D, K, A, Bnd = 3, 6, 0.0, 1.0                            # cubic, six interior knots on [0, 1]
kn = knot_vector(A, Bnd, K, D)
```

```{.python .run #cell-bsplines-equivalence}
def truncated_power(x, interior, d):
    """The basis 1, x, ..., x^d, (x - kappa_1)_+^d, ..., (x - kappa_K)_+^d."""
    x = np.asarray(x, float)
    powers = [x**j for j in range(d + 1)]
    return np.column_stack(powers + [np.maximum(x - k, 0.0) ** d for k in interior])


rng = np.random.default_rng(43043)
n = 300
x = np.sort(rng.uniform(0, 1, n))
y = np.sin(7 * x) + 0.4 * x**2 + 0.3 * rng.normal(size=n)

interior = kn[D + 1:D + 1 + K]
Bx, Tx = bspline_basis(x, kn, D), truncated_power(x, interior, D)
fit_b = Bx @ np.linalg.lstsq(Bx, y, rcond=None)[0]
fit_t = Tx @ np.linalg.lstsq(Tx, y, rcond=None)[0]
print("largest difference between the two fits:", np.max(np.abs(fit_b - fit_t)))
print("condition numbers:", np.linalg.cond(Bx), np.linalg.cond(Tx))
```

## Choosing the knots, and why we shall stop choosing

With a B-spline basis, an unpenalized **regression spline** is ordinary least
squares: \( \hat{\bgamma}=(\B\T\B)^{-1}\B\T\y \), with \( d+1+K \) parameters and
the whole of Parts I–III available for its inference. Nothing is nonparametric
about it. The difficulty is the one @prp-ply-knots describes: the fit depends on
\( K \) and on where the knots are, the dependence is not smooth, and the
selection machinery of [Chapter 29](../ch29-model-selection/index.html) is then
applied to a discrete search whose answer is unstable. Quantile-spaced knots are
better than equally spaced ones when the design is uneven, but no placement rule
removes the problem: the fit is a step function of \( K \).
[Section 43.3](03-penalized-splines.html) takes the opposite route — fix \( K \)
at more knots than the data can support and control the fit with a penalty — so
that a continuous dial takes the place of a discrete one.

## Exercises

### A. Check your understanding

::: {#exr-smo-count}
[A1]

For \( d=3 \) and \( K=6 \) interior knots, write out the clamped knot sequence
@eq-smo-clamped and check that it has \( K+2d+2 \) entries and yields
\( m=K+d+1 \) B-splines. Which B-splines are nonzero at \( x=a \)? At the
midpoint of \( [\kappa_3,\kappa_4] \)?
:::

### B. Practice

::: {#exr-smo-bspline-numeric}
[B1]

Using the cell above, verify numerically for \( d=3 \), \( K=6 \) that
\( \sum_jB_j^3(x)=1 \), that \( B_j^3 \) vanishes outside
\( [\kappa_j,\kappa_{j+4}] \), that at most four B-splines are nonzero at any
point, and that @eq-smo-greville holds. Then check @eq-smo-bderiv against a
central difference.
:::

::: {#exr-smo-bandwidth-matrix}
[B2]

Show that \( \B\T\B \) has \( (\B\T\B)_{jk}=0 \) whenever
\( \lvert j-k\rvert>d \), and deduce that a Cholesky
factorization (@thm-cmp-cholesky) of \( \B\T\B \) costs \( O(nd^2) \) rather than
\( O(nm^2) \) operations. Why does the same statement fail for the truncated power
basis?
:::

::: {.solution}
\( (\B\T\B)_{jk}=\sum_iB_j^d(x_i)B_k^d(x_i) \), and the supports
\( [\kappa_j,\kappa_{j+d+1}) \) and \( [\kappa_k,\kappa_{k+d+1}) \) are disjoint
when \( \lvert j-k\rvert>d \), so every term is zero. The Cholesky factor
inherits the bandwidth: in @eq-cmp-cholesky, if \( r_{ki}=0 \) whenever
\( i-k>d \) for all \( k<j \), then \( r_{ji} \) is a multiple of
\( a_{ji}=0 \) for \( i-j>d \), so by induction the factorization of
@thm-cmp-cholesky visits only \( O(md^2) \) entries rather than \( O(m^3) \), and
forming the band of \( \B\T\B \) costs \( O(nd^2) \). The
truncated power functions have nested, not local, supports: the last one is
nonzero on \( (\kappa_K,b] \) and the first everywhere, so no entry of the
cross-product matrix is structurally zero.
:::

### C. Going deeper

::: {#exr-smo-multiplicity}
[C1]

Suppose an interior knot is repeated \( \nu \) times, \( 1\le\nu\le d+1 \). Show
from @eq-smo-bderiv, by induction on \( \nu \), that the B-splines are
\( C^{d-\nu} \) there, and that \( \nu=d+1 \) allows a jump. How would you use this
to fit a curve that is smooth everywhere except at one known breakpoint, where it
is only continuous?
:::

::: {#exr-smo-marsden}
[C2]

Prove that \( \sum_j\psi_j(y)B_j^d(x)=(y-x)^d \) for all \( x\in[a,b) \) and all
\( y \), where \( \psi_j(y)=(y-\kappa_{j+1})\cdots(y-\kappa_{j+d}) \), for
\( d=0 \) and \( d=1 \). (This is Marsden's identity; the general case is de Boor,
2001, ch. IX. Differentiating it \( d-1 \) times in \( y \) and setting
\( y=x \) recovers @eq-smo-greville.)
:::

::: {.solution}
For \( d=0 \) the claim is \( \sum_jB_j^0(x)=1 \), the empty product being one.
For \( d=1 \), \( \psi_j(y)=y-\kappa_{j+1} \). On \( [\kappa_i,\kappa_{i+1}) \) only
\( B_{i-1}^1 \) and \( B_i^1 \) are nonzero, with values
\( (\kappa_{i+1}-x)/(\kappa_{i+1}-\kappa_i) \) and
\( (x-\kappa_i)/(\kappa_{i+1}-\kappa_i) \). Their combination is
\[
\frac{(y-\kappa_i)(\kappa_{i+1}-x)+(y-\kappa_{i+1})(x-\kappa_i)}{\kappa_{i+1}-\kappa_i}
=y-x ,
\]
after expanding, since the terms in \( \kappa_i\kappa_{i+1} \) cancel.
:::

