# Piecewise polynomials

Changing the basis left the second complaint untouched: polynomials are global,
and a degree-16 fit spends the same sixteen degrees of freedom everywhere, so a
feature at one end is paid for by an oscillation at the other. The remedy is a
different space of functions: cut the range at a few points, fit a low-degree
polynomial on each piece, and join the pieces smoothly. Nothing oscillates, and
the result is still linear in its coefficients.

## Knots and pieces

::: {#def-ply-piecewise}
[Piecewise polynomials and continuity constraints]

Fix an interval \( [a,b] \) and **knots**
\( a=\kappa_0<\kappa_1<\dots<\kappa_K<\kappa_{K+1}=b \), the \( K \) inner ones being
the **interior knots**. A function \( f:[a,b]\to\Real \) is a **piecewise
polynomial of degree \( d \)** if its restriction to each half-open piece
\( [\kappa_k,\kappa_{k+1}) \), \( k=0,\dots,K-1 \), and to \( [\kappa_K,b] \) agrees
with a polynomial of degree at most \( d \); the value at a knot is the limit from
the right, and nothing ties the pieces together. These form a vector space
\( \mathcal{PP}_d(\kappa) \) of dimension \( (K+1)(d+1) \).

For \( -1\le m\le d-1 \), \( f \) satisfies the **continuity constraints of order
\( m \)** at \( \kappa_k \) if
\[
f^{(r)}(\kappa_k^{-})=f^{(r)}(\kappa_k^{+}),\qquad r=0,1,\dots,m,
\]
the case \( m=-1 \) imposing nothing. The **spline space**
\( \mathcal S_d(\kappa) \) is the set of piecewise polynomials of degree \( d \)
satisfying the constraints of order \( d-1 \) at every interior knot; its members
are the **splines of degree \( d \) with knots \( \kappa \)**.
:::

The constraints are linear in the coefficients: writing the unrestricted
piecewise polynomial as \( \sum_{k}\sum_{j\le d}c_{kj}(t-\kappa_k)^j
\mathbf{1}\{t\in[\kappa_k,\kappa_{k+1})\} \), the requirement
\( f^{(r)}(\kappa_{k}^{-})=f^{(r)}(\kappa_{k}^{+}) \) reads
\[
\sum_{j\ge r}\frac{j!}{(j-r)!}c_{k-1,j}(\kappa_k-\kappa_{k-1})^{j-r}=r!\,c_{kr},
\]
one equation for each \( r \) and each interior knot. Fitting under them is
restricted least squares (@prp-glh-restricted-ls), and comparing two levels of
smoothness is a general linear hypothesis (@thm-glh-general-f). Counting
equations gives a candidate dimension, \( (K+1)(d+1) \) coefficients minus
\( (m+1)K \) constraints, which for \( m=d-1 \) is
\[
(K+1)(d+1)-dK=d+1+K .
\]
Only a candidate, since it assumes the constraints independent: the next theorem
proves they are, by exhibiting a basis of exactly that size.

## The spline space and the truncated power basis

Write \( (v)_+=\max(v,0) \), and for \( d\ge1 \) let
\[
(t-\kappa)_+^{d}=\begin{cases}(t-\kappa)^{d}, & t\ge\kappa,\\ 0,& t<\kappa,\end{cases}
\]
the **truncated power function** of degree \( d \) with knot \( \kappa \). For \( d=0 \)
the same formula with the convention \( 0^0=1 \) gives the indicator of
\( [\kappa,\infty) \).

::: {#thm-ply-spline-space}
[Dimension of the spline space and the truncated power basis]

Let \( d\ge0 \) and let \( \kappa_1<\dots<\kappa_K \) be interior knots of \( [a,b] \).

::: {.enumerate options="label=(\alph*)"}
1. \( \mathcal S_d(\kappa) \) is a vector space of dimension \( d+1+K \), and
   \[
   1,\ t,\ \dots,\ t^{d},\ (t-\kappa_1)_+^{d},\ \dots,\ (t-\kappa_K)_+^{d}
   \]{#eq-ply-truncated}

   is a basis of it.

2. More generally, if the constraints imposed at \( \kappa_k \) are of order
   \( m_k \) with \( -1\le m_k\le d-1 \), the resulting space has dimension
   \( d+1+\sum_{k=1}^{K}(d-m_k) \), with basis the powers \( 1,\dots,t^d \) together
   with \( (t-\kappa_k)_+^{r} \) for \( r=m_k+1,\dots,d \) and \( k=1,\dots,K \).
:::

:::

::: {.proof}
It is enough to prove (b); part (a) is the case \( m_k=d-1 \) for every \( k \), for
which the count is \( d+1+K \).

*The space.* Sums and scalar multiples of piecewise polynomials of degree \( d \)
are piecewise polynomials of degree \( d \), and the constraints are linear
equations, so the set is a subspace of the space of functions on \( [a,b] \).

*The listed functions belong to it.* Each power \( t^j \) is a single polynomial,
hence in every such space. Fix \( k \) and \( r>m_k \). The function
\( g(t)=(t-\kappa_k)_+^{r} \) is zero on \( [a,\kappa_k) \) and
\( (t-\kappa_k)^{r} \) on \( [\kappa_k,b] \), a piecewise polynomial of degree
\( r\le d \) whose only possible breakpoint is \( \kappa_k \). For \( s\le r-1 \) its
\( s \)th derivative at \( \kappa_k \) is \( \{r!/(r-s)!\}(t-\kappa_k)^{r-s} \) there,
which is zero, and zero from the left; so \( g \) satisfies the constraints of
order \( r-1 \), in particular those of order \( m_k\le r-1 \). At the other knots
\( g \) is a single polynomial and satisfies any constraint.

*Linear independence.* Suppose
\[
q(t)+\sum_{k=1}^{K}\sum_{r>m_k}c_{kr}(t-\kappa_k)_+^{r}=0
\qquad\text{for all }t\in[a,b],
\]
with \( q\in\mathcal P_d \). On \( [a,\kappa_1) \) every truncated power vanishes, so
\( q \) vanishes on an interval and hence \( q\equiv0 \). Suppose inductively that
\( c_{lr}=0 \) for every \( l<k \). On \( [\kappa_k,\kappa_{k+1}) \) the displayed sum
reduces to \( \sum_{r>m_k}c_{kr}(t-\kappa_k)^{r} \), a polynomial vanishing on an
interval, so all its coefficients are zero. The induction gives \( c_{kr}=0 \)
throughout.

*Spanning.* Let \( f \) belong to the space. On \( [\kappa_0,\kappa_1) \) it agrees
with a polynomial \( q\in\mathcal P_d \); set \( f_1=f-q \), again in the space and
vanishing there. Suppose inductively that \( f_k \) is in the space and vanishes
on \( [a,\kappa_k) \), and on the next piece agrees with
\( r_k\in\mathcal P_d \).
Expanding \( r_k \) in powers of \( t-\kappa_k \),
\[
r_k(t)=\sum_{s=0}^{d}\frac{r_k^{(s)}(\kappa_k)}{s!}(t-\kappa_k)^{s},
\]
and for \( s\le m_k \) the constraint at \( \kappa_k \) forces
\( r_k^{(s)}(\kappa_k)=f_k^{(s)}(\kappa_k^{-})=0 \), because \( f_k \) is identically
zero to the left. Hence \( r_k(t)=\sum_{s>m_k}c_{ks}(t-\kappa_k)^{s} \) for some
coefficients \( c_{ks} \), and
\( f_{k+1}=f_k-\sum_{s>m_k}c_{ks}(t-\kappa_k)_+^{s} \) is in the space and
vanishes on \( [a,\kappa_{k+1}) \): to the left of \( \kappa_k \) both terms vanish,
and between \( \kappa_k \) and \( \kappa_{k+1} \) they cancel. After \( K \) steps
\( f_{K+1} \) vanishes on all of \( [a,b] \), which exhibits \( f \) as a combination of
the listed functions.

*Dimension.* The listed functions are \( d+1+\sum_k(d-m_k) \) in number, and they
are a basis.
:::

::: {.remark}
This settles the counting question: \( (K+1)(d+1) \) minus
\( d+1+\sum_k(d-m_k) \) is \( \sum_k(m_k+1) \), exactly the number of constraint
equations. The constraints are therefore independent, and that a free piecewise
polynomial is in fact a spline is testable with \( \sum_k(m_k+1) \) degrees of
freedom (@def-glh-testable).
:::

## Three familiar special cases

**Degree 0.** \( \mathcal S_0(\kappa) \) is the space of step functions, of
dimension \( K+1 \), and fitting it is the one-way layout with the intervals as
groups (@exm-proj-oneway-M): the fit is *exactly local*, with variance
\( \sigma^2/n_k \) on interval \( k \), at the price of the discontinuity and a
bias of order the interval width.

**Degree 1.** \( \mathcal S_1(\kappa) \), of dimension \( K+2 \), is the space of
continuous broken lines, with basis
\( 1,t,(t-\kappa_1)_+,\dots,(t-\kappa_K)_+ \); the coefficient of
\( (t-\kappa_k)_+ \) is the *change* in slope there, and the two-piece case is the
classical two-phase regression.

**Degree 3.** \( \mathcal S_3(\kappa) \) is the space of **cubic splines**, of
dimension \( K+4 \): value, slope and curvature continuous, only the third
derivative jumping. Cubic is the default, the lowest degree whose joins are
invisible and whose second derivative a smoothness penalty can control.

::: {.warning}
A knot estimated from the data takes the model out of this chapter:
\( (t-\kappa)_+^d \) is not linear in \( \kappa \), the criterion is not even
differentiable in it for \( d=0 \) or \( 1 \), and the usual asymptotics can fail.
This **segmented** or **changepoint** regression has a literature of its own —
Hudson (1966), Feder (1975), Seber and Wild (1989, ch. 9) — which this book does
not develop. The honest options are to fix the knot from outside the data, to
fit many fixed knots and let a penalty
decide ([Chapter 43](../ch43-smoothing/index.html)), or to pay for the
selection (@prp-sel-selection-bias).
:::

## Climbing the smoothness ladder

::: {#exm-ply-nile-ladder}
[Four fits with the same knots]

Take the Nile series again, \( n=100 \) annual flows, with
cubic pieces and interior knots at 1899 — where Cobb (1978) and the analysis of
@exm-ss-nile-polynomials both locate a change — and 1935. Fitting the three
cubics separately uses \( 12 \) parameters and leaves a residual sum of squares
of \( 1394072 \). Continuity at both knots removes two
parameters and raises it to \( 1561452 \); continuity of the
first derivative removes two more (\( 1613581 \)); and the full
cubic spline, with \( 6 \) parameters as @thm-ply-spline-space
predicts, has \( 1781797 \).
[Figure 42.3.1](#fig-ply-ladder) shows the four fits.

Each step down the ladder is a testable linear hypothesis, and against the
unrestricted fit all three are rejected: \( F=5.28 \) on
\( (2,88) \) degrees of freedom for continuity alone
(\( p=0.0068 \)), \( F=3.46 \) on \( (4,88) \)
(\( p=0.0112 \)) for a continuous first derivative, and
\( F=4.08 \) on \( (6,88) \)
(\( p=0.0012 \)) for the cubic spline. The Nile really does
have a discontinuity, and a spline is the wrong space for it.

The right model is the simplest on the ladder: a step function with a single
knot at 1899, whose group means are \( 1097.7 \) and
\( 850.0 \) in units of \( 10^8\,\text{m}^3 \), differing
by \( -247.8 \) with standard error
\( 28.4 \) and \( t=-8.71 \). Its
residual sum of squares, \( 1597457 \), beats the
six-parameter spline and the cubic polynomial's
\( 1909955 \). Two parameters in the right place beat six in
the wrong space.
:::

::: {when-format="html"}
![**Figure 42.3.1.** The Nile flows with cubic pieces cut at 1899 and 1935, under four levels of
continuity; each constraint removes two parameters.](piecewise_nile.svg){#fig-ply-ladder width=100%}
:::

::: {when-format="pdf"}
![The Nile flows with cubic pieces cut at 1899 and 1935, under four levels of
continuity; each constraint removes two parameters.](piecewise_nile.pdf){width=100%}
:::

```{.python .run #cell-piecewise-nile-bases}
import numpy as np
import statsmodels.api as sm

nile = sm.datasets.nile.load_pandas().data
year = nile["year"].to_numpy().astype(float)
y = nile["volume"].to_numpy()
n = len(y)
knots = np.array([1899.0, 1935.0])
d = 3
u = lambda s: 2 * (s - year.min()) / (year.max() - year.min()) - 1


def separate(s):
    """One free cubic per interval: (K+1)(d+1) columns, no constraints at the knots."""
    piece = np.digitize(s, knots)
    return np.column_stack([np.where(piece == k, u(s) ** j, 0.0)
                            for k in range(len(knots) + 1) for j in range(d + 1)])


def truncated(s, m):
    """Continuity of the function and of its first m derivatives at every knot."""
    cols = [u(s) ** j for j in range(d + 1)]
    for k in knots:
        cols += [np.clip(u(s) - u(k), 0.0, None) ** r for r in range(m + 1, d + 1)]
    return np.column_stack(cols)


def fit(X):
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    resid = y - X @ beta
    return np.linalg.matrix_rank(X), resid @ resid, X @ beta


for name, X in [("separate cubics", separate(year)), ("continuous", truncated(year, 0)),
                ("continuous derivative", truncated(year, 1)), ("cubic spline", truncated(year, 2))]:
    p, rss, _ = fit(X)
    print(f"{name:24s} p = {p:2d}   residual sum of squares {rss:10.0f}")
```

```{.python .run #cell-piecewise-nile-ftests}
from scipy import stats

p_full, rss_full, _ = fit(separate(year))
for m, label in ((0, "continuity"), (1, "first derivative"), (2, "second derivative")):
    p_m, rss_m, _ = fit(truncated(year, m))
    q = p_full - p_m
    F = ((rss_m - rss_full) / q) / (rss_full / (n - p_full))
    print(f"up to the {label:17s}: q = {q}, F = {F:5.2f}, p-value {1 - stats.f.cdf(F, q, n - p_full):.4f}")
```

::: {.warning}
The knot at 1899 was not chosen blind: Cobb's analysis and the large quartic
component in @exm-ss-nile-polynomials both point to the end of the nineteenth
century, and both used this series. The \( t \) statistic
\( -8.71 \) therefore overstates the evidence for a break of
unknown date, whose null distribution is that of the largest of ninety-nine
statistics, one for each year one could split at.
:::

::: {.idea}
Degrees of freedom are the currency, and a spline spends them differently. Each
of a polynomial's \( d+1 \) parameters affects the curve everywhere; of a
spline's \( d+1+K \), the \( K \) extra act only to the right of their knots,
adjusting the \( d \)th derivative. Raising \( K \) at fixed low \( d \) is a local
purchase; raising \( d \) is not.
:::

::: {.remark}
[Natural boundary conditions]

Outside \( [\kappa_1,\kappa_K] \) a cubic spline is still a cubic, and cubics
diverge. **Natural** cubic splines add \( f''(a)=f'''(a)=f''(b)=f'''(b)=0 \),
forcing \( f \) to be linear beyond the outer knots and reducing the dimension to
\( K \) (@exr-ply-natural-dimension) — a decision that the curve should flatten
where the data run out, not a discovery about them. It becomes more than a
convention in [Chapter 43](../ch43-smoothing/index.html), where the natural cubic
spline exactly minimizes a penalized criterion (@thm-smo-smoothing-spline).
:::

## Exercises

### A. Check your understanding

::: {#exr-ply-dimension-count}
[A1]

How many parameters has a quadratic spline with 7 interior knots? A cubic spline with
7, one of them a double knot (only \( C^1 \) there)? A piecewise constant function with
7?
:::

::: {.solution}
By @thm-ply-spline-space, \( 2+1+7=10 \), then \( 3+1+7+1=12 \) (the double knot
contributes \( d-m=3-1=2 \) instead of \( 1 \)), and \( 0+1+7=8 \).
:::

::: {#exr-ply-slope-change}
[A2]

In the broken-line model \( f(t)=\beta_0+\beta_1t+\sum_k\theta_k(t-\kappa_k)_+ \), give
the slope of \( f \) on \( (\kappa_j,\kappa_{j+1}) \), the hypothesis that the line is
unbent at \( \kappa_j \), and the degrees of freedom of the test that it is a single
straight line.
:::

::: {#exr-ply-step-locality}
[A3]

Verify that for the piecewise constant fit, changing \( y_i \) by \( \Delta \) moves the
fitted curve by \( \Delta/n_k \) on \( i \)'s own interval and not at all elsewhere.
Contrast with @prp-ply-instability(d).
:::

### B. Practice

::: {#exr-ply-natural-dimension}
[B1]

Show that the natural cubic splines with interior knots \( \kappa_1<\dots<\kappa_K \)
(\( K\ge2 \)) form a space of dimension \( K \), by writing the four boundary conditions
in the truncated power basis and checking that they are independent.
:::

::: {.solution}
Write \( f(t)=\sum_{j\le3}\beta_jt^j+\sum_k\theta_k(t-\kappa_k)_+^3 \), with \( K+4 \)
coefficients. Linearity on \( [a,\kappa_1] \) means \( \beta_2=\beta_3=0 \); linearity on
\( [\kappa_K,b] \) makes the cubic and quadratic coefficients of the last piece vanish,
which after expanding \( (t-\kappa_k)^3 \) reads \( \sum_k\theta_k=0 \) and
\( \sum_k\theta_k\kappa_k=0 \). The four equations involve different coefficients and
are independent, so the dimension is \( K \).
:::

::: {#exr-ply-two-phase}
[B2]

Two straight lines meet at a known \( \gamma \). Write the model in the basis
\( 1,t,(t-\gamma)_+ \), give the least squares estimate of the common value at \( \gamma \)
and its variance, and the \( F \) statistic for equal slopes.
:::

::: {.solution}
With \( \E(Y)=\beta_0+\beta_1t+\theta(t-\gamma)_+ \), the value at \( \gamma \) is the
estimable function \( \beta_0+\beta_1\gamma \), estimated by
\( \hat\beta_0+\gamma\hat\beta_1 \) with variance
\( \sigma^2\mathbf{a}\T(\X\T\X)^{-1}\mathbf{a} \), \( \mathbf{a}=(1,\gamma,0)\T \). Equal slopes
is \( \theta=0 \), tested by the square of the \( t \) statistic for \( \hat\theta \) against
\( F(1,n-3) \) (@thm-glh-t-test).
:::

::: {#exr-ply-derivative-jump}
[B3]

Let \( f\in\mathcal S_d(\kappa) \) have truncated power coefficients \( \theta_k \). Show
that \( f^{(d)} \) is a step function jumping by \( d!\,\theta_k \) at \( \kappa_k \), and
that \( \theta_k=0 \) exactly when \( f \) is a single polynomial across \( \kappa_k \).
:::

::: {.solution}
Differentiating @eq-ply-truncated \( d \) times kills every power below \( t^d \) and
turns \( (t-\kappa_k)_+^d \) into \( d! \) times the indicator of \( [\kappa_k,\infty) \), so
\( f^{(d)}=d!\,\beta_d+d!\sum_k\theta_k\mathbf{1}\{t\ge\kappa_k\} \). If \( \theta_k=0 \)
the pieces on either side of \( \kappa_k \) have the same \( d \)th derivative and agree
to order \( d-1 \) there, hence are the same polynomial; the converse is immediate.
:::

### C. Going deeper

::: {#exr-ply-peano}
[C1]

Show that for \( f\in C^{d+1}[a,b] \),
\[
\begin{aligned}
f(t)&=\sum_{j=0}^{d}\frac{f^{(j)}(a)}{j!}(t-a)^{j}
+\frac{1}{d!}\int_a^{b}f^{(d+1)}(s)\,(t-s)_+^{d}\,ds .
\end{aligned}
\]
Interpret it as saying that the truncated power functions are a continuous basis for
smooth functions, of which @eq-ply-truncated keeps \( K \) members.
:::

::: {#exr-ply-knot-coalescence}
[C2]

Let \( \kappa \) and \( \kappa+\epsilon \) be two knots of a cubic spline. Show that as
\( \epsilon\downarrow0 \) the angle between \( (t-\kappa)_+^3 \) and
\( (t-\kappa-\epsilon)_+^3 \) tends to zero, while the normalized difference
\( \{(t-\kappa)_+^3-(t-\kappa-\epsilon)_+^3\}/(3\epsilon) \) tends to \( (t-\kappa)_+^2 \).
Explain how this both shows that a double knot lowers the smoothness by one order and
warns against placing knots close together in the truncated power basis.
:::

::: {.solution}
For \( t\ge\kappa+\epsilon \),
\( \{(t-\kappa)^3-(t-\kappa-\epsilon)^3\}/(3\epsilon)\to(t-\kappa)^2 \) by the definition
of the derivative; for \( t\le\kappa \) both terms vanish, and the middle interval, of
width \( \epsilon \), contributes nothing. The limit is \( (t-\kappa)_+^2 \), the basis
function @thm-ply-spline-space(b) adds at a knot of multiplicity two, where only
\( C^1 \) continuity is required. Two nearly coincident columns whose difference carries
the information is what makes a design ill-conditioned (@def-mat-condition-number), and
it inflates the coefficients' variances in the sense of @def-col-vif: the fit is fine,
the coefficients are not.
:::
