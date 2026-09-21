# Regression splines and truncated power bases

A spline space is finite-dimensional, so fitting a spline is least squares with
the columns of @eq-ply-truncated. Nothing new is needed to *do* it; what is
needed is judgement about the two things the analyst now chooses, how many knots
and where.

::: {#def-ply-regression-spline}
[Regression spline]

Given a degree \( d \) and fixed interior knots \( \kappa_1<\dots<\kappa_K \) inside
the range of the design, the **regression spline** model is the linear model
with model matrix
\[
\X=\bigl[\bone,\ \x,\ \dots,\ \x^{d},\ (\x-\kappa_1)_+^{d},\ \dots,\
(\x-\kappa_K)_+^{d}\bigr],
\]
of \( d+1+K \) columns, and the **regression spline** is the least squares fit
\( \hat f(t)=\mathbf{x}(t)\T\hbeta \). Other regressors, smooth terms and factors
may be added as further columns.
:::

Being an ordinary linear model, it has \( n-(d+1+K) \) residual degrees of
freedom at full rank; @thm-ci-estimable-interval gives the fitted curve an
interval and @thm-ci-working-hotelling a band, both conditional on the knots;
the leverages average \( (d+1+K)/n \), which is what "\( K \) knots cost \( K \)
degrees of freedom" means; adding a knot is a one-degree-of-freedom nested
test (@thm-glh-f-test) asking, by @exr-ply-derivative-jump, whether the \( d \)th
derivative really jumps there; and choosing \( K \) by GCV (@def-sel-gcv) or
AIC (@prp-sel-aic-bic) belongs to
[Chapter 29](../ch29-model-selection/index.html).

## What the knots control

::: {#prp-ply-knots}
[Number, placement and basis]

Fit a regression spline of degree \( d \) with interior knots
\( \kappa_1<\dots<\kappa_K \) inside the range of the design, and write
\( \M_K \) for the projection onto the column space.

::: {.enumerate options="label=(\alph*)"}
1. *(Number.)* If \( \X \) has full rank, \( \tr\M_K=d+1+K \), so the leverages
   average \( (d+1+K)/n \) and \( \E\norm{(\I-\M_K)\Y}^2=\sigma^2\{n-(d+1+K)\}
   +\norm{(\I-\M_K)\bmu}^2 \). Increasing \( K \) decreases both the squared bias
   and the residual degrees of freedom; the price is the total variance of the
   fitted curve, \( \sigma^2\tr\M_K=\sigma^2(d+1+K) \), which rises by
   \( \sigma^2 \) per knot. Knots buy bias reduction at one degree of freedom
   each.

2. *(Rank.)* If each piece \( [\kappa_k,\kappa_{k+1}) \), \( k=0,\dots,K-1 \), and
   the last piece \( [\kappa_K,b] \) contain at least \( d+1 \) distinct design
   points, then \( \X \) has full column rank. If some knot exceeds
   \( \max_ix_i \), its column is identically zero and \( \X \) does not.

3. *(Placement.)* For \( d=0 \) the variance of the fit on interval \( k \) is
   \( \sigma^2/n_k \), with \( n_k \) the number of design points there, so
   equalizing the counts minimizes the largest pointwise variance; quantile
   knots do that and equally spaced knots need not. For \( d\ge1 \) the effect
   persists in the exact variance
   \( \sigma^2\mathbf{x}(t)\T(\X\T\X)^{-1}\mathbf{x}(t) \), which is large wherever
   an interval is sparsely occupied.

4. *(Basis.)* The columns of @eq-ply-truncated are globally supported and
   increasingly parallel: beyond the last knot every one behaves like \( t^d \),
   the angle between \( (\x-\kappa_j)_+^d \) and \( (\x-\kappa_k)_+^d \) tends to
   zero as \( \kappa_j-\kappa_k\to0 \) (@exr-ply-knot-coalescence), and
   \( \kappa(\X) \) grows with \( K \). A basis of the *same space* whose members
   are each supported on \( d+1 \) consecutive intervals has a **basis condition
   number** — the worst ratio between the size of a coefficient vector and the
   size of the function it builds — bounded by a constant depending only on
   \( d \), whatever the knots (de Boor, 2001). For a design with points in every
   basis function's support, \( \kappa(\X) \) then does not grow with \( K \);
   without one nothing is claimed, and \( \X \) can be singular.
:::

:::

::: {.proof}
(a) \( \M_K \) is an orthogonal projection, so its trace is its
rank (@prp-proj-trace-rank), namely \( d+1+K \) by @thm-ply-spline-space at full
rank; the expected residual sum of squares is the decomposition of
[Section 9.1](../ch09-sums-of-squares/01-total-decomposition.html). Adding a knot
enlarges the column space, so the squared bias cannot increase while
\( \tr(\I-\M) \) falls by one.

(b) Suppose \( \X\mathbf{c}=\bzero \). By @thm-ply-spline-space the function
\( g=\sum_jc_j\phi_j \) built from the same coefficients is a spline vanishing at
every design point. On each piece \( g \) is a polynomial of degree at most
\( d \) with \( d+1 \) distinct roots there, hence identically zero there; so
\( g\equiv0 \) and, the basis being independent, \( \mathbf{c}=\bzero \). For the
second statement, \( (x_i-\kappa)_+^d=0 \) for every \( i \) when
\( \kappa>\max_ix_i \).

(c) For \( d=0 \) the columns are the orthogonal indicators of the intervals, so
\( \X\T\X=\diag(n_0,\dots,n_K) \) and the fit on interval \( k \) is \( \bar y_k \),
with variance \( \sigma^2/n_k \). Over counts summing to \( n \), the largest
\( \sigma^2/n_k \) is smallest when they are as equal as possible, which quantile
knots achieve up to rounding. For \( d\ge1 \) the variance formula is
@thm-ci-estimable-interval, and that it is large over a sparse interval is not
an identity but an observation.

(d) The first statement is @exr-ply-knot-coalescence. The bound on the basis
condition number is quoted from de Boor (2001) and not proved here; that it
controls \( \kappa(\X) \) needs the design to occupy every support, which is the
next paragraph.
:::

Part (c) is demonstrated in @exm-ply-knot-placement and part (d)
in @exm-ply-truncated-conditioning, whose nearly uniform design meets the
proviso; the local basis is built in
[Chapter 43](../ch43-smoothing/index.html).

Part (b) has a sharp form. In the local basis the design matrix has full rank
*if and only if* each basis function's support contains a design point — the
Schoenberg–Whitney condition (Schoenberg and Whitney, 1953), which is also the
proviso in (d). Either way the rule is: no knot where there are no data, and no
two knots with a nearly empty interval between them.

## Six bases for one trend

::: {#exm-ply-co2-spline}
[Trend bases for the carbon dioxide series]

Return to the \( n=521 \) monthly means with the four seasonal
columns fixed, put \( K=6 \) interior knots at the quantiles
\( k/7 \) of the time variable, and compare six trend bases.

| Trend basis | \( p \) | Residual SS | \( \hat\sigma \) | GCV | \( \max h_{ii} \) |
|---|---|---|---|---|---|
| quadratic | 7 | \( 297.03 \) | \( 0.7602 \) | \( 0.5858 \) | \( 0.0256 \) |
| degree 9 | 14 | \( 115.57 \) | \( 0.4774 \) | \( 0.2343 \) | \( 0.2054 \) |
| degree 16 | 21 | \( 92.44 \) | \( 0.4300 \) | \( 0.1926 \) | \( 0.4546 \) |
| piecewise constant | 11 | \( 3081.18 \) | \( 2.4580 \) | \( 6.1718 \) | \( 0.0215 \) |
| piecewise linear | 12 | \( 116.63 \) | \( 0.4787 \) | \( 0.2345 \) | \( 0.0532 \) |
| cubic spline | 14 | \( 108.48 \) | \( 0.4626 \) | \( 0.2199 \) | \( 0.1494 \) |

Read the middle rows against each other. At the same \( p=14 \) the cubic spline
fits better than the degree-9 polynomial (\( 108.48 \) against
\( 115.57 \)) and spreads its influence more evenly (largest
leverage \( 0.1494 \) against \( 0.2054 \));
panel (a) of [Figure 42.4.1](#fig-ply-compare) shows the two parting where the
polynomial has to bend. The step function is a disaster; the broken line is
almost as good as the cubic spline with two fewer parameters, at the cost of
visible kinks (panel (b)).

In-sample, GCV prefers the degree-16 polynomial to everything else
(\( 0.1926 \)) — honest about prediction at the design points
and misleading about the fit, since that model has
\( \max h_{ii}=0.4546 \) and
\( \kappa(\X)=5.39\times 10^{5} \). Among the piecewise bases
the cubic spline is the worst conditioned:
\( \kappa(\X)=7.72\times 10^{3} \) against
\( 9.75\times 10^{1} \).
:::

::: {when-format="html"}
![**Figure 42.4.1.** (a) The degree-9 polynomial trend and the cubic spline with six quantile knots,
both with ten trend parameters, as departures from the quadratic trend. (b) Three fits
on the same knots (vertical lines) over the last fifteen years, seasonally adjusted.](spline_comparison.svg){#fig-ply-compare width=100%}
:::

::: {when-format="pdf"}
![(a) The degree-9 polynomial trend and the cubic spline with six quantile knots,
both with ten trend parameters, as departures from the quadratic trend. (b) Three fits
on the same knots (vertical lines) over the last fifteen years, seasonally adjusted.](spline_comparison.pdf){width=100%}
:::

```{.python .run #cell-spline-comparison-bases}
import numpy as np
import statsmodels.api as sm

co2 = sm.datasets.co2.load_pandas().data["co2"]
monthly = co2.resample("MS").mean().dropna()
y = monthly.to_numpy()
t = (monthly.index.year + (monthly.index.month - 0.5) / 12).to_numpy()
n = len(y)
lo, hi = t.min(), t.max()
u = lambda s: 2 * (s - lo) / (hi - lo) - 1

K = 6
knots = np.quantile(t, np.arange(1, K + 1) / (K + 1))      # quantile knots
print("knots at", np.round(knots, 2))


def polynomial(s, d):
    return np.vander(u(s), d + 1, increasing=True)


def step(s):                                               # piecewise constant, d = 0
    piece = np.digitize(s, knots)
    return np.column_stack([(piece == k).astype(float) for k in range(K + 1)])


def broken_line(s):                                        # continuous piecewise linear
    return np.column_stack([np.ones_like(s), u(s)]
                           + [np.clip(u(s) - u(k), 0.0, None) for k in knots])


def cubic_spline(s):                                       # truncated power basis, d = 3
    return np.column_stack([u(s) ** j for j in range(4)]
                           + [np.clip(u(s) - u(k), 0.0, None) ** 3 for k in knots])


def with_season(f):
    return lambda s: np.column_stack([f(s), np.cos(2 * np.pi * s), np.sin(2 * np.pi * s),
                                      np.cos(4 * np.pi * s), np.sin(4 * np.pi * s)])


bases = {"quadratic": lambda s: polynomial(s, 2), "degree 9": lambda s: polynomial(s, 9),
         "degree 16": lambda s: polynomial(s, 16), "piecewise constant": step,
         "piecewise linear": broken_line, "cubic spline": cubic_spline}
```

```{.python .run #cell-spline-comparison-compare}
for name, f in bases.items():
    X = with_season(f)(t)
    p = np.linalg.matrix_rank(X)
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    resid = y - X @ beta
    rss = resid @ resid
    lev = np.einsum("ij,ij->i", X @ np.linalg.pinv(X.T @ X), X)
    print(f"{name:20s} p = {p:2d}   RSS = {rss:8.2f}   sigma = {np.sqrt(rss / (n - p)):.4f}"
          f"   GCV = {n * rss / (n - p) ** 2:.4f}   max leverage {lev.max():.4f}"
          f"   cond = {np.linalg.cond(X):.1e}")
```

Locality shows up in the same fits: moving the December 2001 value by \( 10 \)
ppm moves the March 1958 trend by \( 0.197 \) ppm under the
degree-12 polynomial but \( 0.004 \) under the cubic spline and
\( 0.001 \) under the step function, where only the shared
seasonal columns keep it from being zero.

## Where to put them

::: {#exm-ply-knot-placement}
[A skewed design]

Take \( n=200 \) design points drawn once, with a fixed seed,
from a skewed distribution on \( [0,1] \), and fit a cubic spline with
\( K=5 \) interior knots. Equally spaced knots leave
\( 67, 71, 35, 21, 5, 1 \) points in the six intervals; quantile
knots leave \( 34, 33, 33, 33, 33, 34 \). Both have full rank —
the truncated powers are globally supported, so a nearly empty interval does not
kill a column — but it does kill precision. The exact pointwise standard
deviation of the fit, in units of \( \sigma \), reaches
\( 7.452 \) with the equally spaced knots against
\( 0.987 \) with the quantile knots; at the median of the
design the quantile knots are slightly *worse*
(\( 0.173 \) against \( 0.123 \)), having
spent four of five knots there.

Quantile knots are a defensible default, not a universal optimum: they make the
worst case tolerable. Where the regression function is known to vary fast in a
sparse region the knots should follow the function — a choice made from outside
the data, or paid for as a selection (@prp-sel-selection-bias).
:::

::: {when-format="html"}
![**Figure 42.4.2.** (a) A skewed design of 200 points with five equally spaced and five quantile
knots. (b) Pointwise standard deviation of the fit, in units of \( \sigma \), log scale.](knot_placement.svg){#fig-ply-knots width=100%}
:::

::: {when-format="pdf"}
![(a) A skewed design of 200 points with five equally spaced and five quantile
knots. (b) Pointwise standard deviation of the fit, in units of \( \sigma \), log scale.](knot_placement.pdf){width=100%}
:::

```{.python .run #cell-knot-placement-design}
import numpy as np

rng = np.random.default_rng(42)
n = 200
x = np.sort(rng.beta(1.5, 6.0, n))                      # a skewed design on [0, 1]
K, d = 5, 3

equal = np.linspace(x.min(), x.max(), K + 2)[1:-1]      # equally spaced knots
quantile = np.quantile(x, np.arange(1, K + 1) / (K + 1))


def spline_design(s, knots):
    """Truncated power basis for degree d with the given interior knots."""
    return np.column_stack([s ** j for j in range(d + 1)]
                           + [np.clip(s - k, 0.0, None) ** d for k in knots])


def counts(knots):
    edges = np.r_[-np.inf, knots, np.inf]
    return [int(((x >= a) & (x < b)).sum()) for a, b in zip(edges[:-1], edges[1:])]


print("equally spaced knots:", np.round(equal, 3), "  points per interval", counts(equal))
print("quantile knots:      ", np.round(quantile, 3), "  points per interval", counts(quantile))
```

```{.python .run #cell-knot-placement-variance}
grid = np.linspace(x.min(), x.max(), 500)
for name, knots in (("equally spaced", equal), ("quantile", quantile)):
    X = spline_design(x, knots)
    G = np.linalg.inv(X.T @ X)
    Xg = spline_design(grid, knots)
    sd = np.sqrt(np.einsum("ij,jk,ik->i", Xg, G, Xg))   # in units of sigma
    print(f"{name:15s} knots: standard deviation of the fit at the median {np.interp(np.median(x), grid, sd):.3f},"
          f" largest {sd.max():.3f} at x = {grid[sd.argmax()]:.3f}")
```

::: {.idea}
Rules of thumb that survive scrutiny, for a fixed-knot regression spline:

- use a low degree, almost always cubic, and control flexibility with \( K \);
- place knots at quantiles unless something outside the data says otherwise,
  keeping several observations between consecutive knots;
- compare a few values of \( K \) by GCV or AIC rather than searching finely;
- report the fit as a curve with a band, never as coefficients, which depend on
  the basis and mean nothing alone.

The middle two leave a discrete, data-dependent choice.
[Chapter 43](../ch43-smoothing/index.html) removes it by taking \( K \) large and
adding a penalty, making flexibility continuous (@def-smo-pspline).
:::

## Why this basis is not the one to compute with

::: {#exm-ply-truncated-conditioning}
[Truncated powers against a local basis]

For the carbon dioxide design rescaled to \( [-1,1] \), take cubic splines with
\( K \) quantile knots and compare the condition number of the truncated power
basis with that of the local (B-spline) basis for the *same* space, verified by
checking that the two sets of columns together have rank \( K+4 \).

| \( K \) | 4 | 8 | 16 | 32 |
|---|---|---|---|---|
| truncated powers | \( 1.75\times 10^{3} \) | \( 2.56\times 10^{4} \) | \( 3.65\times 10^{5} \) | \( 4.97\times 10^{6} \) |
| local basis | \( 5.30 \) | \( 5.18 \) | \( 5.04 \) | \( 5.00 \) |

The truncated power condition number multiplies by about fifteen for every
doubling of \( K \); the local basis does not move. By \( K=32 \) the normal
equations (@thm-proj-normal-equations) in the truncated power basis have lost
about thirteen digits (@eq-cmp-rule), the same fit in the local basis one. This design is
nearly uniform, so every basis function's support is well occupied — the
proviso of part (d).
:::

The remedy is that of [Section 42.2](02-orthogonal-polynomials.html): the space is
fine and the basis is bad. The local basis is built in
[Chapter 43](../ch43-smoothing/index.html) (@def-smo-bspline
and @prp-smo-bspline). The truncated power basis keeps its place in the theory —
@thm-ply-spline-space is easiest to prove with it, and \( \theta_k \) means a
jump in the \( d \)th derivative — but it should not be what a program forms.

## Honest limits of a fitted curve

Three warnings close the chapter, none of them about arithmetic. *Outside the
range there is nothing*: a cubic spline continues as a cubic beyond its last
knot, and refitting to the months up to 1990 and predicting 1991–2001 gives root
mean squared errors of \( 2.71 \) ppm for the broken line,
\( 3.99 \) for the quadratic,
\( 46.11 \) for the cubic spline and
\( 923.15 \) for the degree-9 polynomial — almost the reverse
of the ranking inside the data.

*A curve is a conditional mean, not a mechanism*: dropping linearity removes a
functional-form assumption and nothing else, and @thm-cau-backdoor is still
what licenses a causal reading. And *a flexible fit is no substitute for a good
parameterization*: where the relationship is known to be a power law, the
transformations of [Chapter 22](../ch22-transformations/index.html) give a
two-parameter model that a ten-parameter spline can only approximate and that
extrapolates on a principle (@thm-tr-box-tidwell). A spline is often best used
diagnostically: fit it, look at the shape, then ask whether a parametric form
with meaning would do.

## Exercises

### A. Check your understanding

::: {#exr-ply-knot-df}
[A1]

A cubic regression spline with \( K \) quantile knots is fitted to \( n=200 \)
observations together with three other regressors and an intercept. Give the residual
degrees of freedom and the average leverage as functions of \( K \), and the \( K \) at
which the fit interpolates.
:::

::: {#exr-ply-which-basis}
[A2]

Three analysts fit a cubic spline with the same six knots, in the truncated power
basis, in the local basis, and as free cubics under the continuity constraints. Which
outputs agree exactly, which to rounding, and which disagree?
:::

### B. Practice

::: {#exr-ply-knot-test}
[B1]

Fit cubic splines to the carbon dioxide series with \( K=2,4,\dots,20 \) quantile knots
and plot GCV against \( K \), then the leave-one-out criterion computed by the leverage
shortcut of @thm-sel-loocv. Where do the two disagree, and why should you not report the
\( F \) test comparing the selected \( K \) with \( K=2 \)?
:::

::: {.solution}
The two have the same minimizer here or nearly so; GCV replaces each \( 1-h_{ii} \) by
its average \( 1-p/n \), which matters when the leverages are uneven, as at the
boundary. The \( F \) test is invalid because \( K \) minimized a criterion computed
from the same data: the residual mean square is biased downwards and the statistic is
not \( F \) (@prp-sel-selection-bias).
:::

::: {#exr-ply-band-honesty}
[B2]

Compute the pointwise \( 95\% \) band for the cubic spline fit of
@exm-ply-co2-spline and the simultaneous Working–Hotelling
band (@thm-ci-working-hotelling). By what factor is the simultaneous band wider, and
what does each band cover?
:::

::: {.solution}
The ratio of half-widths is \( \sqrt{pF_{\alpha}(p,n-p)}/t_{\alpha/2}(n-p) \), about two
and a half for \( p=14 \) and \( n=521 \). The pointwise band covers
\( \E(Y\mid x=t) \) at a single pre-chosen \( t \) with probability \( 0.95 \), the
simultaneous band the whole curve. Neither covers anything if the true mean is not in
the spline space: both are conditional on the model, and the model includes the knots.
:::

::: {#exr-ply-add-a-regressor}
[B3]

A model has a cubic spline in \( x \) and a linear term in \( z \). Show that the
coefficient of \( z \) is the slope of the regression of \( \y \) on the part of
\( \mathbf{z} \) orthogonal to the spline space, and say what it would mean for that part
to be nearly zero.
:::

::: {.solution}
This is @thm-proj-fwl with \( \X_2 \) the spline columns. If \( \mathbf{z} \) is nearly a
spline function of \( x \) — because \( z \) is a smooth transformation of it, or because
\( K \) is so large that the spline space reproduces almost anything — then
\( (\I-\M_2)\mathbf{z} \) is short, the coefficient's variance huge, and the two terms
not separately identified: the additive-model version of collinearity, taken up in
[Chapter 44](../ch44-additive-models/index.html) (@prp-add-concurvity).
:::

### C. Going deeper

::: {#exr-ply-equivalent-kernel}
[C1]

For a regression spline the fitted value at \( t \) is \( \sum_i\ell_i(t)y_i \) with
\( \ell_i(t)=\mathbf{x}(t)\T(\X\T\X)^{-1}\mathbf{x}_i \). Show that
\( \sum_i\ell_i(t)=1 \) and \( \sum_i\ell_i(t)(x_i-t)^{r}=0 \) for \( r=1,\dots,d \) at
every \( t \), because \( \mathcal S_d(\kappa) \) contains \( \mathcal P_d \) and a
projection reproduces its own space. Deduce that the fit reproduces polynomials of
degree \( d \) exactly, and interpret \( \sum_i\ell_i(t)^2 \).
:::

::: {#exr-ply-knot-selection-cost}
[C2]

Suppose \( K \) minimizes GCV over \( K\in\{1,\dots,20\} \) and the usual pointwise band
is then reported for the selected fit. Simulate \( 2000 \) data sets from a known smooth
\( f \) with normal errors and estimate the coverage of the nominal \( 95\% \) band at
the quartiles of the design, explaining the pattern in terms of @prp-sel-selection-bias
and of the bias a selected \( K \) leaves behind.
:::
