# Near replicates and nonparametric lack-of-fit tests

Exact replicates are a luxury of designed experiments. With continuous regressors, rows of the model
matrix are rarely identical: stack loss with all three operating variables has \( 21 \) rows and a single
repeated pair, one degree of freedom for pure error. The pure-error test survives in a more general
form: any larger model chosen without looking at the response gives an exact test, and so do regressors
built from the fitted values. What is lost is the guarantee that the larger model contains every
alternative: each test looks in some directions and is blind in others.

## A general recipe

::: {#prp-cor-augmented-test}
[Testing against any larger model]

Let \( \Y\sim\Normal_n(\boldsymbol{\uptheta},\sigma^2\I) \), let \( \X \) have rank \( r \), and let \( \Z \) be a fixed matrix, chosen without
reference to \( \Y \), with \( \C(\X)\subseteq\C(\Z) \) and \( \rank(\Z)=s \), \( r<s<n \). Let \( \M \) and \( \M_Z \) be the two projections and
\[
F=\frac{\bigl[\text{SSE}(\X)-\text{SSE}(\Z)\bigr]/(s-r)}{\text{SSE}(\Z)/(n-s)} .
\]

::: {.enumerate options="label=(\alph*)"}
1. If \( \boldsymbol{\uptheta}\in\C(\X) \), then \( F\sim F(s-r,n-s) \).

2. If \( \boldsymbol{\uptheta}\in\C(\Z) \), then \( F\sim F\bigl(s-r,n-s,\norm{(\I-\M)\boldsymbol{\uptheta}}^2/\sigma^2\bigr) \).

3. In general \( F=\dfrac{U/(s-r)}{V/(n-s)} \) with independent \( U\sim\chi^2(s-r,\gamma_1) \) and \( V\sim\chi^2(n-s,\gamma_2) \), where
   \( \gamma_1=\norm{(\M_Z-\M)\boldsymbol{\uptheta}}^2/\sigma^2 \) and \( \gamma_2=\norm{(\I-\M_Z)\boldsymbol{\uptheta}}^2/\sigma^2 \). For every \( c>0 \),
   \( \Pr(F>c) \) is strictly increasing in \( \gamma_1 \) and strictly decreasing in \( \gamma_2 \).
:::

:::

::: {.proof}
The decomposition \( \I-\M=(\M_Z-\M)+(\I-\M_Z) \) into orthogonal projections of ranks \( s-r \) and \( n-s \) comes from
@thm-proj-nested, and @thm-qf-orthogonal-projections gives (c) apart from the monotonicity. Parts (a) and (b)
are the cases \( \gamma_1=\gamma_2=0 \) and \( \gamma_2=0 \). For the monotonicity, write
\( \Pr(F>c)=\E\bigl[\Pr(U>c'V\mid V)\bigr] \) with \( c'=c(s-r)/(n-s) \). For each value \( V>0 \) the inner probability is
strictly increasing in \( \gamma_1 \) by @prp-qf-ncchisq-monotone. Similarly
\( \Pr(F>c)=\E\bigl[\Pr(V<U/c'\mid U)\bigr] \), and \( \Pr\{\chi^2(n-s,\gamma_2)<u\} \) is strictly decreasing in \( \gamma_2 \) for
\( u>0 \).
:::

The size is exact however \( \Z \) is chosen without using \( \Y \); with random regressors, \( \Z \) may be any
function of \( \X \) (@thm-cor-conditional). The
pure-error test is the case \( \Z= \) cell means. Part (c) describes what happens under the alternatives
that matter, where the true mean lies in neither model. Departures captured by the larger model raise
\( \gamma_1 \) and help; departures outside it raise \( \gamma_2 \), inflate the denominator and hurt. A departure
orthogonal to \( \C(\Z) \) gives power *below* the level (@exr-cor-lof-below-level).

## Near replicates, partitions and the rainbow test

Without exact replicates, the natural substitute is to group rows that are *nearly* identical. A
cluster analysis of the standardized rows of \( \X \) forms \( g \) clusters of near replicates, with indicator
matrix \( \mathbf{D} \). Testing \( \X \) against \( \Z=[\X,\mathbf{D}] \) allows each cluster its own shift
from the fitted plane. This is the near-replicate test studied by Christensen (1989); Shillington (1979) proposed a variant
based on the cluster averages of the regressors, and Christensen (1991) characterized the kinds of lack
of fit each detects. With exact replicates all reduce to the pure-error test.

A second family partitions the cases instead of the regressor space. If the rows are split into two
groups, with model matrices \( \X_1 \) and \( \X_2 \), the model that fits a separate regression to each group has
\( \Z=\begin{psmallmatrix}\X_1&\bzero\\\bzero&\X_2\end{psmallmatrix} \), whose column space contains \( \C(\X) \), and
\( \text{SSE}(\Z) \) is the sum of the two residual sums of squares (@exr-cor-partition-projection). A straight line
fitted to a curved mean is beaten by two lines, one for each half of the range. Utts (1982) proposed an
extreme partition, the **rainbow test**: fit the model to a central subset of \( m \) cases, those with the
smallest leverages, and give each of the remaining \( n-m \) cases a free mean of its own. Then \( \text{SSE}(\Z) \) is the
residual sum of squares of the central fit, and the test asks whether the model fitted to the middle of
the design extrapolates to its edges.

All of these are exact by @prp-cor-augmented-test, and all require a decision made from \( \X \) alone: how many
clusters, where to split, how large a central subset. Choosing the clusters or the split after looking at the
residuals destroys the exactness, because \( \Z \) then depends on \( \Y \).

## Richer mean functions

The most direct alternative is a larger parametric model: add squares and cross-products of the regressors,
or a regression spline. For a single regressor, a cubic spline with interior knots \( \kappa_1,\dots,\kappa_m \) adds the
columns \( x^2 \), \( x^3 \) and \( (x-\kappa_l)_+^3 \), \( l=1,\dots,m \), and can follow any smooth curve if the knots are dense
enough. These are reduced-model \( F \) tests (@thm-glh-f-test), exact by
@prp-cor-augmented-test(a) provided the larger model is fixed before seeing the data. Their power is concentrated on departures that the added columns can represent. A quadratic term is
excellent against smooth curvature and nearly useless against a departure confined to a small part of the
range.

## Regressors built from the fitted values

A classical alternative looks for departures along the fitted values: if the mean is a nonlinear
function of \( \x\T\bbeta \), then powers of \( \hat y \) should help predict \( y \). But \( \hat{\Y}=\M\Y \) depends on the response, so
@prp-cor-augmented-test does not apply. Remarkably, the \( F \) test is exact anyway.

::: {#thm-cor-fitted-regressors}
[Regressors that are functions of the fitted values]

Let \( \Y\sim\Normal_n(\X\bbeta,\sigma^2\I) \) with \( \rank(\X)=r \). Let \( \W=w(\M\Y) \) be an \( n\times q \) matrix computed from the
fitted values alone, with \( \rank[\X,\W]=r+q<n \) with probability one. Then
\[
F=\frac{\bigl[\text{SSE}(\X)-\text{SSE}([\X,\W])\bigr]/q}{\text{SSE}([\X,\W])/(n-r-q)}\ \sim\ F(q,n-r-q).
\]
:::

::: {.proof}
Let \( \tilde{\W}=(\I-\M)\W \), which has rank \( q \), and let \( \bP \) project onto \( \C(\tilde{\W}) \). As in
@lem-proj-fwl-split, \( \C([\X,\W])=\C(\X)\dirsum\C(\tilde{\W}) \) with orthogonal summands, so the projection onto
\( \C([\X,\W]) \) is \( \M+\bP \). With \( \he=(\I-\M)\Y \), and since \( \bP\M=\mathbf{0} \),
\[
\begin{aligned}
\text{SSE}(\X)-\text{SSE}([\X,\W])&=\norm{\bP\Y}^2=\norm{\bP\he}^2,\\
\text{SSE}([\X,\W])&=\norm{\he}^2-\norm{\bP\he}^2 .
\end{aligned}
\]
By @thm-qf-orthogonal-projections, \( \M\Y \) and \( \he \) are independent, and \( \he\sim\Normal_n(\bzero,\sigma^2(\I-\M)) \). So,
given \( \M\Y \), the matrix \( \bP \) is fixed and \( \he \) keeps its distribution, which is that of \( (\I-\M)\be^* \) for
\( \be^*\sim\Normal_n(\bzero,\sigma^2\I) \). Because \( \C(\bP)\subseteq\C(\X)\perpc \), \( \bP(\I-\M)=\bP \), and the two sums of squares
become \( \norm{\bP\be^*}^2 \) and \( \norm{(\I-\M-\bP)\be^*}^2 \). These are quadratic forms in orthogonal projections of
ranks \( q \) and \( n-r-q \), independent \( \sigma^2\chi^2(q) \) and \( \sigma^2\chi^2(n-r-q) \) variables by @thm-qf-orthogonal-projections. The conditional law of \( F \) is therefore \( F(q,n-r-q) \) whatever \( \M\Y \) is, and by
@lem-cor-conditioning it is also the unconditional law.
:::

The theorem is due to Milliken and Graybill (1970). Its best-known case is Ramsey's (1969) **RESET**,
which adds the squares and cubes of the fitted values, entry by entry. Tukey's (1949) one degree of
freedom for nonadditivity, which tests a two-way additive layout against a multiplicative interaction,
can be written in the same form with the squared fitted values. With a single regressor the construction
brings nothing new: when \( \hat{\beta}_1\ne0 \), the columns \( \hat{\y}^2 \) and \( \hat{\y}^3 \) span, together with \( \bone \) and \( \x \), the
same space as \( x^2 \) and \( x^3 \), so RESET is the cubic polynomial test (@exr-cor-reset-cubic). Its value is in
many dimensions, where it spends two degrees of freedom on the one direction along which the fit varies.

## Smoothing-based tests

Instead of committing to a parametric alternative, one can ask whether the residuals contain any smooth
structure. Let \( \bS \) be a fixed smoother matrix, for instance kernel weights that average each residual with
those of its neighbours in the regressor space, and consider
\[
T=\frac{\norm{\bS\he}^2}{\norm{\he}^2} .
\]{#eq-cor-smoother-statistic}

If the model is right, the residuals are noise and smoothing them leaves little; if the mean has structure
the model misses, the smoothed residuals keep it and \( T \) is large.

::: {#prp-cor-smoother-null}
[The smoother statistic is pivotal]

If \( \Y\sim\Normal_n(\X\bbeta,\sigma^2\I) \), the distribution of \( T \) in @eq-cor-smoother-statistic depends only on \( \X \) and \( \bS \),
not on \( \bbeta \) or \( \sigma^2 \). It equals the distribution of \( \norm{\bS(\I-\M)\mathbf{u}}^2/\norm{(\I-\M)\mathbf{u}}^2 \) with
\( \mathbf{u}\sim\Normal_n(\bzero,\I) \).
:::

::: {.proof}
\( \he=(\I-\M)\Y=(\I-\M)\be \), and \( \be=\sigma\mathbf{u} \) with \( \mathbf{u}\sim\Normal_n(\bzero,\I) \). The factor \( \sigma^2 \) cancels in the ratio.
:::

So the null distribution of \( T \) can be simulated as precisely as desired, from the design alone, and the
test has exact size apart from Monte Carlo error. Statistics of this kind were proposed by Azzalini and Bowman
(1993), who compare the residual sums of squares of the parametric fit and of a nonparametric smoother, and
they are developed at length by Hart (1997). The power depends on the amount of smoothing, and a bandwidth
chosen after inspecting the data breaks the exact size. With several regressors,
neighbourhoods become sparse, so in practice one smooths along a single direction, such as the fitted values.
And exactness depends, as for every test in this section, on normal homoscedastic errors.

## Stack loss with all three regressors

::: {#exm-cor-stackloss-full}
[Is the three-variable model adequate?]

Return to the stack loss data with air flow, water temperature and acid concentration, \( n=21 \) and
\( r=4 \). Ward clustering of the standardized rows into \( 8 \) clusters and testing \( \X \) against \( [\X,\mathbf{D}] \) gives
\( F=1.43 \) on \( 7 \) and \( 10 \) degrees of freedom (\( p \)-value \( 0.292 \)). Adding the squares of the three
regressors gives \( F=1.92 \) on \( 3 \) and \( 14 \) (\( p \)-value \( 0.172 \)). The rainbow test with the
\( 11 \) cases of smallest leverage as the central subset gives \( F=1.04 \) on \( 10 \) and \( 7 \)
(\( p \)-value \( 0.497 \)), and RESET gives \( F=1.17 \) on \( 2 \) and \( 15 \) (\( p \)-value
\( 0.336 \)). None of the tests finds lack of fit. The strong lack of fit of the line in air flow alone (@exm-cor-stackloss-lof) was the omitted water temperature, which varies within the groups of equal air
flow. With \( 21 \) cases none of these tests is powerful, and a few unusual days deserve the
scrutiny of [Chapter 20](../ch20-residuals-leverage-influence/index.html).
:::

```{.python .run #cell-lof-alternatives-stackloss}
import numpy as np
import statsmodels.api as sm
from scipy import stats
from scipy.cluster.hierarchy import fcluster, linkage
from statsmodels.stats.diagnostic import linear_reset

def proj(A):
    """Orthogonal projection onto C(A), via an orthonormal basis of the column space."""
    U, s, _ = np.linalg.svd(A, full_matrices=False)
    U = U[:, s > s[0] * 1e-10]
    return U @ U.T

def f_test(y, X, Zbig):
    """F test of C(X) against a larger C(Zbig) ⊇ C(X); returns F, df1, df2, p."""
    n = len(y)
    M, MZ = proj(X), proj(Zbig)
    r, s = round(np.trace(M)), round(np.trace(MZ))
    sse, sse_z = y @ (y - M @ y), y @ (y - MZ @ y)
    F = ((sse - sse_z) / (s - r)) / (sse_z / (n - s))
    return F, s - r, n - s, stats.f.sf(F, s - r, n - s)

data = sm.datasets.stackloss.load_pandas().data
y = data["STACKLOSS"].to_numpy()
R = data[["AIRFLOW", "WATERTEMP", "ACIDCONC"]].to_numpy()
n = len(y)
X = np.column_stack([np.ones(n), R])
fit = sm.OLS(y, X).fit()

# exact replicates: rows of X that occur more than once
_, rows = np.unique(R, axis=0, return_inverse=True)
print("distinct rows:", rows.max() + 1, "of", n)

# near replicates: cluster the standardized regressor rows (Ward linkage), c clusters
Rs = (R - R.mean(axis=0)) / R.std(axis=0)
clusters = fcluster(linkage(Rs, method="ward"), t=8, criterion="maxclust")
D = np.eye(clusters.max())[clusters - 1]              # cluster indicators
tests = {"near replicates": f_test(y, X, np.column_stack([X, D])),
         "added squares": f_test(y, X, np.column_stack([X, R ** 2]))}

# rainbow: refit on the half of the cases with the smallest leverages
h = np.diag(proj(X))
central = np.argsort(h)[: n // 2 + 1]
E = np.eye(n)[:, np.setdiff1d(np.arange(n), central)]  # one free mean per outer case
tests["rainbow"] = f_test(y, X, np.column_stack([X, E]))

# RESET: add squares and cubes of the fitted values
reset = linear_reset(fit, power=3, test_type="fitted", use_f=True)
tests["RESET"] = (reset.fvalue, reset.df_num, reset.df_denom, reset.pvalue)
for name, (F, d1, d2, p) in tests.items():
    print(f"{name:16s} F = {F:6.2f} on ({d1:.0f}, {d2:.0f}) df,  p = {p:.3f}")
```

## Which test?

No lack-of-fit test is best against every alternative: by @prp-cor-augmented-test(c), a test gains
power from departures its larger model can represent and loses power from the others.

::: {#exm-cor-lof-power}
[Two departures, five tests]

Take \( n=50 \) equally spaced values of \( x \) in \( (0,1) \), a straight-line working model, standard normal errors,
and five tests at level \( 0.05 \): adding \( x^2 \); RESET with squared and cubed fitted values; a cubic spline with
knots at \( 0.25 \), \( 0.5 \) and \( 0.75 \); near replicates formed by \( 10 \) clusters of \( 5 \) neighbouring
values; and the smoother statistic @eq-cor-smoother-statistic with Gaussian kernel weights of bandwidth
\( 0.06 \), its critical value simulated from \( 50000 \) null residual vectors. Over \( 4000 \) samples
from the straight line, the rejection rates are \( 0.050 \), \( 0.054 \), \( 0.054 \),
\( 0.052 \) and \( 0.054 \): all five tests are exact.

Against smooth curvature, a mean of \( 1+2x+6(x-0.5)^2 \), the powers are \( 0.874 \) for the
quadratic term, \( 0.786 \) for RESET (which here is the cubic polynomial test), \( 0.606 \) for the
spline, \( 0.422 \) for near replicates and \( 0.723 \) for the smoother. Against a narrow bump
of height \( 2.5 \) near \( x=0.8 \), the order is reversed: the smoother has power \( 0.622 \), the spline
\( 0.501 \), near replicates \( 0.459 \), RESET \( 0.227 \), and the quadratic term
\( 0.029 \), *below* the level: the bump is nearly orthogonal to \( x^2 \) after the line is removed, so it
inflates the denominator of that test without feeding its numerator. [Figure 14.7.1](#fig-cor-lof-power)
shows the power curves.
:::

::: {when-format="html"}
![**Figure 14.7.1.** Power of five exact lack-of-fit tests of a straight line (\( n=50 \), level \( 0.05 \)).
(a) The two departures. (b) Smooth curvature. (c) A local bump, against which the quadratic test has
power below its level.](lof_power.svg){#fig-cor-lof-power width=100%}
:::

::: {when-format="pdf"}
![Power of five exact lack-of-fit tests of a straight line (\( n=50 \), level \( 0.05 \)).
(a) The two departures. (b) Smooth curvature. (c) A local bump, against which the quadratic test has
power below its level.](lof_power.pdf){width=100%}
:::

```{.python .run #cell-lof-alternatives-power}
rng = np.random.default_rng(1471)
n2, reps = 50, 4000
x = (np.arange(n2) + 0.5) / n2                         # fixed design on (0, 1)
X2 = np.column_stack([np.ones(n2), x])
M2 = proj(X2)
knots = (0.25, 0.5, 0.75)
spline = np.column_stack([x ** 2, x ** 3] + [np.clip(x - k, 0, None) ** 3 for k in knots])
cells = np.eye(10)[np.arange(n2) // 5]                 # 10 clusters of 5 neighbouring x
bigger = {"quadratic": np.column_stack([X2, x ** 2]),
          "spline": np.column_stack([X2, spline]),
          "near replicates": np.column_stack([X2, cells])}
proj_big = {k: proj(v) for k, v in bigger.items()}

# kernel smoother of the residuals; T = |S e|^2 / |e|^2, null law free of beta and sigma
bw = 0.06
K = np.exp(-0.5 * ((x[:, None] - x[None, :]) / bw) ** 2)
S = K / K.sum(axis=1, keepdims=True)
null_e = (np.eye(n2) - M2) @ rng.normal(size=(n2, 50_000))
T_null = np.sum((S @ null_e) ** 2, axis=0) / np.sum(null_e ** 2, axis=0)
T_crit = np.quantile(T_null, 0.95)

def powers(mean):
    """Rejection rates at level 0.05 of each test when E(Y) = mean."""
    Y = mean[:, None] + rng.normal(size=(n2, reps))
    e = Y - M2 @ Y
    sse = np.sum(e ** 2, axis=0)
    out = {}
    for name, P in proj_big.items():
        s = round(np.trace(P))
        sse_z = np.sum((Y - P @ Y) ** 2, axis=0)
        F = ((sse - sse_z) / (s - 2)) / (sse_z / (n2 - s))
        out[name] = np.mean(F > stats.f.ppf(0.95, s - 2, n2 - s))
    Fr = np.empty(reps)
    for j in range(reps):                              # RESET: regressors depend on y
        yh = M2 @ Y[:, j]
        W = np.column_stack([yh ** 2, yh ** 3])
        Wt = W - M2 @ W
        g, *_ = np.linalg.lstsq(Wt, e[:, j], rcond=None)
        ss_add = e[:, j] @ (Wt @ g)
        Fr[j] = (ss_add / 2) / ((sse[j] - ss_add) / (n2 - 4))
    out["RESET"] = np.mean(Fr > stats.f.ppf(0.95, 2, n2 - 4))
    T = np.sum((S @ e) ** 2, axis=0) / sse
    out["smoother"] = np.mean(T > T_crit)
    return out

curve = (x - 0.5) ** 2                                 # smooth curvature
bump = np.exp(-0.5 * ((x - 0.8) / 0.05) ** 2)          # a local bump
sizes = powers(1 + 2 * x)
print("sizes:", {k: round(float(v), 3) for k, v in sizes.items()})
reported = {}
for label, shape, a in [("curvature", curve, 6.0), ("bump", bump, 2.5)]:
    reported[label] = powers(1 + 2 * x + a * shape)
    print(f"{label:9s}:", {k: round(float(v), 3) for k, v in reported[label].items()})
```

The lesson is to decide which inadequacy matters before testing: a low-degree polynomial or RESET for
smooth curvature, a flexible alternative or a smoother for local departures. Only the pure-error test of
@thm-cor-lack-of-fit has a larger model containing every mean function of the regressors. Residual plots,
which look in all directions informally, remain indispensable ([Chapter 20](../ch20-residuals-leverage-influence/index.html)).

## Exercises

### A. Check your understanding

::: {#exr-cor-near-df}
[A1]

A regression with an intercept and two regressors has \( n=30 \) cases, whose rows are grouped into \( 8 \) clusters
of near replicates. Assuming no combination of cluster indicators lies in \( \C(\X) \) except the constant,
give the degrees of freedom of the near-replicate test of \( \X \) against \( [\X,\mathbf{D}] \).
:::

### B. Practice

::: {#exr-cor-partition-projection}
[B1]

Let \( \X=\begin{psmallmatrix}\X_1\\\X_2\end{psmallmatrix} \) and \( \Z=\begin{psmallmatrix}\X_1&\bzero\\\bzero&\X_2\end{psmallmatrix} \). Show
that \( \C(\X)\subseteq\C(\Z) \), that the projection onto \( \C(\Z) \) is \( \begin{psmallmatrix}\M_1&\bzero\\\bzero&\M_2\end{psmallmatrix} \), and that
\( \text{SSE}(\Z) \) is the sum of the residual sums of squares of the two separate regressions.
:::

::: {.solution}
\( \X=\Z\begin{psmallmatrix}\I\\\I\end{psmallmatrix} \), so \( \C(\X)\subseteq\C(\Z) \). The block-diagonal matrix of the two projections is symmetric
and idempotent, and its column space is \( \C(\X_1)\times\C(\X_2)=\C(\Z) \), so it is the projection onto \( \C(\Z) \) (@thm-proj-sym-idem). Applying \( \I \) minus it to \( \y=(\y_1\T,\y_2\T)\T \) gives the two residual vectors
\( (\I-\M_1)\y_1 \) and \( (\I-\M_2)\y_2 \), whose squared lengths add.
:::

::: {#exr-cor-rainbow-df}
[B2]

In the rainbow test, order the cases so that the central subset of size \( m \) comes first, and let its model
matrix \( \X_c \) have rank \( r \). Show that \( \Z=\begin{psmallmatrix}\X_c&\bzero\\\bzero&\I_{n-m}\end{psmallmatrix} \) has
\( \C(\X)\subseteq\C(\Z) \), that \( \text{SSE}(\Z) \) is the residual sum of squares of the central fit, and that the test has
\( n-m \) and \( m-r \) degrees of freedom.
:::

::: {#exr-cor-lof-below-level}
[B3]

In @prp-cor-augmented-test, suppose \( \boldsymbol{\uptheta}=\X\bbeta+\boldsymbol{\updelta} \) with \( \boldsymbol{\updelta}\ne\bzero \) orthogonal to \( \C(\Z) \). Show that the
power of the test is strictly less than its size.
:::

::: {.solution}
Here \( (\M_Z-\M)\boldsymbol{\uptheta}=\M_Z\boldsymbol{\updelta}-\M\boldsymbol{\updelta}=\bzero \), since \( \boldsymbol{\updelta}\perp\C(\Z)\supseteq\C(\X) \), so \( \gamma_1=0 \), while
\( \gamma_2=\norm{\boldsymbol{\updelta}}^2/\sigma^2>0 \). By part (c), \( \Pr(F>c) \) is strictly smaller than its value at \( \gamma_2=0 \), which is the
size by part (a).
:::

::: {#exr-cor-reset-cubic}
[B4]

With a single regressor and an intercept, show that if \( \hat{\beta}_1\ne0 \), then
\( \C([\bone,\x,\hat{\y}^2,\hat{\y}^3])=\C([\bone,\x,\x^2,\x^3]) \), powers taken entry by entry. Conclude that RESET is then
the \( F \) test of the cubic polynomial.
:::

### C. Going deeper

::: {#exr-cor-tukey-nonadditivity}
[C1]

In a two-way layout with one observation per cell and the additive model \( \E(y_{ij})=\mu+\alpha_i+\beta_j \), Tukey's
statistic tests the added regressor \( (\hat\alpha_i\hat\beta_j) \). Show that, after the additive model is projected out,
this regressor spans the same space as the entrywise square of the fitted values, and deduce from
@thm-cor-fitted-regressors that Tukey's test has an exact \( F(1,(a-1)(b-1)-1) \) null distribution.
:::

