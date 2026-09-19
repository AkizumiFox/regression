# The multiple linear regression model

A response rarely depends on one regressor alone. The money a country's residents hold
depends on their income and on interest rates. The number of physician visits depends on
health, age and insurance. Even with one measured regressor, a straight line is often
too rigid: the curvature of a trend or the shape of a seasonal cycle needs several
columns to describe. Matrix notation handles all of these cases at once.

## Matrix form

Suppose case \( i \) has regressor values that have been turned into \( p \) numbers
\( x_{i1},\dots,x_{ip} \), and consider
\[
Y_i=\beta_1x_{i1}+\beta_2x_{i2}+\dots+\beta_px_{ip}+\varepsilon_i,\qquad i=1,\dots,n .
\]
Stack the responses, the coefficients and the errors into vectors, and the
regressor values into an \( n\times p \) matrix:
\[
\Y=\begin{pmatrix}Y_1\\\vdots\\Y_n\end{pmatrix},\qquad
\bbeta=\begin{pmatrix}\beta_1\\\vdots\\\beta_p\end{pmatrix},\qquad
\be=\begin{pmatrix}\varepsilon_1\\\vdots\\\varepsilon_n\end{pmatrix},
\]
\[
\X=\begin{pmatrix}x_{11}&\cdots&x_{1p}\\\vdots&&\vdots\\x_{n1}&\cdots&x_{np}\end{pmatrix}
=\begin{pmatrix}\x_{(1)}\T\\\vdots\\\x_{(n)}\T\end{pmatrix}.
\]
The \( n \) equations become the single equation \( \Y=\X\bbeta+\be \). Row \( i \) of \( \X \),
written \( \x_{(i)}\T \), holds the regressor values of case \( i \), its **design point**.
Column \( j \), written \( \x_j \), holds the values of the \( j \)th regressor across the cases.

::: {#def-lm-linear-model}
[The linear model]

The **linear model** for a random vector \( \Y\in\Real^n \) is
\[
\Y=\X\bbeta+\be,\qquad \E(\be)=\bzero,\qquad \Cov(\be)=\sigma^2\I_n,
\]{#eq-lm-model}

where the **model matrix** \( \X \) is a known \( n\times p \) matrix of constants, and the
**coefficient vector** \( \bbeta\in\Real^p \) and the **error variance** \( \sigma^2>0 \) are
unknown parameters. The conditions on \( \be \) are the **second-moment assumptions**.
If in addition \( \be\sim\Normal_n(\bzero,\sigma^2\I_n) \), equivalently
\( \Y\sim\Normal_n(\X\bbeta,\sigma^2\I_n) \), the model is the **normal linear model**.
The model has **full rank** if \( \rank(\X)=p \).
:::

The second-moment assumptions are exactly the layers (L1)–(L3) of
[Section 5.1](01-what-a-model-claims.html). By @thm-rv-linear they say
\[
\E(\Y)=\X\bbeta,\qquad \Cov(\Y)=\sigma^2\I_n .
\]{#eq-lm-model-moments}

The first equation is the whole content of (L1). It says the mean vector
\( \E(\Y) \) is a linear combination of the columns of \( \X \), with unknown coefficients.
The second combines constant variance and uncorrelated errors. When the regressors are
random, both equations are read conditionally on \( \X \), as explained in
[Section 5.1](01-what-a-model-claims.html). The normal version is a model for the full
distribution of \( \Y \), in the sense of [Chapter 3](../ch03-multivariate-normal/index.html).
Since a normal vector with covariance \( \sigma^2\I \) has independent components
(@thm-mvn-independence), the errors are then independent and not merely uncorrelated.

The definition does not require full rank; most results in this chapter assume it, which
forces \( p\le n \), and [Chapter 6](../ch06-projections/index.html) and [Chapter 8](../ch08-estimability/index.html) treat the general case. What must be *known* is
\( \X \), and what must be *linear* is the dependence of \( \E(\Y) \) on \( \bbeta \). Nor does the
model say that the errors are identically distributed, only that they have a common variance.

## The intercept

Almost every model in practice includes a column of ones, \( \x_1=\bone \), so that
\[
\E(Y_i)=\beta_0+\beta_1x_{i1}+\dots+\beta_kx_{ik},\qquad p=k+1 .
\]
In this form it is conventional to number the coefficients from \( 0 \), and to call
\( k \) the number of regressors. The coefficient of the column of ones is the
**intercept**. Including it lets the fitted surface find its own level, rather than being
forced through the origin, and it has two further consequences that
[Section 5.4](04-least-squares.html) proves: the residuals sum to zero, and the slopes
can be computed from centred data. Leaving out the intercept is justified only when
theory requires the mean response to vanish at the origin of every regressor, as in the
regression through the origin at the end of [Section 5.2](02-simple-regression.html).

## Linear in the parameters, not in the regressors

The word *linear* in "linear model" refers to \( \bbeta \). The columns of \( \X \) can
be any known functions of the measured variables. This turns a remarkably wide class of
curved relationships into linear models.

**Transformations.** If the mean response is thought to depend on \( \log z \), the
column \( \log z_i \) enters \( \X \), and \( \E(Y_i)=\beta_0+\beta_1\log z_i \) is linear in
\( (\beta_0,\beta_1) \). The response can be transformed too, as in the money demand
model of [Section 5.4](04-least-squares.html), where \( Y \) is the logarithm of money
holdings. A transformation of the response changes the meaning of the error, and hence
the assumptions; Chapter 22 treats this.

**Polynomials.** A quadratic trend \( \E(Y_i)=\beta_0+\beta_1t_i+\beta_2t_i^2 \) uses the
columns \( \bone \), \( \mathbf{t} \) and \( \mathbf{t}^2 \) (squared entrywise). A polynomial of degree
\( d \) in one variable needs \( d+1 \) columns.

**Periodic terms.** A seasonal cycle of period one year can be written as
\( a\cos(2\pi t+\phi) \) with unknown amplitude \( a \) and phase \( \phi \). This is not
linear in \( \phi \). But
\( a\cos(2\pi t+\phi)=\gamma_1\cos(2\pi t)+\gamma_2\sin(2\pi t) \) with
\( \gamma_1=a\cos\phi \) and \( \gamma_2=-a\sin\phi \), and this *is* linear in
\( (\gamma_1,\gamma_2) \). The amplitude is recovered as \( a=(\gamma_1^2+\gamma_2^2)^{1/2} \).

::: {#exm-lm-co2}
[A curved trend and a seasonal cycle]

The monthly mean concentration of carbon dioxide in the air at Mauna Loa, Hawaii, from
March 1958 to December 2001, gives \( n=521 \) observations in parts per million.
Let \( t \) be the time in years, at the middle of each month, and \( t_c=t-1980 \). The model
\[
\begin{aligned}
\E(Y)={}&\beta_0+\beta_1t_c+\beta_2t_c^2\\
&+\gamma_1\cos 2\pi t+\gamma_2\sin 2\pi t+\gamma_3\cos 4\pi t+\gamma_4\sin 4\pi t
\end{aligned}
\]
has \( p=7 \) columns, and the model matrix has full rank. Its least
squares fit ([Section 5.4](04-least-squares.html)) has
\[
\hat{\beta}_0=337.6163,\quad\hat{\beta}_1=1.3344,\quad\hat{\beta}_2=0.0118,
\]
\[
\hat{\gamma}_1=-1.0413,\quad\hat{\gamma}_2=2.5641,\quad
\hat{\gamma}_3=0.6284,\quad\hat{\gamma}_4=-0.3859 .
\]
The quadratic term says the rate of increase is itself increasing. The slope of the
trend, \( \beta_1+2\beta_2t_c \), is estimated as \( 0.86 \) ppm per year in 1960
and \( 1.81 \) ppm per year in 2000. The annual harmonic has amplitude
\( 2.77 \) ppm and the half-year harmonic \( 0.74 \) ppm. Together
they make a seasonal swing of \( 6.0 \) ppm from the late-spring peak to
the autumn trough, as the vegetation of the northern hemisphere takes up carbon
over the summer. The residual standard deviation (estimated as in
[Section 5.6](06-error-variance.html)) is \( 0.760 \) ppm, against
\( 2.743 \) ppm for a straight line in time alone.
[Figure 5.3.1](#fig-lm-co2) shows the data, the trend and the seasonal cycle.
:::

::: {when-format="html"}
![**Figure 5.3.1.** Monthly carbon dioxide at Mauna Loa, 1958–2001. (a) The data and the
fitted quadratic trend. (b) The data with the fitted trend removed, by calendar month,
and the fitted seasonal cycle made of two harmonics. The model is curved in time and in
the calendar month, but linear in its seven coefficients.](co2_design.svg){#fig-lm-co2 width=100%}
:::

::: {when-format="pdf"}
![Monthly carbon dioxide at Mauna Loa, 1958–2001. (a) The data and the
fitted quadratic trend. (b) The data with the fitted trend removed, by calendar month,
and the fitted seasonal cycle made of two harmonics. The model is curved in time and in
the calendar month, but linear in its seven coefficients.](co2_design.pdf){width=100%}
:::

```{.python .run #cell-co2-design-design}
import numpy as np
import statsmodels.api as sm
co2 = sm.datasets.co2.load_pandas().data["co2"]
monthly = co2.resample("MS").mean().dropna()        # weekly values -> monthly means
y = monthly.to_numpy()
t = (monthly.index.year + (monthly.index.month - 0.5) / 12).to_numpy()
tc = t - 1980.0                                      # years since 1980

X = np.column_stack([
    np.ones_like(tc),                                # intercept
    tc, tc**2,                                       # quadratic trend
    np.cos(2 * np.pi * t), np.sin(2 * np.pi * t),    # annual cycle
    np.cos(4 * np.pi * t), np.sin(4 * np.pi * t),    # half-year harmonic
])
n, p = X.shape
print("n =", n, " p =", p, " rank =", np.linalg.matrix_rank(X))
```

```{.python .run #cell-co2-design-fit}
beta_hat = np.linalg.solve(X.T @ X, X.T @ y)         # normal equations (fine here)
fitted = X @ beta_hat
resid = y - fitted
print("coefficients:", np.round(beta_hat, 4))
print(f"residual standard deviation {np.sqrt(resid @ resid / (n - p)):.3f} ppm")

beta_line = np.linalg.solve(X[:, :2].T @ X[:, :2], X[:, :2].T @ y)
print(f"straight line only: slope {beta_line[1]:.3f} ppm per year")
```

The residuals of this model are strongly correlated from month to month: the
correlation between consecutive residuals is \( 0.90 \), because departures from a
smooth trend persist for many months. So (L3) fails and the standard errors of
[Section 5.5](05-properties.html) would be too optimistic. The *fit* is still a
sensible description. Chapter 21 treats serially correlated errors.

**Models that are not linear.** Some relations cannot be rearranged into linear form.
Exponential decay toward an unknown asymptote,
\( \E(Y)=\beta_0+\beta_1e^{-\beta_2x} \), is nonlinear in \( \beta_2 \), and so is a cycle
of unknown period. Such models are fitted by iterative nonlinear least squares, which
this book does not treat. Others become linear after a transformation: if
\( Y=\alpha x^{\beta}\eta \) with a positive multiplicative error \( \eta \), then
\( \log Y=\log\alpha+\beta\log x+\log\eta \) is a linear model for \( \log Y \), provided the
assumptions are made about \( \log\eta \) rather than \( \eta \).

## Categorical regressors and indicator variables

A regressor that records a category, such as treatment group, region, or a
self-rating on a four-point scale, cannot enter \( \X \) as a single column of codes
\( 1,2,3,4 \). That would force the mean differences between consecutive categories to be
equal, which is an arbitrary consequence of the labels. Instead, each category except
one gets an **indicator variable** (or *dummy variable*): a column that is \( 1 \) for
cases in that category and \( 0 \) otherwise. The omitted category is the **reference
category**.

::: {#exm-lm-health-coding}
[Self-rated health]

The RAND Health Insurance Experiment recorded, for 20,190
person-years, the number of outpatient visits to a physician and the person's own rating
of their health: excellent, good, fair or poor. With excellent as the reference, the model
\[
\E(Y_i)=\beta_0+\beta_1g_i+\beta_2f_i+\beta_3q_i
\]
uses three indicator columns \( \mathbf{g},\mathbf{f},\mathbf{q} \) for good, fair and poor health. For a
person in excellent health the mean is \( \beta_0 \); in good health, \( \beta_0+\beta_1 \); in
fair health, \( \beta_0+\beta_2 \); in poor health, \( \beta_0+\beta_3 \). So \( \beta_1,\beta_2,\beta_3 \)
are the differences in mean visits between each category and the reference.
@prp-lm-indicator-means shows that the least squares estimates are the corresponding
differences of sample means. In these data they are \( 0.268 \),
\( 1.058 \) and \( 3.160 \) visits a year, on top of
\( \hat{\beta}_0=2.634 \) for excellent health.
:::

Why omit a category? With an indicator for *every* category the columns satisfy
\( \mathbf{e}+\mathbf{g}+\mathbf{f}+\mathbf{q}=\bone \), where \( \mathbf{e} \) indicates excellent health. The
model matrix \( [\bone,\mathbf{e},\mathbf{g},\mathbf{f},\mathbf{q}] \) then has five columns and rank four, so it
does not have full rank. Such a model is not wrong. It describes the same set of mean
vectors, but with one more parameter than it can identify, and the least squares
coefficients are no longer unique. This situation is sometimes called the *dummy
variable trap*. It is not a trap at all once the right theory is in place:
[Section 6.4](../ch06-projections/04-least-squares.html) shows that fitted values are
unaffected, and [Chapter 8](../ch08-estimability/index.html) shows which functions of the
coefficients, such as the difference between two categories, are still estimated
uniquely.

```{.python .run #cell-rand-dummies-trap}
import numpy as np
import statsmodels.api as sm
rand = sm.datasets.randhie.load_pandas().data
y = rand["mdvis"].to_numpy(dtype=float)
D = rand[["hlthg", "hlthf", "hlthp"]].to_numpy(dtype=float)   # good, fair, poor
X = np.column_stack([np.ones(len(y)), D])                      # excellent = reference

excellent = 1.0 - D.sum(axis=1)
X_all = np.column_stack([X, excellent])          # an indicator for every level, and 1
print("columns:", X_all.shape[1], " rank:", np.linalg.matrix_rank(X_all))
print("1 - (sum of the four indicators) =", np.abs(X_all[:, 0] - X_all[:, 1:].sum(axis=1)).max())
```

Reference coding is one choice among several. **Effect coding** replaces the indicator
of category \( k \) by a column that is \( 1 \) in category \( k \), \( -1 \) in the reference
category, and \( 0 \) elsewhere. Its coefficients measure deviations of category means from
their unweighted average rather than from a reference. The fitted values are the same
under either coding, because the two model matrices have the same column space
(@prp-lm-reparameterization). Only the meaning of the coefficients changes.

## Interactions

In the model \( \E(Y)=\beta_0+\beta_1x_1+\beta_2x_2 \) the slope in \( x_1 \) is \( \beta_1 \)
whatever the value of \( x_2 \): the effects are **additive**. An **interaction** lets the
slope in one variable depend on another. Adding the product column \( x_1x_2 \),
\[
\E(Y)=\beta_0+\beta_1x_1+\beta_2x_2+\beta_3x_1x_2=(\beta_0+\beta_2x_2)+(\beta_1+\beta_3x_2)x_1 ,
\]{#eq-lm-interaction}

gives a slope in \( x_1 \) of \( \beta_1+\beta_3x_2 \). The model is still linear in the
coefficients. When \( x_2 \) is an indicator for a group, @eq-lm-interaction fits a separate
straight line in each group: intercept \( \beta_0 \) and slope \( \beta_1 \) in the reference
group, intercept \( \beta_0+\beta_2 \) and slope \( \beta_1+\beta_3 \) in the other. When both
variables are categorical, products of indicators give each combination of categories
its own mean, the two-way layout of Chapter 16.

Interactions change what the other coefficients mean. In @eq-lm-interaction, \( \beta_1 \)
is the slope in \( x_1 \) *when \( x_2=0 \)*, which may be far outside the data.
Centring \( x_2 \) before forming the product makes \( \beta_1 \) the slope at the average
\( x_2 \) instead. [Section 5.7](07-interpreting-coefficients.html) returns to this.

## When is the model matrix of full rank?

Full rank means the columns of \( \X \) are linearly independent: no column is an exact
linear combination of the others. By rank–nullity (@thm-mat-rank-nullity), this holds
iff \( \X\bb=\bzero \) only for \( \bb=\bzero \), that is, iff two different coefficient vectors
never give the same mean vector. It fails in three common ways:

1. **Too few cases**: \( n<p \). Chapter 28 deals with this setting.

2. **Redundant coding**: an indicator for every category together with an intercept, as
   above, or the same quantity measured in two units.

3. **Unlucky data**: a design that happens to make columns dependent, such as a
   quadratic in a variable observed at only two distinct values, or an interaction
   between two indicators when some combination of categories never occurs.

Near-failures, where columns are almost dependent, are more common and more
troublesome. They make coefficients imprecise (Chapter 26) and computations inaccurate
([Chapter 10](../ch10-computation/index.html)).

## Exercises

### A. Check your understanding

::: {#exr-lm-write-X}
[A1]

Write out the model matrix for each model, for the data \( x=(1,2,3,4) \) and a group label
\( (A,A,B,B) \), and give its rank: (a) a quadratic in \( x \); (b) separate intercepts for the
two groups and a common slope; (c) separate lines in the two groups, with reference coding;
(d) an intercept, an indicator for \( A \) and an indicator for \( B \).
:::

::: {.solution}
With \( g \) the indicator of group \( B \):
(a) rows \( (1,x,x^2) \): \( (1,1,1),(1,2,4),(1,3,9),(1,4,16) \), rank \( 3 \).
(b) rows \( (1,g,x) \): \( (1,0,1),(1,0,2),(1,1,3),(1,1,4) \), rank \( 3 \).
(c) rows \( (1,g,x,gx) \): \( (1,0,1,0),(1,0,2,0),(1,1,3,3),(1,1,4,4) \), rank \( 4 \); it is
square and nonsingular, so the two lines pass through the data exactly.
(d) rows \( (1,1,0),(1,1,0),(1,0,1),(1,0,1) \), rank \( 2 \): the first column is the sum of the
other two.
:::

::: {#exr-lm-linear-or-not}
[A2]

Which of these mean functions define linear models in the stated parameters (after a
change of parameters if necessary)? (a) \( \beta_0+\beta_1\sqrt{x} \); (b)
\( \beta_0e^{\beta_1x} \); (c) \( \beta_0+\beta_1x+\beta_2\max(x-5,0) \); (d)
\( \beta_0+\beta_1\sin(x-\beta_2) \); (e) \( 1/(\beta_0+\beta_1x) \).
:::

### B. Practice

::: {#exr-lm-effect-coding}
[B1]

For a factor with three levels and an intercept, write the model matrix rows under
reference coding (level 3 as reference) and under effect coding. Find the nonsingular
\( 3\times3 \) matrix \( \mathbf{K} \) with \( \X_{\text{effect}}=\X_{\text{ref}}\mathbf{K} \), and express the
effect-coding coefficients in terms of the three level means \( \mu_1,\mu_2,\mu_3 \).
:::

::: {.solution}
Reference coding has rows \( (1,1,0) \), \( (1,0,1) \), \( (1,0,0) \) for levels 1, 2, 3.
Effect coding has rows \( (1,1,0) \), \( (1,0,1) \), \( (1,-1,-1) \). Each effect-coded row is the
reference-coded row times
\( \mathbf{K}=\begin{psmallmatrix}1&-1&-1\\0&2&1\\0&1&2\end{psmallmatrix} \): for instance
\( (1,0,0)\mathbf{K}=(1,-1,-1) \) and \( (1,1,0)\mathbf{K}=(1,1,0) \). \( \det\mathbf{K}=3\ne0 \). Under effect
coding the level means are \( \mu_1=\alpha_0+\alpha_1 \), \( \mu_2=\alpha_0+\alpha_2 \),
\( \mu_3=\alpha_0-\alpha_1-\alpha_2 \). Adding, \( \alpha_0=(\mu_1+\mu_2+\mu_3)/3 \), and then
\( \alpha_k=\mu_k-\alpha_0 \) for \( k=1,2 \): deviations from the unweighted average of the means.
:::

::: {#exr-lm-interaction-centring}
[B2]

In @eq-lm-interaction, replace \( x_2 \) by \( x_2-c \) in both the main-effect column and the
product. Find the new coefficients in terms of the old, and show that the fitted mean
function is unchanged. For which \( c \) is the new coefficient of \( x_1 \) the slope at the
average value of \( x_2 \)?
:::

### C. Going deeper

::: {#exr-lm-harmonic-orthogonal}
[C1]

Let \( t_i=i/m \) for \( i=0,\dots,m-1 \) (one full period sampled at \( m \) equally spaced
points, with \( m\ge 5 \)). Show that the columns \( \bone,\cos2\pi\mathbf{t},\sin2\pi\mathbf{t},
\cos4\pi\mathbf{t},\sin4\pi\mathbf{t} \) are mutually orthogonal, and find their squared lengths.
Conclude that for a whole number of periods the harmonic coefficients are estimated
separately from each other and from the mean. Why is this not exactly the case in
@exm-lm-co2?
:::
