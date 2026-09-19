# Simple linear regression as a warm-up

With a single regressor \( x \) and an intercept, the linear model is
\[
Y_i=\beta_0+\beta_1x_i+\varepsilon_i,\qquad i=1,\dots,n,
\]{#eq-lm-simple-model}

where the \( x_i \) are known numbers and the errors satisfy the second-moment
assumptions (L1)–(L3) of [Section 5.1](01-what-a-model-claims.html):
\( \E(\varepsilon_i)=0 \), \( \Var(\varepsilon_i)=\sigma^2 \), and
\( \Cov(\varepsilon_i,\varepsilon_j)=0 \) for \( i\ne j \). The parameter \( \beta_1 \) is the
**slope**, the difference in mean response between two values of \( x \) one unit apart.
The **intercept** \( \beta_0 \) is the mean response at \( x=0 \), which may or may not be a
value of \( x \) that makes sense. Everything in this section is a special case of
[Section 5.4](04-least-squares.html). We work it out directly because the formulas are
short, because they show the structure of the general answer, and because they are the
ones used every day.

## The least squares line

For a candidate line \( b_0+b_1x \), the **residual** of case \( i \) is
\( y_i-b_0-b_1x_i \), the vertical distance from the data point to the line. The method of
least squares chooses the line that minimizes the sum of squared residuals,
\[
S(b_0,b_1)=\sum_{i=1}^n(y_i-b_0-b_1x_i)^2 .
\]
Other criteria are possible, for example the sum of absolute residuals (Chapter 45). Squared
error is chosen because it leads to linear equations, because it matches the
conditional mean of @prp-lm-error-decomposition(c), and because its statistical
properties under (L1)–(L3) can be worked out exactly, which is the business of this chapter.

We use the standard abbreviations
\[
\begin{gathered}
\bar{x}=\frac1n\sum_ix_i,\qquad \bar{y}=\frac1n\sum_iy_i,\qquad
S_{xx}=\sum_i(x_i-\bar{x})^2,\\
S_{xy}=\sum_i(x_i-\bar{x})(y_i-\bar{y}),\qquad
S_{yy}=\sum_i(y_i-\bar{y})^2 .
\end{gathered}
\]{#eq-lm-sums}

Two identities are used repeatedly: \( \sum_i(x_i-\bar{x})=0 \), and consequently
\( S_{xy}=\sum_i(x_i-\bar{x})y_i \) and \( S_{xx}=\sum_i(x_i-\bar{x})x_i \). Note that
\( S_{xx}>0 \) unless all the \( x_i \) are equal.

::: {#thm-lm-simple-ls}
[Least squares for a straight line]

Suppose the \( x_i \) are not all equal. Then \( S(b_0,b_1) \) has a unique minimizer
\( (\hat{\beta}_0,\hat{\beta}_1) \), given by
\[
\hat{\beta}_1=\frac{S_{xy}}{S_{xx}},\qquad \hat{\beta}_0=\bar{y}-\hat{\beta}_1\bar{x} .
\]{#eq-lm-simple-estimates}

It is the unique solution of the two **normal equations**
\[
\sum_{i=1}^n(y_i-b_0-b_1x_i)=0,\qquad \sum_{i=1}^n x_i(y_i-b_0-b_1x_i)=0 .
\]{#eq-lm-simple-normal}

Moreover, for all \( b_0,b_1 \),
\[
\begin{aligned}
S(b_0,b_1)={}&S(\hat{\beta}_0,\hat{\beta}_1)+n\bigl[(b_0-\hat{\beta}_0)+(b_1-\hat{\beta}_1)\bar{x}\bigr]^2\\
&+S_{xx}(b_1-\hat{\beta}_1)^2 .
\end{aligned}
\]{#eq-lm-simple-decomposition}

:::

::: {.proof}
\( S \) is a differentiable function of \( (b_0,b_1) \), with
\[
\frac{\partial S}{\partial b_0}=-2\sum_i(y_i-b_0-b_1x_i),\qquad
\frac{\partial S}{\partial b_1}=-2\sum_ix_i(y_i-b_0-b_1x_i).
\]
A minimizer must make both derivatives zero, which gives @eq-lm-simple-normal. The first
normal equation divided by \( n \) says \( \bar{y}-b_0-b_1\bar{x}=0 \), so
\( b_0=\bar{y}-b_1\bar{x} \). Substituting into the second, and using
\( y_i-\bar{y}-b_1(x_i-\bar{x}) \) for the residual,
\[
\sum_ix_i\bigl[(y_i-\bar{y})-b_1(x_i-\bar{x})\bigr]=S_{xy}-b_1S_{xx}=0,
\]
so \( b_1=S_{xy}/S_{xx} \). This is the only stationary point.

A stationary point need not be a minimum, so we check directly. Write
\( \hat{\varepsilon}_i=y_i-\hat{\beta}_0-\hat{\beta}_1x_i \), \( d_0=b_0-\hat{\beta}_0 \) and \( d_1=b_1-\hat{\beta}_1 \).
Then \( y_i-b_0-b_1x_i=\hat{\varepsilon}_i-(d_0+d_1x_i) \), and
\[
S(b_0,b_1)=\sum_i\hat{\varepsilon}_i^2-2d_0\sum_i\hat{\varepsilon}_i-2d_1\sum_ix_i\hat{\varepsilon}_i+\sum_i(d_0+d_1x_i)^2 .
\]
The two middle sums vanish by the normal equations. For the last,
\( d_0+d_1x_i=(d_0+d_1\bar{x})+d_1(x_i-\bar{x}) \), and the cross term again sums to zero, so
\( \sum_i(d_0+d_1x_i)^2=n(d_0+d_1\bar{x})^2+d_1^2S_{xx} \). This proves @eq-lm-simple-decomposition.
Both added terms are nonnegative, and since \( S_{xx}>0 \)
they are both zero only when \( d_1=0 \) and then \( d_0=0 \). So
\( (\hat{\beta}_0,\hat{\beta}_1) \) is the unique global minimizer.
:::

The decomposition @eq-lm-simple-decomposition says more than minimality. It shows the
level sets of \( S \) are ellipses centred at \( (\hat{\beta}_0,\hat{\beta}_1) \), whose shape is set by
\( n \), \( \bar{x} \) and \( S_{xx} \). When \( \bar{x} \) is far from zero the ellipses are long and
thin, tilted along the direction where \( b_0+b_1\bar{x} \) stays constant. Many lines then
fit almost equally well, provided they pass near the centre of the data. This is a first
look at the correlation between the estimated intercept and slope
(@cor-lm-simple-moments).

The **fitted values** are \( \hat{y}_i=\hat{\beta}_0+\hat{\beta}_1x_i \) and the **residuals** are
\( \hat{\varepsilon}_i=y_i-\hat{y}_i \). The minimum \( \text{SSE}=\sum_i\hat{\varepsilon}_i^2 \) is the
**residual sum of squares** (also *error sum of squares*).

::: {#prp-lm-simple-fit}
[Properties of the fitted line]

Under the conditions of @thm-lm-simple-ls:

::: {.enumerate options="label=(\alph*)"}
1. \( \sum_i\hat{\varepsilon}_i=0 \), \( \sum_ix_i\hat{\varepsilon}_i=0 \) and \( \sum_i\hat{y}_i\hat{\varepsilon}_i=0 \);

2. the line passes through \( (\bar{x},\bar{y}) \), and the mean of the fitted values is \( \bar{y} \);

3. with \( \text{SST}=S_{yy} \) and \( \text{SSR}=\sum_i(\hat{y}_i-\bar{y})^2 \),
   \[
   \text{SST}=\text{SSR}+\text{SSE},\qquad \text{SSR}=\hat{\beta}_1^2S_{xx}=\frac{S_{xy}^2}{S_{xx}},
   \]
   and so \( \text{SSE}=S_{yy}-S_{xy}^2/S_{xx} \);

4. if also \( S_{yy}>0 \), the **coefficient of determination** \( R^2=\text{SSR}/\text{SST} \) equals
   \( r^2 \), the square of the sample correlation \( r=S_{xy}/\sqrt{S_{xx}S_{yy}} \).
:::

:::

::: {.proof}
*(a)* The first two are the normal equations @eq-lm-simple-normal. The third follows
because \( \hat{y}_i \) is a combination of \( 1 \) and \( x_i \). *(b)* The first normal equation
is \( \bar{y}=\hat{\beta}_0+\hat{\beta}_1\bar{x} \), and averaging \( \hat{y}_i=y_i-\hat{\varepsilon}_i \) gives \( \bar{y} \).
*(c)* By (b), \( \hat{y}_i-\bar{y}=\hat{\beta}_1(x_i-\bar{x}) \), so \( \text{SSR}=\hat{\beta}_1^2S_{xx}=S_{xy}^2/S_{xx} \).
Also \( y_i-\bar{y}=(\hat{y}_i-\bar{y})+\hat{\varepsilon}_i \), and the cross term
\( \sum_i(\hat{y}_i-\bar{y})\hat{\varepsilon}_i \) vanishes by (a). *(d)* \( R^2=S_{xy}^2/(S_{xx}S_{yy})=r^2 \).
:::

Part (c) is the simplest instance of the analysis of variance: the total variation of
the response about its mean splits into a part along the fitted line and a part around it.
[Chapter 9](../ch09-sums-of-squares/index.html) develops the general decomposition, and
[Section 6.7](../ch06-projections/07-reparameterization.html) interprets \( R^2 \) as a
squared cosine.

## Regression toward the mean

Divide \( S_{xy}/S_{xx} \) above and below by \( n-1 \). With the sample standard deviations
\( s_x=\sqrt{S_{xx}/(n-1)} \) and \( s_y=\sqrt{S_{yy}/(n-1)} \), the slope is
\[
\hat{\beta}_1=r\,\frac{s_y}{s_x},\qquad\text{equivalently}\qquad
\frac{\hat{y}-\bar{y}}{s_y}=r\,\frac{x-\bar{x}}{s_x}.
\]{#eq-lm-toward-mean}

In standard units, the fitted value is the regressor's standardized value multiplied by
\( r \). Since \( \lvert r\rvert\le1 \) (the Cauchy–Schwarz inequality, @prp-mat-cauchy-schwarz),
a case that is two standard deviations above average in \( x \) is predicted to be *less*
than two standard deviations above average in \( y \), unless the points lie exactly on a
line. This is **regression toward the mean**, and it gave the subject its name: Galton
observed it in 1886 for the heights of parents and their adult children.

The effect is purely statistical and involves no mechanism. Regressing \( x \) on \( y \)
instead gives the slope \( r\,s_x/s_y \), which predicts \( x \) from \( y \) with the same
shrinkage toward the mean in the other direction. The two fitted lines differ unless
\( r^2=1 \). They answer different questions: the mean of \( y \) at given \( x \), and the
mean of \( x \) at given \( y \). It is a common error to read regression toward the mean as
evidence that extreme units become more ordinary over time, or that an intervention
given to the most extreme units has worked.

## A worked example

::: {#exm-lm-elnino}
[Forecasting December sea temperature]

The Niño 1+2 index is the mean sea surface temperature in a region of the Pacific off
Peru and Ecuador. It is high during El Niño events. We use the monthly values for the
\( 61 \) years \( 1950 \)–\( 2010 \) and ask how well the August value
predicts the December value four months later. Here \( x \) is the August temperature and
\( y \) the December temperature, both in degrees Celsius. The summary statistics are
\[
\begin{gathered}
\bar{x}=20.843,\qquad \bar{y}=22.693,\\
S_{xx}=77.801,\qquad S_{xy}=56.008,\qquad S_{yy}=70.380,
\end{gathered}
\]
so by @thm-lm-simple-ls
\[
\hat{\beta}_1=\frac{56.008}{77.801}=0.7199,\qquad
\hat{\beta}_0=22.693-0.7199\times20.843=7.6887 .
\]
Each extra degree in August goes with \( 0.7199 \) of a degree more in the following
December, on average. For an August temperature of \( 22 \)°C the fitted December mean is
\( 23.53 \)°C. The intercept, the fitted December temperature after an August at
\( 0 \)°C, is meaningless: no August in the data is below \( 19 \)°C.

The residual sum of squares is \( \text{SSE}=30.061 \), and
\( R^2=0.5729 \), the square of the correlation \( r=0.7569 \). The slope
is less than \( 1 \) partly because \( r<1 \) and partly because December temperatures vary a
little less than August ones: \( s_y=1.083 \) against \( s_x=1.139 \).
:::

::: {when-format="html"}
![**Figure 5.2.1.** December against August sea surface temperature in the Niño 1+2
region, 1950–2010. (a) The data and the least squares line, which passes through the
point of means (cross). (b) Residuals against the regressor. The El Niño year 1997 is
the point at the far right of both panels.](elnino_slr.svg){#fig-lm-elnino width=100%}
:::

::: {when-format="pdf"}
![December against August sea surface temperature in the Niño 1+2
region, 1950–2010. (a) The data and the least squares line, which passes through the
point of means (cross). (b) Residuals against the regressor. The El Niño year 1997 is
the point at the far right of both panels.](elnino_slr.pdf){width=100%}
:::

```{.python .run #cell-elnino-slr-fit}
import numpy as np
import statsmodels.api as sm
sst = sm.datasets.elnino.load_pandas().data
x = sst["AUG"].to_numpy()          # August temperature
y = sst["DEC"].to_numpy()          # December temperature, same year
n = len(y)

xbar, ybar = x.mean(), y.mean()
Sxx = np.sum((x - xbar) ** 2)
Sxy = np.sum((x - xbar) * (y - ybar))
Syy = np.sum((y - ybar) ** 2)
b1 = Sxy / Sxx                     # slope
b0 = ybar - b1 * xbar              # intercept
fitted = b0 + b1 * x
resid = y - fitted
print(f"n = {n}, slope = {b1:.4f}, intercept = {b0:.4f}")
```

```{.python .run #cell-elnino-slr-checks}
print("sum of residuals      ", resid.sum())
print("sum of x * residuals  ", (x * resid).sum())
SSE = resid @ resid
R2 = 1 - SSE / Syy
r = Sxy / np.sqrt(Sxx * Syy)
print(f"R^2 = {R2:.4f}, r^2 = {r**2:.4f}")
print("numpy.polyfit agrees: ", np.polyfit(x, y, 1))
```

The residual plot in [Figure 5.2.1](#fig-lm-elnino)(b) is the first diagnostic to draw
after any fit. The residuals should show no pattern against \( x \). Here there is no
obvious curvature, but the spread is not uniform. The two warmest Augusts, in 1997
and 1983, both fell during strong El Niño events. The 1997 event was still building in
December, while the 1982–83 event was over by then, and the two points sit well away
from the line on opposite sides. The warmest of
all, \( 1997 \) (August \( 24.95 \)°C,
December \( 27.08 \)°C), has a residual of \( 1.430 \)°C. Because it lies
far from \( \bar{x} \), it pulls on the slope: without it the slope would be
\( 0.6214 \). Points with regressor values far from the mean have high
*leverage*, a notion made precise in
[Section 6.8](../ch06-projections/08-leverage.html) and Chapter 20.

## Regression through the origin

Sometimes theory says the mean response is zero when \( x=0 \), and the model
\( Y_i=\beta_1x_i+\varepsilon_i \) is fitted without an intercept. Minimizing
\( \sum_i(y_i-b_1x_i)^2 \) gives one normal equation, \( \sum_ix_i(y_i-b_1x_i)=0 \), with
solution
\[
\tilde{\beta}_1=\frac{\sum_ix_iy_i}{\sum_ix_i^2},
\]
valid whenever some \( x_i\ne0 \). The residuals now satisfy \( \sum_ix_i\hat{\varepsilon}_i=0 \) but
generally *not* \( \sum_i\hat{\varepsilon}_i=0 \), and the decomposition of
@prp-lm-simple-fit(c) about \( \bar{y} \) fails. It holds instead about zero, as
\( \sum_iy_i^2=\sum_i\hat{y}_i^2+\sum_i\hat{\varepsilon}_i^2 \), which is why software reports an
"uncentred" \( R^2 \) for such fits. That number is not comparable with the \( R^2 \) of a
model with an intercept (@exr-lm-origin-r2).

Dropping the intercept is a strong assumption. In @exm-lm-elnino, forcing the line
through the origin gives slope \( 1.0877 \), a fit much worse than the line with an
intercept: the data are nowhere near \( x=0 \), and the constraint has nothing to do with how
ocean temperatures behave.

## Exercises

### A. Check your understanding

::: {#exr-lm-simple-small}
[A1]

For the five points \( (x,y)=(0,1),(1,1),(2,3),(3,4),(4,6) \), compute
\( \bar{x},\bar{y},S_{xx},S_{xy} \), the least squares line, the residuals, and
\( R^2 \). Check the three identities of @prp-lm-simple-fit(a).
:::

::: {.solution}
\( \bar{x}=2 \), \( \bar{y}=3 \), \( S_{xx}=10 \), \( S_{xy}=(-2)(-2)+(-1)(-2)+0+1\cdot1+2\cdot3=13 \).
So \( \hat{\beta}_1=1.3 \) and \( \hat{\beta}_0=3-2.6=0.4 \). The fitted values are
\( 0.4,1.7,3.0,4.3,5.6 \) and the residuals \( 0.6,-0.7,0,-0.3,0.4 \), which sum to zero. Also
\( \sum_ix_i\hat{\varepsilon}_i=0-0.7+0-0.9+1.6=0 \). \( S_{yy}=4+4+0+1+9=18 \), \( \text{SSE}=0.36+0.49+0+0.09+0.16=1.1 \),
and \( R^2=1-1.1/18=0.939 \), which equals \( 13^2/(10\cdot18) \).
:::

::: {#exr-lm-origin-r2}
[A2]

For regression through the origin, show that \( \sum_iy_i^2=\sum_i\hat{y}_i^2+\sum_i\hat{\varepsilon}_i^2 \),
and give a small data set for which \( 1-\sum_i\hat{\varepsilon}_i^2/\sum_iy_i^2 \) is close to \( 1 \)
while \( 1-\sum_i\hat{\varepsilon}_i^2/S_{yy} \) is negative.
:::

### B. Practice

::: {#exr-lm-linear-weights}
[B1]

Show that \( \hat{\beta}_1=\sum_ic_iy_i \) with \( c_i=(x_i-\bar{x})/S_{xx} \), and that
\( \sum_ic_i=0 \), \( \sum_ic_ix_i=1 \) and \( \sum_ic_i^2=1/S_{xx} \). Write \( \hat{\beta}_0 \) as
\( \sum_id_iy_i \) and find \( \sum_id_i \), \( \sum_id_ix_i \) and \( \sum_id_i^2 \).
:::

::: {.solution}
\( S_{xy}=\sum_i(x_i-\bar{x})y_i \), so \( \hat{\beta}_1=\sum_ic_iy_i \). Then
\( \sum_ic_i=0 \), \( \sum_ic_ix_i=\sum_i(x_i-\bar{x})x_i/S_{xx}=1 \), and
\( \sum_ic_i^2=S_{xx}/S_{xx}^2=1/S_{xx} \). For the intercept,
\( \hat{\beta}_0=\bar{y}-\bar{x}\hat{\beta}_1=\sum_id_iy_i \) with \( d_i=1/n-\bar{x}c_i \). Hence
\( \sum_id_i=1 \), \( \sum_id_ix_i=\bar{x}-\bar{x}=0 \), and
\( \sum_id_i^2=\sum_i(1/n^2-2\bar{x}c_i/n+\bar{x}^2c_i^2)=1/n+\bar{x}^2/S_{xx} \). These are the
coefficients that give the moments in @cor-lm-simple-moments.
:::

::: {#exr-lm-weighted-mean-slopes}
[B2]

Suppose all \( x_i \) are distinct. Show that \( \hat{\beta}_1 \) is a weighted average of the
slopes \( (y_j-y_i)/(x_j-x_i) \) of the lines through all pairs of points, with weights
proportional to \( (x_j-x_i)^2 \). *Hint:* show that
\( \sum_{i<j}(x_j-x_i)(y_j-y_i)=nS_{xy} \) and \( \sum_{i<j}(x_j-x_i)^2=nS_{xx} \).
:::

### C. Going deeper

::: {#exr-lm-lad-line}
[C1]

The **least absolute deviations** line minimizes \( \sum_i\lvert y_i-b_0-b_1x_i\rvert \).
Show that some minimizing line passes through at least two of the data points (assume the
\( x_i \) are not all equal). *Hint:* for fixed \( b_1 \) the best \( b_0 \) is a median of the
\( y_i-b_1x_i \), and the objective is piecewise linear and convex in \( b_1 \). Compare
the effect of the \( 1997 \) point in @exm-lm-elnino on this line and on the least squares
line.
:::
