# Cook's distance, DFFITS and DFBETAS

An outlier is a case the model fits badly; an **influential** case is one whose removal would change the conclusions.
A remote case that follows the trend of the others exactly has no influence, and a large residual near the centroid
moves the fit little. Influence measures, all from @thm-res-deletion, quantify different consequences of deleting a case.

## Four measures

::: {#def-res-cooks}
[Cook's distance, DFFITS, DFBETAS]

Let \( h_{ii}<1 \), and let \( c_{jj} \) be the \( j \)th diagonal entry of \( (\X\T\X)^{-1} \), so that \( s^2c_{jj} \) is the estimated
variance of \( \hat{\beta}_j \).

::: {.enumerate options="label=(\alph*)"}
1. **Cook's distance** is
   \[
   D_i=\frac{(\hbeta_{(i)}-\hbeta)\T\X\T\X(\hbeta_{(i)}-\hbeta)}{p\,s^2}.
   \]

2. **DFFITS** is the change in case \( i \)'s own fitted value in units of its standard error,
   \[
   \text{DFFITS}_i=\frac{\hat y_i-\x_{(i)}\T\hbeta_{(i)}}{s_{(i)}\sqrt{h_{ii}}}.
   \]

3. **DFBETAS** is the change in one coefficient in units of its standard error,
   \[
   \text{DFBETAS}_{ij}=\frac{\hat{\beta}_j-\hat{\beta}_{j(i)}}{s_{(i)}\sqrt{c_{jj}}}.
   \]

4. **COVRATIO** compares the generalized variances of the estimates with and without case \( i \),
   \[
   \text{COVRATIO}_i=\frac{\det\bigl(s_{(i)}^2(\X_{(i)}\T\X_{(i)})^{-1}\bigr)}{\det\bigl(s^2(\X\T\X)^{-1}\bigr)}.
   \]
:::

:::

Cook (1977) proposed \( D_i \); DFFITS, DFBETAS and COVRATIO are from Belsley, Kuh and Welsch (1980). DFFITS and DFBETAS
use \( s_{(i)} \), for the same reason as \( t_i \). None requires a refit.

::: {#thm-res-influence}
[Influence identities]

Let \( h_{ii}<1 \) and \( n-p\ge2 \), and let \( \tilde{\x}_j \) and \( \ell_{ij} \) be as in @def-res-partial-leverage.

::: {.enumerate options="label=(\alph*)"}
1. \( D_i=\dfrac{\norm{\hY-\hY_{(i)}}^2}{p\,s^2}=\dfrac{r_i^2}{p}\cdot\dfrac{h_{ii}}{1-h_{ii}} \).

2. \( \text{DFFITS}_i=t_i\sqrt{h_{ii}/(1-h_{ii})} \), and \( D_i=\text{DFFITS}_i^2\,s_{(i)}^2/(p\,s^2) \).

3. \( \text{DFBETAS}_{ij}=t_i\,\dfrac{\tilde x_{ij}}{\norm{\tilde{\x}_j}\sqrt{1-h_{ii}}} \), so
   \( \text{DFBETAS}_{ij}^2=t_i^2\,\ell_{ij}/(1-h_{ii}) \).

4. For every \( \blambda\ne\bzero \),
   \[
   \frac{\bigl(\blambda\T\hbeta-\blambda\T\hbeta_{(i)}\bigr)^2}{s_{(i)}^2\,\blambda\T(\X\T\X)^{-1}\blambda}\le\text{DFFITS}_i^2 ,
   \]
   with equality at \( \blambda=\x_{(i)} \). In particular \( \lvert\text{DFBETAS}_{ij}\rvert\le\lvert\text{DFFITS}_i\rvert \) for every \( j \).

5. \( \text{COVRATIO}_i=\dfrac{1}{1-h_{ii}}\Bigl(\dfrac{n-p}{n-p-1+t_i^2}\Bigr)^{p} \).
:::

:::

::: {.proof}
Put \( \mathbf{d}=\hbeta-\hbeta_{(i)}=(\X\T\X)^{-1}\x_{(i)}g \) with \( g=\hat{\varepsilon}_i/(1-h_{ii}) \) (@thm-res-deletion(b)).

(a) \( \hY-\hY_{(i)}=\X\mathbf{d} \), and \( \norm{\X\mathbf{d}}^2=\mathbf{d}\T\X\T\X\mathbf{d}=g^2\x_{(i)}\T(\X\T\X)^{-1}\x_{(i)}=g^2h_{ii} \).
So \( pD_is^2=\hat{\varepsilon}_i^2h_{ii}/(1-h_{ii})^2 \), and dividing by \( s^2 \) and using \( r_i^2=\hat{\varepsilon}_i^2/\{s^2(1-h_{ii})\} \) gives (a).

(b) The numerator of \( \text{DFFITS}_i \) is \( \x_{(i)}\T\mathbf{d}=h_{ii}g \), so
\( \text{DFFITS}_i=\sqrt{h_{ii}}\,\hat{\varepsilon}_i/\{s_{(i)}(1-h_{ii})\}=t_i\sqrt{h_{ii}/(1-h_{ii})} \). Then
\( \text{DFFITS}_i^2s_{(i)}^2=\hat{\varepsilon}_i^2h_{ii}/(1-h_{ii})^2=pD_is^2 \).

(c) By @thm-proj-fwl(a) the \( j \)th row of \( (\X\T\X)^{-1}\X\T \) is \( \tilde{\x}_j\T/\norm{\tilde{\x}_j}^2 \), so the \( j \)th entry of
\( (\X\T\X)^{-1}\x_{(i)} \), the \( i \)th column of that matrix, is \( \tilde x_{ij}/\norm{\tilde{\x}_j}^2 \). By @eq-proj-partitioned-inverse,
\( c_{jj}=1/\norm{\tilde{\x}_j}^2 \). Hence
\[
\text{DFBETAS}_{ij}=\frac{\tilde x_{ij}\,g\,\norm{\tilde{\x}_j}^{-2}}{s_{(i)}\norm{\tilde{\x}_j}^{-1}}
=\frac{\tilde x_{ij}\,\hat{\varepsilon}_i}{s_{(i)}\norm{\tilde{\x}_j}(1-h_{ii})},
\]
which is the stated form.

(d) \( \blambda\T\mathbf{d}=g\,\blambda\T(\X\T\X)^{-1}\x_{(i)} \). Let \( \mathbf{K}^{1/2} \) be the symmetric square root of \( \mathbf{K}=(\X\T\X)^{-1} \) (@thm-mat-square-root). The
Cauchy–Schwarz inequality (@prp-mat-cauchy-schwarz) for \( \mathbf{K}^{1/2}\blambda \) and \( \mathbf{K}^{1/2}\x_{(i)} \) gives
\( \bigl(\blambda\T(\X\T\X)^{-1}\x_{(i)}\bigr)^2\le\blambda\T(\X\T\X)^{-1}\blambda\cdot h_{ii} \), with equality when \( \blambda \) is a multiple of
\( \x_{(i)} \). So the ratio is at most \( g^2h_{ii}/s_{(i)}^2=\text{DFFITS}_i^2 \). Taking \( \blambda=\vect{e}_j \) gives \( \text{DFBETAS}_{ij}^2 \).

(e) By @exm-mat-deletion, \( \det(\X_{(i)}\T\X_{(i)})=(1-h_{ii})\det(\X\T\X) \), so
\( \text{COVRATIO}_i=(s_{(i)}^2/s^2)^p/(1-h_{ii}) \). From the proof of @thm-res-external-t(c),
\( s_{(i)}^2/s^2=(n-p-r_i^2)/(n-p-1) \), and solving \( t_i^2=r_i^2(n-p-1)/(n-p-r_i^2) \) for \( r_i^2 \) gives
\( n-p-r_i^2=(n-p)(n-p-1)/(n-p-1+t_i^2) \). Substitute.
:::

## What the measures mean

**Cook's distance** has two readings. By (a) it is the squared distance the fitted vector moves when case \( i \) is
deleted, in units of \( p\,s^2 \). By @thm-ci-ellipsoid with \( \bLambda=\I \), the \( 100(1-\alpha)\% \) confidence ellipsoid for \( \bbeta \) is
\( \{\bb:(\bb-\hbeta)\T\X\T\X(\bb-\hbeta)\le p\,s^2F_\alpha(p,n-p)\} \), so \( \hbeta_{(i)} \) lies on the boundary of the ellipsoid of
level \( \Pr\{F(p,n-p)\le D_i\} \).

By the second equality in (a), \( D_i \) is misfit, \( r_i^2 \), times leverage, \( h_{ii}/(1-h_{ii}) \): a case is influential only
if it has both. DFFITS factors the same way with \( t_i \). The
**influence plot** of [Figure 20.4.1](#fig-res-influence) places each case by leverage and studentized residual, with the
curves of constant \( D_i \), \( r^2=pD(1-h)/h \).

**DFBETAS** isolates one coefficient: by (c) the partial leverage \( \ell_{ij} \) replaces the leverage in the product, but the
factor \( 1/(1-h_{ii}) \) still involves the overall leverage. By (d), the largest standardized change in *any* linear
combination \( \blambda\T\hbeta \) is \( \lvert\text{DFFITS}_i\rvert \), so a small DFFITS rules out a large change in every coefficient
and contrast. **COVRATIO** measures precision rather than location: by (e) it exceeds one when \( h_{ii} \) is large and \( t_i \) small
(the case adds precision) and falls below one when \( t_i \) is large (the case inflates \( s^2 \)).

## The null law of Cook's distance and the choice of cut-offs

How large must an influence measure be before it matters? Since \( D_i \) is a multiple of \( r_i^2 \), its law under the
model is known exactly.

::: {#prp-res-cook-null}
[Cook's distance under the model]

Assume @eq-opt-normal-model with \( n-p\ge2 \) and all leverages below one. Then
\( D_i=\dfrac{n-p}{p}\,\dfrac{h_{ii}}{1-h_{ii}}\,B_i \) with \( B_i\sim\mathrm{Beta}\bigl(\tfrac12,\tfrac{n-p-1}2\bigr) \), so
\[
\E(D_i)=\frac{h_{ii}}{p(1-h_{ii})},\qquad \frac1n\sum_{i=1}^n\E(D_i)\ge\frac1{n-p},
\]
with equality in the last bound iff all leverages equal \( p/n \).
:::

::: {.proof}
Combine @thm-res-influence(a) with @prp-res-internal(a), and use \( \E(B_i)=1/(n-p) \). For the bound,
\( u\mapsto u/(1-u) \) is strictly convex on \( [0,1) \), so by Jensen's inequality the average of \( h_{ii}/(1-h_{ii}) \) is at least
\( \bar h/(1-\bar h)=p/(n-p) \), with equality iff the \( h_{ii} \) are equal.
:::

So under the model a typical \( D_i \) is of order \( 1/n \), which explains the popular cut-off \( 4/n \). For the state design of
@exm-res-state-influence the average of \( \E(D_i) \) is \( 0.0240 \), against the bound \( 1/(n-p)=0.0213 \). Comparing \( D_i \) with the
median of \( F(p,n-p) \), which is near one, is far stricter: it flags only cases whose deletion moves the estimate to the edge
of a \( 50\% \) confidence region. The usual cut-offs are these.

| Measure | Cut-off | Reasoning |
|---|---|---|
| \( h_{ii} \) | \( 2p/n \) | twice the average leverage |
| \( D_i \) | \( 4/n \), or the median of \( F(p,n-p) \) | four times the typical value; or the edge of the \( 50\% \) region |
| \( \lvert\text{DFFITS}_i\rvert \) | \( 2\sqrt{p/n} \) | \( \lvert t_i\rvert\approx2 \) at a case of average leverage \( p/n \) |
| \( \lvert\text{DFBETAS}_{ij}\rvert \) | \( 2/\sqrt n \) | \( \lvert t_i\rvert\approx2 \) at a case of average partial leverage \( 1/n \) |
| \( \lvert\text{COVRATIO}_i-1\rvert \) | \( 3p/n \) | first-order expansion of (e) at \( \lvert t_i\rvert\approx2 \) or \( h_{ii}\approx2p/n \) |

All are arbitrary in the same way: each puts a case with studentized residual about two and typical leverage near the
line. Calibrating \( D_i \) exactly by
@prp-res-cook-null would not help: for a fixed design \( D_i \) is a multiple of \( r_i^2 \), so the resulting test is an outlier
test with a case-dependent level, discarding the leverage factor that is the point of an influence measure. Influence is a
consequence to be assessed, not a hypothesis: the right threshold is whether deleting the case changes a conclusion.

::: {#exm-res-state-influence}
[How much does the District of Columbia matter?]

In the murder-rate regression with all \( 51 \) jurisdictions, the District of Columbia has \( D_i=3.430 \), seven times the next
largest (Mississippi, \( 0.481 \)); the median of \( F(4,47) \) is \( 0.851 \), and deleting the District moves \( \hbeta \) to the boundary
of the \( 0.985 \) confidence ellipsoid. Its DFFITS is \( 4.333 \) (cut-off \( 0.560 \)) and its DFBETAS for poverty \( -1.090 \)
(cut-off \( 0.280 \)). Its partial leverage for poverty is only \( 0.032 \) (@exm-res-state-partial), but by @thm-res-influence(c)
it is multiplied by \( t_i^2/(1-h_{ii}) \): \( 4.278^2\times0.032/(1-0.506)\approx1.19 \), so \( \lvert\text{DFBETAS}\rvert\approx1.09 \). Alaska, with the largest
partial leverage (\( 0.141 \)), has \( t_i=-0.760 \) and DFBETAS \( 0.315 \).

The consequences are substantive. With the District, the poverty coefficient is \( 0.1551 \) (standard error \( 0.1059 \),
\( t=1.47 \)); without it, \( 0.2538 \) (standard error \( 0.0934 \), \( t=2.72 \)), the fit of @exm-proj-fwl-crime. One jurisdiction
decides whether poverty has a detectable association with murder rates given the other regressors, and the single-parent
coefficient moves from \( 0.6256 \) to \( 0.3980 \). The cut-offs flag \( 2 \) cases by \( D_i>4/n=0.078 \), \( 3 \) by DFFITS, \( 5 \) by DFBETAS for
poverty and \( 4 \) by leverage above \( 2p/n=0.157 \), but only the District changes a conclusion. Its COVRATIO is \( 0.578 \): it inflates the estimated
generalized variance through its contribution to \( s^2 \).
:::

::: {when-format="html"}
![**Figure 20.4.1.** Influence plot for the state regression: leverage against \( r_i \), area proportional to Cook's
distance, curves \( D_i=0.5 \) (dotted) and \( D_i=1 \) (dashed), and the line \( h_{ii}=2p/n \).](influence_plot.svg){#fig-res-influence width=72%}
:::

::: {when-format="pdf"}
![Influence plot for the state regression: leverage against \( r_i \), area proportional to Cook's
distance, curves \( D_i=0.5 \) (dotted) and \( D_i=1 \) (dashed), and the line \( h_{ii}=2p/n \).](influence_plot.pdf){width=72%}
:::

```{.python .run #cell-state-diagnostics-influence}
import numpy as np
import statsmodels.api as sm
from scipy import stats
data = sm.datasets.statecrime.load_pandas().data           # District of Columbia kept
y = data["murder"].to_numpy()
X = np.column_stack([np.ones(len(y)), data[["poverty", "single", "urban"]]])
names = np.array(data.index)
n, p = X.shape
Q, R = np.linalg.qr(X)
beta = np.linalg.solve(R, Q.T @ y)
e = y - X @ beta                                   # raw residuals
h = np.sum(Q ** 2, axis=1)                         # leverages
s2 = e @ e / (n - p)
r = e / np.sqrt(s2 * (1 - h))                      # internally studentized
s2_del = (e @ e - e ** 2 / (1 - h)) / (n - p - 1)  # s^2 without case i
t = e / np.sqrt(s2_del * (1 - h))                  # externally studentized

cook = r ** 2 * h / (p * (1 - h))                  # Cook's distance
dffits = t * np.sqrt(h / (1 - h))
K = np.linalg.inv(X.T @ X)
C = K @ X.T                                        # row j: weights of y in beta_j
dfbetas = (C * e / (1 - h)).T / np.sqrt(np.outer(s2_del, np.diag(K)))
for i in np.argsort(-cook)[:3]:
    print(f"{names[i]:22s} D = {cook[i]:.3f}  DFFITS = {dffits[i]:6.3f}"
          f"  DFBETAS(poverty) = {dfbetas[i, 1]:6.3f}")
print("conventional cut-offs: D > 4/n =", round(4 / n, 3),
      " |DFFITS| > 2 sqrt(p/n) =", round(2 * np.sqrt(p / n), 3),
      " |DFBETAS| > 2/sqrt(n) =", round(2 / np.sqrt(n), 3))
```

## Exercises

### A. Check your understanding

::: {#exr-res-cook-hand}
[A1]

In a regression with \( p=3 \), case \( i \) has \( r_i=2 \). Compute \( D_i \) when \( h_{ii}=0.1 \) and when \( h_{ii}=0.6 \). Which case would you
look at first, and why might neither be alarming?
:::

::: {.solution}
By @thm-res-influence(a), \( D_i=(4/3)\,h_{ii}/(1-h_{ii}) \): \( (4/3)(1/9)=4/27\approx0.15 \) and \( (4/3)(1.5)=2 \). The second
moves the coefficient vector to about the edge of a high-level confidence region. Neither is alarming by itself if the
change it produces in the estimates of interest is small relative to their uncertainty, or if the case is known to be
correct and the model to hold in its region.
:::

::: {#exr-res-zero-residual}
[A2]

Show that \( D_i=0 \) iff \( \hat{\varepsilon}_i=0 \), whatever the leverage. Give an example of a case with leverage \( 0.9 \) and no influence.
:::

### B. Practice

::: {#exr-res-catcher}
[B1]

Let \( \mathbf{C}=(\X\T\X)^{-1}\X\T \). Show that the sum of squares of the \( j \)th row of \( \mathbf{C} \) is \( c_{jj} \), and deduce that
\( \sum_i\ell_{ij}=1 \) in a second way. Why does this suggest a cut-off of order \( 1/\sqrt n \) for DFBETAS?
:::

::: {.solution}
\( \mathbf{C}\mathbf{C}\T=(\X\T\X)^{-1}\X\T\X(\X\T\X)^{-1}=(\X\T\X)^{-1} \), whose \( j \)th diagonal entry is \( c_{jj} \). The \( j \)th row of \( \mathbf{C} \)
is \( \tilde{\x}_j\T/\norm{\tilde{\x}_j}^2 \) (proof of @thm-res-influence(c)), so \( \sum_i\ell_{ij}=\sum_i\tilde x_{ij}^2/\norm{\tilde{\x}_j}^2=c_{jj}\norm{\tilde{\x}_j}^2=1 \).
The partial leverages for a coefficient average \( 1/n \). With \( \lvert t_i\rvert\approx2 \) and \( h_{ii} \) small,
\( \lvert\text{DFBETAS}_{ij}\rvert\approx2\sqrt{\ell_{ij}}\approx2/\sqrt n \) for a typical case.
:::

::: {#exr-res-simple-dfbetas}
[B2]

For simple regression, show that the DFBETAS of case \( i \) for the slope is \( t_i(x_i-\bar x)/\sqrt{S_{xx}(1-h_{ii})} \). Which cases
can change the slope most?
:::


::: {#exr-res-covratio}
[B3]

Use @thm-res-influence(e) to show that for large \( n \) and small \( p/n \), \( \text{COVRATIO}_i\approx1+h_{ii}-p(t_i^2-1)/n \).
Deduce the cut-off \( \lvert\text{COVRATIO}_i-1\rvert>3p/n \) from the two extreme cases \( h_{ii}\approx0,\ \lvert t_i\rvert=2 \) and
\( h_{ii}=2p/n,\ t_i=0 \).
:::

::: {.solution}
\( 1/(1-h)\approx1+h \), and \( \bigl((n-p)/(n-p-1+t^2)\bigr)^p=\bigl(1+(t^2-1)/(n-p)\bigr)^{-p}\approx1-p(t^2-1)/n \). Multiplying,
\( \text{COVRATIO}\approx1+h-p(t^2-1)/n \) to first order. With \( h\approx0 \) and \( t^2=4 \) this is \( 1-3p/n \): a case of no
leverage and a residual of two standard errors lowers the precision by that much. With \( t=0 \) and \( h=2p/n \) it is
\( 1+3p/n \): a well-fitting case of high leverage raises it by the same amount. Deviations beyond \( 3p/n \) in either
direction are more extreme than both.
:::

### C. Going deeper

::: {#exr-res-cook-bound}
[C1]

How influential can one case be, whatever its response? Fix a design of full rank with
\( h_{ii}<1 \) for all \( i \), and let \( \y \) range over \( \Real^n \).

::: {.enumerate options="label=(\alph*)"}
1. Show that \( \hat{\varepsilon}_i=\bigl\{(\I-\M)\vect{e}_i\bigr\}\T\he \), and deduce from the Cauchy–Schwarz
   inequality (@prp-mat-cauchy-schwarz) that \( r_i^2\le n-p \) for *every* data vector, with
   equality iff \( \he \) is a multiple of \( (\I-\M)\vect{e}_i \). (@prp-res-internal(b) gives the same
   bound from the beta law; this argument needs no model at all.)

2. Deduce \( D_i\le\dfrac{n-p}{p}\cdot\dfrac{h_{ii}}{1-h_{ii}} \), and show that the bound is
   attained: take \( \y=\vect{e}_i \).

3. Conclude that in a design with equal leverages \( h_{ii}=p/n \) no case can have \( D_i>1 \),
   whatever the responses. Using the confidence-ellipsoid reading of \( D_i \), say what this
   means for how far one case can move \( \hbeta \).

4. Contrast a design with a leverage close to one, and explain how (b) makes the
   \( 2p/n \) cut-off for leverage a statement about the *worst case* rather than the observed one.
:::
:::

::: {.solution}
(a) Since \( \I-\M \) is symmetric and idempotent and \( \he=(\I-\M)\y \),
\( \hat{\varepsilon}_i=\vect{e}_i\T(\I-\M)\y=\{(\I-\M)\vect{e}_i\}\T\{(\I-\M)\y\} \). By Cauchy–Schwarz,
\( \hat{\varepsilon}_i^2\le\norm{(\I-\M)\vect{e}_i}^2\norm{\he}^2=(1-h_{ii})\,\text{SSE} \), with equality iff
the two vectors are proportional. Divide by \( s^2(1-h_{ii})=\text{SSE}(1-h_{ii})/(n-p) \) to get
\( r_i^2\le n-p \).

(b) Substitute in @thm-res-influence(a), \( D_i=r_i^2h_{ii}/\{p(1-h_{ii})\} \). For \( \y=\vect{e}_i \)
the residual is \( \he=(\I-\M)\vect{e}_i \), which is proportional to itself, so equality holds in (a)
and \( D_i \) attains the bound.

(c) With \( h_{ii}=p/n \) the bound is \( \dfrac{n-p}{p}\cdot\dfrac{p/n}{1-p/n}=1 \). By the
ellipsoid reading, deleting any one case moves \( \hbeta \) at most to the boundary of the
confidence ellipsoid of level \( \Pr\{F(p,n-p)\le1\} \), which is near \( 0.5 \) for moderate
\( n-p \): in a design that spreads leverage evenly, no single case, however wild its response,
can move the estimate beyond a middling confidence region.

(d) The bound grows without limit as \( h_{ii}\to1 \), and the same case would then also have a
tiny residual, so the *observed* \( D_i \) need not be large. Leverage is the design's exposure to
a bad response: a case with \( h_{ii} \) near one can do unbounded damage if its \( y_i \) is wrong,
and nothing in the data reveals whether it is. That is why leverage is screened on its own,
before any residual is looked at.
:::

::: {#exr-res-cook-dffits-order}
[C2]

Cook's distance and DFFITS differ only in using \( s \) or \( s_{(i)} \), but they can rank cases
differently. Take \( n-p=6 \) and two cases with
\[
h_{11}=\tfrac{10}{19},\quad r_1^2=1,
\qquad
h_{22}=\tfrac1{10},\quad r_2^2=5 .
\]

::: {.enumerate options="label=(\alph*)"}
1. Using @thm-res-external-t(c), show that \( t_1^2=1 \) and \( t_2^2=25 \).

2. Show that \( pD_1=10/9 \) and \( pD_2=5/9 \), while \( \text{DFFITS}_1^2=10/9 \) and
   \( \text{DFFITS}_2^2=25/9 \). So \( D \) ranks case 1 first and DFFITS ranks case 2 first.

3. Check that such a data set is possible: find \( p \) and \( n \) for which the constraints
   \( \sum_ih_{ii}=p \) and \( \sum_i(1-h_{ii})r_i^2=n-p \) (@exr-res-weighted-sum) can be met by
   the remaining cases.

4. Explain the disagreement in one sentence, and say which measure you would use to decide
   whether case 2 is an outlier and which to decide whether it matters.
:::
:::

::: {.solution}
(a) \( t_i^2=r_i^2(n-p-1)/(n-p-r_i^2)=5r_i^2/(6-r_i^2) \), so \( t_1^2=5/5=1 \) and
\( t_2^2=25/1=25 \). The second case has a residual near the algebraic maximum
\( r^2\le n-p=6 \) of @exr-res-cook-bound, so deleting it shrinks \( s^2 \) sharply and \( t_2^2 \)
is much larger than \( r_2^2 \).

(b) With \( g_i=h_{ii}/(1-h_{ii}) \) we have \( g_1=(10/19)/(9/19)=10/9 \) and \( g_2=1/9 \). Then
\( pD_i=r_i^2g_i \) by @thm-res-influence(a), giving \( 10/9 \) and \( 5/9 \); and
\( \text{DFFITS}_i^2=t_i^2g_i \) by @thm-res-influence(b), giving \( 10/9 \) and \( 25/9 \). So
\( D_1=2D_2 \) while \( \text{DFFITS}_2^2=2.5\,\text{DFFITS}_1^2 \).

(c) Take \( p=3 \) and \( n=9 \). The two cases use \( 10/19+1/10=119/190 \) of the leverage budget
\( \sum_ih_{ii}=3 \), leaving \( 451/190 \) for the other seven, each of which can then have
leverage well below one. They use \( 1\cdot\tfrac9{19}+5\cdot\tfrac9{10}=189/38 \) of the budget
\( \sum_i(1-h_{ii})r_i^2=n-p=6 \), leaving \( 39/38 \) for the other seven. Both constraints can
be met.

(d) \( D_i \) measures the move against \( s \), which the case's own residual inflates, while
DFFITS measures it against \( s_{(i)} \), which it does not; so a case with a large residual
looks worse to DFFITS, and a case with high leverage and a modest residual looks worse to
\( D_i \). To ask
whether case 2 is an outlier, use \( t_i \), which has an exact null
law (@thm-res-external-t(b)); to ask whether it matters, use \( D_i \), whose ellipsoid reading is a
statement about the estimate rather than about the error.
:::
