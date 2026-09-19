# Testing lack of fit with replicates

Every test and interval so far has assumed that the mean function is right, that
\( \E(\Y)\in\C(\X) \). When it is wrong, the residual mean square overestimates the error variance:
@exm-rv-rss-bias showed that \( \E(\text{SSE})=\sigma^2(n-r)+\norm{(\I-\M)\boldsymbol{\theta}}^2 \) for a mean vector \( \boldsymbol{\theta} \)
outside the model space. Separating noise from systematic error needs an estimate of \( \sigma^2 \) free of the form of the mean,
and repeated rows provide one: observations at identical settings differ only by noise. The idea, and
the test built on it, go back to Fisher (1922).

## Replicates and the cell-means model

Suppose the model matrix has \( c \) distinct rows \( \x_1\T,\dots,\x_c\T \), the \( i \)th occurring \( n_i\ge1 \) times, with
\( n=\sum_in_i \). Index the observations by the distinct row and the replicate, \( y_{ij} \) with
\( j=1,\dots,n_i \), so that the model is
\[
\E(y_{ij})=\x_i\T\bbeta,\qquad i=1,\dots,c,\quad j=1,\dots,n_i .
\]
Let \( \Z \) be the \( n\times c \) matrix of group indicators, the model matrix of the one-way layout in which each
distinct row has its own mean \( \mu_i \). Call two matrices with \( n \) rows *alike in row structure* as \( \X \) if
their rows are equal whenever the corresponding rows of \( \X \) are equal.

::: {#lem-cor-row-structure}
[The largest model with the same replicates]

::: {.enumerate options="label=(\alph*)"}
1. \( \C(\X)\subseteq\C(\Z) \). More generally, every matrix \( \W \) with \( n \) rows whose rows are constant within the
   groups satisfies \( \C(\W)\subseteq\C(\Z) \).

2. \( \M_Z\y \) replaces each observation by its group mean \( \bar y_i \), and \( \M\y \) is constant within groups.
:::

:::

::: {.proof}
(a) Let \( \bar{\W} \) be the \( c\times s \) matrix whose \( i \)th row is the common row of \( \W \) in group \( i \). Then
\( \W=\Z\bar{\W} \), so \( \C(\W)\subseteq\C(\Z) \). The matrix \( \X \) is a special case.

(b) The first claim is @exm-proj-oneway-M. For the second, \( \M\y\in\C(\X)\subseteq\C(\Z) \), and the vectors of \( \C(\Z) \) are constant within
groups.
:::

So the cell-means model is the largest linear model that treats identical rows identically. Any
mean function of the regressors at all, linear or not, produces a mean vector in \( \C(\Z) \).

::: {#def-cor-pure-error}
[Pure error and lack of fit]

With \( \M_Z \) the projection onto \( \C(\Z) \), the **pure error** and **lack-of-fit** sums of squares are
\[
\begin{aligned}
\text{SSPE}&=\norm{(\I-\M_Z)\y}^2=\sum_{i=1}^c\sum_{j=1}^{n_i}(y_{ij}-\bar y_i)^2,\\
\text{SSLF}&=\norm{(\M_Z-\M)\y}^2=\sum_{i=1}^cn_i(\bar y_i-\hat y_i)^2,
\end{aligned}
\]{#eq-cor-sspe}

with \( n-c \) and \( c-r \) degrees of freedom, where \( \hat y_i=\x_i\T\hbeta \) is the common fitted value of group \( i \)
and \( r=\rank(\X) \).
:::

The two expressions for each sum of squares agree by @lem-cor-row-structure(b): \( (\I-\M_Z)\y \) has entries
\( y_{ij}-\bar y_i \), and \( (\M_Z-\M)\y \) has entries \( \bar y_i-\hat y_i \), repeated \( n_i \) times.

## The test

::: {#thm-cor-lack-of-fit}
[Pure-error lack-of-fit test]

Let \( \Y\sim\Normal_n(\boldsymbol{\theta},\sigma^2\I) \) with \( \boldsymbol{\theta}\in\C(\Z) \), that is, \( \E(y_{ij})=\mu_i \) depends only on the row of
\( \X \); assume \( r<c<n \). Then:

::: {.enumerate options="label=(\alph*)"}
1. \( \text{SSE}=\text{SSLF}+\text{SSPE} \), and \( \text{SSLF} \) and \( \text{SSPE} \) are independent;

2. \( \text{SSPE}/\sigma^2\sim\chi^2(n-c) \) for every \( \mu_1,\dots,\mu_c \);

3. \( \text{SSLF}/\sigma^2\sim\chi^2(c-r,\gamma) \) with
   \[
\gamma=\frac{\norm{(\I-\M)\boldsymbol{\theta}}^2}{\sigma^2}=\frac1{\sigma^2}\sum_{i=1}^cn_i\bigl(\mu_i-\tilde\mu_i\bigr)^2,
\]{#eq-cor-lof-noncentrality}

   where \( \tilde\mu_i \) is the \( i \)th group value of \( \M\boldsymbol{\theta} \), the least squares fit to the true means;

4. the ratio
   \[
F=\frac{\text{SSLF}/(c-r)}{\text{SSPE}/(n-c)}\sim F(c-r,n-c,\gamma),
\]
   and \( \gamma=0 \) iff \( \boldsymbol{\theta}\in\C(\X) \), that is, iff the model is correct. The test that rejects the model
   when \( F>F_\alpha(c-r,n-c) \) has size \( \alpha \), and its power increases strictly with \( \gamma \).
:::

:::

::: {.proof}
By @lem-cor-row-structure, \( \C(\X)\subseteq\C(\Z) \), so by @thm-proj-nested \( \M_Z-\M \) is the projection onto the orthogonal
complement of \( \C(\X) \) in \( \C(\Z) \), of rank \( c-r \), and \( \I-\M=(\M_Z-\M)+(\I-\M_Z) \) splits into two orthogonal
projections. This gives the identity in (a). By @thm-qf-orthogonal-projections, \( \text{SSLF}/\sigma^2 \) and
\( \text{SSPE}/\sigma^2 \) are independent noncentral chi-squared variables with \( c-r \) and \( n-c \) degrees of freedom and
noncentralities \( \norm{(\M_Z-\M)\boldsymbol{\theta}}^2/\sigma^2 \) and \( \norm{(\I-\M_Z)\boldsymbol{\theta}}^2/\sigma^2 \). The second is zero because
\( \boldsymbol{\theta}\in\C(\Z) \), which proves (b). In the first, \( \M_Z\boldsymbol{\theta}=\boldsymbol{\theta} \), so \( (\M_Z-\M)\boldsymbol{\theta}=(\I-\M)\boldsymbol{\theta} \), whose
entries are \( \mu_i-\tilde\mu_i \), repeated \( n_i \) times; this is (c). Part (d) follows from @def-qf-noncentral-f,
the fact that \( (\I-\M)\boldsymbol{\theta}=\bzero \) iff \( \boldsymbol{\theta}\in\C(\X) \), and @thm-qf-f-power.
:::

The test is the \( F \) test of the reduced model \( \C(\X) \) inside the full model \( \C(\Z) \) (@thm-glh-f-test, @thm-qf-nested-f). What makes it special is the full model, which is chosen not by the analyst but by
the replication in the data, and which contains every mean function of the regressors. The
residual sum of squares of the working model is split into two rows of an analysis of variance:

| Source | Degrees of freedom | Sum of squares | Mean square |
|---|---|---|---|
| Lack of fit | \( c-r \) | \( \sum_in_i(\bar y_i-\hat y_i)^2 \) | \( \text{SSLF}/(c-r) \) |
| Pure error | \( n-c \) | \( \sum_i\sum_j(y_{ij}-\bar y_i)^2 \) | \( \text{SSPE}/(n-c) \) |
| Residual | \( n-r \) | \( \text{SSE} \) | \( s^2 \) |

The noncentrality @eq-cor-lof-noncentrality is the weighted squared distance between the true group
means and their best approximation by the model, in units of \( \sigma^2 \). It grows with the replication,
so replicates buy power twice: through \( \gamma \) and through the pure-error degrees of freedom.

::: {#exm-cor-stackloss-lof}
[Stack loss and air flow]

Brownlee's stack loss data record \( n=21 \) days of operation of a plant that oxidizes ammonia. Air flow,
the rate of operation, takes only \( c=7 \) distinct values, so a straight-line regression of stack loss
on air flow can be tested for lack of fit with \( 14 \) degrees of freedom for pure error. The line
has intercept \( -44.13 \) and slope \( 1.020 \), and \( R^2=0.846 \). Its residual sum of squares
\( 319.12 \) splits into \( \text{SSPE}=96.53 \) and \( \text{SSLF}=222.58 \) on \( 5 \) degrees of freedom.
The mean squares are \( 44.52 \) and \( 6.90 \), so \( F=6.46 \) with \( p \)-value
\( 0.0026 \); the \( 5\% \) point of \( F(5,14) \) is \( 2.96 \). The residual mean square of the line,
\( 16.80 \), is more than twice the pure-error estimate of \( \sigma^2 \).

[Figure 14.6.1](#fig-cor-lack-of-fit) shows why. The group means do not lie on any smooth curve: the single
day at air flow \( 70 \) has stack loss \( 15.0 \), below the mean \( 20.6 \) of the five days at air
flow \( 62 \). A parabola does not help, with \( F=7.22 \) on \( 4 \) and \( 14 \) degrees of freedom
(\( p \)-value \( 0.0023 \)). The mean of stack loss is not a function of air flow alone, which points to the
other operating variables, taken up in @exm-cor-stackloss-full.
:::

::: {when-format="html"}
![**Figure 14.6.1.** Stack loss against air flow on \( 21 \) days (replicates spread sideways), the seven
group means (bars), the least squares line and parabola.](lack_of_fit.svg){#fig-cor-lack-of-fit width=68%}
:::

::: {when-format="pdf"}
![Stack loss against air flow on \( 21 \) days (replicates spread sideways), the seven
group means (bars), the least squares line and parabola.](lack_of_fit.pdf){width=68%}
:::

```{.python .run #cell-lack-of-fit-lof}
import numpy as np
import statsmodels.api as sm
from scipy import stats
data = sm.datasets.stackloss.load_pandas().data
y = data["STACKLOSS"].to_numpy()
x = data["AIRFLOW"].to_numpy()
n = len(y)
levels, groups = np.unique(x, return_inverse=True)   # c distinct rows of X
c = len(levels)
Z = np.eye(c)[groups]                                 # cell-means model matrix

def sse(A):
    coef, *_ = np.linalg.lstsq(A, y, rcond=None)
    return np.sum((y - A @ coef) ** 2)

ss_pe = sse(Z)                                        # pure error, n - c df
for name, X in [("line", np.column_stack([np.ones(n), x])),
                ("parabola", np.column_stack([np.ones(n), x, x ** 2]))]:
    r = np.linalg.matrix_rank(X)
    ss_lf = sse(X) - ss_pe                            # lack of fit, c - r df
    F = (ss_lf / (c - r)) / (ss_pe / (n - c))
    print(f"{name:8s}: SSE {sse(X):7.2f} = SSLF {ss_lf:7.2f} ({c - r} df)"
          f" + SSPE {ss_pe:6.2f} ({n - c} df);  F = {F:.2f}, p = {stats.f.sf(F, c - r, n - c):.4f}")
```

## When pure error is not pure

The theorem assumes \( \boldsymbol{\theta}\in\C(\Z) \): the mean depends on the observation only through its row of
\( \X \). If a variable that is not in the model affects the response and varies within the groups of
replicates, as water temperature does among the days with equal air flow in the stack loss data,
then observations in a group do not share a mean. By @thm-ss-expected-mean-squares(b), applied to the
cell-means model,
\[
\E\Bigl(\frac{\text{SSPE}}{n-c}\Bigr)=\sigma^2+\frac{\norm{(\I-\M_Z)\boldsymbol{\theta}}^2}{n-c},
\]
so the "pure" error estimate is inflated, and \( F \) becomes a ratio of two noncentral chi-squared variables,
a *doubly noncentral* \( F \), which rejects less often than it would with a pure estimate of \( \sigma^2 \). Two further cautions: the pooled estimate
assumes a common variance across groups (@exr-cor-lof-hetero), and a non-significant result with few
degrees of freedom is weak evidence that the model is adequate.

::: {.remark}
[Designing for a lack-of-fit test]

The test needs \( c>r \) distinct rows and \( n>c \) observations. The design that estimates a straight line most
precisely puts half the observations at each end of the range; it has \( c=r=2 \) and cannot detect any
curvature. A few intermediate settings buy a check of the model at a small cost in precision. Replicates also cap what any model of the regressors can achieve: \( R^2\le1-\text{SSPE}/\text{SST} \) (@exr-cor-r2-ceiling).
:::

## Exercises

### A. Check your understanding

::: {#exr-cor-lof-table}
[A1]

A straight line is fitted to \( n=24 \) observations at \( c=6 \) distinct values of \( x \), each replicated four
times. The residual sum of squares is \( 50.0 \) and the pure-error sum of squares is \( 36.0 \). Complete the
analysis of variance table for lack of fit and compute \( F \). Is there evidence of lack of fit at level
\( 0.05 \), given that \( F_{0.05}(4,18)=2.93 \)?
:::

### B. Practice

::: {#exr-cor-lof-weighted}
[B1]

Show that the least squares estimate \( \hbeta \) from all \( n \) observations is also the weighted least squares
estimate from the \( c \) group means, with weights \( n_i \), and that \( \text{SSLF} \) is the weighted residual sum of
squares of that fit.
:::

::: {.solution}
By @thm-proj-nested, \( \M\y=\M\M_Z\y \): projecting \( \y \) onto \( \C(\X) \) is the same as projecting the vector of group
means \( \M_Z\y \). For a vector constant within groups, \( \norm{\M_Z\y-\X\bb}^2=\sum_in_i(\bar y_i-\x_i\T\bb)^2 \), so
minimizing it over \( \bb \) is weighted least squares on the \( c \) means with weights \( n_i \). Its minimum is
\( \norm{(\M_Z-\M)\y}^2=\text{SSLF} \).
:::

::: {#exr-cor-r2-ceiling}
[B2]

Suppose \( \bone\in\C(\X) \). Show that every model whose model matrix has the row structure of \( \X \) has
\( R^2\le1-\text{SSPE}/\text{SST} \), with equality for the cell-means model.
:::

::: {.solution}
By @lem-cor-row-structure(a) the column space \( \mathcal V \) of such a model lies in \( \C(\Z) \), so its residual sum of
squares \( \norm{(\I-\bP_{\mathcal V})\y}^2 \) is at least \( \norm{(\I-\M_Z)\y}^2=\text{SSPE} \) (projection onto a smaller space leaves a
longer residual). Hence \( R^2=1-\text{SSE}/\text{SST}\le1-\text{SSPE}/\text{SST} \), with equality when \( \mathcal V=\C(\Z) \).
:::

::: {#exr-cor-hat-rows}
[B3]

Show that rows \( i \) and \( j \) of \( \M \) are equal iff rows \( i \) and \( j \) of \( \X \) are equal. So replicates can be read
off the hat matrix.
:::

### C. Going deeper

::: {#exr-cor-lof-power}
[C1]

A straight line is to be tested for lack of fit at the four values \( x=0,1,2,3 \), each replicated \( m \)
times, when the true mean is \( \beta_0+\beta_1x+\beta_2x^2 \) with \( \beta_2/\sigma=0.5 \). Show that
\( \gamma=m\beta_2^2\norm{(\I-\M_4)\mathbf{q}}^2/\sigma^2 \), where \( \mathbf{q}=(0,1,4,9)\T \) and \( \M_4 \) is the projection for a line
on the four points, and find the smallest \( m \) for which the test at level \( 0.05 \) has power at least
\( 0.8 \).
:::

::: {#exr-cor-lof-hetero}
[C2]

Suppose the groups have different variances \( \sigma_i^2 \). Compute \( \E(\text{SSPE}) \) and \( \E(\text{SSLF}) \) under a correct
mean model, and show by an example with two large groups of very different sizes and variances that
the \( F \) test can then reject far more or far less often than \( \alpha \).
:::
