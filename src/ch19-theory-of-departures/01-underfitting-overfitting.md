# Omitting and adding regressors

The assumption \( \E(\Y)=\X\bbeta \) claims that the columns of \( \X \) are the right ones. This section compares
two ways of getting them wrong. *Underfitting* leaves out regressors that matter. *Overfitting* includes
regressors that do not. This is a first look at the bias-variance trade-off of Part VI.

## The setting

Throughout the section the truth is
\[
\E(\Y)=\X_1\bbeta_1+\X_2\bbeta_2,\qquad \Cov(\Y)=\sigma^2\I_n,
\]
where \( \X=[\X_1,\X_2] \) has full column rank \( p=p_1+p_2<n \). The *long* regression fits both
blocks and gives \( \hbeta=(\hbeta_1\T,\hbeta_2\T)\T \) and \( s^2=\norm{(\I-\M)\Y}^2/(n-p) \). The *short*
regression fits \( \X_1 \) alone and gives
\[
\hbeta_1^{\mathrm{S}}=(\X_1\T\X_1)^{-1}\X_1\T\Y,\qquad s_1^2=\frac{\norm{(\I-\M_1)\Y}^2}{n-p_1},
\]
where the superscript S marks the short fit and \( \M_1 \) projects onto \( \C(\X_1) \). Two quantities recur. The first is
\( \boldsymbol{\Pi}=(\X_1\T\X_1)^{-1}\X_1\T\X_2 \), the coefficients from regressing each column of \( \X_2 \) on \( \X_1 \),
as in @prp-lm-omitted. The second is
\[
\gamma=\frac{\norm{(\I-\M_1)\X_2\bbeta_2}^2}{\sigma^2}
=\frac{\bbeta_2\T\X_2\T(\I-\M_1)\X_2\bbeta_2}{\sigma^2} .
\]{#eq-dep-gamma}

Since \( (\M-\M_1)\X\bbeta=(\I-\M_1)\X\bbeta=(\I-\M_1)\X_2\bbeta_2 \), this is the noncentrality parameter
@eq-glh-noncentrality of the \( F \) test of \( \bbeta_2=\bzero \) in the long model. It measures how far the
true mean lies from the short model's column space, in units of \( \sigma^2 \). It is zero iff
\( \bbeta_2=\bzero \), because \( \X_2\bbeta_2\in\C(\X_1) \) forces \( \bbeta_2=\bzero \) when \( \X \) has full rank.

Both regressions are linked by an exact identity and by one covariance.

::: {#lem-dep-short-long}
[Short and long regressions]

In this setting:

::: {.enumerate options="label=(\alph*)"}
1. \( \hbeta_2=\bigl(\X_2\T(\I-\M_1)\X_2\bigr)^{-1}\X_2\T(\I-\M_1)\Y \) and
   \( \hbeta_1=\hbeta_1^{\mathrm{S}}-\boldsymbol{\Pi}\hbeta_2 \);

2. \( \Cov(\hbeta_1^{\mathrm{S}},\hbeta_2)=\mathbf{0} \) and \( \Cov(\hbeta_2)=\sigma^2\mathbf{C} \), where
   \( \mathbf{C}=\bigl(\X_2\T(\I-\M_1)\X_2\bigr)^{-1} \).
:::

:::

::: {.proof}
(a) The first formula is @thm-proj-fwl(a) with the roles of the two blocks exchanged, and the second is
@prp-lm-omitted(a). (b) Both estimators are linear in \( \Y \). By @thm-rv-linear,
\[
\Cov(\hbeta_1^{\mathrm{S}},\hbeta_2)=\sigma^2(\X_1\T\X_1)^{-1}\X_1\T(\I-\M_1)\X_2\mathbf{C}=\mathbf{0},
\]
because \( \X_1\T(\I-\M_1)=\mathbf{0} \). Similarly
\( \Cov(\hbeta_2)=\sigma^2\mathbf{C}\X_2\T(\I-\M_1)^2\X_2\mathbf{C}=\sigma^2\mathbf{C} \).
:::

## Underfitting

::: {#thm-dep-omitted}
[Omitted regressors]

Fit the short model when the long model is true.

::: {.enumerate options="label=(\alph*)"}
1. \( \E(\hbeta_1^{\mathrm{S}})=\bbeta_1+\boldsymbol{\Pi}\bbeta_2 \) and \( \Cov(\hbeta_1^{\mathrm{S}})=\sigma^2(\X_1\T\X_1)^{-1} \). The bias
   \( \boldsymbol{\Pi}\bbeta_2 \) vanishes iff \( \X_1\T\X_2\bbeta_2=\bzero \), and it vanishes for every \( \bbeta_2 \) iff
   \( \X_1\T\X_2=\mathbf{0} \).

2. The fitted values have mean \( \X\bbeta-(\I-\M_1)\X_2\bbeta_2 \). The residuals have mean
   \( (\I-\M_1)\X_2\bbeta_2 \) and covariance \( \sigma^2(\I-\M_1) \), as if the short model were right.

3. \( \E(s_1^2)=\sigma^2\bigl(1+\gamma/(n-p_1)\bigr) \), so \( s_1^2 \) overestimates \( \sigma^2 \) unless
   \( \bbeta_2=\bzero \).

4. At a new point with regressors \( (\x_{10},\x_{20}) \), the short prediction \( \x_{10}\T\hbeta_1^{\mathrm{S}} \) has
   bias \( (\boldsymbol{\Pi}\T\x_{10}-\x_{20})\T\bbeta_2 \) for the mean \( \x_{10}\T\bbeta_1+\x_{20}\T\bbeta_2 \).
:::

:::

::: {.proof}
(a) The mean is @prp-lm-omitted(b). The covariance is @thm-lm-moments(b) for the model matrix \( \X_1 \),
which needs only \( \Cov(\Y)=\sigma^2\I \). Since \( (\X_1\T\X_1)^{-1} \) is nonsingular,
\( \boldsymbol{\Pi}\bbeta_2=\bzero \) iff \( \X_1\T\X_2\bbeta_2=\bzero \), and this holds for all \( \bbeta_2 \) iff \( \X_1\T\X_2=\mathbf{0} \).
(b) \( \E(\M_1\Y)=\X_1\bbeta_1+\M_1\X_2\bbeta_2=\X\bbeta-(\I-\M_1)\X_2\bbeta_2 \), and
\( \E\bigl((\I-\M_1)\Y\bigr)=(\I-\M_1)\X_2\bbeta_2 \) because \( (\I-\M_1)\X_1=\mathbf{0} \). The covariance is
\( (\I-\M_1)\sigma^2\I(\I-\M_1)=\sigma^2(\I-\M_1) \).
(c) By @eq-rv-rss-mean with \( \M_1 \) of rank \( p_1 \) and mean \( \X\bbeta \),
\( \E\norm{(\I-\M_1)\Y}^2=\sigma^2(n-p_1)+\norm{(\I-\M_1)\X_2\bbeta_2}^2=\sigma^2(n-p_1+\gamma) \).
(d) Take expectations in \( \x_{10}\T\hbeta_1^{\mathrm{S}} \) using (a), and subtract the true mean.
:::

Nothing in the second moments signals the omission: the covariances of \( \hbeta_1^{\mathrm{S}} \) and of the residuals
are exactly what the short model claims. The damage is in the first moments. The estimator is centred at the wrong
place, and the residuals carry the systematic component \( (\I-\M_1)\X_2\bbeta_2 \) that added-variable plots (@thm-proj-fwl) and lack-of-fit tests (@thm-cor-lack-of-fit) look for. Because \( s_1^2 \) is inflated, the
reported variances are on average too large, yet as \( n \) grows they shrink while the bias does not, so
the coverage of intervals for \( \bbeta_1 \) tends to zero. Nothing requires a coefficient to be omitted: for any
true mean \( \X_1\bbeta_1+\boldsymbol{\updelta} \), such as a curved mean fitted by a line, (a)–(c) hold with
\( \X_2\bbeta_2 \) replaced by \( \boldsymbol{\updelta} \) (@prp-lm-misspecified(a)).

## Overfitting

::: {#thm-dep-overfit}
[Irrelevant regressors]

Suppose \( \bbeta_2=\bzero \), so the short model is correct, but the long model is fitted.

::: {.enumerate options="label=(\alph*)"}
1. \( \E(\hbeta_1)=\bbeta_1 \), \( \E(\hbeta_2)=\bzero \) and \( \E(s^2)=\sigma^2 \).

2. \( \Cov(\hbeta_1)=\Cov(\hbeta_1^{\mathrm{S}})+\sigma^2\boldsymbol{\Pi}\mathbf{C}\boldsymbol{\Pi}\T \). Hence
   \( \Var(\mathbf{a}\T\hbeta_1)\ge\Var(\mathbf{a}\T\hbeta_1^{\mathrm{S}}) \) for every \( \mathbf{a}\in\Real^{p_1} \), with equality iff
   \( \boldsymbol{\Pi}\T\mathbf{a}=\bzero \), and \( \Cov(\hbeta_1)=\Cov(\hbeta_1^{\mathrm{S}}) \) iff \( \X_1\T\X_2=\mathbf{0} \).

3. At a new point \( \x_0=(\x_{10}\T,\x_{20}\T)\T \),
   \[
   \begin{aligned}
\Var(\x_0\T\hbeta)={}&\Var(\x_{10}\T\hbeta_1^{\mathrm{S}})\\
&+\sigma^2(\x_{20}-\boldsymbol{\Pi}\T\x_{10})\T\mathbf{C}(\x_{20}-\boldsymbol{\Pi}\T\x_{10}),
\end{aligned}
   \]
   and summed over the design points \( \E\norm{\X\hbeta-\X\bbeta}^2=\sigma^2p \), against \( \sigma^2p_1 \) for the short fit.
:::

:::

::: {.proof}
(a) The long model is correct (with \( \bbeta_2=\bzero \)), so @thm-lm-moments(a) and @thm-lm-sigma2 apply.
(b) By @lem-dep-short-long(a), \( \hbeta_1=\hbeta_1^{\mathrm{S}}-\boldsymbol{\Pi}\hbeta_2 \) with the two terms uncorrelated, so
\( \Cov(\hbeta_1)=\Cov(\hbeta_1^{\mathrm{S}})+\boldsymbol{\Pi}\Cov(\hbeta_2)\boldsymbol{\Pi}\T \), and @lem-dep-short-long(b) gives
\( \Cov(\hbeta_2) \). Because \( \mathbf{C} \) is positive definite, \( \mathbf{a}\T\boldsymbol{\Pi}\mathbf{C}\boldsymbol{\Pi}\T\mathbf{a}=0 \) iff
\( \boldsymbol{\Pi}\T\mathbf{a}=\bzero \), and \( \boldsymbol{\Pi}\mathbf{C}\boldsymbol{\Pi}\T=\mathbf{0} \) iff \( \boldsymbol{\Pi}=\mathbf{0} \), iff \( \X_1\T\X_2=\mathbf{0} \).
(c) Write \( \x_0\T\hbeta=\x_{10}\T\hbeta_1^{\mathrm{S}}+(\x_{20}-\boldsymbol{\Pi}\T\x_{10})\T\hbeta_2 \) using @lem-dep-short-long(a), and
use the zero covariance of (b) of that lemma. For the sum, \( \X\hbeta=\M\Y \) is unbiased for \( \X\bbeta \) with
covariance \( \sigma^2\M \), so \( \E\norm{\M\Y-\X\bbeta}^2=\sigma^2\tr\M=\sigma^2p \). The short fit
\( \M_1\Y \) is also unbiased, and its trace is \( p_1 \).
:::

An irrelevant regressor costs no bias, only variance, and the price grows with its correlation with the relevant
columns: for one extra column the variance of a coefficient is multiplied by \( 1/(1-r^2) \), with \( r \) the partial
correlation of the two regressors (@eq-proj-vif-preview). 

## Bias against variance

When \( \bbeta_2 \) is small but not zero, the short regression is biased but less variable. Which is better?
Compare the two estimators of \( \bbeta \) by their mean squared error matrices, writing
\( \hbeta^{\mathrm{S}}=\bigl((\hbeta_1^{\mathrm{S}})\T,\bzero\T\bigr)\T \) for the short estimator, which sets \( \bbeta_2 \) to zero.

::: {#lem-dep-rank-one}
[A rank-one comparison]

Let \( \mathbf{C} \) be a positive definite \( q\times q \) matrix, \( \bb\in\Real^q \) and \( c>0 \). Then
\( c\,\mathbf{C}-\bb\bb\T \) is nonnegative definite iff \( \bb\T\mathbf{C}^{-1}\bb\le c \).
:::

::: {.proof}
By @lem-ci-cauchy-schwarz with \( \W=\mathbf{C} \) and \( \bu=\bb \), the inequality \( \bb\T\mathbf{C}^{-1}\bb\le c \) holds iff
\( (\mathbf{a}\T\bb)^2\le c\,\mathbf{a}\T\mathbf{C}\mathbf{a} \) for every \( \mathbf{a} \), that is, iff
\( \mathbf{a}\T(c\,\mathbf{C}-\bb\bb\T)\mathbf{a}\ge0 \) for every \( \mathbf{a} \).
:::

::: {#thm-dep-mse}
[When the short regression wins]

In the setting of this section, let \( \mathbf{G}=\begin{pmatrix}-\boldsymbol{\Pi}\\ \I_{p_2}\end{pmatrix} \), a
\( p\times p_2 \) matrix.

::: {.enumerate options="label=(\alph*)"}
1. The mean squared error matrix of \( \hbeta^{\mathrm{S}} \) is
   \[
   \E\bigl[(\hbeta^{\mathrm{S}}-\bbeta)(\hbeta^{\mathrm{S}}-\bbeta)\T\bigr]
   =\Cov(\hbeta)-\mathbf{G}\bigl(\sigma^2\mathbf{C}-\bbeta_2\bbeta_2\T\bigr)\mathbf{G}\T .
   \]

2. \( \E(\blambda\T\hbeta^{\mathrm{S}}-\blambda\T\bbeta)^2\le\Var(\blambda\T\hbeta) \) for every \( \blambda\in\Real^p \) iff
   \( \gamma\le1 \). If \( \gamma>1 \), the reverse strict inequality holds for some \( \blambda \).

3. Over the design points, \( \E\norm{\X\hbeta^{\mathrm{S}}-\X\bbeta}^2=\sigma^2(p_1+\gamma) \), against
   \( \sigma^2p \) for the long fit. The short fit is better there iff \( \gamma<p_2 \).
:::

:::

::: {.proof}
(a) By @lem-dep-short-long(a), \( \hbeta=\hbeta^{\mathrm{S}}+\mathbf{G}\hbeta_2 \), and \( \hbeta^{\mathrm{S}} \) is uncorrelated
with \( \hbeta_2 \) by part (b) of the lemma. Since \( \E(\hbeta)=\bbeta \) and \( \E(\hbeta_2)=\bbeta_2 \),
\( \E(\hbeta^{\mathrm{S}})-\bbeta=-\mathbf{G}\bbeta_2 \). Expand
\( \Cov(\hbeta)=\E\bigl[(\hbeta^{\mathrm{S}}-\bbeta+\mathbf{G}\hbeta_2)(\hbeta^{\mathrm{S}}-\bbeta+\mathbf{G}\hbeta_2)\T\bigr] \). The
cross term is
\( \E\bigl[(\hbeta^{\mathrm{S}}-\bbeta)\hbeta_2\T\bigr]=\Cov(\hbeta^{\mathrm{S}},\hbeta_2)+\bigl(\E\hbeta^{\mathrm{S}}-\bbeta\bigr)\bbeta_2\T=-\mathbf{G}\bbeta_2\bbeta_2\T \),
and \( \E(\hbeta_2\hbeta_2\T)=\sigma^2\mathbf{C}+\bbeta_2\bbeta_2\T \). So
\[
\Cov(\hbeta)=\E\bigl[(\hbeta^{\mathrm{S}}-\bbeta)(\hbeta^{\mathrm{S}}-\bbeta)\T\bigr]
+\mathbf{G}\bigl(\sigma^2\mathbf{C}+\bbeta_2\bbeta_2\T-2\bbeta_2\bbeta_2\T\bigr)\mathbf{G}\T .
\]
(b) Put \( \A=\sigma^2\mathbf{C}-\bbeta_2\bbeta_2\T \). The claimed inequality for all \( \blambda \) says that
\( \mathbf{G}\A\mathbf{G}\T \) is nonnegative definite. The matrix \( \mathbf{L}=[\mathbf{0},\I_{p_2}] \) satisfies \( \mathbf{L}\mathbf{G}=\I \), so
\( \A=\mathbf{L}(\mathbf{G}\A\mathbf{G}\T)\mathbf{L}\T \). Hence \( \mathbf{G}\A\mathbf{G}\T \) is nonnegative definite iff \( \A \) is, and by
@lem-dep-rank-one iff \( \bbeta_2\T\mathbf{C}^{-1}\bbeta_2\le\sigma^2 \), which is \( \gamma\le1 \). If \( \gamma>1 \), the
same lemma gives \( \mathbf{a} \) with \( \mathbf{a}\T\A\mathbf{a}<0 \), and \( \blambda=\mathbf{L}\T\mathbf{a} \) satisfies
\( \mathbf{G}\T\blambda=\mathbf{a} \), so \( \blambda\T\mathbf{G}\A\mathbf{G}\T\blambda<0 \).
(c) \( \X\hbeta^{\mathrm{S}}=\M_1\Y \), whose mean is \( \M_1\X\bbeta \) and covariance \( \sigma^2\M_1 \). So
\( \E\norm{\M_1\Y-\X\bbeta}^2=\sigma^2p_1+\norm{(\I-\M_1)\X\bbeta}^2=\sigma^2(p_1+\gamma) \). The long fit
is @thm-dep-overfit(c), whose proof did not use \( \bbeta_2=\bzero \) for this part.
:::

The threshold in (b) depends only on the noncentrality, not on the design or on \( n \). It covers all linear
functions at once; for a particular function the short estimator can win for larger \( \gamma \), as (c) shows. Since
\( \gamma \) is unknown, neither criterion applies directly. But the \( F \) statistic for \( \bbeta_2=\bzero \) depends on the
parameters only through \( \gamma \) (@thm-glh-f-test), so a small \( F \) is evidence that \( \gamma \) is small. Turning this
into a procedure, and paying for having chosen the model from the data, is the subject of Chapter 29.

::: {#exm-dep-short-wins}
[A biased estimator that wins]

Take \( n=30 \) observations of two centred regressors with unit sample variances and sample
correlation \( r=0.9 \), an intercept, \( \sigma=1 \) and \( \beta_2=0.301 \), chosen so that
\( \gamma=0.5 \). Here \( \boldsymbol{\Pi} \) has slope entry \( r \), so the short slope on \( x_1 \) has bias
\( r\beta_2=0.271 \). Its standard deviation is \( 0.186 \), against
\( 0.426 \) for the long slope, and the ratio of mean squared error to long
variance is \( 0.595 \). A simulation with \( 200000 \) data sets
gives \( 0.594 \).

For two standardized regressors the ratio has a closed form (@exr-dep-two-regressors):
\[
\frac{\E(\hat{\beta}_1^{\mathrm{S}}-\beta_1)^2}{\Var(\hat{\beta}_1)}=1-r^2+r^2\gamma=1+r^2(\gamma-1),
\]{#eq-dep-two-regressor-ratio}

which is \( 1-0.81+0.81\times0.5=0.595 \). The short slope wins iff \( \gamma<1 \), as @thm-dep-mse(b)
requires, and the stakes, gain or loss, are proportional to \( r^2 \). This is panel (a) of
[Figure 19.1.1](#fig-dep-underfitting). The short model's \( s_1^2 \) has mean
\( 1+\gamma/(n-2)=1.0179 \), a bias too small to notice.
:::

```{.python .run #cell-underfitting-design}
import numpy as np
rng = np.random.default_rng(1901)
n, r, sigma = 30, 0.9, 1.0

# two centred, unit-variance regressors with sample correlation exactly r:
# orthonormalize (1, z1, z2) and keep the last two columns, which are centred
Z = np.linalg.qr(np.column_stack([np.ones(n), rng.normal(size=(n, 2))]))[0][:, 1:]
Z *= np.sqrt(n - 1)
x1 = Z[:, 0]
x2 = r * Z[:, 0] + np.sqrt(1 - r**2) * Z[:, 1]
X1 = np.column_stack([np.ones(n), x1])                # the short model
X = np.column_stack([X1, x2])                         # the long model

M1 = X1 @ np.linalg.solve(X1.T @ X1, X1.T)
S22_1 = x2 @ (np.eye(n) - M1) @ x2                    # ||(I - M1) x2||^2
beta2 = np.sqrt(0.5 * sigma**2 / S22_1)               # chosen so that gamma = 0.5
gamma = beta2**2 * S22_1 / sigma**2

var_long = sigma**2 * np.linalg.inv(X.T @ X)[1, 1]    # Var of the long slope on x1
Pi = np.linalg.solve(X1.T @ X1, X1.T @ x2)[1]         # slope of x2 on x1
var_short = sigma**2 * np.linalg.inv(X1.T @ X1)[1, 1]
mse_short = var_short + (Pi * beta2) ** 2             # variance + squared bias
print(f"gamma = {gamma:.3f}, MSE(short)/Var(long) = {mse_short / var_long:.4f}")
print(f"formula 1 - r^2 + r^2 gamma      = {1 - r**2 + r**2 * gamma:.4f}")
```

```{.python .run #cell-underfitting-simulate}
reps = 200_000
beta = np.array([1.0, 1.0, beta2])
Y = X @ beta + sigma * rng.normal(size=(reps, n))
b_long = np.linalg.solve(X.T @ X, X.T @ Y.T)[1]
b_short = np.linalg.solve(X1.T @ X1, X1.T @ Y.T)[1]
print("simulated MSE ratio:", np.mean((b_short - 1) ** 2) / np.mean((b_long - 1) ** 2))
```

::: {#exm-dep-polynomial-risk}
[Choosing a polynomial degree]

Part (c) of @thm-dep-mse applies to any nested sequence of models. Take \( n=40 \) equally spaced points on
\( [-1,1] \), the mean \( \mu(x)=e^x\sin 3x \), which is not a polynomial, and \( \sigma=0.25 \). A
polynomial of degree \( d \) has \( p=d+1 \) columns, and its risk at the design points is
\( \sigma^2(d+1)+\norm{(\I-\M_d)\bmu}^2 \). Panel (b) of [Figure 19.1.1](#fig-dep-underfitting) shows the two
terms. In units of
\( \sigma^2 \) the risk is \( 161.41 \) for a straight line, falls to
\( 5.52 \) at degree \( 4 \), and climbs back to
\( 11.00 \) at degree \( 10 \), where the fit is almost unbiased and pays only
for its variance. @exr-dep-add-one-term shows that adding one term lowers the risk iff that term's own
noncentrality exceeds \( 1 \).
:::

::: {when-format="html"}
![**Figure 19.1.1.** (a) Mean squared error of the short slope relative to the variance of the long
slope, @eq-dep-two-regressor-ratio, for three correlations. All lines cross at \( \gamma=1 \); the point is the
simulation of @exm-dep-short-wins. (b) Squared bias, variance and their sum for polynomial fits of increasing
degree (@exm-dep-polynomial-risk), in units of \( \sigma^2 \) and on a logarithmic scale.](underfitting.svg){#fig-dep-underfitting width=100%}
:::

::: {when-format="pdf"}
![(a) Mean squared error of the short slope relative to the variance of the long
slope, @eq-dep-two-regressor-ratio, for three correlations. All lines cross at \( \gamma=1 \); the point is the
simulation of @exm-dep-short-wins. (b) Squared bias, variance and their sum for polynomial fits of increasing
degree (@exm-dep-polynomial-risk), in units of \( \sigma^2 \) and on a logarithmic scale.](underfitting.pdf){width=100%}
:::

::: {.idea}
Omitting a regressor that matters costs bias; keeping one that does not costs variance, and for \( \gamma<1 \) the
bias is the cheaper. This does not make the biased coefficient *interpretable*: it still estimates \( \bbeta_1+\boldsymbol{\Pi}\bbeta_2 \), and whether that
answers the question is a matter for [Chapter 25](../ch25-causal-interpretation/index.html).
:::

## Exercises

### A. Check your understanding

::: {#exr-dep-line-quadratic}
[A1]

A straight line is fitted at \( x=0,1,2,3 \) when the true mean is \( \beta_0+\beta_1x+\beta_2x^2 \). Find the biases of
the fitted intercept and slope.
:::

::: {#exr-dep-uncorrelated-s2}
[A2]

True or false: if the omitted regressors are orthogonal to the included ones (\( \X_1\T\X_2=\mathbf{0} \)), the
short regression is harmless, including its estimate of \( \sigma^2 \).
:::

### B. Practice

::: {#exr-dep-line-threshold}
[B1]

In @exr-dep-line-quadratic, show that \( \norm{(\I-\M_1)\X_2}^2=4 \). For which values of \( \beta_2/\sigma \) does
the straight line have a smaller mean squared error matrix than the quadratic fit, in the sense of
@thm-dep-mse(b)? For which values is it better at the four design points?
:::

::: {#exr-dep-two-regressors}
[B2]

Derive @eq-dep-two-regressor-ratio. Take an intercept and two centred regressors with
\( \norm{\x_1}^2=\norm{\x_2}^2=S \) and sample correlation \( r \). Show that \( \Var(\hat{\beta}_1)=\sigma^2/\bigl(S(1-r^2)\bigr) \), that
the short slope has variance \( \sigma^2/S \) and bias \( r\beta_2 \), and that \( \gamma=\beta_2^2S(1-r^2)/\sigma^2 \).
:::

::: {.solution}
The long variance is @eq-proj-vif-preview with \( R_1^2=r^2 \). The short slope is \( \x_1\T\Y/S \), with variance
\( \sigma^2/S \) and mean \( \beta_1+\beta_2\x_1\T\x_2/S=\beta_1+r\beta_2 \). The residual of \( \x_2 \) after regression on
\( (\bone,\x_1) \) has squared length \( S(1-r^2) \), which gives \( \gamma \). Then
\( \bigl(\sigma^2/S+r^2\beta_2^2\bigr)S(1-r^2)/\sigma^2=1-r^2+r^2\gamma \).
:::

::: {#exr-dep-add-one-term}
[B3]

Let \( \C(\X_1)\subset\C(\X_2)\subset\cdots \) be nested model spaces with dimensions \( p_1<p_2<\cdots \), and let
\( \bmu=\E(\Y) \) be arbitrary. Show that the risk \( \E\norm{\M_k\Y-\bmu}^2 \) of model \( k \) exceeds that of model \( k+1 \)
iff \( \gamma_k=\norm{(\M_{k+1}-\M_k)\bmu}^2/\sigma^2>p_{k+1}-p_k \).
:::

::: {.solution}
As in @thm-dep-mse(c), the risk of model \( k \) is \( \sigma^2p_k+\norm{(\I-\M_k)\bmu}^2 \). By @thm-proj-nested,
\( \I-\M_k=(\I-\M_{k+1})+(\M_{k+1}-\M_k) \) with orthogonal ranges, so
\( \norm{(\I-\M_k)\bmu}^2=\norm{(\I-\M_{k+1})\bmu}^2+\sigma^2\gamma_k \). The difference of the two risks is
\( \sigma^2\bigl(\gamma_k-(p_{k+1}-p_k)\bigr) \).
:::

### C. Going deeper

::: {#exr-dep-s2-mse}
[C1]

Assume normal errors. Compare \( s_1^2 \) (short) and \( s^2 \) (long) as estimators of \( \sigma^2 \) by mean squared
error, when \( \bbeta_2\neq\bzero \). Show that \( s_1^2 \) is better iff
\( 2(n-p_1)+4\gamma+\gamma^2<2(n-p_1)^2/(n-p) \), and interpret the condition when \( n \) is large.
:::
