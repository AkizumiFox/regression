# Best linear prediction

With random rows, a regression is a rule for predicting the response of a new case from its
regressors. The population theory is in place: the conditional mean \( m(\mathbf{X})=\E(Y\mid\mathbf{X}) \) is the best
predictor (@prp-proj-conditional-expectation), the best *linear* predictor
\( L(Y\mid\mathbf{X})=\alpha^*+\bbeta^{*\top}\mathbf{X} \) is a projection (@thm-proj-blp, @thm-rv-blp), the two coincide
for jointly normal variables (@prp-mvn-best-predictor), and least squares estimates \( L \) consistently (@prp-proj-consistency). This section asks the finite-sample questions: how large is the error of the
estimated predictor on a new case, and how misleading is its error on the data used to fit it?

## Characterizing the best predictor

The next result gives a test for whether a predictor is the best one, and a correlation that
measures how well the best predictor can do.

::: {#prp-cor-best-predictor}
[The best predictor and the correlation ratio]

Let \( \E Y^2<\infty \), \( \Var(Y)>0 \), and \( m(\mathbf{X})=\E(Y\mid\mathbf{X}) \).

::: {.enumerate options="label=(\alph*)"}
1. A predictor \( g(\mathbf{X}) \) with \( \E g(\mathbf{X})^2<\infty \) equals \( m(\mathbf{X}) \) with probability one iff
   \( \E g(\mathbf{X})=\E Y \) and \( \Cov\bigl(h(\mathbf{X}),\,Y-g(\mathbf{X})\bigr)=0 \) for every square-integrable \( h \).

2. For every such \( g \) with \( \Var g(\mathbf{X})>0 \),
   \[
\operatorname{corr}\bigl(Y,g(\mathbf{X})\bigr)^2\le\eta^2=\frac{\Var m(\mathbf{X})}{\Var Y},
\]{#eq-cor-correlation-ratio}

   with equality for \( g=m \) when \( \Var m(\mathbf{X})>0 \). The **correlation ratio** \( \eta^2 \) exceeds the
   squared multiple correlation \( \rho^2_{Y\cdot X} \) of @prp-rv-multiple-correlation by
   \( \E\bigl(m(\mathbf{X})-L(Y\mid\mathbf{X})\bigr)^2/\Var Y \), so \( \eta^2=\rho^2_{Y\cdot X} \) iff \( m \) is affine.
:::

:::

::: {.proof}
(a) If \( g=m \), then \( \E m(\mathbf{X})=\E Y \) by the tower rule, and \( \E[(Y-m)h]=0 \) by
@prp-proj-conditional-expectation(a), so the covariance is zero. Conversely, let \( d=g-m \). Then
\( \E d=0 \) and, taking \( h=d \),
\( \Var d=\Cov(d,Y-m)-\Cov(d,Y-g)=0-0 \), so \( d=0 \) with probability one.

(b) By (a) with \( h=g \), \( \Cov(Y,g)=\Cov(m,g) \). Cauchy–Schwarz for covariances gives
\( \Cov(m,g)^2\le\Var m\,\Var g \), hence
\( \operatorname{corr}(Y,g)^2=\Cov(m,g)^2/(\Var Y\,\Var g)\le\Var m/\Var Y \). For \( g=m \) the covariance is
\( \Var m \) and equality holds. For the last claim write \( L=L(Y\mid\mathbf{X})=L(m\mid\mathbf{X}) \) (@prp-proj-conditional-expectation(b)). The error \( m-L \) has mean zero and is uncorrelated with the
affine function \( L \), so \( \Var m=\Var L+\E(m-L)^2 \). Divide by \( \Var Y \) and use
\( \Var L/\Var Y=\rho^2_{Y\cdot X} \) (@eq-proj-population-pythagoras).
:::

Part (a) is the population basis of residual plots. If a fitted model reproduces the conditional
mean, its residuals are uncorrelated with *every* function of the regressors, not only with the
regressors themselves. A residual plot that shows a trend
against \( x_j^2 \), against a product \( x_jx_k \) or against the fitted values is evidence that the fit is not
\( m \). Sections [14.6](06-lack-of-fit.html) and [14.7](07-near-replicates.html) turn this into tests.
Part (b) says that no transformation of the regressors predicts better, in the correlation sense,
than the conditional mean. The gap \( \eta^2-\rho^2_{Y\cdot X} \) is the share of the variance of \( Y \) that a linear
predictor leaves on the table.

## The error of an estimated predictor

In practice \( \alpha^* \) and \( \bbeta^* \) are replaced by least squares estimates from a training sample
\( (\mathbf{X}_i,Y_i) \), \( i=1,\dots,n \). By centring and @thm-proj-fwl these are the plug-in estimates
\( \hbeta=\bS_{XX}^{-1}\mathbf{s}_{XY} \) and \( \hat{\alpha}=\bar Y-\hbeta\T\bar{\mathbf{X}} \), built from the sample means and
covariances. Their performance is judged on a new case \( (\mathbf{X}_0,Y_0) \) drawn from the same
distribution, independently of the training data.

::: {#prp-cor-prediction-decomposition}
[Excess error of a fitted linear predictor]

Let \( (\mathbf{X}_0,Y_0) \) have finite second moments, \( \Cov(\mathbf{X}_0)=\bSigma_{XX} \), mean \( \bmu_X \), and best linear
predictor \( L_0=\alpha^*+\bbeta^{*\top}\mathbf{X}_0 \) with error variance \( \sigma_*^2=\E(Y_0-L_0)^2 \). Let
\( \hat Y_0=a+\bb\T\mathbf{X}_0 \), where \( (a,\bb) \) depend only on training data independent of
\( (\mathbf{X}_0,Y_0) \). Then, given the training data,
\[
\begin{aligned}
\E\bigl[(Y_0-\hat Y_0)^2\bigm|\text{training data}\bigr]
&=\sigma_*^2+\bigl(\alpha^*-a+(\bbeta^*-\bb)\T\bmu_X\bigr)^2\\
&\quad+(\bbeta^*-\bb)\T\bSigma_{XX}(\bbeta^*-\bb).
\end{aligned}
\]
:::

::: {.proof}
Condition on the training data, so that \( a \) and \( \bb \) are constants. Write
\( Y_0-\hat Y_0=(Y_0-L_0)+(L_0-\hat Y_0) \). The second term is an affine function of \( \mathbf{X}_0 \), and by
@thm-proj-blp the first has mean zero and is uncorrelated with every affine function of \( \mathbf{X}_0 \).
So the cross term vanishes, and the mean square of the affine term is its squared mean plus its
variance \( (\bbeta^*-\bb)\T\bSigma_{XX}(\bbeta^*-\bb) \).
:::

The first term is the error of the best linear predictor, which no amount of data can remove. The
other two are the price of estimation, and by @prp-proj-consistency they vanish as \( n\to\infty \)
whatever the joint distribution. To see how large they are for moderate \( n \), we need the
distribution of the regressors. Under joint normality the calculation can be done exactly. It uses
one new distribution.

::: {#def-cor-wishart}
[Wishart distribution]

Let the rows of the \( m\times k \) matrix \( \Z \) be independent \( \Normal_k(\bzero,\bSigma) \) vectors. The law of
\( \W=\Z\T\Z=\sum_{i=1}^m\bz_i\bz_i\T \) is the **Wishart distribution** \( W_k(m,\bSigma) \) with \( m \) degrees of
freedom. For \( k=1 \), \( W_1(m,\sigma^2) \) is the law of \( \sigma^2\chi^2(m) \).
:::

::: {#lem-cor-centred-wishart}
[Centred cross-products]

Let \( \mathbf{X}_1,\dots,\mathbf{X}_n \) (\( n\ge2 \)) be independent \( \Normal_k(\bmu,\bSigma) \) vectors with mean \( \bar{\mathbf{X}} \), and
let \( \W=\sum_i(\mathbf{X}_i-\bar{\mathbf{X}})(\mathbf{X}_i-\bar{\mathbf{X}})\T \). Then \( \W\sim W_k(n-1,\bSigma) \), and \( \W \) is
independent of \( \bar{\mathbf{X}} \).
:::

::: {.proof}
Let \( \X_1 \) be the \( n\times k \) matrix with rows \( \mathbf{X}_i\T \). Write \( \X_1=\bone\bmu\T+\Z_0\bSigma^{1/2} \), where
\( \bSigma^{1/2} \) is the symmetric square root (@thm-mat-square-root) and \( \Z_0 \) has independent standard normal
entries. Let \( \mathbf{O} \) be an orthogonal \( n\times n \) matrix with first row \( n^{-1/2}\bone\T \), for instance the Helmert
matrix of @exm-mat-helmert, and let \( \mathbf{K} \) consist of its last \( n-1 \) rows. From \( \mathbf{O}\T\mathbf{O}=\I \),
\( \mathbf{K}\T\mathbf{K}=\I-n^{-1}\bone\bone\T \), so \( \W=\X_1\T\mathbf{K}\T\mathbf{K}\X_1=(\mathbf{K}\X_1)\T(\mathbf{K}\X_1) \), and \( \mathbf{K}\bone=\bzero \) gives
\( \mathbf{K}\X_1=(\mathbf{K}\Z_0)\bSigma^{1/2} \).

Each column of \( \mathbf{O}\Z_0 \) is \( \Normal_n(\bzero,\mathbf{O}\mathbf{O}\T)=\Normal_n(\bzero,\I) \) by @thm-mvn-linear, and the
columns are independent, so \( \mathbf{O}\Z_0 \) again has independent standard normal entries. Its first row is
\( \sqrt n\,(\bar{\mathbf{X}}-\bmu)\T\bSigma^{-1/2} \), and its other rows form \( \mathbf{K}\Z_0 \). Hence the \( n-1 \) rows of
\( \mathbf{K}\X_1 \) are independent \( \Normal_k(\bzero,\bSigma) \) vectors, independent of \( \bar{\mathbf{X}} \), and
\( \W\sim W_k(n-1,\bSigma) \).
:::

::: {#lem-cor-inverse-wishart}
[Mean of an inverse Wishart matrix]

Let \( \W\sim W_k(m,\bSigma) \) with \( \bSigma \) positive definite and \( m\ge k+2 \). Then \( \W \) is nonsingular with
probability one, and
\[
\E\bigl(\W^{-1}\bigr)=\frac{\bSigma^{-1}}{m-k-1}.
\]
:::

::: {.proof}
*Reduction.* Write \( \Z=\Z_0\bSigma^{1/2} \) with \( \Z_0 \) as in the previous proof. Then
\( \W=\bSigma^{1/2}\W_0\bSigma^{1/2} \) with \( \W_0=\Z_0\T\Z_0\sim W_k(m,\I) \), and
\( \W^{-1}=\bSigma^{-1/2}\W_0^{-1}\bSigma^{-1/2} \). It suffices to show \( \E(\W_0^{-1})=\I/(m-k-1) \).

*Nonsingularity.* \( \W_0 \) is singular iff some column \( \bz_j \) of \( \Z_0 \) lies in the span of the others.
Given the other columns, that span has dimension at most \( k-1<m \), and \( \bz_j\sim\Normal_m(\bzero,\I) \),
which has a density, falls in it with probability zero.

*Diagonal entries.* Partition \( \Z_0=[\bz_1,\Z_2] \). By @thm-mat-partitioned-inverse, the \( (1,1) \) entry of
\( \W_0^{-1} \) is \( 1/(\bz_1\T\bz_1-\bz_1\T\Z_2(\Z_2\T\Z_2)^{-1}\Z_2\T\bz_1)=1/\bz_1\T(\I-\bP_2)\bz_1 \), where \( \bP_2 \) projects onto
\( \C(\Z_2) \), of dimension \( k-1 \). Given \( \Z_2 \), \( \bz_1 \) is still \( \Normal_m(\bzero,\I) \), so
\( \bz_1\T(\I-\bP_2)\bz_1\sim\chi^2(m-k+1) \) by @thm-qf-chisq(a). This law does not depend on \( \Z_2 \), so it
holds unconditionally (@lem-cor-conditioning), and @lem-cor-inverse-chisq gives mean
\( 1/(m-k-1) \), finite because \( m-k+1>2 \). The same argument applies to every diagonal entry.

*Off-diagonal entries.* Let \( \bD_i \) be the identity with its \( i \)th diagonal entry replaced by \( -1 \).
Then \( \Z_0\bD_i \) has the same law as \( \Z_0 \), so \( \bD_i\W_0^{-1}\bD_i=(\bD_i\W_0\bD_i)^{-1} \) has the same law as
\( \W_0^{-1} \). The mean of an off-diagonal entry \( (i,j) \) therefore equals its own negative, and is zero. It
exists because \( \lvert w^{ij}\rvert\le\tfrac12(w^{ii}+w^{jj}) \) for the entries of the positive definite
\( \W_0^{-1} \).
:::

With these two lemmas the random-design analogue of @thm-lm-moments follows at once.

::: {#cor-cor-slope-covariance}
[Covariance of the slopes under normal regressors]

Let the rows \( (\mathbf{X}_i\T,Y_i) \) be independent draws from a \( (k+1) \)-variate normal distribution with
\( \bSigma_{XX} \) positive definite and conditional variance \( \sigma^2=\Var(Y\mid\mathbf{X}) \), and let \( \hbeta \) be the least
squares slopes from a fit with an intercept. If \( n\ge k+3 \),
\[
\E(\hbeta)=\bbeta,\qquad \Cov(\hbeta)=\frac{\sigma^2}{n-k-2}\,\bSigma_{XX}^{-1}.
\]
:::

::: {.proof}
Given the regressors, the slopes are unbiased with covariance \( \sigma^2\W^{-1} \), where \( \W \) is the
centred cross-product matrix (@eq-proj-partitioned-inverse with \( \X_2=\bone \)). By @lem-cor-centred-wishart
\( \W\sim W_k(n-1,\bSigma_{XX}) \), and by @prp-rv-total-covariance and @lem-cor-inverse-wishart the
unconditional covariance is \( \sigma^2\E(\W^{-1})=\sigma^2\bSigma_{XX}^{-1}/(n-k-2) \).
:::

For \( k=1 \) this is @eq-cor-slope-variance. The factor \( 1/(n-k-2) \), against \( 1/(n-1) \) for a fixed design
with the average spread, grows quickly with \( k \).

::: {#thm-cor-prediction-error}
[Prediction error under normal regressors]

In the setting of @cor-cor-slope-covariance, let \( (\mathbf{X}_0,Y_0) \) be a further independent draw and
\( \hat Y_0=\hat{\alpha}+\hbeta\T\mathbf{X}_0 \) the fitted predictor. If \( n\ge k+3 \), then
\[
\E(Y_0-\hat Y_0)^2=\sigma^2\Bigl(1+\frac1n\Bigr)\frac{n-2}{n-k-2}.
\]{#eq-cor-prediction-error}

:::

::: {.proof}
Write \( \hat Y_0=\bar Y+\hbeta\T(\mathbf{X}_0-\bar{\mathbf{X}}) \), and use the model
\( Y_i=\alpha+\bbeta\T\mathbf{X}_i+e_i \) of @eq-mvn-random-x-model for \( i=0,\dots,n \). Since
\( \bar Y=\alpha+\bbeta\T\bar{\mathbf{X}}+\bar e \),
\[
Y_0-\hat Y_0=e_0-\bar e-(\hbeta-\bbeta)\T(\mathbf{X}_0-\bar{\mathbf{X}}).
\]
Condition on \( \mathbf{X}_0,\dots,\mathbf{X}_n \). Then \( e_0 \), \( \bar e \) and \( \hbeta-\bbeta=\W^{-1}\sum_i(\mathbf{X}_i-\bar{\mathbf{X}})e_i \) are
normal with mean zero, and they are uncorrelated: \( e_0 \) is independent of the training errors, and
\( \Cov(\bar e,\hbeta)=\sigma^2n^{-1}\W^{-1}\sum_i(\mathbf{X}_i-\bar{\mathbf{X}})=\bzero \). So the conditional mean square is
\[
\sigma^2\Bigl[1+\frac1n+(\mathbf{X}_0-\bar{\mathbf{X}})\T\W^{-1}(\mathbf{X}_0-\bar{\mathbf{X}})\Bigr].
\]
Now \( \mathbf{u}=\mathbf{X}_0-\bar{\mathbf{X}}\sim\Normal_k\bigl(\bzero,(1+1/n)\bSigma_{XX}\bigr) \), and \( \mathbf{u} \) is independent of \( \W \),
because \( \mathbf{X}_0 \) is independent of the training sample and \( \bar{\mathbf{X}} \) is independent of \( \W \) (@lem-cor-centred-wishart). Hence
\( \E(\mathbf{u}\T\W^{-1}\mathbf{u})=\tr\bigl[\E(\W^{-1})\E(\mathbf{u}\mathbf{u}\T)\bigr]=(1+1/n)\,k/(n-k-2) \) by @lem-cor-inverse-wishart, and
\[
\E(Y_0-\hat Y_0)^2=\sigma^2\Bigl(1+\frac1n\Bigr)\Bigl(1+\frac{k}{n-k-2}\Bigr)
=\sigma^2\Bigl(1+\frac1n\Bigr)\frac{n-2}{n-k-2}.
\]
:::

The conditional mean square in the proof, \( \sigma^2[1+1/n+(\mathbf{X}_0-\bar{\mathbf{X}})\T\W^{-1}(\mathbf{X}_0-\bar{\mathbf{X}})] \), is the
variance behind the prediction interval of @thm-ci-prediction-interval, so that interval is exact
for random regressors as well (@thm-cor-conditional). What the theorem adds is the average over
new cases and over training designs.

## The optimism of the in-sample error

The residual mean square \( \text{SSE}/n \) is the error the fitted predictor makes on the cases used to
fit it. It is a biased guide to the error on new cases, in both the fixed and the random design.

::: {#prp-cor-optimism}
[In-sample and out-of-sample error]

::: {.enumerate options="label=(\alph*)"}
1. *(Fixed design.)* Let \( \E(\Y)=\X\bbeta \) and \( \Cov(\Y)=\sigma^2\I_n \) with \( \rank(\X)=r \), and let \( \Y^* \) be a
   vector of new responses at the same design, uncorrelated with \( \Y \), with the same mean and
   covariance. Then
   \[
\E\,\frac{\norm{\Y-\X\hbeta}^2}{n}=\sigma^2\,\frac{n-r}{n},\qquad
\E\,\frac{\norm{\Y^*-\X\hbeta}^2}{n}=\sigma^2\,\frac{n+r}{n},
\]
   so the in-sample error understates the error on new responses by \( 2r\sigma^2/n \) on average.

2. *(Normal random design.)* In the setting of @thm-cor-prediction-error, the in-sample error has
   \( \E(\text{SSE}/n)=\sigma^2(n-k-1)/n \), and the error on a new case is @eq-cor-prediction-error.
:::

:::

::: {.proof}
(a) The first identity is \( \E(\text{SSE})=\sigma^2(n-r) \) (@thm-lm-sigma2, or @eq-rv-rss-mean with a mean in
\( \C(\X) \)). For the second, \( \Y^*-\X\hbeta=(\Y^*-\X\bbeta)-\M(\Y-\X\bbeta) \), and the two terms are uncorrelated,
so the expected squared length is \( n\sigma^2+\sigma^2\tr(\M)=\sigma^2(n+r) \) by @thm-rv-quadform-mean.
(b) Given the regressors, \( \E(\text{SSE}\mid\X)=\sigma^2(n-k-1) \), which does not depend on \( \X \); take
expectations.
:::

A new random case brings its own, possibly outlying, regressor values, where the fitted plane is least
reliable, so the random design is worse. [Figure 14.2.1](#fig-cor-prediction-error) shows the three
expected errors for \( n=30 \). With \( k=10 \) regressors the in-sample error is
\( 0.633\sigma^2 \), the error on new responses at the same design is \( 1.367\sigma^2 \), and the error on a new
random case is \( 1.607\sigma^2 \): the in-sample figure flatters the fit by a factor of about
\( 2.5 \). At \( k=20 \) the three values are \( 0.300 \), \( 1.700 \) and \( 3.617 \), in units of
\( \sigma^2 \). The simulated averages over \( 20000 \) training samples per point, for instance
\( 1.604 \) at \( k=10 \) and \( 3.612 \) at \( k=20 \), agree with @eq-cor-prediction-error.
Estimating the out-of-sample error honestly, by penalties such as Mallows' \( C_p \) or by
cross-validation, is the subject of Chapter 29.

::: {when-format="html"}
![**Figure 14.2.1.** Expected squared error of a least squares fit, \( n=30 \): in-sample (bottom), new
responses at the same design (middle) and a new case with normal regressors (top, @eq-cor-prediction-error). Points: simulation.](prediction_error.svg){#fig-cor-prediction-error width=72%}
:::

::: {when-format="pdf"}
![Expected squared error of a least squares fit, \( n=30 \): in-sample (bottom), new
responses at the same design (middle) and a new case with normal regressors (top, @eq-cor-prediction-error). Points: simulation.](prediction_error.pdf){width=72%}
:::

```{.python .run #cell-prediction-error-msep}
import numpy as np
rng = np.random.default_rng(1422)
n, sigma2, reps = 30, 1.0, 20_000

def simulate(k):
    """Average in-sample and out-of-sample squared error over many training samples."""
    ins, out = np.empty(reps), np.empty(reps)
    beta = np.full(k, 0.5)
    for s in range(reps):
        X = np.column_stack([np.ones(n), rng.normal(size=(n, k))])
        y = X @ np.r_[1.0, beta] + rng.normal(size=n)
        b, *_ = np.linalg.lstsq(X, y, rcond=None)
        ins[s] = np.mean((y - X @ b) ** 2)           # SSE / n
        # new point X0 ~ N(0, I), Y0 = 1 + beta^T X0 + e0: exact conditional error
        out[s] = sigma2 + (b[0] - 1.0) ** 2 + np.sum((b[1:] - beta) ** 2)
    return ins.mean(), out.mean()

results = {}
for k in (1, 5, 10, 20):
    ins, out = simulate(k)
    results[k] = (ins, out)
    theory = sigma2 * (1 + 1 / n) * (n - 2) / (n - k - 2)
    print(f"k = {k:2d}: in-sample {ins:.3f} (theory {sigma2 * (n - k - 1) / n:.3f}),"
          f" new point {out:.3f} (theory {theory:.3f})")
```

::: {#exm-cor-rand-prediction}
[Predicting doctor visits]

The RAND Health Insurance Experiment records, for \( 20190 \) person-years, the number of
outpatient visits and characteristics of the insurance plan and of the person. Treat these rows as
the population, and predict \( \log(1+\text{visits}) \) linearly from \( 6 \) of the characteristics
(coinsurance rate, deductible and participation-incentive terms, maximum expenditure, physical
limitation and a chronic-disease index). The population projection has mean squared error
\( 0.6330 \) and explains a fraction \( 0.094 \) of the variance. Draw training samples of
\( 40 \) rows, fit by least squares, and evaluate each fit exactly on the whole population. Over
\( 4000 \) training samples the average in-sample error is \( 0.5182 \) and the average error on
new cases is \( 0.8201 \), against normal-theory values \( 0.5223 \) and \( 0.7705 \) from @prp-cor-optimism(b). The normal-theory value understates the error on new cases: the regressors are far from normal (two are indicators, one skewed), the
linear predictor need not be the conditional mean, and the error variance need not be constant; the
formula assumes none of this. It is a first guide to the optimism, not a substitute for checking it.
:::

```{.python .run #cell-prediction-error-rand}
import statsmodels.api as sm
data = sm.datasets.randhie.load_pandas().data
cols = ["lncoins", "idp", "lpi", "fmde", "physlm", "disea"]
Xpop = np.column_stack([np.ones(len(data)), data[cols].to_numpy()])
ypop = np.log1p(data["mdvis"].to_numpy())            # log(1 + number of visits)
b_pop, *_ = np.linalg.lstsq(Xpop, ypop, rcond=None)  # the population projection
sigma2_pop = np.mean((ypop - Xpop @ b_pop) ** 2)     # its mean squared error

rng = np.random.default_rng(1423)
n_train, draws = 40, 4000
k_rand = len(cols)
ins, out = np.empty(draws), np.empty(draws)
for s in range(draws):
    idx = rng.integers(0, len(ypop), n_train)        # iid draws from the population
    b, *_ = np.linalg.lstsq(Xpop[idx], ypop[idx], rcond=None)
    ins[s] = np.mean((ypop[idx] - Xpop[idx] @ b) ** 2)
    out[s] = np.mean((ypop - Xpop @ b) ** 2)         # exact error over the population
print(f"population projection error {sigma2_pop:.4f}")
print(f"in-sample {ins.mean():.4f}, new case {out.mean():.4f}")
print(f"normal-theory values: {sigma2_pop * (n_train - k_rand - 1) / n_train:.4f},"
      f" {sigma2_pop * (1 + 1 / n_train) * (n_train - 2) / (n_train - k_rand - 2):.4f}")
```

## Best linear unbiased prediction

So far the new response was independent of the training data. When it is not, as for the next
value of a time series or the yield at an unsampled site near sampled ones, the training
residuals carry information about the new error, and a good predictor should use it.

::: {#thm-cor-blup}
[Best linear unbiased predictor]

Let \( \Y=\X\bbeta+\be \) and \( y_0=\x_0\T\bbeta+e_0 \), with \( \E(\be)=\bzero \), \( \E(e_0)=0 \), \( \Cov(\be)=\sigma^2\V \)
for a known positive definite \( \V \), \( \Cov(\be,e_0)=\sigma^2\bv_0 \) and \( \Var(e_0)=\sigma^2v_{00} \). Suppose
\( \x_0\in\C(\X\T) \). A predictor \( \mathbf{a}\T\Y+a_0 \) is **unbiased** if \( \E(\mathbf{a}\T\Y+a_0-y_0)=0 \) for every
\( \bbeta \). Among unbiased linear predictors, the mean squared error \( \E(y_0-\mathbf{a}\T\Y-a_0)^2 \) is minimized
by
\[
\hat y_0=\x_0\T\tilde{\bbeta}+\bv_0\T\V^{-1}(\Y-\X\tilde{\bbeta}),
\]{#eq-cor-blup}

where \( \tilde{\bbeta} \) solves the generalized least squares equations \( \X\T\V^{-1}\X\bb=\X\T\V^{-1}\Y \). With
\( \mathbf{c}=\x_0-\X\T\V^{-1}\bv_0 \), its mean squared error is
\( \sigma^2\bigl[v_{00}-\bv_0\T\V^{-1}\bv_0+\mathbf{c}\T(\X\T\V^{-1}\X)\ginv\mathbf{c}\bigr] \).
:::

::: {.proof}
Put \( \boldsymbol{\delta}=\V^{-1}\bv_0 \) and split the prediction error of any linear predictor:
\[
y_0-\mathbf{a}\T\Y-a_0=\underbrace{\bigl(e_0-\boldsymbol{\delta}\T\be\bigr)}_{u}
+\underbrace{\bigl(\x_0\T\bbeta-\boldsymbol{\delta}\T\X\bbeta-(\mathbf{a}-\boldsymbol{\delta})\T\Y-a_0\bigr)}_{w}.
\]
The term \( u \) is uncorrelated with \( \be \), since \( \Cov(\be,e_0-\boldsymbol{\delta}\T\be)=\sigma^2(\bv_0-\V\boldsymbol{\delta})=\bzero \). The
term \( w \) is a constant plus a linear function of \( \Y \), hence of \( \be \), so \( \E(uw)=0 \) and
\( \E(y_0-\mathbf{a}\T\Y-a_0)^2=\E u^2+\E w^2 \). Only \( \E w^2 \) depends on the predictor.

Unbiasedness means \( \mathbf{a}\T\X\bbeta+a_0=\x_0\T\bbeta \) for all \( \bbeta \), that is, \( a_0=0 \) and \( \X\T\mathbf{a}=\x_0 \).
Equivalently, \( (\mathbf{a}-\boldsymbol{\delta})\T\Y \) is a linear unbiased estimator of \( \mathbf{c}\T\bbeta \), which is estimable
because \( \x_0 \) and \( \X\T\boldsymbol{\delta} \) lie in \( \C(\X\T) \). For such predictors \( \E w=0 \), so \( \E w^2 \) is the variance of
\( (\mathbf{a}-\boldsymbol{\delta})\T\Y \). By @cor-opt-aitken this variance is smallest when \( (\mathbf{a}-\boldsymbol{\delta})\T\Y=\mathbf{c}\T\tilde{\bbeta} \), the
generalized least squares estimator, with variance \( \sigma^2\mathbf{c}\T(\X\T\V^{-1}\X)\ginv\mathbf{c} \). The predictor is then
\( \boldsymbol{\delta}\T\Y+\mathbf{c}\T\tilde{\bbeta}=\x_0\T\tilde{\bbeta}+\boldsymbol{\delta}\T(\Y-\X\tilde{\bbeta}) \), which is @eq-cor-blup. Finally
\( \E u^2=\sigma^2(v_{00}-2\boldsymbol{\delta}\T\bv_0+\boldsymbol{\delta}\T\V\boldsymbol{\delta})=\sigma^2(v_{00}-\bv_0\T\V^{-1}\bv_0) \).
:::

The predictor is the estimated mean at the new point plus a prediction of the new error \( e_0 \) from
the residuals, by the best linear predictor of @thm-rv-blp. In the ordinary model, \( \V=\I \) and \( \bv_0=\bzero \), the correction vanishes and
@eq-cor-blup is the fitted value \( \x_0\T\hbeta \), with mean squared error
\( \sigma^2[v_{00}+\x_0\T(\X\T\X)\ginv\x_0] \). The same theorem, with \( \bv_0 \) built from a spatial covariance
function, is the kriging predictor of geostatistics, and with random effects in place of \( e_0 \) it
gives the predictors of mixed models (Chapter 32).

## Exercises

### A. Check your understanding

::: {#exr-cor-optimism-numbers}
[A1]

With normal regressors, \( n=50 \) and \( k=10 \), compute the expected in-sample error \( \E(\text{SSE}/n) \) and the
expected error on a new case, in units of \( \sigma^2 \). By what factor does the first understate the second?
:::

### B. Practice

::: {#exr-cor-linearized}
[B1]

For any predictor \( g(\mathbf{X}) \) with \( 0<\Var g<\infty \), let \( \tilde g=\E Y+b\,(g-\E g) \) with
\( b=\Cov(Y,g)/\Var g \). Show that \( \tilde g \) predicts at least as well as \( g \) and that
\( \E(Y-\tilde g)^2=\Var(Y)\bigl(1-\operatorname{corr}(Y,g)^2\bigr) \). Conclude that among such *linearized*
predictors, a higher squared correlation with \( Y \) means a smaller mean squared error, and that
@eq-cor-correlation-ratio identifies the best of them.
:::

::: {.solution}
\( \tilde g \) is the best linear predictor of \( Y \) from the single variable \( g \) (@thm-proj-blp with
\( k=1 \)), and \( g \) itself is one such linear function (intercept \( 0 \), slope \( 1 \)), so
\( \E(Y-\tilde g)^2\le\E(Y-g)^2 \). By @eq-proj-population-pythagoras the error of \( \tilde g \) is
\( \Var Y-\Cov(Y,g)^2/\Var g=\Var Y\,(1-\operatorname{corr}(Y,g)^2) \), which is decreasing in the squared correlation. The
largest squared correlation is \( \eta^2 \), attained by \( g=m \), whose linearization is \( m \) itself.
:::

::: {#exr-cor-ar1-forecast}
[B2]

Let the errors of a regression follow a stationary autoregression, \( e_t=\phi e_{t-1}+u_t \), \( \lvert\phi\rvert<1 \),
with white noise \( u_t \), so that \( \Cov(e_s,e_t)\propto\phi^{\lvert s-t\rvert} \). Observe \( t=1,\dots,n \) and predict
\( y_{n+1}=\x_{n+1}\T\bbeta+e_{n+1} \). Show that \( \V^{-1}\bv_0=\phi\,\vect{e}_n \), the last standard basis vector
scaled by \( \phi \), and hence that the BLUP is \( \x_{n+1}\T\tilde{\bbeta}+\phi\,(y_n-\x_n\T\tilde{\bbeta}) \).
:::

::: {.solution}
With \( \V_{st}=\phi^{\lvert s-t\rvert} \) (the common factor cancels in \( \V^{-1}\bv_0 \)), the covariance of \( e_s \) with
\( e_{n+1} \) is \( \phi^{n+1-s}=\phi\cdot\phi^{n-s}=\phi\,\V_{sn} \). So \( \bv_0=\phi\,\V\vect{e}_n \) and \( \V^{-1}\bv_0=\phi\,\vect{e}_n \).
Substituting into @eq-cor-blup, the correction term is \( \phi \) times the last residual
\( y_n-\x_n\T\tilde{\bbeta} \). Only the most recent residual matters, because given \( e_n \) the earlier errors carry
no further linear information about \( e_{n+1} \).
:::

::: {#exr-cor-unbiased-cov}
[B3]

In the setting of @cor-cor-slope-covariance, let \( s^2=\text{SSE}/(n-k-1) \) and let \( \W \) be the centred
cross-product matrix. Show that \( s^2\W^{-1} \), the usual estimate of the covariance of the slopes, is an
unbiased estimator of their *unconditional* covariance \( \sigma^2\bSigma_{XX}^{-1}/(n-k-2) \).
:::

::: {.solution}
Given the regressors, \( (n-k-1)s^2/\sigma^2\sim\chi^2(n-k-1) \) (@thm-opt-sampling(b)), a law free of the
design, so by @lem-cor-conditioning \( s^2 \) is independent of the regressors and hence of \( \W \). Therefore
\( \E(s^2\W^{-1})=\E(s^2)\E(\W^{-1})=\sigma^2\bSigma_{XX}^{-1}/(n-k-2) \) by @lem-cor-centred-wishart and @lem-cor-inverse-wishart.
:::

### C. Going deeper

::: {#exr-cor-omitted-prediction}
[C1]

Let the rows \( (\mathbf{X}_i\T,Y_i) \) be jointly normal with \( k \) regressors, and fit only the first \( k_1<k \) of them.
Show that the expected error on a new case is \( \sigma_1^2(1+1/n)(n-2)/(n-k_1-2) \), where \( \sigma_1^2 \) is the
error variance of the best linear predictor from the first \( k_1 \) regressors. *Hint:* under joint
normality the smaller model is itself a correct normal model with random regressors. For
\( n=30 \), \( k=10 \), \( k_1=4 \) and \( \sigma_1^2=1.2\sigma^2 \), is the smaller model better for prediction? What is
the largest \( k_1 \) for which it is?
:::

::: {.solution}
The pairs \( (\mathbf{X}_{1i},Y_i) \), with \( \mathbf{X}_{1i} \) the first \( k_1 \) regressors, are jointly normal, so given
\( \mathbf{X}_1 \) the response follows a normal linear model with \( k_1 \) regressors and error variance \( \sigma_1^2 \)
(@eq-mvn-random-x-model). The fitted predictor uses only these pairs, so @thm-cor-prediction-error applies with
\( k_1 \) and \( \sigma_1^2 \) in place of \( k \) and \( \sigma^2 \). With the numbers given, the full model has expected
error \( 1.607\sigma^2 \) and the smaller one \( 1.2\cdot(31/30)\cdot28/24\,\sigma^2=1.447\sigma^2 \), so the smaller model
predicts better despite its larger population error. Since the error of the smaller model increases with
\( k_1 \), it wins as long as \( 1.2\cdot28/(28-k_1)<28/18 \), that is \( k_1<6.4 \): at \( k_1=6 \) its error is
\( 1.578\sigma^2 \), at \( k_1=7 \) it is \( 1.653\sigma^2 \). The largest such \( k_1 \) is \( 6 \).
:::

