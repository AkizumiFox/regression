# Mallows' Cp, AIC and BIC

A **selection criterion** assigns a number to each candidate model and picks the smallest. \( C_p \) estimates
the prediction error of [Section 29.1](01-prediction-error.html), AIC a likelihood-based discrepancy, and BIC the
posterior probability of the model. The first two nearly always agree; the third answers a different question, and
the difference shows up in large samples.

## Candidates

The candidates are a finite collection of models \( \C(\X_k) \), \( k\in\mathcal K \), with
\( \rank(\X_k)=p_k \), each fitted by least squares: \( \hat{\bmu}_k=\M_k\Y \), with residual sum of squares
\( \text{SSE}_k=\norm{(\I-\M_k)\Y}^2 \). Typically \( \X_k \) consists of the intercept and a subset of the columns
of a full model matrix \( \X_F \) of rank \( P<n \), and \( s^2=\text{SSE}_F/(n-P) \) is the usual variance estimate
from the full model. The data follow @eq-sel-model. By @thm-sel-optimism(c) the scaled risk of
candidate \( k \) is
\[
J_k=\frac{\E\norm{\hat{\bmu}_k-\bmu}^2}{\sigma^2}=p_k+\gamma_k,\qquad
\gamma_k=\frac{\norm{(\I-\M_k)\bmu}^2}{\sigma^2},
\]{#eq-sel-scaled-risk}

variance plus squared bias. The bias term \( \gamma_k \) is the noncentrality @eq-dep-gamma of the test of
candidate \( k \) against any correct larger model.

## Mallows' Cp

::: {#def-sel-cp}
[Mallows' \( C_p \)]

The \( C_p \) statistic of candidate \( k \) is
\[
C_k=\frac{\text{SSE}_k}{s^2}-n+2p_k .
\]
:::

With \( \sigma^2 \) in place of \( s^2 \), \( \sigma^2(C_k+n)/n \) is the covariance-penalty estimate
@eq-sel-covariance-penalty of the in-sample prediction error of candidate \( k \). The scaling makes \( C_k \)
comparable with \( p_k \).

::: {#prp-sel-cp}
[Properties of \( C_p \)]

::: {.enumerate options="label=(\alph*)"}
1. Under @eq-sel-model, \( C_k^{\sigma}=\text{SSE}_k/\sigma^2-n+2p_k \) is unbiased for \( J_k \).

2. Assume normal errors, \( \bmu\in\C(\X_F) \), \( \C(\X_k)\subseteq\C(\X_F) \) and \( n-P>2 \). Then
   \[
   \E(C_k)=J_k+\frac{2(P-p_k+\gamma_k)}{n-P-2}.
   \]

3. Let \( F_k=[(\text{SSE}_k-\text{SSE}_F)/(P-p_k)]/s^2 \) be the \( F \) statistic for testing candidate \( k \)
   against the full model (\( p_k<P \)). Then
   \[
   C_k=(P-p_k)(F_k-1)+p_k .
   \]
   In particular \( C_F=P \), and \( C_k\le p_k \) iff \( F_k\le1 \).
:::

:::

::: {.proof}
(a) By @thm-sel-optimism(c), \( \E(\text{SSE}_k)=\sigma^2(n-p_k)+\norm{(\I-\M_k)\bmu}^2=\sigma^2(n-p_k+\gamma_k) \).
(b) Write \( \text{SSE}_k=\text{SSE}_F+(\text{SSE}_k-\text{SSE}_F) \). By @thm-glh-f-test(a), applied to \( \C(\X_k)\subseteq\C(\X_F) \),
the two terms are independent, \( \text{SSE}_F/\sigma^2\sim\chi^2(n-P) \), and
\( (\text{SSE}_k-\text{SSE}_F)/\sigma^2 \) is noncentral \( \chi^2 \) with \( P-p_k \) degrees of freedom and noncentrality
\( \norm{(\M_F-\M_k)\bmu}^2/\sigma^2=\gamma_k \), so its mean is \( P-p_k+\gamma_k \). By @lem-cor-inverse-chisq,
\( \E(\sigma^2/s^2)=(n-P)/(n-P-2) \). Hence
\[
\begin{aligned}
\E\Bigl(\frac{\text{SSE}_k}{s^2}\Bigr)&=(n-P)+(P-p_k+\gamma_k)\,\frac{n-P}{n-P-2}\\
&=n-p_k+\gamma_k+\frac{2(P-p_k+\gamma_k)}{n-P-2},
\end{aligned}
\]
and subtracting \( n-2p_k \) gives the result.
(c) By the definition of \( F_k \), \( \text{SSE}_k/s^2=(n-P)+(P-p_k)F_k \). Subtract \( n-2p_k \).
:::

By (b), estimating \( \sigma^2 \) costs little unless the full model uses up most of the data. Mallows (1973) proposed
plotting \( C_k \) against \( p_k \): points near the line \( C_p=p \) have small bias, and the lowest point estimates the best
predictor. By (c), adding \( q \) columns changes \( C_p \) by \( 2q-(\text{SSE}_{\text{small}}-\text{SSE}_{\text{big}})/s^2 \), so
\( C_p \) prefers the larger model iff the \( F \) statistic of the added columns, with the full-model \( s^2 \), exceeds \( 2 \).

## AIC

If \( g \) is the true density and \( f \) a model density, the **Kullback–Leibler divergence** \( \int g\log(g/f)\ge0 \) is zero iff
\( f=g \) (Kullback and Leibler 1951), and up to a term free of \( f \) it is \( \E_g[-\log f(\Y^*)] \) for a new data set
\( \Y^*\sim g \). A fitted model \( f(\cdot\,;\hat{\boldsymbol{\uptheta}}) \) is therefore judged by the **discrepancy**
\[
\Delta=\E_{\Y^*}\bigl[-2\log f(\Y^*;\hat{\boldsymbol{\uptheta}}(\Y))\bigr],
\]
with \( \Y^* \) independent of \( \Y \). The maximized log-likelihood evaluates the same function at the data themselves, and
is optimistic in the same way as the training error.

::: {#def-sel-aic-bic}
[AIC and BIC]

For a model with \( d_k \) free parameters and maximized log-likelihood \( \ell_k \),
\[
\text{AIC}_k=-2\ell_k+2d_k,\qquad \text{BIC}_k=-2\ell_k+d_k\log n .
\]
In the normal linear model with unknown variance, \( d_k=p_k+1 \) and
\( -2\ell_k=n\log(2\pi\,\text{SSE}_k/n)+n \) (@thm-opt-mle), so up to a constant common to all candidates,
\[
\text{AIC}_k=n\log\frac{\text{SSE}_k}n+2(p_k+1),\qquad
\text{BIC}_k=n\log\frac{\text{SSE}_k}n+(p_k+1)\log n .
\]{#eq-sel-aic-normal}

:::

Akaike (1973, 1974) showed that \( 2d \) is the leading term of the optimism of \( -2\ell \) for a correct model
fitted by maximum likelihood, in large samples. In the normal linear model the optimism can be computed
exactly, without asymptotics.

::: {#thm-sel-aicc}
[The expected discrepancy of a normal linear model]

Let \( \Y\sim\Normal_n(\bmu,\sigma^2\I) \) with \( \bmu\in\C(\X_k) \) and \( n>p_k+2 \), and let
\( \hat{\sigma}^2_k=\text{SSE}_k/n \). For the fitted density \( \Normal_n(\hat{\bmu}_k,\hat{\sigma}^2_k\I) \),
\[
\E(\Delta_k)=\E\bigl[n\log(2\pi\hat{\sigma}^2_k)\bigr]+\frac{n(n+p_k)}{n-p_k-2} .
\]
Consequently
\[
\text{AICc}_k=-2\ell_k+2(p_k+1)+\frac{2(p_k+1)(p_k+2)}{n-p_k-2}
\]{#eq-sel-aicc}

is an unbiased estimate of \( \E(\Delta_k) \).
:::

::: {.proof}
For a new data set \( \Y^*=\bmu+\be^* \),
\( -2\log f(\Y^*)=n\log(2\pi\hat{\sigma}^2_k)+\norm{\Y^*-\hat{\bmu}_k}^2/\hat{\sigma}^2_k \), and as in @eq-sel-err-in
\( \E\bigl(\norm{\Y^*-\hat{\bmu}_k}^2\bigm|\Y\bigr)=n\sigma^2+\norm{\hat{\bmu}_k-\bmu}^2 \). So
\[
\Delta_k=n\log(2\pi\hat{\sigma}^2_k)+\frac{n\sigma^2+\norm{\M_k\be}^2}{\hat{\sigma}^2_k},
\]
using \( \hat{\bmu}_k-\bmu=\M_k\be \) when \( \bmu\in\C(\X_k) \). Now \( \norm{\M_k\be}^2/\sigma^2\sim\chi^2(p_k) \) and
\( n\hat{\sigma}_k^2/\sigma^2=\norm{(\I-\M_k)\be}^2/\sigma^2\sim\chi^2(n-p_k) \) are independent (@thm-opt-sampling), and
\( \E(\sigma^2/\hat{\sigma}^2_k)=n/(n-p_k-2) \) by @lem-cor-inverse-chisq. Hence the second term has mean
\( n(n+p_k)/(n-p_k-2) \). For the second statement, \( -2\ell_k=n\log(2\pi\hat{\sigma}^2_k)+n \), and the identity
\( n(n+p)=(n-p-2)\bigl(n+2(p+1)\bigr)+2(p+1)(p+2) \) shows that
\( n(n+p)/(n-p-2)=n+2(p+1)+2(p+1)(p+2)/(n-p-2) \).
:::

The corrected criterion AICc is due to Sugiura (1978) and Hurvich and Tsai (1989). Its extra term vanishes as
\( n\to\infty \) with \( p_k \) fixed, which recovers AIC. It matters when \( p_k \) is not small compared with \( n \).
With \( n=20 \) and \( p_k=5 \), AIC charges \( 2(p_k+1)=12 \), while the exact charge in the proof is
\( 2(p_k+1)\,n/(n-p_k-2)=18.46 \). A simulation of \( 100000 \) data sets gives a mean discrepancy of \( 68.11 \), equal to
the theorem's value to two decimals. AIC underestimates it by \( 6.46 \), and the underestimate favours the larger
models.

```{.python .run #cell-aic-bic-aicc}
import numpy as np
rng = np.random.default_rng(2902)
n_c, p_c, reps_c = 20, 5, 100_000
Xc = np.column_stack([np.ones(n_c), rng.normal(size=(n_c, p_c - 1))])
mu_c = Xc @ np.linspace(1, 2, p_c)                         # the model is correct
H = Xc @ np.linalg.solve(Xc.T @ Xc, Xc.T)
Yc = mu_c + rng.normal(size=(reps_c, n_c))
fit_c = Yc @ H
s2ml = np.mean((Yc - fit_c) ** 2, axis=1)                 # ML variance estimate
# discrepancy: expected -2 log-likelihood of an independent copy Y*, at the fitted values
disc = n_c * np.log(2 * np.pi * s2ml) + (n_c + np.sum((fit_c - mu_c) ** 2, axis=1)) / s2ml
target = np.mean(n_c * np.log(2 * np.pi * s2ml)) + n_c * (n_c + p_c) / (n_c - p_c - 2)
print(f"E(discrepancy) {disc.mean():.2f},  formula {target:.2f}")
print(f"penalty in AIC {2 * (p_c + 1)},  in AICc {2 * (p_c + 1) * n_c / (n_c - p_c - 2):.2f}")
```

For a wrong candidate the penalty is no longer free of \( \bmu \), and \( 2d_k \) is an approximation, justified
when the candidate is close to the truth (Konishi and Kitagawa 2008 discuss corrections). In practice AIC and
\( C_p \) rank candidates almost identically (@exr-sel-aic-cp).

## BIC

Suppose the candidates carry prior probabilities \( \pi_k \), and the parameters \( \boldsymbol{\uptheta}_k \) of model \( k \) a
prior density \( \pi(\boldsymbol{\uptheta}_k\mid k) \). The posterior probability of model \( k \) is proportional to
\( \pi_k\,m_k(\y) \), where
\[
m_k(\y)=\int f_k(\y\mid\boldsymbol{\uptheta}_k)\,\pi(\boldsymbol{\uptheta}_k\mid k)\,d\boldsymbol{\uptheta}_k
\]
is the **marginal likelihood** of the model. For the conjugate normal prior of @thm-opt-bayes-conjugate the integral
has a closed form. BIC is an approximation that needs no particular prior.

::: {.remark}
[Sketch: where \( \log n \) comes from]

This is a heuristic derivation; Schwarz (1978) proved the result for exponential families, and Kass and
Raftery (1995) give conditions and error terms. Let \( \hat{\boldsymbol{\uptheta}} \) maximize \( \ell(\boldsymbol{\uptheta})=\log f_k(\y\mid\boldsymbol{\uptheta}) \), and
let \( \mathbf{J}_n=-\nabla^2\ell(\hat{\boldsymbol{\uptheta}}) \) be the observed information, which for independent observations
grows like \( n \): \( \mathbf{J}_n=n\tilde{\mathbf{J}} \) with \( \tilde{\mathbf{J}} \) of order one. Expanding \( \ell \) to second order
about \( \hat{\boldsymbol{\uptheta}} \) (Laplace's method, Tierney and Kadane 1986),
\[
\begin{aligned}
m_k(\y)&\approx f_k(\y\mid\hat{\boldsymbol{\uptheta}})\,\pi(\hat{\boldsymbol{\uptheta}}\mid k)\int
e^{-(\boldsymbol{\uptheta}-\hat{\boldsymbol{\uptheta}})\T\mathbf{J}_n(\boldsymbol{\uptheta}-\hat{\boldsymbol{\uptheta}})/2}\,d\boldsymbol{\uptheta}\\
&=f_k(\y\mid\hat{\boldsymbol{\uptheta}})\,\pi(\hat{\boldsymbol{\uptheta}}\mid k)\,(2\pi)^{d_k/2}\det(\mathbf{J}_n)^{-1/2}.
\end{aligned}
\]
Since \( \log\det\mathbf{J}_n=d_k\log n+\log\det\tilde{\mathbf{J}} \),
\[
-2\log m_k(\y)=-2\ell_k+d_k\log n+O(1),
\]
where the \( O(1) \) collects the prior density at \( \hat{\boldsymbol{\uptheta}} \), \( \det\tilde{\mathbf{J}} \) and \( 2\pi \), none of which grows
with \( n \). So \( \text{BIC}_k \) approximates \( -2\log m_k \), and with equal prior probabilities the model with the
smallest BIC is approximately the most probable.
:::

The error is of order one, so \( \exp(-\text{BIC}_k/2) \) is a crude estimate of a posterior probability, but the rate
is right: each parameter costs \( \log n \). For a "unit information" prior the error is of order \( n^{-1/2} \) (Kass
and Wasserman 1995).

## Three criteria, one data set, one threshold

Adjusted \( R^2 \) of [Chapter 9](../ch09-sums-of-squares/index.html) is sometimes used as a criterion too.

::: {#prp-sel-adjusted-link}
[Adjusted \( R^2 \) as a criterion]

Maximizing \( \bar{R}^2 \) over candidates that contain the intercept is minimizing \( s_k^2=\text{SSE}_k/(n-p_k) \), and for two
nested candidates the larger has the larger \( \bar{R}^2 \) iff its added columns have \( F>1 \).
:::

::: {.proof}
These are @prp-ss-adjusted-r2(d) and (c).
:::

The threshold \( 1 \) is lower than any of the others, so \( \bar{R}^2 \) chooses the largest models. Since
\( \E(s_k^2)=\sigma^2(1+\gamma_k/(n-p_k)) \) (@thm-dep-omitted(c)), minimizing \( s^2_k \) aims at the smallest *bias*,
not the smallest prediction error.

For two nested candidates, each criterion is a test with its own critical value. Let the larger model
add \( q \) columns, have residual degrees of freedom \( \nu \), and let
\( F=[(\text{SSE}_{\text{small}}-\text{SSE}_{\text{big}})/q]/(\text{SSE}_{\text{big}}/\nu) \), so that
\( \text{SSE}_{\text{small}}/\text{SSE}_{\text{big}}=1+qF/\nu \). The larger model is preferred by

| criterion | when | threshold for \( F \), approximately |
|---|---|---|
| adjusted \( R^2 \) | \( F>1 \) (@prp-sel-adjusted-link) | \( 1 \) |
| \( C_p \) | \( F>2 \), with \( s^2 \) from the full model | \( 2 \) |
| AIC | \( n\log(1+qF/\nu)>2q \) | \( 2 \) |
| BIC | \( n\log(1+qF/\nu)>q\log n \) | \( \log n \) |
| \( F \) test at level \( \alpha \) | \( F>F_\alpha(q,\nu) \) | \( F_\alpha(q,\nu) \) |

The approximations use \( \log(1+x)\approx x \) and \( \nu\approx n \). With \( n=50 \), \( q=1 \) and \( \nu=47 \), the exact
thresholds are \( 1.92 \) for AIC and \( 3.82 \) for BIC (\( \log n=3.91 \)), against \( F_{0.05}(1,47)=4.05 \). As tests
of one column, AIC acts at level \( 0.17 \) and BIC at level \( 0.056 \), close to the conventional \( 5\% \) at this
sample size. BIC's level falls as \( n \) grows, and AIC's does not.

::: {#exm-sel-state-cp}
[Violent crime in the states]

Return to the \( 50 \) states of [Chapter 6](../ch06-projections/index.html) (the District of Columbia left out) and
take as response the logarithm of the violent crime rate in 2009. The candidate regressors are the percentages
of adults with a high-school diploma, of people below the poverty line, of single-parent households, of white
residents and of urban residents, so there are \( 2^5=32 \) candidates, all with an intercept. The full model has
\( P=6 \) and \( s^2=0.1035 \). The eight candidates with the smallest \( C_p \) are

| regressors | \( p \) | \( C_p \) | \( \Delta\text{AIC} \) | \( \Delta\text{BIC} \) | \( \bar{R}^2 \) |
|---|---|---|---|---|---|
| single, urban | 3 | 2.20 | 0.00 | 0.42 | 0.471 |
| poverty, single, urban | 4 | 2.47 | 0.09 | 2.42 | 0.480 |
| single | 2 | 3.54 | 1.49 | 0.00 | 0.445 |
| high school, single, urban | 4 | 3.94 | 1.72 | 4.05 | 0.463 |
| high school, poverty, single, urban | 5 | 4.01 | 1.57 | 5.82 | 0.474 |
| single, white, urban | 4 | 4.04 | 1.82 | 4.16 | 0.462 |
| poverty, single, white, urban | 5 | 4.47 | 2.09 | 6.34 | 0.468 |
| high school, single | 3 | 5.48 | 3.42 | 3.85 | 0.434 |

where \( \Delta \) is the difference from the smallest value over all \( 32 \) candidates. \( C_p \) and AIC choose single
parenthood and urbanization; BIC drops urbanization; \( \bar{R}^2 \) adds poverty. The first three rows are within
two units of each other on every scale, and none of the criteria separates them convincingly. [Figure 29.2.1](#fig-sel-cp)
is Mallows' plot. Every candidate without single parenthood lies far above the line \( C_p=p \): it is badly biased.
:::

::: {when-format="html"}
![**Figure 29.2.1.** \( C_p \) against the number of columns for the \( 32 \) candidate models for log violent
crime in the \( 50 \) states. The dashed line is \( C_p=p \). The three labelled models are the choices of BIC,
of \( C_p \) and AIC, and of adjusted \( R^2 \).](cp_plot.svg){#fig-sel-cp width=72%}
:::

::: {when-format="pdf"}
![\( C_p \) against the number of columns for the \( 32 \) candidate models for log violent
crime in the \( 50 \) states. The dashed line is \( C_p=p \). The three labelled models are the choices of BIC,
of \( C_p \) and AIC, and of adjusted \( R^2 \).](cp_plot.pdf){width=72%}
:::

```{.python .run #cell-state-subsets-subsets}
import itertools
import numpy as np
import statsmodels.api as sm
data = sm.datasets.statecrime.load_pandas().data.drop(index="District of Columbia")
names = ["hs_grad", "poverty", "single", "white", "urban"]
y = np.log(data["violent"].to_numpy())
n = len(y)


def sse(cols):
    X = np.column_stack([np.ones(n)] + [data[c].to_numpy() for c in cols])
    b, *_ = np.linalg.lstsq(X, y, rcond=None)
    return np.sum((y - X @ b) ** 2), X.shape[1]


sse_full, P = sse(names)
s2 = sse_full / (n - P)                                  # sigma^2 estimate from the largest model
sst = np.sum((y - y.mean()) ** 2)
rows = []
for k in range(len(names) + 1):
    for cols in itertools.combinations(names, k):
        e, p = sse(cols)
        rows.append(dict(model=cols, p=p, sse=e,
                         cp=e / s2 - n + 2 * p,
                         aic=n * np.log(e / n) + 2 * (p + 1),
                         bic=n * np.log(e / n) + np.log(n) * (p + 1),
                         adjr2=1 - (e / (n - p)) / (sst / (n - 1))))
for crit in ("cp", "aic", "bic"):
    best = min(rows, key=lambda r: r[crit])
    print(f"{crit:4s} chooses {best['model']}")
print("adj R2 chooses", max(rows, key=lambda r: r["adjr2"])["model"])
```

## Consistency and efficiency

The thresholds explain the large-sample behaviour. The \( F \) statistic of an irrelevant column stays of order one, so
a growing threshold such as \( \log n \) eventually excludes it and a fixed one does not; that of a relevant column grows
like \( n \).

::: {#prp-sel-aic-bic}
[AIC and BIC in large samples]

Let the candidate set \( \mathcal K \) be finite and fixed, and consider a sequence of designs and means with
\( n\to\infty \). Let the errors be independent and identically distributed with mean zero and variance \( \sigma^2 \), and let
\( \sup_n\norm{\bmu}^2/n<\infty \). Call \( k \) *correct* if \( \bmu\in\C(\X_k) \). Assume that there is a correct \( k_0 \) such that
every other correct candidate has \( p_k>p_{k_0} \), and that \( \liminf_n\norm{(\I-\M_k)\bmu}^2/n>0 \) for every
candidate that is not correct. Select \( \hat{k} \) by minimizing \( \text{IC}_k=n\log(\text{SSE}_k/n)+c_np_k \).

::: {.enumerate options="label=(\alph*)"}
1. If \( c_n\to\infty \) and \( c_n/n\to0 \), then \( \Pr(\hat{k}=k_0)\to1 \). In particular BIC (\( c_n=\log n \)) is consistent.

2. If \( c_n=c \) is constant, the errors are normal, and \( k \) is a correct candidate with
   \( \C(\X_{k_0})\subset\C(\X_k) \) and \( q=p_k-p_{k_0} \), then
   \( \Pr(\text{IC}_k<\text{IC}_{k_0})\to\Pr\{\chi^2(q)>cq\}>0 \). In particular AIC (\( c=2 \)) is not consistent when such
   a candidate is available.

3. (*Efficiency, stated with a sketch.*) Suppose instead that no candidate is correct, and that the candidate
   sets grow with \( n \) so that the smallest risk over candidates tends to infinity. Then, under moment conditions on
   the errors and a condition limiting the number of candidates, with \( \sigma^2 \) known or consistently estimated
   in \( C_p \), and for AIC with the largest candidate dimension \( o(n) \), the loss \( \norm{\hat{\bmu}_{\hat{k}}-\bmu}^2 \) of the
   \( C_p \) or AIC choice divided by \( \min_k\norm{\hat{\bmu}_k-\bmu}^2 \) tends to one in probability (Shibata 1981; Li
   1987). A criterion with \( c_n\to\infty \) does not have this property in general.
:::

:::

::: {.proof}
Throughout, \( \text{SSE}_k=\norm{(\I-\M_k)\bmu}^2+2\bmu\T(\I-\M_k)\be+\norm{(\I-\M_k)\be}^2 \). The last term is
\( \norm{\be}^2-\norm{\M_k\be}^2 \), where \( \norm{\be}^2/n\to\sigma^2 \) in probability by the weak law of large numbers and
\( \E\norm{\M_k\be}^2=\sigma^2p_k \), so \( \norm{\M_k\be}^2=O_p(1) \) by Markov's inequality. The middle term has mean
zero and variance \( 4\sigma^2\norm{(\I-\M_k)\bmu}^2=O(n) \), so it is \( O_p(\sqrt n) \) by Chebyshev's inequality.

(a) It suffices to show \( \Pr(\text{IC}_k>\text{IC}_{k_0})\to1 \) for each \( k\ne k_0 \), since \( \mathcal K \) is finite. For
\( k_0 \), \( \text{SSE}_{k_0}=\norm{\be}^2-\norm{\M_{k_0}\be}^2 \), so \( \text{SSE}_{k_0}/n\to\sigma^2 \).

*A correct \( k\ne k_0 \).* Then \( \lvert\text{SSE}_k-\text{SSE}_{k_0}\rvert\le\norm{\M_k\be}^2+\norm{\M_{k_0}\be}^2=O_p(1) \), and
\( n\log(\text{SSE}_k/\text{SSE}_{k_0})=n\log\bigl(1+(\text{SSE}_k-\text{SSE}_{k_0})/\text{SSE}_{k_0}\bigr)=O_p(1) \), because the ratio inside
the logarithm is \( O_p(1/n) \). Hence \( \text{IC}_k-\text{IC}_{k_0}=O_p(1)+c_n(p_k-p_{k_0}) \), which tends to \( +\infty \) in
probability since \( p_k>p_{k_0} \) and \( c_n\to\infty \).

*A candidate that is not correct.* Let \( 0<2b<\liminf_n\norm{(\I-\M_k)\bmu}^2/n \). By the expansion,
\( \text{SSE}_k/n\ge\sigma^2+b \) with probability tending to one, so \( \log(\text{SSE}_k/\text{SSE}_{k_0})\ge\delta \) with probability
tending to one, where \( \delta=\tfrac12\log(1+b/\sigma^2)>0 \). Then
\( \text{IC}_k-\text{IC}_{k_0}\ge n\delta-c_n\lvert p_k-p_{k_0}\rvert \), which tends to \( +\infty \) because \( c_n/n\to0 \).

(b) Here \( D=\text{SSE}_{k_0}-\text{SSE}_k=\norm{(\M_k-\M_{k_0})\be}^2 \), and \( D/\sigma^2\sim\chi^2(q) \) exactly for every \( n \) (@thm-qf-chisq). With \( \text{SSE}_k/n\to\sigma^2 \), Slutsky's theorem gives \( nD/\text{SSE}_k\to\chi^2(q) \) in distribution,
and \( n\log(1+D/\text{SSE}_k)-nD/\text{SSE}_k\to0 \) in probability because \( D/\text{SSE}_k=O_p(1/n) \) and
\( \lvert\log(1+x)-x\rvert\le x^2 \) for \( x\ge0 \). So \( n\log(\text{SSE}_{k_0}/\text{SSE}_k)\to\chi^2(q) \) in distribution. Now
\( \text{IC}_k<\text{IC}_{k_0} \) iff \( n\log(\text{SSE}_{k_0}/\text{SSE}_k)>cq \), and \( cq \) is a continuity point of the limit.

(c) *Sketch.* By @prp-sel-cp(a), \( C_k^\sigma \) is unbiased for the scaled risk of each candidate, and
\( C_k^\sigma-J_k \) and the loss minus the risk are sums of centred quadratic and linear forms in \( \be \). The
conditions make these fluctuations small compared with the risks, uniformly over the candidates, so minimizing
\( C_k \) is asymptotically the same as minimizing the loss. Li (1987) gives such conditions for \( C_p \), cross-validation and
generalized cross-validation: finite moments of order \( 4m \) for the errors and
\( \sum_{k}(n R_k)^{-m}\to0 \) for some \( m \), where \( R_k \) is the risk of candidate \( k \) per observation. AIC differs from
\( C_p \) by terms of order \( p_k^2/n \), which are negligible when the candidate dimensions are \( o(n) \) (Shibata 1981). A growing penalty \( c_n \) biases the choice towards
candidates whose risk is larger by a factor that does not tend to one, as [Figure 29.2.2](#fig-sel-aic-bic)(b) shows.
:::

For AIC the limits in (b) are \( 0.157 \), \( 0.135 \) and \( 0.075 \) for \( q=1,2,5 \): one extra column is kept about one time
in six however large the sample. Nishii (1984) proved (a) and (b) for the normal linear model.

::: {#exm-sel-aic-bic}
[Consistency against efficiency]

Panel (a) of [Figure 29.2.2](#fig-sel-aic-bic) uses six independent normal regressors, three of them with slopes
\( 1 \), \( 0.5 \) and \( 0.25 \) and three with slope zero, and all \( 64 \) subsets as candidates. With \( n=50 \) the true subset is
found by AIC with probability \( 0.28 \) and by BIC with probability \( 0.24 \), because the smallest slope is hard to
detect. With \( n=3200 \) BIC finds it with probability \( 0.99 \), while AIC stays at \( 0.60 \). That is the limit
\( \Pr\{\chi^2(1)\le2\}^3=0.598 \), because in large samples AIC keeps each irrelevant regressor independently when its squared
\( t \) statistic exceeds about \( 2 \).

In panel (b) the truth has \( 1000 \) regressors with slopes \( 2/j \), all nonzero, and the candidates are the nested
models containing the first \( k \) regressors, \( k\le n/4 \). No candidate is correct. The ratio of the loss of the chosen
model to the smallest loss among the candidates falls for AIC, from \( 1.18 \) at \( n=100 \) to \( 1.07 \) at \( n=3200 \),
while for BIC it stays between \( 1.3 \) and \( 1.8 \) at every sample size.
:::

::: {when-format="html"}
![**Figure 29.2.2.** (a) Probability of selecting the true subset among \( 64 \), when a finite true model exists;
dashed: AIC's limit \( 0.598 \). (b) Loss of the selected model relative to the best candidate, when no
candidate is true. \( 1000 \) simulated data sets per point in (a), \( 400 \) in (b).](aic_bic.svg){#fig-sel-aic-bic width=100%}
:::

::: {when-format="pdf"}
![(a) Probability of selecting the true subset among \( 64 \), when a finite true model exists;
dashed: AIC's limit \( 0.598 \). (b) Loss of the selected model relative to the best candidate, when no
candidate is true. \( 1000 \) simulated data sets per point in (a), \( 400 \) in (b).](aic_bic.pdf){width=100%}
:::

```{.python .run #cell-aic-bic-consistency}
import itertools
import numpy as np
rng = np.random.default_rng(2903)
beta = np.array([1.0, 0.5, 0.25, 0.0, 0.0, 0.0])        # three relevant, three irrelevant
true_set = (0, 1, 2)
subsets = [s for k in range(7) for s in itertools.combinations(range(6), k)]


def choose(X, y, penalty):
    """All-subsets choice by n log(SSE/n) + penalty * (number of columns)."""
    n = len(y)
    best, arg = np.inf, None
    for s in subsets:
        Xs = np.column_stack([np.ones(n), X[:, list(s)]])
        b, *_ = np.linalg.lstsq(Xs, y, rcond=None)
        crit = n * np.log(np.sum((y - Xs @ b) ** 2) / n) + penalty * Xs.shape[1]
        if crit < best:
            best, arg = crit, s
    return arg


def hit_rates(n, reps):
    X = rng.normal(size=(n, 6))                           # the design, fixed for this n
    hits = np.zeros(2)
    for _ in range(reps):
        y = 1 + X @ beta + rng.normal(size=n)
        hits += [choose(X, y, 2.0) == true_set, choose(X, y, np.log(n)) == true_set]
    return hits / reps


for n in (50, 200, 800):
    aic, bic = hit_rates(n, 100)
    print(f"n = {n:4d}: AIC finds the true subset {aic:.2f} of the time, BIC {bic:.2f}")
```

The two properties cannot be combined: Yang (2005) showed that no rule can be consistent and also attain the best
rate of prediction risk. The choice between AIC and BIC is a choice of goal: prediction (AIC, \( C_p \),
cross-validation) or identifying a small true set of regressors (BIC).

## Selection in an orthonormal design

When the columns are orthonormal the search disappears, and the criteria can be compared exactly.

::: {#prp-sel-orthonormal}
[Selection as hard thresholding]

Let \( \X \) have orthonormal columns, \( \X\T\X=\I_m \), put \( \Z=\X\T\Y \), and let \( \sigma^2 \) be known. Among all subsets
\( S \) of the columns, \( \text{SSE}_S/\sigma^2+c^2\lvert S\rvert \) is minimized (uniquely, with probability one) by
\( S=\{j:\lvert Z_j\rvert>c\sigma\} \). The selected fit has coefficients \( Z_j1\{\lvert Z_j\rvert>c\sigma\} \).
:::

::: {.proof}
The projection onto the columns in \( S \) is \( \sum_{j\in S}\x_j\x_j\T \), so
\( \text{SSE}_S=\norm{\Y}^2-\sum_{j\in S}Z_j^2 \) and the criterion is
\( \norm{\Y}^2/\sigma^2+\sum_{j\in S}(c^2-Z_j^2/\sigma^2) \). A column lowers it iff \( Z_j^2>c^2\sigma^2 \), and ties have
probability zero.
:::

With known \( \sigma^2 \), \( C_p \) and AIC take \( c=\sqrt2 \) and BIC takes \( c=\sqrt{\log n} \). Testing each coefficient
at level \( 5\% \) and keeping the significant ones (the **pretest estimator** of Bancroft 1944) takes
\( c=1.96 \). This is the hard thresholding of @thm-shr-lasso-orthonormal, and each coordinate can be studied alone.
With \( \sigma=1 \) and \( Z\sim\Normal(\theta,1) \), the risk of \( Z1\{\lvert Z\rvert>c\} \) is
\[
\begin{aligned}
R(\theta;c)&=\theta^2D+1-D+(c-\theta)\phi(c-\theta)+(c+\theta)\phi(c+\theta),\\
D&=\Phi(c-\theta)-\Phi(-c-\theta),
\end{aligned}
\]{#eq-sel-hard-risk}

where \( D \) is the probability that the coordinate is dropped (@exr-sel-hard-risk). Least squares has risk \( 1 \), and
dropping the coordinate always has risk \( \theta^2 \).

[Figure 29.2.3](#fig-sel-pretest-risk) shows the trade. At \( \theta=0 \) the risks are \( 0.572 \) for \( c=\sqrt2 \), \( 0.279 \)
for the pretest and \( 0.075 \) for BIC with \( n=1000 \). (At zero the risk equals the degrees of freedom of
@prp-sel-hard-df: both are \( \E[Z^21\{\lvert Z\rvert>c\}] \).) The curves cross \( 1 \) at \( \theta=0.78 \), \( 0.84 \) and \( 0.92 \)
and peak at \( 1.651 \), \( 2.464 \) and \( 3.956 \), at \( \theta=1.88 \), \( 2.16 \) and \( 2.59 \). A larger threshold gains more
at zero and loses more at moderate \( \theta \), and no threshold, least squares included, is uniformly best. This is
@thm-dep-mse again: dropping a coordinate helps iff \( \theta^2<1 \), and the data cannot tell which side of the boundary
they are on. Shrinkage avoids the bump: the James–Stein estimator of @thm-shr-james-stein beats least squares
everywhere when \( m\ge3 \).

::: {when-format="html"}
![**Figure 29.2.3.** Risk of hard thresholding \( Z1\{\lvert Z\rvert>c\} \), \( Z\sim\Normal(\theta,1) \), for the thresholds
of \( C_p \) and AIC, a \( 5\% \) pretest and BIC with \( n=1000 \). Least squares has risk \( 1 \) (dashed); always
dropping the coordinate has risk \( \theta^2 \) (grey).](pretest_risk.svg){#fig-sel-pretest-risk width=72%}
:::

::: {when-format="pdf"}
![Risk of hard thresholding \( Z1\{\lvert Z\rvert>c\} \), \( Z\sim\Normal(\theta,1) \), for the thresholds
of \( C_p \) and AIC, a \( 5\% \) pretest and BIC with \( n=1000 \). Least squares has risk \( 1 \) (dashed); always
dropping the coordinate has risk \( \theta^2 \) (grey).](pretest_risk.pdf){width=72%}
:::

```{.python .run #cell-thresholding-risk-risk}
import numpy as np
from scipy.stats import norm


def hard_risk(theta, c):
    """E(Z 1{|Z| > c} - theta)^2 for Z ~ N(theta, 1)."""
    D = norm.cdf(c - theta) - norm.cdf(-c - theta)        # Pr(|Z| <= c): coordinate dropped
    return theta**2 * D + 1 - D + (c - theta) * norm.pdf(c - theta) + (c + theta) * norm.pdf(c + theta)


thetas = np.linspace(0, 6, 601)
for label, c in [("C_p / AIC", np.sqrt(2)), ("5% pretest", 1.96), ("BIC, n = 1000", np.sqrt(np.log(1000)))]:
    r = hard_risk(thetas, c)
    print(f"{label:14s} c = {c:.3f}: risk at 0 = {r[0]:.3f}, "
          f"maximum {r.max():.3f} at theta = {thetas[r.argmax()]:.2f}")
```

## Exercises

### A. Check your understanding

::: {#exr-sel-cp-hand}
[A1]

In @exm-sel-state-cp the model with single parenthood and urbanization has \( \text{SSE}=4.782 \). Compute its
\( C_p \) from \( s^2 \), \( n \) and \( p \), and the \( F \) statistic for testing it against the full model.
:::

::: {.solution}
\( C_p=4.782/0.1035-50+6=46.20-44=2.20 \). By @prp-sel-cp(c), \( 2.20=3(F-1)+3 \), so \( F=0.73 \): the three
regressors left out together have an \( F \) statistic below one.
:::

### B. Practice

::: {#exr-sel-aic-cp}
[B1]

Let \( \hat{\sigma}^2 \) be any fixed positive number and write \( \text{SSE}_k/n=\hat{\sigma}^2(1+\delta_k) \). Show that
\[
\text{AIC}_k-\text{AIC}_j=\frac{\text{SSE}_k-\text{SSE}_j}{\hat{\sigma}^2}+2(p_k-p_j)+n\,O(\delta_k^2+\delta_j^2),
\]
and conclude that with \( \hat{\sigma}^2=s^2 \) AIC and \( C_p \) differences agree when every candidate's residual mean
square is close to \( s^2 \).
:::

::: {#exr-sel-hard-risk}
[B2]

Derive @eq-sel-hard-risk. (Split the risk according to whether the coordinate is kept, and use
\( \int_a^bw^2\phi(w)\,dw=\Phi(b)-\Phi(a)-b\phi(b)+a\phi(a) \).)
:::

::: {.solution}
Let \( Z=\theta+W \). When the coordinate is dropped (\( \lvert Z\rvert\le c \), probability \( D \)) the loss is \( \theta^2 \);
when kept it is \( W^2 \). So \( R=\theta^2D+\E[W^21\{\lvert\theta+W\rvert>c\}] \). The event \( \lvert\theta+W\rvert\le c \) is
\( -c-\theta\le W\le c-\theta \), so \( \E[W^21\{\lvert\theta+W\rvert>c\}]=1-\int_{-c-\theta}^{c-\theta}w^2\phi(w)\,dw
=1-D+(c-\theta)\phi(c-\theta)+(c+\theta)\phi(c+\theta) \), using \( \phi(-x)=\phi(x) \).
:::

::: {#exr-sel-aic-threshold}
[B3]

Show that AIC prefers the larger of two nested candidates iff \( F>(\nu/q)(e^{2q/n}-1) \), and BIC iff
\( F>(\nu/q)(n^{q/n}-1) \). Find the limits of both thresholds as \( n\to\infty \) with \( q \) and \( n-\nu \) fixed.
:::

### C. Going deeper

::: {#exr-sel-aic-limit}
[C1]

In panel (a) of @exm-sel-aic-bic, suppose the three regressors with zero slopes are orthogonal to each other and to the
relevant ones. Show that for large \( n \) AIC selects the true subset iff each of them has squared
\( t \) statistic below about \( 2 \), and deduce the limit \( \Pr\{\chi^2(1)\le2\}^3 \).
:::

::: {.solution}
Every candidate that omits a relevant regressor loses with probability tending to one (@prp-sel-aic-bic, proof
of (a), second case, which holds for any \( c_n \) with \( c_n/n\to0 \)). Among candidates containing the relevant block,
orthogonality makes \( \text{SSE}_{k_0\cup A}=\text{SSE}_{k_0}-\sum_{j\in A}\norm{\M_{\{j\}}\be}^2 \), and
\( n\log(\text{SSE}_{k_0}/\text{SSE}_{k_0\cup A})=\sum_{j\in A}T_j^2+o_p(1) \) with \( T_j^2\to\chi^2(1) \) independent. AIC adds \( A \) iff this
exceeds \( 2\lvert A\rvert \), so it chooses \( A=\{j:T_j^2>2\} \) up to \( o_p(1) \), and chooses \( A=\varnothing \) with probability
tending to \( \prod_j\Pr(T_j^2\le2)=\Pr\{\chi^2(1)\le2\}^3 \).
:::

