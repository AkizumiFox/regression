# Sufficiency, completeness and minimum variance

The Gauss–Markov theorem compares least squares only with *linear* unbiased estimators, and
@exm-opt-midrange showed that the restriction matters: with uniform errors a nonlinear estimator
wins easily. This section shows that under normality the restriction can be dropped. In the normal
linear model, \( \blambda\T\hbeta \) has the smallest variance among *all* unbiased estimators of an
estimable \( \blambda\T\bbeta \), and \( s^2 \) has the smallest variance among all unbiased estimators of
\( \sigma^2 \). The argument is the classical one through sufficiency and completeness. It also produces
best unbiased estimators of nonlinear functions, such as \( \sigma \) itself. At the end we compare the
variances with the Cramér–Rao lower bound, which least squares attains for \( \bbeta \) and \( s^2 \) does
not attain for \( \sigma^2 \).

Two facts from mathematical statistics are quoted without proof: the factorization criterion for
sufficiency and the completeness of full-rank exponential families. Both are proved in standard texts
such as Lehmann and Casella (1998) and Lehmann and Romano (2005). Everything built on them is proved here.

## Sufficiency

::: {#def-opt-sufficient}
[Sufficient and complete statistics]

Let \( \Y \) have a distribution \( P_{\bm\theta} \) from a family indexed by \( \bm\theta\in\Theta \), and let
\( \bT=\bT(\Y) \) be a statistic.

::: {.enumerate options="label=(\alph*)"}
1. \( \bT \) is **sufficient** for \( \bm\theta \) if the conditional law of \( \Y \) given \( \bT \)
   is the same for every \( \bm\theta \).

2. \( \bT \) is **complete** if, for every function \( h \) with \( \E_{\bm\theta}\lvert h(\bT)\rvert<\infty \) and
   \( \E_{\bm\theta}h(\bT)=0 \) for all \( \bm\theta\in\Theta \), we have \( \Pr_{\bm\theta}\{h(\bT)=0\}=1 \) for all
   \( \bm\theta \).
:::
:::

A sufficient statistic carries all the information about \( \bm\theta \) in the data. Once \( \bT \) is known,
the rest of \( \Y \) is noise whose distribution does not involve the parameter. Completeness says that
\( \bT \) carries no redundant information either: no nontrivial function of \( \bT \) has expectation zero
for every parameter value. Equivalently, an unbiased estimator that is a function of \( \bT \) is unique.

::: {#lem-opt-factorization}
[Factorization criterion]

Suppose \( \Y \) has a density \( f_{\bm\theta}(\y) \) on \( \Real^n \). Then \( \bT(\Y) \) is sufficient iff there
are functions \( g_{\bm\theta} \) and \( h \), with \( h \) not depending on \( \bm\theta \), such that
\( f_{\bm\theta}(\y)=g_{\bm\theta}(\bT(\y))\,h(\y) \) for all \( \y \) and \( \bm\theta \).
:::

::: {#lem-opt-exponential-family}
[Completeness in exponential families]

Suppose \( \Y \) has a density of the form
\[
f_{\bm\eta}(\y)=c(\bm\eta)\,h(\y)\exp\{\bm\eta\T\bT(\y)\},\qquad \bm\eta\in\mathcal E\subseteq\Real^k,
\]
with \( \bT(\y)\in\Real^k \). If \( \mathcal E \) contains a nonempty open subset of \( \Real^k \), then \( \bT(\Y) \) is
complete and sufficient.
:::

Sufficiency follows from @lem-opt-factorization with \( g_{\bm\eta}(\bT)=c(\bm\eta)e^{\bm\eta\T\bT} \). Completeness
rests on the uniqueness of Laplace transforms: \( \E_{\bm\eta}h(\bT)=0 \) on an open set of \( \bm\eta \) forces a
certain two-sided Laplace transform to vanish on an open set, and hence \( h(\bT)=0 \) almost surely.

For the normal linear model, expand the exponent of the density:
\[
\norm{\y-\X\bbeta}^2=\y\T\y-2\bbeta\T\X\T\y+\bbeta\T\X\T\X\bbeta .
\]
So the density of @eq-opt-normal-model is
\[
f(\y)=(2\pi\sigma^2)^{-n/2}\exp\Bigl\{-\frac{\bbeta\T\X\T\X\bbeta}{2\sigma^2}\Bigr\}
\exp\Bigl\{\frac{\bbeta\T\X\T\y}{\sigma^2}-\frac{\y\T\y}{2\sigma^2}\Bigr\},
\]{#eq-opt-expfam}

and the data enter only through \( \X\T\y \) and \( \y\T\y \).

::: {#prp-opt-complete-sufficient}
[A complete sufficient statistic]

In the normal linear model @eq-opt-normal-model, let \( \Q \) be an \( n\times r \) matrix whose columns are an
orthonormal basis of \( \C(\X) \). Then \( \bT=(\Q\T\Y,\ \Y\T\Y) \) is complete and sufficient for
\( (\bbeta,\sigma^2) \). The statistics \( (\X\T\Y,\Y\T\Y) \) and \( (\M\Y,\text{SSE}) \) are one-to-one
functions of \( \bT \), so they are complete and sufficient too.
:::

::: {.proof}
Since \( \X\bbeta\in\C(\X) \), \( \X\bbeta=\Q\bgamma \) with \( \bgamma=\Q\T\X\bbeta\in\Real^r \), and
\( \bbeta\T\X\T\y=\bgamma\T\Q\T\y \). By @eq-opt-expfam the density has the form of
@lem-opt-exponential-family with \( k=r+1 \),
\[
\bT(\y)=(\Q\T\y,\ \y\T\y),\qquad \bm\eta=\bigl(\bgamma/\sigma^2,\ {-}1/(2\sigma^2)\bigr).
\]
As \( \bbeta \) ranges over \( \Real^p \), \( \bgamma=\Q\T\X\bbeta \) ranges over all of \( \Real^r \): the \( r\times p \) matrix
\( \Q\T\X \) has rank \( r \), because \( \Q\Q\T\X=\M\X=\X \) has rank \( r \). Hence
\( \mathcal E=\Real^r\times(-\infty,0) \), which is open, and @lem-opt-exponential-family applies. For the
equivalent statistics: \( \M\Y=\Q\Q\T\Y \) (@prp-proj-orthonormal-formula) and
\( \text{SSE}=\Y\T\Y-\norm{\Q\T\Y}^2 \) are functions of \( \bT \), and conversely \( \Q\T\Y=\Q\T\M\Y \) and
\( \Y\T\Y=\text{SSE}+\norm{\M\Y}^2 \). Also \( \X\T\Y=\X\T\Q(\Q\T\Y) \), and \( \Q=\X\bm K \) for some \( \bm K \), since the
columns of \( \Q \) lie in \( \C(\X) \), so \( \Q\T\Y=\bm K\T\X\T\Y \). A one-to-one function of a complete sufficient
statistic is again complete and sufficient, since both definitions refer only to the information
(the \( \sigma \)-field) the statistic generates.
:::

The reduction is dramatic. An \( n \)-dimensional data vector is replaced by \( r+1 \) numbers: the fitted
values (equivalently, the least squares estimates of the estimable functions) and the residual sum of squares.
Under normality the residual *vector*, apart from its length, carries no information about the parameters.
Its direction is uniformly distributed over the unit sphere of \( \C(\X)\perpc \) whatever \( (\bbeta,\sigma^2) \)
is (@exr-opt-residual-direction). This is the precise sense in which the residuals are "just noise", and it
is also why they are useful for checking the model. Their pattern can reveal departures from normality
precisely because, under normality, it has a known and parameter-free distribution.

## Rao–Blackwell and Lehmann–Scheffé

::: {#thm-opt-rao-blackwell}
[Rao–Blackwell]

Let \( \bT \) be sufficient for \( \bm\theta \), and let \( \delta(\Y) \) be an unbiased estimator of \( g(\bm\theta) \)
with \( \Var_{\bm\theta}\delta<\infty \) for all \( \bm\theta \). Then \( \delta^*=\E(\delta\mid\bT) \) is a statistic,
it is unbiased for \( g(\bm\theta) \), and
\[
\begin{aligned}
\Var_{\bm\theta}(\delta)&=\Var_{\bm\theta}(\delta^*)+\E_{\bm\theta}\bigl[(\delta-\delta^*)^2\bigr]\\
&\ge\Var_{\bm\theta}(\delta^*)
\end{aligned}
\]
for every \( \bm\theta \), with equality iff \( \delta=\delta^* \) with probability one.
:::

::: {.proof}
The law of \( \Y \) given \( \bT \) does not depend on \( \bm\theta \), so neither does
\( \E(\delta(\Y)\mid\bT) \). It is therefore computable from the data, which makes it a statistic. By the
law of total expectation (@prp-rv-total-covariance), \( \E\delta^*=\E\delta=g(\bm\theta) \). By the law of
total variance in the same result,
\( \Var\delta=\E\{\Var(\delta\mid\bT)\}+\Var(\delta^*) \), and
\( \E\{\Var(\delta\mid\bT)\}=\E\{\E[(\delta-\delta^*)^2\mid\bT]\}=\E(\delta-\delta^*)^2 \). This last term is zero iff
\( \delta=\delta^* \) almost surely.
:::

Conditioning on a sufficient statistic averages out the part of an estimator that responds to
parameter-free noise. It never hurts, and it helps unless the estimator was already a function of
\( \bT \).

::: {#thm-opt-lehmann-scheffe}
[Lehmann–Scheffé]

Let \( \bT \) be complete and sufficient for \( \bm\theta \). If \( \phi(\bT) \) is unbiased for \( g(\bm\theta) \) and has
finite variance, then \( \phi(\bT) \) is a **uniformly minimum variance unbiased estimator** (UMVUE) of
\( g(\bm\theta) \): \( \Var_{\bm\theta}\phi(\bT)\le\Var_{\bm\theta}\delta \) for every unbiased \( \delta \) and every
\( \bm\theta \). It is the only unbiased function of \( \bT \) and, up to sets of probability zero, the only UMVUE.
:::

::: {.proof}
If \( \phi_1(\bT) \) and \( \phi_2(\bT) \) are both unbiased, then \( \E_{\bm\theta}\{\phi_1(\bT)-\phi_2(\bT)\}=0 \) for all
\( \bm\theta \), so \( \phi_1(\bT)=\phi_2(\bT) \) almost surely by completeness. Now let \( \delta \) be any unbiased
estimator with finite variance. By @thm-opt-rao-blackwell, \( \E(\delta\mid\bT) \) is an unbiased function
of \( \bT \), so it equals \( \phi(\bT) \) almost surely, and \( \Var\phi(\bT)\le\Var\delta \), with equality iff
\( \delta=\phi(\bT) \) almost surely. (If \( \delta \) has infinite variance there is nothing to prove.)
:::

## Best unbiased estimators in the normal linear model

::: {#thm-opt-umvue}
[Minimum variance unbiased estimation]

Assume the normal linear model @eq-opt-normal-model with \( n>r \), and let \( s^2=\text{SSE}/(n-r) \).

::: {.enumerate options="label=(\alph*)"}
1. For every \( \blambda\in\C(\X\T) \), \( \blambda\T\hbeta \) is the UMVUE of \( \blambda\T\bbeta \).

2. \( s^2 \) is the UMVUE of \( \sigma^2 \).

3. More generally, every estimator of the form \( \phi(\M\Y,\text{SSE}) \) with finite variance is the UMVUE of
   its expectation. In particular, for \( k=n-r \),
   \[
b_k=\sqrt{\frac2k}\,\frac{\Gamma\bigl(\frac{k+1}2\bigr)}{\Gamma\bigl(\frac k2\bigr)},
\]
   \( s/b_k \) is the UMVUE of \( \sigma \), and \( \x_0\T\hbeta+z\,s/b_k \) is the UMVUE of
   \( \x_0\T\bbeta+z\sigma \) for every constant \( z \) and every \( \x_0\in\C(\X\T) \).
:::
:::

::: {.proof}
By @prp-opt-complete-sufficient, \( (\M\Y,\text{SSE}) \) is complete and sufficient, so by
@thm-opt-lehmann-scheffe it suffices to show that each estimator is an unbiased function of it with finite
variance. (a) Write \( \blambda=\X\T\bm\rho \). Then \( \blambda\T\hbeta=\bm\rho\T\M\Y \) is a function of \( \M\Y \), it is
unbiased by @thm-opt-gauss-markov(a), and its variance is finite. (b) \( s^2 \) is a function of SSE, it is unbiased
because \( \E(\text{SSE})=\sigma^2(n-r) \) (@thm-rv-quadform-mean), and it has finite variance. (c) The first
statement is @thm-opt-lehmann-scheffe again. For \( \sigma \), note that
\( \text{SSE}/\sigma^2=\Y\T(\I-\M)\Y/\sigma^2\sim\chi^2(k) \) by @thm-qf-chisq(a), since \( \I-\M \) is a projection of
rank \( k \) and \( (\I-\M)\X\bbeta=\bzero \). If \( W\sim\chi^2(k) \), with density
\( w^{k/2-1}e^{-w/2}/\{2^{k/2}\Gamma(k/2)\} \), then
\[
\begin{aligned}
\E\sqrt W&=\frac{1}{2^{k/2}\Gamma(k/2)}\int_0^\infty w^{(k+1)/2-1}e^{-w/2}\,dw\\
&=\frac{2^{(k+1)/2}\Gamma\bigl(\frac{k+1}2\bigr)}{2^{k/2}\Gamma\bigl(\frac k2\bigr)}
=\sqrt2\,\frac{\Gamma\bigl(\frac{k+1}2\bigr)}{\Gamma\bigl(\frac k2\bigr)} .
\end{aligned}
\]
Hence \( \E s=\sigma\,\E\sqrt{W/k}=b_k\sigma \), and \( s/b_k \) is unbiased for \( \sigma \), with variance at most
\( \E s^2/b_k^2<\infty \). The last estimator is a sum of two unbiased functions of \( (\M\Y,\text{SSE}) \).
:::

Part (a) is the promised strengthening of Gauss–Markov. Under normality no unbiased estimator of
\( \blambda\T\bbeta \), linear or not, beats least squares. The midrange of @exm-opt-midrange cannot help
when the errors are normal. Part (c) shows the method at work beyond linear estimators. The obvious
estimate \( s \) of \( \sigma \) is biased, because \( \E s<\sqrt{\E s^2} \) by Jensen's inequality, and dividing by
\( b_k \) removes the bias.

::: {#exm-opt-umvue-sigma}
[Estimating \( \sigma \) without bias]

The constant \( b_k=\E(s)/\sigma \) depends only on the residual degrees of freedom. For the stack loss
regression, \( k=17 \) and \( b_{17}=0.9854 \), so the UMVUE of \( \sigma \) is
\( s/b_{17}=1.0148\,s \). The correction is \( 1.5\% \) here. It is large for very small
\( k \), with \( b_2=0.8862 \), and negligible for large \( k \), with
\( b_{200}=0.9988 \).
:::

```{.python .run #cell-umvue-sigma}
import numpy as np
from scipy import special

def mean_s_factor(k):
    """b_k = E(s)/sigma when k s^2 / sigma^2 ~ chi^2(k)."""
    return np.sqrt(2 / k) * np.exp(special.gammaln((k + 1) / 2) - special.gammaln(k / 2))

for k in [2, 5, 10, 17, 50, 200]:
    b = mean_s_factor(k)
    print(f"k = {k:3d}   E(s)/sigma = {b:.4f}   UMVUE of sigma = s / {b:.4f}")
```

## Rao–Blackwellization by reflection

The Rao–Blackwell theorem is not only a device for proofs. In the normal linear model we can carry out the
conditioning explicitly, and it takes every linear unbiased estimator to least squares.

::: {#exm-opt-reflection}
[Conditioning on the sufficient statistic]

Let \( \bm a\T\Y \) be any LUE of an estimable \( \blambda\T\bbeta \), such as the end-point slope of @exm-opt-three-slopes.
Let \( \R=2\M-\I \), the reflection that fixes \( \C(\X) \) and reverses
\( \C(\X)\perpc \). It is symmetric and orthogonal (\( \R^2=\I \)) and \( \R\X=\X \). Consequently:

- \( \R\Y\sim\Normal_n(\R\X\bbeta,\sigma^2\R\R\T)=\Normal_n(\X\bbeta,\sigma^2\I) \) by @thm-mvn-linear, so \( \R\Y \) has the
  same distribution as \( \Y \);
- the sufficient statistic is unchanged: \( \X\T\R\Y=\X\T\Y \) and \( \norm{\R\Y}=\norm{\Y} \).

So the pairs \( (\bm a\T\Y,\bT(\Y)) \) and \( (\bm a\T\R\Y,\bT(\Y)) \) have the same joint distribution, where
\( \bT=(\X\T\Y,\Y\T\Y) \). Hence \( \E(\bm a\T\Y\mid\bT)=\E(\bm a\T\R\Y\mid\bT) \), and averaging,
\[
\E(\bm a\T\Y\mid\bT)=\E\Bigl(\bm a\T\frac{\Y+\R\Y}2\,\Big|\,\bT\Bigr)=\E(\bm a\T\M\Y\mid\bT)=\bm a\T\M\Y=\blambda\T\hbeta ,
\]
because \( \bm a\T\M\Y \) is a function of \( \bT \) and \( \M\bm a=\bm a_* \) (@prp-opt-lue-set). The Rao–Blackwell
improvement of *every* LUE is the least squares estimate. The listing checks the identity
\( \tfrac12(\bm a\T\y+\bm a\T\R\y)=\blambda\T\hbeta \) for simulated data and the end-point slope.
:::

```{.python .run #cell-umvue-reflection}
import numpy as np

x = np.array([1.0, 2.0, 2.5, 3.0, 4.0, 5.0, 5.5, 7.0, 8.0, 9.0, 10.5, 12.0])
n = len(x)
X = np.column_stack([np.ones(n), x])
M = X @ np.linalg.solve(X.T @ X, X.T)
R = 2 * M - np.eye(n)                            # reflection: keeps C(X), flips C(X)-perp

rng = np.random.default_rng(71)
y = X @ np.array([2.0, 0.5]) + rng.normal(size=n)
y_ref = R @ y
print("X^T y unchanged:", np.allclose(X.T @ y_ref, X.T @ y))
print("y^T y unchanged:", np.isclose(y_ref @ y_ref, y @ y))

a_end = np.zeros(n)
a_end[[0, -1]] = [-1.0, 1.0]
a_end /= x[-1] - x[0]                            # the end-point slope estimator
a_ls = X @ np.linalg.solve(X.T @ X, [0.0, 1.0])
print("average of a^T y over the two points:", (a_end @ y + a_end @ y_ref) / 2)
print("least squares slope:                 ", a_ls @ y)
```

The reflection argument uses the normal distribution only through its invariance under orthogonal maps
that fix the mean. It makes concrete what "the residuals are noise" means. Reversing the residual vector
changes neither the likelihood nor the sufficient statistic, so an estimator should not depend on the
residuals' direction. An LUE that uses the residuals is replaced by the average of its value at \( \y \) and at the
reflected data, and that average no longer involves the residuals.

## The Cramér–Rao bound

A different route to optimality bounds the variance of every unbiased estimator from below and checks
whether the bound is attained. We state the bound for the full-rank normal linear model, where every
ingredient can be computed.

Write \( \bm\theta=(\bbeta\T,v)\T \) with \( v=\sigma^2 \). The **score** is the gradient of the
log-likelihood @eq-opt-loglik, evaluated at the random \( \Y \):
\[
\bm S=\begin{pmatrix}\partial\ell/\partial\bbeta\\ \partial\ell/\partial v\end{pmatrix}
=\begin{pmatrix}\X\T\be/v\\[2pt] -\dfrac{n}{2v}+\dfrac{\be\T\be}{2v^2}\end{pmatrix},
\qquad \be=\Y-\X\bbeta .
\]

::: {#prp-opt-information}
[Score and information]

In @eq-opt-normal-model with \( \rank(\X)=p \), \( \E\bm S=\bzero \) and
\[
\bm{\mathcal I}(\bm\theta)=\Cov(\bm S)=\begin{pmatrix}\X\T\X/\sigma^2&\bzero\\ \bzero\T&n/(2\sigma^4)\end{pmatrix}.
\]
:::

::: {.proof}
\( \be\sim\Normal_n(\bzero,v\I) \). Hence \( \E(\X\T\be/v)=\bzero \) and \( \Cov(\X\T\be/v)=\X\T\X/v \) by @thm-rv-linear.
Next, \( \be\T\be/v\sim\chi^2(n) \) by @thm-qf-chisq(a) with \( \bP=\I \), so it has mean \( n \) and variance \( 2n \). So the
last coordinate has mean \( -n/(2v)+nv/(2v^2)=0 \) and variance \( 2n/(4v^2)=n/(2v^2) \). Finally, by
@thm-qf-mean-var(c) with \( \bm K=\X\T \), \( \A=\I \) and \( \bmu=\bzero \), \( \Cov(\X\T\be,\be\T\be)=\bzero \).
:::

::: {#thm-opt-cramer-rao}
[Cramér–Rao bound in the normal linear model]

Assume @eq-opt-normal-model with \( \rank(\X)=p \). Let \( \delta(\Y) \) have finite variance and satisfy
\( \E_{\bm\theta}\delta=g(\bm\theta) \) for all \( \bm\theta \), with \( g \) differentiable. Then
\[
\begin{aligned}
\Var_{\bm\theta}(\delta)&\ge\nabla g(\bm\theta)\T\bm{\mathcal I}(\bm\theta)^{-1}\nabla g(\bm\theta)\\
&=\sigma^2\,\bm g_{\bbeta}\T(\X\T\X)^{-1}\bm g_{\bbeta}+\frac{2\sigma^4}{n}\,g_v^2,
\end{aligned}
\]
where \( \bm g_{\bbeta}=\partial g/\partial\bbeta \) and \( g_v=\partial g/\partial v \).
:::

::: {.proof}
Write \( f_{\bm\theta} \) for the density. Differentiating \( g(\bm\theta)=\int\delta(\y)f_{\bm\theta}(\y)\,d\y \) under the
integral sign, which is justified because the normal family is an exponential family
(Lehmann and Casella 1998, chapter 1), gives
\[
\begin{aligned}
\nabla g(\bm\theta)&=\int\delta(\y)\,\nabla\log f_{\bm\theta}(\y)\,f_{\bm\theta}(\y)\,d\y\\
&=\E(\delta\bm S)=\Cov(\delta,\bm S),
\end{aligned}
\]
the last step because \( \E\bm S=\bzero \). For any constant vector \( \bm c \),
\( 0\le\Var(\delta-\bm c\T\bm S)=\Var\delta-2\bm c\T\nabla g+\bm c\T\bm{\mathcal I}\bm c \). Choosing
\( \bm c=\bm{\mathcal I}^{-1}\nabla g \) gives \( \Var\delta\ge\nabla g\T\bm{\mathcal I}^{-1}\nabla g \). The block-diagonal
form of @prp-opt-information gives the right-hand side.
:::

For \( g(\bm\theta)=\blambda\T\bbeta \) the bound is \( \sigma^2\blambda\T(\X\T\X)^{-1}\blambda \), exactly the variance of
\( \blambda\T\hbeta \). Least squares **attains** the bound, which gives a second proof of @thm-opt-umvue(a) in full
rank. For \( g=\sigma^2 \) the bound is \( 2\sigma^4/n \), while
\[
\Var(s^2)=\frac{\sigma^4}{(n-p)^2}\Var\{\chi^2(n-p)\}=\frac{2\sigma^4}{n-p}>\frac{2\sigma^4}n .
\]
So \( s^2 \) does not attain the bound, although by @thm-opt-umvue(b) no unbiased estimator does better. The
bound is simply not sharp here. Its ratio to the attained minimum, the **efficiency** of \( s^2 \), is
\( (n-p)/n \). That is \( 0.810 \) for the stack loss regression, and it tends to one as \( n\to\infty \)
with \( p \) fixed. The same happens for \( \sigma \): the UMVUE \( s/b_k \) has variance \( 0.0298\,\sigma^2 \)
when \( k=17 \), against the bound \( \sigma^2/(2n)=0.0238\,\sigma^2 \), an efficiency of
\( 0.798 \). The equality condition in the proof explains why. The bound is attained iff
\( \Var(\delta-\bm c\T\bm S)=0 \), that is, iff \( \delta=g(\bm\theta)+\nabla g\T\bm{\mathcal I}^{-1}\bm S \) with probability one.
For \( g=\sigma^2 \) the right-hand side is
\( \sigma^2+(2\sigma^4/n)\{-n/(2\sigma^2)+\be\T\be/(2\sigma^4)\}=\norm{\Y-\X\bbeta}^2/n \). This is not a statistic,
because it involves the unknown \( \bbeta \), so no unbiased estimator of \( \sigma^2 \) attains the bound.

::: {.warning}
The Cramér–Rao bound is a lower bound, not a target. An estimator that fails to attain it may still be
the best unbiased estimator, as \( s^2 \) is. When an estimator does attain it, the bound certifies
optimality without any appeal to completeness. The Lehmann–Scheffé route is stronger in one respect: it
needs no differentiability and covers every estimable function in every rank.
:::

## What changes without normality

Every result of this section uses the normal density. Without it the conclusions can fail. Under uniform
errors the midrange beats \( \bar{Y} \) (@exm-opt-midrange). Under Laplace errors the sample median has
asymptotic variance \( \sigma^2/(2n) \), half that of \( \bar{Y} \). Under the second-moment assumptions alone, the best
that can be said for \( s^2 \) is a result in the spirit of Gauss–Markov. Among quadratic unbiased estimators
\( \Y\T\A\Y \) of \( \sigma^2 \) with \( \A \) nonnegative definite, \( s^2 \) has the smallest variance when the errors have
zero excess kurtosis, or when all leverages are equal. @exr-opt-quadratic-sigma2 outlines this result
of Atiqullah (1962).

## Exercises

### A. Check your understanding

::: {#exr-opt-sufficient-full-rank}
[A1]

Assume \( \rank(\X)=p \). Use @lem-opt-factorization and the identity
\( \norm{\y-\X\bbeta}^2=\text{SSE}+(\hbeta-\bbeta)\T\X\T\X(\hbeta-\bbeta) \) to show directly that \( (\hbeta,\text{SSE}) \) is
sufficient for \( (\bbeta,\sigma^2) \). Is \( \hbeta \) alone sufficient?
:::

::: {#exr-opt-known-variance}
[A2]

Suppose \( \sigma^2=\sigma_0^2 \) is known. Show that \( \X\T\Y \) alone is complete and sufficient for \( \bbeta \), and that
\( \blambda\T\hbeta \) is still the UMVUE of every estimable \( \blambda\T\bbeta \). Why is \( \X\T\Y \) not sufficient when \( \sigma^2 \)
is unknown?
:::

### B. Practice

::: {#exr-opt-umvue-square}
[B1]

For estimable \( \blambda\T\bbeta \), find the UMVUE of \( (\blambda\T\bbeta)^2 \) in the normal linear model. Can it be
negative? Comment.
:::

::: {.solution}
\( \E(\blambda\T\hbeta)^2=(\blambda\T\bbeta)^2+\sigma^2\blambda\T(\X\T\X)\ginv\blambda \) and \( \E s^2=\sigma^2 \), so
\( (\blambda\T\hbeta)^2-s^2\blambda\T(\X\T\X)\ginv\blambda \) is unbiased. It is a function of \( (\M\Y,\text{SSE}) \), so it is the
UMVUE by @thm-opt-umvue(c). It is negative whenever \( \lvert t\rvert<1 \) for the \( t \) statistic of
@prp-opt-profile-t at \( \psi=0 \), which has positive probability. An estimator of a nonnegative quantity
that can be negative is clearly improvable, for example by its positive part. The positive part is biased
but has smaller mean squared error. Unbiasedness can be a costly requirement.
:::

::: {#exr-opt-information}
[B2]

Compute the observed information \( -\nabla^2\ell \) at \( (\bbeta,v) \), take its expectation, and confirm
@prp-opt-information by this second route.
:::

::: {#exr-opt-umvue-quantile}
[B3]

For the stack loss model, the \( 90\% \) quantile of the response at a new operating condition \( \x_0 \) is
\( \x_0\T\bbeta+1.2816\,\sigma \). Write down its UMVUE, and explain why \( \x_0\T\hbeta+1.2816\,s \) is biased. By how much,
in units of \( \sigma \)?
:::

### C. Going deeper

::: {#exr-opt-residual-direction}
[C1]

In the normal linear model with \( n-r\ge2 \), let \( \bU=(\I-\M)\Y/\norm{(\I-\M)\Y} \). Show that \( \bU \) is independent of
\( (\M\Y,\text{SSE}) \), and that its distribution (uniform on the unit sphere of \( \C(\X)\perpc \)) does not depend on
\( (\bbeta,\sigma^2) \). *Hint:* write \( (\I-\M)\Y=\Q_2\Z \) as in @exr-opt-reml, and use the rotation invariance of
\( \Normal(\bzero,\sigma^2\I) \).
:::

::: {.solution}
With \( \Q_2 \) as in @exr-opt-reml, \( (\I-\M)\Y=\Q_2\Z \) with \( \Z\sim\Normal_{n-r}(\bzero,\sigma^2\I) \), and \( \Z \) is independent of
\( \M\Y \) by @thm-mvn-independence, since \( \Cov(\M\Y,\Z)=\sigma^2\M\Q_2=\bm 0 \). Then
\( \text{SSE}=\norm{\Z}^2 \) and \( \bU=\Q_2\Z/\norm{\Z} \). For any orthogonal \( \bm H \), \( \bm H\Z \) has the same law as \( \Z \), so the
law of \( \Z/\norm{\Z} \) is invariant under all rotations. The only such law on the sphere is the uniform one. The
independence of the direction \( \Z/\norm{\Z} \) and the length \( \norm{\Z} \) follows from the polar form of the density
of \( \Z \), which is a function of \( \norm{\bz} \) only. Hence \( \bU \) is independent of \( (\M\Y,\text{SSE}) \), and its law
involves neither \( \bbeta \) nor \( \sigma^2 \).
:::

::: {#exr-opt-quadratic-sigma2}
[C2]

Assume only that the errors are independent with common variance \( \sigma^2 \), zero third moment and fourth
moment \( \mu_4 \). Let \( \Y\T\A\Y \) be an unbiased estimator of \( \sigma^2 \) with \( \A \) symmetric nonnegative definite.
(a) Show that \( \A\X=\bm0 \) and \( \tr\A=1 \). (b) Using @thm-rv-quadform-variance, show that
\( \Var(\Y\T\A\Y)=(\mu_4-3\sigma^4)\sum_ia_{ii}^2+2\sigma^4\tr(\A^2) \). (c) Show that \( \tr(\A^2)\ge1/(n-r) \), with
equality iff \( \A=(\I-\M)/(n-r) \). (d) Conclude that \( s^2 \) has the smallest variance in this class when
\( \mu_4=3\sigma^4 \). *Hint for (c):* \( \A=(\I-\M)\A(\I-\M) \), so the eigenvalues of \( \A \) live on \( \C(\X)\perpc \).
:::
