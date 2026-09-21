# Quasi-likelihood estimating equations

Estimation here starts from @eq-ql-equations and not from a density. There are two useful
ways to read it: as the derivative of a function that is not quite a log-likelihood but
close enough to carry the usual arguments, and as the best member of a large family of
estimating equations. This section develops both, and they meet in the middle.

## The quasi-likelihood function

::: {#def-ql-quasi}
[Quasi-likelihood]

Under @def-ql-variance-function, the **quasi-likelihood** of a mean \( \mu\in\mathcal M \)
given a response \( y \) in the closure of \( \mathcal M \), with weight \( w \), is
\[
Q(\mu;y)=\frac{w}{\phi}\int_{y}^{\mu}\frac{y-t}{V(t)}\,dt ,
\]{#eq-ql-quasi}

assumed finite, the integral being read as a limit at an endpoint of \( \mathcal M \)
where it converges — as it does at \( y=0 \) for \( V(t)=t \) and at \( y\in\{0,1\} \)
for \( V(t)=t(1-t) \). The quasi-likelihood of the model is
\( Q(\bbeta)=\sum_{i=1}^nQ(\mu_i;y_i) \) with \( \mu_i=h(\x_{(i)}\T\bbeta) \) and weight
\( w_i \).
:::

The definition is designed to make one identity true by construction, and everything else
follows from it.

::: {#prp-ql-properties}
[What the quasi-likelihood does and does not share with a log-likelihood]

Write \( Q'=\partial Q/\partial\mu \). Then, for every variance specification:

::: {.enumerate options="label=(\alph*)"}
1. \( Q'(\mu;y)=w(y-\mu)/\{\phi V(\mu)\} \), and \( Q(\mu;y)\le0 \) with equality only at
   \( \mu=y \): the quasi-likelihood is maximized at the observation, as a saturated fit
   would be;

2. \( \E\{Q'(\mu;Y)\}=0 \) and
   \( \Var\{Q'(\mu;Y)\}=-\E\{\partial^2Q/\partial\mu^2\}=w/\{\phi V(\mu)\} \), the two
   Bartlett identities;

3. if \( (V,\phi,w) \) are those of an exponential dispersion family (@def-glm-edf), then
   \( Q(\mu;y)=\ell(\mu;y)-\ell(y;y) \), so \( Q \) is the log-likelihood measured from its
   saturated value, and \( -2\phi\,Q(\bbeta)=D(\y,\bmu) \) is the deviance
   of @eq-glm-deviance.
:::

:::

::: {.proof}
(a) Differentiating @eq-ql-quasi in its upper limit gives the stated \( Q' \). Because
\( V>0 \), the integrand \( (y-t)/V(t) \) is positive for \( t<y \) and negative for
\( t>y \), so for \( \mu>y \) the integral is of a negative function and for
\( \mu<y \) it is minus the integral of a positive one: both give \( Q\le0 \), with
equality only for an empty interval.

(b) \( \E(Y-\mu)=0 \) gives the first identity, and \( \Var(Y)=\phi V(\mu)/w \) gives
\( \Var(Q')=w/\{\phi V(\mu)\} \). Differentiating \( Q' \) again,
\[
\frac{\partial^2Q}{\partial\mu^2}
=-\frac{w}{\phi}\cdot\frac{V(\mu)+(y-\mu)V'(\mu)}{V(\mu)^2},
\]
whose expectation is \( -w/\{\phi V(\mu)\} \), since the term carrying \( y-\mu \) has mean
zero.

(c) For an exponential dispersion family the chain rule used in the proof
of @thm-glm-score gives \( \partial\ell/\partial\mu=w(y-\mu)/\{\phi V(\mu)\} \), which
is \( Q' \); two functions of \( \mu \) with the same derivative differ by a constant, and
both sides vanish at \( \mu=y \). Comparing with @eq-glm-deviance gives
\( D=-2\phi Q \).

:::

The converse of (c) fails: \( \exp Q \) is in general not a density in
\( y \) (@exr-ql-integrability). Part (c) is nevertheless why the objects here look
familiar — everything the exponential family contributed to @def-glm-deviance is already in
the integral @eq-ql-quasi — and for a variance function from no family @eq-ql-quasi still
measures the discrepancy between a fit and an observation, the **quasi-deviance**
\[
D(\y,\bmu)=-2\phi\,Q(\bbeta)=2\sum_{i=1}^nw_i\int_{\mu_i}^{y_i}\frac{y_i-t}{V(t)}\,dt\ \ge 0 ,
\]{#eq-ql-quasi-deviance}

which does not involve \( \phi \) and vanishes only at \( \bmu=\y \).
[Figure 38.2.1](02-estimating-equations.html#fig-ql-shapes) draws one term for four power
variance functions: only \( V\equiv1 \) penalizes symmetrically, and under
\( V(\mu)=\mu \) a fitted mean one unit too small costs \( 2.043 \)
times as much as one unit too large, against \( 4.254 \) for
\( V(\mu)=\mu^2 \) and \( 9.000 \) for \( V(\mu)=\mu^3 \). A
variance function that grows with the mean is a statement about which errors matter.

::: {when-format="html"}
![**Figure 38.2.1.** One term of the quasi-deviance @eq-ql-quasi-deviance at \( y=2 \),
\( w=1 \), against the fitted mean, for four power variance functions. Only the constant
one penalizes symmetrically.](quasi_deviance_shapes.svg){#fig-ql-shapes width=58%}
:::

::: {when-format="pdf"}
![One term of the quasi-deviance @eq-ql-quasi-deviance at \( y=2 \),
\( w=1 \), against the fitted mean, for four power variance functions. Only the constant
one penalizes symmetrically.](quasi_deviance_shapes.pdf){width=58%}
:::

```{.python .run #cell-variance-shapes-qdev}
import numpy as np
from scipy import integrate

def quasi_deviance_term(y, mu, V):
    """2 * integral from mu to y of (y - t)/V(t) dt: one observation's quasi-deviance."""
    value, _ = integrate.quad(lambda t: (y - t) / V(t), mu, y)
    return 2 * value

y0 = 2.0
closed_form = {
    "V = 1": (lambda m: (y0 - m) ** 2, lambda t: 1.0),
    "V = mu": (lambda m: 2 * (y0 * np.log(y0 / m) - (y0 - m)), lambda t: t),
    "V = mu^2": (lambda m: 2 * (-np.log(y0 / m) + (y0 - m) / m), lambda t: t ** 2),
    "V = mu^3": (lambda m: (y0 - m) ** 2 / (y0 * m ** 2), lambda t: t ** 3),
}
for name, (formula, V) in closed_form.items():
    print(f"{name:9s}  at mu = 0.8: {quasi_deviance_term(y0, 0.8, V):.6f}"
          f"   closed form {formula(0.8):.6f}")
```

## The quasi-score

::: {#thm-ql-score}
[The quasi-score]

Let \( \bU(\bbeta)=\partial Q/\partial\bbeta \) be the **quasi-score**, and keep the
notation \( \W=\diag\{w_ih'(\eta_i)^2/V(\mu_i)\} \) and \( \bD=\diag\{h'(\eta_i)\} \)
of @thm-glm-score. Write \( \A_n=\phi^{-1}\X\T\W\X \). Then:

::: {.enumerate options="label=(\alph*)"}
1. \( \bU(\bbeta)=\phi^{-1}\X\T\W\bD^{-1}(\y-\bmu) \), which is @eq-glm-score exactly.
   The estimating equations \( \bU(\hbeta)=\bzero \) are @eq-ql-equations, they do not
   involve \( \phi \), and they are solved by the iteratively reweighted least squares
   of @thm-glm-irls;

2. if the mean model is correct then \( \E\bU(\bbeta^0)=\bzero \), *whatever the true
   variances are*: the quasi-score is an unbiased estimating function;

3. \( -\E\bigl\{\partial\bU/\partial\bbeta\T\bigr\}=\A_n \), again whatever the true
   variances are;

4. if in addition the variance specification @eq-ql-specification is correct, then
   \( \Cov\{\bU(\bbeta^0)\}=\A_n \) as well: sensitivity and variability agree, which is
   the second Bartlett identity for \( Q \);

5. if instead \( \Var(Y_i)=\sigma_i^2 \) for arbitrary \( \sigma_i^2 \), then
   \[
   \Cov\{\bU(\bbeta^0)\}=\B_n
   =\frac1{\phi^2}\sum_{i=1}^n\Bigl\{\frac{w_ih'(\eta_i)}{V(\mu_i)}\Bigr\}^{2}
   \sigma_i^2\,\x_{(i)}\x_{(i)}\T ,
   \]{#eq-ql-meat}

   which equals \( \A_n \) if and only if
   \( \sum_i\{w_ih'_i/V_i\}^2\{\sigma_i^2-\phi V_i/w_i\}\x_{(i)}\x_{(i)}\T=\bzero \).
:::

:::

::: {.proof}
(a) By @prp-ql-properties(a) and the chain rule,
\[
\frac{\partial Q}{\partial\beta_j}
=\sum_{i=1}^nQ'(\mu_i;y_i)\,\frac{d\mu_i}{d\eta_i}\,x_{ij}
=\frac1\phi\sum_{i=1}^n x_{ij}\,\frac{w_i\,h'(\eta_i)}{V(\mu_i)}\,(y_i-\mu_i),
\]
which is the coordinate form of @eq-glm-score. The map from \( (V,h,\bw,\phi) \) to it is
the one computed in @thm-glm-score(a), so the matrix form, the absence of \( \phi \) and
the iteration of @thm-glm-irls transfer unchanged; none of the three used the density.

(b) Every term carries the factor \( y_i-\mu_i \), which has mean zero when
\( \mu_i=h(\x_{(i)}\T\bbeta^0) \). No assumption about \( \Var(Y_i) \) is used.

(c) Write \( a(\eta)=h'(\eta)/V\{h(\eta)\} \), so that
\( \bU=\phi^{-1}\sum_iw_ia(\eta_i)(y_i-\mu_i)\x_{(i)} \). Differentiating,
\[
\frac{\partial\bU}{\partial\bbeta\T}
=\frac1\phi\sum_{i=1}^nw_i\bigl\{a'(\eta_i)(y_i-\mu_i)-a(\eta_i)h'(\eta_i)\bigr\}
\x_{(i)}\x_{(i)}\T ,
\]
the computation of @thm-glm-score(c). Expectation annihilates the first term whatever the
variances, and \( w_ia(\eta_i)h'(\eta_i)=\W_{ii} \) leaves \( -\A_n \).

(d) and (e). The summands of \( \bU \) are independent with mean zero, so
\( \Cov(\bU)=\phi^{-2}\sum_iw_i^2a(\eta_i)^2\Var(Y_i)\x_{(i)}\x_{(i)}\T \), which
is @eq-ql-meat. Substituting \( \sigma_i^2=\phi V(\mu_i)/w_i \) turns the \( i \)th
coefficient into
\( \phi^{-2}\,w_i^2\{h'_i/V_i\}^2\,\phi V_i/w_i=\phi^{-1}w_ih_i'^2/V_i \), giving
\( \A_n \). The stated condition is the difference of the two sums.
:::

::: {.remark}
[Existence and uniqueness]

Nothing above guarantees that @eq-ql-equations has a solution, or only one. Under the
canonical link *of the variance function*, \( g'(\mu)=1/V(\mu) \), the proof of (c) has
\( a\equiv1 \) and \( a'\equiv0 \), so the Jacobian is \( -\A_n \), negative definite
whenever \( \rank(\X)=p \); \( Q \) is then strictly concave and a root, if one exists,
is the unique maximizer, by @thm-glm-concave with \( b'' \) replaced by \( V \). For
other links \( Q \) can have several stationary points, as the warning
after @thm-glm-concave showed.
:::

::: {#exm-ql-strikes-fit}
[Strike durations]

Kennan's data record the length in days of \( 62 \) strikes in US
manufacturing, with a measure of unanticipated industrial production at the time each
strike began. Take a log link and \( V(\mu)=\mu \). The production coefficient is
\( -7.680 \); replacing \( V(\mu)=\mu \) by \( V(\mu)=7.5\mu \)
returns the same number to ten decimal places, since \( \phi \) cancels
from @eq-ql-equations, while \( V(\mu)=\mu^2 \) moves it to
\( -9.353 \), a quadratic variance function downweighting the long
strikes.

The mean duration is \( 42.66 \) days and the sample variance
\( 2102.7 \), nearly fifty times as large: a Poisson likelihood is no
serious description of these data, and the question is how much of the Poisson machinery
survives anyway.
:::

```{.python .run #cell-strikes-fit}
import numpy as np
import statsmodels.api as sm

data = sm.datasets.strikes.load_pandas().data
y = data["duration"].to_numpy(float)          # length of the strike, in days
x = data["iprod"].to_numpy(float)             # unanticipated industrial production
X = np.column_stack([np.ones(len(y)), x])
n, p = X.shape

def quasi_irls(X, y, V, steps=60):
    """Solve the quasi-score equations for a log link and variance function V."""
    beta = np.zeros(X.shape[1])
    beta[0] = np.log(y.mean())
    for _ in range(steps):
        eta = X @ beta
        mu = np.exp(eta)
        W = mu ** 2 / V(mu)                   # working weight (dmu/deta)^2 / V(mu)
        z = eta + (y - mu) / mu               # working response
        beta = np.linalg.solve((X.T * W) @ X, (X.T * W) @ z)
    return beta

beta_lin = quasi_irls(X, y, lambda m: m)              # V(mu) = mu
beta_scaled = quasi_irls(X, y, lambda m: 7.5 * m)     # V(mu) = phi mu: same equations
beta_sq = quasi_irls(X, y, lambda m: m ** 2)          # V(mu) = mu^2: different equations
print(f"V = mu       {beta_lin}")
print(f"V = 7.5 mu   {beta_scaled}")
print(f"V = mu^2     {beta_sq}")
```

## Optimality among linear estimating functions

The second view asks what else might have been used. Any function of the data and the
parameter with mean zero at \( \bbeta^0 \) is an **unbiased estimating function**, and
setting it to zero defines an estimator; the quasi-score is the best of a natural class.

::: {#prp-ql-optimality}
[Godambe optimality]

Assume @eq-ql-specification, and consider **linear estimating functions**
\[
\mathbf{g}(\bbeta)=\sum_{i=1}^n\mathbf{a}_i(\bbeta)\,\bigl\{y_i-\mu_i(\bbeta)\bigr\}
=\mathbf{C}(\bbeta)\,(\y-\bmu),
\]
where \( \mathbf{C}=[\mathbf{a}_1,\dots,\mathbf{a}_n] \) is \( p\times n \), differentiable in
\( \bbeta \), and does not depend on \( \y \). Let
\( \bS=-\E\{\partial\mathbf{g}/\partial\bbeta\T\} \) and suppose \( \bS \) is nonsingular.
Define the **Godambe information**
\[
\boldsymbol{\mathcal J}(\mathbf{g})
=\bS\T\,\Cov\{\mathbf{g}(\bbeta^0)\}^{-1}\,\bS .
\]

::: {.enumerate options="label=(\alph*)"}
1. Every such \( \mathbf{g} \) is unbiased, and \( \bS=\mathbf{C}\bD\X \).

2. \( \boldsymbol{\mathcal J}(\mathbf{g})\preceq\phi^{-1}\X\T\W\X
   =\boldsymbol{\mathcal J}(\bU) \) in the order of nonnegative definiteness, for every
   linear estimating function \( \mathbf{g} \).

3. Equality holds if and only if \( \mathbf{C}\T=\bSigma^{-1}\bD\X\,\mathbf{N} \) for a
   nonsingular \( p\times p \) matrix \( \mathbf{N} \), with \( \bSigma=\Cov(\Y) \): the
   coefficients are those of the quasi-score up to one common nonsingular transformation,
   and \( \mathbf{g} \) then has the same roots as \( \bU \).
:::

:::

::: {.proof}
(a) Unbiasedness follows from \( \E(\y-\bmu)=\bzero \). Differentiating,
\[
\frac{\partial\mathbf{g}}{\partial\bbeta\T}
=\sum_{i=1}^n\frac{\partial\mathbf{a}_i}{\partial\bbeta\T}(y_i-\mu_i)
-\mathbf{C}\,\frac{\partial\bmu}{\partial\bbeta\T},
\]
read entry by entry. Its expectation vanishes because \( \mathbf{a}_i \) does not depend
on \( \y \), and \( \partial\bmu/\partial\bbeta\T=\bD\X \), so
\( \bS=\mathbf{C}\bD\X \).

(b) and (c). Under @eq-ql-specification,
\( \bSigma=\diag\{\phi V(\mu_i)/w_i\} \) is positive definite. Set
\[
\mathbf{F}=\bSigma^{1/2}\mathbf{C}\T,\qquad \mathbf{K}=\bSigma^{-1/2}\bD\X,
\]
both \( n\times p \). Then \( \Cov(\mathbf{g})=\mathbf{C}\bSigma\mathbf{C}\T=\mathbf{F}\T\mathbf{F} \) and
\( \bS=\mathbf{C}\bD\X=\mathbf{F}\T\mathbf{K} \), so
\[
\boldsymbol{\mathcal J}(\mathbf{g})
=\mathbf{K}\T\mathbf{F}(\mathbf{F}\T\mathbf{F})^{-1}\mathbf{F}\T\mathbf{K}
=\mathbf{K}\T\bP\mathbf{K},
\]
where \( \bP=\mathbf{F}(\mathbf{F}\T\mathbf{F})^{-1}\mathbf{F}\T \) is the orthogonal
projection onto \( \C(\mathbf{F}) \) (@thm-proj-M-formula), and
\( \mathbf{F}\T\mathbf{F} \) is nonsingular because \( \bS \) is. As \( \I-\bP \) is
symmetric idempotent,
\( \mathbf{v}\T(\I-\bP)\mathbf{v}=\norm{(\I-\bP)\mathbf{v}}^2\ge0 \), so
\( \bP\preceq\I \) and
\( \boldsymbol{\mathcal J}(\mathbf{g})\preceq\mathbf{K}\T\mathbf{K} \). Finally
\[
\mathbf{K}\T\mathbf{K}=\X\T\bD\bSigma^{-1}\bD\X
=\sum_{i=1}^n\frac{w_i\,h'(\eta_i)^2}{\phi\,V(\mu_i)}\,\x_{(i)}\x_{(i)}\T
=\frac1\phi\X\T\W\X ,
\]
and @thm-ql-score(a),(d) identify this as \( \boldsymbol{\mathcal J}(\bU) \), since for the
quasi-score sensitivity and variability are both \( \A_n \).

Equality in (b) means \( \mathbf{K}\T(\I-\bP)\mathbf{K}=\bzero \), that is
\( \C(\mathbf{K})\subseteq\C(\mathbf{F}) \). Both matrices have \( p \) columns and,
because \( \mathbf{F}\T\mathbf{K} \) is nonsingular, rank \( p \); so their column
spaces coincide and \( \mathbf{K}=\mathbf{F}\mathbf{N} \) for a nonsingular
\( \mathbf{N} \). Written out this is
\( \mathbf{C}\T=\bSigma^{-1}\bD\X\mathbf{N}^{-1} \), whose \( i \)th row is a fixed
multiple of \( w_ih'(\eta_i)\x_{(i)}\T/V(\mu_i) \). Then
\( \mathbf{g}=\mathbf{N}^{-\T}\bU \), because
\( \X\T\bD\bSigma^{-1}(\y-\bmu) \) already carries the factor \( \phi^{-1} \) that
\( \bU \) has; so \( \mathbf{g} \) vanishes exactly where \( \bU \) does.
:::

::: {.warning}
The optimality is *within the class of linear estimating functions*. A known true density
would have a score nonlinear in \( \y \) with a larger Godambe information — the NB1
score @eq-cnt-nb1-score is an instance — so "best linear" is the honest reading, as
in @thm-opt-gauss-markov. Optimality also assumes @eq-ql-specification; under a wrong
variance function the quasi-score stays unbiased but is no longer optimal, and
[Section 38.3](03-sandwich.html) prices the loss.
:::

## Exercises

### A. Check your understanding

::: {#exr-ql-quasi-table}
[A1]

Compute \( \phi w^{-1}Q(\mu;y) \) from @eq-ql-quasi for \( V(\mu)=1 \), \( V(\mu)=\mu \),
\( V(\mu)=\mu(1-\mu) \) and \( V(\mu)=\mu^2 \), and check that in each case it is the
log-likelihood kernel of the corresponding family, up to a function of \( y \) alone.
:::

::: {.solution}
For \( V=1 \), \( \int_y^{\mu}(y-t)\,dt=-(y-\mu)^2/2 \). For \( V(t)=t \) the integrand is
\( y/t-1 \), giving the Poisson kernel \( y\log\mu-\mu \). For \( V(t)=t(1-t) \), partial
fractions give \( y/t-(1-y)/(1-t) \), whose integral is the binomial kernel
\( y\log\mu+(1-y)\log(1-\mu) \). For \( V(t)=t^2 \) the integrand is \( y/t^2-1/t \),
giving the gamma kernel \( -y/\mu-\log\mu \).
:::

::: {#exr-ql-nb-quasi}
[A2]

Show that the NB2 variance function \( V(\mu)=\mu+\mu^2/\kappa \) with \( \kappa \) known
gives \( \phi w^{-1}Q(\mu;y)=y\log\{\mu/(\mu+\kappa)\}+\kappa\log\{\kappa/(\mu+\kappa)\} \)
up to a function of \( y \), which is the negative binomial log-likelihood kernel
of @eq-cnt-nb-pmf. NB2 with \( \kappa \) known is an exponential dispersion family, so
@prp-ql-properties(c) applies; explain why the *quasi-Poisson* specification
\( V(\mu)=\mu \) with \( \phi \) free does not enjoy the same conclusion, and relate
this to the remark of [Section 38.1](01-variance-functions.html).
:::

### B. Practice

::: {#exr-ql-hand-optimality}
[B1]

Take \( n \) independent responses with \( \E Y_i=\mu \) (one unknown parameter, no
regressors) and \( \Var Y_i=\phi V(\mu)/w_i \). Show that the quasi-score estimate
of @thm-ql-score is the weighted mean \( \sum_iw_iy_i/\sum_iw_i \), and verify @prp-ql-optimality directly: among
estimating functions \( \sum_ia_i(y_i-\mu) \) with constants \( a_i \), the Godambe
information is maximized by \( a_i\propto w_i \).
:::

::: {.solution}
Parameterizing by \( \mu \) makes \( h'\equiv1 \), so @eq-ql-equations reads
\( \sum_iw_i(y_i-\mu)/V(\mu)=0 \) and \( \hat\mu=\sum_iw_iy_i/\sum_iw_i \). For a general
\( \mathbf{a} \) the sensitivity is \( S=\sum_ia_i \) and the variability is
\( \phi V(\mu)\sum_ia_i^2/w_i \), so
\[
J=\frac{\bigl(\sum_ia_i\bigr)^2}{\phi V(\mu)\sum_ia_i^2/w_i}
\le\frac{\sum_iw_i}{\phi V(\mu)}
\]
by the Cauchy–Schwarz inequality applied to
\( \sum_ia_i=\sum_i(a_i/\sqrt{w_i})\sqrt{w_i} \), with equality if and only if
\( a_i\propto w_i \). The bound is \( \phi^{-1}\X\T\W\X \) with \( \X=\bone \).
:::

::: {#exr-ql-marginals}
[B2]

Suppose the link satisfies \( g'(\mu)=1/V(\mu) \), the canonical link of the variance
function. Show that @eq-ql-equations becomes
\( \sum_iw_ix_{ij}(y_i-\hat\mu_i)=0 \) for every \( j \), so that the fit reproduces the
weighted column totals of \( \X \) exactly, as in @cor-glm-marginals. Check the condition
for \( V(\mu)=\mu \) with a log link and for \( V(\mu)=\mu(1-\mu) \) with a logit link, and
find the canonical link of \( V(\mu)=\mu^2 \).
:::

::: {.solution}
Since \( h \) is the inverse of \( g \), \( h'(\eta)=1/g'(\mu)=V(\mu) \), so the
coefficient \( w_ih'(\eta_i)/V(\mu_i) \) in @eq-ql-equations is \( w_i \). For
\( V(\mu)=\mu \) the log link has \( g'(\mu)=1/\mu \); for \( V(\mu)=\mu(1-\mu) \) the logit
link has \( g'(\mu)=1/\{\mu(1-\mu)\} \); for \( V(\mu)=\mu^2 \) the condition
\( g'(\mu)=\mu^{-2} \) gives \( g(\mu)=-1/\mu \), the canonical link of the gamma family.
The conclusion is the quasi-likelihood form of @cor-glm-marginals, proved there for
exponential dispersion families but needing only @eq-ql-equations.
:::

### C. Going deeper

::: {#exr-ql-integrability}
[C1]

Take the quasi-binomial specification \( V(t)=t(1-t) \) on \( (0,1) \) with \( w=1 \)
and any \( \phi>0 \). (a) Show that \( y\mapsto\exp\{Q(\mu;y)\} \) is integrable over
\( (0,1) \), so normalizes to a density \( f \). (b) Show that a distribution on
\( [0,1] \) with mean \( m \) has variance at most \( m(1-m) \), with equality only on
\( \{0,1\} \); deduce that \( f \) has variance strictly below \( m(1-m) \) at its own
mean \( m \), never \( \phi\,m(1-m) \) with \( \phi\ge1 \). (c) Show that \( f \) is
symmetric about \( 1/2 \) when \( \mu=1/2 \), and check numerically that at
\( \mu=0.3 \), \( \phi=1 \) it has mean \( 0.4408 \) and
variance \( 0.0687 \), against \( \mu=0.3 \) and
\( \phi\mu(1-\mu)=0.21 \). Conclude that integrability of \( \exp Q \) does not make a
variance specification a family.
:::

::: {.solution}
(a) \( Q\le0 \) by @prp-ql-properties(a), so \( \exp Q\le1 \) on a bounded interval.

(b) For \( Y\in[0,1] \), \( \E\{Y(1-Y)\}\ge0 \) gives \( \E Y^2\le m \) and hence
\( \Var Y\le m-m^2 \); equality forces \( Y(1-Y)=0 \) almost surely, which a density on
the open interval excludes.

(c) By @exr-ql-quasi-table the exponent is
\( \phi^{-1}\{y\log\mu+(1-y)\log(1-\mu)-y\log y-(1-y)\log(1-y)\} \), unchanged by
\( y\mapsto1-y \) when \( \mu=1/2 \). Elsewhere neither moment is reproduced: \( \exp Q \)
normalizes without the result being a distribution with the moments
of @eq-ql-specification.
:::

