# The likelihood equations

In the normal linear model the maximum likelihood estimate was available in closed form,
because the log-likelihood was a quadratic in \( \bbeta \) (@thm-opt-mle). Outside that case
there is no closed form, and everything has to be read off the derivatives of the
log-likelihood. This section computes them. The answer has the same shape as the normal
equations \( \X\T(\y-\X\bbeta)=\bzero \), with two modifications: each observation is
weighted according to its variance, and the residual \( y_i-\mu_i \) is converted from the
mean scale to the linear-predictor scale.

Throughout, the model is @def-glm-model, \( \rank(\X)=p \), and \( \bbeta \) lies in the interior
of the admissible set \( \mathcal B \). Write
\[
\mu_i=h(\eta_i),\qquad \eta_i=\x_{(i)}\T\bbeta,\qquad
\frac{d\mu_i}{d\eta_i}=h'(\eta_i),
\]
and note that \( h'>0 \) everywhere, because \( g \) is strictly monotone (take \( g \) increasing;
a decreasing link only flips signs).

## The score

::: {#thm-glm-score}
[The likelihood equations and the information]

Define the diagonal matrices
\[
\W=\diag\Bigl\{\frac{w_i\,h'(\eta_i)^2}{V(\mu_i)}\Bigr\},\qquad
\bD=\diag\{h'(\eta_i)\} .
\]
Then, with \( \ell(\bbeta) \) the log-likelihood of @eq-glm-joint-loglik:

::: {.enumerate options="label=(\alph*)"}
1. the score is
   \[
   \bU(\bbeta)=\frac{\partial\ell}{\partial\bbeta}
   =\frac1\phi\,\X\T\W\bD^{-1}(\y-\bmu),
   \]{#eq-glm-score}

   so the likelihood equations are
   \( \sum_i x_{ij}\,w_i\,h'(\eta_i)\{y_i-\mu_i\}/V(\mu_i)=0 \) for \( j=1,\dots,p \), and they do
   not involve \( \phi \);

2. \( \E\bU=\bzero \) and the **expected (Fisher) information** is
   \[
   \boldsymbol{\mathcal I}(\bbeta)=\Cov\{\bU(\bbeta)\}=\frac1\phi\,\X\T\W\X ;
   \]{#eq-glm-information}

3. writing \( a(\eta)=h'(\eta)/V\{h(\eta)\} \), the **observed information** is
   \[
   -\frac{\partial^2\ell}{\partial\bbeta\,\partial\bbeta\T}
   =\frac1\phi\,\X\T\bigl(\W-\bD_r\bigr)\X,
   \]
   where \( \bD_r=\diag\{w_i(y_i-\mu_i)\,a'(\eta_i)\} \);

4. under the canonical link \( \bD_r=\bzero \), so the observed and expected information
   coincide, \( \W=\diag\{w_iV(\mu_i)\} \), and @eq-glm-score becomes
   \[
   \bU(\bbeta)=\frac1\phi\,\X\T\diag(\bw)(\y-\bmu).
   \]{#eq-glm-score-canonical}

:::

:::

::: {.proof}
By @eq-glm-joint-loglik, \( \ell=\phi^{-1}\sum_iw_i\{y_i\theta_i-b(\theta_i)\} \) plus a term free
of \( \bbeta \). The chain rule gives
\[
\frac{\partial\theta_i}{\partial\beta_j}
=\frac{d\theta_i}{d\mu_i}\,\frac{d\mu_i}{d\eta_i}\,\frac{\partial\eta_i}{\partial\beta_j}
=\frac{h'(\eta_i)}{V(\mu_i)}\,x_{ij},
\]
using \( d\theta/d\mu=1/V(\mu) \) from @def-glm-variance. Since \( b'(\theta_i)=\mu_i \),
\[
\frac{\partial\ell}{\partial\beta_j}
=\frac1\phi\sum_iw_i\{y_i-b'(\theta_i)\}\frac{\partial\theta_i}{\partial\beta_j}
=\frac1\phi\sum_i x_{ij}\,\frac{w_i\,h'(\eta_i)}{V(\mu_i)}\,(y_i-\mu_i),
\]
which is (a), since \( \W\bD^{-1}=\diag\{w_ih'(\eta_i)/V(\mu_i)\} \). The dispersion appears
only as the overall factor \( \phi^{-1} \), so setting the score to zero gives equations free
of \( \phi \).

(b) \( \E(Y_i-\mu_i)=0 \) gives \( \E\bU=\bzero \). The \( Y_i \) are independent with
\( \Var(Y_i)=\phi V(\mu_i)/w_i \) by @prp-glm-moments, so with the scalars
\( a_i=w_ih'(\eta_i)/V(\mu_i) \),
\[
\Cov(\bU)=\frac1{\phi^2}\sum_ia_i^2\,\Var(Y_i)\,\x_{(i)}\x_{(i)}\T
=\frac1\phi\sum_i\frac{w_ih'(\eta_i)^2}{V(\mu_i)}\,\x_{(i)}\x_{(i)}\T,
\]
which is @eq-glm-information.

(c) With \( a(\eta)=h'(\eta)/V\{h(\eta)\} \) we have
\( \partial\ell/\partial\beta_j=\phi^{-1}\sum_ix_{ij}w_i(y_i-\mu_i)a(\eta_i) \). Differentiating
again, and using \( \partial\mu_i/\partial\beta_k=h'(\eta_i)x_{ik} \),
\[
\frac{\partial^2\ell}{\partial\beta_j\partial\beta_k}
=\frac1\phi\sum_ix_{ij}x_{ik}w_i\bigl\{-h'(\eta_i)a(\eta_i)+(y_i-\mu_i)a'(\eta_i)\bigr\}.
\]
Now \( w_ih'(\eta_i)a(\eta_i)=w_ih'(\eta_i)^2/V(\mu_i) \) is the \( i \)th working weight, which
gives the stated form. Taking expectations kills the second term and recovers (b), a second
proof of @eq-glm-information.

(d) Under the canonical link \( \theta_i=\eta_i \), so \( h=b' \) and
\( h'(\eta)=b''(\theta)=V(\mu) \). Hence \( a(\eta)\equiv1 \), \( a'\equiv0 \), \( \bD_r=\bzero \), and
\( \W=\diag\{w_iV(\mu_i)^2/V(\mu_i)\}=\diag\{w_iV(\mu_i)\} \). Also
\( \W\bD^{-1}=\diag\{w_iV(\mu_i)/V(\mu_i)\}=\diag(\bw) \), which
is @eq-glm-score-canonical.
:::

For the normal family with the identity link and \( w_i=1 \), \( h'\equiv1 \) and
\( V\equiv1 \), so @eq-glm-score reads \( \X\T(\y-\X\bbeta)/\sigma^2 \) and
@eq-glm-information reads \( \X\T\X/\sigma^2 \): the score is the normal equations divided
by the variance, and @prp-opt-information is recovered.

Three features of @eq-glm-score are worth naming. The equations are **weighted** by
\( w_ih'(\eta_i)/V(\mu_i) \), so a more variable observation counts for less — the logic of
weighted least squares (@thm-het-wls(a)), with weights dictated by the family rather than
assumed. They are **free of \( \phi \)**, which multiplies the whole score, so \( \phi \) can
be estimated afterwards without disturbing \( \hbeta \)
([Section 34.5](05-deviance-and-inference.html)), just as least squares does not need
\( \sigma^2 \). And they involve the family **only through \( V \)**: two densities with the
same variance function give the same \( \hbeta \), the seed of the quasi-likelihood of
[Chapter 38](../ch38-quasi-likelihood/index.html) (@def-ql-quasi).

## The fitted marginals

::: {#cor-glm-marginals}
[Canonical links reproduce weighted column totals]

Under the canonical link, any solution \( \hbeta \) of the likelihood equations satisfies
\[
\sum_{i=1}^n w_i\,x_{ij}\,\hat\mu_i=\sum_{i=1}^n w_i\,x_{ij}\,y_i,\qquad j=1,\dots,p .
\]
In particular, if the model contains an intercept, \( \sum_iw_i\hat\mu_i=\sum_iw_iy_i \); and
if a column of \( \X \) is the indicator of a group, the fitted total over that group equals
the observed total.
:::

::: {.proof}
Immediate from @eq-glm-score-canonical: the likelihood equations are
\( \X\T\diag(\bw)(\y-\hat{\bmu})=\bzero \), which is the statement column by column.
:::

This is the exact analogue of the least squares identity \( \X\T\he=\bzero \)
of @thm-proj-ls-projection, and sharper than it looks: not merely that the residuals average
to zero, but that a canonical-link fit *reproduces a table of observed totals exactly*, on
the original scale of the data, whatever else is in the model.

::: {#exm-glm-marginals}
[Fitted votes by party identification]

The 1996 American National Election Studies subset records, for \( 944 \)
respondents, a reported vote (coded \( 1 \) for Dole), a party-identification score from
\( 0 \) (strong Democrat) to \( 6 \) (strong Republican), and age. Fit a logistic model with one
indicator column per party-identification level and age (in decades, centred at 45) as an
eighth column. The age coefficient is \( 0.1044 \) with standard error
\( 0.0730 \).

By @cor-glm-marginals the fitted probabilities must add up, within each of the seven
party-identification groups, to the observed number of Dole voters — though age is in the
model and varies within every group. Among the \( 200 \) strong
Democrats, \( 3 \) reported voting for Dole and the fitted probabilities sum
to \( 3 \) exactly; among the \( 175 \) strong Republicans observed
and fitted totals are both \( 167 \). The eighth column
gives \( \sum_i\text{age}_i\,y_i=\sum_i\text{age}_i\,\hat\mu_i=121.3000 \).

A probit refit destroys the identity: the group totals are then wrong by up
to \( 0.0966 \) and the age total by \( 0.6245 \). The
discrepancies are small, because the fits are similar, but nothing forces them to be zero.
The same contrast holds for the information, which under the logit link agrees with the
expected information to machine precision and under the probit link differs by a
relative \( 0.0056 \).
:::

```{.python .run #cell-score-fit}
import numpy as np
import statsmodels.api as sm

anes = sm.datasets.anes96.load_pandas().data
y = anes["vote"].to_numpy(float)                       # 1 if the respondent voted for Dole
pid = anes["PID"].to_numpy(int)                        # party identification, 0 (Dem) to 6 (Rep)
D = np.eye(7)[pid]                                     # one indicator column per level
age = (anes["age"].to_numpy(float) - 45.0) / 10.0
X = np.column_stack([D, age])                          # no separate intercept: D already sums to 1

logit = sm.GLM(y, X, family=sm.families.Binomial()).fit()
mu = logit.fittedvalues

print("level   n   observed Dole votes   fitted")
for k in range(7):
    inside = pid == k
    print(f"  {k}   {inside.sum():4d}      {y[inside].sum():8.0f}   {mu[inside].sum():12.6f}")
print(f"age column: X^T y = {age @ y:.6f},  X^T mu = {age @ mu:.6f}")
```

::: {.remark}
[Where the marginal property is used]

[Chapter 37](../ch37-counts/index.html) turns @cor-glm-marginals into a tool: for a log-linear model of a contingency
table the canonical link forces the fitted counts to match the observed margins of the terms
in the model, which is what makes log-linear and multinomial fits
agree (@thm-cnt-poisson-multinomial). In the binary case it means those totals carry no
information about goodness of fit — a point @prp-bin-fit takes up.
:::

## Concavity, existence and uniqueness

::: {#thm-glm-concave}
[The canonical-link log-likelihood is concave]

Under the canonical link, \( \ell(\bbeta) \) is concave on the convex set \( \mathcal B \). If in
addition \( \rank(\X)=p \) and \( V(\mu_i)>0 \) for all \( i \), then \( \ell \) is strictly concave, so
it has at most one stationary point, and any stationary point is the unique global maximizer.
:::

::: {.proof}
Under the canonical link \( \theta_i=\eta_i=\x_{(i)}\T\bbeta \), so by @eq-glm-joint-loglik
\[
\ell(\bbeta)=\frac1\phi\sum_{i=1}^n w_i\bigl\{y_i\,\x_{(i)}\T\bbeta-b(\x_{(i)}\T\bbeta)\bigr\}
+\text{const}.
\]
The first term is linear in \( \bbeta \). In the second, \( b \) is convex by @lem-glm-convex and
\( \bbeta\mapsto\x_{(i)}\T\bbeta \) is linear, so \( \bbeta\mapsto b(\x_{(i)}\T\bbeta) \) is convex;
multiplying by \( w_i/\phi>0 \) and summing preserves convexity, and the minus sign makes
\( \ell \) concave. The set \( \mathcal B \) is convex, being an intersection of the preimages of
the interval \( g(\mathcal M) \) under linear maps.

For strictness, @thm-glm-score(d) gives the Hessian
\( -\phi^{-1}\X\T\diag\{w_iV(\mu_i)\}\X \), which for \( \mathbf{v}\ne\bzero \) satisfies
\[
\mathbf{v}\T\Bigl(\frac1\phi\X\T\diag\{w_iV(\mu_i)\}\X\Bigr)\mathbf{v}
=\frac1\phi\sum_iw_iV(\mu_i)(\x_{(i)}\T\mathbf{v})^2>0,
\]
because \( \rank(\X)=p \) forces \( \X\mathbf{v}\ne\bzero \). So the Hessian is negative definite
everywhere on \( \mathcal B \) and \( \ell \) is strictly concave there.
:::

Concavity settles uniqueness but not existence: a strictly concave function on an open set
need not attain its supremum. Both failures happen in practice.

::: {#exm-glm-no-maximum}
[Two ways the maximum fails to exist]

*All counts zero.* For a Poisson log-linear model with an intercept, if every \( y_i=0 \) then
\( \ell(\bbeta)=-\sum_ie^{\x_{(i)}\T\bbeta} \), which increases to its supremum \( 0 \) as
\( \beta_0\to-\infty \) and is never attained.

*Separated binary data.* If there is a vector \( \mathbf{v} \) with \( \x_{(i)}\T\mathbf{v}>0 \)
whenever \( y_i=1 \) and \( \x_{(i)}\T\mathbf{v}<0 \) whenever \( y_i=0 \), then along the ray
\( t\mathbf{v} \) every fitted probability moves towards the observed response, so \( \ell \)
increases without reaching a maximum. [Chapter 35](../ch35-binary-responses/index.html)
proves that this is the only way the binary maximum likelihood estimate can fail to
exist (@thm-bin-separation).
:::

The general statement belongs to exponential-family theory. Under a canonical link
\( \bT=\X\T\diag(\bw)\y \) is the sufficient statistic of
[Section 34.1](01-exponential-dispersion-families.html), and by @cor-glm-marginals the
likelihood equations say \( \X\T\diag(\bw)\hat{\bmu}=\bT \). The attainable mean-value
vectors \( \{\X\T\diag(\bw)\bmu:\ \bmu\in\mathcal M^n\} \) form an open convex set, and the
maximum likelihood estimate exists precisely when \( \bT \) lies inside it rather than on its
boundary (Barndorff-Nielsen 1978, chapter 9; not reproduced here). Both examples above are
boundary cases: all-zero counts put \( \bT \) at the corner \( \bzero \), separation puts it on
a face. For non-canonical links even uniqueness can fail.

::: {.warning}
Concavity is a property of the canonical link, not of generalized linear models. Take the
gamma family with the identity link, where \( \ell=\phi^{-1}\sum_i\{-y_i/\mu_i-\log\mu_i\} \)
with \( \mu_i=\x_{(i)}\T\bbeta \). The second derivative of one term with respect to
\( \mu_i \) is \( (\mu_i-2y_i)/\mu_i^3 \), positive whenever \( \mu_i>2y_i \), so the Hessian is
a difference of nonnegative definite matrices and can be indefinite: an iteration can stop
at a local maximum. Wedderburn (1976) settles existence and uniqueness link by link, cleanly
for the binomial and Poisson with the standard links and not for the gamma with the identity
link.
:::

## Exercises

### A. Check your understanding

::: {#exr-glm-poisson-equations}
[A1]

Write the likelihood equations for a Poisson log-linear model with an intercept. Show that
they are \( \X\T(\y-\bmu)=\bzero \), that \( \sum_i\hat\mu_i=\sum_iy_i \), and that for the
intercept-only model \( \hat\mu_i=\bar y \). (This is @exr-lm-poisson-score, now a special case
of @cor-glm-marginals.)
:::

### B. Practice

::: {#exr-glm-logistic-equations}
[B1]

For grouped binomial data with \( y_i=s_i/m_i \), the logit link and \( w_i=m_i \), write the
likelihood equations in terms of the observed successes \( s_i \). Show that the working weight
is \( \W_{ii}=m_i\mu_i(1-\mu_i) \) and that the information is
\( \X\T\diag\{m_i\mu_i(1-\mu_i)\}\X \). Where have you seen that matrix before?
:::

::: {.solution}
The canonical form @eq-glm-score-canonical with \( w_i=m_i \) and \( \phi=1 \) gives
\( \sum_ix_{ij}m_i(s_i/m_i-\mu_i)=\sum_ix_{ij}(s_i-m_i\mu_i)=0 \): fitted and observed successes
agree in every column of \( \X \). By @thm-glm-score(d), \( \W=\diag\{w_iV(\mu_i)\} \) with
\( V(\mu)=\mu(1-\mu) \), giving \( \W_{ii}=m_i\mu_i(1-\mu_i) \) and the stated information. It is
the weight matrix of the empirical-logit weighted least squares fit
in @exr-tr-empirical-logit, with \( \mu_i \) in place of \( \hat p_i \): that procedure is one step
of the iteration of [Section 34.4](04-irls.html) started from the data.
:::

::: {#exr-glm-canonical-characterization}
[B2]

Show that the observed information equals the expected information *identically in \( \y \)*
if and only if the link is canonical (up to an affine change of \( \eta \)). Assume that no row
of \( \X \) is zero, an observation with \( \x_{(i)}=\bzero \) carrying no information
about \( \bbeta \). Use @thm-glm-score(c).
:::

::: {.solution}
By @thm-glm-score(c) the two agree for all \( \y \) iff \( \X\T\bD_r\X=\bzero \) for all \( \y \),
where \( (\bD_r)_{ii}=w_i(y_i-\mu_i)a'(\eta_i) \). Vary \( y_i \) over two admissible values with
the other responses held fixed: subtracting the two identities leaves
\( \delta_i\,a'(\eta_i)\,\x_{(i)}\x_{(i)}\T=\bzero \) with
\( \delta_i=w_i(y_i-y_i')\ne0 \), and since
\( \x_{(i)}\ne\bzero \) this forces \( a'(\eta_i)=0 \) for every \( i \). (Two admissible values
always exist: the \( y_i \) range over an interval for the normal and gamma families, and over
\( \{0,1/m,\dots,1\} \) or \( \{0,1,2,\dots\} \) for the binomial and Poisson.) So \( a \) is
constant, say \( a\equiv1/\lambda \), that is
\( d\mu/d\eta=V(\mu)/\lambda \). Comparing with the canonical link, for which
\( d\mu/d\theta=V(\mu) \), gives \( \eta=\lambda\theta+\text{const} \).
:::

### C. Going deeper

::: {#exr-glm-concavity-links}
[C1]

For a binary response, \( \ell(\bbeta)=\sum_i\{y_i\log h(\eta_i)+(1-y_i)\log(1-h(\eta_i))\} \).
Show that \( \ell \) is concave for every design and every data set if and only if \( \log h \) and
\( \log(1-h) \) are both concave functions of \( \eta \). Verify that the logit link satisfies
this, and that \( h(\eta)=\Phi(\eta) \) does too, using the fact that the standard normal
density is log-concave.
:::

::: {#exr-glm-offset}
[C2]

An **offset** is a known column added to the linear predictor with coefficient fixed at
\( 1 \): \( \eta_i=o_i+\x_{(i)}\T\bbeta \). Show that @thm-glm-score holds verbatim with
\( \eta_i \) redefined this way, that @cor-glm-marginals is unchanged, and that the offset does
*not* appear among the columns whose totals are reproduced. Where would such a term come
from in a Poisson model for rates? ([Chapter 37](../ch37-counts/index.html) answers this with @prp-cnt-offsets.)
:::
