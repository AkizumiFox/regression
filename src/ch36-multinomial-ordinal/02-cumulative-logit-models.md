# Cumulative logit models

Suppose the \( c \) categories are ordered: none, mild, severe; strongly disagree through
strongly agree; extremely liberal through extremely conservative. The baseline-category
model of [Section 36.1](01-baseline-category-logits.html) fits such data perfectly well
— it fits *any* categorical response — but throws the ordering away. It spends
\( (c-1)p \) coefficients where a model that knows the order can often spend \( p \), and its
fit is unchanged if the category labels are shuffled. That last property is the
diagnosis: a likelihood invariant to permuting the categories cannot be using the order.

## Cutting the scale

The device that uses the ordering is to *cut*: for each \( r \) between \( 1 \) and \( c-1 \),
collapse the response into the binary event \( \{Y_i\le r\} \) and model its probability.
Write
\[
\gamma_{ir}=P(Y_i\le r\mid\x_{(i)})=\pi_{i1}+\dots+\pi_{ir},
\qquad r=1,\dots,c-1,
\]{#eq-mlt-cumulative-prob}

for the **cumulative probabilities**, with \( \gamma_{i0}=0 \) and \( \gamma_{ic}=1 \). The
category probabilities are recovered by differencing,
\( \pi_{ir}=\gamma_{ir}-\gamma_{i,r-1} \), so modelling the \( c-1 \) cumulative
probabilities is modelling the whole distribution.

::: {#def-mlt-cumulative}
[Cumulative link model]

Let \( G \) be a continuous, strictly increasing distribution function on \( \Real \) and let
\( \x_{(i)} \) be a vector of \( p \) regressors *without* an intercept. The **cumulative link
model** is
\[
G^{-1}(\gamma_{ir})=\theta_r-\x_{(i)}\T\bbeta_r ,\qquad r=1,\dots,c-1,
\]{#eq-mlt-cumulative}

with **cutpoints** \( \theta_1,\dots,\theta_{c-1} \) and slope vectors
\( \bbeta_1,\dots,\bbeta_{c-1} \). With \( G=\Lambda \), the standard logistic distribution
function \( \Lambda(z)=e^{z}/(1+e^{z}) \), @eq-mlt-cumulative is the **cumulative logit
model**
\[
\log\frac{P(Y_i\le r\mid\x_{(i)})}{P(Y_i> r\mid\x_{(i)})}=\theta_r-\x_{(i)}\T\bbeta_r ;
\]
with \( G=\Phi \) it is the **cumulative probit model**, and with
\( G(z)=1-\exp(-e^{z}) \) the **cumulative complementary log-log model**. The special case
\( \bbeta_1=\dots=\bbeta_{c-1}=\bbeta \) is the **proportional odds model** when \( G=\Lambda \),
and is the subject of [Section 36.3](03-proportional-odds.html).
:::

Two conventions deserve comment. The intercept is not in \( \x_{(i)} \) because the
cutpoints already play that role: shifting every \( \theta_r \) and a hypothetical intercept
by the same constant would change nothing, so the two cannot both be free. And the minus
sign before \( \x_{(i)}\T\bbeta_r \) is worth the confusion it causes once, because it
makes a *positive* coefficient mean "larger values of this regressor push the response
towards the higher categories". The cutpoints are unconstrained parameters only if the
model is coherent, which is the one real difficulty of this family.

## The ordering constraint

::: {#prp-mlt-ordering}
[Cutpoints must stay in order]

@eq-mlt-cumulative defines a probability distribution over the \( c \) categories at
\( \x_{(i)} \) if and only if
\[
\theta_1-\x_{(i)}\T\bbeta_1<\theta_2-\x_{(i)}\T\bbeta_2<\dots<\theta_{c-1}-\x_{(i)}\T\bbeta_{c-1}.
\]{#eq-mlt-ordering}

If the slopes are common, \( \bbeta_r=\bbeta \) for all \( r \), this reduces to
\( \theta_1<\dots<\theta_{c-1} \), a constraint on the parameters alone. If the slopes
differ, the constraint involves \( \x_{(i)} \), and for any two cuts with
\( \bbeta_r\neq\bbeta_{r+1} \) there are values of \( \x \) at which it fails.
:::

::: {.proof}
The category probabilities are \( \pi_{ir}=\gamma_{ir}-\gamma_{i,r-1} \), so they are all
positive exactly when \( \gamma_{i1}<\dots<\gamma_{i,c-1} \), and they automatically sum to
one. Since \( G \) is strictly increasing, this is @eq-mlt-ordering. With common slopes
the terms \( \x_{(i)}\T\bbeta \) cancel from every inequality. If
\( \mathbf{d}=\bbeta_{r+1}-\bbeta_r\neq\bzero \), the \( r \)th inequality reads
\( \x\T\mathbf{d}<\theta_{r+1}-\theta_r \), a half-space; any \( \x \) far enough in the
direction \( \mathbf{d} \) violates it.
:::

::: {.warning}
[What the constraint costs]

The proportional odds model is safe: the constraint holds for every \( \x \) at once, and
by @exr-mlt-order-free it is satisfied automatically at the maximum likelihood estimate.
A model with cut-specific slopes is *not* a model over the whole of \( \Real^{p} \), only
on the region where @eq-mlt-ordering holds; one fits it, checks that the fitted
probabilities are positive wherever one intends to predict, and reports the check. For a concrete failure take \( c=3 \), one regressor,
\( \theta_1=0 \), \( \beta_1=1 \), \( \theta_2=1 \), \( \beta_2=2 \): at \( x=0 \) the linear predictors
are \( 0 \) and \( 1 \), in order; at \( x=2 \) they are \( -2 \) and \( -3 \), out of order, and the
fitted probability of the middle category is negative.
:::

## Fitting

The cumulative link model is not a generalized linear model in the strict sense of
@def-glm-model: the response is a vector and the parameters are shared across its
\( c-1 \) coordinates in a way no single link expresses. The likelihood machinery survives
intact, because the response is still multinomial and the mean is still a smooth
function of the parameters.

::: {#prp-mlt-cumulative-score}
[Score and expected information for a smooth mean]

Let \( \Y_i\sim\text{Mult}(m_i,\boldsymbol{\uppi}_i) \) independently, \( i=1,\dots,n \), with
\( \boldsymbol{\uppi}_i=\boldsymbol{\uppi}_i(\boldsymbol{\uppsi}) \) a smooth function of a parameter
\( \boldsymbol{\uppsi}\in\Real^{q} \) and every \( \pi_{ir}>0 \). Write
\( \tilde{\boldsymbol{\uppi}}_i \) for the first \( c-1 \) probabilities,
\( \bD_i=\partial\tilde{\boldsymbol{\uppi}}_i/\partial\boldsymbol{\uppsi}\T \) and
\( \V_i=\diag(\tilde{\boldsymbol{\uppi}}_i)-\tilde{\boldsymbol{\uppi}}_i\tilde{\boldsymbol{\uppi}}_i\T \).
Then the score and the expected information are
\[
\bu(\boldsymbol{\uppsi})=\sum_{i=1}^{n}\bD_i\T\V_i^{-1}\bigl(\tilde{\y}_i-m_i\tilde{\boldsymbol{\uppi}}_i\bigr),
\qquad
\mathbf{J}(\boldsymbol{\uppsi})=\sum_{i=1}^{n}m_i\,\bD_i\T\V_i^{-1}\bD_i .
\]{#eq-mlt-score-general}

:::

::: {.proof}
The inverse of \( \V_i \) is \( \diag(1/\tilde{\boldsymbol{\uppi}}_i)+\pi_{ic}^{-1}\bone\bone\T \),
as multiplying out and using \( \bone\T\tilde{\boldsymbol{\uppi}}_i=1-\pi_{ic} \) shows. Also
\( \partial\ell/\partial\boldsymbol{\uppsi}=\sum_i\sum_{r\le c}(y_{ir}/\pi_{ir})\,\partial\pi_{ir}/\partial\boldsymbol{\uppsi} \).
Fix \( i \), drop it from the notation, and substitute the inverse:
\[
\bD\T\V^{-1}(\tilde{\y}-m\tilde{\boldsymbol{\uppi}})
=\sum_{r<c}\frac{\partial\pi_r}{\partial\boldsymbol{\uppsi}}\frac{y_r-m\pi_r}{\pi_r}
 +\Bigl(\sum_{r<c}\frac{\partial\pi_r}{\partial\boldsymbol{\uppsi}}\Bigr)
  \frac{\sum_{s<c}(y_s-m\pi_s)}{\pi_c}.
\]
Because the probabilities sum to one, \( \sum_{r<c}\partial\pi_r/\partial\boldsymbol{\uppsi}
=-\partial\pi_c/\partial\boldsymbol{\uppsi} \), and because the counts sum to \( m \),
\( \sum_{s<c}(y_s-m\pi_s)=-(y_c-m\pi_c) \). The second term is therefore
\( (\partial\pi_c/\partial\boldsymbol{\uppsi})(y_c-m\pi_c)/\pi_c \), and the two terms combine into
\( \sum_{r\le c}(\partial\pi_r/\partial\boldsymbol{\uppsi})(y_r-m\pi_r)/\pi_r \), which equals
\( \sum_{r\le c}(y_r/\pi_r)\partial\pi_r/\partial\boldsymbol{\uppsi} \) since
\( \sum_r\partial\pi_r/\partial\boldsymbol{\uppsi}=\bzero \). Summing over \( i \) gives the score.
For the information, \( \Cov(\tilde{\Y}_i)=m_i\V_i \) by @lem-mlt-multinomial, and the
summands of the score are independent with mean zero, so
\( \Cov(\bu)=\sum_i\bD_i\T\V_i^{-1}(m_i\V_i)\V_i^{-1}\bD_i=\sum_i m_i\bD_i\T\V_i^{-1}\bD_i \).
:::

For @eq-mlt-cumulative the derivatives are immediate: writing
\( \eta_{ir}=\theta_r-\x_{(i)}\T\bbeta_r \) and \( g=G' \) for the density,
\[
\frac{\partial\gamma_{ir}}{\partial\theta_s}=g(\eta_{ir})\,\delta_{rs},
\qquad
\frac{\partial\gamma_{ir}}{\partial\bbeta_s}=-g(\eta_{ir})\,\x_{(i)}\,\delta_{rs},
\]
and \( \bD_i \) is obtained by differencing these rows. Fisher scoring
\( \boldsymbol{\uppsi}\leftarrow\boldsymbol{\uppsi}+\mathbf{J}^{-1}\bu \) converges quickly from ordered
starting cutpoints and zero slopes. Because @eq-mlt-score-general has the form
"derivative of the mean, weighted by the inverse covariance, times the residual", each
step is again a weighted least squares fit — the block-weighted version of @thm-glm-irls
and @prp-res-irls — with \( \V_i \) in the place of scalar weights.

::: {.remark}
[Why the fitting is well behaved]

Suppose the link density \( g=G' \) is log-concave, as it is for the logistic, the normal
and the extreme-value distribution behind the complementary log-log link. Then
\[
G(a)-G(b)=\int g(u)\,\bone\{b<u\le a\}\,du
\]
integrates over \( u \) a function of \( (u,a,b) \) that is log-concave on \( \Real^3 \), being
a log-concave function of \( u \) times the indicator of a convex set; a marginal of a
log-concave function is log-concave, so \( \log\{G(a)-G(b)\} \) is concave on \( \{b<a\} \).
The observation in category \( r \) contributes
\( \log\{G(\eta_{ir})-G(\eta_{i,r-1})\} \) with \( \eta_{ir} \) affine in the parameters, so
the log-likelihood is concave on the region where @eq-mlt-ordering holds at the observed
\( \x_{(i)} \) — an intersection of half-spaces, hence convex. That is Pratt (1981) and
Burridge (1981). Concavity excludes local maxima that are not global, which is what makes
scoring from any ordered start dependable; by itself it gives neither existence nor
uniqueness, and with cut-specific slopes the region it holds on is cut out by the data.
:::

::: {.remark}
[Asymptotics, and what they rest on]

This is not the canonical-link generalized linear model @thm-glm-asymptotics is stated
for, so that theorem cannot be quoted here. What holds is the standard theory of maximum
likelihood in a smooth parametric model: suppose the model is correctly specified and
identifiable, \( c \) and \( p \) are fixed, \( \boldsymbol{\uppsi}^0 \) is interior to the
region where @eq-mlt-ordering holds, \( G \) is twice continuously differentiable with
\( g>0 \), the regressors are bounded, and
\( n^{-1}\mathbf{J}(\boldsymbol{\uppsi})\to\mathbf{J}_\infty \) positive definite uniformly
near \( \boldsymbol{\uppsi}^0 \). Then with probability tending to one the likelihood
equations have a consistent solution, and
\( \sqrt{n}(\hat{\boldsymbol{\uppsi}}-\boldsymbol{\uppsi}^0) \) converges in distribution to
\( \Normal_q(\bzero,\mathbf{J}_\infty^{-1}) \), with the three statistics of
@prp-glm-three-tests behaving as they do there.

We do not prove it; treat it as imported. Its first two steps are those of
@thm-glm-asymptotics — a Lyapunov limit theorem for the score, which
@prp-mlt-cumulative-score exhibits as a sum of independent mean-zero terms, and a uniform
quadratic expansion on compacts — but not its third, which used the global concavity of a
canonical-link generalized linear model. For the classical treatment see van der Vaart
(1998, chapter 5) or Lehmann and Casella (1998, chapter 6). Nothing here is exact in
finite samples.
:::

## Where people place themselves

::: {#exm-mlt-selfplacement}
[A seven-point scale]

The 1996 American National Election Study asked 944 respondents to place
themselves on a seven-point scale from extremely liberal (\( 1 \)) to extremely
conservative (\( 7 \)); the counts are 16, 103,
147, 256, 170,
218 and 34. With the respondent's age in decades, a
seven-point education scale and a twenty-four-point household income scale as
regressors, Fisher scoring for the common-slope cumulative logit model converges in
11 steps to

| | age (decades) | education | income |
|---|---|---|---|
| \( \hat{\beta} \) | 0.0942 | \( -0.1249 \) | 0.0276 |
| standard error | 0.0357 | 0.0395 | 0.0104 |
| \( z \) | 2.637 | \( -3.166 \) | 2.641 |

with cutpoints \( \hat{\theta}=(-3.7744,-1.6378,-0.6159,
0.5580,1.3625,3.6442) \), in order as
@prp-mlt-ordering requires, and maximized log-likelihood
\( -1613.8558 \). Older and richer respondents place themselves
further to the right, more educated respondents further to the left: each decade of age
multiplies the odds of a higher placement by \( 1.0988 \), each step up
the education scale by \( 0.8826 \). Whether one should believe a single
common slope is the question of [Section 36.3](03-proportional-odds.html), and for these
data the answer turns out to be no.
:::

[Figure 36.2.1](02-cumulative-logit-models.html#fig-mlt-ordinal-fit) shows what the fit
says. Panel (a) plots the six fitted cumulative logits against age: exactly parallel
lines, which is what a common slope means. Panel (b) differences them into the seven
category probabilities, whose bands shift smoothly to the right as age increases and can
never cross — the ordering constraint in picture form.

::: {when-format="html"}
![**Figure 36.2.1.** The fitted cumulative logit model for self-placement, as a function
of age at education 4 and income 17. (a) The six cumulative logits are parallel,
because the model gives them a common slope. (b) The same fit as category
probabilities.](ordinal_fit.svg){#fig-mlt-ordinal-fit width=100%}
:::

::: {when-format="pdf"}
![The fitted cumulative logit model for self-placement, as a function of age at
education 4 and income 17. (a) The six cumulative logits are parallel, because the
model gives them a common slope. (b) The same fit as category
probabilities.](ordinal_fit.pdf){width=100%}
:::

```{.python .run #cell-ordinal-data}
import numpy as np
import statsmodels.api as sm
from scipy import stats

NAMES = ["age (decades)", "education", "income"]

anes = sm.datasets.anes96.load_pandas().data
y = anes["selfLR"].astype(int).to_numpy()            # 1 = extremely liberal ... 7 = extremely conservative
X = np.column_stack([anes["age"] / 10.0, anes["educ"], anes["income"]])
n, p = X.shape
c = 7
Y = np.zeros((n, c))
Y[np.arange(n), y - 1] = 1.0

print("category counts:", Y.sum(axis=0))
```

The fitting routine below is @eq-mlt-score-general written out, with any subset of the
regressors allowed cut-specific slopes through the argument `free`. Empty, it is the
proportional odds fit; the other settings are for
[Section 36.3](03-proportional-odds.html).

```{.python .run #cell-ordinal-fit}
def cum_gamma(psi, X, c, free):
    """Cumulative probabilities Lambda(theta_r - x^T beta_r) and their gradients.

    psi = [theta (c-1), beta (p), delta ((c-2) x len(free))]; the columns of X listed
    in `free` get cut-specific slopes beta_r = beta + delta_r, with delta_1 = 0.
    """
    n, p = X.shape
    q, nf, fl = len(psi), len(free), list(free)
    theta, beta = psi[:c - 1], psi[c - 1:c - 1 + p]
    delta = psi[c - 1 + p:].reshape(c - 2, nf)
    G = np.empty((n, c - 1))
    dG = np.zeros((n, c - 1, q))
    for r in range(c - 1):
        beta_r = beta.copy()
        if r >= 1 and nf:
            beta_r[fl] = beta_r[fl] + delta[r - 1]
        g = 1.0 / (1.0 + np.exp(-(theta[r] - X @ beta_r)))
        lam = (g * (1.0 - g))[:, None]              # the logistic density at the cut
        G[:, r] = g
        dG[:, r, r] = lam[:, 0]
        dG[:, r, c - 1:c - 1 + p] = -lam * X
        if r >= 1 and nf:
            off = c - 1 + p + (r - 1) * nf
            dG[:, r, off:off + nf] = -lam * X[:, fl]
    return G, dG


def cum_probs(psi, X, c, free):
    """Category probabilities pi_ir and their derivatives D[i, r, k] = d pi_ir / d psi_k."""
    G, dG = cum_gamma(psi, X, c, free)
    P = np.diff(np.column_stack([np.zeros(len(X)), G, np.ones(len(X))]), axis=1)
    D = np.concatenate([dG[:, :1, :], np.diff(dG, axis=1)], axis=1)
    return P, D


def cum_score(psi, X, Y, c, free):
    """Score and expected information of the multinomial likelihood."""
    P, D = cum_probs(psi, X, c, free)
    Pm = P[:, :c - 1]
    Sigma = np.einsum("ir,rs->irs", Pm, np.eye(c - 1)) - np.einsum("ir,is->irs", Pm, Pm)
    A = np.linalg.solve(Sigma, D)                   # Sigma_i^{-1} D_i, one observation at a time
    U = np.einsum("irk,ir->k", A, Y[:, :c - 1] - Pm)
    J = np.einsum("irk,irl->kl", D, A)
    return U, J


def cum_loglik(psi, X, Y, c, free):
    P, _ = cum_probs(psi, X, c, free)
    return float(np.sum(Y * np.log(np.clip(P, 1e-300, None))))


def cum_fit(X, Y, c, free=(), tol=1e-10, maxit=100):
    """Fisher scoring for the cumulative logit model."""
    n, p = X.shape
    psi = np.zeros((c - 1) + p + (c - 2) * len(free))
    psi[:c - 1] = np.linspace(-1.5, 1.5, c - 1)     # start from ordered cutpoints, no slopes
    for it in range(1, maxit + 1):
        U, J = cum_score(psi, X, Y, c, free)
        step = np.linalg.solve(J, U)
        psi = psi + step
        if np.max(np.abs(step)) < tol:
            break
    return psi, J, it


psi, J, iters = cum_fit(X, Y, c)
theta, beta = psi[:c - 1], psi[c - 1:]
se = np.sqrt(np.diag(np.linalg.inv(J)))
ll_po = cum_loglik(psi, X, Y, c, ())

print("cutpoints", np.round(theta, 4))
for j, name in enumerate(NAMES):
    print(f"{name:>16s}  {beta[j]:8.4f}  ({se[c - 1 + j]:.4f})   z = {beta[j] / se[c - 1 + j]:6.3f}")
print(f"log-likelihood {ll_po:.4f} after {iters} Fisher scoring steps")
```

The script checks the fit, cutpoint by cutpoint, against the `OrderedModel` routine of
statsmodels.

## Choosing the link

The three links of @def-mlt-cumulative are close over the middle of the probability
range and differ in the tails, as in the binary case (@prp-bin-links). Two things are
specific to ordinal data. First, the *direction* of the scale is arbitrary. Relabelling
the categories in reverse order maps the logit and probit models to themselves with
\( \bbeta \) and the cutpoints negated (@exr-mlt-reverse), but turns the complementary
log-log model into a different one, the log-log model; with an asymmetric link the
analysis depends on which end of the scale was called "1", and the choice must be
justified.

Second, the complementary log-log link has an interpretation the others lack. With
\( G(z)=1-\exp(-e^{z}) \) and common slopes, @eq-mlt-cumulative gives
\[
P(Y_i>r\mid\x_{(i)})=\exp\bigl\{-e^{\theta_r}e^{-\x_{(i)}\T\bbeta}\bigr\}
=\bigl[P(Y>r\mid\bzero)\bigr]^{\exp(-\x_{(i)}\T\bbeta)} ,
\]{#eq-mlt-grouped-ph}

so the regressors act multiplicatively on the cumulative hazard of the underlying scale:
the grouped-data version of the proportional hazards model of Cox (1972). If the ordered
categories are intervals of a survival time, reach for it.

## Exercises

### A. Check your understanding

::: {#exr-mlt-reverse}
[A1]

Verify the claim that reversing the category order maps the cumulative logit model to
itself with \( \bbeta \) and the cutpoints negated, and find what it does to the
complementary log-log model.
:::

::: {.solution}
Let \( Y'=c+1-Y \). Then \( \{Y'\le r\}=\{Y\ge c+1-r\}=\{Y>c-r\} \), so
\( \gamma'_{ir}=1-\gamma_{i,c-r} \). Since \( \Lambda^{-1}(1-u)=-\Lambda^{-1}(u) \),
\( \Lambda^{-1}(\gamma'_{ir})=-\theta_{c-r}+\x\T\bbeta \), which is @eq-mlt-cumulative with
cutpoints \( -\theta_{c-r} \) (still increasing in \( r \)) and slope \( -\bbeta \). The same works
for \( \Phi \). For \( G(z)=1-\exp(-e^{z}) \) the identity \( G^{-1}(1-u)=-G^{-1}(u) \) fails,
and \( G^{-1}(1-u)=\log(-\log u) \) defines the log-log link, a different family.
:::

### B. Practice

::: {#exr-mlt-order-free}
[B1]

Show that in the proportional odds model the maximum likelihood estimate automatically
satisfies \( \hat{\theta}_1<\dots<\hat{\theta}_{c-1} \), provided every category is observed,
so the ordering constraint never binds. (Hint: what is the likelihood if two cutpoints
coincide or cross?)
:::

::: {.solution}
The likelihood is defined only where @eq-mlt-ordering holds, and
\( \pi_{ir}\to 0 \) for some \( i \) and \( r \) as the boundary \( \theta_r=\theta_{r+1} \) is
approached. If category \( r \) is observed at least once, some \( y_{ir}>0 \) and
\( \ell\to-\infty \) there. Since \( \ell \) is continuous on the open region and tends to
\( -\infty \) at the part of its boundary where two cutpoints meet, any maximizer lies in
the interior, where the cutpoints are strictly ordered.
:::

::: {#exr-mlt-latent-scale}
[B2]

The cumulative logit model is unchanged if the latent scale is multiplied by a
positive constant \( a \) and shifted by \( b \), provided \( \theta_r\mapsto a\theta_r+b \) and
\( \bbeta\mapsto a\bbeta \). Deduce that only *ratios* of coefficients, and the standardized
distances between cutpoints, are free of the arbitrary scale — and that a coefficient
from a logit fit is not comparable with a coefficient from a probit fit without
rescaling by roughly \( 1.8 \), the standard deviation of the logistic distribution.
:::

::: {#exr-mlt-collapse-fit}
[B3]

Refit @exm-mlt-selfplacement after collapsing the seven categories into three
(\( 1 \)–\( 3 \), \( 4 \), \( 5 \)–\( 7 \)) and again after dichotomizing at the middle category. Compare
the three sets of slopes and their standard errors, and say which comparison is a
check on the model and which is a statement about information.
:::

::: {.solution}
The estimates from the three fits are \( (0.0942,-0.1249,0.0276) \),
\( (0.1350,-0.1148,0.0343) \) and
\( (0.1229,-0.0311,0.0398) \), with standard errors that grow as
categories are merged. That they estimate the *same* parameters is the collapsing
invariance proved in [Section 36.3](03-proportional-odds.html) — a check on the model,
since it holds only if the proportional odds model is true. That the standard
errors grow is a statement about information: collapsing throws data away.
:::

### C. Going deeper

::: {#exr-mlt-interval-censored}
[C1]

Suppose an underlying continuous response \( U_i \) follows the linear model
\( U_i=\x_{(i)}\T\bbeta+\varepsilon_i \) with \( \varepsilon_i \) having distribution function
\( G \), but \( U_i \) is recorded only as the interval \( (\theta_{r-1},\theta_r] \) it falls in,
with *known* cutpoints. Write down the log-likelihood, show it is the cumulative link
model with the cutpoints held fixed, and say what is gained by knowing
them (@exr-mlt-latent-scale).
:::
