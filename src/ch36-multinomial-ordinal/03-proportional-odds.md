# Proportional odds and its checks

Setting \( \bbeta_1=\dots=\bbeta_{c-1}=\bbeta \) in @eq-mlt-cumulative turns \( (c-1)p \)
slopes into \( p \). With the logistic link the result is the model that carries the name
of its most conspicuous consequence.

::: {#thm-mlt-proportional-odds}
[The proportional odds model]

Let
\[
\log\frac{P(Y_i\le r\mid\x_{(i)})}{P(Y_i>r\mid\x_{(i)})}=\theta_r-\x_{(i)}\T\bbeta,
\qquad r=1,\dots,c-1 .
\]{#eq-mlt-po}

Then:

::: {.enumerate options="label=(\alph*)"}
1. *(latent variable)* @eq-mlt-po holds exactly when there is a latent response
      \( U_i=\x_{(i)}\T\bbeta+\varepsilon_i \), with \( \varepsilon_i \) independent of
      \( \x_{(i)} \) and standard logistic, and \( Y_i=r \) precisely when
      \( \theta_{r-1}<U_i\le\theta_r \), where
      \( -\infty=\theta_0<\theta_1<\dots<\theta_{c-1}<\theta_c=\infty \). Replacing the logistic
      distribution by any continuous, strictly increasing \( G \) gives the general cumulative
      link model of @def-mlt-cumulative with common slopes;

2. *(collapsing)* let \( 0=a_0<a_1<\dots<a_k=c \) be integers and define the coarsened
      response \( Y_i^{*}=j \) when \( a_{j-1}<Y_i\le a_j \). Then \( Y_i^{*} \) follows the
      proportional odds model with the *same* \( \bbeta \) and cutpoints
      \( \theta^{*}_j=\theta_{a_j} \), \( j=1,\dots,k-1 \). In particular dichotomizing at any cut
      \( r \) gives a binary logistic model with slope \( \bbeta \);

3. *(coherence)* the ordering constraint @eq-mlt-ordering reduces to
      \( \theta_1<\dots<\theta_{c-1} \) and so holds at every \( \x \) at once;

4. *(interpretation)* for any two regressor vectors and every \( r \),
      \[
      \frac{\text{odds}(Y>r\mid\x_2)}{\text{odds}(Y>r\mid\x_1)}
        =\exp\{(\x_2-\x_1)\T\bbeta\} ,
      \]{#eq-mlt-common-or}

      a single odds ratio that does not depend on \( r \); and if
      \( (\x_2-\x_1)\T\bbeta>0 \) then the response at \( \x_2 \) is stochastically larger than
      the response at \( \x_1 \).
:::

:::

::: {.proof}
(a) If such a latent variable exists, then
\( P(Y_i\le r)=P(U_i\le\theta_r)=P(\varepsilon_i\le\theta_r-\x_{(i)}\T\bbeta)
=\Lambda(\theta_r-\x_{(i)}\T\bbeta) \), which is @eq-mlt-po after taking logits, since
\( \log\{\Lambda(z)/(1-\Lambda(z))\}=z \). Conversely, given @eq-mlt-po, define \( U_i \) by the
displayed construction with \( \varepsilon_i \) standard logistic; the computation just made
shows that the induced distribution of \( Y_i \) is the one specified by the model, and the
cutpoints are increasing by (c).

(b) The events are nested: \( \{Y^{*}_i\le j\}=\{Y_i\le a_j\} \), directly from the
definition. Hence \( P(Y^{*}_i\le j\mid\x_{(i)})=\gamma_{i,a_j} \), and applying the logit
to both sides,
\( \log\{\gamma_{i,a_j}/(1-\gamma_{i,a_j})\}=\theta_{a_j}-\x_{(i)}\T\bbeta \). This is
@eq-mlt-po for \( Y^{*} \) with cutpoints \( \theta_{a_j} \), which are increasing because the
\( \theta_r \) are, and with \( \bbeta \) untouched. Taking \( k=2 \), \( a_1=r \) gives the binary
statement.

(c) Substitute \( \bbeta_r=\bbeta \) into @eq-mlt-ordering: the terms
\( \x_{(i)}\T\bbeta \) cancel, leaving \( \theta_1<\dots<\theta_{c-1} \), a condition on the
parameters alone. This is the case already recorded in @prp-mlt-ordering.

(d) From @eq-mlt-po, \( \text{odds}(Y\le r\mid\x)=\exp(\theta_r-\x\T\bbeta) \), so
\( \text{odds}(Y>r\mid\x)=\exp(\x\T\bbeta-\theta_r) \) and the ratio is
\( \exp\{(\x_2-\x_1)\T\bbeta\} \); the cutpoint cancels, which is the whole content
of @eq-mlt-common-or. If \( (\x_2-\x_1)\T\bbeta>0 \) the ratio exceeds one for every \( r \), so
\( P(Y>r\mid\x_2)>P(Y>r\mid\x_1) \) for every \( r \): stochastic dominance.
:::

::: {.idea}
Part (b) is why the model is worth having: the slope means the same thing however finely
the scale was cut, so studies using three categories and seven are comparable. It is
also the model's exposed flank, since *every* dichotomization must give that slope.
:::

::: {.warning}
[What proportional odds does not say]

@eq-mlt-common-or is about *cumulative* odds. It does not say that a unit of \( x \) has
the same effect on every category probability: differencing
\( \pi_{ir}=G(\eta_{ir})-G(\eta_{i,r-1}) \) gives
\( \partial\pi_{ir}/\partial\x=\{g(\eta_{i,r-1})-g(\eta_{ir})\}\,\bbeta \), whose sign changes
as \( r \) moves across the distribution. Nor does it say that the latent variable of part
(a) exists in any physical sense: the construction is a *device for generating the
model*, always available once the model holds, and it carries no empirical content.
Claims that the latent scale is "really" attitude, utility or severity are claims about
the world, and reading \( \bbeta \) causally needs the whole apparatus of
[Chapter 25](../ch25-causal-interpretation/index.html).
:::

## Checking the assumption

The restriction is \( (c-2)p \) linear constraints on the unrestricted cumulative model, so
the three statistics of @prp-glm-three-tests are all available. The likelihood ratio and
Wald tests need the unrestricted model, which by @prp-mlt-ordering may not be fittable;
the score test needs only the restricted fit, which always is.

::: {#prp-mlt-checks}
[Testing and diagnosing the proportional odds restriction]

Let \( A\subseteq\{1,\dots,p\} \) index the regressors whose slopes are allowed to vary,
and embed @eq-mlt-po in the cumulative logit model
\( \text{logit}\,\gamma_{ir}=\theta_r-\x_{(i)}\T(\bbeta+\boldsymbol{\updelta}_r) \), where
\( \boldsymbol{\updelta}_1=\bzero \) and \( \boldsymbol{\updelta}_r \) is supported on \( A \) for
\( r\ge 2 \). Write \( \boldsymbol{\uppsi}=(\boldsymbol{\uptheta},\bbeta,\boldsymbol{\updelta}) \) and let
\( \hat{\boldsymbol{\uppsi}}_0 \) be the maximum likelihood estimate under
\( \boldsymbol{\updelta}=\bzero \). Then, under the regularity conditions recorded in
[Section 36.2](02-cumulative-logit-models.html) and if @eq-mlt-po holds:

::: {.enumerate options="label=(\alph*)"}
1. the score statistic
      \( S_{\mathrm{sc}}=\bu(\hat{\boldsymbol{\uppsi}}_0)\T\mathbf{J}(\hat{\boldsymbol{\uppsi}}_0)^{-1}
      \bu(\hat{\boldsymbol{\uppsi}}_0) \), with \( \bu \) and \( \mathbf{J} \) from
      @eq-mlt-score-general, is computable from the restricted fit alone, equals
      \( \bu_{\boldsymbol{\updelta}}\T[\mathbf{J}^{-1}]_{\boldsymbol{\updelta}\boldsymbol{\updelta}}
      \bu_{\boldsymbol{\updelta}} \), and converges in distribution to \( \chi^2 \) with
      \( (c-2)|A| \) degrees of freedom;

2. the likelihood ratio and Wald statistics for the same hypothesis have the same
      limit and are asymptotically equivalent to \( S_{\mathrm{sc}} \);

3. the \( c-1 \) separate binary logistic regressions of \( \bone\{Y_i\le r\} \) on
      \( \x_{(i)} \) estimate the same \( \bbeta \), so plotting their slopes with confidence
      intervals shows where any departure lies. The \( c-1 \) fits are dependent, so the plot
      is a diagnostic and not a test.
:::

:::

::: {.proof}
(a) By definition of a maximum, the components of \( \bu(\hat{\boldsymbol{\uppsi}}_0) \) in the
\( (\boldsymbol{\uptheta},\bbeta) \) block vanish: those are the likelihood equations of the
restricted model, which is the full model with \( \boldsymbol{\updelta} \) held at zero. Writing
\( \bu=(\bzero\T,\bu_{\boldsymbol{\updelta}}\T)\T \) and partitioning \( \mathbf{J}^{-1} \)
conformably gives the stated reduction. That
\( \bu\T\mathbf{J}^{-1}\bu \) has the \( \chi^2 \) limit with degrees of freedom the number of
restrictions is the score-test half of @prp-glm-three-tests, which for a generalized
linear model is proved there and holds here by the same argument applied to the
likelihood of @prp-mlt-cumulative-score, under the conditions imported in
[Section 36.2](02-cumulative-logit-models.html); the restrictions number \( (c-2)|A| \), one
vector \( \boldsymbol{\updelta}_r \) of length \( |A| \) for each of the cuts \( r=2,\dots,c-1 \).

(b) The asymptotic equivalence of the three statistics is @prp-glm-three-tests, and it
carries over under the same imported conditions; the hypothesis is a linear restriction
on a parameter interior to the parameter space.

(c) Dichotomizing at cut \( r \) gives, by @thm-mlt-proportional-odds(b), the binary
logistic model \( \text{logit}\,P(Y\le r)=\theta_r-\x\T\bbeta \), with intercept
\( \theta_r \) and the same \( \bbeta \). Each fit is therefore a correctly
specified binary logistic regression under the null and its estimate is consistent for
\( \bbeta \) by @thm-glm-asymptotics. The fits use overlapping data, so their errors are
correlated, which is why they cannot simply be pooled or differenced into a test.
:::

Brant (1990) turns (c) into a test by estimating the joint covariance of the \( c-1 \)
slope estimates — they solve stacked estimating equations, so a sandwich-type covariance
in the spirit of @thm-het-sandwich is available — and forming a Wald statistic for their
equality. It tests the same null hypothesis and has the same chi-squared limit. It is
*not* the score statistic, it is built from different estimating equations, and the two
can differ appreciably in small samples.

::: {.warning}
[The test is not the question]

With \( n \) in the thousands the score test rejects departures too small to matter; with
\( n \) in the dozens it misses departures that would change every conclusion. The plot in
(c) answers the useful question — *how much* do the slopes move, and where. Report both.
:::

## Self-placement, checked

::: {#exm-mlt-po-check}
[The proportional odds model does not hold]

The common odds ratio of @thm-mlt-proportional-odds(d) is worth having only if it is
there. For @exm-mlt-selfplacement, let all three slopes vary across the six cuts. The
score statistic of @prp-mlt-checks is 67.667 on 15 degrees
of freedom, \( p=1.16\times 10^{-8} \); the likelihood ratio
statistic against the fitted unrestricted model is 67.420, agreeing to
within half a percent as @prp-mlt-checks(b) predicts. One regressor at a time, on
\( c-2=5 \) degrees of freedom:

| regressor | score statistic | \( p \) | likelihood ratio |
|---|---|---|---|
| age | 17.542 | \( 3.58\times 10^{-3} \) | 18.367 |
| education | 41.258 | \( 8.32\times 10^{-8} \) | 40.190 |
| income | 20.056 | \( 1.22\times 10^{-3} \) | 20.155 |

Education is the main offender. [Figure 36.3.1](03-proportional-odds.html#fig-mlt-cuts)
shows why: its slope is near \( -0.3 \) at the three lowest cuts, near \( 0 \) at cuts
\( 4 \) and \( 5 \), and back near \( -0.3 \) at cut \( 6 \). Education sharply separates the
liberal end of the scale from everything above it, says almost nothing about position
within the upper middle, and then separates the extreme conservative category again.
That is a fact about the data, and a single number cannot express it.
:::

::: {when-format="html"}
![**Figure 36.3.1.** Slopes from the six separate binary logistic fits of
\( \bone\{Y\le r\} \) on the three regressors, with 95% intervals; the horizontal line and
band are the proportional odds estimate and its interval. Under proportional odds the
six points scatter about the line (@prp-mlt-checks). Education does
not.](ordinal_cuts.svg){#fig-mlt-cuts width=100%}
:::

::: {when-format="pdf"}
![Slopes from the six separate binary logistic fits of \( \bone\{Y\le r\} \) on the three
regressors, with 95% intervals; the horizontal line and band are the proportional odds
estimate and its interval. Under proportional odds the six points scatter about the
line (@prp-mlt-checks). Education does not.](ordinal_cuts.pdf){width=100%}
:::

The cells on this page continue those of
[Section 36.2](02-cumulative-logit-models.html); this first one repeats the data and the
fitting routine so that the page runs on its own.

```{.python .run #cell-ordinal-setup}
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


def cum_gamma(psi, X, c, free):
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
    G, dG = cum_gamma(psi, X, c, free)
    P = np.diff(np.column_stack([np.zeros(len(X)), G, np.ones(len(X))]), axis=1)
    D = np.concatenate([dG[:, :1, :], np.diff(dG, axis=1)], axis=1)
    return P, D


def cum_score(psi, X, Y, c, free):
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
print(f"log-likelihood {ll_po:.4f} after {iters} Fisher scoring steps")
```

```{.python .run #cell-ordinal-cuts}
# the same slopes estimated one cut at a time: c - 1 separate binary logistic fits
Xc = sm.add_constant(X)
cut_beta = np.empty((c - 1, p))
cut_se = np.empty((c - 1, p))
for r in range(1, c):
    fit_r = sm.Logit((y <= r).astype(float), Xc).fit(disp=0)
    cut_beta[r - 1] = -np.asarray(fit_r.params)[1:]    # logit P(Y <= r) = theta_r - x^T beta_r
    cut_se[r - 1] = np.asarray(fit_r.bse)[1:]

for r in range(c - 1):
    print(f"cut {r + 1}: " + "  ".join(f"{b:7.4f} ({s:.4f})" for b, s in zip(cut_beta[r], cut_se[r])))
```

The score test needs no second fit: evaluate the score and information of the *larger*
model at the restricted estimate, with the extra coefficients zero.

```{.python .run #cell-ordinal-score}
# score test of proportional odds: all three slopes free to vary across cuts
free = (0, 1, 2)
psi_free = np.zeros((c - 1) + p + (c - 2) * len(free))
psi_free[:c - 1 + p] = psi                          # the proportional odds fit, with delta = 0
U, J_free = cum_score(psi_free, X, Y, c, free)
score_stat = float(U @ np.linalg.solve(J_free, U))
df = (c - 2) * len(free)

psi_np, _, _ = cum_fit(X, Y, c, free=free)
ll_np = cum_loglik(psi_np, X, Y, c, free)
lr_stat = 2 * (ll_np - ll_po)

print(f"score {score_stat:.3f} on {df} df, p = {stats.chi2.sf(score_stat, df):.3g}")
print(f"LR    {lr_stat:.3f} on {df} df, p = {stats.chi2.sf(lr_stat, df):.3g}")
```

## When it fails

Rejecting proportional odds is information, not a disaster. The options, roughly in
order of how much they give up.

**Partial proportional odds.** Free the slopes of the offending regressors only, as in
@prp-mlt-checks, keeping the restriction for the rest: the model of Peterson and
Harrell (1990). It costs the coherence guarantee, so @prp-mlt-ordering must now be
checked at every \( \x \) of interest. Freeing education alone raises the log-likelihood to
\( -1593.7610 \), a likelihood ratio statistic of 40.190
on five degrees of freedom (\( p=1.37\times 10^{-7} \)), and lowers Akaike's
criterion (@def-sel-aic-bic) from 3245.71 to
3215.52; the fully unrestricted model reaches 3208.29
with ten more parameters.
The fitted education slopes by cut are \( -0.2683 \), \( -0.2988 \),
\( -0.2523 \), \( -0.0206 \), \( -0.0498 \) and
\( -0.2970 \), and every fitted probability stays positive.

```{.python .run #cell-ordinal-partial}
# partial proportional odds: education alone gets cut-specific slopes
psi_pp, J_pp, _ = cum_fit(X, Y, c, free=(1,))
ll_pp = cum_loglik(psi_pp, X, Y, c, (1,))
educ_slopes = np.concatenate([[psi_pp[c - 1 + 1]], psi_pp[c - 1 + 1] + psi_pp[c - 1 + p:]])

aic = lambda ll, k: -2 * ll + 2 * k
print("education slope by cut:", np.round(educ_slopes, 4))
print("AIC  proportional odds %.2f | partial %.2f | unrestricted %.2f"
      % (aic(ll_po, c - 1 + p), aic(ll_pp, len(psi_pp)), aic(ll_np, len(psi_np))))
```

**A different link.** Non-proportionality is sometimes link misspecification in
disguise: the three links disagree in the tails (@prp-bin-links), so data can satisfy a
cumulative probit or complementary log-log model with common slopes while violating the
logit version. Trying the other two costs nothing.

**A scale effect.** Let the latent variable of @thm-mlt-proportional-odds(a) have a
dispersion that depends on the regressors, \( \log\sigma(\x)=\x\T\bgamma \). The cumulative
logits become \( (\theta_r-\x\T\bbeta)e^{-\x\T\bgamma} \), non-parallel in a particular
way: they fan out rather than crossing arbitrarily (@exr-mlt-po-scale). This costs \( p \)
parameters instead of \( (c-2)p \), and when it fits, "the effect is the same but the
responses are more spread out here" is a far more useful statement than a table of
slopes.

**A different family.** Data that violate proportional odds badly may satisfy a
sequential model ([Section 36.4](04-sequential-models.html)) with common slopes; failing
that, the baseline-category model abandons the ordering but is never wrong about it.

**Report the cuts.** Sometimes no single number describes the effect, and the honest
answer is @prp-mlt-checks(c). What is *not* a remedy is keeping the proportional odds
fit and repairing its standard errors with a sandwich estimator: @thm-het-sandwich fixes
the variance of an estimator whose target is still well defined, and here the target
itself has evaporated.

## Exercises

### A. Check your understanding

::: {#exr-mlt-po-df}
[A1]

A proportional odds model has \( c=4 \) and \( p=3 \). How many parameters does it have?
How many does the model with cut-specific slopes have, and what are the degrees of
freedom of the score test of @prp-mlt-checks with \( A=\{1,2,3\} \) and with \( A=\{2\} \)?
:::

::: {#exr-mlt-po-interpret}
[A2]

In a proportional odds fit the coefficient of a treatment indicator is
\( \hat{\beta}=0.7 \). Write down three correct sentences and one common incorrect one.
:::

::: {.solution}
Correct: the odds that a treated subject exceeds any fixed level are \( e^{0.7}\approx 2 \)
times the odds for a control; the same factor applies at every level; the treated
response is stochastically larger. Incorrect: "treatment doubles the probability of each
category" — the effect on a category probability is a difference of two logistic
densities times \( \hat{\beta} \), and changes sign across the scale.
:::

### B. Practice

::: {#exr-mlt-collapse-numeric}
[B1]

Using @thm-mlt-proportional-odds(b), predict what happens to the estimated slopes when
the seven categories of @exm-mlt-selfplacement are collapsed to three, and compare with
the actual refit in @exr-mlt-collapse-fit. Does the comparison support or undermine the
proportional odds assumption here?
:::

::: {.solution}
Under the model the collapsed fit estimates the same \( \bbeta \), so the two should differ
only by sampling error, the collapsed one being less precise. The education slope moves
from \( -0.1249 \) to \( -0.1148 \), small
against its standard error of \( 0.0416 \), and the age slope from
\( 0.0942 \) to \( 0.1350 \), about one standard
error: mild evidence at worst. Collapsing to three categories hides the failure that the
cut-by-cut plot makes obvious, so a check should use all the cuts.
:::

::: {#exr-mlt-po-scale}
[B2]

Derive the cumulative logits of the scale-effect model
\( U_i=\x_{(i)}\T\bbeta+\sigma(\x_{(i)})\varepsilon_i \), \( \log\sigma(\x)=\x\T\bgamma \),
with \( \varepsilon_i \) standard logistic. Show that it contains @eq-mlt-po as the case
\( \bgamma=\bzero \), and that its cumulative logits never cross as functions of \( \x \).
:::

::: {.solution}
\( P(Y_i\le r)=P(\varepsilon_i\le(\theta_r-\x\T\bbeta)e^{-\x\T\bgamma}) \), so
\( \text{logit}\,\gamma_{ir}=(\theta_r-\x\T\bbeta)e^{-\x\T\bgamma} \), which is @eq-mlt-po
when \( \bgamma=\bzero \). For two cuts \( r<s \) the difference of the cumulative logits is
\( (\theta_r-\theta_s)e^{-\x\T\bgamma} \), which has the sign of \( \theta_r-\theta_s \) for every
\( \x \), since the exponential is positive: the curves fan out but never cross, so the
ordering constraint of @prp-mlt-ordering holds everywhere.
:::

### C. Going deeper

::: {#exr-mlt-po-monotone}
[C1]

Show that if the proportional odds model holds and \( \beta_j>0 \), then
\( \E\{h(Y)\mid\x\} \) is increasing in \( x_j \) for *every* nondecreasing function \( h \).
Deduce that the sign of \( \beta_j \) is a statement that survives any recoding of the
categories into numbers, as long as the recoding respects the order.
:::

::: {.solution}
By @thm-mlt-proportional-odds(d), increasing \( x_j \) makes the response stochastically
larger. Stochastic dominance is equivalent to \( \E\{h(Y)\mid\x_2\}\ge\E\{h(Y)\mid\x_1\} \)
for every nondecreasing \( h \): write
\( \E h(Y)=h(1)+\sum_{r=1}^{c-1}\{h(r+1)-h(r)\}P(Y>r) \), a sum of nonnegative
coefficients times the survivor probabilities, each of which increases. Since any
order-respecting scoring of the categories is such an \( h \), the sign of \( \beta_j \)
cannot be reversed by rescoring.
:::

::: {#exr-mlt-po-rank}
[C2]

Let @thm-mlt-proportional-odds hold with a single binary regressor, let \( Y_1 \) and
\( Y_2 \) be independent responses from the two groups, and let \( U_1,U_2 \) be the latent
variables of part (a) of that theorem.

::: {.enumerate options="label=(\alph*)"}
1. Show that \( P(U_2>U_1)=\E\{\Lambda(\beta-\varepsilon)\} \) with \( \varepsilon \) standard
     logistic, so that it is a function of \( \beta \) alone and strictly increasing in it.

2. Show that the concordance probability of the *observed* responses,
     \( P(Y_2>Y_1)+\tfrac12P(Y_2=Y_1) \), is increasing in \( \beta \) for fixed cutpoints.
     (Hint: write it as \( \E\{h(Y_2)\} \) for a nondecreasing \( h \) that depends only on the
     distribution of \( Y_1 \).) Then evaluate it at \( \beta=1 \) for \( c=3 \) with cutpoints
     \( (-0.5,0.5) \) and again with \( (-3,3) \), and conclude that, unlike the odds ratio
     \( e^{\beta} \), it is not free of the cutpoints.

3. The sample version of (b) is the Mann–Whitney statistic divided by \( n_1n_2 \), with
     ties counted half. Say in what sense the proportional odds model is a regression
     version of the Wilcoxon–Mann–Whitney test, and in what sense it is not.
:::

:::

::: {.solution}
(a) \( U_2-U_1=\beta+\varepsilon_2-\varepsilon_1 \), so
\( P(U_2>U_1)=P(\varepsilon_1<\beta+\varepsilon_2)=\E\{\Lambda(\beta+\varepsilon_2)\} \), which
is \( \E\{\Lambda(\beta-\varepsilon)\} \) by symmetry of the logistic law. No cutpoint
appears, and \( \Lambda \) is strictly increasing, so this increases strictly with
\( \beta \); at \( \beta=1 \) it is 0.6613.

(b) Conditioning on \( Y_1 \) gives \( \E\{h(Y_2)\} \) with
\( h(b)=P(Y_1<b)+\tfrac12P(Y_1=b) \), nondecreasing in \( b \) and depending only on the
distribution of \( Y_1 \) — which does not involve \( \beta \), group one being the
reference. Raising \( \beta \) makes \( Y_2 \) stochastically larger by
@thm-mlt-proportional-odds(d), and \( \E\{h(Y_2)\} \) then increases, by the equivalence
used in @exr-mlt-po-monotone. It is not a function of \( \beta \) alone: at \( \beta=1 \)
the cutpoints \( (-0.5,0.5) \) give 0.6370 and \( (-3,3) \) give
0.5482, while \( e^{\beta} \) is the same in both. Widely spaced cutpoints
push mass into the middle category, so ties are common and the concordance probability is
pulled towards \( \tfrac12 \).

(c) It is a regression version in that the two-group comparison is the same
stochastic-ordering comparison, extended to several regressors and to adjustment. It is
not, in two ways: the rank test is distribution free where the model imposes proportional
odds at every cut, a real restriction by @exm-mlt-po-check; and what the model reports is
\( \beta \), not the concordance probability, which by (b) cannot be recovered from
\( \beta \) without the cutpoints.
:::
