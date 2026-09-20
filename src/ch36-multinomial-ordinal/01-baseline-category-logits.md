# Baseline-category logits

Let the response of observation \( i \) fall in one of \( c \) categories, recorded as a
vector of counts \( \Y_i=(Y_{i1},\dots,Y_{ic})\T \) with \( \sum_r Y_{ir}=m_i \). When each row
is a single subject, \( m_i=1 \) and \( \Y_i \) is an indicator vector; when rows are covariate
patterns, \( m_i \) counts the subjects sharing that pattern. As with grouped and ungrouped
binary data in [Chapter 35](../ch35-binary-responses/index.html), the two cases are the
same model with different bookkeeping, and the distinction matters only for the deviance
as a goodness-of-fit statistic.

## The multinomial distribution is an exponential family

[Chapter 34](../ch34-exponential-families-glm/index.html) was built on a one-dimensional
exponential dispersion family (@def-glm-edf). A categorical response is a *vector* with
\( c-1 \) free coordinates, and the structure survives the move, with the natural parameter
and the cumulant function becoming vector-valued.

::: {#lem-mlt-multinomial}
[The multinomial as an exponential family]

Let \( \Y\sim\text{Mult}(m,\boldsymbol{\uppi}) \) with \( \pi_r>0 \) for every \( r \) and
\( \sum_{r=1}^{c}\pi_r=1 \). Write \( \tilde{\Y}=(Y_1,\dots,Y_{c-1})\T \) for the first \( c-1 \)
counts and
\[
\eta_r=\log\frac{\pi_r}{\pi_c},\qquad r=1,\dots,c-1 .
\]

Then the probability mass function is
\[
f(\y)=h(\y)\,\exp\bigl\{\tilde{\y}\T\boldsymbol{\upeta}-m\,b(\boldsymbol{\upeta})\bigr\},
\qquad
b(\boldsymbol{\upeta})=\log\Bigl(1+\sum_{s=1}^{c-1}e^{\eta_s}\Bigr),
\]{#eq-mlt-expfam}

with \( h(\y)=m!/(y_1!\cdots y_c!) \). The map \( \boldsymbol{\upeta}\mapsto\boldsymbol{\uppi} \) is a
bijection from \( \Real^{c-1} \) onto the interior of the probability simplex, with inverse
\[
\pi_r=\frac{e^{\eta_r}}{1+\sum_{s<c}e^{\eta_s}}\ (r<c),
\qquad
\pi_c=\frac{1}{1+\sum_{s<c}e^{\eta_s}},
\]{#eq-mlt-inverse-link}

and the first two moments are
\[
\E(\tilde{\Y})=m\,\nabla b(\boldsymbol{\upeta})=m\tilde{\boldsymbol{\uppi}},
\qquad
\Cov(\tilde{\Y})=m\,\nabla^2 b(\boldsymbol{\upeta})
 =m\bigl(\diag(\tilde{\boldsymbol{\uppi}})-\tilde{\boldsymbol{\uppi}}\tilde{\boldsymbol{\uppi}}\T\bigr),
\]{#eq-mlt-moments}

where \( \tilde{\boldsymbol{\uppi}}=(\pi_1,\dots,\pi_{c-1})\T \).
:::

::: {.proof}
Since \( \sum_r y_r=m \),
\[
\prod_{r=1}^{c}\pi_r^{y_r}
=\pi_c^{m}\prod_{r<c}\Bigl(\frac{\pi_r}{\pi_c}\Bigr)^{y_r}
=\exp\Bigl\{\sum_{r<c}y_r\eta_r+m\log\pi_c\Bigr\},
\]
and \( \log\pi_c=-b(\boldsymbol{\upeta}) \) by the definition of \( \eta_r \), which
gives @eq-mlt-expfam. Solving \( \eta_r=\log\pi_r-\log\pi_c \) with
\( \sum_r\pi_r=1 \) gives @eq-mlt-inverse-link, and every \( \boldsymbol{\uppi} \) with positive entries arises from
exactly one \( \boldsymbol{\upeta} \), so the map is a bijection.

For the moments, sum @eq-mlt-expfam over the support. By the multinomial theorem,
\[
M(\boldsymbol{\upeta})=\sum_{\y}h(\y)\exp(\tilde{\y}\T\boldsymbol{\upeta})
=\Bigl(1+\sum_{s<c}e^{\eta_s}\Bigr)^{m}=e^{m b(\boldsymbol{\upeta})},
\]
a finite sum, which may therefore be differentiated term by term without any regularity
condition. Since \( \log M=m\,b \) is the cumulant generating function of \( \tilde{\Y} \)
shifted to \( \boldsymbol{\upeta} \), its first two derivatives are
\( \E(\tilde{\Y})=m\nabla b \) and \( \Cov(\tilde{\Y})=m\nabla^2 b \). Finally
\[
\frac{\partial b}{\partial \eta_r}=\frac{e^{\eta_r}}{1+\sum_{s<c}e^{\eta_s}}=\pi_r,
\qquad
\frac{\partial^2 b}{\partial \eta_r\,\partial\eta_s}=\pi_r(\delta_{rs}-\pi_s),
\]
where \( \delta_{rs} \) is one if \( r=s \) and zero otherwise.
:::

::: {.remark}
[What is new and what is not]

@eq-mlt-expfam has the shape of @def-glm-edf with the scalar natural parameter replaced
by a vector and the dispersion \( \phi/w_i \) by \( 1/m \): the number of trials plays the role
of a prior weight. The identity "mean is the gradient of the cumulant function, variance
its second derivative" is @prp-glm-moments, proved there by the same argument. New is
only that \( \nabla^2 b \) is a matrix, and a singular one once the baseline coordinate is
put back: \( \diag(\boldsymbol{\uppi})-\boldsymbol{\uppi}\boldsymbol{\uppi}\T \) kills \( \bone \),
because the \( c \) counts sum to the fixed total \( m \).
:::

## The model

::: {#def-mlt-baseline}
[Baseline-category logit model]

Let \( \Y_1,\dots,\Y_n \) be independent with \( \Y_i\sim\text{Mult}(m_i,\boldsymbol{\uppi}_i) \), and
let \( \x_{(i)}\in\Real^{p} \) be the \( i \)th row of a model matrix \( \X \) whose first column is
\( \bone \). The **baseline-category logit model** with baseline \( c \) is
\[
\log\frac{\pi_{ir}}{\pi_{ic}}=\x_{(i)}\T\bbeta_r,\qquad r=1,\dots,c-1,
\]{#eq-mlt-baseline}

with unknown \( \bbeta_1,\dots,\bbeta_{c-1}\in\Real^{p} \). Equivalently, writing
\( \bbeta_c=\bzero \),
\( \pi_{ir}=\exp(\x_{(i)}\T\bbeta_r)/\sum_{s=1}^{c}\exp(\x_{(i)}\T\bbeta_s) \) for every \( r \).
:::

This is a generalized linear model in the sense of @def-glm-model, with the multinomial
random component of @lem-mlt-multinomial and the *canonical* link, since
@eq-mlt-baseline sets \( \eta_{ir} \) itself equal to a linear predictor. The one departure
from [Chapter 34](../ch34-exponential-families-glm/index.html) is that each observation
carries a vector of natural parameters, so the model's coefficient vector is the stack
\( \bbeta=(\bbeta_1\T,\dots,\bbeta_{c-1}\T)\T \) of length \( (c-1)p \). Every regressor buys
\( c-1 \) coefficients: the price of treating the categories as unordered names.

## The likelihood equations

::: {#thm-mlt-score}
[Score, information and concavity]

Under @def-mlt-baseline the log-likelihood, as a function of
\( \bbeta=(\bbeta_1\T,\dots,\bbeta_{c-1}\T)\T \), is
\[
\ell(\bbeta)=\sum_{i=1}^{n}\Bigl\{\sum_{r<c}y_{ir}\,\x_{(i)}\T\bbeta_r
 -m_i\,b(\boldsymbol{\upeta}_i)\Bigr\},\qquad \eta_{ir}=\x_{(i)}\T\bbeta_r .
\]

Then:

::: {.enumerate options="label=(\alph*)"}
1. the score has blocks
      \( \partial\ell/\partial\bbeta_r=\sum_{i}(y_{ir}-m_i\pi_{ir})\,\x_{(i)}
         =\X\T(\y_r-\bmu_r) \), where \( \bmu_r \) is the vector with entries \( m_i\pi_{ir} \);

2. the observed and expected information coincide and have blocks
      \( \X\T\W_{rs}\X \) with \( \W_{rs}=\diag\{m_i\pi_{ir}(\delta_{rs}-\pi_{is})\}_{i=1}^n \);

3. \( \ell \) is concave, and strictly concave when \( \rank(\X)=p \); in that case it has at
      most one maximizer, and the maximizer, if it exists, is the unique solution of the
      likelihood equations \( \X\T(\y_r-\bmu_r)=\bzero \) for \( r=1,\dots,c-1 \).
:::

:::

::: {.proof}
(a) Differentiate: \( \partial b(\boldsymbol{\upeta}_i)/\partial\eta_{ir}=\pi_{ir} \) by
@lem-mlt-multinomial, and \( \partial\eta_{ir}/\partial\bbeta_r=\x_{(i)} \).

(b) Differentiating again, and using
\( \partial^2 b/\partial\eta_{ir}\partial\eta_{is}=\pi_{ir}(\delta_{rs}-\pi_{is}) \),
\[
-\frac{\partial^2\ell}{\partial\bbeta_r\,\partial\bbeta_s\T}
=\sum_i m_i\,\pi_{ir}(\delta_{rs}-\pi_{is})\,\x_{(i)}\x_{(i)}\T .
\]
The data enter only through the score, so the second derivative is nonrandom and
observed information equals expected information. This is the canonical-link
simplification of @thm-glm-score.

(c) Take \( \mathbf{a}=(\mathbf{a}_1\T,\dots,\mathbf{a}_{c-1}\T)\T \) and put
\( t_{ir}=\x_{(i)}\T\mathbf{a}_r \) for \( r<c \) and \( t_{ic}=0 \). By (b),
\[
\mathbf{a}\T\bigl(-\nabla^2\ell\bigr)\mathbf{a}
=\sum_i m_i\Bigl\{\sum_{r=1}^{c}\pi_{ir}t_{ir}^2-\Bigl(\sum_{r=1}^{c}\pi_{ir}t_{ir}\Bigr)^{2}\Bigr\},
\]
the inner brace being the variance of the random variable that takes the value \( t_{ir} \)
with probability \( \pi_{ir} \). A variance is nonnegative, so \( \ell \) is concave. It
vanishes only if \( t_{i1}=\dots=t_{ic} \), and since \( t_{ic}=0 \) this forces
\( \x_{(i)}\T\mathbf{a}_r=0 \) for every \( i \) and every \( r<c \), that is \( \X\mathbf{a}_r=\bzero \).
When \( \rank(\X)=p \) this gives \( \mathbf{a}=\bzero \), so the Hessian is negative definite
everywhere and \( \ell \) is strictly concave. A strictly concave function has at most one
stationary point, and a stationary point is its maximum.
:::

::: {.idea}
The likelihood equations say \( \X\T\y_r=\X\T\hat{\bmu}_r \) for each non-baseline category:
*fitted counts reproduce observed counts in every direction spanned by a regressor*, so
with an intercept each category's fitted total is its observed total. This is the
canonical-link moment matching of @thm-glm-score, one copy per logit.
:::

Because the information is the expected information, Newton–Raphson and Fisher scoring
are the same algorithm, and by @thm-glm-irls each step is a weighted least squares fit,
with block weights
\( \W_i=m_i(\diag(\tilde{\boldsymbol{\uppi}}_i)-\tilde{\boldsymbol{\uppi}}_i\tilde{\boldsymbol{\uppi}}_i\T) \)
and working response
\( \boldsymbol{\upeta}_i+\W_i^{-1}(\y_i-m_i\boldsymbol{\uppi}_i)_{1:c-1} \). The only change from
[Chapter 34](../ch34-exponential-families-glm/index.html) is that the weights are
\( (c-1)\times(c-1) \) blocks rather than scalars, so each step is a generalized rather
than a diagonally weighted least squares problem (@thm-het-wls).

::: {.warning}
[Existence]

Concavity guarantees uniqueness, not existence. As for binary data, the maximum is
attained unless the categories can be separated. Take \( m_i=1 \) and \( \rank(\X)=p \) as in
@thm-mlt-score(c), with observation \( i \) in
category \( y_i \). If there are vectors \( \mathbf{a}_1,\dots,\mathbf{a}_c \), not all equal, with
\( \x_{(i)}\T(\mathbf{a}_{y_i}-\mathbf{a}_r)\ge 0 \) for every \( i \) and \( r \), then differentiating
\( \ell \) along that direction gives
\( \sum_i\sum_r\pi_{ir}\,\x_{(i)}\T(\mathbf{a}_{y_i}-\mathbf{a}_r)\ge 0 \) at *every* point of
the ray, the fitted probabilities being positive throughout. By @thm-mlt-score(c) the
restriction of \( \ell \) to the ray is strictly concave, so its derivative is strictly
decreasing; were it zero somewhere it would be negative beyond, contradicting the
display. It is therefore positive throughout, and the supremum is approached only at
infinity. (Without the rank condition the derivative can vanish along the whole ray,
which forces every \( \x_{(i)} \) orthogonal to every difference
\( \mathbf{a}_s-\mathbf{a}_t \); the maximum is then attained, but not unique.) The binary case,
with its proof and its remedies, is @thm-bin-separation; the multinomial statement is
Albert and Anderson (1984). With many categories and small cells separation is common,
and the remedies are the same: a penalty as in @def-shr-ridge, or the equivalent proper
prior of @prp-reg-map — see [Chapter 27](../ch27-shrinkage/index.html) and
[Chapter 30](../ch30-regularization-boosting/index.html).
:::

## Reading the coefficients

Two facts about @eq-mlt-baseline are easy to state and easy to misuse.

::: {#prp-mlt-invariance}
[Changing the baseline]

Let \( \hbeta_1,\dots,\hbeta_{c-1} \) maximize the likelihood of @def-mlt-baseline, and set
\( \hbeta_c=\bzero \). Fix any \( s\in\{1,\dots,c\} \). Then:

::: {.enumerate options="label=(\alph*)"}
1. for every pair of categories \( r,r' \),
     \( \log(\pi_{ir}/\pi_{ir'})=\x_{(i)}\T(\bbeta_r-\bbeta_{r'}) \); in particular a model with
     baseline \( s \) is the same model, with coefficients \( \bbeta_r^{(s)}=\bbeta_r-\bbeta_s \);

2. refitting the model with baseline \( s \) gives
     \( \hbeta_r^{(s)}=\hbeta_r-\hbeta_s \), the same maximized log-likelihood and the same
     fitted probabilities;

3. the estimated covariance matrix of \( \hbeta^{(s)} \) obtained from the refitted model
     equals \( \bL\,\hat{\Cov}(\hbeta)\,\bL\T \), where \( \bL \) is the linear map of (a); so
     standard errors of pairwise contrasts do not depend on which baseline was used.
:::

:::

::: {.proof}
(a) Subtract two instances of @eq-mlt-baseline; the baseline cancels because
\( \log(\pi_{ir}/\pi_{ir'})=\log(\pi_{ir}/\pi_{ic})-\log(\pi_{ir'}/\pi_{ic}) \).

(b) The map \( \bL:\bbeta\mapsto\bbeta^{(s)} \) defined by \( \bbeta_r^{(s)}=\bbeta_r-\bbeta_s \)
for \( r\neq s \) is linear, and it is invertible: from \( \bbeta^{(s)} \) one recovers
\( \bbeta_r=\bbeta^{(s)}_r-\bbeta^{(s)}_c \) with \( \bbeta^{(s)}_s=\bzero \). By (a) the two
parameterizations give exactly the same family of probability vectors
\( \boldsymbol{\uppi}_i \), so the two likelihoods are the same function composed with \( \bL \).
A reparameterization by a bijection carries maximizers to maximizers and leaves the
maximum value and every fitted quantity unchanged.

(c) Under a linear reparameterization \( \boldsymbol{\uppsi}=\bL\bbeta \) the information
transforms as \( \mathbf{J}^{(s)}=\bL^{-\top}\mathbf{J}\bL^{-1} \), by the chain rule applied to
the score. Hence \( (\mathbf{J}^{(s)})^{-1}=\bL\mathbf{J}^{-1}\bL\T \), which is exactly the delta-method covariance of
\( \bL\hbeta \) — and the delta method is exact here because \( \bL \) is linear.
:::

So the coefficients are log odds ratios: \( \hat{\beta}_{rj} \) is the estimated change in
\( \log(\pi_r/\pi_c) \) per unit of regressor \( j \), and differences of coefficients are log
odds ratios between two non-baseline categories. The baseline is a choice of *display*,
not of model. What is not invariant is the reading of a single coefficient, and here is
the trap.

::: {#prp-mlt-derivative}
[A coefficient is not a derivative of a probability]

In @def-mlt-baseline, with \( \bbeta_c=\bzero \) and \( x_j \) a continuous regressor,
\[
\frac{\partial \pi_{ir}}{\partial x_{ij}}
=\pi_{ir}\Bigl(\beta_{rj}-\sum_{s=1}^{c}\pi_{is}\beta_{sj}\Bigr).
\]{#eq-mlt-derivative}

So the \( c \) derivatives sum to zero, the derivative for category \( r \) is positive exactly
when \( \beta_{rj} \) exceeds the probability-weighted average
\( \bar{\beta}_j(\x_{(i)})=\sum_s\pi_{is}\beta_{sj} \), and the baseline, whose coefficient is
zero by construction, has a strictly positive derivative whenever
\( \bar{\beta}_j(\x_{(i)})<0 \).
:::

::: {.proof}
Write \( \pi_{ir}=e^{\eta_{ir}}/\sum_s e^{\eta_{is}} \) with
\( \partial\eta_{is}/\partial x_{ij}=\beta_{sj} \). The quotient rule gives
\[
\frac{\partial\pi_{ir}}{\partial x_{ij}}
=\pi_{ir}\beta_{rj}-\pi_{ir}\frac{\sum_s e^{\eta_{is}}\beta_{sj}}{\sum_s e^{\eta_{is}}}
=\pi_{ir}\Bigl(\beta_{rj}-\sum_s\pi_{is}\beta_{sj}\Bigr).
\]
Summing over \( r \) gives \( \bar{\beta}_j-\bar{\beta}_j=0 \). The remaining statements read off
the sign of the bracket, with \( \beta_{cj}=0 \).
:::

Since \( \bar{\beta}_j(\x) \) moves with \( \x \), a fitted probability can rise over one range
of the regressor and fall over another, whatever the sign of its coefficient. Read
coefficients as odds ratios; read effects on probabilities off the fitted curves.

## Travel mode choice

::: {#exm-mlt-mode}
[Choosing how to travel]

Two hundred and ten travellers on the Sydney–Melbourne corridor each chose one of four
modes: air (58), train (63), bus (30) or car (59). The regressors
are the household income in units of $10,000 and the size of the travelling party.
Take car as the baseline. Fitting @def-mlt-baseline by Newton–Raphson from
\( \bbeta=\bzero \) takes 7 steps to full double precision, and gives the
following coefficients, with standard errors from the inverse information in
parentheses.

| logit | intercept | income (per $10,000) | party size |
|---|---|---|---|
| air vs car | 0.9435 (0.5498) | 0.0354 (0.1030) | \( -0.6006 \) (0.1992) |
| train vs car | 2.4938 (0.5357) | \( -0.5731 \) (0.1184) | \( -0.3098 \) (0.1956) |
| bus vs car | 1.9780 (0.6717) | \( -0.3033 \) (0.1322) | \( -0.9404 \) (0.3245) |

Income does not distinguish air from car (\( z=0.34 \)) but moves people strongly away
from the ground alternatives: each extra $10,000 of household income multiplies the odds
of train rather than car by \( e^{-0.5731}=0.5638 \). Read through
@prp-mlt-invariance, the same fit raises the odds of air over train by
\( e^{0.6085} \) per $10,000. Party size pushes travellers towards the
car, most sharply away from the bus. The likelihood ratio statistic against the
intercept-only model is \( 2(-253.3408+283.7588)=60.836 \) on six
degrees of freedom.
:::

```{.python .run #cell-mode-choice-data}
import numpy as np
import statsmodels.api as sm

MODES = ["air", "train", "bus", "car"]          # car is category 4, the baseline

raw = sm.datasets.modechoice.load_pandas().data
chosen = raw[raw["choice"] == 1].sort_values("individual").reset_index(drop=True)
y = chosen["mode"].astype(int).to_numpy()       # 1 air, 2 train, 3 bus, 4 car
n, c = len(y), 4
X = np.column_stack([np.ones(n), chosen["hinc"] / 10.0, chosen["psize"]])
Y = np.zeros((n, c))
Y[np.arange(n), y - 1] = 1.0                    # indicator matrix, baseline last

print("chosen mode counts:", Y.sum(axis=0))
```

The fit itself is @thm-mlt-score written out: build the probabilities, form the score
and the block information, solve, repeat.

```{.python .run #cell-mode-choice-fit}
def probs(X, B):
    """Fitted probabilities of the baseline-category logit model; B is p x (c-1)."""
    eta = np.column_stack([X @ B, np.zeros(len(X))])     # baseline has eta = 0
    eta = eta - eta.max(axis=1, keepdims=True)           # stabilize the exponentials
    E = np.exp(eta)
    return E / E.sum(axis=1, keepdims=True)


def loglik(X, Y, B):
    return float(np.sum(Y * np.log(probs(X, B))))


def mnlogit(X, Y, tol=1e-10, maxit=50):
    """Newton-Raphson (= Fisher scoring, the link is canonical) for baseline logits."""
    n, p = X.shape
    c = Y.shape[1]
    B = np.zeros((p, c - 1))
    for it in range(1, maxit + 1):
        P = probs(X, B)
        U = (X.T @ (Y[:, :c - 1] - P[:, :c - 1])).T.ravel()       # score, stacked by category
        J = np.zeros(((c - 1) * p, (c - 1) * p))                  # expected = observed information
        for r in range(c - 1):
            for s in range(c - 1):
                w = P[:, r] * ((r == s) - P[:, s])
                J[r * p:(r + 1) * p, s * p:(s + 1) * p] = X.T @ (w[:, None] * X)
        step = np.linalg.solve(J, U)
        B = B + step.reshape(c - 1, p).T
        if np.max(np.abs(step)) < tol:
            break
    return B, J, it


B, J, iters = mnlogit(X, Y)
se = np.sqrt(np.diag(np.linalg.inv(J))).reshape(c - 1, X.shape[1]).T

for r in range(c - 1):
    print(f"{MODES[r]:>6s} vs car:  " + "  ".join(f"{b:7.4f} ({s:.4f})" for b, s in zip(B[:, r], se[:, r])))
print(f"log-likelihood {loglik(X, Y, B):.4f} after {iters} Newton steps")
```

The script behind this section checks the output against the `MNLogit` routine of
statsmodels, coefficient by coefficient. It also carries out the substitution of
@prp-mlt-invariance, refitting with train as the baseline and confirming that
subtracting the train column from the fitted coefficient matrix reproduces the refit and
leaves every fitted probability unchanged.

```{.python .run #cell-mode-choice-baseline}
# changing the baseline: subtract the train column (category 2) from every column
Btrain = B - B[:, [1]]
Btrain = np.delete(np.column_stack([Btrain, -B[:, 1]]), 1, axis=1)      # drop train, add car
Ytrain = Y[:, [0, 2, 3, 1]]                                            # put train last
Bcheck, _, _ = mnlogit(X, Ytrain)

print("air vs train, refitted :", np.round(Bcheck[:, 0], 6))
print("air vs train, subtracted:", np.round(Btrain[:, 0], 6))
print("fitted probabilities identical:", np.allclose(probs(X, Bcheck)[:, [0, 3, 1, 2]], probs(X, B)))
```

[Figure 36.1.1](01-baseline-category-logits.html#fig-mlt-mode) plots the fitted probabilities
against income, and the bus curve in the left panel is @prp-mlt-derivative made visible:
the coefficient of income in the bus logit is \( -0.3033 \), yet the fitted
probability of the bus *rises* over the first twenty thousand dollars before falling. At
an income of $10,000 with a party of one the fitted probabilities are
\( (0.1528,0.5243,0.2182,0.1047) \) for air, train, bus and car, with
derivatives per $10,000 of \( (0.0606,-0.1111,
0.0126,0.0378) \); the car, whose coefficient is fixed at zero, has the
second largest of the four, behind air alone.

::: {when-format="html"}
![**Figure 36.1.1.** Fitted mode probabilities from the baseline-category logit model,
against household income, for a traveller alone (left) and in a party of three (right).
Coefficients are log odds ratios against the car; the curves show what those
coefficients do to the probabilities themselves.](mode_probabilities.svg){#fig-mlt-mode width=100%}
:::

::: {when-format="pdf"}
![Fitted mode probabilities from the baseline-category logit model, against household
income, for a traveller alone (left) and in a party of three (right). Coefficients are
log odds ratios against the car; the curves show what those coefficients do to the
probabilities themselves.](mode_probabilities.pdf){width=100%}
:::

```{.python .run #cell-mode-choice-derivative}
# the derivative of a fitted probability is not the coefficient
Bfull = np.column_stack([B, np.zeros(X.shape[1])])          # append the baseline's column of zeros
x0 = np.array([1.0, 1.0, 1.0])                              # income $10,000, party of one
p0 = probs(x0[None, :], B)[0]
dp0 = p0 * (Bfull[1] - p0 @ Bfull[1])                       # d pi_r / d (income in $10,000s)

print("fitted probabilities", np.round(p0, 4))
print("derivatives         ", np.round(dp0, 4))
print("income coefficients ", np.round(Bfull[1], 4))
```

::: {.remark}
[Fitting one logit at a time]

It is tempting to estimate \( \bbeta_r \) by discarding every observation outside
\( \{r,c\} \) and running an ordinary binary logistic regression. By @eq-mlt-baseline,
\[
P(Y_i=r\mid Y_i\in\{r,c\},\x_{(i)})
=\frac{\pi_{ir}}{\pi_{ir}+\pi_{ic}}
=\frac{e^{\x_{(i)}\T\bbeta_r}}{1+e^{\x_{(i)}\T\bbeta_r}},
\]
so the conditional model is exactly the binary logistic model of @def-bin-logistic with
the same \( \bbeta_r \), and the individualized fit is consistent. It is not efficient, and
the \( c-1 \) fits are correlated in a way the separate outputs do not report. Begg and
Gray (1984) quantified the loss, small when the baseline is a large category and severe
otherwise.
:::

## Inference

The argument of @thm-glm-asymptotics carries over verbatim in block notation, although
the theorem there is stated for a scalar response: by @thm-mlt-score the score is again a
sum of independent mean-zero vectors, the information is again nonrandom, and the global
concavity its proof leans on is @thm-mlt-score(c). Under the analogues of its conditions
— categories and regressors fixed, sample size growing, the leverages of the stacked
design tending to zero, the true parameter interior — the
maximum likelihood estimate is consistent and asymptotically normal with covariance the
inverse information; nested models are compared by the analysis of deviance
of @thm-glm-deviance.

Two cautions are specific to categorical responses. The deviance is *not* a
goodness-of-fit statistic for ungrouped data: when every \( m_i=1 \) the saturated model
has as many parameters as observations and the \( \chi^2 \) approximation fails, as
@thm-glm-deviance says in general and @prp-bin-fit for binary data. And a regressor is
not one coefficient here but \( c-1 \), so "does income matter?" has \( c-1 \) degrees of freedom;
testing the coefficients separately and reporting the smallest \( p \)-value is the
multiplicity error of [Chapter 13](../ch13-multiplicity/index.html) in a new costume.

## Exercises

### A. Check your understanding

::: {#exr-mlt-count-parameters}
[A1]

A response has \( c=5 \) categories and the model matrix has \( p=4 \) columns including the
intercept. How many free coefficients has the baseline-category logit model? How many
has the saturated model for \( n \) ungrouped observations? How many degrees of freedom
does the test that one regressor may be dropped have?
:::

::: {.solution}
The model has \( (c-1)p=16 \) free coefficients. The saturated model for ungrouped data
gives each observation its own probability vector, \( n(c-1) \) free parameters. Dropping
one regressor removes one coefficient from each of the \( c-1 \) logits, so the test has
\( c-1=4 \) degrees of freedom.
:::

### B. Practice

::: {#exr-mlt-pairwise}
[B1]

From the fitted coefficients of @exm-mlt-mode, compute the estimated log odds ratio of
bus against train per $10,000 of income, and the estimated odds ratio of air against
bus per extra person in the party. State which of these depend on the choice of
baseline.
:::

::: {.solution}
By @prp-mlt-invariance the bus-against-train income coefficient is
\( -0.3033-(-0.5731)=0.2698 \), so the odds of bus rather than train rise by
a factor \( e^{0.2698} \) per $10,000. For party size, air against bus is
\( -0.6006-(-0.9404)=0.3398 \), an odds ratio of \( e^{0.3398} \). Neither
depends on the baseline: by @prp-mlt-invariance(a) a difference of coefficients is the
pairwise log odds ratio \( \log(\pi_{ir}/\pi_{ir'}) \) per unit of the regressor, which the
model specifies without reference to any baseline.
:::

::: {#exr-mlt-information-singular}
[B2]

Show directly that the \( c\times c \) matrix
\( \bSigma(\boldsymbol{\uppi})=\diag(\boldsymbol{\uppi})-\boldsymbol{\uppi}\boldsymbol{\uppi}\T \) is
nonnegative definite with \( \bSigma\bone=\bzero \), and that its restriction to the first
\( c-1 \) coordinates is positive definite when every \( \pi_r>0 \).
:::

::: {.solution}
For any \( \mathbf{t} \), \( \mathbf{t}\T\bSigma\mathbf{t}=\sum_r\pi_rt_r^2-(\sum_r\pi_rt_r)^2 \) is the
variance of a random variable taking the value \( t_r \) with probability \( \pi_r \), hence
nonnegative, and it vanishes iff all \( t_r \) with \( \pi_r>0 \) are equal. Taking
\( \mathbf{t}=\bone \) gives \( \bSigma\bone=\bzero \). If \( \mathbf{t} \) has \( t_c=0 \) and
\( \mathbf{t}\T\bSigma\mathbf{t}=0 \), then all \( t_r \) are equal to \( t_c=0 \); so the leading
\( (c-1)\times(c-1) \) block is positive definite.
:::

### C. Going deeper

::: {#exr-mlt-adjacent}
[C1]

The **adjacent-category logit model** for an ordered response specifies
\( \log(\pi_{i,r+1}/\pi_{ir})=\theta_r+\x_{(i)}\T\bbeta \) for \( r=1,\dots,c-1 \), with
\( \x_{(i)} \) not containing an intercept. Show that it is the baseline-category model
@def-mlt-baseline with slope vectors \( \bbeta_r=(c-r)\,\bbeta^{*} \) for a single vector
\( \bbeta^{*} \), plus a free intercept for each \( r \): a submodel obtained by \( (c-2)p \)
linear restrictions, fitted with the same software after recoding the model matrix.
:::

::: {.solution}
Summing the adjacent logits from \( r \) to \( c-1 \) gives
\( \log(\pi_{ic}/\pi_{ir})=\sum_{s=r}^{c-1}\theta_s+(c-r)\,\x_{(i)}\T\bbeta \), so
\( \log(\pi_{ir}/\pi_{ic})=-\sum_{s\ge r}\theta_s-(c-r)\x_{(i)}\T\bbeta \). This is
@eq-mlt-baseline with a free intercept for each \( r \) and slope vector proportional to
\( \bbeta^{*}=-\bbeta \) with the known factor \( c-r \). Fitting is a baseline-category fit on
the model matrix whose slope block for logit \( r \) is \( (c-r)\X \).
:::

::: {#exr-mlt-iia}
[C2]

Show that in @def-mlt-baseline the ratio \( \pi_{ir}/\pi_{ir'} \) does not involve any
third category — the *independence of irrelevant alternatives*. Construct a
three-category example in which this is substantively absurd, and say what the
assumption forces about the fitted probabilities when a new category is added that is
nearly a copy of an existing one.
:::

