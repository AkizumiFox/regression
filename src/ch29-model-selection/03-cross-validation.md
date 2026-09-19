# Cross-validation

The criteria of [Section 29.2](02-cp-aic-bic.html) correct the training error by a formula that rests on assumptions.
Cross-validation instead predicts cases that the fit has not seen. It needs neither \( \sigma^2 \) nor degrees of freedom,
and applies to any procedure, including one that selects a model. Its price is computation, a bias that is easy to
describe, and a variance that is not.

## K-fold and leave-one-out

::: {#def-sel-cv}
[Cross-validation]

Let the cases \( (\x_{(i)},y_i) \), \( i=1,\dots,n \), be divided into \( K \) disjoint **folds** \( I_1,\dots,I_K \). For each
\( k \), let \( \hat{f}^{(-k)} \) be the prediction rule obtained by applying the procedure to the cases outside \( I_k \). The
**\( K \)-fold cross-validation** estimate of prediction error is
\[
\text{CV}_K=\frac1n\sum_{k=1}^K\sum_{i\in I_k}\bigl(y_i-\hat{f}^{(-k)}(\x_{(i)})\bigr)^2 .
\]
With \( K=n \), each fold is a single case, and \( \text{CV}_n \) is **leave-one-out** cross-validation. With \( K=2 \) and
one fold used only for fitting, the estimate is a **holdout** or **validation-set** estimate.
:::

The idea is old (Stone 1974; Geisser 1975; Allen 1974 for least squares). What cross-validation estimates is
clearest when the cases are a random sample.

::: {#prp-sel-cv-target}
[What K-fold cross-validation estimates]

Let the cases be independent draws of \( (\mathbf{X},Y) \), let the folds be fixed in advance with equal sizes \( n/K \),
and let \( \text{Err}(m)=\E\bigl(Y_0-\hat{f}_m(\mathbf{X}_0)\bigr)^2 \) be the expected squared error of the procedure fitted to
\( m \) cases, for an independent new case \( (\mathbf{X}_0,Y_0) \). Then
\[
\E(\text{CV}_K)=\text{Err}(n-n/K).
\]
:::

::: {.proof}
For \( i\in I_k \), the rule \( \hat{f}^{(-k)} \) is the procedure applied to \( n-n/K \) independent cases, and case \( i \) is
independent of them and has the same distribution as \( (\mathbf{X}_0,Y_0) \). So each term of the double sum has
expectation \( \text{Err}(n-n/K) \).
:::

The expectation is over training sets as well as new cases. Cross-validation therefore estimates the expected
error of the *procedure* at a slightly smaller sample size, not the error of the particular fit computed from all
\( n \) cases. For least squares with \( k \) normal regressors and a correct model, @thm-cor-prediction-error gives
\( \text{Err}(m)=\sigma^2(1+1/m)(m-2)/(m-k-2) \), so the bias of \( \text{CV}_K \) can be computed exactly.

::: {#exm-sel-cv-bias}
[Bias and spread of K-fold cross-validation]

Take \( n=40 \) cases with \( k=10 \) independent standard normal regressors, an intercept, slopes \( 0.5 \) and
\( \sigma=1 \), and fit by least squares. The target is \( \text{Err}(40)=1.391 \). Over \( 4000 \) training sets, the
averages of \( \text{CV}_K \) for \( K=2,5,10,40 \) are \( 2.378 \), \( 1.551 \), \( 1.461 \) and \( 1.408 \), against
\( \text{Err}(n-n/K) \) of \( 2.363 \), \( 1.547 \), \( 1.456 \) and \( 1.406 \) from @prp-sel-cv-target. Fitting to half the data
inflates the error by \( 70\% \) here, because ten slopes are estimated from twenty cases. Their standard
deviations across training sets are \( 1.062 \), \( 0.477 \), \( 0.411 \) and \( 0.378 \). For this stable procedure
leave-one-out has both the smallest bias and the smallest variance.

The error of the fitted model itself, \( \sigma^2+\norm{\hat{\bbeta}-\bbeta}^2 \) for a new standard normal case, varies too:
its mean is \( 1.395 \) and its standard deviation \( 0.203 \). Panel (b) of [Figure 29.3.1](#fig-sel-cv-bias) plots
\( \text{CV}_{10} \) against it. The correlation is \( -0.005 \) (\( -0.014 \) for leave-one-out). Cross-validation reports how well
least squares does on average with forty such cases, not how well this fit will do.
:::

::: {when-format="html"}
![**Figure 29.3.1.** (a) Mean and standard deviation of \( \text{CV}_K \) over \( 4000 \) training sets
(\( n=40 \), ten regressors), with the expected error at sample size \( n-n/K \) (crosses) and at \( n \) (dashed).
(b) \( \text{CV}_{10} \) against the prediction error of the model actually fitted.](cv_bias_variance.svg){#fig-sel-cv-bias width=100%}
:::

::: {when-format="pdf"}
![(a) Mean and standard deviation of \( \text{CV}_K \) over \( 4000 \) training sets
(\( n=40 \), ten regressors), with the expected error at sample size \( n-n/K \) (crosses) and at \( n \) (dashed).
(b) \( \text{CV}_{10} \) against the prediction error of the model actually fitted.](cv_bias_variance.pdf){width=100%}
:::

```{.python .run #cell-cross-validation-kfold}
import numpy as np
rng = np.random.default_rng(2905)
n_s, k_s, sigma = 40, 10, 1.0
beta = np.r_[1.0, np.full(k_s, 0.5)]                         # intercept and ten slopes


def err_theory(m):
    """Expected error of a least squares fit to m cases (Chapter 14, normal regressors)."""
    return sigma**2 * (1 + 1 / m) * (m - 2) / (m - k_s - 2)


def cv_kfold(X, y, K):
    folds = np.arange(len(y)) % K                            # equal folds (cases are exchangeable)
    sq = np.empty(len(y))
    for f in range(K):
        tr, te = folds != f, folds == f
        b, *_ = np.linalg.lstsq(X[tr], y[tr], rcond=None)
        sq[te] = (y[te] - X[te] @ b) ** 2
    return sq.mean()


Ks = [2, 5, 10, n_s]


def simulate(reps):
    cv = np.empty((reps, len(Ks)))
    err_cond = np.empty(reps)                                 # error of the fit at hand
    for r in range(reps):
        Xr = np.column_stack([np.ones(n_s), rng.normal(size=(n_s, k_s))])
        yr = Xr @ beta + sigma * rng.normal(size=n_s)
        b, *_ = np.linalg.lstsq(Xr, yr, rcond=None)
        err_cond[r] = sigma**2 + np.sum((b - beta) ** 2)     # new case X0 ~ N(0, I)
        cv[r] = [cv_kfold(Xr, yr, K) for K in Ks]
    return cv, err_cond


cv, err_cond = simulate(300)                                 # the book's figures use 4000
for j, K in enumerate(Ks):
    print(f"K = {K:2d}: mean {cv[:, j].mean():.3f} (theory {err_theory(n_s - n_s // K):.3f}), sd {cv[:, j].std():.3f}")
print(f"error of the fit at hand: mean {err_cond.mean():.3f} (theory {err_theory(n_s):.3f})")
print("correlation of CV_10 with it:", np.corrcoef(cv[:, 2], err_cond)[0, 1])
```

Bates, Hastie and Tibshirani (2024) study this phenomenon for linear models. The spread falling with \( K \) is typical of
least squares, whose fits change little when a few cases are removed; for procedures that jump, such as subset
selection, leave-one-out can be the more variable choice (Breiman and Spector 1992). Five or ten folds are the usual compromise.

## The leverage shortcut

Leave-one-out looks expensive: \( n \) refits. For least squares and for ridge regression it costs one fit. The
reason is a small fixed-point argument.

::: {#thm-sel-loocv}
[Leave-one-out without refitting]

Let \( \boldsymbol{\Omega} \) be a \( p\times p \) nonnegative definite matrix, let \( \hat{\bbeta} \) minimize
\( \norm{\y-\X\bb}^2+\bb\T\boldsymbol{\Omega}\bb \), and let \( \hat{\bbeta}_{(i)} \) minimize the same criterion with case \( i \) removed.
Assume that \( \mathbf{A}_{(i)}=\X_{(i)}\T\X_{(i)}+\boldsymbol{\Omega} \) is nonsingular for every \( i \), where \( \X_{(i)} \) is \( \X \) without
row \( i \). Let \( \bS=\X(\X\T\X+\boldsymbol{\Omega})^{-1}\X\T \) and \( \hat{\y}=\bS\y \). Then:

::: {.enumerate options="label=(\alph*)"}
1. \( s_{ii}<1 \) for every \( i \);

2. \( y_i-\x_{(i)}\T\hat{\bbeta}_{(i)}=(y_i-\hat{y}_i)/(1-s_{ii}) \);

3. consequently
   \[
   \text{CV}_n=\frac1n\sum_{i=1}^n\Bigl(\frac{y_i-\hat{y}_i}{1-s_{ii}}\Bigr)^2 .
   \]{#eq-sel-loo-shortcut}
:::

For least squares (\( \boldsymbol{\Omega}=\mathbf{0} \)) this is @thm-res-deletion(d), and \( \text{CV}_n=\text{PRESS}/n \) with PRESS
as in @eq-res-press.
:::

::: {.proof}
(a) Write \( \x=\x_{(i)} \) and \( \mathbf{A}=\mathbf{A}_{(i)} \), which is positive definite because it is nonnegative definite and
nonsingular. Then \( \X\T\X+\boldsymbol{\Omega}=\mathbf{A}+\x\x\T \), and \( (\mathbf{A}+\x\x\T)\mathbf{A}^{-1}\x=\x(1+a) \) with
\( a=\x\T\mathbf{A}^{-1}\x\ge0 \). Hence \( (\mathbf{A}+\x\x\T)^{-1}\x=\mathbf{A}^{-1}\x/(1+a) \) and
\( s_{ii}=\x\T(\mathbf{A}+\x\x\T)^{-1}\x=a/(1+a)<1 \).

(b) Let \( \tilde{y}_i=\x_{(i)}\T\hat{\bbeta}_{(i)} \), the leave-one-out prediction, and let \( \tilde{\y} \) be \( \y \) with \( y_i \) replaced by
\( \tilde{y}_i \). For any \( \bb \), the full criterion evaluated at \( \tilde{\y} \) is the leave-one-out criterion plus
\( (\tilde{y}_i-\x_{(i)}\T\bb)^2\ge0 \). At \( \bb=\hat{\bbeta}_{(i)} \) the extra term is zero and the leave-one-out criterion is at its
minimum. So \( \hat{\bbeta}_{(i)} \) minimizes the full criterion for the data \( \tilde{\y} \), and since that minimizer is unique
(\( \X\T\X+\boldsymbol{\Omega}=\mathbf{A}+\x\x\T \) is positive definite), \( \X\hat{\bbeta}_{(i)}=\bS\tilde{\y} \). Its \( i \)th entry reads
\[
\tilde{y}_i=(\bS\tilde{\y})_i=(\bS\y)_i-s_{ii}(y_i-\tilde{y}_i)=\hat{y}_i-s_{ii}(y_i-\tilde{y}_i).
\]
Rearranging, \( (1-s_{ii})(y_i-\tilde{y}_i)=y_i-\hat{y}_i \), and (a) allows division by \( 1-s_{ii} \).

(c) follows from (b), and for \( \boldsymbol{\Omega}=\mathbf{0} \), \( s_{ii}=h_{ii} \) and (b) is @thm-res-deletion(d).
:::

The argument uses only that deleting a case gives the same answer as replacing its response by its own prediction.
It covers the ridge estimators of [Chapter 27](../ch27-shrinkage/index.html) and the smoothing splines of Chapter 43,
but not the lasso or subset selection, where deleting a case can change which coefficients are nonzero.

[Section 27.2](../ch27-shrinkage/02-ridge.html) already used the leave-one-out shortcut and GCV to choose the ridge
\( \lambda \); the theorem proves the shortcut for every quadratic penalty. Replacing each \( s_{ii} \) by its average
\( \tr(\bS)/n \) gives a criterion that needs only the trace.

::: {#def-sel-gcv}
[Generalized cross-validation]

For a linear smoother \( \hat{\y}=\bS\y \) with \( \tr(\bS)<n \),
\[
\text{GCV}=\frac{n^{-1}\norm{\y-\hat{\y}}^2}{\bigl(1-\tr(\bS)/n\bigr)^2} .
\]
:::

Craven and Wahba (1979) introduced GCV for smoothing splines and Golub, Heath and Wahba (1979) for ridge regression,
where it is leave-one-out after a rotation that equalizes the diagonal of \( \bS \). Since \( (1-x)^{-2}=1+2x+O(x^2) \),
\( \text{GCV}\approx\overline{\text{err}}\,(1+2\,\text{df}/n) \), the covariance penalty @eq-sel-covariance-penalty with
\( \sigma^2 \) replaced by \( \overline{\text{err}} \). For least squares \( \text{GCV}=ns^2/(n-p) \), and for a correct model (\( \bmu\in\C(\X) \)), by @prp-res-press,
\( \E(\text{PRESS})/n\ge\sigma^2n/(n-p)=\E(\text{GCV}) \): high-leverage cases inflate leave-one-out.

For the state regression of @exm-sel-state-cp with all five regressors, \( \text{PRESS}/n=0.1426 \) and
\( \text{GCV}=0.1176 \). For ridge regression on the standardized regressors with \( \lambda=5 \) (intercept unpenalized),
\( \tr(\bS)=4.870 \), and the shortcut gives \( \text{CV}_n=0.1309 \), equal to the brute-force value from fifty refits to
rounding error, against \( \text{GCV}=0.1140 \).

```{.python .run #cell-cross-validation-shortcut}
import numpy as np
import statsmodels.api as sm
data = sm.datasets.statecrime.load_pandas().data.drop(index="District of Columbia")
names = ["hs_grad", "poverty", "single", "white", "urban"]
y = np.log(data["violent"].to_numpy())
Zs = data[names].to_numpy()
Zs = (Zs - Zs.mean(axis=0)) / Zs.std(axis=0)               # standardized regressors
X = np.column_stack([np.ones(len(y)), Zs])
n, p = X.shape
Omega = np.diag([0.0] + [1.0] * (p - 1))                    # ridge penalty, intercept left free


def ridge_fit(X, y, lam):
    return np.linalg.solve(X.T @ X + lam * Omega, X.T @ y)


lam = 5.0
S = X @ np.linalg.solve(X.T @ X + lam * Omega, X.T)          # the smoother matrix
resid = y - S @ y
loo_short = resid / (1 - np.diag(S))                         # the shortcut
loo_brute = np.array([y[i] - X[i] @ ridge_fit(np.delete(X, i, 0), np.delete(y, i), lam)
                      for i in range(n)])                    # n refits
print("largest discrepancy:", np.max(np.abs(loo_short - loo_brute)))
cv_loo = np.mean(loo_short**2)
gcv = np.mean(resid**2) / (1 - np.trace(S) / n) ** 2
print(f"CV_n = {cv_loo:.5f}, GCV = {gcv:.5f}, tr(S) = {np.trace(S):.3f}")
```

::: {.remark}
[Sketch: leave-one-out behaves like AIC]

If every leverage is close to \( p/n \), then \( (1-h_{ii})^{-2}\approx1+2p/n \) and
\( \text{CV}_n\approx(\text{SSE}/n)(1+2p/n) \), so \( n\log\text{CV}_n\approx n\log(\text{SSE}/n)+2p \), the AIC of
@eq-sel-aic-normal up to a constant. Stone (1977) made this equivalence precise. Leave-one-out therefore shares AIC's
large-sample behaviour: efficient when no candidate is correct (Li 1987), not consistent when one is. Shao (1993) showed
that consistency requires leaving out a fraction of the cases that tends to one, the opposite of what prediction
suggests; \( K \)-fold with fixed \( K \) is not consistent either.
:::

## Choosing a model by cross-validation

To choose among candidates, take the smallest \( \text{CV}_K \). The fold errors \( \text{CV}^{(k)} \) give a rough standard
error \( \text{sd}(\text{CV}^{(1)},\dots,\text{CV}^{(K)})/\sqrt K \), which treats the folds as independent although their
training sets overlap; no unbiased estimate of the variance of \( \text{CV}_K \) exists (Bengio and Grandvalet 2004). The
**one-standard-error rule** (Breiman, Friedman, Olshen and Stone 1984) takes the simplest candidate within one
standard error of the minimum, trading a little estimated error for a smaller model.

::: {#exm-sel-state-cv}
[Cross-validating the state regressions]

For each size from \( 0 \) to \( 5 \), take the best subset of the state regressors by residual sum of squares and
estimate its prediction error by ten-fold cross-validation, with one random split into folds of five states.
[Figure 29.3.2](#fig-sel-cv-state) shows the result. The minimum, \( 0.1136 \) with standard error \( 0.0198 \), is at
two regressors (single parenthood and urbanization), the choice of \( C_p \) and AIC. Single parenthood alone has
\( 0.1162 \), well within one standard error, so the one-standard-error rule chooses it, the choice of BIC. The full
model has \( 0.1351 \) and the intercept alone \( 0.1993 \). Leave-one-out over all \( 32 \) candidates gives
\( 0.1129 \) for single parenthood and urbanization and \( 0.1130 \) for single parenthood alone, a tie, and \( 0.1184 \)
for poverty with single parenthood and urbanization.
:::

::: {when-format="html"}
![**Figure 29.3.2.** Ten-fold cross-validation error of the best subset of each size for log violent crime in the
\( 50 \) states, with one-standard-error bars. Dashed: the minimum plus one standard
error.](cv_state.svg){#fig-sel-cv-state width=70%}
:::

::: {when-format="pdf"}
![Ten-fold cross-validation error of the best subset of each size for log violent crime in the
\( 50 \) states, with one-standard-error bars. Dashed: the minimum plus one standard
error.](cv_state.pdf){width=70%}
:::

```{.python .run #cell-state-cv-cv}
import itertools
import numpy as np
import statsmodels.api as sm
data = sm.datasets.statecrime.load_pandas().data.drop(index="District of Columbia")
names = ["hs_grad", "poverty", "single", "white", "urban"]
y = np.log(data["violent"].to_numpy())
n = len(y)


def design(cols):
    return np.column_stack([np.ones(n)] + [data[c].to_numpy() for c in cols])


def sse(cols):
    X = design(cols)
    b, *_ = np.linalg.lstsq(X, y, rcond=None)
    return np.sum((y - X @ b) ** 2)


# best subset of each size, by residual sum of squares
best = [min(itertools.combinations(names, k), key=sse) for k in range(len(names) + 1)]

rng = np.random.default_rng(2906)
K = 10
folds = rng.permutation(np.arange(n) % K)                   # a random split into 10 folds of 5
cv, se = [], []
for cols in best:
    X = design(cols)
    fold_err = []
    for f in range(K):
        tr, te = folds != f, folds == f
        b, *_ = np.linalg.lstsq(X[tr], y[tr], rcond=None)
        fold_err.append(np.mean((y[te] - X[te] @ b) ** 2))
    cv.append(np.mean(fold_err))
    se.append(np.std(fold_err, ddof=1) / np.sqrt(K))        # the usual (optimistic) standard error
cv, se = np.array(cv), np.array(se)
k_min = int(np.argmin(cv))
k_1se = int(np.min(np.nonzero(cv <= cv[k_min] + se[k_min])[0]))
for k, cols in enumerate(best):
    print(f"{k} regressors {cols}: CV {cv[k]:.4f} (se {se[k]:.4f})")
print("minimum CV at size", k_min, "; one-standard-error rule chooses size", k_1se)
```

The example has a flaw: the best subset of each size was found from all fifty states, and only its fitting was
cross-validated.

## Cross-validating the whole procedure

Cross-validation estimates the error of whatever is repeated inside the folds. Every decision that looked at the
responses, such as screening regressors, choosing a subset, choosing \( \lambda \) or deleting outliers, must be repeated
inside each fold. Otherwise the held-out cases have already influenced the fit that predicts them.

::: {#exm-sel-nested}
[Screening, then cross-validating]

Take \( n=50 \) cases, \( 200 \) candidate regressors and a response that is independent of all of them, with
\( \sigma=1 \). The procedure keeps the \( 5 \) regressors most correlated with the response and fits least squares on
them. Every procedure has expected prediction error at least \( \sigma^2=1 \) here. Over \( 200 \) data sets:

- ten-fold cross-validation of the least squares fit, after screening on all the data, averages \( 0.720 \), and
  is below \( 1 \) in a fraction \( 0.93 \) of the data sets;
- ten-fold cross-validation that repeats the screening inside every fold averages \( 1.427 \);
- the actual prediction error of the screened and fitted model averages \( 1.372 \).

The wrong version reports that noise predicts the response. The right one is slightly pessimistic, because it
uses \( 45 \) cases instead of \( 50 \) (@prp-sel-cv-target).
:::

```{.python .run #cell-nested-cv-nested}
import numpy as np
rng = np.random.default_rng(2907)
n, m, keep, K = 50, 200, 5, 10


def screen(X, y):
    """Indices of the `keep` columns most correlated with y."""
    Xc, yc = X - X.mean(axis=0), y - y.mean()
    corr = np.abs(Xc.T @ yc) / np.sqrt(np.sum(Xc**2, axis=0) * (yc @ yc))
    return np.argsort(-corr)[:keep]


def fit(X, y):
    A = np.column_stack([np.ones(len(y)), X])
    return np.linalg.lstsq(A, y, rcond=None)[0]


def predict(b, X):
    return b[0] + X @ b[1:]


def one_data_set():
    X = rng.normal(size=(n, m))
    y = rng.normal(size=n)                                   # unrelated to every column
    folds = rng.permutation(np.arange(n) % K)
    cols = screen(X, y)                                      # screening on ALL the data
    wrong = right = 0.0
    for f in range(K):
        tr, te = folds != f, folds == f
        b = fit(X[tr][:, cols], y[tr])                       # wrong: screening done outside
        wrong += np.sum((y[te] - predict(b, X[te][:, cols])) ** 2)
        cols_f = screen(X[tr], y[tr])                        # right: screen again in the fold
        b = fit(X[tr][:, cols_f], y[tr])
        right += np.sum((y[te] - predict(b, X[te][:, cols_f])) ** 2)
    b = fit(X[:, cols], y)
    truth = 1 + b[0] ** 2 + np.sum(b[1:] ** 2)               # error on a new case
    return wrong / n, right / n, truth


res = np.array([one_data_set() for _ in range(200)])
print("average CV, screening outside the folds:", res[:, 0].mean().round(3))
print("average CV, screening inside the folds: ", res[:, 1].mean().round(3))
print("average true error of the fitted model: ", res[:, 2].mean().round(3))
```

When cross-validation is used both to choose among candidates and to report the error of the choice, the same problem
appears in a milder form. The smallest of many cross-validation estimates is biased downwards as an estimate of the error
of the model it picks ([Section 29.5](05-selection-bias.html)). The remedy is **nested cross-validation**: an outer loop
holds out folds for assessment, and inside each outer training set an inner cross-validation makes the choice (Varma
and Simon 2006). The outer estimate then describes the whole procedure, choice included, at sample size \( n-n/K \).

## Exercises

### A. Check your understanding

::: {#exr-sel-holdout}
[A1]

In @exm-sel-cv-bias, compute \( \text{Err}(20) \) from @thm-cor-prediction-error and confirm the value quoted for
\( K=2 \). Why is a holdout estimate that fits on half the data especially biased when the number of regressors is
large?
:::

::: {.solution}
\( \text{Err}(20)=(1+1/20)(18/8)=1.05\times2.25=2.3625 \). The factor \( (m-2)/(m-k-2) \) blows up as the training size
\( m \) approaches \( k+2 \). With ten regressors, halving \( m \) from \( 40 \) to \( 20 \) moves it from \( 38/28 \) to \( 18/8 \).
:::

### B. Practice

::: {#exr-sel-loo-mean}
[B1]

For the intercept-only model, show that the leave-one-out residual is \( n(y_i-\bar{y})/(n-1) \) and that
\( \text{CV}_n=ns^2/(n-1) \), where \( s^2 \) is the sample variance. Check this against @thm-sel-loocv.
:::

::: {.solution}
Without case \( i \) the mean is \( \bar{y}_{(i)}=(n\bar{y}-y_i)/(n-1) \), so
\( y_i-\bar{y}_{(i)}=(ny_i-n\bar{y})/(n-1) \). Squaring and averaging,
\( \text{CV}_n=n^{-1}\{n/(n-1)\}^2\sum_i(y_i-\bar{y})^2=n\sum_i(y_i-\bar{y})^2/(n-1)^2=ns^2/(n-1) \). In @thm-sel-loocv,
\( s_{ii}=h_{ii}=1/n \) and \( 1/(1-1/n)=n/(n-1) \).
:::

::: {#exr-sel-ridge-leverage}
[B2]

In @thm-sel-loocv, suppose \( \X \) has full column rank. Show that \( s_{ii}\le h_{ii} \), where \( h_{ii} \) is the least squares
leverage, so that the shortcut inflates ridge residuals less than least squares residuals.
:::

::: {.solution}
\( \X\T\X+\boldsymbol{\Omega}-\X\T\X=\boldsymbol{\Omega} \) is nonnegative definite, so \( (\X\T\X+\boldsymbol{\Omega})^{-1}\le(\X\T\X)^{-1} \) in the
nonnegative definite ordering (for positive definite \( \A\ge\B \), \( \B^{-1}-\A^{-1}=\B^{-1/2}(\I-\mathbf{C}^{-1})\B^{-1/2} \) with
\( \mathbf{C}=\B^{-1/2}\A\B^{-1/2}\ge\I \)). Evaluating both quadratic forms at \( \x_{(i)} \) gives \( s_{ii}\le h_{ii} \).
:::

### C. Going deeper

::: {#exr-sel-unequal-folds}
[C1]

Generalize @prp-sel-cv-target to folds of unequal sizes \( n_1,\dots,n_K \), and to leave-one-out for a design that is
fixed rather than random. In the fixed-design case, what does \( \E(\text{CV}_n) \) equal for least squares, and how does it
compare with the in-sample prediction error of @thm-sel-optimism?
:::

::: {#exr-sel-stone}
[C2]

Suppose \( \max_i\lvert h_{ii}-p/n\rvert\le\epsilon\,p/n \). Show that
\[
\frac{\text{SSE}}n\Bigl(1-\frac pn(1+\epsilon)\Bigr)^{-2}\ \ge\ \text{CV}_n\ \ge\ \frac{\text{SSE}}n\Bigl(1-\frac pn(1-\epsilon)\Bigr)^{-2},
\]
and deduce that \( n\log\text{CV}_n=n\log(\text{SSE}/n)+2p+O(p\epsilon+p^2/n) \). Which designs satisfy the hypothesis with
small \( \epsilon \)?
:::

