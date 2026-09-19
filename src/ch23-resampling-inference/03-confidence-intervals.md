# Bootstrap confidence intervals

A bootstrap distribution gives a confidence interval through its quantiles, in several ways. They agree to first order, which
is all the consistency theorems guarantee, and differ in how fast their coverage approaches the nominal level.

Let \( \text{se} \) be a standard error for \( \hat{\theta}=\mathbf{c}\T\hbeta \), classical or sandwich, and let \( \hat{\theta}^* \) and
\( \text{se}^* \) be computed from a resample in the same way. Write \( \hat{q}(\gamma) \) and \( \hat{t}(\gamma) \) for the \( \gamma \)
quantiles of the bootstrap distributions of \( \hat{\theta}^* \) and of the **bootstrap \( t \) statistic**
\( T^*=(\hat{\theta}^*-\hat{\theta})/\text{se}^* \).

::: {#def-bs-intervals}
[Percentile, basic and studentized intervals]

The nominal \( 1-\alpha \) bootstrap intervals for \( \theta \) are:

::: {.enumerate options="label=(\alph*)"}
1. the **percentile** interval \( \bigl[\hat{q}(\alpha/2),\ \hat{q}(1-\alpha/2)\bigr] \);

2. the **basic** interval \( \bigl[2\hat{\theta}-\hat{q}(1-\alpha/2),\ 2\hat{\theta}-\hat{q}(\alpha/2)\bigr] \);

3. the **studentized** (or bootstrap-\( t \)) interval
   \( \bigl[\hat{\theta}-\hat{t}(1-\alpha/2)\,\text{se},\ \hat{\theta}-\hat{t}(\alpha/2)\,\text{se}\bigr] \).
:::

:::

The basic interval translates the bootstrap idea directly: if \( \hat{\theta}-\theta \) had the distribution of
\( \hat{\theta}^*-\hat{\theta} \), then \( \hat{q}(\alpha/2)-\hat{\theta}\le\hat{\theta}-\theta\le\hat{q}(1-\alpha/2)-\hat{\theta} \) would have
probability \( 1-\alpha \), and solving for \( \theta \) gives (b). The studentized interval does the same with the root
\( (\hat{\theta}-\theta)/\text{se} \), as the \( t \) interval of @cor-opt-t does. The percentile interval looks naive, but if some
increasing \( \phi \) makes \( \phi(\hat{\theta})-\phi(\theta) \) a pivot symmetric about zero, it is exact without knowing \( \phi \) (@exr-bs-percentile-transform). It is also **transformation respecting**, because quantiles commute with increasing functions;
the basic and studentized intervals are not.

## First-order correctness

The consistency theorems translate into coverage through a lemma about quantiles.

::: {#lem-bs-quantiles}
[Consistent bootstrap quantiles]

Let \( R_n \) be random variables with distribution functions \( H_n \), and let \( \hat{H}_n \) be random distribution functions with
\( \sup_x\lvert\hat{H}_n(x)-H_n(x)\rvert\to0 \) in probability. Suppose \( H_n(x)\to H(x) \) for every \( x \), where \( H \) is continuous
and strictly increasing on \( \Real \). For \( 0<\gamma<1 \), let \( \hat{\xi}_n(\gamma)=\inf\{x:\hat{H}_n(x)\ge\gamma\} \) and
\( \xi(\gamma)=H^{-1}(\gamma) \). Then \( \hat{\xi}_n(\gamma)\to\xi(\gamma) \) in probability, and
\( \Pr\bigl(R_n\le\hat{\xi}_n(\gamma)\bigr)\to\gamma \).
:::

::: {.proof}
By Pólya's theorem \( \sup_x\lvert H_n(x)-H(x)\rvert\to0 \), so \( \delta_n=\sup_x\lvert\hat{H}_n(x)-H(x)\rvert\to0 \) in probability. Fix
\( \epsilon>0 \), write \( \xi=\xi(\gamma) \), and let
\( \kappa=\min\{\gamma-H(\xi-\epsilon),\ H(\xi+\epsilon)-\gamma\}>0 \). On the event \( \delta_n<\kappa \),
\( \hat{H}_n(\xi-\epsilon)<\gamma \) and \( \hat{H}_n(\xi+\epsilon)>\gamma \), so \( \xi-\epsilon\le\hat{\xi}_n(\gamma)\le\xi+\epsilon \). This proves
the first claim. For the second,
\[
\begin{aligned}
\Pr\bigl(R_n\le\hat{\xi}_n\bigr)&\ge H_n(\xi-\epsilon)-\Pr\bigl(\lvert\hat{\xi}_n-\xi\rvert>\epsilon\bigr),\\
\Pr\bigl(R_n\le\hat{\xi}_n\bigr)&\le H_n(\xi+\epsilon)+\Pr\bigl(\lvert\hat{\xi}_n-\xi\rvert>\epsilon\bigr).
\end{aligned}
\]
Let \( n\to\infty \) and then \( \epsilon\to0 \), using the continuity of \( H \) at \( \xi \).
:::

::: {#cor-bs-interval-consistency}
[Bootstrap intervals are first-order correct]

Under the conditions of @thm-bs-residual, the percentile and basic intervals built from the residual bootstrap have
coverage probability tending to \( 1-\alpha \).
:::

::: {.proof}
Take \( R_n=(\hat{\theta}-\theta)/\norm{\mathbf{a}_n} \), with \( \hat{H}_n \) the bootstrap distribution of
\( (\hat{\theta}^*-\hat{\theta})/\norm{\mathbf{a}_n} \) and \( H(x)=\Phi(x/\sigma) \). @thm-bs-residual gives the hypotheses of @lem-bs-quantiles. The basic interval covers \( \theta \) if and only if
\( \hat{q}(\alpha/2)-\hat{\theta}\le\hat{\theta}-\theta\le\hat{q}(1-\alpha/2)-\hat{\theta} \), that is,
\( \hat{\xi}_n(\alpha/2)\le R_n\le\hat{\xi}_n(1-\alpha/2) \) with \( \hat{\xi}_n \) the quantiles of \( \hat{H}_n \). By the lemma the
probability tends to \( (1-\alpha/2)-\alpha/2 \), the lower endpoint contributing no mass in the limit because \( H \) is
continuous. The percentile interval covers \( \theta \) if and only if \( -\hat{\xi}_n(1-\alpha/2)\le R_n\le-\hat{\xi}_n(\alpha/2) \).
Since \( H \) is symmetric, \( -\xi(1-\alpha/2)=\xi(\alpha/2) \), and the same argument applies.
:::

The same argument applies to the wild bootstrap under the conditions of @thm-bs-wild-consistency. Since \( \tau_n \)
need not converge, take \( R_n=(\hat{\theta}-\theta)/(\norm{\mathbf{a}_n}\tau_n) \) and \( H=\Phi \); dividing by the constant
\( \tau_n \) leaves the intervals unchanged. It applies also to studentized intervals whenever the standard error is consistent. None of this separates the methods from one another or from the normal interval
\( \hat{\theta}\pm z_{\alpha/2}\,\text{se} \).

## Pivots and higher-order accuracy

The case for studentizing is clearest when the root is an exact pivot.

::: {#prp-bs-parametric}
[The parametric bootstrap in the normal model]

In the normal linear model with \( \rank(\X)=p<n \), let the **parametric bootstrap** draw
\( \Y^*\sim\Normal_n(\X\hbeta,\,s^2\I) \). Then for every data set:

::: {.enumerate options="label=(\alph*)"}
1. \( T^*=(\mathbf{c}\T\hbeta^*-\mathbf{c}\T\hbeta)/(s^*\norm{\mathbf{a}})\sim t(n-p) \) under \( {\Pr}_* \), so the studentized interval is
   the \( t \) interval of @cor-opt-t, apart from Monte Carlo error;

2. \( \mathbf{c}\T\hbeta^*\sim\Normal(\mathbf{c}\T\hbeta,\,s^2\norm{\mathbf{a}}^2) \) under \( {\Pr}_* \), so the percentile and basic intervals are
   \( \hat{\theta}\pm z_{\alpha/2}\,s\norm{\mathbf{a}} \), whose coverage is
   \( \Pr\bigl(\lvert T\rvert\le z_{\alpha/2}\bigr)<1-\alpha \) with \( T\sim t(n-p) \).
:::

:::

::: {.proof}
Given the data, \( \Y^* \) follows the normal linear model with coefficient vector \( \hbeta \) and variance \( s^2 \). (a) is @cor-opt-t
applied to that model, with true value \( \mathbf{c}\T\hbeta \). (b) is @thm-opt-sampling(a) for the same model. The
interval \( \hat{\theta}\pm z_{\alpha/2}s\norm{\mathbf{a}} \) covers \( \theta \) if and only if \( \lvert T\rvert\le z_{\alpha/2} \),
where \( T=(\hat{\theta}-\theta)/(s\norm{\mathbf{a}})\sim t(n-p) \) by @cor-opt-t, and \( t(n-p) \) has heavier tails than
\( \Normal(0,1) \).
:::

The studentized root has the same distribution in both worlds, so the bootstrap reproduces it exactly. The unstudentized root
depends on \( \sigma \), and the bootstrap substitutes \( s \), losing the spread due to estimating \( \sigma \). The nonparametric
bootstrap is never exact, but the same mechanism operates: an **asymptotic pivot** is reproduced more accurately than a root whose
distribution depends strongly on the unknown law. Hall (1988) made this precise with Edgeworth expansions. Under smoothness
conditions the one-sided studentized interval has coverage error of order \( n^{-1} \), the one-sided percentile and basic
intervals of order \( n^{-1/2} \). We do not prove these rates (see Hall 1992). The price of studentizing is a standard error
recomputed in each resample, and instability when it is poorly estimated.

**BCa.** Efron (1987) corrected the percentile interval to the same order of accuracy without a standard error. The
**bias-corrected and accelerated** (BCa) interval uses the percentile interval at adjusted levels
\[
\gamma_k=\Phi\Bigl(\hat{z}_0+\frac{\hat{z}_0+z_k}{1-\hat{a}(\hat{z}_0+z_k)}\Bigr),\qquad z_1=z_{(\alpha/2)},\ z_2=z_{(1-\alpha/2)},
\]
where \( z_{(\gamma)}=\Phi^{-1}(\gamma) \), the **bias correction** \( \hat{z}_0=\Phi^{-1}\bigl({\Pr}_*(\hat{\theta}^*<\hat{\theta})\bigr) \)
measures median bias, and the **acceleration** \( \hat{a} \), the rate at which the standard error changes with \( \theta \) on a
normalizing scale, is usually estimated from the jackknife values \( \hat{\theta}_{(i)} \), computed without case \( i \):
\[
\hat{a}=\frac{\sum_id_i^3}{6\bigl(\sum_id_i^2\bigr)^{3/2}},\qquad d_i=\bar{\theta}_{(\cdot)}-\hat{\theta}_{(i)} .
\]
In regression no refitting is needed: \( \hat{\theta}_{(i)}=\hat{\theta}-a_i\hat{\varepsilon}_i/(1-h_{ii}) \) by @prp-cmp-loo(a) (@exr-bs-jackknife-loo). With \( \hat{z}_0=\hat{a}=0 \) BCa is the percentile interval. Its derivation and second-order accuracy are
in Efron (1987) and DiCiccio and Efron (1996), for resampling cases, to which the jackknife formula above belongs. Pairing it
with the residual bootstrap, as below, is a common heuristic that Efron's proof does not cover.

## An example and a simulation

::: {#exm-bs-crime-intervals}
[Intervals for the poverty coefficient]

Return to the state murder rates of @exm-proj-fwl-crime, with the District of Columbia omitted, and the coefficient
\( 0.2538 \) of poverty in the regression on poverty, single parenthood and urbanization. Seven \( 95\% \) intervals
from \( B=9999 \) resamples:

| method | interval | length |
|---|---|---|
| \( t \), classical standard error | \( (0.066,\ 0.442) \) | \( 0.376 \) |
| residual bootstrap, percentile | \( (0.064,\ 0.437) \) | \( 0.374 \) |
| residual bootstrap, basic | \( (0.070,\ 0.444) \) | \( 0.374 \) |
| residual bootstrap, studentized | \( (0.071,\ 0.449) \) | \( 0.378 \) |
| residual bootstrap, BCa | \( (0.040,\ 0.418) \) | \( 0.378 \) |
| case bootstrap, percentile | \( (0.080,\ 0.428) \) | \( 0.348 \) |
| wild bootstrap-\( t \), HC2 | \( (0.079,\ 0.428) \) | \( 0.348 \) |

The residual bootstrap used leverage-adjusted residuals. The first four intervals agree closely: with fifty states and a
well-behaved fit, the classical interval is hard to improve on. BCa is shifted down; its bias correction
\( \hat{z}_0=-0.016 \) is negligible, but the jackknife acceleration \( \hat{a}=-0.0458 \) moves both endpoints. The shift comes from that heuristic pairing, not from
BCa as such, and opposes the small upward shift of the studentized interval. The case and wild bootstraps give
shorter intervals, because the HC2 standard error \( 0.0869 \) is smaller than the classical
\( 0.0934 \) for these data. All seven
exclude zero.
:::

```{.python .run #cell-crime-intervals-intervals}
import numpy as np
import statsmodels.api as sm
from scipy import stats
data = sm.datasets.statecrime.load_pandas().data
data = data.drop(index="District of Columbia")

rng = np.random.default_rng(2303)
B = 9999
y = data["murder"].to_numpy()
X = np.column_stack([np.ones(len(y)), data[["poverty", "single", "urban"]]])
n, p = X.shape
XtX_inv = np.linalg.inv(X.T @ X)
A = XtX_inv @ X.T
a, c = A[1], np.sqrt(XtX_inv[1, 1])               # the poverty coefficient is a @ y
h = np.sum((X @ XtX_inv) * X, axis=1)
fit = X @ (A @ y)
e = y - fit
b = a @ y
s = np.sqrt(e @ e / (n - p))
q = [0.025, 0.975]
out = {"t": b + np.array([-1, 1]) * stats.t.ppf(0.975, n - p) * s * c}

r = e / np.sqrt(1 - h)                            # residual bootstrap, modified residuals
r -= r.mean()
Ystar = fit + r[rng.integers(0, n, size=(B, n))]
bstar = Ystar @ a
sstar = np.sqrt(np.sum((Ystar - (Ystar @ A.T) @ X.T) ** 2, axis=1) / (n - p))
tstar = (bstar - b) / (sstar * c)
lo, hi = np.quantile(bstar, q)
out["percentile"] = np.array([lo, hi])
out["basic"] = np.array([2 * b - hi, 2 * b - lo])
out["studentized"] = b - np.quantile(tstar, q[::-1]) * s * c

z0 = stats.norm.ppf(np.mean(bstar < b))           # BCa
loo = b - a * e / (1 - h)                         # case-jackknife values (heuristic with residual resampling)
dev = loo.mean() - loo
acc = np.sum(dev ** 3) / (6 * np.sum(dev ** 2) ** 1.5)
zs = z0 + stats.norm.ppf(q)
out["BCa"] = np.quantile(bstar, stats.norm.cdf(z0 + zs / (1 - acc * zs)))

idx = rng.integers(0, n, size=(B, n))             # case bootstrap
bcase = np.array([np.linalg.lstsq(X[i], y[i], rcond=None)[0][1] for i in idx])
out["case"] = np.quantile(bcase, q)

se2 = np.sqrt(np.sum(a ** 2 * e ** 2 / (1 - h)))  # wild bootstrap-t with HC2
Ystar = fit + (e / np.sqrt(1 - h)) * rng.choice([-1.0, 1.0], size=(B, n))
Estar = Ystar - (Ystar @ A.T) @ X.T
wstar = (Ystar @ a - b) / np.sqrt((Estar ** 2 / (1 - h)) @ a ** 2)
out["wild-t"] = b - np.quantile(wstar, q[::-1]) * se2

for m, (l, u) in out.items():
    print(f"{m:12s} ({l:.3f}, {u:.3f})   length {u - l:.3f}")
```

One data set cannot show which interval is right. [Figure 23.3.1](#fig-bs-coverage) reports coverages for the slope of a
line, from \( 2000 \) simulated data sets with \( B=999 \) resamples each, and Monte Carlo standard error about
\( 0.005 \). The fixed design is skewed, with largest leverage \( 0.433 \) at \( n=20 \) and \( 0.309 \) at
\( n=40 \). The errors are either centred exponential with constant variance, or normal with standard deviation proportional to
\( x \), so that the largest variances sit at the high-leverage end.

::: {when-format="html"}
![**Figure 23.3.1.** Coverage of eight nominal 95% intervals for a slope. Circles: skewed errors, constant variance. Squares:
normal errors with variance increasing in \( x \). The band is two Monte Carlo standard errors around 0.95.](coverage.svg){#fig-bs-coverage width=100%}
:::

::: {when-format="pdf"}
![Coverage of eight nominal 95% intervals for a slope. Circles: skewed errors, constant variance. Squares:
normal errors with variance increasing in \( x \). The band is two Monte Carlo standard errors around 0.95.](coverage.pdf){width=100%}
:::

Three lessons stand out.

1. *With constant variance, studentization is what matters.* At \( n=20 \) the classical \( t \) interval (\( 0.954 \))
   and the studentized residual bootstrap (\( 0.953 \)) are on target despite skewed errors. The percentile
   (\( 0.938 \)) and BCa (\( 0.930 \)) intervals are a little short, as @prp-bs-parametric predicts for intervals that
   ignore the variability of \( s \). The robust methods (case \( 0.912 \), wild-\( t \) \( 0.918 \)) pay for
   estimating a variance that did not need estimating.

2. *With heteroscedastic errors the classical \( t \) interval and the four residual-bootstrap intervals all fail, and more
   data do not help:* all five cover
   between \( 0.71 \) and \( 0.78 \) at both sample sizes (@exr-bs-heteroscedastic-residual).

3. *The robust methods improve with \( n \) but start below nominal.* The case bootstrap, the HC2 \( t \) interval and the wild
   bootstrap-\( t \) cover \( 0.899 \), \( 0.852 \) and \( 0.868 \) at \( n=20 \), and \( 0.914 \),
   \( 0.893 \) and \( 0.922 \) at \( n=40 \). A high-leverage point with a large error variance is hard for all of
   them, since one residual carries much of the information about the variance that matters (@exm-bs-engel-wild).

```{.python .run #cell-coverage-small}
import numpy as np
from scipy import stats
level = 0.95
lo_q, hi_q = (1 - level) / 2, (1 + level) / 2
beta = np.array([1.0, 1.0])
methods = ["t", "percentile", "basic", "studentized", "BCa", "case", "HC2 t", "wild-t"]

def make_design(n):
    """A skewed fixed design and the quantities every data set reuses."""
    x = np.exp(0.5 * stats.norm.ppf((np.arange(1, n + 1) - 0.5) / n))
    X = np.column_stack([np.ones(n), x])
    XtX_inv = np.linalg.inv(X.T @ X)
    A = XtX_inv @ X.T                             # beta_hat = A y
    h = np.sum((X @ XtX_inv) * X, axis=1)         # leverages
    sd = x / np.sqrt(np.mean(x ** 2))             # heteroscedastic sd, mean square one
    return dict(n=n, x=x, X=X, A=A, h=h, a=A[1], c=np.sqrt(XtX_inv[1, 1]), sd=sd)

def errors(rng, law, d):
    if law == "skewed":                           # centred exponential, constant variance
        return rng.exponential(size=d["n"]) - 1.0
    return d["sd"] * rng.normal(size=d["n"])      # normal, sd proportional to x

def one_data_set(rng, law, d, B):
    """Return, for each method, whether its 95% interval covers the true slope."""
    n, x, X, A, h, a, c = (d[k] for k in ("n", "x", "X", "A", "h", "a", "c"))
    y = X @ beta + errors(rng, law, d)
    b = a @ y
    e = y - X @ (A @ y)
    s = np.sqrt(e @ e / (n - 2))
    cover = {"t": abs(b - beta[1]) <= stats.t.ppf(hi_q, n - 2) * s * c}
    # residual bootstrap from leverage-adjusted, centred residuals
    r = e / np.sqrt(1 - h)
    r -= r.mean()
    Ystar = X @ (A @ y) + r[rng.integers(0, n, size=(B, n))]
    bstar = Ystar @ a
    Estar = Ystar - (Ystar @ A.T) @ X.T
    tstar = (bstar - b) / (np.sqrt(np.sum(Estar ** 2, axis=1) / (n - 2)) * c)
    q_lo, q_hi = np.quantile(bstar, [lo_q, hi_q])
    t_lo, t_hi = np.quantile(tstar, [lo_q, hi_q])
    cover["percentile"] = q_lo <= beta[1] <= q_hi
    cover["basic"] = 2 * b - q_hi <= beta[1] <= 2 * b - q_lo
    cover["studentized"] = b - t_hi * s * c <= beta[1] <= b - t_lo * s * c
    # BCa: residual-bootstrap bias correction, case-jackknife acceleration (a heuristic pairing)
    z0 = stats.norm.ppf(np.clip(np.mean(bstar < b), 1 / B, 1 - 1 / B))
    loo = b - a * e / (1 - h)                     # slopes with one case deleted
    dev = loo.mean() - loo
    acc = np.sum(dev ** 3) / (6 * np.sum(dev ** 2) ** 1.5)
    zs = z0 + stats.norm.ppf([lo_q, hi_q])
    bca_lo, bca_hi = np.quantile(bstar, stats.norm.cdf(z0 + zs / (1 - acc * zs)))
    cover["BCa"] = bca_lo <= beta[1] <= bca_hi
    # case bootstrap percentile
    idx = rng.integers(0, n, size=(B, n))
    xs, ys = x[idx], y[idx]
    xc = xs - xs.mean(axis=1, keepdims=True)
    sxx = np.sum(xc ** 2, axis=1)
    ok = sxx > 1e-12                              # drop the (rare) resamples with one distinct x
    bcase = np.sum(xc[ok] * ys[ok], axis=1) / sxx[ok]
    cl, ch = np.quantile(bcase, [lo_q, hi_q])
    cover["case"] = cl <= beta[1] <= ch
    # HC2 t interval, and the wild bootstrap-t studentized by HC2
    se2 = np.sqrt(np.sum(a ** 2 * e ** 2 / (1 - h)))
    cover["HC2 t"] = abs(b - beta[1]) <= stats.t.ppf(hi_q, n - 2) * se2
    v = rng.choice([-1.0, 1.0], size=(B, n))
    Ystar = X @ (A @ y) + (e / np.sqrt(1 - h)) * v
    Estar = Ystar - (Ystar @ A.T) @ X.T
    wstar = (Ystar @ a - b) / np.sqrt((Estar ** 2 / (1 - h)) @ a ** 2)
    w_lo, w_hi = np.quantile(wstar, [lo_q, hi_q])
    cover["wild-t"] = b - w_hi * se2 <= beta[1] <= b - w_lo * se2
    return cover

def coverage(law, n, reps, B, seed):
    rng = np.random.default_rng(seed)
    d = make_design(n)
    hits = {m: 0 for m in methods}
    for _ in range(reps):
        for m, ok in one_data_set(rng, law, d, B).items():
            hits[m] += ok
    return {m: hits[m] / reps for m in methods}

for law in ["skewed", "heteroscedastic"]:
    cov = coverage(law, n=20, reps=200, B=499, seed=2331)
    print(f"{law:16s}" + "  ".join(f"{m} {v:.2f}" for m, v in cov.items()))
```

With \( 200 \) data sets the cell's coverages have Monte Carlo standard errors of about
\( \sqrt{0.95\times0.05/200}=0.015 \).

## Bootstrap tests of a linear hypothesis

A test of \( H:\bLambda\T\bbeta=\mathbf{d} \) can take its critical value from the bootstrap instead of the \( F \) table. In each
residual-bootstrap resample from the full model, compute \( F^*_H \) of @eq-glh-general-F for the hypothesis
\( \bLambda\T\bbeta=\bLambda\T\hbeta \), which is true in the bootstrap world, and report the proportion of resamples with
\( F^*_H\ge F_H \). This is the testing counterpart of the studentized interval (@prp-ci-duality). Alternatively, resample the
residuals of the restricted fit, so that the null hypothesis holds in the bootstrap world; that is the bootstrap analogue of the
Freedman–Lane scheme of [Section 23.5](05-freedman-lane.html).

## Exercises

### A. Check your understanding

::: {#exr-bs-basic-vs-percentile}
[A1]

Show that the percentile and basic intervals coincide when the bootstrap distribution of \( \hat{\theta}^* \) is symmetric about
\( \hat{\theta} \). With \( \hat{\theta}=2 \), \( \hat{q}(0.025)=1.5 \) and \( \hat{q}(0.975)=3.5 \), compute both intervals and explain why they
lean in opposite directions.
:::

::: {.solution}
Symmetry about \( \hat{\theta} \) means \( \hat{q}(\gamma)-\hat{\theta}=\hat{\theta}-\hat{q}(1-\gamma) \), so
\( 2\hat{\theta}-\hat{q}(1-\alpha/2)=\hat{q}(\alpha/2) \), and similarly for the other endpoint. In the example the percentile interval is
\( [1.5,3.5] \) and the basic interval is \( [4-3.5,\ 4-1.5]=[0.5,2.5] \). The bootstrap distribution has a long right tail.
The percentile interval copies it. The basic interval reflects it, reasoning that if \( \hat{\theta} \) tends to fall far above
\( \theta \) with some probability, then \( \theta \) may lie far below \( \hat{\theta} \).
:::

### B. Practice

::: {#exr-bs-percentile-transform}
[B1]

Suppose there is a strictly increasing \( \phi \) such that \( \phi(\hat{\theta})-\phi(\theta) \) has a continuous distribution \( K \),
symmetric about zero and the same for every \( \theta \), and suppose the bootstrap reproduces it exactly:
\( \phi(\hat{\theta}^*)-\phi(\hat{\theta}) \) has distribution \( K \) under \( {\Pr}_* \). Show that the percentile interval has coverage exactly
\( 1-\alpha \).
:::

::: {.solution}
Let \( k_\gamma \) be the \( \gamma \) quantile of \( K \), so \( k_{\alpha/2}=-k_{1-\alpha/2} \). Under \( {\Pr}_* \),
\( \phi(\hat{\theta}^*)=\phi(\hat{\theta})+\) a \( K \)-variable, so the \( \gamma \) quantile of \( \hat{\theta}^* \) is
\( \hat{q}(\gamma)=\phi^{-1}\bigl(\phi(\hat{\theta})+k_\gamma\bigr) \). Hence \( \hat{q}(\alpha/2)\le\theta\le\hat{q}(1-\alpha/2) \) if and only if
\( k_{\alpha/2}\le\phi(\theta)-\phi(\hat{\theta})\le k_{1-\alpha/2} \), that is,
\( -k_{1-\alpha/2}\le\phi(\hat{\theta})-\phi(\theta)\le-k_{\alpha/2} \). By symmetry this is
\( k_{\alpha/2}\le\phi(\hat{\theta})-\phi(\theta)\le k_{1-\alpha/2} \), which has probability \( 1-\alpha \).
:::

::: {#exr-bs-jackknife-loo}
[B2]

Show that the jackknife values used in the BCa acceleration satisfy \( \hat{\theta}_{(i)}=\hat{\theta}-a_i\hat{\varepsilon}_i/(1-h_{ii}) \), where
\( \mathbf{a}=\X(\X\T\X)^{-1}\mathbf{c} \). Deduce that the jackknife variance estimate
\( \frac{n-1}{n}\sum_i(\hat{\theta}_{(i)}-\bar{\theta}_{(\cdot)})^2 \) is close to the HC3 variance \( \sum_ia_i^2\hat{\varepsilon}_i^2/(1-h_{ii})^2 \).
([Chapter 21](../ch21-nonnormality-heteroscedasticity-serial/index.html) shows that the uncentred sum
\( \sum_i(\hat{\theta}_{(i)}-\hat{\theta})^2 \) equals it exactly.)
:::

### C. Going deeper

::: {#exr-bs-iterated}
[C1]

The coverage of a bootstrap interval can itself be estimated by a second level of bootstrapping: treat each first-level resample
as data, and ask how often its interval covers \( \hat{\theta} \). Describe how this can be used to recalibrate the nominal level, and
count the least squares fits needed with \( B_1=B_2=999 \). How does the design-fixed residual bootstrap reduce the cost?
:::
