# The likelihood ratio test and its relatives

The \( F \) test was motivated in [Section 11.1](01-reduced-models.html) by a heuristic: compare the
extra sum of squares with the noise level. General statistical theory offers other routes to a test,
chiefly the likelihood ratio, the Wald statistic and the score statistic. In the normal linear model
all three routes lead back to \( F \). They differ only when they are calibrated by their common large-sample
\( \chi^2 \) approximation instead of the exact distribution, and then they can disagree about the
same data. This section derives all three, shows how different the approximate versions can be, and
then proves that the \( F \) test is the best test among those that respect the symmetries of
[Section 11.2](02-comparing-projections.html).

## The likelihood ratio

For a hypothesis \( H_0 \) about the parameter of a family of densities \( L(\boldsymbol{\theta};\y) \), the
**likelihood ratio** is
\[
\Lambda(\y)=\frac{\sup_{\boldsymbol{\theta}\in H_0}L(\boldsymbol{\theta};\y)}{\sup_{\boldsymbol{\theta}}L(\boldsymbol{\theta};\y)} ,
\]
and the likelihood ratio test rejects when \( \Lambda \) is small, that is, when the best explanation the
hypothesis can offer is much worse than the best explanation overall. In the linear model both suprema
are available in closed form from @thm-opt-mle, and @exr-opt-lr-nested already carried out the
computation for reduced models. The theorem records it for testable hypotheses, including the
calibration.

::: {#thm-glh-lrt}
[The \( F \) test is the likelihood ratio test]

Assume @eq-opt-normal-model with \( r<n \), and let \( H:\bLambda\T\bbeta=\mathbf{d} \) be testable with
\( q=\rank(\bLambda) \). Let \( \text{SSE}_H=\min_{\bLambda\T\bb=\mathbf{d}}\norm{\y-\X\bb}^2 \) and suppose
\( \text{SSE}>0 \), which holds with probability one. Then:

::: {.enumerate options="label=(\alph*)"}
1. under \( H \) the likelihood is maximized at any minimizer \( \hbeta_H \) of \( \norm{\y-\X\bb}^2 \) subject
   to \( \bLambda\T\bb=\mathbf{d} \) (for \( \bLambda \) of full column rank, the estimate of
   @prp-glh-restricted-ls; compare @thm-ss-restricted), with \( \tilde{\sigma}^2_H=\text{SSE}_H/n \);

2. the likelihood ratio is
   \[
\Lambda=\Bigl(\frac{\text{SSE}}{\text{SSE}_H}\Bigr)^{n/2}=\Bigl(1+\frac{q\,F_H}{n-r}\Bigr)^{-n/2},
\]{#eq-glh-lr}

   a strictly decreasing function of the statistic \( F_H \) of @eq-glh-general-F;

3. the likelihood ratio test of size \( \alpha \) rejects exactly when \( F_H>F_\alpha(q,n-r) \).
:::

The same holds for a reduced model \( \E(\Y)\in\C(\X_0) \), with \( q=r-r_0 \) and the statistic of @thm-glh-f-test.
:::

::: {.proof}
(a) By @thm-glh-general-f(b), \( H \) holds iff \( \E(\Y)\in\X\bb_0+\mathcal S_0 \). Maximizing the likelihood
over \( H \) is therefore maximizing it in the linear model \( \Y-\X\bb_0\sim\Normal_n(\mathbf{m},\sigma^2\I) \),
\( \mathbf{m}\in\mathcal S_0 \). By @thm-opt-mle applied to that model (with any matrix whose columns span
\( \mathcal S_0 \)), the maximizing mean is the least squares fit and the maximizing variance is its
residual sum of squares divided by \( n \). That fit is \( \bP_{\mathcal S_0}(\y-\X\bb_0) \), so the
maximizing mean of \( \Y \) is \( \X\bb_0+\bP_{\mathcal S_0}(\y-\X\bb_0) \), the point of
\( \X\bb_0+\mathcal S_0 \) nearest to \( \y \). As shown in the proof of @thm-glh-general-f(b), no rank
condition on \( \bLambda \) is needed for this. That point is \( \X\bb \) exactly for the minimizers \( \bb \) of
\( \norm{\y-\X\bb}^2 \) subject to \( \bLambda\T\bb=\mathbf{d} \), so it is \( \X\hbeta_H \), and its residual sum of
squares is \( \text{SSE}_H \). This needs \( \text{SSE}_H>0 \), which follows from
\( \text{SSE}_H\ge\text{SSE}>0 \).

(b) By @thm-opt-mle the two maxima are \( (2\pi e\,\text{SSE}/n)^{-n/2} \) and
\( (2\pi e\,\text{SSE}_H/n)^{-n/2} \). Their ratio is \( (\text{SSE}/\text{SSE}_H)^{n/2} \), and
\( \text{SSE}_H=\text{SSE}+\text{SS}_H=\text{SSE}\{1+qF_H/(n-r)\} \) by @thm-glh-general-f(b).

(c) By (b), \( \Lambda\le k \) iff \( F_H\ge c \) for a constant \( c \) determined by \( k \). By
@thm-glh-general-f(c) the null distribution of \( F_H \) is \( F(q,n-r) \) at every point of \( H \), so size
\( \alpha \) requires \( c=F_\alpha(q,n-r) \).
:::

The large-sample theory of likelihood ratios (Wilks's theorem) would calibrate
\( -2\log\Lambda=n\log\{1+qF_H/(n-r)\} \) by the \( \chi^2(q) \) distribution. That calibration is only
approximate. In the normal linear model the exact distribution is available through \( F_H \) and
should be used.

## Wald and score statistics

Two other large-sample recipes are in wide use. The **Wald statistic** measures the estimated departure
\( \bLambda\T\hbeta-\mathbf{d} \) against its estimated covariance, both computed at the unrestricted
maximum likelihood estimate. The **score statistic**, also called the Lagrange multiplier statistic,
measures the gradient of the log-likelihood at the *restricted* maximum likelihood estimate against the
information there. In the linear model, with \( \bLambda \) of full column rank \( q \) and the maximum
likelihood variance estimates \( \hat{\sigma}^2=\text{SSE}/n \) and \( \tilde{\sigma}^2_H=\text{SSE}_H/n \), they are
\[
\text{Wald}=\frac{(\bLambda\T\hbeta-\mathbf{d})\T(\bLambda\T\G\bLambda)^{-1}(\bLambda\T\hbeta-\mathbf{d})}{\hat{\sigma}^2}
\]
and
\[
\text{Score}=\mathbf{U}\T\Bigl(\frac{\X\T\X}{\tilde{\sigma}^2_H}\Bigr)\ginv\mathbf{U}=\tilde{\sigma}^2_H\,\mathbf{U}\T(\X\T\X)\ginv\mathbf{U},
\]
where \( \mathbf{U}=\X\T(\y-\X\hbeta_H)/\tilde{\sigma}^2_H \). Here \( \mathbf{U} \) is the gradient of the log-likelihood @eq-opt-normal-model with respect to \( \bbeta \) at
\( (\hbeta_H,\tilde{\sigma}^2_H) \), and \( \X\T\X/\tilde{\sigma}^2_H \) is the corresponding block of the
Fisher information. The gradient with respect to \( \sigma^2 \) vanishes at the restricted maximum, and the
information matrix is block diagonal, because \( \E\{\X\T(\Y-\X\bbeta)\}=\bzero \). So the \( \sigma^2 \)
coordinates contribute nothing, and the formula is the full score statistic. The generalized inverse
is harmless because \( \mathbf{U}\in\C(\X\T)=\C(\X\T\X) \).

::: {#prp-glh-trinity}
[Wald, likelihood ratio and score]

In the setting above, let \( u=\text{SS}_H/\text{SSE}=qF_H/(n-r) \). Then
\[
\text{Wald}=n\,u,\qquad -2\log\Lambda=n\log(1+u),\qquad \text{Score}=\frac{n\,u}{1+u}.
\]
Each is a strictly increasing function of \( F_H \), so the three tests, calibrated exactly, coincide with
the \( F \) test. For every sample with \( \text{SS}_H>0 \),
\[
\text{Wald}>-2\log\Lambda>\text{Score}.
\]
If each is compared with the upper \( \alpha \) point \( c \) of \( \chi^2(q) \), the exact sizes are
\[
\Pr\Bigl\{F(q,n-r)>\frac{n-r}{q}\,u_c\Bigr\},\qquad
u_c=\frac cn,\quad e^{c/n}-1,\quad \frac{c}{n-c}\ \ (n>c)
\]
for the Wald, likelihood ratio and score tests respectively.
:::

::: {.proof}
With \( \hat{\sigma}^2=\text{SSE}/n \), \( \text{Wald}=n\,\text{SS}_H/\text{SSE}=nu \) by @eq-glh-general-F. For the
score, \( \X\T\y=\X\T\X\hbeta \), so \( \X\T(\y-\X\hbeta_H)=\X\T\X(\hbeta-\hbeta_H) \) and
\[
\mathbf{U}\T(\X\T\X)\ginv\mathbf{U}\cdot\tilde{\sigma}^2_H
=\frac{(\hbeta-\hbeta_H)\T\X\T\X(\X\T\X)\ginv\X\T\X(\hbeta-\hbeta_H)}{\tilde{\sigma}^2_H}
=\frac{\norm{\X(\hbeta-\hbeta_H)}^2}{\tilde{\sigma}^2_H}.
\]
By the proof of @prp-glh-restricted-ls, \( \norm{\X(\hbeta-\hbeta_H)}^2=\text{SS}_H \), so
\( \text{Score}=n\,\text{SS}_H/\text{SSE}_H=nu/(1+u) \). The likelihood ratio is @eq-glh-lr. All three
are increasing in \( u \), which is proportional to \( F_H \). For \( u>0 \),
\( \log(1+u)=\int_0^u(1+t)^{-1}\,dt \) lies strictly between \( u/(1+u) \) and \( u \), since the integrand lies
strictly between \( (1+u)^{-1} \) and \( 1 \) on \( (0,u) \). Finally, each statistic exceeds \( c \) iff \( u>u_c \),
with \( u_c \) as stated (for the score, \( nu/(1+u)>c \) iff \( u(n-c)>c \), and it never exceeds \( c \) when
\( n\le c \)), and \( u>u_c \) iff \( F_H>(n-r)u_c/q \).
:::

The inequality \( \text{Wald}\ge\text{LR}\ge\text{Score} \) is a general feature of linear
restrictions in normal regression, and it has a practical consequence: with \( \chi^2 \) critical values,
the Wald test rejects most often and the score test least often, *on the same data*.

::: {#exm-glh-trinity}
[Three verdicts from one data set]

For the region test of @exm-glh-region, \( n=50 \), \( q=3 \) and \( n-r=43 \). The three
statistics are

| Statistic | Value | p-value from \( \chi^2(3) \) |
|:---|:---:|:---:|
| Wald | \( 8.512 \) | \( 0.0365 \) |
| Likelihood ratio, \( -2\log\Lambda \) | \( 7.861 \) | \( 0.0490 \) |
| Score | \( 7.274 \) | \( 0.0637 \) |
| \( F=2.440 \), exact | | \( 0.0773 \) |

With the \( 5\% \) critical value \( 7.815 \) of \( \chi^2(3) \), the Wald test rejects, the
likelihood ratio test rejects barely, and the score test does not. The exact test, which all three
approximate, does not reject. The disagreement is not a paradox. The approximate tests are different
tests, and their true sizes differ from \( 0.05 \). With \( q=3 \) and \( r=7 \), the formulas of
@prp-glh-trinity give

| \( n \) | \( 15 \) | \( 30 \) | \( 50 \) | \( 100 \) | \( 400 \) |
|:---|:---:|:---:|:---:|:---:|:---:|
| Wald | \( 0.315 \) | \( 0.143 \) | \( 0.097 \) | \( 0.071 \) | \( 0.055 \) |
| Likelihood ratio | \( 0.221 \) | \( 0.106 \) | \( 0.079 \) | \( 0.063 \) | \( 0.053 \) |
| Score | \( 0.102 \) | \( 0.069 \) | \( 0.060 \) | \( 0.055 \) | \( 0.051 \) |

At the sample size of the example the nominal \( 5\% \) Wald test has true size close to \( 10\% \).
All three sizes approach \( 0.05 \) as \( n \) grows, but slowly when \( r \) is not small compared with \( n \).
:::

```{.python .run #cell-likelihood-ratio-trinity}
import numpy as np
import statsmodels.api as sm
from scipy import stats

data = sm.datasets.statecrime.load_pandas().data
data.index = data.index.str.strip()
data = data.drop(index="District of Columbia")
codes = "SWWSWWNSSSWWMMMMSSNSNMMSMWMWNNWNSMMSWNNSMSSWNSWSMW"
region = np.array(list(codes))
y = data["murder"].to_numpy()
n = len(y)
covariates = np.column_stack([data["poverty"], data["single"], data["urban"]])
X0 = np.column_stack([np.ones(n), covariates])
D = np.column_stack([(region == g).astype(float) for g in "NSW"])
X = np.column_stack([X0, D])
r, q = 7, 3

def fit(Z):
    b, *_ = np.linalg.lstsq(Z, y, rcond=None)
    return b, np.sum((y - Z @ b) ** 2)

b, sse = fit(X)
b0, sse0 = fit(X0)
F = (sse0 - sse) / q / (sse / (n - r))

u = (sse0 - sse) / sse                          # the ratio behind every statistic
wald = n * u                                    # n (SSE0 - SSE) / SSE
lr = n * np.log1p(u)                            # n log(SSE0 / SSE)
score = n * u / (1 + u)                         # n (SSE0 - SSE) / SSE0
crit = stats.chi2.ppf(0.95, q)
for name, stat in (("Wald", wald), ("LR", lr), ("score", score)):
    print(f"{name:5s} {stat:.3f}  chi2 p = {stats.chi2.sf(stat, q):.4f}")
print(f"F     {F:.3f}  exact p = {stats.f.sf(F, q, n - r):.4f}")
```

The lesson for other models is double. In the normal linear model none of the three approximations is
needed, because the exact \( F \) test is available. In models without an exact theory,
generalized linear models among them (Chapter 34 onwards), the three statistics are the standard tools. Their
disagreement in small samples is a warning that the \( \chi^2 \) approximation is poor. The linear model
suggests a remedy: dividing by \( q \) and referring to an \( F \) distribution with finite denominator degrees
of freedom often calibrates better.

## Invariance and optimality

Is the \( F \) test a *good* test? For a single constraint and a one-sided alternative, the \( t \) test of
[Section 11.6](06-coefficients.html) has strong optimality properties. For \( q\ge2 \) no test is
uniformly most powerful: the most powerful test against alternatives in one direction of the test
space rejects for large values of the estimated departure *in that direction*, and it is a different
test for each direction (@exr-glh-no-ump). Some restriction on the class of tests is needed. The
natural one is invariance. By @prp-glh-invariance, the problem does not change if the data are rescaled,
shifted within the reduced model, or rotated within the test space and within the residual space. A
test whose verdict changes under such operations would be favouring particular directions for no
reason given by the problem. We show that among tests that do not, the \( F \) test is the most powerful,
simultaneously against every alternative.

We need one property of the noncentral \( F \) family.

::: {#lem-glh-mlr}
[Monotone likelihood ratio of the noncentral \( F \)]

For \( q,\nu\ge1 \) let \( f_\gamma \) be the density of \( F(q,\nu,\gamma) \). For every \( \gamma>0 \), the ratio
\( f_\gamma(x)/f_0(x) \) is a strictly increasing function of \( x>0 \).
:::

::: {.proof}
Let \( V\sim\chi^2(\nu) \), and for \( k=0,1,\dots \) let \( h_k \) be the density of \( (U_k/q)/(V/\nu) \), where
\( U_k\sim\chi^2(q+2k) \) is independent of \( V \). By @exr-qf-beta,
\( B=U_k/(U_k+V)\sim\mathrm{Beta}(a_k,b) \) with \( a_k=q/2+k \) and \( b=\nu/2 \). The map
\( B\mapsto x=\nu B/\{q(1-B)\} \) is increasing, with inverse \( B=t(x):=qx/(\nu+qx) \), and
\( 1-t(x)=\nu/(\nu+qx) \). Changing variables,
\[
h_k(x)=\frac{t(x)^{a_k-1}\{1-t(x)\}^{b-1}}{\mathrm{B}(a_k,b)}\cdot\frac{q\nu}{(\nu+qx)^2},
\qquad\text{so}\qquad
\frac{h_k(x)}{h_0(x)}=\frac{\mathrm{B}(q/2,\nu/2)}{\mathrm{B}(q/2+k,\nu/2)}\,t(x)^k=c_k\,t(x)^k ,
\]
where \( \mathrm B \) is the beta function and \( c_k>0 \). By @thm-qf-ncchisq(d), \( U\sim\chi^2(q,\gamma) \) has the
law of \( U_K \) with \( K\sim\text{Poisson}(\gamma/2) \) independent of the \( U_k \) and \( V \). Conditioning on \( K \),
\( F=(U/q)/(V/\nu) \) has density \( f_\gamma(x)=\sum_{k\ge0}p_kh_k(x) \) with
\( p_k=e^{-\gamma/2}(\gamma/2)^k/k! \), and \( f_0=h_0 \). Hence
\[
\frac{f_\gamma(x)}{f_0(x)}=\sum_{k\ge0}p_kc_k\,t(x)^k .
\]
Every coefficient is positive when \( \gamma>0 \), so the series is strictly increasing in \( t \), and \( t(x) \) is
strictly increasing in \( x \).
:::

A test is a measurable function \( \phi \) of the data with values in \( [0,1] \), the probability of
rejecting. It is **invariant** if \( \phi(g(\y))=\phi(\y) \) for all \( \y \) and all \( g \) in the set \( \mathcal G \) of
maps
\[
g(\y)=a\mathbf{H}\y+\bv,\qquad a>0,\quad \bv\in\C(\X_0),\quad
\mathbf{H}\text{ orthogonal with }\mathbf{H}\C(\X_0)=\C(\X_0),\ \mathbf{H}\C(\X)=\C(\X).
\]
For a testable hypothesis with \( \mathbf{d}\ne\bzero \), the same applies to the data \( \y-\X\bb_0 \), with
\( \C(\X_0) \) replaced by \( \mathcal S_0 \).

::: {#thm-glh-ump-invariant}
[The \( F \) test is uniformly most powerful invariant]

Assume @eq-opt-normal-model and the setting of @def-glh-reduced-model. Let \( \phi \) be an invariant test
with \( \E\phi(\Y)\le\alpha \) at every point of \( H_0 \). Then at every point \( (\bbeta,\sigma^2) \) of the full model,
\[
\E_{\bbeta,\sigma^2}\,\phi(\Y)\le\Pr\{F(q,n-r,\gamma)>F_\alpha(q,n-r)\},
\]
with \( \gamma \) given by @eq-glh-noncentrality. The right side is the power of the \( F \) test, which is itself
invariant and of size \( \alpha \).
:::

::: {.proof}
Use the canonical coordinates \( \bz_k=\Q_k\T\y \) of @prp-glh-canonical, and let
\( \mathcal Y=\{\y:\bz_2\ne\bzero\} \), a set of probability one on which \( F \) is defined.

*Step 1: an invariant test depends on \( \y \) only through \( F \).* Let \( \y,\y'\in\mathcal Y \) with
\( F(\y)=F(\y') \). Put \( a=\norm{\bz_2'}/\norm{\bz_2}>0 \). Equality of the \( F \) values gives
\( \norm{\bz_1'}=a\norm{\bz_1} \). Choose orthogonal matrices \( \mathbf{O}_1 \) (\( q\times q \)) and \( \mathbf{O}_2 \)
(\( (n-r)\times(n-r) \)) with \( \mathbf{O}_1(a\bz_1)=\bz_1' \) and \( \mathbf{O}_2(a\bz_2)=\bz_2' \). This is possible because the
vectors in each pair have equal lengths, and if \( \bz_1=\bzero \) then \( \bz_1'=\bzero \) and \( \mathbf{O}_1=\I \) will do. Put
\[
\mathbf{H}=\Q_0\Q_0\T+\Q_1\mathbf{O}_1\Q_1\T+\Q_2\mathbf{O}_2\Q_2\T,\qquad \bv=\Q_0(\bz_0'-a\bz_0).
\]
\( \mathbf{H} \) is orthogonal, maps \( \C(\X_0)=\C(\Q_0) \) and \( \C(\X)=\C([\Q_0,\Q_1]) \) onto themselves, and
\( \bv\in\C(\X_0) \). The canonical coordinates of \( a\mathbf{H}\y+\bv \) are \( a\bz_0+\bz_0'-a\bz_0=\bz_0' \),
\( a\mathbf{O}_1\bz_1=\bz_1' \) and \( a\mathbf{O}_2\bz_2=\bz_2' \), so \( a\mathbf{H}\y+\bv=\y' \) and \( \phi(\y')=\phi(\y) \). Now fix, for
each \( f\ge0 \), the point \( \y_f=\sqrt{qf/(n-r)}\,\Q_1\mathbf{e}_1+\Q_2\mathbf{e}_1 \), where \( \mathbf{e}_1 \) denotes a first coordinate
vector, and put \( \psi(f)=\phi(\y_f) \). Then \( F(\y_f)=f \), so \( \phi(\y)=\psi(F(\y)) \) for every \( \y\in\mathcal Y \).
The map \( f\mapsto\y_f \) is continuous, so \( \psi \) is measurable.

*Step 2: reduction to one dimension.* By Step 1 and @thm-glh-f-test(b),
\( \E\phi(\Y)=\E\psi(F)=\int_0^\infty\psi(x)f_\gamma(x)\,dx \), with \( f_\gamma \) the density of \( F(q,n-r,\gamma) \). The
hypothesis is \( \gamma=0 \) (@thm-glh-f-test(c)), so the size condition says \( \int\psi f_0\le\alpha \).

*Step 3: Neyman–Pearson.* Fix \( \gamma>0 \), let \( c=F_\alpha(q,n-r) \), \( \phi^*(x)=1\{x>c\} \) and
\( k=f_\gamma(c)/f_0(c) \). By @lem-glh-mlr, \( f_\gamma(x)-kf_0(x) \) is positive for \( x>c \) and negative for
\( x<c \). Since \( \phi^*-\psi\ge0 \) where \( x>c \) and \( \phi^*-\psi\le0 \) where \( x<c \), the product
\( (\phi^*-\psi)(f_\gamma-kf_0) \) is nonnegative everywhere. Integrating,
\[
\int(\phi^*-\psi)f_\gamma\ \ge\ k\int(\phi^*-\psi)f_0=k\Bigl(\alpha-\int\psi f_0\Bigr)\ \ge\ 0 .
\]
So \( \E\phi(\Y)=\int\psi f_\gamma\le\int\phi^*f_\gamma \), which is the power of the \( F \) test. For
\( \gamma=0 \) the inequality is the size condition. The \( F \) test is invariant by @prp-glh-invariance.
:::

The theorem says more than "the \( F \) test is good". It says that once one agrees that no direction
in the test space and no scale deserves special treatment, the choice of test is forced: every
invariant test is a function of \( F \) (Step 1), and among functions of \( F \) the upper tail is best
against every alternative at once. The likelihood ratio principle arrives at the same test without any
appeal to optimality, which is reassuring but not a proof. In other models the likelihood ratio test
need not be optimal in any sense.

::: {.remark}
[What invariance costs]

Invariance is a restriction, and the \( F \) test is not the most powerful test against every
alternative. If one contrast in the test space is of special interest, the \( t \) test for that contrast
is more powerful against departures in its direction and less powerful against others
([Section 11.6](06-coefficients.html)). A Bayesian analysis with a prior concentrated on some
directions behaves similarly. Invariance expresses the absence of such preferences. It is appropriate
for an omnibus question like "do the regions differ?", and inappropriate when the scientific question
singles out a direction.
:::

## Exercises

### A. Check your understanding

::: {#exr-glh-lr-numbers}
[A1]

From \( F=2.440 \), \( q=3 \), \( n=50 \) and \( r=7 \), compute the three statistics of
@prp-glh-trinity and check them against @exm-glh-trinity.
:::

### B. Practice

::: {#exr-glh-wilks}
[B1]

Fix \( q \) and \( r \). Show that under \( H \), as \( n\to\infty \), \( qF_H \) converges in distribution to
\( \chi^2(q) \), and that the differences \( \text{Wald}-qF_H \), \( -2\log\Lambda-qF_H \) and \( \text{Score}-qF_H \) converge
to zero in probability. Deduce that the sizes in @exm-glh-trinity tend to \( 0.05 \).
:::

::: {.solution}
Write \( qF_H=(\text{SS}_H/\sigma^2)/(s^2/\sigma^2) \). The numerator is \( \chi^2(q) \) for every \( n \), and
\( s^2/\sigma^2\to1 \) in probability (@exr-glh-known-sigma), so \( qF_H\to\chi^2(q) \) in distribution. With
\( u=qF_H/(n-r) \), \( \text{Wald}-qF_H=qF_H\{n/(n-r)-1\}=qF_H\,r/(n-r)\to0 \) in probability. For the others,
\( nu=O_p(1) \), so \( u\to0 \), and \( n\log(1+u)-nu=nu\{\log(1+u)/u-1\} \) and \( nu/(1+u)-nu=-nu^2/(1+u) \) both
tend to zero. Each statistic therefore converges in distribution to \( \chi^2(q) \), whose upper \( 5\% \) point
is the critical value used, so the sizes converge to \( 0.05 \).
:::

::: {#exr-glh-nr2}
[B2]

Let \( \hat{\be}_H=\y-\X\hbeta_H \) be the residual vector of the restricted fit. Show that the score
statistic equals \( n\norm{\M\hat{\be}_H}^2/\norm{\hat{\be}_H}^2 \), which is \( n \) times the uncentred \( R^2 \) from
regressing \( \hat{\be}_H \) on \( \X \). Explain why this "\( nR^2 \)" recipe needs only the restricted fit.
:::

::: {.solution}
\( \M\hat{\be}_H=\M\y-\X\hbeta_H=\X(\hbeta-\hbeta_H) \), whose squared length is \( \text{SS}_H \) (proof of @prp-glh-trinity), and \( \norm{\hat{\be}_H}^2=\text{SSE}_H \). So
\( n\norm{\M\hat{\be}_H}^2/\norm{\hat{\be}_H}^2=n\,\text{SS}_H/\text{SSE}_H=\text{Score} \). The uncentred \( R^2 \) of a
regression of \( \hat{\be}_H \) on \( \X \) is \( \norm{\M\hat{\be}_H}^2/\norm{\hat{\be}_H}^2 \). Computing it requires the restricted
residuals and one auxiliary regression, but not \( \hbeta \) or the unrestricted model's residual sum of
squares as such. This is why score tests are popular when the full model is hard to fit, as in tests
for omitted nonlinear terms or for heteroscedasticity (Chapter 21).
:::

::: {#exr-glh-known-sigma-lr}
[B3]

Suppose \( \sigma^2=\sigma_0^2 \) is known. Show that the likelihood ratio test of a reduced model rejects for
large \( (\text{SSE}_0-\text{SSE})/\sigma_0^2 \), the statistic of @exr-glh-known-sigma, and that
\( -2\log\Lambda \) is then exactly \( \chi^2(q) \) under \( H_0 \).
:::

### C. Going deeper

::: {#exr-glh-no-ump}
[C1]

In the canonical form with \( \sigma^2 \) known and \( q\ge2 \), consider testing \( \boldsymbol{\eta}_1=\bzero \) against the
single alternative \( \boldsymbol{\eta}_1=\mathbf{a}\ne\bzero \). Show by the Neyman–Pearson lemma that the most powerful
level-\( \alpha \) test rejects when \( \mathbf{a}\T\bz_1>\sigma\norm{\mathbf{a}}z_\alpha \), where \( z_\alpha \) is the upper \( \alpha \) point of
\( \Normal(0,1) \). Show that this test has power below \( \alpha \) against \( -\mathbf{a} \), and conclude that no
uniformly most powerful test of \( \boldsymbol{\eta}_1=\bzero \) exists.
:::

::: {.solution}
The likelihood ratio of \( \Normal(\mathbf{a},\sigma^2\I) \) to \( \Normal(\bzero,\sigma^2\I) \) at \( \bz_1 \) is
\( \exp\{(\mathbf{a}\T\bz_1-\norm{\mathbf{a}}^2/2)/\sigma^2\} \), increasing in \( \mathbf{a}\T\bz_1 \). Under the null,
\( \mathbf{a}\T\bz_1\sim\Normal(0,\sigma^2\norm{\mathbf{a}}^2) \), which gives the critical value, and the Neyman–Pearson
argument of Step 3 of the proof of @thm-glh-ump-invariant shows that this test is most powerful. Against
\( -\mathbf{a} \), \( \mathbf{a}\T\bz_1\sim\Normal(-\norm{\mathbf{a}}^2,\sigma^2\norm{\mathbf{a}}^2) \), so the power is
\( \Pr\{Z>z_\alpha+\norm{\mathbf{a}}/\sigma\}<\alpha \). A uniformly most powerful test would have to be most powerful
against both \( \mathbf{a} \) and \( -\mathbf{a} \). The most powerful tests against these two alternatives are
essentially unique and different, so no single test can be both.
:::

::: {#exr-glh-mlr-general}
[C2]

Extend @lem-glh-mlr: show that for \( 0\le\gamma_1<\gamma_2 \) the ratio \( f_{\gamma_2}(x)/f_{\gamma_1}(x) \) is increasing in
\( x \). (Write both densities as power series in \( t(x) \) with coefficients \( a_k(\gamma) \), and show that
\( a_k(\gamma_2)/a_k(\gamma_1) \) increases with \( k \).) Deduce that the \( F \) test is also uniformly most powerful
invariant for testing \( \gamma\le\gamma_1 \) against \( \gamma>\gamma_1 \), with a noncentral critical value. Such tests
appear in equivalence testing.
:::
