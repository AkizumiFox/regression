# Potential outcomes and the regression coefficient

Whether a treatment works needs much less than a model of the whole mechanism. Each unit has two possible responses,
with and without the treatment, and the effect on that unit is their difference. This language of potential outcomes,
due to Neyman (1923) for field trials and extended to observational studies by Rubin (1974), needs no model for the
responses, and under randomization it gives exact results. The regression coefficient of a treatment indicator is the natural estimator; its
usual standard error is wrong in a specific way.

## Potential outcomes

::: {#def-cau-potential-outcomes}
[Potential outcomes and average effects]

Consider \( n \) units and a binary treatment. Unit \( i \) has two **potential outcomes**: \( y_i(1) \), the response
it would show if treated, and \( y_i(0) \), the response it would show if not. With \( T_i\in\{0,1\} \) the treatment
it actually receives, the observed response is
\[
Y_i=T_iy_i(1)+(1-T_i)y_i(0).
\]{#eq-cau-consistency}

The **unit effect** is \( \tau_i=y_i(1)-y_i(0) \), and the **average treatment effect** of the \( n \) units is
\[
\tau=\frac1n\sum_{i=1}^n\tau_i=\bar{y}(1)-\bar{y}(0).
\]
When the units are instead drawn at random from a population, so that \( (Y(0),Y(1),T) \) is a random vector, the
population **average treatment effect** is \( \E\bigl(Y(1)-Y(0)\bigr) \) and the **average effect on the
treated** is \( \E\bigl(Y(1)-Y(0)\mid T=1\bigr) \).
:::

Equation @eq-cau-consistency assumes that unit \( i \)'s response depends only on its own treatment (**no
interference**, which fails for vaccines) and that there is one version of each treatment: Rubin's *stable unit
treatment value assumption*. In a structural causal model, \( Y(t) \) is the value of \( Y \), as a function of the
disturbances, in the model modified by \( \operatorname{do}(T=t) \), so \( \E\,Y(t)=\E\bigl(Y\mid\operatorname{do}(T=t)\bigr) \).
On the event \( T=t \), solving the original and the modified assignments in order gives every variable the same value, so
\( Y=Y(t) \) there, which is @eq-cau-consistency.

For each unit one potential outcome is observed and the other is **counterfactual**, so no unit effect is ever
observed (Holland's (1986) fundamental problem of causal inference). Averages can still be estimated from *different*
units, if treated and untreated units are comparable.

::: {#prp-cau-selection-bias}
[Selection bias]

Let \( (Y(0),Y(1),T) \) be random with finite means and \( 0<\Pr(T=1)<1 \). Then
\[
\E(Y\mid T=1)-\E(Y\mid T=0)
=\E\bigl(Y(1)-Y(0)\mid T=1\bigr)+\Bigl[\E\bigl(Y(0)\mid T=1\bigr)-\E\bigl(Y(0)\mid T=0\bigr)\Bigr].
\]{#eq-cau-selection}

:::

::: {.proof}
By @eq-cau-consistency, \( \E(Y\mid T=1)=\E(Y(1)\mid T=1) \) and \( \E(Y\mid T=0)=\E(Y(0)\mid T=0) \). Add and subtract
\( \E(Y(0)\mid T=1) \).
:::

The observed difference is the effect on the treated plus a **selection bias**, the difference between treated and
untreated units in the response they would have shown untreated. In @exm-cau-visits heavy users of the clinic would have
been less healthy than light users even without their visits. The data on \( (T,Y) \) cannot reveal the bias, because
\( Y(0) \) is never observed for the treated.

## Randomization

A **completely randomized experiment** treats \( n_1 \) of the \( n \) units, every subset of size \( n_1 \) being equally
likely, and leaves the other \( n_0=n-n_1 \) as controls. In Neyman's framework the potential outcomes are fixed numbers
and the only randomness is the assignment: no linearity, normality or common variance is assumed, and the units need
not be a sample from anything.

Write \( \bar{Y}_1 \) and \( \bar{Y}_0 \) for the mean observed responses of the treated and control units, \( s_1^2 \)
and \( s_0^2 \) for their sample variances, and
\[
S_t^2=\frac1{n-1}\sum_{i=1}^n\bigl(y_i(t)-\bar{y}(t)\bigr)^2,\qquad
S_\tau^2=\frac1{n-1}\sum_{i=1}^n(\tau_i-\tau)^2
\]
for the variances of the potential outcomes and of the unit effects over the \( n \) units. The treated units are a
simple random sample from the \( n \) units, and the moments of such samples are as follows (compare @lem-dsn-srs).

::: {#lem-cau-srs}
[Means of a simple random sample]

Let \( u_1,\dots,u_N \) and \( v_1,\dots,v_N \) be fixed numbers with means \( \bar{u} \), \( \bar{v} \) and
\( S_{uv}=\sum_a(u_a-\bar{u})(v_a-\bar{v})/(N-1) \). If \( \mathcal S \) is a simple random sample of \( m \) of the indices,
\( 1\le m\le N \), with sample means \( \bar{u}_{\mathcal S} \) and \( \bar{v}_{\mathcal S} \), then
\[
\E\,\bar{u}_{\mathcal S}=\bar{u},\qquad
\Cov(\bar{u}_{\mathcal S},\bar{v}_{\mathcal S})=\Bigl(\frac1m-\frac1N\Bigr)S_{uv}.
\]
:::

::: {.proof}
Let \( W_a \) indicate that \( a\in\mathcal S \). Then \( \E W_a=m/N \), and for \( a\ne c \),
\( \E W_aW_c=m(m-1)/\{N(N-1)\} \). Hence \( \Var W_a=q \) and
\( \Cov(W_a,W_c)=m(m-1)/\{N(N-1)\}-m^2/N^2=-q/(N-1) \), where \( q=(m/N)(1-m/N) \). Put \( \tilde{u}_a=u_a-\bar{u} \) and
\( \tilde{v}_a=v_a-\bar{v} \). Since \( \sum_aW_a=m \), \( \bar{u}_{\mathcal S}-\bar{u}=m^{-1}\sum_aW_a\tilde{u}_a \), whose mean
is \( N^{-1}\sum_a\tilde{u}_a=0 \). For the covariance,
\[
m^2\Cov(\bar{u}_{\mathcal S},\bar{v}_{\mathcal S})=\sum_{a,c}\Cov(W_a,W_c)\,\tilde{u}_a\tilde{v}_c
=q\Bigl(\sum_a\tilde{u}_a\tilde{v}_a-\frac1{N-1}\sum_{a\ne c}\tilde{u}_a\tilde{v}_c\Bigr).
\]
Because \( \sum_a\tilde{u}_a=0 \), \( \sum_{a\ne c}\tilde{u}_a\tilde{v}_c=-\sum_a\tilde{u}_a\tilde{v}_a \), and the bracket is
\( \sum_a\tilde{u}_a\tilde{v}_a\,N/(N-1)=NS_{uv} \). So the covariance is \( qNS_{uv}/m^2=(1/m-1/N)S_{uv} \).
:::

::: {#thm-cau-randomized}
[Randomization and the difference in means]

In a completely randomized experiment with \( n_1\ge1 \) and \( n_0\ge1 \), let \( \hat{\tau}=\bar{Y}_1-\bar{Y}_0 \).

::: {.enumerate options="label=(\alph*)"}
1. \( \hat{\tau} \) is unbiased for the average treatment effect: \( \E\hat{\tau}=\tau \).

2. Its variance over the randomization is
   \[
   \Var(\hat{\tau})=\frac{S_1^2}{n_1}+\frac{S_0^2}{n_0}-\frac{S_\tau^2}{n}.
   \]{#eq-cau-neyman-variance}

3. If \( n_1,n_0\ge2 \), the estimator \( \hat{V}=s_1^2/n_1+s_0^2/n_0 \) has
   \( \E\hat{V}=\Var(\hat{\tau})+S_\tau^2/n\ge\Var(\hat{\tau}) \), with equality iff the unit effects are all equal.

4. \( \hat{\tau} \) is the least squares coefficient of \( T \) in the regression of the observed responses on an
   intercept and the treatment indicator.

5. If instead the units are independent copies of \( (Y(0),Y(1),T) \) with \( T \) independent of \( (Y(0),Y(1)) \),
   \( 0<\Pr(T=1)<1 \) and finite means, then \( \E(Y\mid T=1)-\E(Y\mid T=0)=\E\bigl(Y(1)-Y(0)\bigr) \), and the
   average effect on the treated equals the average treatment effect.
:::

:::

::: {.proof}
Let \( \mathcal S \) be the set of treated units, a simple random sample of size \( n_1 \). Then \( \bar{Y}_1 \) is the
sample mean of the \( y_i(1) \) over \( \mathcal S \), and \( \bar{Y}_0 \) is the mean of the \( y_i(0) \) over the complement.

(a) The complement of \( \mathcal S \) is a simple random sample of size \( n_0 \). By @lem-cau-srs,
\( \E\bar{Y}_1=\bar{y}(1) \) and \( \E\bar{Y}_0=\bar{y}(0) \).

(b) Let \( \bar{y}_{\mathcal S}(0) \) be the mean of the \( y_i(0) \) over \( \mathcal S \). The control mean is
\( \bar{Y}_0=\bigl(n\bar{y}(0)-n_1\bar{y}_{\mathcal S}(0)\bigr)/n_0 \), so
\[
\hat{\tau}=\bar{u}_{\mathcal S}-\frac{n}{n_0}\bar{y}(0),\qquad u_i=y_i(1)+\frac{n_1}{n_0}y_i(0).
\]{#eq-cau-u-representation}

By @lem-cau-srs, \( \Var(\hat{\tau})=(1/n_1-1/n)S_u^2=\{n_0/(nn_1)\}S_u^2 \). Expanding,
\( S_u^2=S_1^2+(n_1/n_0)^2S_0^2+2(n_1/n_0)S_{10} \), where \( S_{10} \) is the covariance of the two potential outcomes over
the units, and \( S_\tau^2=S_1^2+S_0^2-2S_{10} \). Eliminating \( S_{10} \),
\[
S_u^2=\frac{n}{n_0}S_1^2+\frac{nn_1}{n_0^2}S_0^2-\frac{n_1}{n_0}S_\tau^2 ,
\]
and multiplying by \( n_0/(nn_1) \) gives @eq-cau-neyman-variance.

(c) The treated responses are a simple random sample of \( n_1\ge2 \) of the numbers \( y_i(1) \), so
\( \E s_1^2=S_1^2 \) by @lem-dsn-srs; likewise \( \E s_0^2=S_0^2 \). Hence \( \E\hat{V}=S_1^2/n_1+S_0^2/n_0 \), which
exceeds @eq-cau-neyman-variance by \( S_\tau^2/n \). This is zero iff \( \tau_i=\tau \) for every \( i \).

(d) This is @prp-lm-indicator-means with two groups and control as the reference.

(e) By @eq-cau-consistency and independence, \( \E(Y\mid T=1)=\E(Y(1)\mid T=1)=\E\,Y(1) \), and likewise for
\( T=0 \). The same independence gives \( \E(Y(1)-Y(0)\mid T=1)=\E(Y(1)-Y(0)) \).
:::

Part (e) shows randomization removing the selection bias of @eq-cau-selection; part (a) needs no sampling model at all. The last term of (b) cannot be estimated, because it depends
on how \( y_i(1) \) is paired with \( y_i(0) \) within units, which no unit reveals. Part (c) accepts a conservative answer,
exact when the treatment adds the same amount to every unit.

::: {#exm-cau-twenty}
[All randomizations of twenty units]

A synthetic population of \( n=20 \) units has \( \tau=2.100 \), \( S_1^2=12.129 \),
\( S_0^2=5.678 \) and \( S_\tau^2=1.698 \); the effect grows with \( y_i(0) \). With \( n_1=8 \)
there are \( 125970 \) assignments, few enough to enumerate. Over all of them \( \hat{\tau} \) has mean
\( 2.100 \) and variance \( 1.9043 \), which is @eq-cau-neyman-variance to every digit. Neyman's
\( \hat{V} \) has mean \( 1.9892 \), too large by \( S_\tau^2/n=0.0849 \). The pooled least squares
variance estimate averages \( 1.7055 \), below the true variance, because it gives the more variable treated
arm the weight of the smaller arm (@exr-cau-pooled-se). [Figure 25.2.1](#fig-cau-randomization) shows the
distribution of \( \hat{\tau} \).
:::

```{.python .run #cell-potential-outcomes-population}
import itertools
import numpy as np

rng = np.random.default_rng(2502)
n, n1 = 20, 8
n0 = n - n1
y0 = np.round(rng.normal(10, 3, size=n), 1)                      # potential outcome, control
y1 = np.round(y0 + 2 + 0.5 * (y0 - 10) + rng.normal(0, 1, n), 1)  # potential outcome, treated
tau = np.mean(y1 - y0)                                            # the sample average effect
S1, S0, St = np.var(y1, ddof=1), np.var(y0, ddof=1), np.var(y1 - y0, ddof=1)
print(f"tau = {tau:.3f};  S1^2 = {S1:.3f}, S0^2 = {S0:.3f}, S_tau^2 = {St:.3f}")
```

```{.python .run #cell-potential-outcomes-enumerate}
idx = np.array(list(itertools.combinations(range(n), n1)))       # every possible treated set
treated = np.zeros((len(idx), n), dtype=bool)
treated[np.arange(len(idx))[:, None], idx] = True                 # one row per assignment
m1, m0 = (treated @ y1) / n1, (~treated @ y0) / n0                # arm means, one row per assignment
diff = m1 - m0                                                    # difference in means
s1 = (treated @ y1 ** 2 - n1 * m1 ** 2) / (n1 - 1)                # within-arm variances
s0 = (~treated @ y0 ** 2 - n0 * m0 ** 2) / (n0 - 1)
V_hat = s1 / n1 + s0 / n0                                         # Neyman's variance estimate
V_neyman = S1 / n1 + S0 / n0 - St / n
print(f"{len(diff)} assignments")
print(f"mean of difference in means {diff.mean():.3f}   (tau = {tau:.3f})")
print(f"variance over assignments   {diff.var():.4f}  (Neyman's formula {V_neyman:.4f})")
print(f"mean of variance estimate   {V_hat.mean():.4f}  (excess S_tau^2/n = {St / n:.4f})")
```

::: {when-format="html"}
![**Figure 25.2.1.** The randomization distribution of the difference in means in @exm-cau-twenty, over all
\( 125970 \) assignments of \( 8 \) of the \( 20 \) units to treatment. The dashed line is \( \tau \); the curve is
the normal density with Neyman's variance @eq-cau-neyman-variance.](randomization_distribution.svg){#fig-cau-randomization width=70%}
:::

::: {when-format="pdf"}
![The randomization distribution of the difference in means in @exm-cau-twenty, over all
\( 125970 \) assignments of \( 8 \) of the \( 20 \) units to treatment. The dashed line is \( \tau \); the curve is
the normal density with Neyman's variance @eq-cau-neyman-variance.](randomization_distribution.pdf){width=70%}
:::

Finite-population central limit theorems (Hájek 1960; Li and Ding 2017), whose proofs are beyond this chapter, make
\( (\hat{\tau}-\tau)/\sqrt{\Var\hat{\tau}} \) asymptotically standard normal under conditions on the potential outcomes, so
\( \hat{\tau}\pm1.96\sqrt{\hat{V}} \) has coverage at least \( 95\% \) in the limit. In the example, with twenty units, it covers \( \tau \)
for a fraction \( 0.931 \) of the assignments.

## What the regression coefficient estimates

By part (d), the least squares coefficient of \( T \) is unbiased for the average effect because of the design, not
because the model \( \E(Y_i)=\beta_0+\beta_1T_i \) is true: its "errors" are not random draws and have different variances
in the two arms. What the model gets wrong is the standard error. The pooled estimate of
\( \sigma^2(1/n_1+1/n_0) \) (@thm-lm-moments) can be too small or too large when \( n_1\ne n_0 \) (@exr-cau-pooled-se),
while Neyman's \( \hat{V} \) is exactly the HC2 sandwich estimate of
[Chapter 21](../ch21-nonnormality-heteroscedasticity-serial/index.html) (@exr-cau-hc2, @thm-het-sandwich). Covariates
measured before randomization can be added to this regression to gain precision;
[Section 25.6](06-randomized-adjustment.html) shows how, and finds a surprise.

::: {.remark}
[Fisher's sharp null]

Fisher (1935) tested the *sharp* null \( y_i(1)=y_i(0) \) for every unit, under which the randomization distribution of
any statistic is known exactly. [Chapter 23](../ch23-resampling-inference/index.html) develops these permutation tests (@thm-bs-permutation).
:::

## Exercises

### A. Check your understanding

::: {#exr-cau-four-units}
[A1]

Four units have potential outcomes \( \bigl(y_i(0),y_i(1)\bigr)=(1,4) \), \( (3,3) \), \( (5,9) \), \( (7,8) \). Two are
treated by complete randomization. List the six possible values of \( \hat{\tau} \), and verify
@thm-cau-randomized(a) and (b) directly.
:::

::: {.solution}
The unit effects are \( 3,0,4,1 \), so \( \tau=2 \). For the treated pairs \( \{1,2\},\{1,3\},\{1,4\},\{2,3\},\{2,4\},\{3,4\} \)
the treated means are \( 3.5,6.5,6,6,5.5,8.5 \) and the control means \( 6,5,4,4,3,2 \), giving
\( \hat{\tau}=-2.5,1.5,2,2,2.5,6.5 \). Their mean is \( 12/6=2=\tau \). Their squared deviations from \( 2 \) sum to
\( 20.25+0.25+0+0+0.25+20.25=41 \), so \( \Var\hat{\tau}=41/6 \). The formula: \( S_1^2=26/3 \), \( S_0^2=20/3 \),
\( S_\tau^2=10/3 \), and \( (26/3)/2+(20/3)/2-(10/3)/4=46/6-5/6=41/6 \).
:::

### B. Practice

::: {#exr-cau-pooled-se}
[B1]

Let \( s_p^2=\{(n_1-1)s_1^2+(n_0-1)s_0^2\}/(n-2) \) and \( \hat{V}_{\mathrm{OLS}}=s_p^2(1/n_1+1/n_0) \), the least squares
variance estimate for the coefficient of \( T \). Show that \( \hat{V}_{\mathrm{OLS}}=\hat{V} \) when \( n_1=n_0 \), and that in
general
\[
\E\hat{V}_{\mathrm{OLS}}-\E\hat{V}=\frac{(n-1)(n_1-n_0)(S_1^2-S_0^2)}{n_1n_0(n-2)} .
\]
Conclude that when the smaller arm has the more variable potential outcomes, the least squares variance estimate is
smaller on average than Neyman's, and it can be smaller than the true variance.
:::

::: {.solution}
If \( n_1=n_0=m \), then \( s_p^2=(s_1^2+s_0^2)/2 \) and \( \hat{V}_{\mathrm{OLS}}=s_p^2\cdot2/m=\hat{V} \). In general
\( \E s_t^2=S_t^2 \) (@lem-dsn-srs), so \( \E\hat{V}_{\mathrm{OLS}}=n\{(n_1-1)S_1^2+(n_0-1)S_0^2\}/\{n_1n_0(n-2)\} \). The
coefficient of \( S_1^2 \) in \( \E\hat{V}_{\mathrm{OLS}}-\E\hat{V} \) is
\[
\frac{n(n_1-1)-n_0(n-2)}{n_1n_0(n-2)}=\frac{n(n_1-n_0)-(n_1-n_0)}{n_1n_0(n-2)}=\frac{(n-1)(n_1-n_0)}{n_1n_0(n-2)},
\]
using \( n=n_1+n_0 \); by symmetry the coefficient of \( S_0^2 \) is its negative. If \( n_1<n_0 \) and \( S_1^2>S_0^2 \), the
difference is negative. Since \( \E\hat{V} \) exceeds \( \Var\hat{\tau} \) only by \( S_\tau^2/n \), the least squares estimate
underestimates the variance on average whenever the negative difference is larger than that. In @exm-cau-twenty,
\( \E\hat{V}_{\mathrm{OLS}}=1.7055 \) while \( \Var\hat{\tau}=1.9043 \).
:::

::: {#exr-cau-hc2}
[B2]

In the regression of the observed \( Y_i \) on \( \bone \) and the treatment indicator, the leverages are
\( h_{ii}=1/n_1 \) for treated units and \( 1/n_0 \) for controls (@prp-proj-leverage). The HC2 estimate of
\( \Cov(\hbeta) \) is \( (\X\T\X)^{-1}\X\T\diag\bigl(\hat{e}_i^2/(1-h_{ii})\bigr)\X(\X\T\X)^{-1} \). Show that its entry for the
coefficient of \( T \) is exactly Neyman's \( \hat{V} \).
:::

::: {.solution}
Reparameterize with the two arm indicators as columns, which changes neither the fit nor the estimate of the difference
of the arm means (@prp-lm-reparameterization). Then \( \X\T\X=\diag(n_1,n_0) \), the residuals are \( Y_i-\bar{Y}_1 \) and
\( Y_i-\bar{Y}_0 \), and the sandwich is diagonal with entries
\( n_1^{-2}\sum_{\text{treated}}\hat{e}_i^2/(1-1/n_1)=s_1^2/n_1 \) and similarly \( s_0^2/n_0 \). The variance of the difference
of the two coefficients is their sum, \( \hat{V} \).
:::

### C. Going deeper

::: {#exr-cau-sharper-bound}
[C1]

Show that \( S_\tau^2\ge(S_1-S_0)^2 \), and hence that
\[
\Var(\hat{\tau})\le\frac{S_1^2}{n_1}+\frac{S_0^2}{n_0}-\frac{(S_1-S_0)^2}{n},
\]
a bound that is attained when \( y_i(1) \) is an increasing affine function of \( y_i(0) \). Propose a conservative variance
estimator that is never larger than \( \hat{V} \).
:::
