# Complete cases and quick repairs

Confronted with holes, software does one of two things without being asked: it
deletes the rows that have any, or it fills them and proceeds as though nothing had
happened. Both are defensible in circumstances that can be stated precisely, and
disastrous outside them. This section states them.

## When deleting rows is safe

Write \( C_i=1 \) when unit \( i \)'s record is complete in every variable the
analysis uses. A **complete-case analysis** applies the analysis to the subsample
\( \{i:C_i=1\} \), and its behaviour is governed by one thing: the conditional
distribution of \( \mathbf{d}_{(i)} \) given \( C_i=1 \).

::: {#prp-mis-complete-case}
[Validity of a complete-case analysis]

Let the rows \( (\mathbf{d}_{(i)},C_i) \) be independent and identically
distributed, with \( \Pr(C_i=1)>0 \).

::: {.enumerate options="label=(\alph*)"}
1. *(General condition.)* Let a procedure estimate a functional \( T \)
   consistently at every law in a class containing both \( f \), the law of
   \( \mathbf{d}_{(i)} \), and \( f_{C=1} \), its conditional law given
   \( C_i=1 \). Applied to the complete cases the procedure stays consistent for
   \( T(f) \) if and only if \( T(f_{C=1})=T(f) \). Under MCAR
   \( f_{C=1}=f \), and every such functional is preserved.

2. *(Regression with selection on the covariates.)* Partition
   \( \mathbf{d}_{(i)}=(y_i,\x_{(i)}) \) and suppose
   \[
   \Pr(C_i=1\mid y_i,\x_{(i)})=\Pr(C_i=1\mid \x_{(i)})\quad\text{almost surely}.
   \]{#eq-mis-selection-on-x}

   Then the conditional law of \( y_i \) given \( \x_{(i)} \) is the same among the
   complete cases as in the population. So if a model specifies that conditional law
   correctly, the complete-case least squares estimate in a linear model is unbiased
   for \( \bbeta \) conditionally on the selected covariates; and for a generalized
   linear model (@def-glm-model) the complete-case maximum likelihood estimate is
   consistent and asymptotically normal, with the usual model-based covariance of
   @thm-glm-asymptotics computed from the complete cases alone, provided their number
   grows and their covariate second moments converge to a nonsingular limit.

3. *(Otherwise.)* If @eq-mis-selection-on-x fails, write
   \( \pi(y,\x)=\Pr(C=1\mid y,\x) \). The complete-case regression function is
   \[
   \E(y\mid\x,C=1)=\E(y\mid\x)
   +\frac{\Cov\{y,\pi(y,\x)\mid\x\}}{\E\{\pi(y,\x)\mid\x\}} ,
   \]{#eq-mis-cc-bias}

   so the complete cases follow the intended regression function exactly when
   \( \pi \) is conditionally uncorrelated with the response at every \( \x \).
:::

:::

::: {.proof}
(a) The complete cases are an independent sample from \( f_{C=1} \), and their
number grows without bound since \( \Pr(C=1)>0 \). By the assumed consistency of the
procedure at \( f_{C=1} \), it converges to \( T(f_{C=1}) \), which is the target
\( T(f) \) if and only if the two agree. Under MCAR, \( C \) is independent of
\( \mathbf{d} \), hence \( f_{C=1}=f \).

(b) By Bayes's rule,
\[
f(y\mid\x,C=1)=\frac{\Pr(C=1\mid y,\x)\,f(y\mid\x)}{\Pr(C=1\mid\x)} ,
\]{#eq-mis-bayes-selection}

and under @eq-mis-selection-on-x the two selection probabilities cancel, leaving
\( f(y\mid\x) \). The complete cases are therefore a sample from a joint law whose
conditional of \( y \) given \( \x \) is the modelled one, with only the marginal
law of \( \x \) distorted. Every likelihood built from \( f(y\mid\x;\bbeta) \)
alone conditions on the covariates and is unaffected by that distortion.

(c) Multiply @eq-mis-bayes-selection by \( y \) and integrate:
\[
\begin{aligned}
\E(y\mid\x,C=1)&=\frac{\E\{y\,\pi(y,\x)\mid\x\}}{\E\{\pi(y,\x)\mid\x\}}\\
&=\E(y\mid\x)+\frac{\Cov\{y,\pi(y,\x)\mid\x\}}{\E\{\pi(y,\x)\mid\x\}} ,
\end{aligned}
\]
using \( \E(y\pi)=\E(y)\E(\pi)+\Cov(y,\pi) \) throughout at fixed \( \x \).
:::

Part (b) is the most useful positive result in the subject and is routinely
overlooked. *Missingness may depend on the covariates as strongly as it likes, in any
functional form, without harming a correctly specified regression.* The complete
cases are a biased sample of units — their covariate distribution is wrong, so any
marginal summary of them is wrong — and yet the fitted conditional mean is right.
This is selection on a regressor, as in
[Chapter 25](../ch25-causal-interpretation/index.html); the binary-response version
is @exr-bin-selection-on-x.

::: {.warning}
Part (b) needs the conditional *model* to be right. If the true regression is curved
and the fitted one is a straight line, selecting on \( \x \) changes which part of the
covariate range dominates the fit, and the complete-case slope moves with it.
Selection on the covariates is harmless for a correct model and can be very harmful
for a wrong one.
:::

::: {#exm-mis-four-mechanisms}
[Four mechanisms, one regression]

Take \( y_i=1+0.8x_i+0.9w_i+\varepsilon_i \) with \( x_i,w_i,\varepsilon_i \)
independent standard normal and \( n=300 \), and lose the response under four
mechanisms: completely at random; at random on \( x \); at random on \( x \) and
\( w \); and not at random, on \( y \).

[**Table 41.2.1.** Mean complete-case estimate of \( \beta_1=0.8 \) over
\( 2000 \) replicates, for two analyses under four
mechanisms.]{#tab-mis-four}

| mechanism | complete cases | fit \( y\sim x \) | fit \( y\sim x+w \) |
|---|---|---|---|
| MCAR | \( 194 \) | \( 0.801 \) | \( 0.801 \) |
| MAR on \( x \) | \( 181 \) | \( 0.804 \) | \( 0.802 \) |
| MAR on \( x \) and \( w \) | \( 178 \) | \( 0.637 \) | \( 0.800 \) |
| MNAR on \( y \) | \( 178 \) | \( 0.593 \) | \( 0.669 \) |

The first two rows are @prp-mis-complete-case(b): whatever the missingness does with
\( x \), the regression of \( y \) on \( x \) is untouched. The third is the warning
attached to it. The mechanism is MAR, using only recorded variables, yet the short
regression omitting \( w \) is badly biased, because selecting on \( w \) selects on
the part of \( y \) that the short model calls error; adding \( w \) repairs it
completely. The fourth row cannot be repaired by any choice of regressors, because
the selection is on \( y \) itself.
:::

::: {when-format="html"}
![**Figure 41.2.1.** Complete-case estimates of \( \beta_1=0.8 \) over
\( 2000 \) replicates under four mechanisms, for the regression of
\( y \) on \( x \) alone and on \( x \) and \( w \). The horizontal line is the
true value.](complete_case_bias.svg){#fig-mis-cc width=100%}
:::

::: {when-format="pdf"}
![Complete-case estimates of \( \beta_1=0.8 \) over
\( 2000 \) replicates under four mechanisms, for the regression of
\( y \) on \( x \) alone and on \( x \) and \( w \). The horizontal line is the
true value.](complete_case_bias.pdf){width=100%}
:::

```{.python .run #cell-complete-case-scenarios}
import numpy as np

rng = np.random.default_rng(4101)
n, reps = 300, 2000
beta = np.array([1.0, 0.8, 0.9])                      # intercept, coefficient of x, of w


def one_run(mechanism, rng):
    x, w = rng.normal(size=n), rng.normal(size=n)
    y = beta[0] + beta[1] * x + beta[2] * w + rng.normal(size=n)
    if mechanism == "MCAR":
        eta = np.full(n, -0.6)
    elif mechanism == "MAR on x":
        eta = -0.6 + 1.6 * x                          # depends on a regressor only
    elif mechanism == "MAR on x and w":
        eta = -0.6 + 1.2 * x + 1.5 * w                # depends on observed variables only
    else:                                             # MNAR
        eta = -0.6 + 1.2 * (y - beta[0])              # depends on the value itself
    missing = rng.uniform(size=n) < 1 / (1 + np.exp(-eta))
    obs = ~missing
    short = np.column_stack([np.ones(obs.sum()), x[obs]])            # y on x
    long = np.column_stack([np.ones(obs.sum()), x[obs], w[obs]])     # y on x and w
    b_short, *_ = np.linalg.lstsq(short, y[obs], rcond=None)
    b_long, *_ = np.linalg.lstsq(long, y[obs], rcond=None)
    return b_short[1], b_long[1], obs.sum()


names = ["MCAR", "MAR on x", "MAR on x and w", "MNAR"]
out = {m: np.array([one_run(m, rng) for _ in range(reps)]) for m in names}
for m in names:
    s, l, k = out[m].T
    print(f"{m:16s} complete cases {k.mean():5.0f}   y ~ x: {s.mean():.4f}"
          f"   y ~ x + w: {l.mean():.4f}")
```

::: {.idea}
Ask what the missingness depends on, then put those variables into the analysis. If
everything the mechanism uses is in the model and the model is right, the complete
cases are enough. If the mechanism uses the response, or a variable the model leaves
out, they are not.
:::

## Available cases

A tempting middle course computes each summary from whichever records carry the
variables it needs: a mean from everyone who reported that variable, a covariance
from everyone who reported both. This **available-case** analysis, or pairwise
deletion, uses more data, at the price of a covariance matrix whose entries come from
different subsamples and need not be jointly realizable.

::: {#exm-mis-pairwise}
[A correlation matrix with a negative eigenvalue]

Twelve records on three variables, in three patterns of four: \( d_3 \) absent in
the first, \( d_2 \) in the second, \( d_1 \) in the third, with numbers chosen so
that the available pairs give \( r_{12}=0.9 \), \( r_{13}=0.9 \) and
\( r_{23}=-0.9 \). Each is an ordinary correlation on its own four points; together
they are impossible, since by @prp-rv-three-correlations the third must lie between
\( 0.81-\sqrt{0.19\times0.19}=0.62 \) and \( 0.81+0.19=1 \). The matrix has
smallest eigenvalue \( -0.800 \), so any procedure needing
it to be nonnegative definite — a Cholesky factor, a GLS weight matrix, a normal
likelihood — fails, even though the matrix itself is perfectly invertible.
:::

```{.python .run #cell-complete-case-pairwise}
nan = np.nan
u = np.array([-1.5, -0.5, 0.5, 1.5])
d = np.array([0.5, -1.5, 1.5, -0.5])                  # orthogonal to u, same length
pair = np.column_stack([u, 0.9 * u + np.sqrt(1 - 0.81) * d])             # correlation +0.9
flip = np.column_stack([u, -(0.9 * u + np.sqrt(1 - 0.81) * d)])          # correlation -0.9
D = np.vstack([
    np.column_stack([pair, np.full(4, nan)]),         # pattern A: z3 absent
    np.column_stack([pair[:, 0], np.full(4, nan), pair[:, 1]]),          # pattern B: z2 absent
    np.column_stack([np.full(4, nan), flip]),         # pattern C: z1 absent
])
R = np.ones((3, 3))
for i in range(3):
    for j in range(i + 1, 3):
        both = ~np.isnan(D[:, i]) & ~np.isnan(D[:, j])
        R[i, j] = R[j, i] = np.corrcoef(D[both, i], D[both, j])[0, 1]
print("available-case correlation matrix\n", R)
print("eigenvalues", np.round(np.linalg.eigvalsh(R), 3))
```

Under MCAR each entry is consistent and the difficulty is a finite-sample one; under
MAR the entries are not consistent at all. Available-case analysis has uses in
exploratory work and none in inference.

## Single imputation

The other reflex is to fill the holes and carry on. **Mean imputation** replaces every
missing entry of a variable by the mean of its recorded entries, leaving a pile of
entries exactly at the sample mean and a spread that is too small.

::: {#prp-mis-mean-imputation}
[Mean imputation deflates the variance]

Let \( z_1,\dots,z_n \) be independent and identically distributed with finite
variance \( \sigma^2 \), of which \( k \) are recorded and \( n-k \) missing, the
mechanism being completely at random; let \( \tilde z_i \) equal \( z_i \) when recorded and
\( \bar z_{\mathrm{obs}} \) otherwise. Then the sample variance of the filled
column satisfies
\[
\E\Bigl\{\frac1{n-1}\sum_{i=1}^n(\tilde z_i-\bar{\tilde z})^2\ \Big|\ k\Bigr\}
=\frac{k-1}{n-1}\,\sigma^2 ,
\]{#eq-mis-deflation}

so the estimate of \( \sigma^2 \) is too small by the factor
\( (k-1)/(n-1)\approx k/n \), the proportion recorded.
:::

::: {.proof}
The filled column's mean is \( \bar{\tilde z}=\bar z_{\mathrm{obs}} \), because the
\( n-k \) inserted values are all equal to it. The imputed entries therefore
contribute nothing to \( \sum_i(\tilde z_i-\bar{\tilde z})^2 \), which reduces to
\( \sum_{i\ \mathrm{obs}}(z_i-\bar z_{\mathrm{obs}})^2 \). Given \( k \), and under
MCAR, that sum has expectation \( (k-1)\sigma^2 \). Dividing by \( n-1 \)
gives @eq-mis-deflation.
:::

With \( n=300 \) and \( 30\% \) of values missing the factor is
\( 0.699 \), and the simulation gives
\( 0.700 \) as the average filled-in sample variance of a
standard normal column. Covariances with fully recorded variables are deflated by the
same mechanism, so correlations shrink towards zero and coefficients with them.

**Regression imputation** replaces a missing \( z_i \) by its predicted value from a
regression on the recorded variables. It fixes the mean but not the spread, since the
imputed points lie exactly on the fitted surface. In the same simulation, filling
\( 30\% \) of a variable whose correlation with the predictor is
\( 0.6 \) raises the sample correlation to
\( 0.667 \) and lowers the standard deviation to
\( 0.899 \). Adding a draw from the estimated residual
distribution — **stochastic regression imputation** — repairs both, and is the first
honest imputation on this list. What it does not repair is the *uncertainty*: the
analysis treats \( n \) values as measured when \( n-k \) were invented, and every
standard error is too small. [Section 41.5](05-multiple-imputation.html) measures the
damage — a nominal \( 95\% \) interval covering
\( 75.9\% \) of the time in the example there — and repairs it.

::: {.warning}
Never treat a single set of imputed values as data. The estimate may be fine; the
standard error, the \( t \) statistic and the interval are not, because the imputed
values are correlated with the imputation model's own estimates. The same warning
appeared in
[Section 18.6](../ch18-covariance-and-design/06-missing-observations.html) for the
filled-in plot of a designed experiment, where @thm-dsn-missing(e) gives the exact
correction available in that special case.
:::

## Weighting

Survey statisticians treat unit nonresponse differently. If unit \( i \) is
complete with probability \( \pi_i=\Pr(C_i=1\mid\mathbf{d}_{(i)}) \), then for any
estimating function \( \mathbf{g} \) with \( \E\{\mathbf{g}(\mathbf{d};\bbeta)\}=\bzero \) at
the true \( \bbeta \), the **inverse-probability weighted** equations
\[
\sum_{i=1}^n \frac{C_i}{\pi_i}\,\mathbf{g}(\mathbf{d}_{(i)};\bbeta)=\bzero
\]{#eq-mis-weighted}

are unbiased, because \( \E(C_i/\pi_i\mid\mathbf{d}_{(i)})=1 \). Under MAR,
\( \pi_i \) is a function of recorded variables and can be estimated, typically by a
logistic regression of \( C_i \) on them, and the estimator is then consistent when
that model is right. This is Horvitz and Thompson's (1952) idea transplanted from
sampling to nonresponse, and its weakness is visible in the formula: a unit with a
small estimated \( \pi_i \) receives a huge weight, so the estimator is least stable
where the data are thinnest. Augmenting @eq-mis-weighted with a term built from a
model for \( \E(\mathbf{g}\mid\text{recorded}) \) gives the **doubly robust**
estimators of Robins, Rotnitzky and Zhao (1994), consistent if either model is right.
Weighting is the natural language for unit nonresponse, where there is nothing to
impute; likelihood and imputation are the language for item nonresponse, where there
is.

## Exercises

### A. Check your understanding

::: {#exr-mis-cc-mean}
[A1]

Under the mechanism of @prp-mis-complete-case(b), show that the complete-case
sample mean of \( y \) is generally not consistent for \( \E(y) \), even though the
complete-case regression coefficients are consistent. Reconcile the two statements.
:::

::: {.solution}
Selection changes the law of \( \x \), and \( \E(y)=\E\{\E(y\mid\x)\} \) averages
the conditional mean over that law; the complete cases average it over the selected
law instead. The regression coefficients describe \( \E(y\mid\x) \), which selection
leaves alone, so both statements hold at once: the conditional model is right and
the marginal average of it is wrong.
:::

### B. Practice

::: {#exr-mis-cc-bias-normal}
[B1]

In @eq-mis-cc-bias take \( y\mid\x\sim\Normal(\x\T\bbeta,\sigma^2) \) and
\( \pi(y,\x)=\Phi(a+by) \). Show that
\[
\E(y\mid\x,C=1)=\x\T\bbeta+\frac{b\sigma^{2}}{\sqrt{1+b^{2}\sigma^{2}}}\,
\frac{\phi(c)}{\Phi(c)},\qquad
c=\frac{a+b\,\x\T\bbeta}{\sqrt{1+b^{2}\sigma^{2}}},
\]
and deduce that the complete-case regression function is nonlinear in \( \x \)
whenever \( b\ne0 \).
:::

::: {.solution}
With \( y=\x\T\bbeta+\sigma e \), \( e \) standard normal,
\( \E\{\Phi(a+by)\}=\Phi(c) \) by the standard convolution identity, and
\( \E\{e\,\Phi(a+by)\} = b\sigma\phi(c)/\sqrt{1+b^2\sigma^2} \), obtained by writing
\( \Phi(a+by)=\Pr(U\le a+by) \) for an independent standard normal \( U \) and
integrating by parts (equivalently, by Stein's identity
\( \E\{e\,h(e)\}=\E\{h'(e)\} \)). Dividing the second by the first and substituting
into @eq-mis-cc-bias gives the display. The added term depends on \( \x \) through
\( c \), and the ratio \( \phi/\Phi \) is not affine, so the complete-case mean is
not linear in \( \x \); in particular no single shift of the intercept absorbs it.
:::

::: {#exr-mis-deflation-check}
[B2]

Verify @eq-mis-deflation by simulation for \( n=300 \) and three missing fractions.
Then show that the sample covariance of the filled column with a fully recorded
variable \( w \) carries the *same* factor \( (k-1)/(n-1) \), and explain why the
fraction recorded \( k/n \) is only the limiting form of it.
:::

::: {.solution}
The imputed entries satisfy \( \tilde z_i-\bar{\tilde z}=0 \), so
\( \sum_i(\tilde z_i-\bar{\tilde z})(w_i-\bar w)
=\sum_{\mathrm{obs}}(z_i-\bar z_{\mathrm{obs}})w_i \). Given \( k \) and under MCAR
its expectation is \( k\sigma_{zw}-\sigma_{zw}=(k-1)\sigma_{zw} \), the second term
coming from \( \E(\bar z_{\mathrm{obs}}\sum_{\mathrm{obs}}w_i)=\sigma_{zw} \).
Dividing by \( n-1 \) gives \( \{(k-1)/(n-1)\}\sigma_{zw} \), the factor
of @eq-mis-deflation. The two candidates differ only through the divisors used for
the sample moments, and \( (k-1)/(n-1)\to k/n \) as both grow: \( k/n \) is the
population-level statement, \( (k-1)/(n-1) \) the exact finite-sample one, and there
is no second phenomenon hiding in the gap.
:::

::: {#exr-mis-weighted-unbiased}
[B3]

Prove that @eq-mis-weighted is an unbiased estimating equation when \( \pi_i>0 \) for
every unit, and show by example that the conclusion fails if \( \pi_i=0 \) for a
set of covariate values with positive probability. What does that failure
correspond to in @prp-mis-complete-case?
:::

### C. Going deeper

::: {#exr-mis-cc-nonlinear}
[C1]

Let \( \E(y\mid x)=\beta_0+\beta_1x+\beta_2x^2 \) with \( \beta_2\ne0 \), and let
the analyst fit a straight line. Show that even under @eq-mis-selection-on-x the
complete-case slope generally differs from the full-data slope, and identify the
functional of the covariate distribution that the straight-line slope estimates
(see @thm-rv-blp and @thm-proj-blp). Under what selection would the two agree?
:::

::: {.solution}
The straight-line least squares coefficient converges to the best linear predictor
coefficient \( \Cov(y,x)/\Var(x) \) computed under whatever law \( x \) has. With
\( \E(y\mid x)=\beta_0+\beta_1x+\beta_2x^2 \) this is
\( \beta_1+\beta_2\Cov(x^2,x)/\Var(x) \), which depends on the third and second
central moments of \( x \). Selection changes those moments, hence the limit, unless
it leaves \( \Cov(x^2,x)/\Var(x) \) unchanged — for instance under MCAR, or under a
selection probability that is an even function of \( x \) about its mean when the
law of \( x \) is symmetric.
:::

