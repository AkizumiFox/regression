# Principal component regression

Ridge regression keeps every principal direction and shrinks each by a factor between zero and
one. Principal component regression makes the factors zero or one: it keeps the directions in
which the regressors vary most and discards the rest. Principal components (Hotelling 1933) are an old
remedy for collinearity (Massy 1965), and @prp-cmp-tsvd computed the bias and variance of
regression on the leading ones in uncentred form. This section proves when it beats least squares
and makes precise the warning that goes with it.

## Definition

Keep the conventions of [Section 27.2](02-ridge.html): \( \X \) is \( n\times p \), centred, of full column
rank, with its columns on a common scale, and \( \X=\bU\bD\V\T \) with \( d_1\ge\dots\ge d_p>0 \).

::: {#def-shr-pcr}
[Principal component regression]

The **principal components** of the regressors are the vectors \( \bz_j=\X\bv_j=d_j\bu_j \),
\( j=1,\dots,p \). Their sample variances are \( d_j^2/(n-1) \), they are mutually orthogonal, and
\( \bv_j \) is the \( j \)th **loading vector**. For \( 0\le k\le p \), the **principal component
regression** estimate with \( k \) components is
\[
\hbeta^{\mathrm{PC}}_k=\sum_{j\le k}\hat{\alpha}_j\bv_j=\V_k\V_k\T\hbeta,\qquad \V_k=[\bv_1,\dots,\bv_k],
\]
obtained by regressing \( \y \) on \( \bz_1,\dots,\bz_k \) and transforming the coefficients back.
:::

Because the \( \bz_j \) are orthogonal, the coefficient of \( \bz_j \) in the regression on
\( \bz_1,\dots,\bz_k \) is \( \bz_j\T\y/\norm{\bz_j}^2=c_j/d_j=\hat{\alpha}_j \), the same for every \( k\ge j \);
the regression on the leading components is least squares with the trailing canonical
coordinates set to zero. The shrinkage factors are \( f_j=1 \) for \( j\le k \) and \( f_j=0 \) for
\( j>k \).

::: {#prp-shr-pcr}
[Bias and variance of principal component regression]

Let \( \E(\Y)=\X\bbeta \), \( \Cov(\Y)=\sigma^2\I \), \( \boldsymbol{\upalpha}=\V\T\bbeta \) and
\( \tau_j^2=d_j^2\alpha_j^2/\sigma^2 \).

::: {.enumerate options="label=(\alph*)"}
1. \( \hbeta^{\mathrm{PC}}_k \) is the least squares estimate subject to the \( p-k \) restrictions
   \( \bv_j\T\bb=0 \), \( j>k \). Its fitted vector is \( \M_k\y \), where \( \M_k=\sum_{j\le k}\bu_j\bu_j\T \) is
   the orthogonal projection onto \( \spn(\bz_1,\dots,\bz_k) \).

2. \( \E\hbeta^{\mathrm{PC}}_k-\bbeta=-\sum_{j>k}\alpha_j\bv_j \) and
   \( \Cov(\hbeta^{\mathrm{PC}}_k)=\sigma^2\sum_{j\le k}d_j^{-2}\bv_j\bv_j\T \). The total mean squared
   error is \( \sum_{j>k}\alpha_j^2+\sigma^2\sum_{j\le k}d_j^{-2} \), and the mean squared error of the
   fitted values is \( \E\norm{\X\hbeta^{\mathrm{PC}}_k-\X\bbeta}^2=k\sigma^2+\sum_{j>k}d_j^2\alpha_j^2 \).

3. \( \hbeta^{\mathrm{PC}}_k \) is at least as good as \( \hbeta \) in the matrix sense iff
   \( \sum_{j>k}\tau_j^2\le1 \).

4. Among all estimators \( \tilde{\bbeta}_{\mathbf{f}} \) with every \( f_j\in\{0,1\} \), both mean squared errors in
   (b) are minimized by keeping exactly the directions with \( \tau_j^2>1 \) (ties may go either way).
   Principal component regression with \( k \) components is this best selection iff
   \( \tau_j^2\ge1 \) for \( j\le k \) and \( \tau_j^2\le1 \) for \( j>k \).
:::
:::

::: {.proof}
(a) The restrictions say that \( \bb\in\spn(\bv_1,\dots,\bv_k) \), that is, \( \bb=\V_k\mathbf{g} \). Then
\( \X\bb=\bU_k\bD_k\mathbf{g} \) with \( \bU_k,\bD_k \) the leading blocks, and
\( \norm{\y-\bU_k\bD_k\mathbf{g}}^2 \) is minimized at \( \bD_k\mathbf{g}=\bU_k\T\y \), so
\( \mathbf{g}=(\hat{\alpha}_1,\dots,\hat{\alpha}_k)\T \) and the fit is \( \bU_k\bU_k\T\y=\M_k\y \).
(b) The moments are @prp-cmp-tsvd, and the mean squared errors are @prp-shr-canonical(a) and (b)
with the factors \( 0 \) and \( 1 \).
(c) By @thm-shr-bias-variance(c) with the saving \( \boldsymbol{\Delta}_{\A}=\sigma^2\sum_{j>k}d_j^{-2}\bv_j\bv_j\T \), which is
nonnegative definite with column space \( \spn(\bv_{k+1},\dots,\bv_p) \), and bias
\( \boldsymbol{\updelta}=-\sum_{j>k}\alpha_j\bv_j \), which lies in that column space. The Moore–Penrose inverse
is \( \boldsymbol{\Delta}_{\A}^+=\sigma^{-2}\sum_{j>k}d_j^2\bv_j\bv_j\T \), so \( \boldsymbol{\updelta}\T\boldsymbol{\Delta}_{\A}^+\boldsymbol{\updelta}=\sum_{j>k}d_j^2\alpha_j^2/\sigma^2 \).
(d) The losses are sums over directions, and by @prp-shr-canonical(d) the \( j \)th term is smaller
with \( f_j=0 \) than with \( f_j=1 \) iff \( \tau_j^2<1 \).
:::

By part (c), discarding directions improves the estimate of *every* linear combination exactly
when the total signal they carry, in noise units, is at most one. Part (b) shows that the two losses charge a dropped direction differently: for
\( \bbeta \) it costs \( \alpha_j^2 \) and saves \( \sigma^2/d_j^2 \), which is large when \( d_j \) is small; for the
mean it costs \( d_j^2\alpha_j^2 \) and saves only \( \sigma^2 \). Collinearity is a problem for coefficients, not
for fitted values (@prp-col-prediction).

## Low-variance components need not be irrelevant

Part (d) is the crux. The right directions to drop are those with small \( \tau_j^2=d_j^2\alpha_j^2/\sigma^2 \);
principal component regression drops those with small \( d_j^2 \). The two orderings agree only if
the \( \alpha_j \) do not grow too fast as \( d_j \) shrinks, and nothing in the model says they should: the
components are computed from \( \X \) alone. Jolliffe (1982) gave published examples in which
low-variance components were important predictors, and Hadi and Ling (1998) showed how a component
with a tiny variance can carry most of a regression (@exr-shr-hidden-component).

The signal-to-noise ratios can be estimated: under normal errors \( c_j^2/\sigma^2-1 \) is unbiased for
\( \tau_j^2 \), and the component \( t \) statistics \( t_j=c_j/s \) estimate \( \tau_j \) (@exr-shr-component-t).
Part (d) then suggests keeping the components with \( t_j^2>2 \). Such a supervised selection is not a
fixed projection, and the leave-one-out and standard-error formulas do not apply to it without
adjustment ([Chapter 29](../ch29-model-selection/index.html)).

::: {#exm-shr-longley-pcr}
[Principal components of Longley's regressors]

The six standardized Longley regressors have principal components whose shares of the total
variance are \( 76.723\% \), \( 19.589\% \), \( 3.390\% \), \( 0.249\% \), \( 0.043\% \) and \( 0.006\% \), so the usual
rules of thumb keep two or three. The component \( t \) statistics are \( 42.66 \),
\( 5.40 \), \( 10.66 \), \( -0.55 \), \( 3.96 \) and \( 1.72 \) ([Figure 27.3.1](#fig-shr-pcr)(a)). The fifth
component, with less than a two-thousandth of the variance, is nearly as strongly related to
employment as the second, while the fourth, with more variance than the fifth and sixth
together, is irrelevant.

Because principal component regression fits a projection (@prp-shr-pcr(a)), the leave-one-out
prediction errors follow from @thm-res-deletion with the leverages
\( 1/n+\sum_{j\le k}u_{ij}^2 \). The root mean squared leave-one-out errors for \( k=0,1,\dots,6 \) are
\( 3627 \), \( 1109 \), \( 1067 \), \( 506 \), \( 606 \), \( 435 \) and \( 425 \), in thousands. Adding the fourth
component makes predictions *worse*, and adding the fifth makes them much better; the best
number of leading components is \( 6 \), which is least squares itself. The rule \( t_j^2>2 \) keeps
components \( 1,2,3,5,6 \). Applied once to all sixteen years and then treated as fixed, it gives
\( 392 \), but that figure ignores the selection. Redoing the selection within each fold (*nested*
leave-one-out) gives \( 444 \), worse than least squares. Searching all \( 64 \) subsets of components
picks \( 1,2,3,5 \) with \( 389 \), optimistic for the same reason. With three leading
components the coefficients are \( 1023 \), \( 1260 \), \( -1085 \), \( -417 \), \( 1070 \) and \( 966 \): all six
regressors are retained, with the GNP and YEAR coefficients now positive and comparable.
:::

::: {when-format="html"}
![**Figure 27.3.1.** Principal components of Longley's regressors. (a) Share of the regressors'
variance (log scale, left) and the absolute component \( t \) statistic (right; dashed line at
\( \sqrt2 \)). (b) Leave-one-out root mean squared prediction error of principal component
regression with the leading \( k \) components, of ridge regression against its effective degrees of
freedom, and of the components with \( t_j^2>2 \), each for a fixed choice of components or \( \lambda \).](pcr_longley.svg){#fig-shr-pcr width=100%}
:::

::: {when-format="pdf"}
![Principal components of Longley's regressors. (a) Share of the regressors'
variance (log scale, left) and the absolute component \( t \) statistic (right; dashed line at
\( \sqrt2 \)). (b) Leave-one-out root mean squared prediction error of principal component
regression with the leading \( k \) components, of ridge regression against its effective degrees of
freedom, and of the components with \( t_j^2>2 \), each for a fixed choice of components or \( \lambda \).](pcr_longley.pdf){width=100%}
:::

```{.python .run #cell-pcr-pcr}
import numpy as np
import statsmodels.api as sm
data = sm.datasets.longley.load_pandas().data
names = ["GNPDEFL", "GNP", "UNEMP", "ARMED", "POP", "YEAR"]
Z = data[names].to_numpy()
y = data["TOTEMP"].to_numpy()
n, p = Z.shape
X = (Z - Z.mean(axis=0)) / Z.std(axis=0, ddof=1)
yc = y - y.mean()
U, d, Vt = np.linalg.svd(X, full_matrices=False)
c = U.T @ yc                                   # c_j = u_j^T y, the fitted coordinate of component j

def fit_components(keep):
    """Least squares on the chosen component scores: coefficients and leave-one-out error."""
    keep = list(keep)
    beta = Vt[keep].T @ (c[keep] / d[keep])
    resid = yc - X @ beta
    h = 1 / n + np.sum(U[:, keep] ** 2, axis=1)            # leverages, intercept included
    return beta, np.mean((resid / (1 - h)) ** 2)

resid_ls = yc - U @ c
s = np.sqrt(resid_ls @ resid_ls / (n - p - 1))
t = c / s                                        # t statistic of each component
share = d ** 2 / np.sum(d ** 2)                  # share of the regressors' variance
for k in range(p + 1):
    beta_k, loo_k = fit_components(range(k))
    print(f"k = {k}  LOO = {loo_k:9.0f}   coefficients", np.round(beta_k, 0))
print("variance shares:", np.round(share, 5))
print("component t statistics:", np.round(t, 2))
```

## Ridge and principal components compared

Both methods order the directions by \( d_j \) alone, ridge softly and principal components hard,
and neither dominates. If the trailing canonical
coefficients are zero, principal components with the right \( k \) is close to the oracle, while ridge
must leave variance in the trailing directions or bias the leading ones; if all the \( \alpha_j^2 \) are equal, ridge with \( \lambda=\sigma^2/\alpha_j^2 \) *is* the oracle and no
zero–one choice matches it (@exr-shr-pcr-vs-ridge). In the simulations of Frank and Friedman (1993)
ridge was generally at least as good as principal components. For Longley's data neither
method gains over least squares (leading components never do, and tuned ridge has nested
leave-one-out error \( 441 \) against \( 425 \)), because the fifth direction matters and both treat it
as close to noise. Partial least squares, which builds its directions from \( \y \) as well as \( \X \),
is one response (@exr-shr-pls).

## Exercises

### A. Check your understanding

::: {#exr-shr-pcr-not-selection}
[A1]

Show that for \( k<p \) the estimate \( \hbeta^{\mathrm{PC}}_k \) generally has no zero coordinates, so that
principal component regression reduces dimension without selecting variables. When does it set
the coefficient of regressor \( j \) to zero?
:::

::: {#exr-shr-pcr-projection}
[A2]

Show that the fitted values of principal component regression with \( k \) components, with an
unpenalized intercept, are \( (n^{-1}\bone\bone\T+\M_k)\y \), that this matrix is an orthogonal projection of
rank \( k+1 \), and that its leverages are \( 1/n+\sum_{j\le k}u_{ij}^2 \). Why does @thm-res-deletion apply
to principal component regression but not to ridge regression?
:::

### B. Practice

::: {#exr-shr-component-t}
[B1]

Under the normal linear model with \( \rank(\X)=p<n-1 \) (centred data), show that the
\( c_j=\bu_j\T\Y \) are independent \( \Normal(d_j\alpha_j,\sigma^2) \), independent of the residual sum of
squares, and that \( t_j=c_j/s \) has a noncentral \( t \) distribution with \( n-p-1 \) degrees of freedom and
noncentrality \( \tau_j \). Show that \( c_j^2/\sigma^2-1 \) is unbiased for \( \tau_j^2 \).
:::

::: {.solution}
\( \bU\T\Y \) is normal with mean \( \bU\T\X\bbeta=\bD\V\T\bbeta=\bD\boldsymbol{\upalpha} \) and covariance
\( \sigma^2\bU\T\bU=\sigma^2\I \), so the \( c_j \) are independent with the stated laws. They are functions of
the fitted values \( \M\Y=\bU\bU\T\Y \), which are independent of the residual sum of squares (@thm-opt-sampling(c)); with the intercept removed by centring the residual degrees of freedom are
\( n-p-1 \). Hence \( t_j=(c_j/\sigma)/(s/\sigma) \) is a normal with mean \( \tau_j \) (up to sign) divided by an
independent \( \sqrt{\chi^2(n-p-1)/(n-p-1)} \), which is the noncentral \( t \) (@def-qf-noncentral-t). Finally
\( \E(c_j^2/\sigma^2)=\Var(c_j/\sigma)+\{\E(c_j/\sigma)\}^2=1+\tau_j^2 \).
:::

::: {#exr-shr-hidden-component}
[B2]

Let \( x_1,x_2 \) be standardized regressors with correlation \( \rho=0.99 \) and let
\( \E(Y)=x_1-x_2 \). Find the canonical coefficients and show that \( \alpha_1=0 \). Compute \( \tau_1^2 \) and
\( \tau_2^2 \) with \( n=100 \), unit-variance columns and \( \sigma=1 \). Which principal component
regression is best, and what does a one-component principal component regression estimate?
:::

::: {.solution}
With \( \X\T\X=(n-1)\mathbf{R} \), the eigenvectors are \( \bv_1=(1,1)\T/\sqrt2 \) and \( \bv_2=(1,-1)\T/\sqrt2 \),
with \( d_1^2=(n-1)(1+\rho) \) and \( d_2^2=(n-1)(1-\rho) \). For \( \bbeta=(1,-1)\T \),
\( \alpha_1=\bv_1\T\bbeta=0 \) and \( \alpha_2=\sqrt2 \). So \( \tau_1^2=0 \) and
\( \tau_2^2=(n-1)(1-\rho)\cdot2/\sigma^2=99\times0.01\times2=1.98 \). By @prp-shr-pcr(d), the best zero–one
selection keeps only the *second* component. Principal component regression with one component
keeps the first, estimates \( \alpha_1=0 \) and has expected value \( \bzero \): it discards the whole
regression.
:::

::: {#exr-shr-pcr-vs-ridge}
[B3]

Consider a design with two principal directions. (i) If \( \alpha_2=0\ne\alpha_1 \), show that the
oracle factors of @prp-shr-canonical(c) are \( (\tau_1^2/(1+\tau_1^2),0) \), that no ridge estimate attains
them, and that principal component regression with \( k=1 \) exceeds the oracle's total mean squared
error by only \( (\sigma^2/d_1^2)/(1+\tau_1^2) \). (ii) If \( \alpha_1^2=\alpha_2^2>0 \), show that ridge with
\( \lambda=\sigma^2/\alpha_1^2 \) attains the oracle, and that every principal component regression is strictly worse.
:::

::: {.solution}
Part (i): \( \tau_2^2=0 \) gives \( f_2^*=0 \). Ridge has \( f_2(\lambda)=d_2^2/(d_2^2+\lambda)>0 \) for every finite
\( \lambda \), so it never attains \( f_2^*=0 \). Principal components with \( k=1 \) has factors \( (1,0) \). By
@prp-shr-canonical(a) its first term is \( \sigma^2/d_1^2 \) against the oracle's
\( (\sigma^2/d_1^2)\tau_1^2/(1+\tau_1^2) \), and the second terms are both zero; the difference is
\( (\sigma^2/d_1^2)/(1+\tau_1^2) \). (ii) With \( \alpha_1^2=\alpha_2^2=a \) and \( \lambda=\sigma^2/a \), ridge has factors
\( d_j^2/(d_j^2+\sigma^2/a)=d_j^2a/(d_j^2a+\sigma^2)=\tau_j^2/(1+\tau_j^2) \), the oracle. Each term of
@prp-shr-canonical(a) is a strictly convex quadratic in \( f_j \), so the oracle is its unique minimizer, and
since \( 0<\tau_j^2/(1+\tau_j^2)<1 \), every choice of factors in \( \{0,1\} \) is strictly worse.
:::

### C. Going deeper

::: {#exr-shr-pls}
[C1]

The first direction of **partial least squares** is \( \bw_1\propto\X\T\y \), the direction maximizing
the sample covariance of \( \X\bw \) with \( \y \) over unit vectors \( \bw \). Show that
\( \bw_1\propto\sum_jd_jc_j\bv_j \), and compare with the first principal direction \( \bv_1 \). Show that
the one-dimensional fit along \( \bw_1 \) never ignores a direction \( j \) with \( c_j\ne0 \).
:::

