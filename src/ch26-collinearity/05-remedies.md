# What to do about collinearity

Refitting the same data in another way cannot supply information the design lacks, and every genuine
remedy does one of four things: it
gets more information of the right kind, changes the question to one the data can answer, adds
information from outside the data, or accepts bias in exchange for variance. Before choosing among
them, decide whether anything needs to be done at all.

## Does it matter?

By @thm-col-variance(c) and @prp-col-prediction, collinearity is harmless for quantities with no
component along the weak directions: fitted values, predictions that follow the pattern of the data, and
combinations such as the sum in @exm-col-two-regressors. If those are the targets, report them and stop.

## Collect better data

The only remedy that keeps the question is more information along the weak directions. The next
proposition shows how much a single well-placed observation can supply.

::: {#prp-col-new-observation}
[One new observation]

Let \( \A=\X\T\X \) be positive definite with eigenvalues \( \lambda_1\ge\dots\ge\lambda_p \) and orthonormal
eigenvectors \( \bv_1,\dots,\bv_p \), \( p\ge2 \), and let a case with regressor vector \( \x \) be added, so that \( \A \)
becomes \( \A+\x\x\T \). Write \( \G=\A^{-1} \).

::: {.enumerate options="label=(\alph*)"}
1. For every \( \mathbf{a} \), the variance of \( \mathbf{a}\T\hbeta \) falls by \( \sigma^2(\mathbf{a}\T\G\x)^2/(1+\x\T\G\x) \); in particular it never
   increases.

2. If \( \norm{\x}\le c \), then \( \lambda_{\min}(\A+\x\x\T)\le\min(\lambda_p+c^2,\ \lambda_{p-1}) \), with equality when \( \x=c\bv_p \).

3. Among all \( \x \) with \( \norm{\x}=c \), the choice \( \x=c\bv_p \) minimizes the largest variance of \( \mathbf{a}\T\hbeta \) over unit
   vectors \( \mathbf{a} \).
:::

:::

::: {.proof}
(a) By @eq-mat-sherman-morrison, \( (\A+\x\x\T)^{-1}=\G-\G\x\x\T\G/(1+\x\T\G\x) \).

(b) By @thm-mat-extremal-rayleigh, \( \lambda_{\min}(\A+\x\x\T)\le\bv_p\T(\A+\x\x\T)\bv_p=\lambda_p+(\x\T\bv_p)^2\le\lambda_p+c^2 \). For the
second bound, the plane \( \spn\{\bv_{p-1},\bv_p\} \) contains a unit vector \( \bb \) orthogonal to \( \x \), since the linear map
\( \bb\mapsto\x\T\bb \) from the plane to \( \Real \) has a nontrivial null space. Then
\( \lambda_{\min}(\A+\x\x\T)\le\bb\T\A\bb\le\lambda_{p-1} \). If \( \x=c\bv_p \), then \( \A+\x\x\T \) has the same eigenvectors as \( \A \) and
eigenvalues \( \lambda_1,\dots,\lambda_{p-1},\lambda_p+c^2 \), whose minimum is \( \min(\lambda_p+c^2,\lambda_{p-1}) \).

(c) By @thm-col-variance(b), the largest variance over unit vectors is \( \sigma^2/\lambda_{\min} \), and (b) says that
\( \x=c\bv_p \) maximizes \( \lambda_{\min} \).
:::

The new case should lie along the weakest direction, where the existing data are thin. Such a point
breaks the pattern of the regressors, and in observational data it may not exist; in an experiment, or
when the sample can be chosen, the proposition says what to look for.

::: {#exm-col-new-point}
[The value of one observation]

In the design of @exm-col-two-regressors, add one case at the centred values
\( (x_1,x_2)=(1.077,-1.044) \), the point at distance \( 1.5 \) from the centroid along the weak slope direction.
Each coordinate is well inside the range of its regressor; only the combination is new. The slope
estimates are governed by the slope block with the intercept partialled out (@thm-proj-fwl), and the
new case, which has an intercept entry, changes that block from \( \mathbf{S} \) to
\( \mathbf{S}+\tfrac{n}{n+1}\x\x\T \), where \( \x \) holds its two slope entries: the columns must be centred again about
the new mean, and that costs the factor \( n/(n+1) \). Applied to this block,
@prp-col-new-observation(b) predicts that the smaller eigenvalue rises from \( 0.172 \) to
\( 0.172+1.5^2\cdot30/31=2.349 \), which is what the script finds. The standard deviation of \( \hat{\beta}_1 \) falls from \( 1.734\,\sigma \) to
\( 0.476\,\sigma \). To achieve the same reduction with more data of the old pattern, the whole sample would have
to be repeated \( 13.3 \) times, about \( 370 \) further cases in place of one. The script also checks, on the
centred block, that no other direction does better.
:::

```{.python .run #cell-remedies-design}
import numpy as np
rng = np.random.default_rng(2601)                     # the design of Section 26.1
n = 30
x1 = rng.normal(size=n)
x2 = x1 + 0.15 * rng.normal(size=n)
X = np.column_stack([np.ones(n), x1 - x1.mean(), x2 - x2.mean()])
G = np.linalg.inv(X.T @ X)
lam, V = np.linalg.eigh(X.T @ X)                      # ascending eigenvalues
v_weak = V[:, 0]                                      # the weak direction (it has no intercept part)
```

```{.python .run #cell-remedies-new-point}
c = 1.5
x_new = np.r_[1.0, c * v_weak[1:] / np.linalg.norm(v_weak[1:])]   # breaks the pattern
print("new point (centred x1, x2):", np.round(x_new[1:], 3))
G_new = np.linalg.inv(X.T @ X + np.outer(x_new, x_new))
print("sd of b1 before / after:", round(np.sqrt(G[1, 1]), 3), round(np.sqrt(G_new[1, 1]), 3))
S_old = X[:, 1:].T @ X[:, 1:]                         # slope block: the intercept is orthogonal
xs = x_new[1:]
S_new = S_old + n / (n + 1) * np.outer(xs, xs)        # centred slope block after recentring
print("smallest slope eigenvalue before / after:", round(np.linalg.eigvalsh(S_old)[0], 3),
      round(np.linalg.eigvalsh(S_new)[0], 3))
# repeating the old pattern k times multiplies X'X by k: Var(b1) falls by the factor 1/k
k_equiv = G[1, 1] / G_new[1, 1]
print("equivalent number of copies of the whole data set:", round(k_equiv, 1))
```

## Change the question

A reparameterization \( \X\bbeta=(\X\mathbf{T})(\mathbf{T}^{-1}\bbeta) \) leaves the least squares estimate of every linear
function of \( \bbeta \) unchanged (@thm-proj-reparam), so it cannot make \( \beta_1-\beta_2 \) better determined. It can make
the reported parameters the well-determined ones. Replacing \( x_1 \) and \( x_2 \) by their average and their
difference gives a coefficient for the average that estimates \( \beta_1+\beta_2 \) precisely, and leaves all the
imprecision in the coefficient of the difference. Centring polynomial terms (@exm-col-quadratic) is the same move.

Two common changes of model are really linear restrictions on \( \bbeta \), and should be judged as such.
*Dropping a regressor* imposes \( \beta_2=0 \). By @thm-dep-omitted the short slope estimates \( \beta_1+\pi\beta_2 \), where
\( \pi \) is the slope of \( x_2 \) on \( x_1 \); in @exm-col-two-regressors \( \pi=1.0261 \), and the short slope has standard
deviation \( 0.171\,\sigma \) against \( 1.734\,\sigma \). It is precise because it answers a different question, giving \( x_1 \)
the credit for both regressors: with true slopes \( (1,1) \) it centres at about \( 2.03 \). @thm-dep-mse gives the
condition under which it has the smaller mean squared error. *Combining regressors into an index*
imposes a restriction such as \( \beta_1=\beta_2 \). It concerns the weak direction \( \beta_1-\beta_2 \), whose estimate has
standard deviation \( 3.412\,\sigma \) in the example, so an \( F \) test of it has almost no power: the restriction
is cheap in variance and nearly untestable, for the same reason.

## Add information from outside the data

An exact restriction known from theory is fitted by restricted least squares (@prp-glh-restricted-ls).
An approximate one can be written as extra observations \( \mathbf{r}=\bL\bbeta+\mathbf{u} \) and combined with the data
by generalized least squares, the mixed estimation of Theil and Goldberger (1961) (@exr-col-mixed). The Bayesian version is a prior: with the conjugate prior of @thm-opt-bayes-conjugate, under which
\( \bbeta \) given \( \sigma^2 \) has covariance \( \sigma^2\V_0 \), the posterior mean solves equations with matrix
\( \V_0^{-1}+\X\T\X \). A prior informative along the weak eigenvector raises the smallest eigenvalue as a new observation
does, and the answer along that direction then comes mostly from the prior.

## Accept bias

The last option gives up unbiasedness. Principal component regression drops the canonical coordinates
with small eigenvalues, the truncated SVD estimator of @prp-cmp-tsvd; ridge regression shrinks every
canonical coordinate towards zero, the weak ones most. [Chapter 27](../ch27-shrinkage/index.html)
develops both (@thm-shr-ridge and @def-shr-pcr). Neither creates information: they trade variance along
the weak directions for bias along the same directions.

## What not to do

- Do not drop regressors because their variance inflation factors exceed a threshold; the fit may
  hardly change while the meaning of every remaining coefficient does.
- Do not expect centring or rescaling to cure collinearity that is not created by the parameterization.
- Do not let an automatic search choose among collinear regressors: noise decides which enters
  ([Chapter 29](../ch29-model-selection/index.html)).

## Exercises

### A. Check your understanding

::: {#exr-col-reparam-cannot}
[A1]

Explain, using @thm-col-variance, why no reparameterization of @exm-col-two-regressors can estimate
\( \beta_1-\beta_2 \) with a standard deviation smaller than \( 3.412\,\sigma \), although the coefficient of the new
difference regressor may have any standard deviation you like.
:::

::: {#exr-col-short-credit}
[A2]

In @exm-col-two-regressors, what does the short regression of \( y \) on \( x_1 \) alone estimate if the true slopes are
\( (2,0) \)? If they are \( (0,2) \)? What does this say about using the short regression to decide which regressor
matters?
:::

### B. Practice

::: {#exr-col-mixed}
[B1]

Suppose that, besides \( \Y=\X\bbeta+\be \) with \( \Cov(\be)=\sigma^2\I \), there is independent information
\( \mathbf{r}=\bL\bbeta+\mathbf{u} \) with \( \E(\mathbf{u})=\bzero \) and \( \Cov(\mathbf{u})=\sigma^2\boldsymbol{\Omega} \), \( \boldsymbol{\Omega} \) known and positive
definite. Show that the generalized least squares estimator from the stacked data is
\[
\tilde{\bbeta}=(\X\T\X+\bL\T\boldsymbol{\Omega}^{-1}\bL)^{-1}(\X\T\Y+\bL\T\boldsymbol{\Omega}^{-1}\mathbf{r}),
\]
that it is unbiased with covariance \( \sigma^2(\X\T\X+\bL\T\boldsymbol{\Omega}^{-1}\bL)^{-1} \), and that it equals the posterior mean of
@thm-opt-bayes-conjugate for a suitable prior. With \( \bL=\bv_p\T \), how does the smallest eigenvalue change?
:::

::: {.solution}
Stack \( (\Y\T,\mathbf{r}\T)\T=(\X\T,\bL\T)\T\bbeta+(\be\T,\mathbf{u}\T)\T \), with error covariance
\( \sigma^2\diag(\I,\boldsymbol{\Omega}) \). Generalized least squares (@cor-opt-aitken) gives the stated estimator, and it is unbiased with
covariance \( \sigma^2 \) times the inverse of the weighted cross-product matrix. If \( \bL \) is square and nonsingular, the prior
\( \bbeta\sim\Normal(\bL^{-1}\mathbf{r},\sigma^2\bL^{-1}\boldsymbol{\Omega}(\bL^{-1})\T) \) has precision \( \bL\T\boldsymbol{\Omega}^{-1}\bL/\sigma^2 \) and mean
\( \bL^{-1}\mathbf{r} \), and the posterior mean of @thm-opt-bayes-conjugate is the same expression. With \( \bL=\bv_p\T \) and
\( \boldsymbol{\Omega}=\omega^2 \), the matrix becomes \( \X\T\X+\bv_p\bv_p\T/\omega^2 \), whose smallest eigenvalue is
\( \min(\lambda_p+1/\omega^2,\lambda_{p-1}) \): exactly the effect of a new observation of length \( 1/\omega \) along \( \bv_p \) (@prp-col-new-observation(b)).
:::

::: {#exr-col-composite}
[B2]

Let \( x_1,x_2 \) be centred with \( \norm{\x_1}^2=S_1 \), \( \norm{\x_2}^2=S_2 \) and \( \x_1\T\x_2=S_{12} \), and fit the intercept and the
single composite \( \mathbf{s}=\x_1+\x_2 \). Show that its coefficient has expectation
\( \beta_1+(\beta_2-\beta_1)(S_2+S_{12})/(S_1+S_2+2S_{12}) \). When is this \( (\beta_1+\beta_2)/2 \)?
:::

::: {.solution}
The coefficient is \( \mathbf{s}\T\Y/\mathbf{s}\T\mathbf{s} \), because \( \mathbf{s} \) is centred and so orthogonal to \( \bone \). Its expectation is
\( \mathbf{s}\T(\x_1\beta_1+\x_2\beta_2)/\mathbf{s}\T\mathbf{s} \), with \( \mathbf{s}\T\x_1=S_1+S_{12} \), \( \mathbf{s}\T\x_2=S_2+S_{12} \) and
\( \mathbf{s}\T\mathbf{s}=S_1+S_2+2S_{12} \). Writing \( \beta_1(S_1+S_{12})=\beta_1(\mathbf{s}\T\mathbf{s})-\beta_1(S_2+S_{12}) \) gives the formula. It is
\( (\beta_1+\beta_2)/2 \) for all \( \bbeta \) iff \( S_2+S_{12}=(S_1+S_2+2S_{12})/2 \), that is, iff \( S_1=S_2 \): the composite estimates the
average slope when the two regressors have equal spread, whatever their correlation.
:::

::: {#exr-col-best-direction}
[B3]

For \( p=2 \) and \( \A=\diag(\lambda_1,\lambda_2) \) with \( \lambda_1>\lambda_2 \), let \( \x=c(\cos\theta,\sin\theta)\T \). Compute
\( \lambda_{\min}(\A+\x\x\T) \) and show directly that it is maximized at \( \theta=\pi/2 \).
:::

::: {.solution}
The trace is \( T=\lambda_1+\lambda_2+c^2 \), independent of \( \theta \), and the determinant is
\( D=(\lambda_1+c^2\cos^2\theta)(\lambda_2+c^2\sin^2\theta)-c^4\cos^2\theta\sin^2\theta=\lambda_1\lambda_2+c^2(\lambda_1\sin^2\theta+\lambda_2\cos^2\theta) \).
The smaller eigenvalue is \( \bigl(T-\sqrt{T^2-4D}\bigr)/2 \), which increases with \( D \) for fixed \( T \). Since \( \lambda_1>\lambda_2 \), \( D \)
is largest when \( \sin^2\theta=1 \).
:::

### C. Going deeper

::: {#exr-col-two-points}
[C1]

@prp-col-new-observation(c) says where to put *one* new case. Repeating that advice for two cases is not
optimal. Let \( p=2 \), \( \A=\diag(\lambda_1,\lambda_2) \) with \( \lambda_1>\lambda_2>0 \), and add two
rows \( \x \) and \( \tilde{\x} \), each of length \( c \), so that \( \A \) becomes
\( \A+\x\x\T+\tilde{\x}\tilde{\x}\T \).

::: {.enumerate options="label=(\alph*)"}
1. Show that both cases along the weak direction give
   \( \lambda_{\min}=\min(\lambda_1,\ \lambda_2+2c^2) \), and that one along each eigenvector gives
   \( \lambda_2+c^2 \).

2. Show that for *every* such pair,
   \[
   \lambda_{\min}\le\min\Bigl\{\lambda_2+2c^2,\ \tfrac12(\lambda_1+\lambda_2)+c^2\Bigr\},
   \]
   using \( \bv_2\T(\cdot)\bv_2 \) for the first bound and the trace for the second.

3. Show that the bound is attained. When \( 2c^2\le\lambda_1-\lambda_2 \) it is attained by (a)'s first
   choice. When \( 2c^2>\lambda_1-\lambda_2 \) it is attained by
   \[
   \x=\bigl(\sqrt{a_1/2},\ \sqrt{a_2/2}\bigr)\T,\qquad
   \tilde{\x}=\bigl(\sqrt{a_1/2},\ -\sqrt{a_2/2}\bigr)\T,
   \qquad a_{1,2}=c^2\mp\tfrac12(\lambda_1-\lambda_2),
   \]
   neither of which lies along an eigenvector.

4. Take \( \lambda_1=2 \), \( \lambda_2=1 \) and \( c=3/2 \). Compute \( \lambda_{\min} \) for the three
   designs of (a) and (c), and compare.
:::

:::

::: {.solution}
(a) If \( \x=\tilde{\x}=c\bv_2 \), the sum \( \A+2c^2\bv_2\bv_2\T \) is diagonal with entries
\( \lambda_1 \) and \( \lambda_2+2c^2 \). If \( \x=c\bv_1 \) and \( \tilde{\x}=c\bv_2 \), it is diagonal with
entries \( \lambda_1+c^2 \) and \( \lambda_2+c^2 \), whose minimum is \( \lambda_2+c^2 \) because
\( \lambda_1>\lambda_2 \).

(b) Write \( \M=\x\x\T+\tilde{\x}\tilde{\x}\T \), which is nonnegative definite with
\( \tr\M=2c^2 \). By @thm-mat-extremal-rayleigh,
\( \lambda_{\min}(\A+\M)\le\bv_2\T(\A+\M)\bv_2=\lambda_2+\bv_2\T\M\bv_2\le\lambda_2+\tr\M=\lambda_2+2c^2 \),
since the eigenvalues of \( \M \) are nonnegative and sum to \( \tr\M \). And the smaller of two eigenvalues
is at most their average, so \( \lambda_{\min}(\A+\M)\le\tfrac12\tr(\A+\M)=\tfrac12(\lambda_1+\lambda_2)+c^2 \).

(c) For the first case the design of (a) gives exactly \( \lambda_2+2c^2 \), which is then the smaller of the
two bounds and is attained. For the second, \( a_1 \) and \( a_2 \) are nonnegative because
\( 2c^2>\lambda_1-\lambda_2 \), they sum to \( 2c^2 \), and
\( \x\x\T+\tilde{\x}\tilde{\x}\T=\diag(a_1,a_2) \), the off-diagonal terms cancelling. Also
\( \norm{\x}^2=\norm{\tilde{\x}}^2=(a_1+a_2)/2=c^2 \). Then
\( \A+\M=\diag(\lambda_1+a_1,\ \lambda_2+a_2) \), and both entries equal
\( \tfrac12(\lambda_1+\lambda_2)+c^2 \), so \( \lambda_{\min} \) attains the second bound.

(d) Here \( c^2=9/4 \) and \( \lambda_1-\lambda_2=1<2c^2 \). Both cases along \( \bv_2 \) give
\( \min(2,\ 1+\tfrac92)=2 \), no better than doing nothing to the strong direction. One along each eigenvector
gives \( 1+\tfrac94=\tfrac{13}4=3.25 \). The optimal pair, with \( a_1=\tfrac74 \) and \( a_2=\tfrac{11}4 \),
gives \( \tfrac32+\tfrac94=\tfrac{15}4=3.75 \).

The greedy rule fails because after the first case the weakest direction is no longer \( \bv_2 \), and a plan
made for two cases should anticipate that. Choosing several cases at once is a design problem — here the
criterion \( \lambda_{\min} \), called E-optimality — and one-at-a-time optimization solves it only when the
budget is too small to level the two eigenvalues.
:::

::: {#exr-col-prior-hurts}
[C2]

Outside information helps only if it is good enough. In the setting of @exr-col-mixed take \( \bL=\bv_p\T \)
and \( \boldsymbol{\Omega}=\omega^2 \), and suppose the source is *wrong*: \( \E(r)=\gamma_p+\delta \) with
\( \gamma_p=\bv_p\T\bbeta \) and \( \delta\ne0 \), while \( \Var(r)=\sigma^2\omega^2 \) as claimed, and
\( r \) is independent of \( \Y \). Write \( \kappa=1/\omega^2 \).

::: {.enumerate options="label=(\alph*)"}
1. Show that \( \bv_\ell\T\tilde{\bbeta}=\hat\gamma_\ell \) for \( \ell\ne p \), while
   \( \tilde\gamma_p=(\lambda_p\hat\gamma_p+\kappa r)/(\lambda_p+\kappa) \). Only the weak coordinate moves.

2. Show that \( \tilde\gamma_p \) has bias \( \kappa\delta/(\lambda_p+\kappa) \) and variance
   \( \sigma^2/(\lambda_p+\kappa) \), so that
   \[
   \E(\tilde\gamma_p-\gamma_p)^2=\frac{\sigma^2}{\lambda_p+\kappa}
   +\frac{\kappa^2\delta^2}{(\lambda_p+\kappa)^2}.
   \]

3. Deduce that the outside information lowers the mean squared error if and only if
   \[
   \delta^2<\sigma^2\omega^2+\frac{\sigma^2}{\lambda_p}=\Var(r)+\Var(\hat\gamma_p):
   \]
   its error must be smaller than the two standard errors combined in quadrature.

4. Interpret: what happens when the source is honest, in the sense that \( \delta \) is itself drawn with
   mean zero and variance \( \sigma^2\omega^2 \)? What happens when it is confident and wrong, that is when
   \( \omega^2 \) is small and \( \delta \) is not?
:::

:::

::: {.solution}
(a) By @exr-col-mixed, \( \tilde{\bbeta}=(\X\T\X+\kappa\bv_p\bv_p\T)^{-1}(\X\T\Y+\kappa\bv_pr) \). The matrix
has the eigenvectors of \( \X\T\X \), with eigenvalues \( \lambda_\ell \) for \( \ell\ne p \) and
\( \lambda_p+\kappa \), so its inverse is
\( \sum_{\ell\ne p}\lambda_\ell^{-1}\bv_\ell\bv_\ell\T+(\lambda_p+\kappa)^{-1}\bv_p\bv_p\T \). Since
\( \bv_\ell\T\X\T\Y=\lambda_\ell\hat\gamma_\ell \) (@thm-col-variance(a)) and \( \bv_\ell\T\bv_p=0 \) for
\( \ell\ne p \), the stated coordinates follow.

(b) \( \E\tilde\gamma_p=(\lambda_p\gamma_p+\kappa(\gamma_p+\delta))/(\lambda_p+\kappa) \), whose difference
from \( \gamma_p \) is \( \kappa\delta/(\lambda_p+\kappa) \). Since \( r \) and \( \hat\gamma_p \) are
independent with variances \( \sigma^2/\kappa \) and \( \sigma^2/\lambda_p \) (@thm-col-variance(a)),
\[
\Var(\tilde\gamma_p)=\frac{\lambda_p^2\sigma^2/\lambda_p+\kappa^2\sigma^2/\kappa}{(\lambda_p+\kappa)^2}
=\frac{\sigma^2}{\lambda_p+\kappa}.
\]
Add the squared bias.

(c) The inequality \( \E(\tilde\gamma_p-\gamma_p)^2<\sigma^2/\lambda_p \) reads
\[
\frac{\kappa^2\delta^2}{(\lambda_p+\kappa)^2}<\sigma^2\Bigl(\frac1{\lambda_p}-\frac1{\lambda_p+\kappa}\Bigr)
=\frac{\sigma^2\kappa}{\lambda_p(\lambda_p+\kappa)} .
\]
Multiply by \( (\lambda_p+\kappa)^2/\kappa \) and divide by \( \kappa \):
\( \delta^2<\sigma^2(\lambda_p+\kappa)/(\kappa\lambda_p)=\sigma^2/\kappa+\sigma^2/\lambda_p \), which is the
stated bound because \( \sigma^2/\kappa=\sigma^2\omega^2=\Var(r) \).

(d) If the source is honest then \( \E\delta^2=\sigma^2\omega^2 \), which is strictly less than the bound, so
on average the outside information helps, and it helps most where the data are weakest, \( \lambda_p \) small.
If the source is confident and wrong — \( \omega^2 \) small, so that \( \kappa \) is large and the answer
along \( \bv_p \) comes almost entirely from outside — the bound \( \sigma^2\omega^2+\sigma^2/\lambda_p \) is
itself small, and a modest \( \delta \) makes the mean squared error worse than plain least squares, which at
least was unbiased. A stated prior variance is a claim about \( \delta \), and the inequality says the claim
has to be roughly right, not merely convenient.
:::
