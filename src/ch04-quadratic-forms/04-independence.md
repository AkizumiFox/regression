# Independence of quadratic forms and linear forms

An \( F \) or \( t \) statistic is a ratio, and its distribution is known only if numerator and
denominator are independent. For normal vectors, uncorrelated linear functions are
independent (@thm-mvn-independence). Quadratic forms are not linear, but the
following observation brings them within reach of that result.

::: {.idea}
A quadratic form is a function of a linear form: for any generalized inverse \( \A\ginv \) of
\( \A \),
\[
\Y\T\A\Y=\Y\T\A\A\ginv\A\Y=(\A\Y)\T\A\ginv(\A\Y).
\]
So \( \Y\T\A\Y \) is independent of anything that is independent of the vector \( \A\Y \).
:::

## The independence theorems

::: {#thm-qf-indep-linear}
[Linear and quadratic forms]

Let \( \Y\sim\Normal_n(\bmu,\bSigma) \), let \( \B \) be \( m\times n \) and \( \A \) symmetric \( n\times n \).
If \( \B\bSigma\A=\mathbf{0} \), then \( \B\Y \) and \( \Y\T\A\Y \) are independent. Conversely, if \( \B\Y \) and
\( \Y\T\A\Y \) are independent for every \( \bmu\in\Real^n \), with \( \bSigma \) fixed, then
\( \B\bSigma\A=\mathbf{0} \).
:::

::: {.proof}
The vector \( (\B\Y,\A\Y) \) is a linear function of \( \Y \), hence jointly normal
(@thm-mvn-linear), and its cross-covariance is \( \Cov(\B\Y,\A\Y)=\B\bSigma\A\T=\B\bSigma\A \)
(@thm-rv-linear). If this is zero, \( \B\Y \) and \( \A\Y \) are independent by
@thm-mvn-independence, and so are \( \B\Y \) and the function \( (\A\Y)\T\A\ginv(\A\Y) \) of
\( \A\Y \), which is \( \Y\T\A\Y \) by the key idea above. Conversely, independent variables with
finite variances are uncorrelated, so @thm-qf-mean-var(c) gives
\( 2\B\bSigma\A\bmu=\bzero \) for every \( \bmu \), that is, \( \B\bSigma\A=\mathbf{0} \).
:::

::: {#thm-qf-indep-quadratic}
[Two quadratic forms]

Let \( \Y\sim\Normal_n(\bmu,\bSigma) \) and let \( \A,\B \) be symmetric. If \( \A\bSigma\B=\mathbf{0} \),
equivalently \( \B\bSigma\A=\mathbf{0} \), then \( \Y\T\A\Y \) and \( \Y\T\B\Y \) are independent.
:::

::: {.proof}
The two conditions are transposes of each other. As before, \( \A\Y \) and \( \B\Y \) are jointly
normal with \( \Cov(\A\Y,\B\Y)=\A\bSigma\B=\mathbf{0} \), hence independent, and each form is a
function of one of them.
:::

The condition is sufficient for every \( \bmu \) and every \( \bSigma \). When \( \bSigma \) is singular one
might hope that the weaker requirement \( \bSigma\A\bSigma\B\bSigma=\mathbf{0} \) suffices, since only
directions in \( \C(\bSigma) \) carry noise. For nonnegative definite \( \A \) and \( \B \) the two requirements
turn out to be equivalent (@exr-qf-singular-independence). In the other direction, the
condition is necessary under mild assumptions, a result usually credited to
Craig (1943). Its general proof is surprisingly delicate, and the history of
incorrect published proofs is recounted by Driscoll and Gundberg (1986). The case that
matters most for sums of squares has a short proof.

::: {#prp-qf-craig-nnd}
[Necessity for nonnegative definite forms]

Let \( \Y\sim\Normal_n(\bzero,\I) \) and let \( \A,\B \) be nonnegative definite. If \( \Y\T\A\Y \) and
\( \Y\T\B\Y \) are independent, then \( \A\B=\mathbf{0} \).
:::

::: {.proof}
Independence makes the covariance zero, and by @thm-qf-mean-var(d) the covariance is
\( 2\tr(\A\B) \). With the square roots of @thm-mat-square-root,
\[
\tr(\A\B)=\tr(\A^{1/2}\A^{1/2}\B^{1/2}\B^{1/2})=\tr\bigl\{(\A^{1/2}\B^{1/2})\T(\A^{1/2}\B^{1/2})\bigr\},
\]
the sum of squares of the entries of \( \A^{1/2}\B^{1/2} \). So \( \A^{1/2}\B^{1/2}=\mathbf{0} \), and
\( \A\B=\A^{1/2}(\A^{1/2}\B^{1/2})\B^{1/2}=\mathbf{0} \).
:::

For nonnegative definite forms in a standard normal vector, then, zero correlation and
independence coincide, just as for jointly normal linear forms. This fails for general
symmetric matrices (@exr-qf-uncorrelated-dependent).

## Orthogonal projections of a spherical normal vector

The situation in regression is the tidiest one: a spherical normal vector and several
projections onto mutually orthogonal subspaces.

::: {#thm-qf-orthogonal-projections}
[Independent projections]

Let \( \Y\sim\Normal_n(\bmu,\sigma^2\I) \), and let \( \bP_1,\dots,\bP_k \) be symmetric idempotent
matrices with \( \bP_i\bP_j=\mathbf{0} \) for \( i\ne j \), and ranks \( r_1,\dots,r_k \). Then
\( \bP_1\Y,\dots,\bP_k\Y \) are mutually independent normal vectors,
\( \bP_i\Y\sim\Normal_n(\bP_i\bmu,\sigma^2\bP_i) \), and the sums of squares
\[
\frac{\norm{\bP_i\Y}^2}{\sigma^2}\sim\chi^2\Bigl(r_i,\frac{\norm{\bP_i\bmu}^2}{\sigma^2}\Bigr),
\qquad i=1,\dots,k,
\]
are mutually independent.
:::

::: {.proof}
Stack the vectors \( \bP_i\Y \) into one \( nk \)-vector. It is a linear function of \( \Y \), so it is
normal, and the covariance between blocks \( i \) and \( j \) is
\( \bP_i(\sigma^2\I)\bP_j\T=\sigma^2\bP_i\bP_j \). This is \( \mathbf{0} \) for \( i\ne j \) and \( \sigma^2\bP_i \)
for \( i=j \). By the block form of @thm-mvn-independence the blocks are mutually
independent, and functions of them, the squared lengths, are then mutually independent.
Their laws are @thm-qf-chisq(a).
:::

In geometric language, which [Chapter 6](../ch06-projections/index.html) develops: the coordinates of a spherical
normal vector in mutually orthogonal subspaces are independent. It is the reason the rows
of an analysis of variance table can be treated separately.

::: {#cor-qf-sample-variance}
[Sample mean and sample variance]

Let \( Y_1,\dots,Y_n \) (\( n\ge2 \)) be independent \( \Normal(\mu,\sigma^2) \), with sample mean \( \bar{Y} \)
and sample variance \( S^2=(n-1)^{-1}\sum_i(Y_i-\bar{Y})^2 \). Then

::: {.enumerate options="label=(\alph*)"}
1. \( \bar{Y} \) and \( S^2 \) are independent, and \( \bar{Y}\sim\Normal(\mu,\sigma^2/n) \);

2. \( (n-1)S^2/\sigma^2\sim\chi^2(n-1) \);

3. for any \( \mu_0 \), \( T=\sqrt n(\bar{Y}-\mu_0)/S\sim t(n-1,\delta) \) with
           \( \delta=\sqrt n(\mu-\mu_0)/\sigma \). In particular \( T\sim t(n-1) \) when \( \mu=\mu_0 \).
:::

:::

::: {.proof}
Here \( \Y\sim\Normal_n(\mu\bone,\sigma^2\I) \). Let \( \bP_1=n^{-1}\bone\bone\T \) and
\( \mathbf{C}=\I-\bP_1 \). Both are symmetric and idempotent, \( \bP_1\mathbf{C}=\mathbf{0} \), and
\( \rank(\mathbf{C})=\tr(\mathbf{C})=n-1 \). Since \( \bP_1\Y=\bar{Y}\bone \) and \( \norm{\mathbf{C}\Y}^2=(n-1)S^2 \),
@thm-qf-orthogonal-projections gives (a), which @exm-mvn-mean-deviations obtained
directly, and, because \( \mathbf{C}\mu\bone=\bzero \), also (b).
For (c), \( X=\sqrt n(\bar{Y}-\mu_0)/\sigma\sim\Normal(\delta,1) \) and \( V=(n-1)S^2/\sigma^2\sim\chi^2(n-1) \)
are independent, and \( T=X/\sqrt{V/(n-1)} \). Apply @def-qf-noncentral-t.
:::

The independence in (a) characterizes the normal distribution: for \( n\ge2 \) independent
identically distributed observations with finite variance, \( \bar{Y} \) and \( S^2 \) are independent
only under normality (Lukacs 1942). We will not need this, but it shows that normality is used essentially here.
It is not a technicality that more careful arguments could remove.

## Preview: the \( F \) statistic for nested models

The main use of these results comes in [Chapter 11](../ch11-general-linear-hypothesis/index.html). Its core fits in
a few lines.

::: {#thm-qf-nested-f}
[\( F \) statistic for nested subspaces]

Let \( \Y\sim\Normal_n(\bmu,\sigma^2\I) \), and let \( \bP_0 \) and \( \bP \) be symmetric idempotent
matrices of ranks \( r_0<r<n \) with \( \bP\bP_0=\bP_0 \). Suppose that \( \bP\bmu=\bmu \). Then
\[
F=\frac{\Y\T(\bP-\bP_0)\Y/(r-r_0)}{\Y\T(\I-\bP)\Y/(n-r)}
\sim F\Bigl(r-r_0,\;n-r,\;\frac{\norm{(\bP-\bP_0)\bmu}^2}{\sigma^2}\Bigr),
\]
and the noncentrality is zero iff \( \bP_0\bmu=\bmu \).
:::

::: {.proof}
Transposing \( \bP\bP_0=\bP_0 \) gives \( \bP_0\bP=\bP_0 \). Hence \( (\bP-\bP_0)^2=\bP-2\bP_0+\bP_0=\bP-\bP_0 \),
so \( \bP-\bP_0 \) is symmetric idempotent with rank \( \tr(\bP)-\tr(\bP_0)=r-r_0 \). Also
\( (\bP-\bP_0)(\I-\bP)=\bP-\bP-\bP_0+\bP_0=\mathbf{0} \), and \( \I-\bP \) is symmetric idempotent of rank
\( n-r \). By @thm-qf-orthogonal-projections the numerator and denominator sums of squares
are independent, with \( \sigma^{-2} \) times them distributed as
\( \chi^2(r-r_0,\norm{(\bP-\bP_0)\bmu}^2/\sigma^2) \) and \( \chi^2(n-r,\norm{(\I-\bP)\bmu}^2/\sigma^2) \).
The second noncentrality is \( 0 \) since \( \bP\bmu=\bmu \). The unknown \( \sigma^2 \) cancels in the
ratio, and @def-qf-noncentral-f applies. Finally, \( (\bP-\bP_0)\bmu=\bmu-\bP_0\bmu \).
:::

In a linear model, \( \bP=\M \) projects onto the column space of the full model matrix \( \X \) and
\( \bP_0=\Mo \) onto that of a submodel \( \X_0 \). The numerator is the drop in residual sum of squares
from the submodel to the full model ([Chapter 6](../ch06-projections/index.html)). The theorem says three things. The
test that rejects for large \( F \) has an exactly known null distribution. Its power depends on
the unknown parameters only through
\( \gamma=\norm{(\M-\Mo)\X\bbeta}^2/\sigma^2 \). And, by @thm-qf-f-power, that power increases
with \( \gamma \).

::: {#exm-qf-curvature}
[Detecting curvature]

Let \( x_1,\dots,x_n \), \( n=15 \), be equally spaced on \( [0,1] \), and suppose
\( Y_i=1+0.5x_i+\beta_2x_i^2+\varepsilon_i \) with independent \( \Normal(0,1) \) errors. To test for
curvature, compare the straight line (\( \X_0=[\bone,\x] \)) with the quadratic
(\( \X=[\bone,\x,\x^2] \)). The \( F \) statistic has \( 1 \) and \( 12 \) degrees of freedom,
and the \( 5\% \) critical value is \( 4.747 \). For \( \beta_2=8.0 \), the
listing finds \( \gamma=6.873 \) and power \( 0.673 \) from the noncentral
\( F \), and the rejection rate in \( 100{,}000 \) simulated data sets is \( 0.674 \).

The noncentrality has a clear meaning. It is \( \beta_2^2 \) times the squared length of the part
of \( \x^2 \) orthogonal to \( \bone \) and \( \x \), divided by \( \sigma^2 \). It grows roughly in proportion to \( n \) for
a fixed design pattern. Keeping the same design on \( [0,1] \), the smallest \( n \) with power at
least \( 0.8 \) is \( 21 \), where \( \gamma=8.97 \).
[Figure 4.4.1](04-independence.html#fig-qf-f-power) places the example on the power curve and shows the loss of power
when the same noncentrality is spread over more numerator degrees of freedom \( r \). At
\( \gamma=10 \) the power is \( 0.827 \), \( 0.606 \) and
\( 0.421 \) for \( r=1,3,6 \).
:::

```{.python .run #cell-f-power-curvature}
import numpy as np
from scipy import stats
rng = np.random.default_rng(1104)
alpha, s = 0.05, 12
def power(gamma, q, s, alpha=0.05):
    crit = stats.f.ppf(1 - alpha, q, s)
    return stats.ncf.sf(crit, q, s, np.maximum(gamma, 1e-12))
gammas = np.linspace(0, 30, 301)
gamma_ref, rs = 10.0, (1, 3, 6)                        # r = numerator degrees of freedom
curves = {q: power(gammas, q, s) for q in rs}

n, sigma, beta0, beta1, beta2, reps = 15, 1.0, 1.0, 0.5, 8.0, 100_000
x = np.linspace(0, 1, n)
X0 = np.column_stack([np.ones(n), x])                  # straight line
X = np.column_stack([X0, x**2])                        # adds curvature

def proj(Z):
    Q, _ = np.linalg.qr(Z)
    return Q @ Q.T

M, M0 = proj(X), proj(X0)
mean = X @ np.array([beta0, beta1, beta2])
gamma = mean @ (M - M0) @ mean / sigma**2              # noncentrality
q, df_resid = 1, n - 3
crit = stats.f.ppf(0.95, q, df_resid)
exact_power = stats.ncf.sf(crit, q, df_resid, gamma)

Y = mean + sigma * rng.standard_normal((reps, n))
num = np.einsum("ij,jk,ik->i", Y, M - M0, Y) / q
den = np.einsum("ij,jk,ik->i", Y, np.eye(n) - M, Y) / df_resid
print(gamma, exact_power, np.mean(num / den > crit))
```

::: {when-format="html"}
![**Figure 4.4.1.** Power \( \Pr\{F(r,12,\gamma)>F_{0.05}(r,12)\} \) of the level-\( 0.05 \) \( F \) test as a
function of the noncentrality. Every curve starts at \( 0.05 \) and increases to \( 1 \)
(@thm-qf-f-power). For a fixed \( \gamma \), power falls as the number \( r \) of numerator
degrees of freedom grows. The dot is the curvature test of @exm-qf-curvature.](f_power.svg){#fig-qf-f-power width=66%}
:::

::: {when-format="pdf"}
![Power \( \Pr\{F(r,12,\gamma)>F_{0.05}(r,12)\} \) of the level-\( 0.05 \) \( F \) test as a
function of the noncentrality. Every curve starts at \( 0.05 \) and increases to \( 1 \)
(@thm-qf-f-power). For a fixed \( \gamma \), power falls as the number \( r \) of numerator
degrees of freedom grows. The dot is the curvature test of @exm-qf-curvature.](f_power.pdf){width=66%}
:::

## Exercises

### A. Check your understanding

::: {#exr-qf-drift-independence}
[A1]

In @exm-qf-drift let \( \Y\sim\Normal_n(\mu\bone,\sigma^2\I) \) with \( n\ge3 \). Use
@thm-qf-mean-var(d) to show that \( \Cov(S^2,Q)=\Var(S^2) \), and deduce that \( S^2 \) and \( Q \) are
dependent with correlation \( \{2(n-1)/(3n-4)\}^{1/2} \). Which linear forms \( \mathbf{a}\T\Y \) are independent of
\( Q \) for every mean vector \( \bmu\in\Real^n \)?
:::

::: {.solution}
With \( \bSigma=\sigma^2\I \), part (d) of
@thm-qf-mean-var gives
\( \Cov(S^2,Q)=2\sigma^4\tr(\mathbf{C}\bD\T\bD)/\{2(n-1)^2\} \). Since \( \bD\bone=\bzero \), \( \mathbf{C}\bD\T\bD=\bD\T\bD \),
whose trace is \( 2(n-1) \), so \( \Cov(S^2,Q)=2\sigma^4/(n-1)=\Var(S^2) \). The mean terms vanish because
\( \mathbf{C}\bmu=\bD\bmu=\bzero \). Using \( \Var(Q) \) from @exm-qf-drift-variance, the correlation is
\( \{\Var(S^2)/\Var(Q)\}^{1/2}=\{2(n-1)/(3n-4)\}^{1/2} \), which is nonzero, so the forms are dependent.
By @thm-qf-indep-linear, \( \mathbf{a}\T\Y \) is independent of \( Q \) for every \( \bmu \) iff
\( \mathbf{a}\T\bD\T\bD=\bzero\T \). Multiplying by \( \mathbf{a} \) gives \( \norm{\bD\mathbf{a}}^2=0 \), so \( \bD\mathbf{a}=\bzero \) and
\( \mathbf{a} \) is a multiple of \( \bone \); conversely such \( \mathbf{a} \) work. The only such linear forms are the
multiples of \( \bar{Y} \).
:::

### B. Practice

::: {#exr-qf-uncorrelated-dependent}
[B1]

Let \( \Y\sim\Normal_2(\bzero,\I) \). Show that \( Y_1^2-Y_2^2 \) and \( 2Y_1Y_2 \) are uncorrelated but not
independent. Why does @prp-qf-craig-nnd not apply?
:::

::: {.solution}
With \( \A=\diag(1,-1) \) and
\( \B=\begin{psmallmatrix}0&1\\1&0\end{psmallmatrix} \), \( \tr(\A\B)=0 \), so the covariance
\( 2\tr(\A\B) \) is zero. In polar coordinates \( Y_1=R\cos\Theta \), \( Y_2=R\sin\Theta \), with
\( R^2\sim\chi^2(2) \) independent of \( \Theta \) uniform on \( [0,2\pi) \), the two forms are
\( R^2\cos2\Theta \) and \( R^2\sin2\Theta \). Now \( \E(R^4)=8 \) and \( \E(R^8)=384 \), so
\( \E[(R^2\cos2\Theta)^2(R^2\sin2\Theta)^2]=384/8=48 \), but the product of the expected squares is
\( (8/2)^2=16 \). The forms are dependent. @prp-qf-craig-nnd needs nonnegative definite
matrices, and \( \A \) and \( \B \) are indefinite.
:::

::: {#exr-qf-singular-independence}
[B2]

Let \( \Y\sim\Normal_n(\bmu,\bSigma) \) and let \( \A,\B \) be nonnegative definite with
\( \bSigma\A\bSigma\B\bSigma=\mathbf{0} \). Show that \( \A\bSigma\B=\mathbf{0} \), so that \( \Y\T\A\Y \) and \( \Y\T\B\Y \) are
independent. Show by an example that the implication can fail when \( \A \) is not nonnegative definite.
:::

::: {.solution}
Write \( \A=\R\R\T \) and \( \B=\bS\bS\T \), for instance
with symmetric square roots, and \( \bSigma=\bL\bL\T \) with \( \bL \) of full column rank. Cancelling \( \bL \) on
both sides of \( \bL(\bL\T\A\bL)(\bL\T\B\bL)\bL\T=\mathbf{0} \) gives \( \bL\T\A\bL\bL\T\B\bL=\mathbf{0} \). The trace
of this matrix is
\[
\tr(\bL\T\R\R\T\bL\bL\T\bS\bS\T\bL)=\tr\{(\R\T\bSigma\bS)(\R\T\bSigma\bS)\T\}=\norm{\R\T\bSigma\bS}_F^2,
\]
so \( \R\T\bSigma\bS=\mathbf{0} \). Hence \( \A\bSigma\B=\R(\R\T\bSigma\bS)\bS\T=\mathbf{0} \), and
@thm-qf-indep-quadratic applies. For nonnegative definite matrices the two conditions are
therefore equivalent, since \( \A\bSigma\B=\mathbf{0} \) trivially implies \( \bSigma\A\bSigma\B\bSigma=\mathbf{0} \).
Without nonnegative definiteness the argument breaks at the trace step. Take \( \bSigma=\diag(1,0) \),
\( \A=\begin{psmallmatrix}0&1\\1&0\end{psmallmatrix} \) and \( \B=\I_2 \). Then \( \bSigma\A\bSigma=\mathbf{0} \), so
\( \bSigma\A\bSigma\B\bSigma=\mathbf{0} \), but \( \A\bSigma\B=\begin{psmallmatrix}0&0\\1&0\end{psmallmatrix}\ne\mathbf{0} \).
With \( \bmu=(0,m)\T \) and \( m\ne0 \) we have \( Y_2=m \) with probability one, so the forms are \( 2mY_1 \) and
\( Y_1^2+m^2 \). The second is a nonconstant function of the first, so they are dependent.
:::

::: {#exr-qf-regression-t}
[B3]

Let \( \Y\sim\Normal_n(\X\bbeta,\sigma^2\I) \) with \( \X \) of full column rank \( p<n \), let
\( \hbeta=(\X\T\X)^{-1}\X\T\Y \), \( \M=\X(\X\T\X)^{-1}\X\T \) and \( s^2=\norm{(\I-\M)\Y}^2/(n-p) \). With \( v_{jj} \) the
\( j \)th diagonal entry of \( (\X\T\X)^{-1} \), show that
\( (\hat{\beta}_j-b)/(s\sqrt{v_{jj}})\sim t\bigl(n-p,(\beta_j-b)/(\sigma\sqrt{v_{jj}})\bigr) \).
:::

::: {.solution}
Here \( \M \) is symmetric and idempotent with \( \tr\M=\tr\{(\X\T\X)^{-1}\X\T\X\}=p \),
and \( \M\X=\X \). So \( (n-p)s^2/\sigma^2=\norm{(\I-\M)\Y}^2/\sigma^2\sim\chi^2(n-p) \) by
@thm-qf-chisq(a), because \( (\I-\M)\X\bbeta=\bzero \). With \( \B=(\X\T\X)^{-1}\X\T \),
\( \B(\sigma^2\I)(\I-\M)=\sigma^2(\X\T\X)^{-1}(\X\T-\X\T\M)=\mathbf{0} \), so \( \hbeta=\B\Y \) is independent of
\( s^2 \) (@thm-qf-indep-linear). Since \( \hbeta\sim\Normal_p(\bbeta,\sigma^2(\X\T\X)^{-1}) \), the
variable \( X=(\hat{\beta}_j-b)/(\sigma\sqrt{v_{jj}}) \) is \( \Normal\bigl((\beta_j-b)/(\sigma\sqrt{v_{jj}}),1\bigr) \).
The statistic is \( X/\sqrt{V/(n-p)} \) with \( V=(n-p)s^2/\sigma^2 \), and @def-qf-noncentral-t applies.
:::

::: {#exr-qf-f-limit}
[B4]

Show that \( F_\alpha(r,s)\to\chi^2_\alpha(r)/r \) as \( s\to\infty \), where \( \chi^2_\alpha(r) \) is the upper \( \alpha \)
point of \( \chi^2(r) \), and that the power \( \Pr\{F(r,s,\gamma)>F_\alpha(r,s)\} \) tends to
\( \Pr\{\chi^2(r,\gamma)>\chi^2_\alpha(r)\} \), the power of the test that knows \( \sigma^2 \).
:::

### C. Going deeper

::: {#exr-qf-shared-denominator}
[C1]

Keep the notation of @thm-qf-nested-f and assume the null hypothesis \( \bP_0\bmu=\bmu \). Put
\( a=r-r_0 \), \( b=n-r \), \( U=\Y\T(\bP-\bP_0)\Y \), \( V=\Y\T(\I-\bP)\Y \), and consider
\[
R=\frac{\Y\T(\bP-\bP_0)\Y}{\Y\T(\I-\bP_0)\Y}=\frac{U}{U+V},
\]
the share of the reduced model's residual sum of squares that the full model removes.

::: {.enumerate options="label=(\alph*)"}
1. Show that the numerator and the denominator of \( R \) are *not* independent, and that their
           correlation is \( \sqrt{a/(a+b)}=\sqrt{(r-r_0)/(n-r_0)} \).

2. Show nevertheless that \( R\sim\mathrm{Beta}(a/2,b/2) \) exactly, that \( R \) is independent of
           \( \Y\T(\I-\bP_0)\Y \), and that \( F=(b/a)R/(1-R) \), so the two tests are the same test.

3. Which hypothesis of @thm-qf-indep-quadratic fails for the pair
           \( \{\Y\T(\bP-\bP_0)\Y,\ \Y\T(\I-\bP_0)\Y\} \), and why does its failure not damage the
           exact null distribution?
:::

:::

::: {.solution}
(a) The matrices \( \bP-\bP_0 \) and \( \I-\bP \) are symmetric idempotent with
zero product (proof of @thm-qf-nested-f), so \( U \) and \( V \) are independent with
\( U/\sigma^2\sim\chi^2(a) \) and \( V/\sigma^2\sim\chi^2(b) \) by @thm-qf-orthogonal-projections; the
noncentralities vanish because \( \bP\bmu=\bmu \) and \( \bP_0\bmu=\bmu \). The denominator is
\( U+V \), so \( \Cov(U,U+V)=\Var(U) \) and
\[
\operatorname{corr}(U,U+V)=\frac{\Var(U)}{\sqrt{\Var(U)\Var(U+V)}}
=\sqrt{\frac{2a\sigma^4}{2(a+b)\sigma^4}}=\sqrt{\frac{a}{a+b}} ,
\]
which is positive: a large \( U \) inflates the denominator too.

(b) This is @exr-qf-beta applied to the independent pair \( (U/\sigma^2,V/\sigma^2) \): the ratio
\( U/(U+V) \) is \( \mathrm{Beta}(a/2,b/2) \) and is independent of \( U+V \), and
\( U/V=R/(1-R) \) gives \( F=(U/a)/(V/b)=(b/a)R/(1-R) \), a strictly increasing function of \( R \).
So rejecting for large \( R \) and rejecting for large \( F \) are the same rule at the same level,
and \( \sigma^2 \) has cancelled in both.

(c) Independence would require \( (\bP-\bP_0)(\I-\bP_0)=\mathbf{0} \), whereas
\( (\bP-\bP_0)(\I-\bP_0)=\bP-\bP_0 \ne\mathbf{0} \). What the theorems supply is a decomposition of
the *denominator* into two independent pieces, one of which is the numerator. Dependence
between numerator and denominator is not a problem in itself; what an exact distribution
needs is that the pair be a function of independent chi-squared variables, which it is.
:::

::: {#exr-qf-gls-independence}
[C2]

Let \( \Y\sim\Normal_n(\X\bbeta,\sigma^2\V) \) with \( \V \) positive definite and \( \X \) of full column
rank \( p<n \). Write \( \mathbf{K}=(\X\T\V^{-1}\X)^{-1} \),
\[
\tilde{\bbeta}=\mathbf{K}\X\T\V^{-1}\Y,\qquad
\A=\V^{-1}-\V^{-1}\X\mathbf{K}\X\T\V^{-1},\qquad Q=\Y\T\A\Y .
\]

::: {.enumerate options="label=(\alph*)"}
1. Show that \( Q=(\Y-\X\tilde{\bbeta})\T\V^{-1}(\Y-\X\tilde{\bbeta}) \), that
           \( \tilde{\bbeta}\sim\Normal_p(\bbeta,\sigma^2\mathbf{K}) \) and \( Q/\sigma^2\sim\chi^2(n-p) \), and
           that \( \tilde{\bbeta} \) and \( Q \) are independent. Deduce an exact \( t \) distribution for
           \( \blambda\T\tilde{\bbeta} \) studentized by \( Q/(n-p) \).

2. Now take the ordinary least squares \( \hbeta=(\X\T\X)^{-1}\X\T\Y \) and
           \( \text{RSS}=\Y\T(\I-\M)\Y \) with \( \M=\X(\X\T\X)^{-1}\X\T \). Show that they are
           independent for every \( \bbeta \) if and only if \( \C(\V\X)\subseteq\C(\X) \).

3. Check that the condition holds for equicorrelated errors, \( \V=(1-\rho)\I+\rho\bone\bone\T \),
           in any model whose columns include \( \bone \), and fails for \( \X=\bone \), \( n=3 \) and
           \( \V=\diag(1,1,4) \).
:::

:::

::: {.solution}
(a) Since \( \Y-\X\tilde{\bbeta}=(\I-\X\mathbf{K}\X\T\V^{-1})\Y \), expanding
\( (\I-\X\mathbf{K}\X\T\V^{-1})\T\V^{-1}(\I-\X\mathbf{K}\X\T\V^{-1}) \) gives
\( \V^{-1}-2\V^{-1}\X\mathbf{K}\X\T\V^{-1}+\V^{-1}\X\mathbf{K}(\X\T\V^{-1}\X)\mathbf{K}\X\T\V^{-1}=\A \). By
@thm-mvn-linear, \( \tilde{\bbeta} \) is normal with mean \( \mathbf{K}\X\T\V^{-1}\X\bbeta=\bbeta \) and
covariance \( \sigma^2\mathbf{K}\X\T\V^{-1}\V\V^{-1}\X\mathbf{K}=\sigma^2\mathbf{K} \). Next, with \( \bSigma=\sigma^2\V \), the
matrix \( (\A/\sigma^2)\bSigma=\A\V \) is idempotent, because
\( \A\V\A=\A \) by the same expansion, and \( \tr(\A\V)=n-\tr(\mathbf{K}\X\T\V^{-1}\X)=n-p \); also
\( \A\X\bbeta=\bzero \) because \( \A\X=\V^{-1}\X-\V^{-1}\X=\mathbf{0} \). So
\( Q/\sigma^2\sim\chi^2(n-p) \) by @thm-qf-chisq(b). For independence take
\( \B=\mathbf{K}\X\T\V^{-1} \) in @thm-qf-indep-linear: \( \B\bSigma\A=\sigma^2\mathbf{K}\X\T\A=\mathbf{0} \), again
because \( \A\X=\mathbf{0} \). Hence
\( \blambda\T(\tilde{\bbeta}-\bbeta)/\sqrt{\blambda\T\mathbf{K}\blambda\,Q/(n-p)}\sim t(n-p) \), exactly as
in @exr-qf-regression-t.

(b) With \( \B=(\X\T\X)^{-1}\X\T \) and \( \A=\I-\M \), @thm-qf-indep-linear says that
independence for every mean vector \( \X\bbeta \) holds iff
\( \B\bSigma\A=\sigma^2(\X\T\X)^{-1}\X\T\V(\I-\M)=\mathbf{0} \), that is, iff \( \X\T\V(\I-\M)=\mathbf{0} \).
Transposing, this says \( (\I-\M)\V\X=\mathbf{0} \), that is, \( \V\X=\M\V\X \), which is
\( \C(\V\X)\subseteq\C(\X) \). (Necessity is for all \( \bbeta \); for one particular \( \bbeta \) the
covariance may vanish by accident.)

(c) For the equicorrelated \( \V \) and a model with \( \bone\in\C(\X) \),
\( \V\X=(1-\rho)\X+\rho\bone(\bone\T\X) \), and each column of \( \bone\bone\T\X \) is a multiple of
\( \bone \), so \( \C(\V\X)\subseteq\C(\X) \): ordinary least squares keeps its exact \( t \) and \( F \)
theory. For \( \X=\bone \) and \( \V=\diag(1,1,4) \), \( \V\bone=(1,1,4)\T \) is not a multiple of
\( \bone \), so \( \bar{Y} \) and the residual sum of squares are dependent, and the usual \( t \)
statistic has no \( t \) distribution. The condition \( \C(\V\X)\subseteq\C(\X) \) is the same one
under which ordinary and generalized least squares coincide
([Chapter 31](../ch31-general-gauss-markov/index.html)).
:::
