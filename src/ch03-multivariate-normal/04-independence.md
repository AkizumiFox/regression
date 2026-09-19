# Independence and zero covariance

For any random vectors with finite second moments, independence implies zero
cross-covariance, because the expectation of a product of functions of independent
variables factors. The converse fails in general, as @exm-mvn-sign-flip and @exm-mvn-mixture
show. For jointly normal vectors it holds. This is the single most used property of
the normal distribution in linear-model theory. Every independence statement
behind a \( t \) or \( F \) statistic reduces to it.

::: {#thm-mvn-independence}
[Independence]

Let \( \Y\sim\Normal_n(\bmu,\bSigma) \) be partitioned into blocks
\( \Y=(\Y_1\T,\dots,\Y_k\T)\T \), with \( \bmu \) and \( \bSigma=(\bSigma_{ij}) \) partitioned to
match.

::: {.enumerate options="label=(\alph*)"}
1. The blocks \( \Y_1,\dots,\Y_k \) are mutually independent iff
           \( \bSigma_{ij}=\Cov(\Y_i,\Y_j)=\bzero \) for all \( i\ne j \).

2. In particular, \( \Y_1 \) and \( \Y_2 \) are independent iff \( \Cov(\Y_1,\Y_2)=\bzero \),
           and the scalar entries \( Y_1,\dots,Y_n \) are mutually independent iff
           \( \bSigma \) is diagonal.
:::

:::

::: {.proof}
If the blocks are independent, every cross-covariance is zero. Conversely, suppose
\( \bSigma_{ij}=\bzero \) for \( i\ne j \), and partition \( \mathbf{t} \) to match. Then
\[
\mathbf{t}\T\bmu+\tfrac12\mathbf{t}\T\bSigma\mathbf{t}
=\sum_{i=1}^k\Bigl(\mathbf{t}_i\T\bmu_i+\tfrac12\mathbf{t}_i\T\bSigma_{ii}\mathbf{t}_i\Bigr),
\]
because all the cross terms \( \mathbf{t}_i\T\bSigma_{ij}\mathbf{t}_j \) vanish. Exponentiating,
@thm-mvn-mgf gives \( M_{\Y}(\mathbf{t})=\prod_iM_i(\mathbf{t}_i) \), where \( M_i \) is the
moment generating function of \( \Normal(\bmu_i,\bSigma_{ii}) \). By
@cor-mvn-marginals that is the distribution of \( \Y_i \), so \( M_i=M_{\Y_i} \).
The factorization property of [Section 3.2](02-mgf-density.html) gives mutual independence. Part (b) is the cases
\( k=2 \) and \( k=n \).
:::

A consequence worth stating explicitly: for jointly normal blocks, *pairwise*
independence implies *mutual* independence, because condition (a) involves
only pairs. For general random vectors this is false, as the example of three random
signs in @exm-rv-pairwise shows. Joint normality of the whole vector matters here:
joint normality of each pair is not enough (@exr-mvn-pairwise-normal).

The theorem is usually applied to two linear functions of the same normal vector.

::: {#cor-mvn-AY-BY}
Let \( \Y\sim\Normal_n(\bmu,\bSigma) \), and let \( \A \) and \( \B \) be matrices with \( n \)
columns. Then \( \A\Y \) and \( \B\Y \) are independent iff \( \A\bSigma\B\T=\bzero \). In that
case \( g(\A\Y) \) and \( h(\B\Y) \) are independent for all measurable functions \( g \) and \( h \).
:::

::: {.proof}
The stacked vector \( \bigl((\A\Y)\T,(\B\Y)\T\bigr)\T \) is a linear function of \( \Y \),
so it is normal by @thm-mvn-linear. Its off-diagonal covariance block is
\( \Cov(\A\Y,\B\Y)=\A\bSigma\B\T \) (@thm-rv-linear). Apply @thm-mvn-independence. Functions of independent
random vectors are independent.
:::

The condition \( \A\bSigma\B\T=\bzero \) involves only second moments, but the
conclusion is full independence. In particular, any quadratic form in \( \B\Y \), such
as \( \norm{\B\Y}^2 \), is independent of \( \A\Y \). [Chapter 4](../ch04-quadratic-forms/index.html) extends this to
quadratic forms that are not written as squared lengths (@thm-qf-indep-linear).

::: {#exm-mvn-mean-deviations}
[Sample mean and deviations]

Let \( \Y\sim\Normal_n(\mu\bone,\sigma^2\I) \), so that \( Y_1,\dots,Y_n \) is a random sample
from \( \Normal(\mu,\sigma^2) \). The mean is \( \bar{Y}=n^{-1}\bone\T\Y \), and the vector of
deviations is \( \mathbf{d}=(\I-n^{-1}\bone\bone\T)\Y \). As in @exm-rv-centering,
\[
\Cov(\bar{Y},\mathbf{d})=\sigma^2n^{-1}\bone\T(\I-n^{-1}\bone\bone\T)=\bzero\T ,
\]
By @cor-mvn-AY-BY, \( \bar{Y} \) is independent of \( \mathbf{d} \), and so of every function
of the deviations: the sample variance \( S^2=\norm{\mathbf{d}}^2/(n-1) \), the range, the
sample skewness. The distribution of \( (n-1)S^2/\sigma^2 \) is found in @cor-qf-sample-variance.
The computation used \( \Cov(\Y)=\sigma^2\I \) only through
\( \bSigma\bone\in\spn(\bone) \), which is the precise condition for the conclusion
(@exr-mvn-equicorrelated-mean).
:::

::: {#exm-mvn-ls-independence}
[Coefficients and residuals]

In the setting of @exm-mvn-ls-estimator,
\[
\Cov\bigl(\hbeta,\he\bigr)=(\X\T\X)^{-1}\X\T(\sigma^2\I)(\I-\M)
=\sigma^2(\X\T\X)^{-1}\bigl\{(\I-\M)\X\bigr\}\T=\bzero ,
\]
since \( (\I-\M)\X=\bzero \). Hence \( \hbeta \) is independent of the residual vector, and so of
the residual sum of squares \( \norm{\he}^2 \) and of the estimate \( s^2=\norm{\he}^2/(n-p) \).
The same argument shows that the fitted values \( \hY \) and the residuals are independent.
This independence is what allows a \( t \) statistic
\( (\blambda\T\hbeta-\blambda\T\bbeta)/\{s^2\blambda\T(\X\T\X)^{-1}\blambda\}^{1/2} \)
to have a \( t \) distribution (@thm-glh-t-test).
:::

::: {.remark}
[What joint normality buys]

The two examples use only a covariance calculation. Without joint normality, zero
covariance between \( \hbeta \) and \( \he \) still holds under the second-moment
assumptions of @def-lm-linear-model, but independence does not follow,
and neither do exact \( t \) and \( F \) distributions. Chapter 21
examines what survives approximately.
:::

## Exercises

### A. Check your understanding

::: {#exr-mvn-sum-difference}
[A1]

Let \( Y_1,Y_2 \) be independent \( \Normal(\mu,\sigma^2) \). Show that \( Y_1+Y_2 \) and \( Y_1-Y_2 \) are
independent. What if the variances differ?
:::

::: {#exr-mvn-projection-independence}
[A2]

Let \( \Y\sim\Normal_n(\bzero,\I) \) and let \( \bu \) be a unit vector. Show that \( \bu\T\Y \) and
\( \norm{\Y-(\bu\T\Y)\bu}^2 \) are independent.
:::

### B. Practice

::: {#exr-mvn-equicorrelated-mean}
[B1]

Let \( \Y\sim\Normal_n(\mu\bone,\bSigma) \) with \( \bSigma \) arbitrary. Show that \( \bar{Y} \) is
independent of the deviations \( Y_i-\bar{Y} \) iff \( \bSigma\bone\in\spn(\bone) \), that is, iff all
row sums of \( \bSigma \) are equal. Check the condition for the equicorrelation matrix
\( \sigma^2\{(1-\rho)\I+\rho\bone\bone\T\} \), and for a two-group matrix: the first \( m_1 \) and the
last \( m_2 \) observations form groups, all variances are \( 1 \), the correlation is \( \rho \) within a
group and \( \tau \) between groups.
:::

::: {.solution}
\( \bar{Y}=n^{-1}\bone\T\Y \) and the deviations are
\( (\I-n^{-1}\bone\bone\T)\Y \), so by @cor-mvn-AY-BY they are independent iff
\( n^{-1}\bone\T\bSigma(\I-n^{-1}\bone\bone\T)=\bzero\T \), that is,
\( (\I-n^{-1}\bone\bone\T)\bSigma\bone=\bzero \). The null space of \( \I-n^{-1}\bone\bone\T \) is
\( \spn(\bone) \), and \( \bSigma\bone \) is the vector of row sums. For the equicorrelation matrix,
\( \bSigma\bone=\sigma^2(1-\rho+n\rho)\bone \), so independence holds for every admissible \( \rho \).
In the two-group matrix, a row in the first group sums to \( 1+(m_1-1)\rho+m_2\tau \) and a row in
the second to \( 1+(m_2-1)\rho+m_1\tau \). The difference is \( (m_1-m_2)(\rho-\tau) \), so independence
holds iff the groups have equal size or \( \rho=\tau \), and the latter is the equicorrelated case.
:::

::: {#exr-mvn-rank-deficient-independence}
[B2]

Let \( \Y\sim\Normal_n(\X\bbeta,\sigma^2\I) \) with \( \rank(\X)=r<p \), let \( \G \) be a generalized inverse of
\( \X\T\X \) and \( \M=\X\G\X\T \). Show that the least squares estimate \( \G\X\T\Y \) and the residual
vector \( (\I-\M)\Y \) are independent, and that every estimable function \( \blambda\T\bbeta \) with
\( \blambda\in\C(\X\T) \) has an estimate independent of the residual sum of squares.
:::

### C. Going deeper

::: {#exr-mvn-pairwise-normal}
[C1]

Let \( \varphi \) be the \( \Normal(0,1) \) density, \( 0<\varepsilon\le e^{3/2} \), and
\[
f(x,y,z)=\varphi(x)\varphi(y)\varphi(z)\bigl\{1+\varepsilon\,xyz\,e^{-(x^2+y^2+z^2)/2}\bigr\}.
\]
Show that \( f \) is a probability density, and that under \( f \) every pair of coordinates is
\( \Normal_2(\bzero,\I) \). So the three variables are pairwise jointly normal and pairwise
independent. Show that they are not mutually independent and that \( (X,Y,Z) \) is not jointly
normal. Which step of @thm-mvn-independence fails?
:::

::: {.solution}
Since \( \lvert t\rvert e^{-t^2/2}\le e^{-1/2} \) for all
\( t \), the perturbation satisfies \( \lvert\varepsilon\,xyz\,e^{-(x^2+y^2+z^2)/2}\rvert\le\varepsilon e^{-3/2}\le1 \),
so \( f\ge0 \). The function \( z\mapsto z\varphi(z)e^{-z^2/2} \) is odd and integrable, so integrating
\( f \) over \( z \) removes the perturbation and leaves \( \varphi(x)\varphi(y) \). Integrating further
shows that \( f \) integrates to \( 1 \), and by symmetry each pair has density
\( \varphi\otimes\varphi \), the \( \Normal_2(\bzero,\I) \) density. Mutual independence would force
\( f=\varphi(x)\varphi(y)\varphi(z) \) almost everywhere, but the perturbation is nonzero off the
coordinate planes. Directly,
\( \E(XYZ)=\varepsilon\bigl(\int t^2\varphi(t)e^{-t^2/2}\,dt\bigr)^3>0 \), whereas independence with
mean-zero margins gives \( 0 \). If \( (X,Y,Z) \) were jointly normal, its covariance matrix would be \( \I \)
(unit variances from the margins, zero covariances from the pairs), and by
@thm-mvn-independence the coordinates would be independent, which they are not. The
theorem's hypothesis is joint normality of the whole vector. Here only the pairs are normal,
and the moment generating function of the triple does not factor.
:::
