# Noncentral chi-squared, t and F distributions

A regression analysis reports a handful of numbers whose sampling behaviour
decides everything: the residual sum of squares, the reduction in it when a
term is added, the ratio of the two. Each is a function of the response vector
of the form
\[
\Y\T\A\Y=\sum_{i=1}^n\sum_{j=1}^n a_{ij}Y_iY_j ,
\]
a **quadratic form** in \( \Y \). The sample variance is one, with
\( \A=(n-1)^{-1}(\I-n^{-1}\bone\bone\T) \). So is the squared length
\( \norm{\bP\Y}^2=\Y\T\bP\Y \) of any symmetric idempotent image of \( \Y \), and this
is how the sums of squares of [Chapter 6](../ch06-projections/index.html) will arise.

Since \( \Y\T\A\Y \) is a scalar, it equals its own transpose \( \Y\T\A\T\Y \), and so
also \( \Y\T\tfrac12(\A+\A\T)\Y \). Replacing \( \A \) by its symmetric part changes
nothing, and *throughout this chapter the matrix of a quadratic form is
symmetric*. Several formulas below, the variance formula among them, are false
for a non-symmetric \( \A \) (@exr-qf-symmetrize).

When \( \Y\sim\Normal_n(\bzero,\I) \) and \( \A \) is symmetric idempotent, \( \Y\T\A\Y \) has a
chi-squared distribution, which is the null distribution of a test. Under an
alternative the mean is not zero, and a slightly larger family appears. We
start with that family.

## Noncentral chi-squared

Recall that \( \chi^2(r) \), the chi-squared distribution with \( r \) degrees of
freedom, is the law of \( Z_1^2+\dots+Z_r^2 \) for independent standard normal
\( Z_i \). It is the gamma distribution with shape \( r/2 \) and scale \( 2 \), with density
\[
g_r(x)=\frac{x^{r/2-1}e^{-x/2}}{2^{r/2}\,\Gamma(r/2)},\qquad x>0,
\]{#eq-qf-chisq-density}

mean \( r \), variance \( 2r \) and moment generating function \( (1-2t)^{-r/2} \) for
\( t<1/2 \). By convention \( \chi^2(0) \) is the point mass at zero. Everything about
the noncentral version follows from one computation.

::: {#lem-qf-shifted-square}
Let \( W\sim\Normal(0,1) \) and \( \nu\in\Real \). For \( t<1/2 \),
\[
\E\bigl[e^{t(W+\nu)^2}\bigr]=(1-2t)^{-1/2}\exp\Bigl(\frac{\nu^2t}{1-2t}\Bigr).
\]
:::

::: {.proof}
With \( u=w+\nu \) the expectation is
\( (2\pi)^{-1/2}\int\exp\{tu^2-(u-\nu)^2/2\}\,du \). Completing the square in the
exponent,
\[
tu^2-\tfrac12(u-\nu)^2
=-\frac{1-2t}{2}\Bigl(u-\frac{\nu}{1-2t}\Bigr)^2+\frac{\nu^2}{2}\Bigl(\frac{1}{1-2t}-1\Bigr).
\]
The last term equals \( \nu^2t/(1-2t) \). The remaining Gaussian integral is
\( (2\pi)^{1/2}(1-2t)^{-1/2} \) when \( 1-2t>0 \).
:::

::: {#def-qf-noncentral-chisq}
[Noncentral chi-squared]

Let \( \Z\sim\Normal_r(\bmu,\I_r) \). The distribution of \( \Z\T\Z=\sum_{i}Z_i^2 \) is the
**noncentral chi-squared distribution** with \( r \) degrees of freedom and
**noncentrality** \( \gamma=\bmu\T\bmu=\norm{\bmu}^2 \), written \( \chi^2(r,\gamma) \).
:::

The definition hides a claim: the law of \( \norm{\Z}^2 \) depends on the mean
vector \( \bmu \) only through its length. @thm-qf-ncchisq proves it. Intuitively,
a rotation \( \Q\Z \) of a spherical normal is again spherical normal, with mean \( \Q\bmu \),
and it has the same length as \( \Z \), so every mean vector of a given length gives
the same answer.

::: {.remark}
[Convention]

In this book the noncentrality is \( \gamma=\bmu\T\bmu \), with no factor \( \tfrac12 \).
Then \( \E\,\chi^2(r,\gamma)=r+\gamma \), which is easy to remember, and \( \gamma \)
agrees with the parameter `nc` of `scipy.stats.ncx2`. Several texts,
including Rencher and Schaalje (2008) and Christensen (2020), use
\( \lambda=\bmu\T\bmu/2 \) instead. When reading them, replace their \( \lambda \) by \( \gamma/2 \).
The same convention applies to the noncentral \( F \) distribution below.
:::

::: {#thm-qf-ncchisq}
[Properties of the noncentral chi-squared]

Let \( U\sim\chi^2(r,\gamma) \).

::: {.enumerate options="label=(\alph*)"}
1. The moment generating function is
     \[
M_U(t)=(1-2t)^{-r/2}\exp\Bigl(\frac{\gamma t}{1-2t}\Bigr),\qquad t<\tfrac12 ,
\]{#eq-qf-ncchisq-mgf}

     so the law depends on \( \bmu \) only through \( \gamma \), and \( \chi^2(r,0)=\chi^2(r) \).

2. \( \E(U)=r+\gamma \) and \( \Var(U)=2r+4\gamma \).

3. If \( U_1,\dots,U_k \) are independent with \( U_i\sim\chi^2(r_i,\gamma_i) \), then
           \( \sum_iU_i\sim\chi^2\bigl(\sum_ir_i,\sum_i\gamma_i\bigr) \).

4. *(Poisson mixture.)* Let \( K\sim\text{Poisson}(\gamma/2) \) and, given
           \( K=k \), let \( V\sim\chi^2(r+2k) \). Then \( V\sim\chi^2(r,\gamma) \). Consequently
           \( U \) has density
           \[
f_{r,\gamma}(x)=\sum_{k=0}^\infty e^{-\gamma/2}\frac{(\gamma/2)^k}{k!}\,g_{r+2k}(x),
          \qquad x>0 .
\]{#eq-qf-ncchisq-density}

:::

:::

::: {.proof}
(a) With \( \Z\sim\Normal_r(\bmu,\I) \), the \( Z_i \) are independent, and
\( Z_i=W_i+\mu_i \) with \( W_i\sim\Normal(0,1) \). By @lem-qf-shifted-square,
\[
\E\,e^{t\Z\T\Z}=\prod_{i=1}^r(1-2t)^{-1/2}\exp\Bigl(\frac{\mu_i^2t}{1-2t}\Bigr)
=(1-2t)^{-r/2}\exp\Bigl(\frac{\gamma t}{1-2t}\Bigr).
\]
A moment generating function that is finite on an open interval around zero
determines the distribution, so the law depends only on \( (r,\gamma) \). Setting
\( \gamma=0 \) gives the central mgf.

(b) The cumulant generating function is \( K(t)=\log M_U(t)=-\tfrac r2\log(1-2t)+\gamma t/(1-2t) \).
Differentiating, \( K'(t)=r(1-2t)^{-1}+\gamma(1-2t)^{-2} \) and
\( K''(t)=2r(1-2t)^{-2}+4\gamma(1-2t)^{-3} \). The mean and variance are \( K'(0) \) and \( K''(0) \).

(c) By independence the mgf of the sum is the product of the mgfs
@eq-qf-ncchisq-mgf, which has the same form with \( \sum r_i \) and \( \sum\gamma_i \).

(d) Conditioning on \( K \),
\[
\E\,e^{tV}=\sum_{k=0}^\infty e^{-\gamma/2}\frac{(\gamma/2)^k}{k!}(1-2t)^{-(r+2k)/2}
=(1-2t)^{-r/2}e^{-\gamma/2}\exp\Bigl(\frac{\gamma/2}{1-2t}\Bigr),
\]
and \( -\tfrac{\gamma}{2}+\tfrac{\gamma}{2(1-2t)}=\tfrac{\gamma t}{1-2t} \). This is @eq-qf-ncchisq-mgf.
The density of a mixture is the mixture of the densities.
:::

The mixture in (d) is more than a formula. It says that noncentrality acts like a
random number of extra degrees of freedom, two at a time, with \( \gamma/2 \) of them
on average. Since each extra pair of degrees of freedom pushes the distribution to
the right, a larger \( \gamma \) should give a stochastically larger variable. That is
the next result, and it is what makes tests based on these distributions sensible.

::: {#exm-qf-ncx-four}
[A noncentral chi-squared with four degrees of freedom]

Take \( r=4 \) and \( \gamma=6 \). The Poisson\( (3) \)
weights of the central densities \( g_{4},g_{6},g_{8},g_{10} \) in
@eq-qf-ncchisq-density are \( 0.050 \), \( 0.149 \),
\( 0.224 \) and \( 0.224 \), and the mean and variance are
\( 10 \) and \( 32 \). The listing checks the mixture
against a library implementation of the density. The script also squares
normal vectors with a mean of length \( \sqrt{6} \) in an arbitrary direction, and their
average squared length is \( 10.002 \). A test that rejects when a
\( \chi^2(4) \) statistic exceeds its upper \( 5\% \) point \( 9.488 \) rejects
with probability \( 0.470 \) when the statistic is in fact
\( \chi^2(4,6) \). [Figure 4.1.1](01-noncentral.html#fig-qf-ncx) shows how the density moves right and spreads out
as \( \gamma \) grows.
:::

```{.python .run #cell-noncentral-chisq-mixture}
import numpy as np
from scipy import stats
rng = np.random.default_rng(20240401)

def ncx2_pdf_mixture(x, r, gamma, terms=200):
    """Density of chi^2(r, gamma) as a Poisson(gamma/2) mixture of central densities."""
    k = np.arange(terms)
    weights = stats.poisson.pmf(k, gamma / 2)            # P(K = k)
    dens = stats.chi2.pdf(np.asarray(x)[:, None], r + 2 * k)
    return dens @ weights

r, gamma = 4, 6.0
x = np.linspace(0.05, 40, 400)
print(np.max(np.abs(ncx2_pdf_mixture(x, r, gamma) - stats.ncx2.pdf(x, r, gamma))))
```

::: {when-format="html"}
![**Figure 4.1.1.** Densities of \( \chi^2(4,\gamma) \). Dotted ticks mark the means
\( 4+\gamma \). The variance \( 2r+4\gamma \) grows with the noncentrality, so the distributions
flatten as they move right.](noncentral_chisq.svg){#fig-qf-ncx width=66%}
:::

::: {when-format="pdf"}
![Densities of \( \chi^2(4,\gamma) \). Dotted ticks mark the means
\( 4+\gamma \). The variance \( 2r+4\gamma \) grows with the noncentrality, so the distributions
flatten as they move right.](noncentral_chisq.pdf){width=66%}
:::

::: {#prp-qf-ncchisq-monotone}
[Monotonicity in the noncentrality]

For fixed \( r\ge1 \) and \( x>0 \), the tail probability \( \gamma\mapsto\Pr\{\chi^2(r,\gamma)>x\} \)
is continuous and strictly increasing on \( [0,\infty) \), and tends to \( 1 \) as \( \gamma\to\infty \).
:::

::: {.proof}
Let \( h(k)=\Pr\{\chi^2(r+2k)>x\} \) and \( p_k(\lambda)=e^{-\lambda}\lambda^k/k! \). By
@thm-qf-ncchisq(d), the tail probability is \( \phi(\lambda)=\sum_kp_k(\lambda)h(k) \)
with \( \lambda=\gamma/2 \). First, \( h \) is strictly increasing: if \( V\sim\chi^2(r+2k) \) and
\( E\sim\chi^2(2) \) are independent, then \( V+E\sim\chi^2(r+2k+2) \) and
\( h(k+1)-h(k)=\Pr(V\le x<V+E)>0 \), since \( (V,E) \) has a positive joint density on
\( (0,\infty)^2 \). Next, \( p_k'(\lambda)=p_{k-1}(\lambda)-p_k(\lambda) \) (with \( p_{-1}=0 \)).
Because \( 0\le h\le1 \), the series of derivatives converges uniformly on bounded
intervals, so term-by-term differentiation is valid and
\[
\phi'(\lambda)=\sum_{k\ge0}\bigl(p_{k-1}(\lambda)-p_k(\lambda)\bigr)h(k)
=\sum_{k\ge0}p_k(\lambda)\bigl(h(k+1)-h(k)\bigr)>0 .
\]
For the limit, fix \( \varepsilon>0 \). Since \( \chi^2(r+2m) \) has mean \( r+2m \) and
standard deviation \( \sqrt{2r+4m} \), Chebyshev's inequality gives \( h(m)>1-\varepsilon \)
for some \( m \). Then \( \phi(\lambda)\ge h(m)\Pr(K\ge m) \) for \( K\sim\text{Poisson}(\lambda) \),
and \( \Pr(K\ge m)\to1 \) as \( \lambda\to\infty \).
:::

## Noncentral \( F \) and \( t \)

Test statistics divide one quadratic form by another, to remove the unknown
scale. The two ratio distributions that result are defined next.

::: {#def-qf-noncentral-f}
[Noncentral \( F \)]

Let \( U\sim\chi^2(r,\gamma) \) and \( V\sim\chi^2(s) \) be independent, \( r,s\ge1 \). The
law of
\[
F=\frac{U/r}{V/s}
\]
is the **noncentral \( F \) distribution** \( F(r,s,\gamma) \), with numerator and
denominator degrees of freedom \( r \) and \( s \) and noncentrality \( \gamma \). The central
\( F(r,s)=F(r,s,0) \) has upper \( \alpha \) point \( F_\alpha(r,s) \).
:::

::: {#def-qf-noncentral-t}
[Noncentral \( t \)]

Let \( X\sim\Normal(\delta,1) \) and \( V\sim\chi^2(s) \) be independent. The law of
\( T=X/\sqrt{V/s} \) is the **noncentral \( t \) distribution** \( t(s,\delta) \), and
\( t(s)=t(s,0) \).
:::

Only the numerator is allowed to be noncentral. A noncentral denominator gives a
*doubly noncentral* \( F \), which appears when the model used for the error sum of
squares is itself wrong (@exr-qf-sse-bias). Two relations are immediate.
First, \( X^2\sim\chi^2(1,\delta^2) \), so \( T^2\sim F(1,s,\delta^2) \): a two-sided \( t \)
test is an \( F \) test with one numerator degree of freedom. Second, since
\( \E(V^{-1})=1/(s-2) \) for \( s>2 \), independence gives
\[
\E(F)=\frac{r+\gamma}{r}\cdot\frac{s}{s-2},\qquad s>2 .
\]

The central \( F \) is the null distribution of the tests of
[Chapter 11](../ch11-general-linear-hypothesis/index.html), and \( F(r,s,\gamma) \) is their distribution
under an alternative. The power of such a test is therefore
\( \Pr\{F(r,s,\gamma)>F_\alpha(r,s)\} \), and the following fact is what guarantees
that it behaves as a test should.

::: {#thm-qf-f-power}
[Power of the \( F \) test increases with \( \gamma \)]

For fixed \( r,s\ge1 \) and \( c>0 \), \( \beta(\gamma)=\Pr\{F(r,s,\gamma)>c\} \) is continuous
and strictly increasing in \( \gamma\ge0 \), and \( \beta(\gamma)\to1 \) as
\( \gamma\to\infty \). In particular, the level-\( \alpha \) test that rejects when
\( F>F_\alpha(r,s) \) has power greater than \( \alpha \) at every \( \gamma>0 \).
:::

::: {.proof}
With \( U,V \) as in @def-qf-noncentral-f, \( F>c \) iff \( U>crV/s \). Condition on \( V \),
which is independent of \( U \):
\[
\beta(\gamma)=\E\bigl[\psi_\gamma(crV/s)\bigr],\qquad
\psi_\gamma(x)=\Pr\{\chi^2(r,\gamma)>x\}.
\]
Since \( V>0 \) with probability one, @prp-qf-ncchisq-monotone gives, for
\( \gamma<\gamma' \), \( \psi_\gamma(crV/s)<\psi_{\gamma'}(crV/s) \) almost surely, so
\( \beta(\gamma)<\beta(\gamma') \). Continuity and the limit follow from the same
proposition by dominated convergence, since \( 0\le\psi_\gamma\le1 \).
:::

Two further monotonicity facts hold for fixed \( \gamma>0 \) and fixed level. The power
decreases as the numerator degrees of freedom \( r \) grow, and increases as the
denominator degrees of freedom \( s \) grow (Ghosh 1973). The first says
that spreading a fixed amount of signal over more dimensions makes it harder to
detect. The second says that estimating \( \sigma^2 \) better helps. We use both only
qualitatively, and [Figure 4.4.1](04-independence.html#fig-qf-f-power) in [Section 4.4](04-independence.html) shows
the first. Because \( T^2\sim F(1,s,\delta^2) \), @thm-qf-f-power also shows that
the power of the two-sided \( t \) test increases with \( \lvert\delta\rvert \)
(@exr-qf-t-square).

## Exercises

### A. Check your understanding

::: {#exr-qf-chisq-two}
[A1]

Show that \( \chi^2(2) \) is the exponential distribution with mean \( 2 \). Use
@thm-qf-ncchisq(d) to show that
\[
\Pr\{\chi^2(2,\gamma)>x\}=\sum_{k\ge0}e^{-\gamma/2}\frac{(\gamma/2)^k}{k!}\;e^{-x/2}\sum_{j=0}^k\frac{(x/2)^j}{j!} .
\]
:::

::: {#exr-qf-t-square}
[A2]

Show that \( T\sim t(s,\delta) \) implies \( T^2\sim F(1,s,\delta^2) \). Deduce that the power of the
two-sided level-\( \alpha \) \( t \) test depends on \( \delta \) only through \( \lvert\delta\rvert \) and increases
with it. Show directly, by conditioning on the denominator, that the power of the one-sided test
that rejects for \( T>c \) increases with \( \delta \).
:::

::: {.solution}
\( T^2=X^2/(V/s) \) with \( X^2\sim\chi^2(1,\delta^2) \) independent of
\( V \), which is @def-qf-noncentral-f. The two-sided test rejects when
\( \lvert T\rvert>t_{\alpha/2}(s) \), that is, when \( T^2>F_\alpha(1,s) \). Its power is therefore
\( \Pr\{F(1,s,\delta^2)>F_\alpha(1,s)\} \), a strictly increasing function of \( \delta^2 \) by @thm-qf-f-power.
For the one-sided test, \( \Pr(T>c)=\E\,\Pr\{X>c\sqrt{V/s}\mid V\}
=\E\,\Phi(\delta-c\sqrt{V/s}) \), and \( \Phi \) is strictly increasing.
:::

### B. Practice

::: {#exr-qf-cumulants}
[B1]

Show that the \( k \)th cumulant of \( \chi^2(r,\gamma) \) is \( 2^{k-1}(k-1)!\,(r+k\gamma) \). Deduce that its
skewness tends to zero like \( 3/\sqrt{\gamma} \) as \( \gamma\to\infty \) with \( r \) fixed.
:::

::: {.solution}
Using \( -\tfrac12\log(1-2t)=\sum_{k\ge1}(2t)^k/(2k) \) and
\( t/(1-2t)=\sum_{k\ge1}2^{k-1}t^k \), the cumulant generating function is
\[
K(t)=\sum_{k\ge1}t^k\,2^{k-1}\Bigl(\frac rk+\gamma\Bigr)=\sum_{k\ge1}\frac{t^k}{k!}\,2^{k-1}(k-1)!\,(r+k\gamma).
\]
The \( k \)th cumulant is \( k! \) times the coefficient of \( t^k \). In particular \( \kappa_2=2(r+2\gamma) \)
and \( \kappa_3=8(r+3\gamma) \), so the skewness is
\( \kappa_3/\kappa_2^{3/2}=8(r+3\gamma)/\{2(r+2\gamma)\}^{3/2} \). As \( \gamma\to\infty \) this behaves like
\( 24\gamma/(4\gamma)^{3/2}=3/\sqrt\gamma \).
:::

::: {#exr-qf-one-df-density}
[B2]

Show that \( \chi^2(1,\gamma) \) has density \( (2\pi x)^{-1/2}e^{-(x+\gamma)/2}\cosh\sqrt{\gamma x} \) for
\( x>0 \). Verify that expanding \( \cosh \) in a power series reproduces @eq-qf-ncchisq-density.
*Hint:* \( \Gamma(k+\tfrac12)=\sqrt\pi\,(2k)!/(4^kk!) \).
:::

::: {.solution}
Let \( X=(W+\sqrt\gamma)^2 \) with \( W\sim\Normal(0,1) \). For
\( x>0 \), \( X\le x \) iff \( -\sqrt x-\sqrt\gamma\le W\le\sqrt x-\sqrt\gamma \). Differentiating with \( \phi \)
the standard normal density,
\[
f(x)=\frac{\phi(\sqrt x-\sqrt\gamma)+\phi(\sqrt x+\sqrt\gamma)}{2\sqrt x}
=\frac{e^{-(x+\gamma)/2}}{\sqrt{2\pi x}}\cdot\frac{e^{\sqrt{\gamma x}}+e^{-\sqrt{\gamma x}}}{2}.
\]
Expanding, \( f(x)=(2\pi x)^{-1/2}e^{-(x+\gamma)/2}\sum_k(\gamma x)^k/(2k)! \). The \( k \)th term of
@eq-qf-ncchisq-density is
\[
\frac{e^{-\gamma/2}(\gamma/2)^k\,x^{k-1/2}e^{-x/2}}{k!\,2^{k+1/2}\,\Gamma(k+\frac12)} .
\]
By the hint the denominator is \( k!\,2^{k+1/2}\sqrt\pi\,(2k)!/(4^kk!)=\sqrt{2\pi}\,2^{-k}(2k)! \), and the
\( k \)th term becomes \( (2\pi x)^{-1/2}e^{-(x+\gamma)/2}(\gamma x)^k/(2k)! \), as required.
:::

::: {#exr-qf-normal-limit}
[B3]

Let \( U_\gamma\sim\chi^2(r,\gamma) \) with \( r \) fixed. Show that
\( (U_\gamma-r-\gamma)/\sqrt{2r+4\gamma} \) converges in distribution to \( \Normal(0,1) \) as \( \gamma\to\infty \).
*Hint:* \( U_\gamma \) has the law of \( (W_1+\sqrt\gamma)^2+W_2^2+\dots+W_r^2 \) with
\( W_1,\dots,W_r \) independent standard normal.
:::

::: {#exr-qf-beta}
[B4]

Let \( U\sim\chi^2(a) \) and \( V\sim\chi^2(b) \) be independent. Show that \( B=U/(U+V) \) has the
\( \mathrm{Beta}(a/2,b/2) \) distribution and is independent of \( U+V \), and that
\( F=(U/a)/(V/b)=(b/a)\,B/(1-B) \). Deduce the null distribution of \( R^2 \) asked for in
[Chapter 6](../ch06-projections/index.html).
:::

::: {.solution}
The joint density of \( (U,V) \) is proportional to
\( u^{a/2-1}v^{b/2-1}e^{-(u+v)/2} \). Change variables to \( w=u/(u+v)\in(0,1) \) and \( s=u+v>0 \), so that
\( u=ws \), \( v=(1-w)s \) and the Jacobian is \( s \). The density of \( (W,S) \) is proportional to
\[
w^{a/2-1}(1-w)^{b/2-1}\cdot s^{(a+b)/2-1}e^{-s/2},
\]
a product of a function of \( w \) and a function of \( s \). So \( B \) and \( U+V \) are independent,
\( B\sim\mathrm{Beta}(a/2,b/2) \) and \( U+V\sim\chi^2(a+b) \). Also \( U/V=B/(1-B) \). For \( R^2 \) in a model
with intercept, \( R^2=\text{SSR}/(\text{SSR}+\text{SSE}) \), where \( \text{SSR}=\norm{(\M-\bP_1)\Y}^2 \)
and \( \text{SSE}=\norm{(\I-\M)\Y}^2 \) come from orthogonal projections of ranks \( r-1 \) and \( n-r \). If
\( \bmu=\mu\bone \), both noncentralities vanish, and @thm-qf-orthogonal-projections with
\( a=r-1 \), \( b=n-r \) gives \( R^2\sim\mathrm{Beta}\bigl((r-1)/2,(n-r)/2\bigr) \).
:::

::: {#exr-qf-mlr}
[B5]

Show that \( f_{r,\gamma}(x)/g_r(x) \) is strictly increasing in \( x>0 \) for every \( \gamma>0 \). Conclude
that among tests based on one observation \( U\sim\chi^2(r,\gamma) \) of \( \gamma=0 \) against a fixed
\( \gamma>0 \), rejecting for large \( U \) is most powerful.
:::

### C. Going deeper

::: {#exr-qf-weighted-sum}
[C1]

Let \( W_1,W_2 \) be independent \( \Normal(0,1) \) and \( U=W_1^2+\sigma^2W_2^2 \) with \( \sigma^2>0 \). A
weighted sum of squares like this is what a general quadratic form in a normal vector turns
out to be ([Section 4.2](02-moments.html)), so it is worth knowing that it is almost never
noncentral chi-squared.

::: {.enumerate options="label=(\alph*)"}
1. Show that matching the mean and variance of \( U \) to those of \( \chi^2(r,\gamma) \) forces
           \( \gamma=\sigma^4-\sigma^2 \) and \( r=1+2\sigma^2-\sigma^4 \). Deduce that no \( \chi^2(r,\gamma) \)
           has the first two moments of \( U \) when \( \sigma^2<1 \).

2. Show more: if \( cU\sim\chi^2(r,\gamma) \) for some \( c>0 \), \( r\ge1 \) and \( \gamma\ge0 \), then
           \( c=\sigma^2=1 \). *Hint:* by @exr-qf-cumulants the \( k \)th cumulant of \( \chi^2(r,\gamma) \) is
           \( 2^{k-1}(k-1)!\,(r+k\gamma) \), and that of \( cU \) is \( 2^{k-1}(k-1)!\,(c^k+(c\sigma^2)^k) \).
           Use \( k=1,2,3 \).

3. Check that \( \sigma^2=1 \) does give a member of the family, and say which one.
:::

:::

::: {.solution}
(a) By @thm-qf-ncchisq(b), matching gives \( r+\gamma=\E U=1+\sigma^2 \) and
\( 2r+4\gamma=\Var U=2(1+\sigma^4) \), since \( \Var(W_1^2)=2 \) and \( \Var(\sigma^2W_2^2)=2\sigma^4 \). The
second equation is \( r+2\gamma=1+\sigma^4 \); subtracting the first gives
\( \gamma=\sigma^4-\sigma^2=\sigma^2(\sigma^2-1) \) and then \( r=1+2\sigma^2-\sigma^4 \). A noncentrality must
be nonnegative, so \( \sigma^2<1 \) is impossible.

(b) Write \( e=c \) and \( d=c\sigma^2 \), both positive. The cumulants of \( cU \) are
\( 2^{k-1}(k-1)!(e^k+d^k) \), because cumulants of a sum of independent variables add and those
of \( aW^2 \) are \( a^k \) times those of \( \chi^2(1) \). Equating with \( 2^{k-1}(k-1)!(r+k\gamma) \) for
\( k=1,2,3 \) gives
\[
e+d=r+\gamma,\qquad e^2+d^2=r+2\gamma,\qquad e^3+d^3=r+3\gamma .
\]
Each consecutive difference equals \( \gamma \), so \( e^2-e+d^2-d=e^3-e^2+d^3-d^2 \), that is
\( e(e-1)^2+d(d-1)^2=0 \). Both terms are nonnegative, so \( e=d=1 \), which is \( c=1 \) and
\( \sigma^2=1 \).

(c) With \( \sigma^2=1 \), \( U=W_1^2+W_2^2\sim\chi^2(2) \), and the three equations are indeed
satisfied by \( c=1 \), \( r=2 \), \( \gamma=0 \). The conclusion is that the noncentral chi-squared
family is not closed under weighting: only equal weights give a member of it, which is why
[Section 4.3](03-chisq.html) will ask for an idempotent matrix.
:::

::: {#exr-qf-local-power}
[C2]

@prp-qf-ncchisq-monotone says the power of a chi-squared test increases with \( \gamma \). How
fast?

::: {.enumerate options="label=(\alph*)"}
1. For \( z>0 \), show that
           \( \Pr\{\chi^2(1,\gamma)>z^2\}=\Phi(-z+\sqrt\gamma)+\Phi(-z-\sqrt\gamma) \), where \( \Phi \) is
           the standard normal distribution function. This is the power of the two-sided test
           based on a single \( \Normal(\sqrt{\gamma},1) \) observation.

2. Show that the right side equals \( \alpha+z\varphi(z)\,\gamma+o(\gamma) \) as \( \gamma\downarrow0 \),
           where \( \alpha=2\Phi(-z) \) and \( \varphi \) is the standard normal density. The power rises
           *linearly* in \( \gamma \), with no first-order term in \( \sqrt\gamma \).

3. For general \( r \), show from the Poisson mixture in the proof of @prp-qf-ncchisq-monotone
           that
           \[
\frac{d}{d\gamma}\Pr\{\chi^2(r,\gamma)>x\}\Big|_{\gamma=0}
          =\tfrac12\bigl[\Pr\{\chi^2(r+2)>x\}-\Pr\{\chi^2(r)>x\}\bigr].
\]
           Check that this agrees with (b) when \( r=1 \) and \( x=z^2 \), given that
           \( \Pr\{\chi^2(3)>x\}=\Pr\{\chi^2(1)>x\}+2\sqrt x\,\varphi(\sqrt x) \).
:::

:::

::: {.solution}
(a) By @def-qf-noncentral-chisq, \( X^2\sim\chi^2(1,\gamma) \) when
\( X\sim\Normal(\sqrt\gamma,1) \). Then \( \Pr(X^2>z^2)=\Pr(X>z)+\Pr(X<-z) \), and
\( \Pr(X>z)=1-\Phi(z-\sqrt\gamma)=\Phi(-z+\sqrt\gamma) \), \( \Pr(X<-z)=\Phi(-z-\sqrt\gamma) \).

(b) Put \( \delta=\sqrt\gamma \) and \( \beta(\delta)=\Phi(-z+\delta)+\Phi(-z-\delta) \). Then
\( \beta(0)=2\Phi(-z)=\alpha \), \( \beta'(\delta)=\varphi(-z+\delta)-\varphi(-z-\delta) \), which
vanishes at \( \delta=0 \), and \( \beta''(\delta)=\varphi'(-z+\delta)+\varphi'(-z-\delta) \). Since
\( \varphi'(u)=-u\varphi(u) \), \( \beta''(0)=2z\varphi(z) \). Taylor's theorem gives
\( \beta(\delta)=\alpha+z\varphi(z)\delta^2+o(\delta^2) \), which is the claim because
\( \delta^2=\gamma \). The even symmetry in \( \delta \) is what kills the \( \sqrt\gamma \) term: a
two-sided test cannot tell the sign of the shift.

(c) In the notation of that proof, with its \( \phi \) written \( \Psi \) here to keep \( \varphi \)
for the normal density, the tail probability is \( \Psi(\lambda)=\sum_kp_k(\lambda)h(k) \) with
\( \lambda=\gamma/2 \) and \( h(k)=\Pr\{\chi^2(r+2k)>x\} \), and
\( \Psi'(\lambda)=\sum_kp_k(\lambda)\{h(k+1)-h(k)\} \). At \( \lambda=0 \) only \( p_0=1 \) survives, so
\( \Psi'(0)=h(1)-h(0) \), and \( d\gamma=2\,d\lambda \) gives the stated derivative. For \( r=1 \) and
\( x=z^2 \) the bracket is \( 2\sqrt x\varphi(\sqrt x)=2z\varphi(z) \), and half of it is
\( z\varphi(z) \), matching (b). The same calculation for a general \( r \) shows how much a design
must raise \( \gamma \) to buy a given amount of power near the null.
:::
