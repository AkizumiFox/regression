# Moment generating functions and independence

Means and covariance matrices summarize a distribution. They do not determine it, and they
cannot certify independence. For that we need a transform of the whole distribution. This
section collects the facts that [Chapter 3](../ch03-multivariate-normal/index.html) and [Chapter 4](../ch04-quadratic-forms/index.html) use.

::: {#def-rv-mgf}
[Moment generating function]

The **moment generating function** of a \( p\times1 \) random vector \( \Y \) is
\[
M_{\Y}(\mathbf{t})=\E\exp(\mathbf{t}\T\Y),\qquad \mathbf{t}\in\Real^p ,
\]
where the expectation may be \( +\infty \). We say that \( M_{\Y} \) *exists* if it is finite for
all \( \mathbf{t} \) in some neighbourhood \( \{\norm{\mathbf{t}}<h\} \) of the origin, \( h>0 \).
:::

::: {#thm-rv-mgf}
[Properties of moment generating functions]

Suppose \( M_{\Y} \) exists.

::: {.enumerate options="label=(\alph*)"}
1. *Uniqueness.* If \( M_{\mathbf{V}} \) also exists and \( M_{\Y}=M_{\mathbf{V}} \) on a neighbourhood of
           \( \bzero \), then \( \Y \) and \( \mathbf{V} \) have the same distribution.

2. *Moments.* All moments of \( \Y \) are finite. The gradient of \( M_{\Y} \) at \( \bzero \) is
           \( \E(\Y) \), and its Hessian at \( \bzero \) is \( \E(\Y\Y\T) \). The Hessian of \( \log M_{\Y} \) at
           \( \bzero \) is \( \Cov(\Y) \).

3. *Affine maps.* For a constant \( k\times p \) matrix \( \A \) and \( k\times1 \) vector \( \bb \),
           \( M_{\A\Y+\bb}(\mathbf{s})=\exp(\mathbf{s}\T\bb)\,M_{\Y}(\A\T\mathbf{s}) \).

4. *Independence.* Partition \( \Y=(\Y_1\T,\Y_2\T)\T \) and \( \mathbf{t}=(\mathbf{t}_1\T,\mathbf{t}_2\T)\T \)
           conformably. Then \( \Y_1 \) and \( \Y_2 \) are independent iff
           \( M_{\Y}(\mathbf{t})=M_{\Y_1}(\mathbf{t}_1)\,M_{\Y_2}(\mathbf{t}_2) \) for all \( \mathbf{t} \) in a neighbourhood of
           \( \bzero \).
:::

:::

::: {.proof}
For (a), the scalar case is a theorem of analysis: a moment generating function that is
finite near the origin extends analytically, which pins down the characteristic function and
hence the distribution (Billingsley 1995). For vectors, fix \( \mathbf{a} \). For
small \( \lvert s\rvert \), \( M_{\mathbf{a}\T\Y}(s)=M_{\Y}(s\mathbf{a})=M_{\mathbf{V}}(s\mathbf{a})=M_{\mathbf{a}\T\mathbf{V}}(s) \), so
\( \mathbf{a}\T\Y \) and \( \mathbf{a}\T\mathbf{V} \) have the same distribution by the scalar case, and
@prp-rv-cramer-wold below finishes the argument.
For (b), finiteness near \( \bzero \) implies \( \E\exp(h'\lvert Y_i\rvert)<\infty \) for some \( h'>0 \)
and each \( i \), so every polynomial in the components is integrable and differentiation under
the expectation is justified. Then
\( \partial M/\partial t_i=\E\bigl(Y_ie^{\mathbf{t}\T\Y}\bigr) \) and
\( \partial^2M/\partial t_i\partial t_j=\E\bigl(Y_iY_je^{\mathbf{t}\T\Y}\bigr) \). At \( \mathbf{t}=\bzero \) these
are \( \E(Y_i) \) and \( \E(Y_iY_j) \). For \( K=\log M \), the chain rule gives
\( \nabla^2K=M^{-1}\nabla^2M-M^{-2}(\nabla M)(\nabla M)\T \), which at \( \bzero \), where \( M=1 \), is
\( \E(\Y\Y\T)-\bmu\bmu\T=\Cov(\Y) \). For (c),
\( \E\exp\bigl(\mathbf{s}\T(\A\Y+\bb)\bigr)=e^{\mathbf{s}\T\bb}\,\E\exp\bigl((\A\T\mathbf{s})\T\Y\bigr) \).
For (d), independence gives
\( \E\bigl(e^{\mathbf{t}_1\T\Y_1}e^{\mathbf{t}_2\T\Y_2}\bigr)=\E\bigl(e^{\mathbf{t}_1\T\Y_1}\bigr)\E\bigl(e^{\mathbf{t}_2\T\Y_2}\bigr) \).
Conversely, let \( \mathbf{V}_1 \) and \( \mathbf{V}_2 \) be independent with the distributions of \( \Y_1 \) and \( \Y_2 \).
The vector \( (\mathbf{V}_1\T,\mathbf{V}_2\T)\T \) has moment generating function
\( M_{\Y_1}(\mathbf{t}_1)M_{\Y_2}(\mathbf{t}_2) \), which by hypothesis equals \( M_{\Y} \) near \( \bzero \). By (a) it
has the same joint distribution as \( \Y \), so \( \Y_1 \) and \( \Y_2 \) are independent.
:::

The next fact explains why [Chapter 3](../ch03-multivariate-normal/index.html) can define the multivariate normal
through linear combinations.

::: {#prp-rv-cramer-wold}
[Distributions are determined by linear combinations]

If \( \mathbf{a}\T\Y \) and \( \mathbf{a}\T\mathbf{V} \) have the same distribution for every \( \mathbf{a}\in\Real^p \), then \( \Y \)
and \( \mathbf{V} \) have the same distribution.
:::

::: {.proof}
The characteristic function \( \phi_{\Y}(\mathbf{t})=\E\exp(i\mathbf{t}\T\Y) \) exists for every random vector
and determines its distribution (Billingsley 1995). Since
\( \phi_{\Y}(\mathbf{t})=\phi_{\mathbf{t}\T\Y}(1)=\phi_{\mathbf{t}\T\mathbf{V}}(1)=\phi_{\mathbf{V}}(\mathbf{t}) \) for every \( \mathbf{t} \), the
two distributions coincide.
:::

If \( \Y_1 \) and \( \Y_2 \) are independent, so are \( g(\Y_1) \) and \( h(\Y_2) \) for measurable \( g \) and \( h \).
Independence of statistics is usually proved this way. For normal vectors,
@thm-mvn-independence reduces the independence of two blocks to a zero
cross-covariance, and @thm-qf-indep-linear and @thm-qf-indep-quadratic build on that.

Independence of *several* vectors means that the joint distribution factors into all the
marginals, not merely that each pair is independent.

::: {#exm-rv-pairwise}
[Pairwise but not mutually independent]

Let \( S_1 \) and \( S_2 \) be independent, each equal to \( \pm1 \) with probability \( \tfrac12 \), and put
\( S_3=S_1S_2 \). Each \( S_j \) is \( \pm1 \) with probability \( \tfrac12 \). Each pair is independent: for
instance, \( \Pr(S_1=s,S_3=u)=\Pr(S_1=s,S_2=su)=\tfrac14 \) for all signs \( s,u \). So
\( \Cov(\mathbf{S})=\I_3 \). But the three are not mutually independent, since \( S_3 \) is a function of the
other two. A third-order moment shows the failure directly:
\( \E(S_1S_2S_3)=1 \), whereas mutual independence would force it to be \( 0 \). No condition on the
covariance matrix can detect this, because covariances involve only pairs.
:::

## Exercises

### A. Check your understanding

::: {#exr-rv-mgf-multinomial}
[A1]

Find the moment generating function of the multinomial count vector of @exm-rv-multinomial,
and recover its mean and covariance matrix from @thm-rv-mgf(b).
:::

### B. Practice

::: {#exr-rv-pairwise-many}
[B1]

Let \( S_1,\dots,S_k \) be independent random signs. For each nonempty \( A\subseteq\{1,\dots,k\} \) put
\( T_A=\prod_{i\in A}S_i \). Show that the \( 2^k-1 \) variables \( T_A \) are pairwise independent with
covariance matrix \( \I \), and find a triple of them that is not mutually independent.
:::

### C. Going deeper

::: {#exr-rv-lognormal-moments}
[C1]

Let \( \varphi \) be the \( \Normal(0,1) \) density and, for \( \lvert a\rvert\le1 \), let
\[
f_a(x)=\frac{\varphi(\log x)}{x}\bigl\{1+a\sin(2\pi\log x)\bigr\},\qquad x>0 .
\]

::: {.enumerate options="label=(\alph*)"}
1. Show that \( f_a \) is a probability density for every such \( a \), and that \( f_0 \) is the
           lognormal density, the density of \( e^{W} \) with \( W\sim\Normal(0,1) \).

2. Show that a random variable \( X_a \) with density \( f_a \) has \( \E(X_a^k)=e^{k^2/2} \) for every
           integer \( k\ge0 \), the same for all \( a \). *Hint:* substitute \( u=\log x \) and use
           \( e^{ku}\varphi(u)=e^{k^2/2}\varphi(u-k) \), then the symmetry of \( \varphi \).

3. Show that \( \E\,e^{tX_0}=\infty \) for every \( t>0 \). Conclude that the uniqueness statement
           @thm-rv-mgf(a) cannot be weakened to “all moments are finite and agree”: these
           distributions differ, yet no moment, and hence no moment-based argument, can separate
           them.
:::

:::

::: {.solution}
(a) Since \( \lvert\sin\rvert\le1 \) and \( \lvert a\rvert\le1 \), the brace is
nonnegative, so \( f_a\ge0 \). The substitution \( u=\log x \), \( dx=x\,du \), turns
\( \int_0^\infty f_a(x)\,dx \) into \( \int_{-\infty}^{\infty}\varphi(u)\{1+a\sin(2\pi u)\}\,du=1 \),
because \( \varphi \) is a density and \( u\mapsto\varphi(u)\sin(2\pi u) \) is odd and integrable. The
same substitution identifies \( f_0 \) as the density of \( e^W \).

(b) With the same substitution,
\[
\E(X_a^k)=\int_{-\infty}^{\infty}e^{ku}\varphi(u)\bigl\{1+a\sin(2\pi u)\bigr\}du
=e^{k^2/2}\Bigl[1+a\int_{-\infty}^{\infty}\varphi(u-k)\sin(2\pi u)\,du\Bigr].
\]
The remaining integral is \( \E\sin\{2\pi(W+k)\} \) with \( W\sim\Normal(0,1) \), which expands as
\( \sin(2\pi k)\,\E\cos(2\pi W)+\cos(2\pi k)\,\E\sin(2\pi W) \). For integer \( k \) the first
coefficient \( \sin(2\pi k) \) is zero, and \( \E\sin(2\pi W)=0 \) because \( W \) is symmetric. So
\( \E(X_a^k)=e^{k^2/2} \) whatever \( a \) is.

(c) For \( t>0 \), \( \E e^{tX_0}=\int_{-\infty}^\infty\exp(te^u)\varphi(u)\,du \), and
\( te^u-u^2/2\to\infty \) as \( u\to\infty \), so the integrand does not even tend to zero: the
integral diverges. Hence @thm-rv-mgf does not apply to the lognormal law, and indeed it
cannot: the perturbed laws have exactly the same moments of every order but different
densities, so they are genuinely different distributions.
:::

::: {#exr-rv-joint-mgf-exists}
[C2]

::: {.enumerate options="label=(\alph*)"}
1. Suppose each component of \( \Y \) has a moment generating function that is finite on
           \( \lvert s\rvert<h \). Show that \( M_{\Y}(\mathbf{t}) \) is finite whenever
           \( \max_i\lvert t_i\rvert<h/p \). *Hint:* convexity of the exponential gives
           \( \exp(\sum_it_iY_i)\le p^{-1}\sum_i\exp(pt_iY_i) \). So the hypothesis of @thm-rv-mgf is a
           statement about the components, but its conclusions are not.

2. Let \( Z\sim\Normal(0,1) \), let \( S \) be an independent random sign, and put \( W=SZ \). Show
           that
           \[
M_{(Z,W)}(\mathbf{t})=\tfrac12\exp\bigl\{\tfrac12(t_1+t_2)^2\bigr\}+\tfrac12\exp\bigl\{\tfrac12(t_1-t_2)^2\bigr\}.
\]
           Deduce that \( Z \) and \( W \) are each \( \Normal(0,1) \), that \( \Cov(Z,W)=0 \), and that they are
           *not* independent.

3. Compute \( \E(Z^2W^2) \) and compare it with the value it would have under independence.
:::

:::

::: {.solution}
(a) Each \( pt_i \) lies in \( (-h,h) \), so
\( \E\exp(pt_iY_i)=M_{Y_i}(pt_i)<\infty \). Averaging the hint over \( i \) and taking expectations
gives \( M_{\Y}(\mathbf{t})\le p^{-1}\sum_iM_{Y_i}(pt_i)<\infty \). The hint itself is the convexity
inequality \( \exp\bigl(p^{-1}\sum_i s_i\bigr)\le p^{-1}\sum_i e^{s_i} \) with \( s_i=pt_iY_i \).

(b) Conditioning on \( S \),
\( \E\exp(t_1Z+t_2SZ)=\tfrac12\E\exp\{(t_1+t_2)Z\}+\tfrac12\E\exp\{(t_1-t_2)Z\} \), and each term is
a normal moment generating function, which gives the display. Setting \( t_2=0 \) leaves
\( e^{t_1^2/2} \), so \( Z\sim\Normal(0,1) \); setting \( t_1=0 \) leaves the same, so \( W\sim\Normal(0,1) \).
Expanding to second order, \( M_{(Z,W)}(\mathbf{t})=1+\tfrac12(t_1^2+t_2^2)+O(\norm{\mathbf{t}}^4) \) with no
\( t_1t_2 \) term, so \( \Cov(Z,W)=0 \) by @thm-rv-mgf(b). Were they independent, the joint moment
generating function would be \( \exp\{\tfrac12(t_1^2+t_2^2)\} \) by @thm-rv-mgf(d). Along
\( \mathbf{t}=(s,s)\T \) the display equals \( \tfrac12(e^{2s^2}+1)=1+s^2+s^4+O(s^6) \), whereas
\( \exp(s^2)=1+s^2+\tfrac12s^4+O(s^6) \). The two disagree in every neighbourhood of the origin,
so \( Z \) and \( W \) are dependent, although both are normal and uncorrelated.

(c) \( W^2=Z^2 \), so \( \E(Z^2W^2)=\E(Z^4)=3 \), whereas independence would give
\( \E(Z^2)\E(W^2)=1 \). This is the failure of @exm-rv-circle in a normal disguise; joint
normality, defined in [Chapter 3](../ch03-multivariate-normal/index.html), is exactly what rules it
out.
:::
