# Moment generating functions and independence

Means and covariance matrices summarize a distribution. They do not determine it, and they
cannot certify independence. For that we need a transform of the whole distribution. This
section collects the facts that [Chapter 3](../ch03-multivariate-normal/index.html) and [Chapter 4](../ch04-quadratic-forms/index.html) use.

::: {#def-rv-mgf}
[Moment generating function]

The **moment generating function** of a \( p\times1 \) random vector \( \Y \) is
\[
M_{\Y}(\bm t)=\E\exp(\bm t\T\Y),\qquad \bm t\in\Real^p ,
\]
where the expectation may be \( +\infty \). We say that \( M_{\Y} \) *exists* if it is finite for
all \( \bm t \) in some neighbourhood \( \{\norm{\bm t}<h\} \) of the origin, \( h>0 \).
:::

::: {#thm-rv-mgf}
[Properties of moment generating functions]

Suppose \( M_{\Y} \) exists.

::: {.enumerate options="label=(\alph*)"}
1. *Uniqueness.* If \( M_{\bm V} \) also exists and \( M_{\Y}=M_{\bm V} \) on a neighbourhood of
           \( \bzero \), then \( \Y \) and \( \bm V \) have the same distribution.

2. *Moments.* All moments of \( \Y \) are finite. The gradient of \( M_{\Y} \) at \( \bzero \) is
           \( \E(\Y) \), and its Hessian at \( \bzero \) is \( \E(\Y\Y\T) \). The Hessian of \( \log M_{\Y} \) at
           \( \bzero \) is \( \Cov(\Y) \).

3. *Affine maps.* For a constant \( k\times p \) matrix \( \A \) and \( k\times1 \) vector \( \bb \),
           \( M_{\A\Y+\bb}(\bm s)=\exp(\bm s\T\bb)\,M_{\Y}(\A\T\bm s) \).

4. *Independence.* Partition \( \Y=(\Y_1\T,\Y_2\T)\T \) and \( \bm t=(\bm t_1\T,\bm t_2\T)\T \)
           conformably. Then \( \Y_1 \) and \( \Y_2 \) are independent iff
           \( M_{\Y}(\bm t)=M_{\Y_1}(\bm t_1)\,M_{\Y_2}(\bm t_2) \) for all \( \bm t \) in a neighbourhood of
           \( \bzero \).
:::

:::

::: {.proof}
For (a), the scalar case is a theorem of analysis: a moment generating function that is
finite near the origin extends analytically, which pins down the characteristic function and
hence the distribution (Billingsley 1995). For vectors, fix \( \bm a \). For
small \( \lvert s\rvert \), \( M_{\bm a\T\Y}(s)=M_{\Y}(s\bm a)=M_{\bm V}(s\bm a)=M_{\bm a\T\bm V}(s) \), so
\( \bm a\T\Y \) and \( \bm a\T\bm V \) have the same distribution by the scalar case, and
@prp-rv-cramer-wold below finishes the argument.
For (b), finiteness near \( \bzero \) implies \( \E\exp(h'\lvert Y_i\rvert)<\infty \) for some \( h'>0 \)
and each \( i \), so every polynomial in the components is integrable and differentiation under
the expectation is justified. Then
\( \partial M/\partial t_i=\E\bigl(Y_ie^{\bm t\T\Y}\bigr) \) and
\( \partial^2M/\partial t_i\partial t_j=\E\bigl(Y_iY_je^{\bm t\T\Y}\bigr) \). At \( \bm t=\bzero \) these
are \( \E(Y_i) \) and \( \E(Y_iY_j) \). For \( K=\log M \), the chain rule gives
\( \nabla^2K=M^{-1}\nabla^2M-M^{-2}(\nabla M)(\nabla M)\T \), which at \( \bzero \), where \( M=1 \), is
\( \E(\Y\Y\T)-\bmu\bmu\T=\Cov(\Y) \). For (c),
\( \E\exp\bigl(\bm s\T(\A\Y+\bb)\bigr)=e^{\bm s\T\bb}\,\E\exp\bigl((\A\T\bm s)\T\Y\bigr) \).
For (d), independence gives
\( \E\bigl(e^{\bm t_1\T\Y_1}e^{\bm t_2\T\Y_2}\bigr)=\E\bigl(e^{\bm t_1\T\Y_1}\bigr)\E\bigl(e^{\bm t_2\T\Y_2}\bigr) \).
Conversely, let \( \bm V_1 \) and \( \bm V_2 \) be independent with the distributions of \( \Y_1 \) and \( \Y_2 \).
The vector \( (\bm V_1\T,\bm V_2\T)\T \) has moment generating function
\( M_{\Y_1}(\bm t_1)M_{\Y_2}(\bm t_2) \), which by hypothesis equals \( M_{\Y} \) near \( \bzero \). By (a) it
has the same joint distribution as \( \Y \), so \( \Y_1 \) and \( \Y_2 \) are independent.
:::

The next fact explains why [Chapter 3](../ch03-multivariate-normal/index.html) can define the multivariate normal
through linear combinations.

::: {#prp-rv-cramer-wold}
[Distributions are determined by linear combinations]

If \( \bm a\T\Y \) and \( \bm a\T\bm V \) have the same distribution for every \( \bm a\in\Real^p \), then \( \Y \)
and \( \bm V \) have the same distribution.
:::

::: {.proof}
The characteristic function \( \phi_{\Y}(\bm t)=\E\exp(i\bm t\T\Y) \) exists for every random vector
and determines its distribution (Billingsley 1995). Since
\( \phi_{\Y}(\bm t)=\phi_{\bm t\T\Y}(1)=\phi_{\bm t\T\bm V}(1)=\phi_{\bm V}(\bm t) \) for every \( \bm t \), the
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
\( \Cov(\bm S)=\I_3 \). But the three are not mutually independent, since \( S_3 \) is a function of the
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
