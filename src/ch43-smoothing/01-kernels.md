# Kernel and local polynomial regression

Drop the linear predictor. Keep one covariate, keep additive errors, and assume
nothing about the shape of the regression function except that it is smooth:
\[
y_i = f(x_i)+\varepsilon_i,\qquad
\E(\varepsilon_i)=0,\quad \Var(\varepsilon_i)=\sigma^2,\quad i=1,\dots,n,
\]{#eq-smo-model}

with independent errors and design points \( x_1,\dots,x_n \) in an interval
\( [a,b] \). The unknown is now a *function*, not a vector of coefficients, and
there is no finite-dimensional parameter to estimate consistently. What replaces
it is an assumption of *local* regularity: \( f \) does not change much over a
short interval, so observations with covariate values near \( x \) carry
information about \( f(x) \).

That sentence is the whole idea of this section. Fix a target point \( x \), take
the observations whose covariates lie in a window around it, fit something simple
to them by least squares, and read off the fitted value at \( x \); move \( x \)
across the interval and the fitted values trace a curve. Everything else is
bookkeeping: which observations count as "near", how much they count, and what
"something simple" means.

## Local averages

::: {#exm-smo-running-mean}
[The running mean]

The crudest smoother averages the responses in a window: fix a half-width
\( h \) and set
\( \hat f(x)=\bigl\lvert N(x)\bigr\rvert^{-1}\sum_{i\in N(x)}y_i \) with
\( N(x)=\{i:\lvert x_i-x\rvert\le h\} \). Three things go wrong at once. The curve
is a step function, since \( \hat f \) changes only when a point enters or leaves
the window; every point in the window counts the same, at the centre or at the
edge; and near \( x=a \) the window is one-sided, so the average estimates
\( f \) somewhere to the right of \( a \).
:::

Weighting repairs the first two defects: count a point by a weight that decreases
smoothly with distance and vanishes outside the window. The third is deeper, and
is repaired by fitting a line rather than a constant. Both repairs fit in one
definition.

::: {#def-smo-kernel}
[Kernel weights and the local polynomial estimator]

A **kernel** is a bounded function \( K \) with \( \int K(u)\,du=1 \) and
\( K(u)=K(-u) \), here always nonnegative and supported on \( [-1,1] \). Its
moments are \( \mu_j=\int u^jK(u)\,du \) and
\( \nu_j=\int u^jK(u)^2\,du \), with \( R(K)=\nu_0 \). For a **bandwidth**
\( h>0 \) write \( K_h(u)=h^{-1}K(u/h) \).

Fix a degree \( q\ge0 \). The **local polynomial estimator of degree \( q \)** at
the point \( x \) is \( \hat f(x)=\hat b_0 \), where
\( (\hat b_0,\dots,\hat b_q) \) minimizes
\[
\sum_{i=1}^n K_h(x_i-x)\Bigl\{y_i-\sum_{j=0}^{q}b_j(x_i-x)^j\Bigr\}^2 .
\]{#eq-smo-localpoly}

For \( q=0 \) this is the **local constant** or **Nadaraya–Watson** estimator
\[
\hat f(x)=\frac{\sum_iK_h(x_i-x)\,y_i}{\sum_iK_h(x_i-x)} ,
\]{#eq-smo-nw}

and for \( q=1 \) the **local linear** estimator.
:::

Three kernels are used throughout: the uniform kernel
\( K(u)=\tfrac12 \) on \( [-1,1] \), which gives the running mean back; the
**Epanechnikov** kernel \( K(u)=\tfrac34(1-u^2) \) on \( [-1,1] \), with
\( \mu_2=1/5 \) and \( R(K)=3/5 \); and the standard normal density, whose
unbounded support does not literally satisfy the conditions below, though the
conclusions are the same. The kernel matters far less than \( h \), as
@exr-smo-efficiency quantifies.

@eq-smo-localpoly is a weighted least squares
problem (@thm-het-wls) with model matrix and weight matrix
\[
\bU_x=\begin{pmatrix}
1 & x_1-x & \cdots & (x_1-x)^q\\
\vdots & \vdots & & \vdots\\
1 & x_n-x & \cdots & (x_n-x)^q
\end{pmatrix},
\qquad
\W_x=\diag\{K_h(x_i-x)\}_{i=1}^{n} .
\]

Nothing in it is new; what is new is that it is solved afresh at every \( x \),
and that the weights, not a global model, decide which observations matter.

## The estimator is linear, and it reproduces polynomials

::: {#prp-smo-equivalent-kernel}
[Equivalent kernel weights]

Suppose \( \bU_x\T\W_x\bU_x \) is nonsingular. Then
\( \hat f(x)=\sum_{i}\ell_i(x)\,y_i \) with
\[
\ell_i(x)=\vect{e}_1\T(\bU_x\T\W_x\bU_x)^{-1}\bU_x\T\W_x\vect{e}_i ,
\]{#eq-smo-equivalent}

where \( \vect{e}_1 \) is the first standard basis vector, and

::: {.enumerate options="label=(\alph*)"}
1. \( \sum_i\ell_i(x)(x_i-x)^j=\delta_{j0} \) for \( j=0,\dots,q \), where
   \( \delta_{j0} \) is one if \( j=0 \) and zero otherwise;

2. consequently \( \hat f \) reproduces polynomials: if \( f \) is a polynomial of
   degree at most \( q \) then \( \sum_i\ell_i(x)f(x_i)=f(x) \) for every \( x \),
   whatever the design, the kernel and the bandwidth;

3. the weights do not depend on \( \y \), so
   \( \Var\{\hat f(x)\}=\sigma^2\sum_i\ell_i(x)^2 \).
:::

:::

::: {.proof}
The weighted normal equations (@thm-proj-normal-equations) are
\( \bU_x\T\W_x\bU_x\hat{\bb}=\bU_x\T\W_x\y \), and \( \hat f(x)=\vect{e}_1\T\hat{\bb} \),
which is @eq-smo-equivalent. For (a), the \( j \)th column of \( \bU_x \) has
entries \( (x_i-x)^j \), so
\( \sum_i\ell_i(x)(x_i-x)^j=\vect{e}_1\T(\bU_x\T\W_x\bU_x)^{-1}\bU_x\T\W_x\bU_x\vect{e}_{j+1}
=\vect{e}_1\T\vect{e}_{j+1} \). Part (b) follows by writing
\( f(x_i)=\sum_{j\le q}c_j(x_i-x)^j \) with \( c_0=f(x) \), and (c) because
\( \Cov(\Y)=\sigma^2\I \).
:::

Part (a) is the engine of everything that follows. The weights sum to one, so
the estimator has no bias against a constant; for \( q\ge1 \) they also
annihilate \( (x_i-x) \), so it has no bias against a linear trend either. That
is exactly what the running mean lacked.

The weight profile \( u\mapsto\ell_i(x) \), with \( u=(x_i-x)/h \), is the
**equivalent kernel** of the estimator. For \( q=1 \) at an interior point it is
asymptotically \( K(u)/\{nhg(x)\} \): a local linear fit *is* a kernel average,
with the kernel unchanged. For \( q=3 \) it is
\( (\mu_4-\mu_2u^2)K(u)/(\mu_4-\mu_2^2) \), which takes negative values — the
price of killing the \( h^2 \) bias term.

## Bias and variance

Four assumptions make the asymptotics exact rather than heuristic. A random
design would work too, at the cost of replacing every statement by one that holds
in probability; the proof below is for a *regular* design, where the argument is a
clean Riemann-sum estimate.

**(A) Design.** \( x_i=G^{-1}\{(i-\tfrac12)/n\} \) for \( i=1,\dots,n \), where
\( G \) is a distribution function on \( [a,b] \) whose density \( g \) is
continuously differentiable and bounded away from zero.

**(B) Kernel.** \( K \) is as in @def-smo-kernel, of bounded variation, with
\( \mu_2>0 \).

**(C) Function.** \( f \) is twice continuously differentiable on \( [a,b] \).

**(D) Bandwidth.** \( h=h_n\to0 \) and \( nh\to\infty \).

Assumption (A) holds exactly for an equally spaced or quantile-spaced design. The
next lemma, the only analysis in this chapter, says that a kernel-weighted sum
over such a design is an integral plus a bounded error.

::: {#lem-smo-moment-sums}
[Kernel sums are integrals]

Let \( \varphi \) be of bounded variation on \( \Real \) with support in
\( [-1,1] \) and total variation \( V(\varphi) \). Under (A), for every \( x \) and
every \( h>0 \),
\[
\Bigl\lvert\ \sum_{i=1}^n\varphi\Bigl(\frac{x_i-x}{h}\Bigr)
-nh\int_{(a-x)/h}^{(b-x)/h}\varphi(u)\,g(x+hu)\,du\ \Bigr\rvert\le V(\varphi) .
\]{#eq-smo-riemann}

:::

::: {.proof}
Put \( \Psi(v)=\varphi\{(G^{-1}(v)-x)/h\} \) for \( v\in[0,1] \). Since
\( G^{-1} \) is increasing, \( \Psi \) has the same total variation as
\( t\mapsto\varphi\{(t-x)/h\} \), namely \( V(\varphi) \). The sum on the left is
\( \sum_{i}\Psi\{(i-\tfrac12)/n\} \), a midpoint sum. Write \( I_i \) for the
\( i \)th subinterval of \( [0,1] \) of length \( 1/n \) and \( V_i \) for the
variation of \( \Psi \) on \( I_i \). For \( v\in I_i \),
\( \lvert\Psi\{(i-\tfrac12)/n\}-\Psi(v)\rvert\le V_i \), so
\( \lvert n^{-1}\Psi\{(i-\tfrac12)/n\}-\int_{I_i}\Psi\rvert\le V_i/n \).
Summing over \( i \) and using \( \sum_iV_i=V(\varphi) \) gives
\( \lvert\sum_i\Psi\{(i-\tfrac12)/n\}-n\int_0^1\Psi\rvert\le V(\varphi) \).
Finally substitute \( t=G^{-1}(v) \), so that \( dv=g(t)\,dt \), and then
\( t=x+hu \).
:::

Applied with \( \varphi(u)=K(u)u^j \) and with \( \varphi=K^2u^j \), the lemma
turns every sum that appears below into a moment of the kernel. Write
\[
S_{n,j}=\sum_{i=1}^nK_h(x_i-x)(x_i-x)^j ,
\]
and let \( \A=(\mu_{j+k})_{j,k=0}^{q} \) and
\( \A^{*}=(\nu_{j+k})_{j,k=0}^{q} \) be the kernel's moment matrices, which are
positive definite because \( \mu_{j+k} \) and \( \nu_{j+k} \) are the Gram matrices
of \( 1,u,\dots,u^q \) in \( L^2(K\,du) \) and \( L^2(K^2du) \).

::: {#thm-smo-local-poly}
[Bias and variance of a local polynomial estimate]

Assume (A)–(D) and let \( x \) satisfy \( [x-h,x+h]\subset(a,b) \).

::: {.enumerate options="label=(\alph*)"}
1. *(Local linear, \( q=1 \).)*
   \[
   \E\{\hat f(x)\}-f(x)=\tfrac12\mu_2f''(x)h^2+o(h^2),
   \qquad
   \Var\{\hat f(x)\}=\frac{\sigma^2R(K)}{n h\,g(x)}\{1+o(1)\} .
   \]

2. *(Local constant, \( q=0 \).)* The variance is the same to first order, and
   \[
   \E\{\hat f(x)\}-f(x)
   =\mu_2h^2\Bigl\{\tfrac12f''(x)+\frac{f'(x)g'(x)}{g(x)}\Bigr\}+o(h^2) .
   \]{#eq-smo-design-bias}

3. *(At the boundary.)* Let \( x=a+ch \) with \( c\in[0,1) \) fixed, and write
   \( \mu_{j,c}=\int_{-c}^{1}u^jK(u)\,du \). Then the local constant estimator has
   bias \( f'(a)\,h\,\mu_{1,c}/\mu_{0,c}+O(h^2) \), of exact order \( h \) unless
   \( f'(a)=0 \), while the local linear estimator still has bias \( O(h^2) \).

4. *(Optimal bandwidth.)* For the local linear estimator with
   \( f''(x)\ne0 \), the pointwise mean squared error
   \( \tfrac14\mu_2^2f''(x)^2h^4+\sigma^2R(K)/\{nhg(x)\} \) is minimized at
   \[
   h^{*}(x)=\Bigl\{\frac{\sigma^2R(K)}{g(x)\,\mu_2^2f''(x)^2}\Bigr\}^{1/5}n^{-1/5},
   \]{#eq-smo-hopt}

   and the minimum value is
   \( \tfrac54\{\mu_2^2f''(x)^2\}^{1/5}\{\sigma^2R(K)/g(x)\}^{4/5}n^{-4/5} \).
:::

:::

::: {.proof}
Write \( u_i=(x_i-x)/h \), let \( \bU \) have rows \( (1,u_i,\dots,u_i^q) \) and
\( \W=\W_x \). Because \( \bU_x=\bU\,\diag(1,h,\dots,h^q) \) and the first diagonal
entry is one, @eq-smo-equivalent gives
\( \ell_i(x)=\vect{e}_1\T(\bU\T\W\bU)^{-1}\bU\T\W\vect{e}_i \).

*Moment sums.* By @lem-smo-moment-sums with \( \varphi(u)=K(u)u^j \), and because
\( [x-h,x+h]\subset(a,b) \) makes the range of integration all of \( [-1,1] \),
\[
S_{n,j}=h^{j-1}\sum_iK(u_i)u_i^j
=nh^{j}\Bigl\{\int K(u)u^jg(x+hu)\,du+O\bigl((nh)^{-1}\bigr)\Bigr\} .
\]
Since \( g\in C^1 \), \( g(x+hu)=g(x)+hug'(x+\theta hu) \) with
\( \theta\in(0,1) \), and \( g' \) is uniformly continuous, so
\[
S_{n,j}=nh^{j}\bigl\{g(x)\mu_j+hg'(x)\mu_{j+1}+o(h)+O((nh)^{-1})\bigr\} .
\]{#eq-smo-moment-sum}

Because \( (\bU\T\W\bU)_{jk}=\sum_iK_h(x_i-x)u_i^{j+k}=h^{-(j+k)}S_{n,j+k} \), this
says \( \bU\T\W\bU=n\{g(x)\A+O(h)+O((nh)^{-1})\} \). In the same way, with
\( \varphi=K^2u^j \),
\( \bU\T\W^2\bU=(n/h)\{g(x)\A^{*}+O(h)+O((nh)^{-1})\} \). Assumption (D) makes both
error terms vanish, so \( \bU\T\W\bU \) is eventually nonsingular and
@prp-smo-equivalent-kernel applies.

*Variance.* By @prp-smo-equivalent-kernel(c),
\[
\Var\{\hat f(x)\}=\sigma^2\,\vect{e}_1\T(\bU\T\W\bU)^{-1}(\bU\T\W^2\bU)(\bU\T\W\bU)^{-1}\vect{e}_1
=\frac{\sigma^2\,\vect{e}_1\T\A^{-1}\A^{*}\A^{-1}\vect{e}_1}{nh\,g(x)}\{1+o(1)\} .
\]
For \( q=0 \), \( \A=(1) \) and \( \A^{*}=(R(K)) \). For \( q=1 \) the symmetry of
\( K \) makes \( \mu_1=0 \), so \( \A=\diag(1,\mu_2) \), \( \A^{-1}\vect{e}_1=\vect{e}_1 \)
and the quadratic form is \( \nu_0=R(K) \). Both give the stated variance.

*Bias, \( q=1 \).* By @prp-smo-equivalent-kernel(a) the weights annihilate
constants and linear terms, so with a second-order Taylor expansion of \( f \)
about \( x \), \( f(x_i)=f(x)+f'(x)(x_i-x)+\tfrac12f''(\xi_i)(x_i-x)^2 \) for some
\( \xi_i \) between \( x \) and \( x_i \),
\[
\E\{\hat f(x)\}-f(x)
=\tfrac12f''(x)\sum_i\ell_i(x)(x_i-x)^2
+\tfrac12\sum_i\ell_i(x)\{f''(\xi_i)-f''(x)\}(x_i-x)^2 .
\]
For the first sum, \( \sum_i\ell_i(x)(x_i-x)^2
=h^2\vect{e}_1\T(\bU\T\W\bU)^{-1}\bU\T\W\bu^{(2)} \) with \( \bu^{(2)}_i=u_i^2 \), and
\( (\bU\T\W\bu^{(2)})_j=\sum_iK_h(x_i-x)u_i^{j+2}=h^{-(j+2)}S_{n,j+2} \),
which by @eq-smo-moment-sum equals \( n\{g(x)\mu_{j+2}+O(h)+O((nh)^{-1})\} \). Hence
\[
\sum_i\ell_i(x)(x_i-x)^2
=h^2\bigl[\vect{e}_1\T\A^{-1}(\mu_2,\mu_3)\T+O(h)+O((nh)^{-1})\bigr]
=h^2\{\mu_2+o(1)\},
\]
using \( \mu_3=0 \). For the second sum, only \( \lvert x_i-x\rvert\le h \)
contributes, so \( \lvert f''(\xi_i)-f''(x)\rvert\le\omega(h) \), the modulus of
continuity of \( f'' \), which tends to zero. Since
\( (\bU\T\W\bU)^{-1}=O(n^{-1}) \) entrywise and
\( \sum_iK_h(x_i-x)\lvert u_i\rvert^{j}=O(n) \) by @lem-smo-moment-sums, we get
\( \sum_i\lvert\ell_i(x)\rvert(x_i-x)^2=O(h^2) \), and the second sum is
\( O\{\omega(h)h^2\}=o(h^2) \). This proves (a).

*Bias, \( q=0 \).* Now only constants are annihilated, so a term in
\( f'(x) \) survives:
\[
\E\{\hat f(x)\}-f(x)=f'(x)\frac{S_{n,1}}{S_{n,0}}
+\tfrac12f''(x)\frac{S_{n,2}}{S_{n,0}}+o(h^2) ,
\]
the remainder being controlled as above. By @eq-smo-moment-sum and \( \mu_1=0 \),
\( S_{n,0}=n\{g(x)+o(h)\} \) and
\( S_{n,1}=nh^2\{g'(x)\mu_2+o(1)\} \), so
\( S_{n,1}/S_{n,0}=h^2\mu_2g'(x)/g(x)+o(h^2) \), while
\( S_{n,2}/S_{n,0}=h^2\mu_2+o(h^2) \). This is @eq-smo-design-bias.

*Boundary.* For \( x=a+ch \) the range of integration in
@lem-smo-moment-sums is \( [-c,(b-x)/h] \), which for large \( n \) is
\( [-c,1] \). Every moment \( \mu_j \) above is therefore replaced by
\( \mu_{j,c} \) and the moment matrix by
\( \A_c=(\mu_{j+k,c}) \), still positive definite for every \( c\in[0,1) \) because the
monomials remain linearly independent in \( L^2(K\,du) \) on \( [-c,1] \). For
\( q=0 \) the leading term is \( f'(a)hS_{n,1}/(hS_{n,0})
=f'(a)h\,\mu_{1,c}/\mu_{0,c}+O(h^2) \), and \( \mu_{1,c}\neq0 \) whenever
\( c<1 \), since the truncated kernel is no longer symmetric. For \( q=1 \) the
identity \( \sum_i\ell_i(x)(x_i-x)=0 \) of @prp-smo-equivalent-kernel(a) holds
exactly, *at every point, boundary or not*, so the first-order term is absent and
the same computation as in (a) leaves
\( \tfrac12f''(a)h^2\vect{e}_1\T\A_c^{-1}(\mu_{2,c},\mu_{3,c})\T+o(h^2) \).

*Optimal bandwidth.* Write the mean squared error as \( Ah^4+B/(nh) \) with
\( A=\tfrac14\mu_2^2f''(x)^2 \) and \( B=\sigma^2R(K)/g(x) \). Differentiating,
\( 4Ah^3=B/(nh^2) \), that is \( h^5=B/(4An) \), which is @eq-smo-hopt; at that
\( h \) the two terms are in the ratio \( 1:4 \), so the minimum is
\( 5Ah^4=5A^{1/5}\{B/(4n)\}^{4/5} \), the stated value.
:::

::: {.idea}
Local linear beats local constant twice, and for the same reason both times: the
weights of a fit of degree \( q\ge1 \) reproduce linear functions *exactly*,
whatever the design and wherever the target point. The local linear estimator is
therefore blind to how the design points are distributed (no term in
\( g'/g \)) and to the one-sidedness of a window at the boundary (no term of
order \( h \)). The local constant estimator is not, and neither extra bias is
removed by more data, only by a smaller \( h \).
:::

Part (d) is the rate that every nonparametric method in this chapter obeys. A
parametric estimator has variance of order \( n^{-1} \) and no bias; here only
the \( O(nh) \) points in the window contribute, so the variance is of order
\( (nh)^{-1} \) against a bias of order \( h^2 \). Balancing the two forces
\( h\asymp n^{-1/5} \) and leaves a mean squared error of order \( n^{-4/5} \),
strictly slower than \( n^{-1} \). The optimal bandwidth shrinks slowly: ten
times the data calls for a window only
\( 1.585 \) times narrower.

::: {.warning}
@eq-smo-hopt is not a recipe. It involves the unknown \( f''(x) \), the unknown
\( \sigma^2 \) and the design density, and it is a *pointwise* optimum, so the
ideal bandwidth differs at every \( x \). Any formula that gives a bandwidth
without an estimate of \( f'' \) is smuggling one in.
:::

## What the numbers look like

::: {#exm-smo-bandwidth}
[Bias, variance and the optimal bandwidth]

Take \( n=200 \) equally spaced design points on \( [0,1] \)
(so \( g\equiv1 \)), \( \sigma=0.25 \) and
\[
f(x)=2x+e^{-25(x-0.35)^2}-\tfrac34 e^{-50(x-0.80)^2},
\]
a straight line with a bump and a narrower dip. At the top of the bump,
\( x=0.35 \), \( f''(x)=-50.06 \), and with the
Epanechnikov kernel @eq-smo-hopt gives
\( h^{*}=0.0715 \). Averaging \( 2000 \) data sets at
each of \( 29 \) bandwidths
([Figure 43.1.2](01-kernels.html#fig-smo-mse)), the bias at \( h^{*} \) is
\( 0.98 \) times the asymptotic value
\( \tfrac12\mu_2f''(x)h^2 \) and the variance
\( 1.02 \) times \( \sigma^2R(K)/(nh) \), while the
empirically best bandwidth is \( 0.0703 \). For a
first-order expansion at \( n=200 \) this is closer agreement than one has any
right to expect.

[Figure 43.1.1](01-kernels.html#fig-smo-boundary) shows fits at three bandwidths.
At the left endpoint the two estimators part company. With
\( h=0.08 \), the local constant estimator has bias
\( 0.0938 \) there and the local linear estimator
\( -0.0058 \), sixteen times smaller. Halving the
bandwidth divides the first by \( 2.10 \) and the
second by \( 4.38 \): order \( h \) against order
\( h^2 \), as @thm-smo-local-poly(c) says.
:::

::: {when-format="html"}
![**Figure 43.1.1.** (a) Local linear fits to one data set at three bandwidths: too
small, about right, too large. (b) Bias near the left endpoint over 2000
replicates. Inside \( [0,h] \) the local constant window is one-sided, with bias
of order \( h \); the local linear bias stays of order \( h^2
\) (@thm-smo-local-poly).](local_linear_boundary.svg){#fig-smo-boundary width=100%}
:::

::: {when-format="pdf"}
![(a) Local linear fits to one data set at three bandwidths: too
small, about right, too large. (b) Bias near the left endpoint over 2000
replicates. Inside \( [0,h] \) the local constant window is one-sided, with bias
of order \( h \); the local linear bias stays of order \( h^2
\) (@thm-smo-local-poly).](local_linear_boundary.pdf){width=100%}
:::

```{.python .run #cell-local-polynomial-fit}
import numpy as np

A1, C1, M1 = 1.0, 25.0, 0.35          # the bump
A2, C2, M2 = -0.75, 50.0, 0.80        # the dip


def f_true(x):
    """The mean function: a straight line plus a bump and a narrower dip."""
    return 2 * x + A1 * np.exp(-C1 * (x - M1) ** 2) + A2 * np.exp(-C2 * (x - M2) ** 2)


def epanechnikov(u):
    """K(u) = (3/4)(1 - u^2) on [-1, 1], zero outside."""
    return np.where(np.abs(u) < 1, 0.75 * (1 - u**2), 0.0)


def local_poly_weights(x0, x, h, degree):
    """The weights l_i(x0) of the local polynomial fit: f-hat(x0) = sum_i l_i(x0) y_i."""
    w = epanechnikov((x - x0) / h) / h                 # kernel weights K_h(x_i - x0)
    U = np.vander(x - x0, degree + 1, increasing=True)  # columns (x_i - x0)^j
    A = U.T @ (w[:, None] * U)                          # the weighted moment matrix
    e1 = np.zeros(degree + 1)
    e1[0] = 1.0
    return np.linalg.solve(A, e1) @ U.T * w             # row vector of weights


def local_poly(grid, x, y, h, degree):
    """Local polynomial estimate of degree `degree` on a grid of target points."""
    return np.array([local_poly_weights(x0, x, h, degree) @ y for x0 in grid])


n, sigma = 200, 0.25
x = (np.arange(1, n + 1) - 0.5) / n          # the regular design of assumption (A)
rng = np.random.default_rng(20250921)
y = f_true(x) + sigma * rng.normal(size=n)


def endpoint_bias(h, degree):
    w = local_poly_weights(0.0, x, h, degree)
    return w @ f_true(x) - f_true(0.0)


print("bias at the left endpoint, local constant:", round(endpoint_bias(0.08, 0), 4))
print("bias at the left endpoint, local linear:  ", round(endpoint_bias(0.08, 1), 4))
```

::: {when-format="html"}
![**Figure 43.1.2.** Squared bias, variance and mean squared error of the local
linear estimate at \( x=0.35 \) against the bandwidth, from 2000 data sets
(solid) and from @thm-smo-local-poly(a) (dotted); the vertical line is
\( h^{*} \). Logarithmic scales.](bias_variance_bandwidth.svg){#fig-smo-mse width=58%}
:::

::: {when-format="pdf"}
![Squared bias, variance and mean squared error of the local
linear estimate at \( x=0.35 \) against the bandwidth, from 2000 data sets
(solid) and from @thm-smo-local-poly(a) (dotted); the vertical line is
\( h^{*} \). Logarithmic scales.](bias_variance_bandwidth.pdf){width=58%}
:::

## Choosing the bandwidth in practice, and what kernels cannot do

Two families of rules are used. A **plug-in** rule estimates the unknowns in
@eq-smo-hopt — usually an integrated version of it, with \( f''(x)^2 \) replaced
by \( \int f''(t)^2g(t)\,dt \) — from a pilot fit, and iterates; Ruppert,
Sheather and Wand (1995) is the standard version. A **cross-validatory** rule
estimates prediction error directly, and since the local polynomial estimator is
linear in \( \y \), everything in [Section 43.5](05-choosing-lambda.html) applies
to it with \( \tr(\bS_h) \) for \( \tr(\bS_\lambda) \).

::: {.remark}
[Loess]

The most widely used kernel smoother in practice is Cleveland's **loess**
(Cleveland 1979), which makes two changes. The bandwidth is a
*nearest-neighbour* one: \( h(x) \) is the distance from \( x \) to its
\( \lceil\alpha n\rceil \)th nearest design point, so the window holds a fixed
fraction of the data and widens where the design is sparse. And the fit is
iterated with bisquare weights on the residuals, the reweighting
of @prp-res-irls, which makes it resistant to a few wild responses. Neither
change affects @thm-smo-local-poly at an interior point where \( g \) is
continuous, since \( h(x)\approx\alpha/\{2g(x)\} \) is then a fixed bandwidth
divided by \( g(x) \).
:::

There are two reasons this chapter does not stop here. The first is that a
kernel estimate is the solution of no global optimization problem: there is no
criterion whose minimizer it is and no coefficient vector to report, so there is
no obvious way to add a second covariate, a non-normal response or a random
effect. The second is dimension. With \( d \) covariates the window is a
\( d \)-dimensional ball holding of order \( nh^d \) points, and the optimal rate
degrades to \( n^{-4/(4+d)} \) — the *curse of dimensionality*. At \( d=5 \) it is
\( n^{-4/9} \), so matching what \( 100 \) observations give in one dimension
needs of the order of \( 10^5 \) of them;
[Chapter 44](../ch44-additive-models/index.html) escapes by assuming additivity,
not by smoothing harder. The rest of this chapter writes the estimate as the
minimizer of a penalized criterion, which keeps a coefficient vector, a linear
operator and a criterion, all of which generalize.

## Exercises

### A. Check your understanding

::: {#exr-smo-weights-sum}
[A1]

Show directly from @eq-smo-nw that the Nadaraya–Watson weights are nonnegative and
sum to one, and deduce that \( \min_iy_i\le\hat f(x)\le\max_iy_i \). Show by
example that the local linear estimator has no such property.
:::

::: {.solution}
The weights \( K_h(x_i-x)/\sum_jK_h(x_j-x) \) are nonnegative because \( K\ge0 \)
and sum to one by construction, so \( \hat f(x) \) is a convex combination of the
\( y_i \). For the local linear estimator take \( n=2 \), \( x_1=0 \),
\( x_2=1 \), both inside the window at \( x=2 \) with equal kernel weights: the
line through the two points is extrapolated, giving
\( \hat f(2)=2y_2-y_1 \), outside \( [y_1,y_2] \) when \( y_2>y_1 \).
:::

### B. Practice

::: {#exr-smo-uniform-kernel}
[B1]

Compute \( \mu_2 \) and \( R(K) \) for the uniform kernel
\( K(u)=\tfrac12 \) on \( [-1,1] \) and for the Epanechnikov
kernel, and verify the values \( \mu_2=1/5 \), \( R(K)=3/5 \) quoted
after @def-smo-kernel.
:::

::: {#exr-smo-efficiency}
[B2]

By @thm-smo-local-poly(d), the smallest achievable mean squared error is
proportional to \( \{\mu_2^2R(K)^4\}^{1/5} \). Compute this quantity for the
uniform, Epanechnikov and normal kernels (for the normal,
\( \mu_2=1 \) and \( R(K)=1/(2\sqrt{\pi}) \)), and express each as a ratio to the
Epanechnikov value. Conclude that the choice of kernel changes the attainable
error by a few per cent, whereas @eq-smo-hopt shows that getting \( h \) wrong by a
factor of two can change it by much more.
:::

::: {.solution}
Write \( C(K)=\{\mu_2^2R(K)^4\}^{1/5} \). Epanechnikov:
\( \{(1/5)^2(3/5)^4\}^{1/5}=(0.04\times0.1296)^{1/5}=0.3491 \). Uniform:
\( \{(1/3)^2(1/2)^4\}^{1/5}=(0.11111\times0.0625)^{1/5}=0.3701 \), a ratio of
\( 1.0602 \). Normal: \( R(K)=0.2821 \), so
\( (0.2821^4)^{1/5}=0.3633 \), a ratio of \( 1.0408 \). Six per cent and four per
cent. By contrast, replacing \( h^{*} \) by \( 2h^{*} \) in
\( Ah^4+B/(nh) \) multiplies the first term by 16 and divides the second by 2, and
with the terms in the ratio \( 1:4 \) at the optimum this raises the criterion by a
factor \( (16+2)/5=3.6 \).
:::

::: {#exr-smo-linear-exact}
[B3]

Suppose \( f \) is exactly linear, \( f(x)=\alpha+\beta x \). Show that the local
linear estimator is unbiased at every \( x \) and for every bandwidth, and that
the Nadaraya–Watson estimator is unbiased at \( x \) only if
\( \sum_iK_h(x_i-x)(x_i-x)=0 \). Relate the second condition
to @eq-smo-design-bias.
:::

### C. Going deeper

::: {#exr-smo-higher-order}
[C1]

Suppose \( f \) is four times and \( g \) three times continuously
differentiable, and let \( K \) be symmetric with \( \mu_2=0 \) and
\( \mu_4\ne0 \) — a *fourth-order* kernel, which must take negative values, and
which violates assumption (B).

::: {.enumerate options="label=(\alph*)"}
1. Carry @eq-smo-moment-sum two terms further, to
   \( S_{n,j}=nh^{j}\sum_{r=0}^{3}h^{r}g^{(r)}(x)\mu_{j+r}/r!+o(nh^{j+3}) \), and
   redo the local constant (\( q=0 \)) bias calculation in the proof
   of @thm-smo-local-poly to get
   \[
   \E\{\hat f(x)\}-f(x)=\mu_4h^4\Bigl\{\frac{f'g'''}{6g}+\frac{f''g''}{4g}
   +\frac{f'''g'}{6g}+\frac{f^{(4)}}{24}\Bigr\}(x)+o(h^4) .
   \]
   Deduce that the optimal bandwidth is of order \( n^{-1/9} \) and the mean
   squared error of order \( n^{-8/9} \).

2. Show that for \( q=1 \) the limiting moment matrix is
   \( \A=\diag(\mu_0,\mu_2)=\diag(1,0) \), which is singular, so the expansion in
   the proof does not apply. Writing \( a=S_{n,0} \), \( b=h^{-1}S_{n,1} \),
   \( c=h^{-2}S_{n,2} \) and \( u_i=(x_i-x)/h \), show from @eq-smo-equivalent that
   \( \ell_i(x)=K_h(x_i-x)(c-bu_i)/(ac-b^2) \), that \( b/c=O(h) \), and hence that
   the local linear weights collapse to the local constant ones.

3. Why is a fourth-order kernel rarely worth using? Consider the sign of
   \( \hat f \) when every \( y_i\ge0 \), the size of the constants, and (b).
:::
:::

::: {.solution}
(a) Symmetry gives \( \mu_1=\mu_3=\mu_5=0 \), and \( \mu_2=0 \) by assumption, so
the expansion leaves \( S_{n,0}=ng(x)\{1+O(h^4)\} \) and, to relative order
\( o(1) \), \( S_{n,1}=nh^4\mu_4g'''(x)/6 \),
\( S_{n,2}=nh^4\mu_4g''(x)/2 \), \( S_{n,3}=nh^4\mu_4g'(x) \) and
\( S_{n,4}=nh^4\mu_4g(x) \). Substituting these in
\( \sum_{j=1}^{4}f^{(j)}(x)S_{n,j}/(j!\,S_{n,0}) \) gives the display: every term
of order \( h \), \( h^2 \) and \( h^3 \) has vanished. Balancing the squared bias
\( h^8 \) against the variance \( (nh)^{-1} \) gives \( h\asymp n^{-1/9} \) and
mean squared error of order \( n^{-8/9} \).

(b) \( \A=(\mu_{j+k})_{j,k=0}^{1} \) has entries \( \mu_0=1 \), \( \mu_1=0 \),
\( \mu_2=0 \); it is singular, and the step of the proof that inverts it fails.
For the finite-sample weights, \( \bU\T\W\bU=\begin{pmatrix}a&b\\b&c\end{pmatrix} \)
in the notation of the proof, whose inverse has first row
\( (c,-b)/(ac-b^2) \), and \( (\bU\T\W)_{ji}=K_h(x_i-x)u_i^{j} \); multiplying gives
the stated \( \ell_i(x) \). By (a), \( b=nh^3\mu_4g'''/6 \) and
\( c=nh^2\mu_4g''/2 \) to leading order, so \( b/c=O(h) \) and
\( \ell_i(x)=K_h(x_i-x)/a+O(h) \), the Nadaraya–Watson weights. In floating point
the matrix is worse than singular: it is ill conditioned with an entry of
indefinite sign.

(c) The constant now involves \( f^{(4)} \), far harder to estimate than
\( f'' \); the gain over \( n^{-4/5} \) shows only at sample sizes where the two
rates differ appreciably; negative weights can make \( \hat f \) negative where
the data are not, which is unacceptable for a variance or a probability; and by
(b) the fit loses the boundary and design-bias advantages that were the reason
for \( q=1 \).
:::

::: {#exr-smo-random-design}
[C2]

Replace assumption (A) by: \( x_1,\dots,x_n \) are independent with density
\( g \), continuous and positive at \( x \). Show that
\( \E\{S_{n,j}\}=nh^j\{g(x)\mu_j+O(h)\} \) and
\( \Var(S_{n,j})=O(nh^{2j-1}) \), so that
\( (nh^j)^{-1}S_{n,j}\to g(x)\mu_j \) in probability. Explain why this gives the
conditional bias and variance of @thm-smo-local-poly with \( o_p \) in place of
\( o \), and why the argument needs more care than the fixed-design one when the
bias term of order \( h^2 \) is required (compare the size of the stochastic error
\( (nh^3)^{-1/2} \) with \( h \) at \( h\asymp n^{-1/5} \)).
:::

::: {.solution}
\( \E\{S_{n,j}\}=n\E\{K_h(x_1-x)(x_1-x)^j\}
=nh^j\int K(u)u^jg(x+hu)\,du=nh^j\{g(x)\mu_j+O(h)\} \) by continuity of \( g \),
and \( \Var(S_{n,j})\le n\E\{K_h(x_1-x)^2(x_1-x)^{2j}\}=O(nh^{2j-1}) \). The
standard deviation of \( (nh^j)^{-1}S_{n,j} \) is therefore \( O\{(nh)^{-1/2}\} \),
which tends to zero under (D), and Slutsky's theorem transfers every limit in the
proof. The difficulty is that the stochastic error in the bias is of order
\( h^2(nh)^{-1/2} \), which at \( h\asymp n^{-1/5} \) is \( n^{-4/5} \), the same
order as the squared bias. A statement about the *conditional* bias given the
design therefore needs the sharper expansions of Fan and Gijbels (1996, ch. 3).
:::

