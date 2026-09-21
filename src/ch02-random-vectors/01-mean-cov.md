# Mean vectors and covariance matrices

Random vectors enter regression in two different ways. The first is the response
vector. In the model
\[
\Y=\X\bbeta+\be,
\]
\( \Y=(Y_1,\dots,Y_n)\T \) has one entry per observational unit, and so does the
unobservable error \( \be \). The assumptions that do the real work in
[Chapter 5](../ch05-model-and-least-squares/index.html), [Chapter 6](../ch06-projections/index.html) and [Chapter 7](../ch07-optimality/index.html) are not about the shape
of any distribution. They are statements about moments: \( \E(\be)=\bzero \), and the errors have
equal variances and are uncorrelated. The second way is the measurement vector. When
several quantities are recorded on the same unit, such as income, age and years of
schooling, the unit contributes a vector whose components are usually correlated and
measured in different units. Random regressors, errors-in-variables models and
multivariate responses are all built from such vectors.

Both uses need the same toolkit: expectation and covariance as operations on vectors and
matrices, their behaviour under linear maps, and the expected values of quadratic forms.
Throughout, a random vector *has finite second moments* if every component has finite
variance.

## Random matrices and their expectation

::: {#def-rv-expectation}
[Expectation of a random matrix]

A **random matrix** \( \mathbf{Z}=(Z_{ij}) \) is an \( m\times k \) array of random variables. If
every \( Z_{ij} \) has a finite expectation, the **expectation** of \( \mathbf{Z} \) is the
\( m\times k \) matrix of expectations,
\[
\E(\mathbf{Z})=\bigl(\E(Z_{ij})\bigr).
\]
A random vector is the case \( k=1 \).
:::

The definition works entry by entry, so every rule that is linear in the entries carries
over at once.

::: {#prp-rv-expectation-linear}
[Linearity of expectation]

Let \( \mathbf{Z} \) and \( \mathbf{Z}_2 \) be \( m\times k \) random matrices with integrable entries, and let
\( \A \) (\( l\times m \)), \( \B \) (\( k\times s \)) and \( \mathbf{C} \) (\( l\times s \)) be constant matrices. Then

::: {.enumerate options="label=(\alph*)"}
1. \( \E(\A\mathbf{Z}\B+\mathbf{C})=\A\,\E(\mathbf{Z})\,\B+\mathbf{C} \);

2. \( \E(\mathbf{Z}+\mathbf{Z}_2)=\E(\mathbf{Z})+\E(\mathbf{Z}_2) \) and \( \E(\mathbf{Z}\T)=\E(\mathbf{Z})\T \);

3. if \( m=k \), then \( \E\tr(\mathbf{Z})=\tr\E(\mathbf{Z}) \).
:::

:::

::: {.proof}
The \( (i,j) \) entry of \( \A\mathbf{Z}\B+\mathbf{C} \) is
\( \sum_{r}\sum_{t}a_{ir}Z_{rt}b_{tj}+c_{ij} \), a finite linear combination of integrable
random variables with constant coefficients. Its expectation is
\( \sum_r\sum_t a_{ir}\E(Z_{rt})b_{tj}+c_{ij} \), which is the \( (i,j) \) entry of
\( \A\E(\mathbf{Z})\B+\mathbf{C} \). Parts (b) and (c) are the same argument applied to \( Z_{ij}+Z_{2,ij} \),
to \( Z_{ji} \), and to \( \sum_iZ_{ii} \).
:::

Two warnings go with this proposition. The multipliers must be constant: if \( \A \) is random,
\( \E(\A\mathbf{Z}) \) is generally not \( \E(\A)\E(\mathbf{Z}) \). And the expectation of a product of
random matrices is generally not the product of expectations. It is when every entry of
the first is independent of every entry of the second, because then
\( \E(Z_{ir}W_{rj})=\E(Z_{ir})\E(W_{rj}) \) for each term.

## Mean vector and covariance matrix

::: {#def-rv-mean-cov}
[Mean vector and covariance matrices]

Let \( \Y \) be a \( p\times1 \) random vector with finite second moments. Its **mean vector**
is \( \bmu=\E(\Y) \), and its **covariance matrix** is the \( p\times p \) matrix
\[
\Cov(\Y)=\E\bigl[(\Y-\bmu)(\Y-\bmu)\T\bigr],
\]
whose \( (i,j) \) entry is \( \sigma_{ij}=\Cov(Y_i,Y_j) \) and whose diagonal holds the
variances \( \sigma_{ii}=\Var(Y_i) \). If \( \bU \) (\( p\times1 \)) and \( \V \) (\( q\times1 \)) both have
finite second moments, their **cross-covariance matrix** is the \( p\times q \) matrix
\[
\Cov(\bU,\V)=\E\bigl[(\bU-\E\bU)(\V-\E\V)\T\bigr],
\]
with \( (i,j) \) entry \( \Cov(U_i,V_j) \).
:::

The expectations exist: \( \lvert U_iV_j\rvert\le\tfrac12(U_i^2+V_j^2) \), so finite second
moments make every product integrable. We often write \( \bSigma \) for \( \Cov(\Y) \), and
\( \Cov(\Y)=\Cov(\Y,\Y) \). Some books write \( \Var(\Y) \) or \( D(\Y) \) for the covariance matrix.
We keep \( \Var \) for scalars. The symbol \( \Cov \) then has three uses: \( \Cov(Y_i,Y_j) \) is a
number, \( \Cov(\Y) \) is a symmetric square matrix, and \( \Cov(\bU,\V) \) is a rectangular
matrix that is usually not symmetric even when \( p=q \).

::: {#prp-rv-cov-basic}
[Elementary properties]

Let \( \Y \), \( \bU \), \( \V \) have finite second moments, with \( \bmu=\E(\Y) \) and
\( \bSigma=\Cov(\Y) \), and let \( \mathbf{a} \), \( \bb \) be constant vectors of matching sizes. Then:

::: {.enumerate options="label=(\alph*)"}
1. \( \bSigma \) is symmetric;

2. \( \Cov(\bU,\V)=\E(\bU\V\T)-\E(\bU)\E(\V)\T \), and in particular
           \( \bSigma=\E(\Y\Y\T)-\bmu\bmu\T \);

3. \( \Cov(\V,\bU)=\Cov(\bU,\V)\T \);

4. \( \Cov(\bU+\mathbf{a},\V+\bb)=\Cov(\bU,\V) \);

5. \( \E\bigl[(\Y-\mathbf{a})(\Y-\mathbf{a})\T\bigr]=\bSigma+(\bmu-\mathbf{a})(\bmu-\mathbf{a})\T \), and hence
           \( \E\norm{\Y-\mathbf{a}}^2=\tr(\bSigma)+\norm{\bmu-\mathbf{a}}^2 \).
:::

:::

::: {.proof}
(a) and (c) hold because \( (\mathbf{u}\mathbf{v}\T)\T=\mathbf{v}\mathbf{u}\T \) and expectation commutes with
transposition. For (b), expand
\( (\bU-\bmu_U)(\V-\bmu_V)\T=\bU\V\T-\bU\bmu_V\T-\bmu_U\V\T+\bmu_U\bmu_V\T \) and take
expectations with @prp-rv-expectation-linear. The two middle terms each give
\( -\bmu_U\bmu_V\T \). (d) holds because adding a constant does not change the centred
vectors. For (e), write \( \Y-\mathbf{a}=(\Y-\bmu)+(\bmu-\mathbf{a}) \) and expand. The cross terms
contain the factor \( \E(\Y-\bmu)=\bzero \) and vanish. The trace identity follows from
\( \norm{\mathbf{v}}^2=\tr(\mathbf{v}\mathbf{v}\T) \) and @prp-rv-expectation-linear(c).
:::

Part (e) is the vector form of “mean squared error equals variance plus squared
bias”. It shows that the mean is the constant vector closest to \( \Y \) in expected squared
distance, and that the minimal value is \( \tr(\bSigma)=\sum_i\Var(Y_i) \). The matrix statement is stronger than the trace statement.
It says that the second-moment matrix about any point exceeds \( \bSigma \) by a
nonnegative definite matrix of rank at most one.

## Uncorrelated is not independent

Random vectors \( \bU \) and \( \V \) are **uncorrelated** if \( \Cov(\bU,\V)=\bzero \). If \( \bU \) and
\( \V \) are independent, then \( \E(U_iV_j)=\E(U_i)\E(V_j) \) for all \( i,j \), so they are
uncorrelated. The converse is false, and it fails in the most innocent-looking cases.

::: {#exm-rv-circle}
[A point on a circle]

Let \( \Theta \) be uniform on \( [0,2\pi) \) and \( \Y=(\cos\Theta,\sin\Theta)\T \). Then
\( \E(\cos\Theta)=\E(\sin\Theta)=0 \), \( \E(\cos^2\Theta)=\E(\sin^2\Theta)=\tfrac12 \), and
\( \E(\cos\Theta\sin\Theta)=\tfrac12\E(\sin2\Theta)=0 \). So \( \bmu=\bzero \) and
\( \Cov(\Y)=\tfrac12\I_2 \), exactly the moments of two independent coordinates. Yet
\( Y_1^2+Y_2^2=1 \): once \( Y_1 \) is known, \( Y_2 \) is known up to sign. Covariance measures
*linear* association only, and the dependence here is entirely nonlinear.
:::

For jointly normal vectors, zero covariance and independence coincide
(@thm-mvn-independence). That equivalence rests on joint normality, not on normal
marginals, and it is one of the main reasons normal theory is so tractable.

**Scalar summaries.**  The **total variance** \( \tr(\bSigma) \) ignores correlation; the
**generalized variance** \( \det(\bSigma) \) does not. If \( \lambda_1,\dots,\lambda_p \) are the eigenvalues of
\( \bSigma \) (@thm-mat-spectral), these are \( \sum_i\lambda_i \) and \( \prod_i\lambda_i \), and both are unchanged
by a rotation \( \Y\mapsto\Q\Y \). The determinant can be tiny even when every variance is large, if one
eigenvalue is near zero. @thm-rv-cov-nnd shows that a zero eigenvalue means some linear
combination of \( \Y \) does not vary at all.

## Exercises

### A. Check your understanding

::: {#exr-rv-expected-gram}
[A1]

Let \( \mathbf{Z} \) be an \( n\times p \) random matrix whose entries are independent, with mean \( 0 \) and
variance \( \sigma^2 \). Show that \( \E(\mathbf{Z}\T\mathbf{Z})=n\sigma^2\I_p \) and \( \E(\mathbf{Z}\mathbf{Z}\T)=p\sigma^2\I_n \). For a
constant \( p\times p \) matrix \( \B \), find \( \E(\mathbf{Z}\B\mathbf{Z}\T) \).
:::

### B. Practice

::: {#exr-rv-sign-product}
[B1]

Let \( X \) have a distribution symmetric about \( 0 \) with \( \Var(X)=1 \) and \( \E X^4<\infty \), let \( S \) be
independent of \( X \) with \( \Pr(S=1)=\Pr(S=-1)=\tfrac12 \), and put \( \Y=(X,SX)\T \). Show that the two
components have the same distribution and are uncorrelated, and that
\( \Cov(Y_1^2,Y_2^2)=\Var(X^2) \). Deduce that \( Y_1 \) and \( Y_2 \) are independent iff \( X \) is itself a
random sign. ([Chapter 3](../ch03-multivariate-normal/index.html) takes \( X \) standard normal.)
:::

::: {.solution}
Since \( S \) is independent of \( X \) and \( X \) is symmetric,
\( \Pr(SX\le y)=\tfrac12\Pr(X\le y)+\tfrac12\Pr(-X\le y)=\Pr(X\le y) \). By independence,
\( \Cov(X,SX)=\E(S)\E(X^2)=0 \). Both squares equal \( X^2 \), so \( \Cov(Y_1^2,Y_2^2)=\Var(X^2) \). Functions of
independent variables are independent, hence uncorrelated, so \( \Var(X^2)>0 \) makes \( Y_1 \) and \( Y_2 \)
dependent. If \( \Var(X^2)=0 \), then \( X^2=\E X^2=1 \) with probability one, and by symmetry \( X=\pm1 \)
with probability \( \tfrac12 \) each. Then
\( \Pr(Y_1=a,Y_2=b)=\Pr(X=a)\Pr(S=ab)=\tfrac14 \) for all signs \( a,b \), which is the product of the
marginal probabilities, as in @exm-rv-pairwise.
:::

::: {#exr-rv-best-constant-matrix}
[B2]

Show that \( \E\bigl[(\Y-\mathbf{a})(\Y-\mathbf{a})\T\bigr]-\Cov(\Y) \) is nonnegative definite for every constant
\( \mathbf{a} \), with equality iff \( \mathbf{a}=\bmu \). Deduce that \( \mathbf{a}=\bmu \) minimizes
\( \E\bigl[(\Y-\mathbf{a})\T\mathbf{N}(\Y-\mathbf{a})\bigr] \) for every positive definite \( \mathbf{N} \).
:::

### C. Going deeper

::: {#exr-rv-moment-support}
[C1]

Let \( \bmu\in\Real^p \) and let \( \bSigma \) be nonnegative definite with spectral decomposition
\( \bSigma=\sum_{i=1}^p\lambda_i\boldsymbol{\gamma}_i\boldsymbol{\gamma}_i\T \) (@thm-mat-spectral).

::: {.enumerate options="label=(\alph*)"}
1. Let \( \Y \) take the \( 2p \) values \( \bmu\pm\sqrt{p\lambda_i}\,\boldsymbol{\gamma}_i \), each with probability
           \( 1/(2p) \). Show that \( \E(\Y)=\bmu \) and \( \Cov(\Y)=\bSigma \). So a mean vector and a covariance
           matrix say nothing at all about the support: some vector with these moments lives on
           \( 2p \) points, while, for positive definite \( \bSigma \), another has a density that is
           positive everywhere (@thm-rv-cov-nnd(b)).

2. Show that if a random vector takes only \( k \) values then \( \rank\Cov(\Y)\le k-1 \). Hence at
           least \( p+1 \) atoms are needed when \( \bSigma \) is positive definite.

3. Show that \( p+1 \) atoms always suffice. *Hint:* let \( \A \) be the last \( p \) rows of the
           \( (p+1)\times(p+1) \) Helmert matrix of @exm-mat-helmert, so that \( \A\A\T=\I_p \) and
           \( \A\bone=\bzero \), and put the mass \( 1/(p+1) \) on each column of \( \sqrt{p+1}\,\A \).
:::

:::

::: {.solution}
(a) The \( 2p \) values come in pairs symmetric about \( \bmu \), so
\( \E(\Y)=\bmu \). Each pair contributes
\( 2\cdot\frac{1}{2p}\,p\lambda_i\boldsymbol{\gamma}_i\boldsymbol{\gamma}_i\T=\lambda_i\boldsymbol{\gamma}_i\boldsymbol{\gamma}_i\T \) to
\( \E[(\Y-\bmu)(\Y-\bmu)\T] \), and the sum over \( i \) is \( \bSigma \). (Atoms with \( \lambda_i=0 \) coincide with
\( \bmu \) and may be merged.)

(b) If \( \Y \) takes the values \( \y_1,\dots,\y_k \), then \( \Y-\bmu \) lies with probability one in the
span of \( \y_1-\bmu,\dots,\y_k-\bmu \), which has dimension at most \( k-1 \) because these \( k \) vectors
sum to \( \bzero \) after weighting by their probabilities. By @thm-rv-cov-nnd(d),
\( \C(\bSigma) \) is contained in that span, so \( \rank(\bSigma)\le k-1 \). A positive definite \( \bSigma \)
has rank \( p \), forcing \( k\ge p+1 \).

(c) Let \( \bv_i=\sqrt{p+1}\,\A\mathbf{e}_i \) for \( i=1,\dots,p+1 \) and let \( \Y=\bmu+\bSigma^{1/2}\bv_I \),
where \( I \) is uniform on \( \{1,\dots,p+1\} \). Then
\( \sum_i\bv_i=\sqrt{p+1}\,\A\bone=\bzero \), so \( \E(\Y)=\bmu \), and
\[
\frac{1}{p+1}\sum_{i=1}^{p+1}\bv_i\bv_i\T=\frac{1}{p+1}(p+1)\,\A\Bigl(\sum_i\mathbf{e}_i\mathbf{e}_i\T\Bigr)\A\T=\A\A\T=\I_p ,
\]
so \( \Cov(\Y)=\bSigma^{1/2}\I_p\bSigma^{1/2}=\bSigma \) by @thm-rv-linear(b). The \( p+1 \) atoms are the
vertices of a regular simplex centred at \( \bmu \).
:::

::: {#exr-rv-isotropic}
[C2]

Call \( \Y \) *isotropic in law* if \( \mathbf{a}\T\Y \) has the same distribution for every unit vector
\( \mathbf{a} \).

::: {.enumerate options="label=(\alph*)"}
1. Show that an isotropic \( \Y \) with finite second moments has \( \bmu=\bzero \) and
           \( \bSigma=c\,\I_p \) with \( c=\Var(\mathbf{a}\T\Y)\ge0 \). *Hint:* a symmetric matrix whose quadratic
           form vanishes identically is \( \bzero \).

2. Show that the converse fails: for \( \Y \) uniform on the four points \( (\pm\sqrt2,0)\T \),
           \( (0,\pm\sqrt2)\T \) we have \( \Cov(\Y)=\I_2 \), yet \( \mathbf{e}_1\T\Y \) and
           \( (\mathbf{e}_1+\mathbf{e}_2)\T\Y/\sqrt2 \) have different distributions.

3. Check that the vector of @exm-rv-circle *is* isotropic, and deduce that isotropy in law
           does not identify a distribution either.
:::

:::

::: {.solution}
(a) For a unit \( \mathbf{a} \), \( \E(\mathbf{a}\T\Y)=\mathbf{a}\T\bmu \) is the same number for
\( \mathbf{a} \) and \( -\mathbf{a} \), so \( \mathbf{a}\T\bmu=-\mathbf{a}\T\bmu=0 \) for every \( \mathbf{a} \), giving \( \bmu=\bzero \).
Likewise \( \mathbf{a}\T\bSigma\mathbf{a}=c \) for every unit \( \mathbf{a} \), and by homogeneity
\( \mathbf{a}\T\bSigma\mathbf{a}=c\norm{\mathbf{a}}^2 \) for every \( \mathbf{a} \). The symmetric matrix \( \bSigma-c\I \) then has
\( \mathbf{a}\T(\bSigma-c\I)\mathbf{a}=0 \) for all \( \mathbf{a} \); polarization,
\( 2\mathbf{a}\T(\bSigma-c\I)\bb=(\mathbf{a}+\bb)\T(\bSigma-c\I)(\mathbf{a}+\bb)-\mathbf{a}\T(\bSigma-c\I)\mathbf{a}-\bb\T(\bSigma-c\I)\bb \),
makes every entry zero.

(b) The four points are symmetric about the origin, so \( \bmu=\bzero \), and
\( \E(\Y\Y\T)=\tfrac12\cdot2\,\mathbf{e}_1\mathbf{e}_1\T+\tfrac12\cdot2\,\mathbf{e}_2\mathbf{e}_2\T=\I_2 \). But
\( Y_1 \) takes the values \( \pm\sqrt2 \) with probability \( \tfrac14 \) each and \( 0 \) with probability
\( \tfrac12 \), whereas \( (Y_1+Y_2)/\sqrt2 \) takes the values \( \pm1 \) with probability \( \tfrac12 \) each.
Equal covariance matrices, different projections.

(c) There \( \mathbf{a}\T\Y=\cos(\Theta-\psi) \) for the angle \( \psi \) of \( \mathbf{a} \), and \( \Theta-\psi \) is
again uniform on a full period, so the law of \( \mathbf{a}\T\Y \) does not depend on \( \mathbf{a} \); its
variance is \( \tfrac12 \), matching \( \Cov(\Y)=\tfrac12\I_2 \). A spherical normal vector with
covariance \( \tfrac12\I_2 \) is isotropic too, and the two laws are different.
:::
