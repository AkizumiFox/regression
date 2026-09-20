# Singular covariance matrices

A covariance matrix is nonnegative definite, not necessarily positive definite.
When it is singular, some linear combination of the observations has zero
variance and is therefore a *constant*: the data satisfy an exact linear
relation. This is not a pathology to be regularized away — it happens whenever a linear
operation has destroyed information — and it changes the theory, because some
functions of \( \bbeta \) become known without error.

## Where a singular covariance comes from

Three mechanisms cover most cases.

**Constraints built into the data.** Percentages adjusted to add to a hundred,
effects reported under a sum-to-zero side condition
([Chapter 8](../ch08-estimability/index.html)), residuals from an earlier fit: the
reported vector satisfies an exact linear identity, and if \( \mathbf{a}\T\Y \) is
constant then \( \V\mathbf{a}=\bzero \) by @thm-rv-cov-nnd(c).

**A deterministic component.** Some observations may be free of measurement
error — a boundary value known exactly, a fitted value substituted for a missing
one — and the corresponding rows and columns of \( \V \) vanish. So is
\( \V=\Z\G\Z\T \), the covariance of a random-effects model \( \Y=\X\bbeta+\Z\bu \)
with no separate measurement error, as soon as \( \Z \) has fewer columns than
rows; [Chapter 32](../ch32-linear-mixed-models/index.html) studies that family.

**Linear filtering with redundancy.** Replacing uncorrelated observations by more
linear combinations than the rank of the transformation — differences together
with their total, overlapping averages together with the underlying ones —
produces a singular covariance. Redundancy is what matters, not filtering:
\( n-1 \) successive differences of \( n \) uncorrelated observations have a
nonsingular tridiagonal covariance.

## The consistency condition

::: {.idea}
A singular \( \V \) confines the *errors* to a subspace: \( \Pr(\be\in\C(\V))=1 \).
With \( \E(\Y)=\X\bbeta\in\C(\X) \), this confines the data to
\( \C(\X:\V)=\C(\X)+\C(\V) \) and \( \X\bbeta \) to the affine set \( \y+\C(\V) \).
Estimation becomes a constrained problem: find the point of \( \C(\X) \) nearest
\( \y \) *among those the data allow*.
:::

Write \( \C(\X:\V) \) for \( \C(\X)+\C(\V) \), the column space of \( [\X\ \V] \). The
following small lemma is used repeatedly.

::: {#lem-ggm-colspace-sum}
If \( \A \) and \( \B \) are nonnegative definite of the same order, then
\( \C(\A+\B)=\C(\A)+\C(\B) \).
:::

::: {.proof}
\( (\A+\B)\bw=\A\bw+\B\bw \) shows \( \C(\A+\B)\subseteq\C(\A)+\C(\B) \). For the
reverse, take \( \bw\perp\C(\A+\B) \). Then \( \bw\T\A\bw+\bw\T\B\bw=0 \), and both terms
are nonnegative, so both vanish. By @thm-mat-square-root,
\( \norm{\A^{1/2}\bw}^2=\bw\T\A\bw=0 \), so \( \A\bw=\A^{1/2}\A^{1/2}\bw=\bzero \), and
likewise \( \B\bw=\bzero \). Hence \( \bw\perp\C(\A) \) and \( \bw\perp\C(\B) \), so
\( \C(\A+\B)\perpc\subseteq\bigl(\C(\A)+\C(\B)\bigr)\perpc \), which gives the
inclusion.
:::

## The general theorem

Fix an orthogonal splitting adapted to \( \V \): let \( \boldsymbol{\Gamma} \) be \( n\times q \) and
\( \mathbf{N} \) be \( n\times(n-q) \) with orthonormal columns, \( \C(\boldsymbol{\Gamma})=\C(\V) \) and
\( \C(\mathbf{N})=\C(\V)\perpc \), where \( q=\rank(\V) \). Put \( \bLambda=\boldsymbol{\Gamma}\T\V\boldsymbol{\Gamma} \), which is
\( q\times q \) positive definite, and let \( \V^{+} \) be the Moore–Penrose inverse (@def-mat-moore-penrose),
so that \( \V^{+}=\boldsymbol{\Gamma}\bLambda^{-1}\boldsymbol{\Gamma}\T \). Write
\( r=\rank(\X) \) and \( d=\dim\bigl(\C(\X)\cap\C(\V)\bigr) \), and let \( \bU \) be any
\( n\times d \) matrix of full column rank with \( \C(\bU)=\C(\X)\cap\C(\V) \).

::: {#thm-ggm-singular}
[The general Gauss–Markov model]

Assume @eq-ggm-model with \( \V \) nonnegative definite of rank \( q \).

::: {.enumerate options="label=(\alph*)"}
1. *(Consistency.)* \( \mathbf{N}\T\Y=\mathbf{N}\T\X\bbeta \) with probability one; consequently
   \( \Pr\bigl(\Y\in\C(\X:\V)\bigr)=1 \) and \( \Pr\bigl(\X\bbeta\in\Y+\C(\V)\bigr)=1 \).

2. *(Estimability and perfect estimation.)* \( \blambda\T\bbeta \) has a linear
   unbiased estimator iff \( \blambda\in\C(\X\T) \), exactly as when \( \V=\I \). It has a
   linear unbiased estimator of variance zero iff \( \blambda\in\C(\X\T\mathbf{N}) \), a
   subspace of dimension \( r-d \), and then \( \blambda\T\bbeta=\bb\T\mathbf{N}\T\Y \) with
   probability one for any \( \bb \) with \( \X\T\mathbf{N}\bb=\blambda \).

3. *(Rao's unified formula.)* Let \( \bT=\V+\X\X\T \). Then \( \C(\bT)=\C(\X:\V) \),
   \( \bT\bT^{+}\X=\X \), and
   \[
\A=\X(\X\T\bT^{+}\X)\ginv\X\T\bT^{+}
\]{#eq-ggm-rao}

   satisfies \( \A\X=\X \) and \( \A\V(\I-\M)=\mathbf{0} \). Hence \( \A\Y \) is a BLUE of
   \( \X\bbeta \).

4. *(Uniqueness.)* Every estimable \( \blambda\T\bbeta \) has a BLUE, and any two
   BLUEs of it agree for every \( \y\in\C(\X:\V) \), hence with probability one.

5. *(Projection form.)* Let \( \y\in\C(\X:\V) \). The set
   \( \mathcal F(\y)=\{\bb:\y-\X\bb\in\C(\V)\} \) is nonempty, and
   \[
\tilde{\bbeta}\in\argmin_{\bb\in\mathcal F(\y)}\ (\y-\X\bb)\T\V^{+}(\y-\X\bb)
\]{#eq-ggm-projection-form}

   determines \( \X\tilde{\bbeta} \) uniquely; the map \( \y\mapsto\X\tilde{\bbeta} \) is
   linear and \( \X\tilde{\bbeta} \) is a BLUE of \( \X\bbeta \). Explicitly, with
   \( \mathbf{K}=\mathbf{N}\T\X \) and any generalized inverse \( \mathbf{K}\ginv \),
   \[
\X\tilde{\bbeta}=\X\bb_*+\bU(\bU\T\V^{+}\bU)^{-1}\bU\T\V^{+}(\y-\X\bb_*),
\qquad \bb_*=\mathbf{K}\ginv\mathbf{N}\T\y .
\]{#eq-ggm-two-step}

6. *(Error sum of squares.)* Put \( \he=\Y-\X\tilde{\bbeta} \) and
   \( \text{SSE}_{\V}=\he\T\V^{+}\he \). Then \( \E(\text{SSE}_{\V})=\sigma^2(q-d) \),
   and \( q-d=\rank(\X:\V)-\rank(\X) \). If \( \Y \) is normal, then
   \( \text{SSE}_{\V}/\sigma^2\sim\chi^2(q-d) \), independently of
   \( \X\tilde{\bbeta} \).
:::
:::

::: {.proof}
(a) By @thm-rv-cov-nnd(d), \( \Pr(\be\in\C(\sigma^2\V))=1 \), and \( \C(\V)=\C(\mathbf{N})\perpc \),
so \( \mathbf{N}\T\be=\bzero \) almost surely, giving \( \mathbf{N}\T\Y=\mathbf{N}\T\X\bbeta \). Then
\( \Y=\X\bbeta+\be\in\C(\X)+\C(\V) \) and \( \X\bbeta=\Y-\be\in\Y+\C(\V) \).

(b) \( \mathbf{c}\T\Y \) is unbiased for \( \blambda\T\bbeta \) for every \( \bbeta \) iff
\( \X\T\mathbf{c}=\blambda \), which is possible iff \( \blambda\in\C(\X\T) \). Its variance is
\( \sigma^2\mathbf{c}\T\V\mathbf{c} \), which vanishes iff \( \V\mathbf{c}=\bzero \) (as in @lem-ggm-colspace-sum),
that is iff \( \mathbf{c}\in\C(\V)\perpc=\C(\mathbf{N}) \), say
\( \mathbf{c}=\mathbf{N}\bb \). Then \( \blambda=\X\T\mathbf{N}\bb \) and \( \mathbf{c}\T\Y=\bb\T\mathbf{N}\T\Y \), which
equals \( \blambda\T\bbeta \) almost surely by (a). For the dimension,
\( \rank(\X\T\mathbf{N})=\rank(\mathbf{N})-\dim\{\bb:\mathbf{N}\bb\perp\C(\X)\} \) by
@thm-mat-rank-nullity, and \( \{\mathbf{N}\bb:\mathbf{N}\bb\perp\C(\X)\}=\C(\V)\perpc\cap\C(\X)\perpc
=\bigl(\C(\X:\V)\bigr)\perpc \), of dimension \( n-\rank(\X:\V) \). Since
\( \rank(\X:\V)=r+q-d \), we get \( \rank(\X\T\mathbf{N})=(n-q)-(n-r-q+d)=r-d \).

(c) By @lem-ggm-colspace-sum and \( \C(\X\X\T)=\C(\X) \) (@cor-proj-gram-colspace),
\( \C(\bT)=\C(\V)+\C(\X)=\C(\X:\V) \). Since \( \C(\X)\subseteq\C(\bT) \), the system
\( \bT\mathbf{H}=\X \) is consistent (@thm-mat-consistency), so \( \X=\bT\mathbf{H} \) for some
\( \mathbf{H} \) and \( \bT\bT^{+}\X=\bT\bT^{+}\bT\mathbf{H}=\bT\mathbf{H}=\X \).

Write \( \bT^{+}=\bL\bL\T \), possible because \( \bT^{+} \) is symmetric nonnegative
definite, and set \( \X_*=\bL\T\X \), so that \( \X\T\bT^{+}\X=\X_*\T\X_* \). If \( \X_*\bw=\bzero \)
then \( \X\bw\in\Null(\bL\T)=\Null(\bT^{+})=\C(\bT)\perpc \), while
\( \X\bw\in\C(\X)\subseteq\C(\bT) \); hence \( \X\bw=\bzero \). Every column of
\( \I-(\X_*\T\X_*)\ginv\X_*\T\X_* \) lies in \( \Null(\X_*\T\X_*)=\Null(\X_*) \), so
\( \X\bigl(\I-(\X_*\T\X_*)\ginv\X_*\T\X_*\bigr)=\mathbf{0} \), which is
\( \A\X=\X(\X_*\T\X_*)\ginv\X_*\T\bL\T\X=\X(\X_*\T\X_*)\ginv\X_*\T\X_*=\X \).

For the second condition, let \( \mathbf{d} \) satisfy \( \X\T\mathbf{d}=\bzero \). Using
\( \V=\bT-\X\X\T \) and the symmetry of \( \bT \) and \( \bT^{+} \),
\[
\X\T\bT^{+}\V\mathbf{d}=\X\T\bT^{+}\bT\mathbf{d}-\X\T\bT^{+}\X\X\T\mathbf{d}
=(\bT\bT^{+}\X)\T\mathbf{d}=\X\T\mathbf{d}=\bzero,
\]
so \( \A\V\mathbf{d}=\bzero \). As \( \mathbf{d} \) ranges over \( \C(\X)\perpc \) this says
\( \A\V(\I-\M)=\mathbf{0} \). Now @lem-ggm-blue-criterion(c) makes \( \A\Y \) a BLUE of
\( \X\bbeta \).

(d) Existence is (c): for estimable \( \blambda\T\bbeta=\boldsymbol{\uprho}\T\X\bbeta \), the statistic
\( \boldsymbol{\uprho}\T\A\Y \) is a BLUE. Uniqueness is @lem-ggm-blue-criterion(b).

(e) If \( \y=\X\bb_0+\bv \) with \( \bv\in\C(\V) \) then \( \bb_0\in\mathcal F(\y) \), so
\( \mathcal F(\y)\ne\emptyset \); and \( \bb\in\mathcal F(\y) \) iff \( \mathbf{K}\bb=\mathbf{N}\T\y \),
a consistent system, whose solutions are \( \mathbf{K}\ginv\mathbf{N}\T\y \) plus \( \Null(\mathbf{K}) \).
Since \( \X\Null(\mathbf{K})=\{\X\mathbf{h}:\mathbf{N}\T\X\mathbf{h}=\bzero\}=\C(\X)\cap\C(\V)=\C(\bU) \), the
fitted vectors available are \( \X\bb=\X\bb_*+\bU\mathbf{c} \) with
\( \bb_*=\mathbf{K}\ginv\mathbf{N}\T\y \) and \( \mathbf{c}\in\Real^d \) free. Also
\( \y-\X\bb_*\in\C(\V) \), and on \( \C(\V) \) the form \( \mathbf{a}\T\V^{+}\mathbf{a} \) is positive
definite; so minimizing over \( \mathbf{c} \) is an ordinary generalized least squares
problem in the \( \V^{+} \) inner product with model matrix \( \bU \) of full column
rank, whose unique solution is
\( \hat{\mathbf{c}}=(\bU\T\V^{+}\bU)^{-1}\bU\T\V^{+}(\y-\X\bb_*) \). This gives
@eq-ggm-two-step, shows that \( \X\tilde{\bbeta} \) is unique and that the map is
linear, say \( \X\tilde{\bbeta}=\A_*\y \).

It remains to check the two conditions of @lem-ggm-blue-criterion(c) for
\( \A_* \). If \( \y=\X\bbeta \) then \( \mathbf{K}\bb_*=\mathbf{N}\T\X\bbeta \) and
\( \mathbf{N}\T(\y-\X\bb_*)=\bzero \), so \( \y-\X\bb_*\in\C(\X)\cap\C(\V)=\C(\bU) \); the
\( \V^{+} \)-projection onto \( \C(\bU) \) leaves it fixed, and
\( \A_*\X\bbeta=\X\bb_*+(\y-\X\bb_*)=\X\bbeta \). Thus \( \A_*\X=\X \). If \( \y=\V\mathbf{d} \)
with \( \X\T\mathbf{d}=\bzero \), then \( \mathbf{N}\T\V\mathbf{d}=\bzero \) so we may take \( \bb_*=\bzero \),
and \( \bU\T\V^{+}\V\mathbf{d}=\bU\T\boldsymbol{\Gamma}\boldsymbol{\Gamma}\T\mathbf{d}=\bU\T\mathbf{d}=\bzero \), using
\( \V^{+}\V=\boldsymbol{\Gamma}\boldsymbol{\Gamma}\T \), \( \C(\bU)\subseteq\C(\V) \) and
\( \C(\bU)\subseteq\C(\X) \). So \( \A_*\V\mathbf{d}=\bzero \), that is
\( \A_*\V(\I-\M)=\mathbf{0} \).

(f) Put \( \Z_*=\bLambda^{-1/2}\boldsymbol{\Gamma}\T(\I-\X\mathbf{K}\ginv\mathbf{N}\T)\Y \), a random vector in
\( \Real^q \). Because \( \V\mathbf{N}=\mathbf{0} \),
\( (\I-\X\mathbf{K}\ginv\mathbf{N}\T)\V(\I-\X\mathbf{K}\ginv\mathbf{N}\T)\T=\V \), so
\( \Cov(\Z_*)=\sigma^2\bLambda^{-1/2}\boldsymbol{\Gamma}\T\V\boldsymbol{\Gamma}\bLambda^{-1/2}=\sigma^2\I_q \). Its
mean is \( \bLambda^{-1/2}\boldsymbol{\Gamma}\T\X(\I-\mathbf{K}\ginv\mathbf{K})\bbeta \), which lies in
\( \bLambda^{-1/2}\boldsymbol{\Gamma}\T\C(\bU) \), a subspace of dimension \( d \) (the map
\( \bLambda^{-1/2}\boldsymbol{\Gamma}\T \) is injective on \( \C(\V)\supseteq\C(\bU) \)). By the
computation in (e), \( \text{SSE}_{\V}=\norm{(\I-\M_d)\Z_*}^2 \) where \( \M_d \) is the
orthogonal projection of \( \Real^q \) onto that \( d \)-dimensional subspace, and
\( \X\tilde{\bbeta} \) is a fixed vector plus a function of \( \M_d\Z_* \), since
\( \mathbf{N}\T\Y \) is constant by (a). Now @thm-rv-quadform-mean gives
\( \E(\text{SSE}_{\V})=\sigma^2\tr(\I-\M_d)=\sigma^2(q-d) \) with no distributional
assumption, and under normality @thm-opt-sampling applied to \( \Z_* \) gives the
\( \chi^2(q-d) \) law and the independence. Finally \( \rank(\X:\V)=r+q-d \) by the dimension formula for a sum of
subspaces, so \( q-d=\rank(\X:\V)-r \).
:::

The singular model is the nonsingular model with a deterministic layer on top.
Part (b) isolates that layer, \( r-d \) linear functions of \( \bbeta \) known exactly
once the data are seen; part (e) pins down a particular solution with them and
then does generalized least squares inside \( \C(\X)\cap\C(\V) \), the part of the
model space the noise can still reach; part (f) counts what is left. If
\( \C(\X)\subseteq\C(\V) \) then \( d=r \), nothing is known exactly, and
@eq-ggm-two-step is generalized least squares with \( \V^{+} \) for \( \V^{-1} \) on
\( q-r \) degrees of freedom; if \( \V \) is nonsingular then \( q=n \), \( d=r \), and
everything collapses to @thm-ggm-inference.

::: {when-format="html"}
![**Figure 31.3.1.** The geometry of a singular model with \( n=3 \). The plane is
\( \C(\X) \); \( \C(\V) \) is horizontal, so \( \C(\X)\cap\C(\V) \) is the lower line. The
data \( \y \) determine \( \mathbf{N}\T\X\bbeta \) exactly, which confines \( \X\bbeta \) to the
dashed line \( \C(\X)\cap(\y+\C(\V)) \); the estimate is the point of that line nearest
\( \y \) in the \( \V^{+} \) geometry, and the residual lies in
\( \C(\V) \).](singular_geometry.svg){#fig-ggm-singular width=66%}
:::

::: {when-format="pdf"}
![The geometry of a singular model with \( n=3 \). The plane is
\( \C(\X) \); \( \C(\V) \) is horizontal, so \( \C(\X)\cap\C(\V) \) is the lower line. The
data \( \y \) determine \( \mathbf{N}\T\X\bbeta \) exactly, which confines \( \X\bbeta \) to the
dashed line \( \C(\X)\cap(\y+\C(\V)) \); the estimate is the point of that line nearest
\( \y \) in the \( \V^{+} \) geometry, and the residual lies in
\( \C(\V) \).](singular_geometry.pdf){width=66%}
:::

## A worked example

::: {#exm-ggm-effects}
[Regression on estimated treatment effects]

A balanced one-way layout with \( a \) treatments and \( m \) observations each is
analysed, and the estimated effects \( \hat\alpha_1,\dots,\hat\alpha_a \) are reported
under the side condition \( \sum_k\hat\alpha_k=0 \) (@thm-est-side-conditions).
A second analyst, with access only to these numbers,
regresses them on the dose \( x_k \) applied to treatment \( k \), fitting
\( \E(\hat\alpha_k)=\beta_0+\beta_1x_k \) with an intercept.

The reported vector has \( \Cov(\hat{\boldsymbol{\upalpha}})=\tau^2\V \) with
\( \V=\I_a-\bone\bone\T/a \) and \( \tau^2=\sigma^2/m \), because
\( \hat\alpha_k=\bar Y_k-\bar Y \) and the \( \bar Y_k \) are uncorrelated with
variance \( \tau^2 \). Here \( q=a-1 \), \( \C(\V)=\bone\perpc \), and
\( \mathbf{N}=\bone/\sqrt a \). With \( \bone\in\C(\X) \) we get \( d=1 \) and \( r=2 \), so
@thm-ggm-singular(b) gives exactly \( r-d=1 \) perfectly known function:
\( \X\T\mathbf{N}\propto(a,\sum_kx_k)\T \), so
\[
a\beta_0+\Bigl(\sum_kx_k\Bigr)\beta_1=\bone\T\X\bbeta=\bone\T\hat{\boldsymbol{\upalpha}}=0
\]
with probability one. In words: the intercept is determined by the slope,
\( \beta_0=-\bar x\beta_1 \), and the data say so without error. Fitting an
intercept was a mistake, and the model itself reveals it.

Everything else is unchanged. Because \( \V\X=\X-\bone(\bone\T\X)/a \) has columns
in \( \C(\X) \), @thm-ggm-ols-blue applies and ordinary least squares is already
best: \( \hat\beta_1 \) is the usual centred slope and
\( \hat\beta_0=-\bar x\hat\beta_1 \). Since \( \V \) is idempotent, \( \V^{+}=\V \) and
\( \text{SSE}_{\V}=\he\T\he \) is the usual residual sum of squares, on
\( q-d=a-2 \) degrees of freedom. What is *not* unchanged is the covariance
matrix. Because \( \bone\T\hat{\boldsymbol{\upalpha}}=0 \) exactly,
\[
\begin{aligned}
\Var(\hat\beta_1)&=\frac{\tau^2}{S_{xx}},\qquad
\Var(\hat\beta_0)=\bar x^2\,\Var(\hat\beta_1),\\
\Cov(\hat\beta_0,\hat\beta_1)&=-\bar x\,\Var(\hat\beta_1),
\end{aligned}
\]
whereas the ordinary formula \( s^2(\X\T\X)^{-1} \) reports
\( s^2(1/a+\bar x^2/S_{xx}) \) for the intercept. With \( a=5 \) doses
\( 0,1,2,4,8 \), \( m=6 \) replicates and \( \sigma=1.2 \), so \( \tau^2=0.24 \),
\( \bar x=3.0 \) and \( S_{xx}=40.0 \): the true variances are
\( 0.00600 \) for the slope and \( 0.0540 \) for the intercept, while the
ordinary formula reports \( 0.1020 \) on average for the intercept, too
large by a factor of \( 1.89 \). A simulation of \( 200{,}000 \) layouts
reproduces \( 0.00600 \) and \( 0.0540 \).
:::

```{.python .run #cell-singular-effects-setup}
import numpy as np

dose = np.array([0.0, 1.0, 2.0, 4.0, 8.0])           # a = 5 treatments
a, reps, sigma = len(dose), 6, 1.2
X = np.column_stack([np.ones(a), dose])
V = np.eye(a) - np.ones((a, a)) / a                  # Cov(effects) = (sigma^2/reps) V
xbar = dose.mean()
Sxx = float(((dose - xbar) ** 2).sum())

Vplus = V                                            # V is symmetric idempotent, so V+ = V
N = np.ones((a, 1)) / np.sqrt(a)                     # an orthonormal basis of C(V) perp
print("rank(V) =", np.linalg.matrix_rank(V), " rank([X V]) =",
      np.linalg.matrix_rank(np.column_stack([X, V])), " rank(X) =", np.linalg.matrix_rank(X))
print("the perfectly known function:", (X.T @ N).ravel())
```

```{.python .run #cell-singular-effects-simulate}
rng = np.random.default_rng(3107)
trials = 200_000
cell = rng.normal(scale=sigma / np.sqrt(reps), size=(trials, a))     # treatment means
cell += 0.35 * (dose - xbar)                                         # the true effects
effects = cell - cell.mean(axis=1, keepdims=True)                    # sum to zero exactly
B = effects @ np.linalg.solve(X.T @ X, X.T).T                        # ordinary least squares

tau2 = sigma**2 / reps
print("Var(slope):     %.5f simulated, %.5f from tau^2/Sxx" % (B[:, 1].var(), tau2 / Sxx))
print("Var(intercept): %.5f simulated, %.5f from xbar^2 tau^2/Sxx"
      % (B[:, 0].var(), xbar**2 * tau2 / Sxx))
print("what OLS reports for the intercept: %.5f" % (tau2 * (1 / a + xbar**2 / Sxx)))
```

The script behind @exm-ggm-effects also computes the estimate four ways — by
@eq-ggm-projection-form, by the closed form @eq-ggm-two-step, by @eq-ggm-rao with
\( \bT=\V+\X\X\T \), and by ordinary least squares — and checks that all four
agree, as @thm-ggm-singular(d) requires, and that @eq-ggm-two-step gives the same
answer for two different generalized inverses \( \mathbf{K}\ginv \).

::: {.warning}
Software that inverts \( \V \) will fail here, and software that silently replaces
\( \V^{-1} \) by a pseudo-inverse solves the wrong problem: it minimizes the
\( \V^{+} \) criterion over all \( \bb \) instead of over \( \mathcal F(\y) \), which is
what @thm-ggm-singular(e) is for, so the point estimate can be wrong. In
@exm-ggm-effects the formula \( (\X\T\V^{+}\X)\ginv\X\T\V^{+}\y \) always returns
\( \hat\beta_0=0 \), contradicting \( \bone\T\hat{\boldsymbol{\upalpha}}=0 \). The smallest
instance has two observations, the second measured without error:
\( \V=\diag(1,0) \), \( \X=\bone_2 \), so \( \beta=Y_2 \) exactly, while the formula
reports \( Y_1 \). The degrees of freedom are wrong too: \( n-r \) is not \( q-d \). Here \( a-2=3 \)
happens to agree with \( n-r \), but with a \( \V \) of lower rank it would not; always
compute \( \rank(\X:\V)-\rank(\X) \).
:::

## Testing

A reduced model \( \C(\X_0)\subseteq\C(\X) \) is tested as in @thm-ggm-inference(e),
with two adjustments. The degrees of freedom are those of @thm-ggm-singular(f)
for each model, \( \nu=\rank(\X:\V)-\rank(\X) \) and
\( \nu_0=\rank(\X_0:\V)-\rank(\X_0) \), so the numerator carries \( \nu_0-\nu \). And
the reduced model can be *refuted with certainty*: if \( \Y\notin\C(\X_0:\V) \) it
is impossible and the test rejects at once. Otherwise @eq-ggm-f has the
\( F(\nu_0-\nu,\nu,\gamma) \) distribution; @exr-ggm-singular-f works through the
argument.

## Exercises

### A. Check your understanding

::: {#exr-ggm-zero-variance}
[A1]

Show that if \( \V \) is singular with \( \V\mathbf{a}=\bzero \), then \( \mathbf{a}\T\Y \) is the same
number for every possible data vector, and identify that number in terms of
\( \bbeta \). What does this say about \( \mathbf{a}\T\X\bbeta \) when \( \X\T\mathbf{a}=\bzero \)?
:::

::: {#exr-ggm-rank-formula}
[A2]

Verify \( \rank(\X:\V)=\rank(\X)+\rank(\V)-\dim(\C(\X)\cap\C(\V)) \) and compute both
sides for @exm-ggm-effects.
:::

::: {.solution}
This is the dimension formula for the sum of two subspaces,
\( \dim(\mathcal S+\mathcal T)=\dim\mathcal S+\dim\mathcal T-\dim(\mathcal S\cap\mathcal T) \).
In @exm-ggm-effects, \( \C(\X)+\C(\V)=\Real^a \) because \( \bone\in\C(\X) \) and
\( \C(\V)=\bone\perpc \); so the left side is \( a \), and the right side is
\( 2+(a-1)-1=a \).
:::

### B. Practice

::: {#exr-ggm-shares}
[B1]

Three fractions of a sample are measured and then rescaled so that they add to
one. Model \( \E(Y_k)=\mu+\delta_k \) with \( \sum_k\delta_k=0 \) known, and suppose
\( \Cov(\Y)=\sigma^2(\I_3-\bone\bone\T/3) \). Which function of \( (\mu,\delta_1,\delta_2) \)
is known exactly? Find the variance of \( Y_1-Y_2 \) as an estimate of
\( \delta_1-\delta_2 \), and explain why no other linear unbiased estimate exists.
:::

::: {.solution}
Here \( \mathbf{N}=\bone/\sqrt3 \) and \( \X\T\mathbf{N}\propto\X\T\bone \); with the parameterization
\( \bbeta=(\mu,\delta_1,\delta_2)\T \) and \( \delta_3=-\delta_1-\delta_2 \), we get
\( \X\T\bone=(3,0,0)\T \), so \( 3\mu \) is known exactly: \( \mu=\bone\T\Y/3=1/3 \), the
rescaling constant. Here \( \rank(\X)=3=n \), so \( \X\T \) is nonsingular and the
unbiasedness equation \( \X\T\mathbf{c}=\blambda \) has exactly one solution: every estimable
function has a single linear unbiased estimate, which is therefore its BLUE. For
\( \delta_1-\delta_2 \) that solution is \( \mathbf{c}=(1,-1,0)\T \), and
\( \Var(Y_1-Y_2)=\sigma^2\mathbf{c}\T\V\mathbf{c}=\sigma^2(2/3+2/3+2/3)=2\sigma^2 \).
:::

::: {#exr-ggm-deterministic-row}
[B2]

Suppose the first observation is measured without error, so
\( \V=\begin{pmatrix}0&\bzero\T\\\bzero&\V_{22}\end{pmatrix} \) with \( \V_{22} \)
positive definite, and \( \X=\begin{pmatrix}\x_{(1)}\T\\\X_2\end{pmatrix} \). Show
that \( \x_{(1)}\T\bbeta \) is known exactly, and that the BLUE of any estimable
\( \blambda\T\bbeta \) is the generalized least squares estimate from the remaining
\( n-1 \) observations *subject to* the constraint \( \x_{(1)}\T\bb=y_1 \).
:::

::: {.solution}
\( \mathbf{N}=\mathbf{e}_1 \) spans \( \C(\V)\perpc \), so \( \mathbf{N}\T\Y=Y_1=\x_{(1)}\T\bbeta \) with
probability one. In @eq-ggm-projection-form the feasible set is
\( \{\bb:\y-\X\bb\in\C(\V)\}=\{\bb:\x_{(1)}\T\bb=y_1\} \), and on \( \C(\V) \) the form
\( \mathbf{a}\T\V^{+}\mathbf{a} \) is \( \mathbf{a}_2\T\V_{22}^{-1}\mathbf{a}_2 \); so the criterion is the
generalized least squares criterion of the last \( n-1 \) observations, minimized
subject to the constraint. This is restricted least squares (@thm-ss-restricted)
in the \( \V_{22}^{-1} \) geometry.
:::

::: {#exr-ggm-two-step-check}
[B3]

Verify directly from @eq-ggm-two-step that \( \X\tilde{\bbeta} \) does not depend on
which generalized inverse \( \mathbf{K}\ginv \) is used, given \( \y\in\C(\X:\V) \).
:::

### C. Going deeper

::: {#exr-ggm-singular-f}
[C1]

In the setting of @thm-ggm-singular, let \( \C(\X_0)\subseteq\C(\X) \) with
\( \nu_0=\rank(\X_0:\V)-\rank(\X_0) \) and \( \nu=\rank(\X:\V)-\rank(\X) \), and assume
\( \Y \) is normal and \( \X\bbeta\in\C(\X_0:\V) \). Working in the coordinates
\( \Z_* \) of the proof of @thm-ggm-singular(f), show that
\( \text{SSE}_{\V,0}-\text{SSE}_{\V} \) and \( \text{SSE}_{\V} \) are independent with
\( (\text{SSE}_{\V,0}-\text{SSE}_{\V})/\sigma^2\sim\chi^2(\nu_0-\nu,\gamma) \), and
identify \( \gamma \). Why is the assumption \( \X\bbeta\in\C(\X_0:\V) \) needed?
:::

::: {.solution}
Because \( \X\bbeta\in\C(\X_0:\V) \) makes \( \mathbf{N}\T\Y\in\C(\mathbf{N}\T\X_0) \), the same
\( \bb_* \) is feasible in both models. In the coordinates
\( \Z_*\sim\Normal_q(\bmu_*,\sigma^2\I) \) the two sums of squares are
\( \norm{(\I-\M_d)\Z_*}^2 \) and \( \norm{(\I-\M_{d_0})\Z_*}^2 \) for nested subspaces of
dimensions \( d_0\le d \), with \( q-d=\nu \) and \( q-d_0=\nu_0 \). Now @thm-glh-f-test
applied to \( \Z_* \) gives independence, the \( \chi^2 \) laws and
\( \gamma=\norm{(\M_d-\M_{d_0})\bmu_*}^2/\sigma^2 \). If \( \X\bbeta\notin\C(\X_0:\V) \) the
reduced model assigns probability zero to the observed \( \Y \), so it is
contradicted outright and no distribution theory is needed.
:::

::: {#exr-ggm-lsce}
[C2]

Call \( \X\tilde{\bbeta} \) a *least squares consistent estimate* when it minimizes
the ordinary sum of squares \( \norm{\y-\X\bb}^2 \) over \( \mathcal F(\y) \); this is
the criterion used by Christensen (2020). Show that it agrees with @eq-ggm-projection-form when
\( \V \) is idempotent, and that in general the two differ. Show also that when
\( \C(\V\bU)\subseteq\C(\bU) \) the least squares consistent estimate is a BLUE.
:::
