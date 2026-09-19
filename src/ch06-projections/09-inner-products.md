# Other inner products

Everything so far has used the Euclidean inner product \( \bu\T\bv \). That choice treats
the \( n \) coordinates of observation space symmetrically: an error of a given size counts
the same in every observation, and errors in different observations are unrelated. When
observations have different variances or are correlated, Euclidean distance is the wrong
way to measure how far \( \y \) is from the model space. This section shows that the whole
theory carries over to any inner product, and that the result is generalized least
squares.

## Projection in a general inner product

Let \( \A \) be an \( n\times n \) symmetric positive definite matrix, and define
\[
\inner{\bu}{\bv}_{\A}=\bu\T\A\bv,\qquad \norm{\bu}_{\A}=(\bu\T\A\bu)^{1/2}.
\]
This is an inner product: bilinear, symmetric, and positive on nonzero vectors. Say
\( \bu \) and \( \bv \) are *\( \A \)-orthogonal* if \( \bu\T\A\bv=0 \). Nothing in the proofs of
@prp-proj-gram-schmidt, @thm-proj-direct-sum or
@thm-proj-projection-theorem used anything about \( \bu\T\bv \) except these three
properties. So those results hold verbatim with \( \inner{\cdot}{\cdot}_{\A} \) in place of
the Euclidean inner product. In particular, every \( \y \) has a unique
\( \A \)-nearest point in a subspace \( \mathcal S \), characterized by \( \A \)-orthogonality of
the error to \( \mathcal S \). The map sending \( \y \) to it is linear, and we write it
\( \bP_{\mathcal S,\A} \).

::: {#thm-proj-A-projection}
Let \( \A \) be positive definite and \( \X \) any \( n\times p \) matrix. The \( \A \)-orthogonal
projection onto \( \C(\X) \) is
\[
\bP_{\X,\A}=\X(\X\T\A\X)\ginv\X\T\A,
\]
for any generalized inverse. It satisfies \( \bP_{\X,\A}^2=\bP_{\X,\A} \) and
\( \C(\bP_{\X,\A})=\C(\X) \), and the matrix \( \A\bP_{\X,\A} \) is symmetric. The matrix
\( \bP_{\X,\A} \) itself is symmetric iff \( \A\C(\X)\subseteq\C(\X) \), and in general it is
not.
:::

::: {.proof}
Let \( \A=\bL\bL\T \) with \( \bL \) nonsingular (for instance the Cholesky factor), and put
\( \X_*=\bL\T\X \) and \( \y_*=\bL\T\y \). Then \( \norm{\y-\X\bb}_{\A}^2=\norm{\y_*-\X_*\bb}^2 \), so
\( \bb \) minimizes the \( \A \)-distance iff it is an ordinary least squares estimate for
\( (\y_*,\X_*) \). By @thm-proj-M-formula and @thm-proj-normal-equations, one such
estimate is \( (\X_*\T\X_*)\ginv\X_*\T\y_*=(\X\T\A\X)\ginv\X\T\A\y \), and
\( \X_*(\X_*\T\X_*)\ginv\X_*\T \) is invariant to the choice of generalized inverse. Multiplying
on the left by \( \bL^{-\top} \) and on the right by \( \bL\T \) shows that
\( \bP_{\X,\A}=\bL^{-\top}\M_*\bL\T \), where \( \M_* \) projects onto \( \C(\X_*) \). Idempotence and
\( \C(\bP_{\X,\A})=\bL^{-\top}\C(\X_*)=\C(\X) \) follow, and
\( \A\bP_{\X,\A}=\bL\M_*\bL\T \) is symmetric. For the last claim, see
@exr-proj-kruskal.
:::

The proof also gives the general principle: *\( \A \)-geometry is Euclidean
geometry after the change of variables \( \y\mapsto\bL\T\y \).* Any fact about orthogonal
projections has an \( \A \)-version, obtained by transporting it through \( \bL \).

Symmetry is the Euclidean signature of an orthogonal projection. The
\( \A \)-version of @thm-proj-sym-idem replaces it with self-adjointness.

::: {#prp-proj-A-selfadjoint}
An idempotent matrix \( \bP \) is the \( \A \)-orthogonal projection onto \( \C(\bP) \) iff \( \A\bP \)
is symmetric, that is, iff \( \inner{\bP\bu}{\bv}_{\A}=\inner{\bu}{\bP\bv}_{\A} \) for all \( \bu,\bv \).
:::

::: {.proof}
For idempotent \( \bP \), the residual \( (\I-\bP)\y \) is \( \A \)-orthogonal to \( \C(\bP) \) for
every \( \y \) iff \( \bP\T\A(\I-\bP)=\bzero \), that is, iff \( \bP\T\A=\bP\T\A\bP \). If \( \A\bP \)
is symmetric, then \( \bP\T\A=\A\bP \), and
\( \bP\T\A\bP=\A\bP\bP=\A\bP=\bP\T\A \). Conversely, \( \bP\T\A=\bP\T\A\bP \) has a symmetric
right-hand side, so \( \bP\T\A \) and hence its transpose \( \A\bP \) are symmetric.
:::

## Generalized least squares

Suppose \( \Cov(\be)=\sigma^2\V \) with \( \V \) known and positive definite. Observations with
large variance should count for less, and correlated observations partly repeat each
other. The natural distance is the one that makes the errors look uncorrelated:
\( \norm{\cdot}_{\V^{-1}} \). If \( \V=\bL\bL\T \), then \( \bL^{-1}\be \) has covariance \( \sigma^2\I \),
and \( \norm{\bu}_{\V^{-1}}=\norm{\bL^{-1}\bu} \).

::: {#def-proj-gls}
[Generalized least squares]

A **generalized least squares** (GLS) estimate minimizes
\( (\y-\X\bb)\T\V^{-1}(\y-\X\bb) \). Equivalently, \( \X\hbeta_{\text{GLS}}=\bP_{\X,\V^{-1}}\y \), and
\( \hbeta_{\text{GLS}} \) solves the generalized normal equations
\( \X\T\V^{-1}\X\bb=\X\T\V^{-1}\y \).
:::

@thm-proj-A-projection with \( \A=\V^{-1} \) is the whole of the algebra.
Computationally, GLS is ordinary least squares applied to the *whitened* data
\( (\bL^{-1}\y,\bL^{-1}\X) \). Aitken (1935) showed that GLS gives the best linear
unbiased estimator under this covariance, which is the Gauss–Markov theorem transported through
\( \bL \). Chapter 31 proves this, and handles
singular \( \V \) and the more realistic case where \( \V \) must be estimated.

When is the ordinary least squares fit already the GLS fit? Geometrically, the two
projections onto \( \C(\X) \) coincide iff the Euclidean and the \( \V^{-1} \) notions of
orthogonality agree on this particular subspace.

::: {#thm-proj-kruskal}
[Kruskal]

For positive definite \( \V \), \( \M\y=\bP_{\X,\V^{-1}}\y \) for every \( \y \) iff
\( \C(\V\X)\subseteq\C(\X) \).
:::

::: {.proof}
The two projections have the same range. They agree iff they have the same null
space, since an idempotent matrix is determined by its range and null space
(@prp-proj-oblique). The null space of \( \M \) is \( \C(\X)\perpc=\Null(\X\T) \). The
null space of \( \bP_{\X,\V^{-1}} \) is the set of \( \bu \) with \( \X\T\V^{-1}\bu=\bzero \), that is,
\( \V\,\Null(\X\T) \). So the projections agree iff \( \V\,\Null(\X\T)=\Null(\X\T) \), and since
\( \V \) is nonsingular and the spaces have equal dimension, iff \( \V\,\Null(\X\T)\subseteq\Null(\X\T) \).
By symmetry of \( \V \), this says \( \X\T\V\bu=\bzero \) whenever \( \X\T\bu=\bzero \), that is,
\( \C(\V\X)\perp\C(\X)\perpc \), which is \( \C(\V\X)\subseteq\C(\X) \).
:::

The condition holds, for example, for equicorrelated errors
\( \V=(1-\rho)\I+\rho\bone\bone\T \) whenever \( \bone\in\C(\X) \)
(@exr-proj-equicorrelation). It also underlies the analysis of balanced
split-plot designs in Chapter 33. The result is due to
Kruskal (1968), and the coordinate-free argument is his.

::: {#exm-proj-ar1}
[Autocorrelated errors]

Take eight equally spaced observations with a linear trend and error covariance
\( v_{st}=0.7^{\lvert s-t\rvert} \). The \( \V^{-1} \)-projection onto \( \C(\X) \) is idempotent and
self-adjoint for \( \inner{\cdot}{\cdot}_{\V^{-1}} \) but not symmetric: its largest
asymmetry \( \max\lvert p_{st}-p_{ts}\rvert \) is \( 0.444 \). The ordinary and
generalized slope estimates are \( 0.648 \) and
\( 0.650 \). The listing checks each claim, including that the GLS fit is
ordinary least squares on whitened data.
:::

```{.python .run #cell-other-inner-products-vproj}
import numpy as np
rng = np.random.default_rng(9)
n = 8
X = np.column_stack([np.ones(n), np.arange(n, dtype=float)])
y = X @ np.array([1.0, 0.5]) + rng.normal(size=n)
t = np.arange(n)
V = 0.7 ** np.abs(t[:, None] - t[None, :])             # AR(1)-type covariance

Vinv = np.linalg.inv(V)
P_V = X @ np.linalg.solve(X.T @ Vinv @ X, X.T @ Vinv)    # V^{-1}-orthogonal projection onto C(X)

print("idempotent:", np.allclose(P_V @ P_V, P_V))
print("symmetric: ", np.allclose(P_V, P_V.T))                   # False: oblique in the usual geometry
print("self-adjoint for <u,v> = u'V^{-1}v:", np.allclose(Vinv @ P_V, (Vinv @ P_V).T))

# GLS = ordinary least squares after whitening with a square root of V^{-1}
L = np.linalg.cholesky(V)                               # V = L L'
Xw, yw = np.linalg.solve(L, X), np.linalg.solve(L, y)
beta_whitened, *_ = np.linalg.lstsq(Xw, yw, rcond=None)
beta_gls = np.linalg.solve(X.T @ Vinv @ X, X.T @ Vinv @ y)
print("GLS via P_V equals OLS on whitened data:", np.allclose(beta_gls, beta_whitened))
```

## Weighted least squares

The simplest non-Euclidean geometry is a diagonal one. If the errors are uncorrelated
with variances \( \sigma^2/w_i \), then \( \V^{-1}=\bm W=\diag(w_1,\dots,w_n) \), and GLS
minimizes \( \sum_i w_i(y_i-\x_{(i)}\T\bb)^2 \). Observations believed to be noisier get
less weight in the distance. Whitening multiplies row \( i \) of \( \X \) and \( \y \) by
\( \sqrt{w_i} \). The fitted values are \( \bP_{\X,\bm W}\y \), and the natural leverages are
the diagonal entries of \( \bP_{\X,\bm W} \). These equal the diagonal of the orthogonal
projection onto \( \C(\bm W^{1/2}\X) \) and still sum to \( \rank(\X) \)
(@exr-proj-wls-leverage).

::: {#exm-proj-engel}
[Engel's food expenditure data]

Engel's 1857 survey of \( 235 \) Belgian working-class households records annual income and food
expenditure. The spread of expenditure grows with income, as [Figure 6.9.1](09-inner-products.html#fig-proj-wls)(a)
shows. Weighting by \( w_i=1/\text{income}_i^2 \) (standard deviation proportional to income)
changes the fitted slope from \( 0.4852 \) to \( 0.5740 \).
The change in geometry shows up more clearly in the leverages, panel (b). Under
ordinary least squares the richest households have the most leverage, with a maximum of
\( 0.255 \). In the weighted geometry their large variance shrinks their
influence, and no household has leverage above \( 0.037 \). Whether this
is the *right* geometry depends on whether the variance model holds.
Chapter 21 discusses how to check it, and how to
get valid standard errors when it fails.
:::

::: {when-format="html"}
![**Figure 6.9.1.** Engel's data. (a) Ordinary and weighted least squares fits, with weights
proportional to \( 1/\text{income}^2 \). (b) Leverages under each geometry. The weighted
projection spreads influence away from the high-variance, high-income households.](weighted_ls.svg){#fig-proj-wls width=100%}
:::

::: {when-format="pdf"}
![Engel's data. (a) Ordinary and weighted least squares fits, with weights
proportional to \( 1/\text{income}^2 \). (b) Leverages under each geometry. The weighted
projection spreads influence away from the high-variance, high-income households.](weighted_ls.pdf){width=100%}
:::

```{.python .run #cell-weighted-ls-wls}
import numpy as np
import statsmodels.api as sm
df = sm.datasets.engel.load_pandas().data
x, y = df["income"].to_numpy(), df["foodexp"].to_numpy()
X = np.column_stack([np.ones(len(x)), x])

w = 1.0 / x**2                                   # working model: sd proportional to income
sw = np.sqrt(w)
beta_ols, *_ = np.linalg.lstsq(X, y, rcond=None)
beta_wls, *_ = np.linalg.lstsq(X * sw[:, None], y * sw, rcond=None)   # OLS on (W^1/2 X, W^1/2 y)

# leverages in the weighted geometry: squared row norms of an orthonormal basis of C(W^1/2 X)
Qw, _ = np.linalg.qr(X * sw[:, None])
h_wls = np.sum(Qw**2, axis=1)
Q, _ = np.linalg.qr(X)
h_ols = np.sum(Q**2, axis=1)
print("OLS:", beta_ols, " WLS:", beta_wls)
print("max leverage OLS %.3f, WLS %.3f" % (h_ols.max(), h_wls.max()))
```

## Oblique projections

The matrices \( \bP_{\X,\A} \) are idempotent without being symmetric. In the Euclidean
geometry of \( \Real^n \) they are **oblique projections**. Idempotent matrices in general
have a clean description.

::: {#prp-proj-oblique}
Let \( \bP \) be idempotent. Then \( \Real^n=\C(\bP)\dirsum\Null(\bP) \), and \( \bP \) maps each
\( \y=\bu+\bv \) (\( \bu\in\C(\bP) \), \( \bv\in\Null(\bP) \)) to \( \bu \). Conversely, for any
decomposition \( \Real^n=\mathcal S\dirsum\mathcal T \) there is exactly one idempotent \( \bP \) with
\( \C(\bP)=\mathcal S \) and \( \Null(\bP)=\mathcal T \), the *projection onto \( \mathcal S \) along
\( \mathcal T \)*. It is an orthogonal projection iff \( \mathcal T=\mathcal S\perpc \), that is, iff
\( \bP \) is symmetric.
:::

::: {.proof}
Write \( \y=\bP\y+(\I-\bP)\y \). The first term lies in \( \C(\bP) \), and
\( \bP(\I-\bP)\y=\bzero \) puts the second in \( \Null(\bP) \). If \( \bu=\bP\bw\in\Null(\bP) \)
then \( \bu=\bP^2\bw=\bP\bu=\bzero \), so the sum is direct. For \( \bu\in\C(\bP) \),
\( \bP\bu=\bu \) (write \( \bu=\bP\bw \)), so \( \bP(\bu+\bv)=\bu \). Conversely, given
\( \mathcal S\dirsum\mathcal T \), the map \( \bu+\bv\mapsto\bu \) is well defined and linear by
uniqueness of the decomposition, and it is idempotent with the required range and null
space. It is unique because its values are forced on all of \( \Real^n \). If
\( \mathcal T=\mathcal S\perpc \), the map is the orthogonal projection, which is symmetric. If
\( \bP \) is symmetric, then \( \Null(\bP)=\Null(\bP\T)=\C(\bP)\perpc \) by
@lem-proj-null-colspace.
:::

::: {when-format="html"}
![**Figure 6.9.2.** Projections onto a line \( \mathcal S \) in \( \Real^2 \). (a) The orthogonal projection:
the residual is perpendicular to \( \mathcal S \), and \( \bP\y \) is the nearest point. (b) An
oblique projection moves \( \y \) along a fixed direction that is not perpendicular to
\( \mathcal S \). The result is idempotent but is not the nearest point, and the
residual is longer.](oblique_projection.svg){#fig-proj-oblique width=85%}
:::

::: {when-format="pdf"}
![Projections onto a line \( \mathcal S \) in \( \Real^2 \). (a) The orthogonal projection:
the residual is perpendicular to \( \mathcal S \), and \( \bP\y \) is the nearest point. (b) An
oblique projection moves \( \y \) along a fixed direction that is not perpendicular to
\( \mathcal S \). The result is idempotent but is not the nearest point, and the
residual is longer.](oblique_projection.pdf){width=85%}
:::

[Figure 6.9.2](09-inner-products.html#fig-proj-oblique) shows the difference. The nearest-point property belongs
to the orthogonal projection alone. An oblique projection can move a vector
*farther* away: its operator norm exceeds one unless it is orthogonal
(@exr-proj-oblique-norm). The GLS projection is still a nearest-point map,
but in the \( \V^{-1} \) norm, not the Euclidean one. Which geometry is “right”
is not a mathematical question. It is decided by the covariance of the errors.

## Exercises

### A. Check your understanding

::: {#exr-proj-equicorrelation}
[A1]

Let \( \V=(1-\rho)\I+\rho\bone\bone\T \) with \( -1/(n-1)<\rho<1 \). Show that \( \V \) is positive
definite and that \( \C(\V\X)\subseteq\C(\X) \) whenever \( \bone\in\C(\X) \). Conclude that
ordinary least squares is generalized least squares for equicorrelated errors in any
model with an intercept.
:::

::: {.solution}
\( \V\bone=(1+(n-1)\rho)\bone \) and
\( \V\bu=(1-\rho)\bu \) for \( \bu\perp\bone \). So the eigenvalues are \( 1+(n-1)\rho \) and
\( 1-\rho \), both positive exactly when \( -1/(n-1)<\rho<1 \). Then
\( \V\X=(1-\rho)\X+\rho\bone(\bone\T\X) \), whose columns lie in \( \C(\X) \) when \( \bone \) does.
@thm-proj-kruskal applies.
:::

### B. Practice

::: {#exr-proj-kruskal}
[B1]

Complete the proof of @thm-proj-A-projection. Show that \( \bP_{\X,\A} \) is symmetric
iff it equals \( \M \), and deduce from @thm-proj-kruskal (with \( \V=\A^{-1} \)) that this
happens iff \( \A\C(\X)\subseteq\C(\X) \).
:::

::: {#exr-proj-oblique-norm}
[B2]

Let \( \bP\neq\bzero \) be idempotent. Show that the operator norm
\( \norm{\bP}_2=\max_{\norm{\bu}=1}\norm{\bP\bu} \) is at least \( 1 \), with equality iff \( \bP \)
is symmetric. *Hint:* if \( \bP \) is not symmetric, find \( \bu\in\C(\bP) \) and
\( \bv\in\Null(\bP) \) with \( \bu\T\bv\ne0 \), and consider \( \bu+t\bv \) for small \( t \).
:::

::: {.solution}
For \( \bu\in\C(\bP) \), \( \bP\bu=\bu \), so
\( \norm{\bP}_2\ge1 \). If \( \bP \) is symmetric, \( \norm{\bP\bu}\le\norm{\bu} \) by
@prp-proj-trace-rank(d). Suppose \( \bP \) is not symmetric. Then
\( \Null(\bP)\neq\C(\bP)\perpc \) by @prp-proj-oblique. The two spaces have the same
dimension, so some \( \bv\in\Null(\bP) \) is not orthogonal to \( \C(\bP) \). Choose a unit
\( \bu\in\C(\bP) \) with \( \bu\T\bv>0 \). For \( t>0 \), \( \bP(\bu-t\bv)=\bu \), while
\( \norm{\bu-t\bv}^2=1-2t\,\bu\T\bv+t^2\norm{\bv}^2<1 \) for small \( t \). So \( \bP \) stretches
\( \bu-t\bv \), and \( \norm{\bP}_2>1 \).
:::

::: {#exr-proj-wls-leverage}
[B3]

For weighted least squares, \( \V=\diag(1/w_1,\dots,1/w_n) \) with \( w_i>0 \). Show that the
diagonal entries of \( \bP_{\X,\V^{-1}} \) equal those of the orthogonal projection onto
\( \C(\bm W^{1/2}\X) \), where \( \bm W=\diag(w_i) \), and that they sum to \( \rank(\X) \). Interpret
them as leverages.
:::
