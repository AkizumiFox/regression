# Penalized splines

The last section ended with a dilemma. A regression spline with few knots is
too rigid to follow the data; one with many knots overfits; and the number of
knots is a discrete parameter that cannot be tuned gently. The resolution is to
stop treating the knots as the control: take *too many* of them — enough that the
basis could interpolate — and shrink the coefficients towards a smooth
configuration with a quadratic penalty. This is @def-reg-penalized with a
particular penalty, and by @thm-shr-ridge it is ridge regression in the B-spline
basis.

## The difference penalty

If the basis is rich, neighbouring B-splines overlap heavily and neighbouring
coefficients describe nearly the same part of the curve. A wiggly fit is one
whose coefficient sequence \( \gamma_1,\dots,\gamma_m \) jumps about, a smooth
fit one whose coefficients change slowly with \( j \), so penalizing the
differences of the coefficients penalizes roughness — with no integrals to
evaluate, since the penalty is a quadratic form whose matrix is known in
advance.

Write \( \Delta\gamma_j=\gamma_j-\gamma_{j-1} \) and
\( \Delta^k=\Delta(\Delta^{k-1}) \), so that
\( \Delta^2\gamma_j=\gamma_j-2\gamma_{j-1}+\gamma_{j-2} \). Let \( \bD_k \) be the
\( (m-k)\times m \) matrix with
\( (\bD_k\bgamma)_j=\Delta^k\gamma_{j+k} \); for \( k=1 \) and \( k=2 \) its rows are
\( (\dots,-1,1,\dots) \) and \( (\dots,1,-2,1,\dots) \).

::: {#def-smo-pspline}
[Penalized spline]

Let \( \B \) be the B-spline basis matrix of @def-smo-bspline for \( K \) interior
knots, equally spaced, with \( K \) chosen generously (say \( 20 \) to \( 40 \), and
in any case no more than \( n \)). For \( \lambda\ge0 \) and an integer
\( k\ge1 \), the **penalized spline** or **P-spline** estimate is the function
\( \hat f=\sum_j\hat\gamma_jB_j^d \) whose coefficients minimize
\[
\norm{\y-\B\bgamma}^2+\lambda\,\norm{\bD_k\bgamma}^2
=\norm{\y-\B\bgamma}^2+\lambda\,\bgamma\T\bP\bgamma,
\qquad \bP=\bD_k\T\bD_k .
\]{#eq-smo-pspline}

\( \lambda \) is the **smoothing parameter** and \( k \) the **order** of the
penalty.
:::

This is the construction of Eilers and Marx (1996), who gave it the name
P-spline. In the scaling of @def-reg-penalized the criterion carries a factor
\( \tfrac12 \) and \( P(\bgamma)=\tfrac12\bgamma\T\bP\bgamma \); nothing depends on
the convention except the numerical value of \( \lambda \).

Three things have changed since [Section 43.2](02-b-splines.html). The number of
knots is a computational decision, not a modelling one: take enough. The estimate
is no longer a projection, so it has bias even when \( f \) lies in the spline
space. And the dial is continuous.

## The fit, and what the penalty does to it

::: {#lem-smo-diagonalize}
[The two matrices diagonalize together]

Suppose \( \B \) has full column rank \( m \), and let \( \bP \) be nonnegative
definite. There is a nonsingular \( \bT \) with
\( \bT\T\B\T\B\,\bT=\I \) and \( \bT\T\bP\bT=\diag(e_1,\dots,e_m) \), where
\( 0\le e_1\le\dots\le e_m \) and exactly \( \rank(\bP) \) of the \( e_j \) are
positive. The columns \( \bu_1,\dots,\bu_m \) of \( \B\bT \) are orthonormal.
:::

::: {.proof}
\( \B\T\B \) is positive definite, so @cor-mat-generalized-rayleigh(b) supplies a
nonsingular \( \bT \) with \( \bT\T\B\T\B\bT=\I \) and \( \bT\T\bP\bT \) diagonal.
Its entries are nonnegative because \( \bP \) is nonnegative definite, and the
number of nonzero ones is \( \rank(\bT\T\bP\bT)=\rank(\bP) \); permuting the
columns of \( \bT \) puts them in increasing order and preserves both
identities. Orthonormality of
the columns of \( \B\bT \) is the statement \( (\B\bT)\T(\B\bT)=\I \).
:::

The vectors \( \bu_j \) are the **Demmler–Reinsch basis** of the smoother. Ordered
by \( e_j \), they run from smooth to rough: those with \( e_j=0 \) are the shapes
the penalty does not see at all.

::: {#prp-smo-pspline-ridge}
[A P-spline is a ridge regression]

Let \( \B \) have full column rank and \( \lambda>0 \), and write
\( \bT,e_j,\bu_j \) as in @lem-smo-diagonalize.

::: {.enumerate options="label=(\alph*)"}
1. \( \hat{\bgamma}=(\B\T\B+\lambda\bP)^{-1}\B\T\y \) is the unique minimizer of
   @eq-smo-pspline, and the fitted vector is \( \hat{\mathbf f}=\bS_\lambda\y \) with
   \[
   \bS_\lambda=\B(\B\T\B+\lambda\bP)^{-1}\B\T
   =\sum_{j=1}^m\frac{1}{1+\lambda e_j}\,\bu_j\bu_j\T .
   \]{#eq-smo-smoother}

2. \( \bS_\lambda \) is symmetric with eigenvalues \( (1+\lambda e_j)^{-1}\in(0,1] \)
   and the same eigenvectors for every \( \lambda \). Its **effective degrees of
   freedom** — the count of @thm-shr-ridge(b), now for a roughness penalty —
   \[
   \mathrm{df}(\lambda)=\tr(\bS_\lambda)=\sum_{j=1}^m\frac{1}{1+\lambda e_j}
   \]{#eq-smo-df}

   decrease strictly and continuously from \( m \) at \( \lambda=0 \) to
   \( m-\rank(\bP) \) as \( \lambda\to\infty \). For a \( k \)th-order difference
   penalty, \( \rank(\bP)=m-k \), so the limit is \( k \).

3. As \( \lambda\to\infty \), \( \bS_\lambda \) tends to the orthogonal projection
   onto \( \C\{\B\,\Null(\bP)\} \), a space of dimension \( k \).
:::

:::

::: {.proof}
(a) The criterion is \( \bgamma\T(\B\T\B+\lambda\bP)\bgamma-2\bgamma\T\B\T\y+\y\T\y \),
and \( \B\T\B+\lambda\bP \) is positive definite because \( \B\T\B \) is and
\( \bP \) is nonnegative definite; @prp-mat-quadratic-min gives the unique
minimizer. In the coordinates \( \bgamma=\bT\boldsymbol{\upalpha} \) the criterion
becomes
\( \norm{\y-\B\bT\boldsymbol{\upalpha}}^2+\lambda\sum_je_j\alpha_j^2 \), and because
the columns of \( \B\bT \) are orthonormal this separates:
\( \norm{\y}^2-\sum_j\{2\alpha_j\bu_j\T\y-\alpha_j^2(1+\lambda e_j)\} \) plus the
part of \( \y \) orthogonal to all the \( \bu_j \). Minimizing term by term gives
\( \hat\alpha_j=\bu_j\T\y/(1+\lambda e_j) \) and
\( \hat{\mathbf f}=\sum_j\hat\alpha_j\bu_j \), which is @eq-smo-smoother.

(b) The representation @eq-smo-smoother is a spectral decomposition, since the
\( \bu_j \) are orthonormal; the trace is the sum of the eigenvalues. Each term
\( (1+\lambda e_j)^{-1} \) is constant in \( \lambda \) if \( e_j=0 \) and strictly
decreasing to zero otherwise, and there are \( \rank(\bP) \) of the latter. For
the difference penalty, \( \bD_k \) has \( m-k \) rows and full row rank (its
leading \( (m-k)\times(m-k) \) submatrix is triangular with \( \pm1 \) on the
diagonal), so \( \rank(\bP)=\rank(\bD_k)=m-k \).

(c) By @eq-smo-smoother, \( \bS_\lambda\to\sum_{j:e_j=0}\bu_j\bu_j\T \), the
orthogonal projection onto the span of those \( \bu_j \). Those columns of
\( \B\bT \) are \( \B \) times the columns of \( \bT \) with \( \bT\T\bP\bT \)-entry
zero, and \( \bT\T\bP\bT\vect{e}_j=0 \) with \( \bT \) nonsingular means
\( \bT\vect{e}_j\in\Null(\bP) \); counting dimensions gives all of
\( \Null(\bP) \).
:::

::: {.idea}
Compare @thm-shr-ridge(b). Ridge shrinks the coordinate along the \( j \)th
singular direction by \( d_j^2/(d_j^2+\lambda) \), and the small singular values
are shrunk most. A P-spline shrinks the coordinate along \( \bu_j \) by
\( 1/(1+\lambda e_j) \), and the *rough* directions are shrunk most. The
difference is only in what the geometry calls "small": ridge measures size by the
design, a smoother measures it by roughness. Both count what is left with a
trace.
:::

## What survives as the penalty grows

Part (c) says the limit is a projection onto a \( k \)-dimensional space. Which
one?

::: {#prp-smo-penalty-limit}
[The limiting fit]

Let \( \bD_k \) be the \( k \)th difference matrix.

::: {.enumerate options="label=(\alph*)"}
1. \( \Null(\bD_k)=\{\bgamma:\gamma_j=p(j)\text{ for a polynomial }p\text{ of degree}<k\} \),
   of dimension \( k \).

2. For \( k=1 \), the limiting fit of @prp-smo-pspline-ridge(c) is the constant
   \( \bar y \), whatever the knots.

3. For \( k=2 \) and equally spaced Greville abscissae (@prp-smo-greville) — the
   uniform knot sequence of [Section 43.2](02-b-splines.html) — the limiting fit
   is the ordinary least squares line.
:::

:::

::: {.proof}
(a) If \( \gamma_j=p(j) \) with \( \deg p<k \) then \( \Delta \) lowers the degree by
one, so \( \Delta^k\gamma=0 \). Conversely the set of such \( \bgamma \) has
dimension \( k \) (the polynomials \( 1,j,\dots,j^{k-1} \) restricted to
\( j=1,\dots,m \) are independent as long as \( m\ge k \)), and
\( \dim\Null(\bD_k)=m-\rank(\bD_k)=k \) by the proof
of @prp-smo-pspline-ridge(b). Two subspaces of the same dimension, one inside the
other, coincide.

(b) For \( k=1 \), \( \Null(\bD_1) \) is spanned by \( \bone \), and
\( \B\bone \) is the vector of values of \( \sum_jB_j^d \), which is \( \bone \)
by @prp-smo-bspline(c). The projection onto \( \spn\{\bone\} \) sends \( \y \) to
\( \bar y\bone \).

(c) For \( k=2 \), \( \Null(\bD_2) \) is spanned by \( \bone \) and
\( \mathbf{j}=(1,\dots,m)\T \). If the Greville abscissae satisfy
\( \xi_j=\xi_1+(j-1)\delta \), then \( \mathbf{j}=\{\bone(\delta-\xi_1)+\boldsymbol{\upxi}\}/\delta \)
with \( \boldsymbol{\upxi}=(\xi_1,\dots,\xi_m)\T \), so
\( \spn\{\bone,\mathbf{j}\}=\spn\{\bone,\boldsymbol{\upxi}\} \). By @prp-smo-bspline(c)
and @eq-smo-greville, \( \B\bone \) is the vector of ones and
\( \B\boldsymbol{\upxi} \) is the vector \( (x_1,\dots,x_n)\T \). So
\( \C\{\B\Null(\bP)\} \) is the column space of \( [\bone,\x] \), and the projection
onto it is the least squares line.
:::

::: {.remark}
[Why the knot sequence matters here and nowhere else]

Part (c) is exactly true for the uniform sequence and only approximately true for
the clamped sequence @eq-smo-clamped, whose Greville abscissae bunch up near the
endpoints. On the Nile data of @exm-smo-knot-insensitivity, against a response
range of \( 914 \), the uniform fit at
\( \lambda=10^{12} \) departs from the least squares line by
\( 3.2\times 10^{-2} \) — the numerical remainder of a
penalty that is large but finite — while the clamped fit departs by
\( 11.4 \), some
\( 1.2 \) per cent of that range: visually
almost nothing, but not the line. Since the two sequences span the same
functions, the choice affects only what "the penalty is zero" means. For
\( k\ge3 \) the limiting space is a \( k \)-dimensional space of splines
containing the straight lines, and consists exactly of the polynomials of degree
\( k-1 \) when the knots are equally spaced (Eilers and Marx, 1996; de Boor,
2001, ch. IX), a fact we use but do not prove.
:::

## Knots stop mattering

::: {#exm-smo-knot-insensitivity}
[Nile flow, 1871 to 1970]

The annual flow of the Nile at Aswan was recorded for
\( 100 \) years from \( 1871 \) to
\( 1970 \). Rescale the year to \( [0,1] \) and fit cubic
P-splines with a second-order penalty. Reaching
\( 3 \), \( 8 \) and \( 20 \) effective degrees of freedom needs
\( \lambda=9.1\times 10^{2} \),
\( 4.7\times 10^{0} \) and
\( 1.5\times 10^{-2} \): five orders of magnitude for
a threefold change in df, which is why \( \lambda \) is searched on a logarithmic
scale and reported as df.

Now hold df at \( 8 \) and change the number of knots.
With \( K=10,20,40 \) and \( 80 \) — the last giving
\( 84 \) coefficients for \( 100 \) observations — the four
curves differ nowhere by more than
\( 16.60 \), or
\( 1.82 \) per cent of the response range, and the
\( K=40 \) and \( K=80 \) curves by at most
\( 1.373 \). Without the penalty the same knot counts
give curves differing by up to \( 258.2 \), fifteen
times as much.

[Figure 43.3.1](03-penalized-splines.html#fig-smo-pspline) draws both
comparisons. The fitted curve descends steeply around 1900, and it is worth
saying what a smoother does and does not know: the historical record is of an
abrupt change, and a smoother with a second-order penalty cannot represent a
jump. It renders one as a steep smooth descent and reports no lack of
fit (@thm-cor-lack-of-fit would, if the design had replicates). Smoothness is an
assumption, and like every assumption in this book it buys precision by excluding
something.
:::

::: {when-format="html"}
![**Figure 43.3.1.** Annual Nile flow. (a) Cubic P-splines with a second-order
penalty, labelled by effective degrees of freedom @eq-smo-df. (b) The same fit at
df \( =8 \) with \( K=10,20,40 \) and \( 80 \) interior knots: four curves are
drawn and one is visible.](pspline_nile.svg){#fig-smo-pspline width=100%}
:::

::: {when-format="pdf"}
![Annual Nile flow. (a) Cubic P-splines with a second-order
penalty, labelled by effective degrees of freedom @eq-smo-df. (b) The same fit at
df \( =8 \) with \( K=10,20,40 \) and \( 80 \) interior knots: four curves are
drawn and one is visible.](pspline_nile.pdf){width=100%}
:::

```{.python .run #cell-pspline-pspline}
import numpy as np
import statsmodels.api as sm


def uniform_knots(a, b, K, d):
    """K interior knots equally spaced in (a, b), continued d further steps at each end."""
    delta = (b - a) / (K + 1)
    return a + delta * np.arange(-d, K + 2 + d)


def bspline_basis(x, kn, d):
    """The B-splines of degree d on the knot vector kn, by the de Boor recurrence."""
    x = np.atleast_1d(np.asarray(x, float))
    B = np.array([(x >= kn[j]) & (x < kn[j + 1]) for j in range(len(kn) - 1)], float).T
    last = np.max(np.nonzero(kn[:-1] < kn[1:])[0])       # close the rightmost interval
    B[x == kn[-1], last] = 1.0
    for deg in range(1, d + 1):
        C = np.zeros((len(x), B.shape[1] - 1))
        for j in range(C.shape[1]):
            if kn[j + deg] > kn[j]:
                C[:, j] += (x - kn[j]) / (kn[j + deg] - kn[j]) * B[:, j]
            if kn[j + deg + 1] > kn[j + 1]:
                C[:, j] += (kn[j + deg + 1] - x) / (kn[j + deg + 1] - kn[j + 1]) * B[:, j + 1]
        B = C
    return B


def difference_matrix(m, k):
    """The k-th difference operator on m coefficients, an (m-k) x m matrix."""
    return np.diff(np.eye(m), n=k, axis=0)


def pspline(x, y, lam, K=20, d=3, k=2, a=None, b=None):
    """Penalized spline fit: minimize ||y - B gamma||^2 + lam ||Delta^k gamma||^2."""
    a = x.min() if a is None else a
    b = x.max() if b is None else b
    kn = uniform_knots(a, b, K, d)
    B = bspline_basis(x, kn, d)
    Dk = difference_matrix(B.shape[1], k)
    A = B.T @ B + lam * Dk.T @ Dk
    gamma = np.linalg.solve(A, B.T @ y)
    df = np.trace(np.linalg.solve(A, B.T @ B))          # trace of the smoother matrix
    return gamma, df, kn, d


def evaluate(gamma, kn, d, grid):
    return bspline_basis(grid, kn, d) @ gamma


data = sm.datasets.nile.load_pandas().data
year = data["year"].to_numpy()
y = data["volume"].to_numpy()
n = len(y)
x = (year - year.min()) / (year.max() - year.min())      # the covariate, rescaled to [0, 1]
grid = np.linspace(0, 1, 401)

gamma_big, df_big, kn, d = pspline(x, y, 1e12)
line = np.polyval(np.polyfit(x, y, 1), grid)
print("df at lambda = 1e12:", round(df_big, 6))
print("largest gap from the least squares line:", np.max(np.abs(evaluate(gamma_big, kn, d, grid) - line)))
```

## Choices that remain

**How many knots.** Enough that adding more changes nothing, which
@exm-smo-knot-insensitivity suggests is a few dozen; Ruppert's rule of thumb is
\( K=\min(n/4,35) \), and by @exr-smo-bandwidth-matrix generosity is cheap.

**The degree.** Cubic (\( d=3 \)) throughout; higher degrees buy nothing once a
penalty is present.

**The order of the penalty.** \( k=2 \) is the default, because its null space is
the straight lines: when the data say there is no curvature, the fit degrades to
the linear model of Part II rather than to a horizontal line (\( k=1 \)) or a
parabola (\( k=3 \)). The order also fixes the behaviour *beyond* the data, where
the penalty is all there is: a \( k \)th-order penalty extrapolates as a
polynomial of degree \( k-1 \), one more reason to distrust an extrapolated
smooth (@prp-ply-instability).

**The smoothing parameter.** [Section 43.5](05-choosing-lambda.html).

::: {.remark}
[Derivative penalties]

A different penalty measures roughness by an integral,
\( \int_a^b\{f''(t)\}^2dt=\bgamma\T\boldsymbol{\Omega}\bgamma \) with
\( \Omega_{jk}=\int B_j''B_k'' \), computable exactly from @eq-smo-bderiv and
Gauss–Legendre quadrature and banded for the same reason \( \B\T\B \) is. This is
the O'Sullivan (1986) penalty, and it makes the P-spline a finite-dimensional
version of the smoothing spline of the next section; Wand and Ormerod (2008)
show that for equally spaced knots @def-smo-pspline's penalty is a difference
approximation to it, which is why the two fits look identical. The difference
penalty is simpler to state and easier to generalize — nothing in @eq-smo-pspline
mentions an integral — and is the one used in
[Chapter 44](../ch44-additive-models/index.html).
:::

## Exercises

### A. Check your understanding

::: {#exr-smo-df-monotone}
[A1]

Using @eq-smo-df, sketch \( \mathrm{df}(\lambda) \) against \( \log\lambda \) for
\( m=24 \) and a second-order penalty. What are the two horizontal asymptotes,
and why is \( \log\lambda \) the right scale for the horizontal axis?
:::

### B. Practice

::: {#exr-smo-df-not-integer}
[B1]

A model with \( 8 \) parameters uses \( 8 \) degrees of freedom; a P-spline at
\( \mathrm{df}(\lambda)=8 \) uses \( 24 \) parameters. Explain,
using @eq-smo-smoother, in what sense the two are comparable, and show that
\( \tr(\bS_\lambda\bS_\lambda\T)\le\tr(\bS_\lambda) \) with equality only when
\( \bS_\lambda \) is a projection. (Both traces are used as degrees of freedom in
the literature; see [Section 43.5](05-choosing-lambda.html).)
:::

::: {.solution}
By @eq-smo-smoother, \( \bS_\lambda \) has eigenvalues
\( s_j=(1+\lambda e_j)^{-1}\in(0,1] \), so \( \tr(\bS_\lambda)=\sum_js_j \) and
\( \tr(\bS_\lambda\bS_\lambda\T)=\sum_js_j^2\le\sum_js_j \), with equality iff
every \( s_j\in\{0,1\} \), that is iff \( \bS_\lambda \) is a
projection (@thm-proj-sym-idem). A least squares fit on \( p \) columns has
\( p \) eigenvalues equal to one; a smoother spreads the same total over more
directions, using each only partly. This is the covariance-penalty degrees of
freedom @def-sel-df applied to a linear rule.
:::

::: {#exr-smo-ridge-equivalence}
[B2]

Show that the P-spline criterion @eq-smo-pspline is @thm-sel-loocv's criterion
\( \norm{\y-\X\bb}^2+\bb\T\boldsymbol{\Omega}\bb \) with \( \X=\B \) and
\( \boldsymbol{\Omega}=\lambda\bP \), and conclude that the leave-one-out shortcut
\( (y_i-\hat f(x_i))/(1-s_{ii}) \) is exact for a P-spline. What changes if the
knots are chosen by looking at the data first?
:::

### C. Going deeper

::: {#exr-smo-eigenvectors}
[C1]

Compute the Demmler–Reinsch basis of @lem-smo-diagonalize for a cubic P-spline
with \( K=20 \), \( k=2 \), on \( 100 \) equally spaced design points, and plot
\( \bu_1,\bu_2,\bu_5,\bu_{10} \) against \( x \). Confirm that \( e_1=e_2=0 \) and
that the corresponding \( \bu_j \) span the constant and linear functions, and that
\( \bu_j \) crosses zero more often as \( e_j \) grows. Relate the picture to
@eq-smo-smoother: which shapes survive a large \( \lambda \)?
:::

::: {#exr-smo-rank-deficient}
[C2]

@prp-smo-pspline-ridge assumes \( \B \) has full column rank, which fails if some
knot interval contains no design point. Show that the minimizer of
@eq-smo-pspline is still unique provided
\( \Null(\B)\cap\Null(\bP)=\{\bzero\} \), and give an example with \( K \) large
and a gap in the design where this fails for \( k=1 \) but holds for \( k=2 \)…
or show that it cannot fail for \( k=2 \) either. (Hint: @thm-reg-existence(c),
and the fact that \( \Null(\bP)=\Null(\bD_k) \).)
:::

::: {.solution}
\( \B\T\B+\lambda\bP \) is nonnegative definite, and
\( \bgamma\T(\B\T\B+\lambda\bP)\bgamma=0 \) forces \( \B\bgamma=\bzero \) and
\( \bD_k\bgamma=\bzero \); so the matrix is positive definite, and the minimizer
unique, exactly when the two null spaces meet only at zero. This is
@thm-reg-existence(c) with a strictly convex penalty on the relevant subspace.
For \( k=2 \), \( \bgamma\in\Null(\bD_2) \) means \( \gamma_j=\alpha+\beta j \), and
by @prp-smo-penalty-limit(c) the corresponding function is
\( \alpha'+\beta'x \) with \( \beta'\propto\beta \); such a function vanishes at all
\( n\ge2 \) distinct design points only if \( \alpha'=\beta'=0 \). So the
intersection is trivial for \( k=2 \) whatever the gaps, and the fit is always
unique. For \( k=1 \) the same argument with constants gives the same conclusion.
The moral is that a penalty repairs rank deficiency exactly on its own null space
— which is why a generous \( K \) is safe.
:::

::: {#exr-smo-shrinkage-bias}
[C3]

Assume \( \E(\Y)=\mathbf f \) and \( \Cov(\Y)=\sigma^2\I \).
Using @eq-smo-smoother, show that
\[
\E\norm{\hat{\mathbf f}-\mathbf f}^2
=\sum_{j}\Bigl(\frac{\lambda e_j}{1+\lambda e_j}\Bigr)^2(\bu_j\T\mathbf f)^2
+\sigma^2\sum_j\frac{1}{(1+\lambda e_j)^2}
+\norm{(\I-\bU\bU\T)\mathbf f}^2,
\]
where \( \bU=[\bu_1,\dots,\bu_m] \). Identify the three terms, and say which of
them no choice of \( \lambda \) can reduce.
:::

::: {.solution}
Write \( \mathbf f=\bU\bU\T\mathbf f+(\I-\bU\bU\T)\mathbf f \). Since
\( \bS_\lambda=\bU\diag(s_j)\bU\T \) with \( s_j=(1+\lambda e_j)^{-1} \),
\( \hat{\mathbf f}-\mathbf f=\bU\{\diag(s_j)\bU\T(\mathbf f+\be)-\bU\T\mathbf f\}
-(\I-\bU\bU\T)\mathbf f \), and the three pieces are orthogonal in expectation.
The first sum is the squared shrinkage bias,
\( (s_j-1)\bu_j\T\mathbf f=-\lambda e_j\bu_j\T\mathbf f/(1+\lambda e_j) \); the
second the variance \( \sigma^2\sum_js_j^2 \); the third the *approximation*
error, the part of \( \mathbf f \) outside the spline space, which does not
involve \( \lambda \). Only a richer basis reduces it, which is the argument for
a generous \( K \).
:::
