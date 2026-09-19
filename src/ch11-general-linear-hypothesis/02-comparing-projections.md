# The F statistic as a comparison of projections

The formula @eq-glh-F mentions two residual sums of squares, but the test is really about three
subspaces and the way the data vector sits among them. Seen this way, the \( F \) statistic is a
function of a single angle, its rejection region is a cone, and its null distribution survives
far beyond normal errors. The same view gives a canonical form in which every linear hypothesis
looks the same. [Section 11.4](04-likelihood-ratio.html) needs that form to prove that the
\( F \) test is optimal.

## Three orthogonal pieces

By @thm-proj-nested, nested model spaces split \( \Real^n \) into three mutually orthogonal
subspaces:
\[
\Real^n=\underbrace{\C(\X_0)}_{\dim r_0}\;\dirsum\;
\underbrace{\C(\X)\cap\C(\X_0)\perpc}_{\dim q=r-r_0}\;\dirsum\;
\underbrace{\C(\X)\perpc}_{\dim n-r},
\qquad
\I=\Mo+(\M-\Mo)+(\I-\M).
\]{#eq-glh-three-spaces}

The middle subspace, \( \C(\M-\Mo) \), is the **test space** of the hypothesis. It consists of the
directions in which the full model can move the mean and the reduced model cannot. Every data vector
splits accordingly,
\[
\y=\Mo\y+(\M-\Mo)\y+(\I-\M)\y ,
\]
and so does the mean vector, \( \X\bbeta=\Mo\X\bbeta+(\M-\Mo)\X\bbeta+\bzero \). The hypothesis
says that the middle component of the *mean* is zero. The middle component of the *data*,
\( (\M-\Mo)\Y \), is its unbiased estimate, with covariance \( \sigma^2(\M-\Mo) \). The residual
component \( (\I-\M)\Y \) has mean zero whether or not the hypothesis holds, and it serves as a
yardstick for noise. The first component \( \Mo\Y \) estimates the part of the mean that both models
allow. It carries no information about the hypothesis, and the \( F \) statistic ignores it.

What the statistic does use is \( (\I-\Mo)\y \), the residual vector of the reduced model. It is
split into a piece in the test space and a piece in the residual space, and
\[
F=\frac{\norm{(\M-\Mo)\y}^2/q}{\norm{(\I-\M)\y}^2/(n-r)}
\]
compares their squared lengths per dimension. Under the hypothesis both pieces are pure noise, and
each dimension carries on average \( \sigma^2 \). Under the alternative the first piece also carries
the signal \( (\M-\Mo)\X\bbeta \).

The test space can be described by a spanning set. Since \( \M-\Mo=(\I-\Mo)\M \) by
@thm-proj-nested(a), and \( \C(\M)=\C(\X) \),
\[
\C(\M-\Mo)=\C\bigl((\I-\Mo)\X\bigr).
\]
The test space is spanned by the columns of \( \X \) after the reduced model has been regressed out
of them. When \( \X=[\X_0,\X_2] \) adds a block \( \X_2 \), this is \( \C((\I-\Mo)\X_2) \), the space of @lem-ss-residualized. For a hypothesis on coefficients it is the constraint space
\( \C(\M\bT) \) of @thm-proj-constraint-space, as [Section 11.3](03-testable-hypotheses.html)
shows.

## One angle

Because the two pieces of \( (\I-\Mo)\y \) are orthogonal, their lengths are the legs of a right
triangle whose hypotenuse is the reduced model's residual vector.

::: {#prp-glh-angle}
[\( F \) as a function of an angle]

Let \( \bw=(\I-\Mo)\y\ne\bzero \) and let \( \theta\in[0,\pi/2] \) be the angle between \( \bw \) and the
test space, defined by \( \cos\theta=\norm{(\M-\Mo)\y}/\norm{\bw} \). Then
\[
\norm{(\M-\Mo)\y}^2=\norm{\bw}^2\cos^2\theta,\qquad
\text{SSE}=\norm{\bw}^2\sin^2\theta,\qquad
F=\frac{n-r}{q}\cot^2\theta .
\]
The partial coefficient of determination of @def-ss-partial-r2 is \( \cos^2\theta \), and the
\( F \) test rejects exactly when \( \theta<\theta_\alpha \), where
\( \tan^2\theta_\alpha=(n-r)/\{q\,F_\alpha(q,n-r)\} \).
:::

::: {.proof}
\( \bw=(\M-\Mo)\y+(\I-\M)\y \) with orthogonal summands, since
\( (\M-\Mo)(\I-\M)=\bzero \). Pythagoras gives
\( \norm{(\I-\M)\y}^2=\norm{\bw}^2-\norm{\bw}^2\cos^2\theta=\norm{\bw}^2\sin^2\theta \). Since
\( (\M-\Mo)\y \) is the projection of \( \bw \) onto the test space, \( \theta \) is also the angle
between \( \bw \) and its projection, the smallest angle between \( \bw \) and any vector of the
test space. The ratio of the two squared lengths is \( \cot^2\theta \), which gives \( F \). The partial
coefficient of determination is \( (\text{SSE}_0-\text{SSE})/\text{SSE}_0=\norm{(\M-\Mo)\y}^2/\norm{\bw}^2=\cos^2\theta \).
Finally \( \cot^2 \) is strictly decreasing on \( (0,\pi/2] \), so \( F>F_\alpha \) iff
\( \cot^2\theta>qF_\alpha/(n-r) \) iff \( \theta<\theta_\alpha \).
:::

So the \( F \) test asks a geometric question: *does the reduced model's residual vector point
suspiciously close to the test space?* Under the hypothesis, \( \bw \) is noise spread over
\( n-r_0 \) dimensions, only \( q \) of which belong to the test space, so it typically points well away
from it. The rejection region \( \{\y:\theta<\theta_\alpha\} \) is a circular cone around the test
space inside \( \C(\X_0)\perpc \), extended by arbitrary components in \( \C(\X_0) \). That this is the
same as the partial \( R^2 \) test of @thm-ss-partial-r2(b) is no coincidence. The partial \( R^2 \) is a
squared cosine, just as @thm-proj-r2-cosine showed for the ordinary \( R^2 \).

::: {#exm-glh-region-angle}
[The region test as an angle]

In @exm-glh-region, the reduced model's residual vector has a component of length
\( 3.851 \) in the three-dimensional test space and a component of length \( 9.334 \) in
the \( 43 \)-dimensional residual space. The angle between the residual vector and the test space is
\( \theta=67.6^\circ \), and the partial \( R^2 \) is
\( \cos^2\theta=0.145 \). The \( 5\% \) test rejects when \( \theta<66.1^\circ \), so the
data fall just outside the cone.

[Figure 11.2.1](#fig-glh-cone) plots the two lengths for the data and for \( 4000 \)
simulated data sets. With the reduced model true, \( 0.053 \) of them fall in the rejection
region. With the Northeast lowered by two murders per 100,000, and all else as fitted, the fraction is
\( 0.774 \).
:::

::: {when-format="html"}
![**Figure 11.2.1.** The region test of @exm-glh-region in two coordinates: the length of
\( (\M-\Mo)\y \) in the test space and the length of the residual vector \( (\I-\M)\y \). The \( F \)
statistic depends only on the angle of the point, and the test rejects inside the shaded wedge.
The grey points are simulated with the reduced model true, the green points with the Northeast
lowered by two murders per 100,000. The observed data are the black point, just above the
wedge.](test_cone.svg){#fig-glh-cone width=70%}
:::

::: {when-format="pdf"}
![The region test of @exm-glh-region in two coordinates: the length of
\( (\M-\Mo)\y \) in the test space and the length of the residual vector \( (\I-\M)\y \). The \( F \)
statistic depends only on the angle of the point, and the test rejects inside the shaded wedge.
The grey points are simulated with the reduced model true, the green points with the Northeast
lowered by two murders per 100,000. The observed data are the black point, just above the
wedge.](test_cone.pdf){width=70%}
:::

```{.python .run #cell-test-geometry-angle}
import numpy as np
import statsmodels.api as sm
from scipy import stats

data = sm.datasets.statecrime.load_pandas().data
data.index = data.index.str.strip()
data = data.drop(index="District of Columbia")
codes = "SWWSWWNSSSWWMMMMSSNSNMMSMWMWNNWNSMMSWNNSMSSWNSWSMW"
region = np.array(list(codes))
y = data["murder"].to_numpy()
n = len(y)
covariates = np.column_stack([data["poverty"], data["single"], data["urban"]])
X0 = np.column_stack([np.ones(n), covariates])
D = np.column_stack([(region == g).astype(float) for g in "NSW"])
X = np.column_stack([X0, D])

def proj(Z):
    """Orthogonal projection onto C(Z), from an orthonormal basis."""
    Q, _ = np.linalg.qr(Z)
    return Q @ Q.T

M, M0 = proj(X), proj(X0)
q, nu = 3, n - 7                               # r - r0 and n - r
u = (M - M0) @ y                               # component in the test space
e = y - M @ y                                  # residual vector
w = y - M0 @ y                                 # what the reduced model leaves unexplained
theta = np.arccos(np.linalg.norm(u) / np.linalg.norm(w))
F = (nu / q) / np.tan(theta) ** 2
theta_crit = np.arctan(np.sqrt(nu / (q * stats.f.ppf(0.95, q, nu))))
print(f"angle = {np.degrees(theta):.1f} degrees, F = {F:.3f}")
print(f"reject when the angle is below {np.degrees(theta_crit):.1f} degrees")
```

## What the statistic ignores

Since \( F \) depends on \( \y \) only through the angle \( \theta \), many transformations of the data
leave it unchanged. They form the symmetry group of the testing problem, and
[Section 11.4](04-likelihood-ratio.html) uses them to characterize the \( F \) test.

::: {#prp-glh-invariance}
[Invariance of the \( F \) statistic]

The statistic @eq-glh-F depends on \( \X \) and \( \X_0 \) only through \( \C(\X) \) and \( \C(\X_0) \).
Where it is defined, that is, when \( (\I-\M)\y\ne\bzero \), it is unchanged when \( \y \) is replaced by
\[
a\,\mathbf{H}\y+\bv,\qquad a\ne0,\quad \bv\in\C(\X_0),
\]
where \( \mathbf{H} \) is any orthogonal matrix that maps \( \C(\X_0) \) onto itself and \( \C(\X) \) onto itself.
:::

::: {.proof}
The first claim holds because \( \M \) and \( \Mo \) are determined by the column spaces. For the second,
the transformations can be applied one at a time. Adding \( \bv \) changes nothing, because \( \M-\Mo \)
and \( \I-\M \) annihilate \( \C(\X_0) \). The factor \( a^2 \) cancels in the ratio. For \( \mathbf{H} \): the
matrix \( \mathbf{H}\M\mathbf{H}\T \) is symmetric and idempotent with column space
\( \mathbf{H}\C(\X)=\C(\X) \), so it equals \( \M \) (@thm-proj-sym-idem). Hence \( \mathbf{H}\M=\M\mathbf{H} \), and
likewise \( \mathbf{H}\Mo=\Mo\mathbf{H} \). Then
\( \norm{(\M-\Mo)\mathbf{H}\y}=\norm{\mathbf{H}(\M-\Mo)\y}=\norm{(\M-\Mo)\y} \), because \( \mathbf{H} \) preserves
lengths, and the same holds for \( \I-\M \).
:::

The orthogonal matrices allowed here are exactly those that rotate the test space within itself
and the residual space within itself, and fix \( \C(\X_0) \) as a set. The invariance says that no
direction within the test space is favoured. A departure of given size \( \norm{(\M-\Mo)\X\bbeta} \)
is equally detectable in every direction. When some directions matter more than others, for
instance when one contrast is of primary interest, a test directed at them can have more power
there. It pays for this with less power elsewhere ([Section 11.6](06-coefficients.html)).

## Canonical form

Orthonormal bases adapted to @eq-glh-three-spaces turn every linear hypothesis into the same
problem about independent normal coordinates.

::: {#prp-glh-canonical}
[Canonical form of a linear hypothesis]

Let \( \Q=[\Q_0,\Q_1,\Q_2] \) be an orthogonal \( n\times n \) matrix whose blocks of \( r_0 \), \( q \) and
\( n-r \) columns are orthonormal bases of \( \C(\X_0) \), of the test space and of \( \C(\X)\perpc \). Put
\( \bz_k=\Q_k\T\Y \) and \( \boldsymbol{\upeta}_k=\Q_k\T\X\bbeta \), \( k=0,1,2 \). Under @eq-opt-normal-model:

::: {.enumerate options="label=(\alph*)"}
1. \( \bz_0,\bz_1,\bz_2 \) are independent, \( \bz_k\sim\Normal(\boldsymbol{\upeta}_k,\sigma^2\I) \), and \( \boldsymbol{\upeta}_2=\bzero \);

2. \( \norm{\bz_1}^2=\text{SSE}_0-\text{SSE} \) and \( \norm{\bz_2}^2=\text{SSE} \), so
   \( F=\{\norm{\bz_1}^2/q\}/\{\norm{\bz_2}^2/(n-r)\} \);

3. \( \norm{\boldsymbol{\upeta}_1}^2=\sigma^2\gamma \), and \( H_0 \) holds iff \( \boldsymbol{\upeta}_1=\bzero \).
:::

:::

::: {.proof}
(a) \( \Q\T\Y\sim\Normal_n(\Q\T\X\bbeta,\sigma^2\Q\T\Q)=\Normal_n(\Q\T\X\bbeta,\sigma^2\I) \) by @thm-mvn-linear. Its blocks are uncorrelated, hence independent (@thm-mvn-independence). The columns of
\( \Q_2 \) are orthogonal to \( \C(\X) \), so \( \boldsymbol{\upeta}_2=\Q_2\T\X\bbeta=\bzero \). (b) By
@prp-proj-orthonormal-formula, \( \Q_1\Q_1\T=\M-\Mo \) and \( \Q_2\Q_2\T=\I-\M \), so
\( \norm{\bz_1}^2=\Y\T\Q_1\Q_1\T\Y=\Y\T(\M-\Mo)\Y \) and similarly for \( \bz_2 \). (c) In the same way
\( \norm{\boldsymbol{\upeta}_1}^2=\norm{(\M-\Mo)\X\bbeta}^2 \), and @thm-glh-f-test(c) applies.
:::

In canonical form the problem is stripped to its essentials. One observes a \( q \)-vector \( \bz_1 \),
normal with unknown mean \( \boldsymbol{\upeta}_1 \) and covariance \( \sigma^2\I \), and wishes to test
\( \boldsymbol{\upeta}_1=\bzero \). An independent vector \( \bz_2 \) of pure noise, of dimension \( n-r \), measures
\( \sigma^2 \). The vector \( \bz_0 \) has an unknown mean that neither hypothesis restricts, so it is of no
use. Every hypothesis of this chapter, whether about coefficients, contrasts or whole models, is an
instance of this problem with particular \( q \) and \( n-r \).

## Beyond normal errors

The angle picture shows why normality is less important for the null distribution than it first
appears. Under \( H_0 \), \( (\I-\Mo)\Y=(\I-\Mo)\be \), so \( F \) is a function of the error vector alone,
and it is a function of its *direction* only: replacing \( \be \) by \( c\be \), \( c>0 \), leaves every
angle unchanged. Any error law under which the direction of \( \be \) behaves as it does for normal
errors therefore gives the same null distribution. We call \( \mathbf{U} \) **uniformly distributed on the
unit sphere** if it has the distribution of \( \Z/\norm{\Z} \) for \( \Z\sim\Normal_n(\bzero,\I) \).

::: {#prp-glh-spherical}
[Spherical errors]

Let \( \Y=\X\bbeta+\be \) with \( \X\bbeta\in\C(\X_0) \), where \( \be=R\,\mathbf{U} \) for a random variable
\( R>0 \) and a vector \( \mathbf{U} \) uniformly distributed on the unit sphere. The variables \( R \) and
\( \mathbf{U} \) need not be independent. Then \( F\sim F(r-r_0,n-r) \).
:::

::: {.proof}
Write \( g(\bv)=\{\norm{(\M-\Mo)\bv}^2/q\}/\{\norm{(\I-\M)\bv}^2/(n-r)\} \) for vectors with
\( (\I-\M)\bv\neq\bzero \). Under \( H_0 \), \( (\M-\Mo)\Y=(\M-\Mo)\be \) and \( (\I-\M)\Y=(\I-\M)\be \),
since both projections annihilate \( \X\bbeta\in\C(\X_0) \). So \( F=g(\be) \). Because
\( g(c\bv)=g(\bv) \) for \( c>0 \), \( F=g(R\mathbf{U})=g(\mathbf{U}) \). The law of \( g(\mathbf{U}) \) is the law of
\( g(\Z/\norm{\Z})=g(\Z) \) with \( \Z\sim\Normal_n(\bzero,\I) \). By @thm-glh-f-test(d) with
\( \sigma=1 \), that law is \( F(r-r_0,n-r) \). (The events \( \Z=\bzero \) and \( (\I-\M)\Z=\bzero \) have
probability zero.)
:::

The class of error laws covered is large. If \( \Z\sim\Normal_n(\bzero,\I) \) and \( S>0 \) is any
random scale independent of \( \Z \), or not, then \( \be=S\Z=(S\norm{\Z})(\Z/\norm{\Z}) \) qualifies.
Taking \( S^{-2} \) proportional to a \( \chi^2(\nu) \) variable gives the spherical multivariate
\( t \) distribution. Its components are uncorrelated and heavy-tailed but *not* independent: one
common scale inflates or shrinks them all. The \( F \) test is exactly valid for all such laws. What
the proposition does not cover are independent non-normal errors, whose direction is not uniform.
For them the null distribution is only approximately \( F \).

::: {#exm-glh-spherical}
[Four error laws]

Keep the design of @exm-glh-region and simulate \( 20000 \) data sets with the reduced
model true, each with errors from one of four laws: independent normal; spherical \( t \) with three
degrees of freedom; independent \( t \) with three degrees of freedom; independent centred
exponential. The spherical errors are uncorrelated, but their sizes move together: the rank
correlation between the absolute values of two components is \( 0.20 \). The
rejection rates of the \( 5\% \) test are

| Errors | normal | spherical \( t_3 \) | independent \( t_3 \) | exponential |
|:---|:---:|:---:|:---:|:---:|
| Rejection rate | \( 0.0502 \) | \( 0.0502 \) | \( 0.0444 \) | \( 0.0474 \) |

with Monte Carlo standard error about \( 0.0015 \). The first two agree exactly, not just
approximately: the spherical errors were made by rescaling the normal ones, and each rescaled data
set gives the same \( F \) as its normal original. The two independent non-normal laws give
rates slightly below \( 0.05 \). For this design, with \( 43 \) residual degrees of freedom and no
extreme leverage, the independent \( t_3 \) errors make the test slightly conservative, the exponential
rate is within Monte Carlo error of \( 0.05 \), and both are close to the nominal level.
:::

```{.python .run #cell-null-distribution-spherical}
import numpy as np
import statsmodels.api as sm
from scipy import stats

rng = np.random.default_rng(1103)

data = sm.datasets.statecrime.load_pandas().data
data.index = data.index.str.strip()
data = data.drop(index="District of Columbia")
codes = "SWWSWWNSSSWWMMMMSSNSNMMSMWMWNNWNSMMSWNNSMSSWNSWSMW"
region = np.array(list(codes))
n = len(data)
covariates = np.column_stack([data["poverty"], data["single"], data["urban"]])
X0 = np.column_stack([np.ones(n), covariates])
D = np.column_stack([(region == g).astype(float) for g in "NSW"])
X = np.column_stack([X0, D])

def proj(Z):
    Q, _ = np.linalg.qr(Z)
    return Q @ Q.T

M, M0 = proj(X), proj(X0)
q, nu = 3, n - 7
crit = stats.f.ppf(0.95, q, nu)

def F_stat(E):
    """F statistic for each row of E, used as a data vector under the reduced model."""
    num = np.sum((E @ (M - M0)) ** 2, axis=1) / q
    den = np.sum((E - E @ M) ** 2, axis=1) / nu
    return num / den

reps = 20000
Z = rng.standard_normal((reps, n))
scale = 1 / np.sqrt(rng.chisquare(3, size=(reps, 1)) / 3)       # one scale per data set
errors = {
    "normal": Z,
    "spherical t3": Z * scale,                                    # dependent, heavy-tailed
    "independent t3": rng.standard_t(3, size=(reps, n)),
    "exponential": rng.exponential(size=(reps, n)) - 1.0,         # skewed
}
for name, E in errors.items():
    print(f"{name:15s} rejection rate {np.mean(F_stat(E) > crit):.4f}")
```

::: {.remark}
[Only the null distribution]

@prp-glh-spherical concerns the level of the test. Under an alternative, \( F \) depends on the
error length \( R \) as well as the direction, and the power under spherical \( t \) errors differs from
the power under normal errors with the same covariance. Nor does the proposition say that the
\( F \) test is a *good* test for heavy-tailed errors. Its power can be poor compared with robust
methods. Those questions belong to Chapter 19 and Chapter 21.
:::

## Exercises

### A. Check your understanding

::: {#exr-glh-angle-numbers}
[A1]

Use @prp-glh-angle to express \( \cos^2\theta \) in terms of \( F \), \( q \) and \( n-r \). Check the value
\( 0.145 \) of @exm-glh-region-angle from \( F=2.440 \), \( q=3 \) and
\( n-r=43 \), and compute \( \theta_\alpha \) from \( F_{0.05}(3,43)=2.822 \).
:::

### B. Practice

::: {#exr-glh-beta-null}
[B1]

Show that under \( H_0 \) the squared cosine \( \cos^2\theta \) of @prp-glh-angle has the
\( \mathrm{Beta}\bigl(q/2,(n-r)/2\bigr) \) distribution, and that it is independent of
\( \norm{(\I-\Mo)\Y}^2 \). What is its expectation?
:::

::: {.solution}
By @prp-glh-canonical, \( \cos^2\theta=U/(U+V) \) with \( U=\norm{\bz_1}^2/\sigma^2\sim\chi^2(q) \) and
\( V=\norm{\bz_2}^2/\sigma^2\sim\chi^2(n-r) \) independent under \( H_0 \), and
\( \norm{(\I-\Mo)\Y}^2=\sigma^2(U+V) \). @exr-qf-beta gives both the Beta law and the independence
of \( U/(U+V) \) from \( U+V \). The mean of \( \mathrm{Beta}(a,b) \) is \( a/(a+b) \), here
\( q/(n-r_0) \): under the hypothesis the residual vector of the reduced model puts, on average, a
share of its squared length in the test space equal to that space's share of the dimensions.
:::

::: {#exr-glh-test-space}
[B2]

Prove that \( \C(\M-\Mo)=\C\bigl((\I-\Mo)\X\bigr) \), and that for \( \X=[\X_0,\X_2] \) this equals
\( \C\bigl((\I-\Mo)\X_2\bigr) \). Show that \( \dim\C((\I-\Mo)\X_2)=\rank(\X)-\rank(\X_0) \).
:::

::: {.solution}
\( \M-\Mo=(\I-\Mo)\M \) by @thm-proj-nested(a), so
\( \C(\M-\Mo)\subseteq\C((\I-\Mo)\M)=(\I-\Mo)\C(\M)=(\I-\Mo)\C(\X)=\C((\I-\Mo)\X) \). Conversely each
column \( (\I-\Mo)\x_j=\x_j-\Mo\x_j \) lies in \( \C(\X) \) and is orthogonal to \( \C(\X_0) \), so it lies in
\( \C(\X)\cap\C(\X_0)\perpc=\C(\M-\Mo) \). For \( \X=[\X_0,\X_2] \), \( (\I-\Mo)\X=[\bzero,(\I-\Mo)\X_2] \).
The dimension is \( \rank(\M-\Mo)=r-r_0 \) (@thm-proj-nested(c)).
:::

::: {#exr-glh-spherical-examples}
[B3]

Which of the following error vectors satisfy the hypothesis of @prp-glh-spherical?
(i) \( \be=\sigma\Z \) with \( \sigma \) fixed; (ii) \( \be=S\Z \) with \( S \) uniform on \( (1,2) \) and independent of
\( \Z \); (iii) \( \be \) with independent Laplace components; (iv) \( \be=\mathbf{D}\Z \) with \( \mathbf{D} \) a fixed
diagonal matrix with unequal entries; (v) \( \be=\Z/\norm{\Z}^2 \).
:::

### C. Going deeper

::: {#exr-glh-maximal-invariant}
[C1]

In the canonical form, let two data vectors \( \y \) and \( \y' \) with \( \bz_2,\bz_2'\ne\bzero \) have the same
value of \( F \). Construct \( a>0 \), \( \bv\in\C(\X_0) \) and an orthogonal \( \mathbf{H} \) of the kind allowed in
@prp-glh-invariance such that \( \y'=a\mathbf{H}\y+\bv \). (So \( F \) is a *maximal invariant*: any statistic
unchanged by these transformations is a function of \( F \).)
:::

::: {#exr-glh-non-spherical-size}
[C2]

Take a one-way layout with two groups of sizes \( n_1=5 \) and \( n_2=20 \), independent normal errors with
standard deviations \( 3 \) in the small group and \( 1 \) in the large one, and equal means. Compute, by
simulation or numerical integration, the true level of the nominal \( 5\% \) \( F \) test of equal means.
Explain the direction of the error using the angle picture: which way does the unequal variance
tilt the residual vector?
:::
