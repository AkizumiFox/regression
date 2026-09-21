# Repeated measures

A repeated-measures study observes the same subject on several occasions: a
clustered design whose cluster is the subject and whose "treatments" are the
occasions, laid out like the split-plot design of
[Section 33.1](01-split-plot.html), and for many years analysed that way.

The analogy has one weak joint. What makes the split-plot covariance compound
symmetric is the *randomization* of subplot treatments within the whole plot, and
occasions cannot be randomized: measurements close together in time are more alike
than distant ones, so compound symmetry here is usually false. What is remarkable
is that the univariate \( F \) test survives a much weaker condition, and that the
condition can be written down exactly.

## The univariate analysis

Let \( \y_1,\dots,\y_N \) be independent \( \Normal_m(\bmu,\bSigma) \) vectors, one per
subject, with \( \bSigma \) positive definite and unknown, and write \( y_{ij} \) for
occasion \( j \) of subject \( i \). The hypothesis of interest is
\[
H:\ \mu_1=\dots=\mu_m,\qquad\text{that is,}\qquad \bmu\in\spn\{\bone\} .
\]

Laid out with subjects as rows and occasions as columns the data have one
observation per cell, and the classical analysis is the additive two-way analysis
of variance of [Section 16.1](../ch16-multiway-layouts/01-additive.html), with the
subject-by-occasion interaction as error. With \( \bar{y}_{i\cdot} \),
\( \bar{y}_{\cdot j} \) and \( \bar{y} \) the subject, occasion and grand means,
\[
\text{SS}_T=N\sum_{j=1}^m(\bar{y}_{\cdot j}-\bar{y})^2,\qquad
\text{SS}_E=\sum_{i=1}^N\sum_{j=1}^m(y_{ij}-\bar{y}_{i\cdot}-\bar{y}_{\cdot j}+\bar{y})^2,
\]{#eq-cls-rm-ss}

and the statistic is
\( F=\{\text{SS}_T/(m-1)\}/\{\text{SS}_E/((N-1)(m-1))\} \). It is the subplot
stratum of a split-plot design with one whole plot per subject, and under
compound symmetry it is exact (@thm-cls-sphericity(b) below,
with @exr-cls-cs-implies-sphericity). The question is how much of that exactness
survives a general \( \bSigma \).

Everything reduces to the within-subject contrasts. Let \( \bU \) be
\( m\times(m-1) \) with orthonormal columns orthogonal to \( \bone \), so that
\( \bU\bU\T=\I-\bar{\mathbf{J}}_m \), and set \( \bz_i=\bU\T\y_i \): independent
\( \Normal_{m-1}(\bU\T\bmu,\boldsymbol{\Psi}) \) with
\( \boldsymbol{\Psi}=\bU\T\bSigma\bU \), with \( H \) saying exactly that
\( \bU\T\bmu=\bzero \). The subject means carry no information about \( H \).

## Sphericity

::: {#def-cls-sphericity}
[Sphericity and Box's measure]

Let \( \bSigma \) be an \( m\times m \) positive definite matrix and \( \bU \) as above.
\( \bSigma \) satisfies **sphericity** if \( \bU\T\bSigma\bU=\lambda\I_{m-1} \) for some
\( \lambda>0 \). **Box's measure of sphericity** is
\[
\epsilon(\bSigma)=\frac{\{\tr\boldsymbol{\Psi}\}^2}{(m-1)\tr(\boldsymbol{\Psi}^2)},
\qquad \boldsymbol{\Psi}=\bU\T\bSigma\bU .
\]{#eq-cls-epsilon}

Neither the condition nor \( \epsilon \) depends on the choice of \( \bU \): two such
matrices satisfy \( \tilde{\bU}=\bU\Q \) with \( \Q \) orthogonal, so
\( \tilde{\boldsymbol{\Psi}}=\Q\T\boldsymbol{\Psi}\Q \) has the same trace, the same
trace of square, and is a multiple of the identity exactly when
\( \boldsymbol{\Psi} \) is.
:::

::: {#lem-cls-hf-form}
[Equivalent forms of sphericity]

For a positive definite \( m\times m \) matrix \( \bSigma \) the following are
equivalent:

::: {.enumerate options="label=(\roman*)"}
1. \( \bU\T\bSigma\bU=\lambda\I_{m-1} \);

2. \( \Var(Y_j-Y_k)=2\lambda \) for every \( j\ne k \);

3. \( \sigma_{jk}=a_j+a_k+\lambda\delta_{jk} \) for some \( a_1,\dots,a_m \), where
   \( \delta_{jk} \) is one if \( j=k \) and zero otherwise;

4. \( \bSigma=\lambda\I+\bone\mathbf{a}\T+\mathbf{a}\bone\T \) for some
   \( \mathbf{a}\in\Real^m \).
:::

Moreover \( 1/(m-1)\le\epsilon(\bSigma)\le1 \), with \( \epsilon=1 \) iff sphericity
holds, and compound symmetry implies sphericity but not conversely.
:::

::: {.proof}
Write \( \mathbf{C}=\I-\bar{\mathbf{J}}_m=\bU\bU\T \). Since \( \bU\T\bU=\I_{m-1} \),
condition (i) holds iff \( \mathbf{C}\bSigma\mathbf{C}=\lambda\mathbf{C} \).

(i) \( \Rightarrow \) (ii): \( \vect{e}_j-\vect{e}_k \) is orthogonal to \( \bone \),
hence fixed by \( \mathbf{C} \), so \( \Var(Y_j-Y_k) \) equals
\( (\vect{e}_j-\vect{e}_k)\T\mathbf{C}\bSigma\mathbf{C}(\vect{e}_j-\vect{e}_k) \), which
is \( \lambda\norm{\vect{e}_j-\vect{e}_k}^2=2\lambda \).

(ii) \( \Rightarrow \) (iii): \( \Var(Y_j-Y_k)=\sigma_{jj}+\sigma_{kk}-2\sigma_{jk}=2\lambda \)
gives \( \sigma_{jk}=(\sigma_{jj}-\lambda)/2+(\sigma_{kk}-\lambda)/2 \) for \( j\ne k \).
Put \( a_j=(\sigma_{jj}-\lambda)/2 \); then \( \sigma_{jk}=a_j+a_k \) off the diagonal and
\( \sigma_{jj}=2a_j+\lambda \) on it, which is (iii).

(iii) \( \Leftrightarrow \) (iv) is the same statement in matrix form.

(iv) \( \Rightarrow \) (i): \( \mathbf{C}\bone=\bzero \), so
\( \mathbf{C}\bSigma\mathbf{C}=\lambda\mathbf{C} \).

For the bounds, let \( \psi_1,\dots,\psi_{m-1}>0 \) be the eigenvalues of
\( \boldsymbol{\Psi} \). Then
\( \epsilon=(\sum_l\psi_l)^2/\{(m-1)\sum_l\psi_l^2\} \). The Cauchy–Schwarz inequality
gives \( (\sum_l\psi_l)^2\le(m-1)\sum_l\psi_l^2 \) with equality iff all \( \psi_l \) are
equal, so \( \epsilon\le1 \) with equality iff \( \boldsymbol{\Psi}=\lambda\I \); and
\( \sum_l\psi_l^2\le(\sum_l\psi_l)^2 \) for nonnegative numbers gives
\( \epsilon\ge1/(m-1) \). Compound symmetry is (iv) with all \( a_j \) equal; taking the
\( a_j \) unequal gives sphericity without compound symmetry, and such a \( \bSigma \)
is positive definite as soon as \( \lambda \) is large enough.
:::

Condition (ii) is the form to remember and to check: *all pairwise differences of
the repeated measurements have the same variance*. Compound symmetry says more
than the \( F \) test needs.

::: {#thm-cls-sphericity}
[Sphericity is necessary and sufficient for the exact F test]

Let \( \y_1,\dots,\y_N \) be independent \( \Normal_m(\bmu,\bSigma) \) with \( \bSigma \)
positive definite, \( m\ge2 \), \( N\ge2 \), and let \( F \) be the statistic built
from @eq-cls-rm-ss. Put \( \boldsymbol{\Psi}=\bU\T\bSigma\bU \),
\( \mathbf{q}=\sqrt N\,\bar{\bz} \) and \( \W=\sum_i(\bz_i-\bar{\bz})(\bz_i-\bar{\bz})\T \).

::: {.enumerate options="label=(\alph*)"}
1. \( \text{SS}_T=\norm{\mathbf{q}}^2 \) and \( \text{SS}_E=\tr(\W) \), and under \( H \) the
   vector \( \mathbf{q}\sim\Normal_{m-1}(\bzero,\boldsymbol{\Psi}) \) is independent of
   \( \W\sim W_{m-1}(N-1,\boldsymbol{\Psi}) \). Hence
   \( F=N(N-1)\norm{\bar{\bz}}^2/\tr(\W) \).

2. If \( \bSigma \) satisfies sphericity, then under \( H \)
   \( F\sim F(m-1,(N-1)(m-1)) \) exactly, for every \( N\ge2 \).

3. Conversely, if the null distribution of \( F \) is \( F(m-1,(N-1)(m-1)) \) for every
   \( N\ge2 \), then \( \bSigma \) satisfies sphericity.
:::

:::

::: {.proof}
(a) Since \( \mathbf{C}=\bU\bU\T \) and \( \bU \) has orthonormal columns,
\( \sum_j(\bar{y}_{\cdot j}-\bar{y})^2=\norm{\mathbf{C}\bar{\y}}^2=\norm{\bar{\bz}}^2 \),
so \( \text{SS}_T=N\norm{\bar{\bz}}^2=\norm{\mathbf{q}}^2 \); likewise subject \( i \)'s
contribution to \( \text{SS}_E \) is \( \norm{\bz_i-\bar{\bz}}^2 \), summing to
\( \tr(\W) \). Under \( H \) the \( \bz_i \) are independent
\( \Normal_{m-1}(\bzero,\boldsymbol{\Psi}) \), so @lem-cor-centred-wishart gives the
Wishart law of \( \W \) and its independence of \( \bar{\bz} \). Finally
\( F=\{\norm{\mathbf{q}}^2/(m-1)\}\{(N-1)(m-1)/\tr\W\} \).

(b) If \( \boldsymbol{\Psi}=\lambda\I_{m-1} \), then \( \mathbf{q}/\sqrt\lambda \) is
standard normal, so \( \norm{\mathbf{q}}^2/\lambda\sim\chi^2(m-1) \). And
\( \W/\lambda\sim W_{m-1}(N-1,\I) \), which is the matrix of sums of squares and
cross products of \( (N-1)(m-1) \) independent standard normal variables, so
\( \tr(\W)/\lambda\sim\chi^2((N-1)(m-1)) \). The two are independent, and \( \lambda \)
cancels from the ratio, leaving the ratio of two independent chi-squared variables
each divided by its degrees of freedom (@def-qf-noncentral-f).

(c) Write \( F_N \) for the statistic with \( N \) subjects. By (a),
\[
F_N=\norm{\mathbf{q}}^2\cdot\frac{N-1}{\tr(\W)},
\]
where \( \mathbf{q}\sim\Normal_{m-1}(\bzero,\boldsymbol{\Psi}) \) has the *same*
distribution for every \( N \) and is independent of \( \W \). Since
\( \W/(N-1) \) is an average of \( N-1 \) independent \( W_{m-1}(1,\boldsymbol{\Psi}) \)
matrices, the weak law of large numbers gives \( \tr(\W)/(N-1)\to\tr\boldsymbol{\Psi} \)
in probability, so by Slutsky's lemma
\( F_N\to\norm{\mathbf{q}}^2/\tr\boldsymbol{\Psi} \) in distribution. On the other hand
\( F(m-1,\nu)\to\chi^2(m-1)/(m-1) \) in distribution as \( \nu\to\infty \), because its
denominator chi-squared divided by \( \nu \) tends to one in probability. If the
exactness holds for every \( N \), the two limits agree:
\[
\norm{\mathbf{q}}^2\ \overset{d}{=}\ \frac{\tr\boldsymbol{\Psi}}{m-1}\,\chi^2(m-1).
\]
Now match variances. By @thm-qf-mean-var,
\( \Var(\norm{\mathbf{q}}^2)=2\tr(\boldsymbol{\Psi}^2) \), while the right-hand side has
variance \( 2\{\tr\boldsymbol{\Psi}\}^2/(m-1) \). Equality says
\( (m-1)\tr(\boldsymbol{\Psi}^2)=\{\tr\boldsymbol{\Psi}\}^2 \), that is
\( \epsilon(\bSigma)=1 \), and @lem-cls-hf-form makes that sphericity.
:::

::: {.idea}
The repeated-measures \( F \) test does not require the occasions to be
exchangeable, only that the \( m-1 \) orthonormalized within-subject contrasts be
uncorrelated with equal variances — that the contrast space look spherical, which
series of measurements in time almost never do.
:::

## Two corrections and one alternative

When sphericity fails, \( \norm{\mathbf{q}}^2 \) and \( \tr(\W) \) are still weighted
sums of chi-squared variables, and matching their first two moments to scaled
chi-squared variables gives a usable approximation (Box 1954).

::: {#prp-cls-box}
[The moment-matched approximation]

Under the assumptions of @thm-cls-sphericity, with \( H \) true, write
\( \epsilon=\epsilon(\bSigma) \), \( d_1=m-1 \) and \( d_2=(N-1)(m-1) \). Then

::: {.enumerate options="label=(\alph*)"}
1. \( \E\,\text{SS}_T=\tr\boldsymbol{\Psi} \),
   \( \Var\,\text{SS}_T=2\tr(\boldsymbol{\Psi}^2) \),
   \( \E\,\text{SS}_E=(N-1)\tr\boldsymbol{\Psi} \) and
   \( \Var\,\text{SS}_E=2(N-1)\tr(\boldsymbol{\Psi}^2) \).

2. The scaled chi-squared variable \( g\chi^2(h) \) matching the first two moments of
   \( \text{SS}_T \) has \( h=\epsilon d_1 \) and
   \( g=\tr(\boldsymbol{\Psi}^2)/\tr\boldsymbol{\Psi} \); the one matching \( \text{SS}_E \)
   has \( h'=\epsilon d_2 \) and the *same* \( g \). Since \( \text{SS}_T \) and
   \( \text{SS}_E \) are independent by @thm-cls-sphericity(a), replacing each by its
   matched scaled chi-squared gives the approximation
   \( F\approx F(\epsilon d_1,\epsilon d_2) \), which is exact when \( \epsilon=1 \).

:::

:::

::: {.proof}
(a) The first two are @thm-qf-mean-var with \( \A=\I \) and mean zero. Writing
\( \W=\sum_{l=1}^{N-1}\bw_l\bw_l\T \) with \( \bw_l \) independent
\( \Normal_{m-1}(\bzero,\boldsymbol{\Psi}) \) (@lem-cor-centred-wishart),
\( \tr\W=\sum_l\norm{\bw_l}^2 \) is a sum of \( N-1 \) independent copies of
\( \norm{\mathbf{q}}^2 \), so its mean and variance are \( N-1 \) times theirs.

(b) Matching \( gh=\tr\boldsymbol{\Psi} \) and \( 2g^2h=2\tr(\boldsymbol{\Psi}^2) \)
gives \( g=\tr(\boldsymbol{\Psi}^2)/\tr\boldsymbol{\Psi} \) and
\( h=\{\tr\boldsymbol{\Psi}\}^2/\tr(\boldsymbol{\Psi}^2)=\epsilon d_1 \). For
\( \text{SS}_E \) both moments are multiplied by \( N-1 \), so \( g \) is unchanged and
\( h'=(N-1)h=\epsilon d_2 \). Since \( d_2/d_1=h'/h \), the ratio
\( (\text{SS}_T/d_1)/(\text{SS}_E/d_2) \) becomes
\( \{\chi^2(h)/h\}/\{\chi^2(h')/h'\} \), an \( F(h,h') \) variable, the two
chi-squared variables being independent by @thm-cls-sphericity(a); at
\( \epsilon=1 \) this is the exact result of @thm-cls-sphericity(b).


:::

By @lem-cls-hf-form, \( \epsilon\ge1/(m-1) \), so the approximating distribution
is never more dispersed than \( F(d_1/(m-1),d_2/(m-1))=F(1,N-1) \), the largest
that @prp-cls-box(b) can ask for. Referring \( F \) to \( F(1,N-1) \) is therefore
conservative within Box's approximation, and Geisser and Greenhouse (1958) proved
it exactly conservative, for every \( \bSigma \) and \( N \). This is the
**lower-bound** correction, turned into a recommendation by Greenhouse and
Geisser (1959).

Part (b) becomes a procedure once \( \epsilon \) is estimated from the sample
covariance matrix \( \bS \) of the \( \y_i \).

- The **Greenhouse–Geisser** correction refers \( F \) to
  \( F(\hat{\epsilon}d_1,\hat{\epsilon}d_2) \) with the plug-in
  \( \hat{\epsilon}=\epsilon(\bS) \), which is biased downward, so the test is
  conservative.
- The **Huynh–Feldt** correction removes the leading bias. For a single group of
  \( N \) subjects it uses
  \[
\tilde{\epsilon}=\min\Bigl\{1,\ \frac{N(m-1)\hat{\epsilon}-2}
   {(m-1)\bigl\{(N-1)-(m-1)\hat{\epsilon}\bigr\}}\Bigr\},
\]{#eq-cls-hf-epsilon}

  which exceeds \( \hat{\epsilon} \) and gives a less conservative test; its own
  bias is upward when \( \epsilon \) is small.
- The **lower-bound** correction takes \( \epsilon=1/(m-1) \): safe, and usually
  too conservative to be useful.

A third option avoids the issue. Since the \( \bz_i \) are independent
\( \Normal_{m-1}(\bU\T\bmu,\boldsymbol{\Psi}) \) with \( \boldsymbol{\Psi} \)
unrestricted, \( H \) says a multivariate normal mean is zero, and the classical
test is Hotelling's.

::: {#prp-cls-multivariate}
[The multivariate test]

Under the assumptions of @thm-cls-sphericity with \( N\ge m \), let
\( \bS_z=\W/(N-1) \) and
\[
T^2=N\,\bar{\bz}\T\bS_z^{-1}\bar{\bz}.
\]
Then \( T^2 \) does not depend on the choice of \( \bU \), nor on replacing \( \bU \) by
\( \bU\mathbf{K} \) for any nonsingular \( (m-1)\times(m-1) \) matrix \( \mathbf{K} \); and
under \( H \),
\[
\frac{N-m+1}{(N-1)(m-1)}\,T^2\sim F(m-1,\,N-m+1),
\]
whatever \( \bSigma \) is. The distributional statement is Hotelling's and is not
proved here; see Anderson (2003, §5.2).
:::

::: {.proof}
Only the invariance is proved. Replacing \( \bU \) by \( \bU\mathbf{K} \) replaces
\( \bar{\bz} \) by \( \mathbf{K}\T\bar{\bz} \) and \( \bS_z \) by
\( \mathbf{K}\T\bS_z\mathbf{K} \), leaving \( T^2 \) unchanged; two orthonormal bases
differ by an orthogonal \( \mathbf{K} \). The requirement \( N\ge m \) makes \( \bS_z \),
a Wishart matrix on \( N-1\ge m-1 \) degrees of freedom, invertible with probability
one.
:::

Exactness under any \( \bSigma \) sounds decisive, but it costs power: the
multivariate test estimates a whole covariance matrix and has \( N-m+1 \)
denominator degrees of freedom rather than \( (N-1)(m-1) \). When sphericity
nearly holds the univariate test is much more powerful, when it fails badly the
two are comparable, and when \( N<m \) the multivariate test does not exist.

## An example, four ways

::: {#exm-cls-grip}
[Grip strength over six visits]

Twenty-four patients have grip strength measured, in kilograms, at six weekly
visits, simulated from \( \Normal_6(\bmu,\bSigma) \) with \( \bSigma \) first-order
autoregressive, \( \sigma_{jk}=4\times0.6^{\lvert j-k\rvert} \), and a mean drifting
slowly upward; the observed occasion means are \( 47.818 \),
\( 48.099 \), \( 49.056 \), \( 49.177 \),
\( 49.561 \) and \( 49.263 \) kg. The true Box measure is
\( \epsilon=0.700 \): adjacent visits correlate \( 0.6 \) and visits five weeks
apart only \( 0.6^5 \), so adjacent differences are much less variable than
distant ones. The four analyses of \( H \) give

| Analysis | statistic | degrees of freedom | \( p \) |
|:---|---:|:---|---:|
| univariate, uncorrected | 4.770 | 5, 115 | 0.0005 |
| Greenhouse–Geisser, \( \hat{\epsilon}=0.645 \) | 4.770 | 3.23, 74.2 | 0.0035 |
| Huynh–Feldt, \( \tilde{\epsilon}=0.763 \) | 4.770 | 3.81, 87.7 | 0.0019 |
| multivariate | 3.605 | 5, 19 | 0.0184 |

All four reject, but the evidence differs by a factor of nearly forty between
first and last. Distrust the uncorrected test: its nominal degrees of freedom are
too large by the factor \( \epsilon \).
:::

```{.python .run #cell-repeated-measures-contrast}
import numpy as np
from scipy import stats

N, m = 24, 6
mu_true = np.array([48.0, 48.5, 48.9, 49.2, 49.4, 49.5])   # kg, a slow improvement
sigma, rho_true = 2.0, 0.6
lags = np.abs(np.subtract.outer(np.arange(m), np.arange(m)))
Sigma_true = sigma**2 * rho_true**lags               # AR(1) within a patient
rng = np.random.default_rng(3303)
Yobs = np.round(mu_true + rng.multivariate_normal(np.zeros(m), Sigma_true, N), 2)

def contrast_basis(m):
    """Orthonormal basis U of the space of contrasts among m repeated measures."""
    Q, _ = np.linalg.qr(np.column_stack([np.ones(m), np.eye(m)[:, 1:]]))
    return Q[:, 1:]                                  # columns orthonormal and orthogonal to 1

def epsilon(Psi):
    """Box's measure of departure from sphericity, 1/(m-1) <= eps <= 1."""
    k = Psi.shape[0]
    return np.trace(Psi)**2 / (k * np.sum(Psi * Psi))

U = contrast_basis(m)
Z = Yobs @ U                                         # N x (m-1) contrasts, one row per patient
zbar = Z.mean(axis=0)
Sz = np.cov(Z, rowvar=False, ddof=1)                 # Wishart, (N-1) degrees of freedom

F_uncorrected = N * (zbar @ zbar) / np.trace(Sz)
eps_hat = epsilon(Sz)
eps_hf = min(1.0, (N * (m - 1) * eps_hat - 2)
             / ((m - 1) * ((N - 1) - (m - 1) * eps_hat)))
d1, d2 = m - 1, (N - 1) * (m - 1)
p_un = stats.f.sf(F_uncorrected, d1, d2)
p_gg = stats.f.sf(F_uncorrected, d1 * eps_hat, d2 * eps_hat)
p_hf = stats.f.sf(F_uncorrected, d1 * eps_hf, d2 * eps_hf)
T2 = N * zbar @ np.linalg.solve(Sz, zbar)            # Hotelling's statistic
F_mv = (N - m + 1) / ((N - 1) * (m - 1)) * T2
p_mv = stats.f.sf(F_mv, m - 1, N - m + 1)
print(f"F = {F_uncorrected:.3f} on ({d1}, {d2}) df,  p = {p_un:.4f}")
print(f"Greenhouse-Geisser eps = {eps_hat:.3f}, p = {p_gg:.4f}")
print(f"Huynh-Feldt eps = {eps_hf:.3f}, p = {p_hf:.4f}")
print(f"multivariate F({m - 1}, {N - m + 1}) = {F_mv:.3f}, p = {p_mv:.4f}")
```

How much does the failure of sphericity cost? [Figure 33.3.1](#fig-cls-sphericity)
answers with a simulation under \( H \). Panel (a) shows the null distribution of
the uncorrected \( F \): more dispersed than \( F(5,115) \), with a heavier right
tail, and tracking @prp-cls-box closely. Panel (b) traces the level of a nominal
\( 5\% \) test as the autoregressive parameter grows. At \( \rho=0 \) it is
\( 0.050 \), as @thm-cls-sphericity(b) requires; by \( \rho=0.8 \),
where \( \epsilon=0.562 \), it is \( 0.085 \). The Greenhouse–Geisser
test is conservative throughout (\( 0.038 \) to \( 0.049 \)), the
Huynh–Feldt test close to nominal (\( 0.049 \) to
\( 0.056 \)), and the multivariate test exact at \( 0.049 \)
for every \( \rho \). A spherical covariance of the form @lem-cls-hf-form(iv) with
unequal \( a_j \), far from compound symmetric, gives the uncorrected test a level
of \( 0.050 \): sphericity, not compound symmetry, is what matters.

::: {when-format="html"}
![**Figure 33.3.1.** (a) The null distribution of the uncorrected \( F \) statistic
with autoregressive errors, against the nominal \( F(5,115) \) density and the
moment-matched \( F(\epsilon d_1,\epsilon d_2) \) of @prp-cls-box, with the nominal
\( 5\% \) point marked. (b) Simulated level of four tests as \( \rho \) grows.](sphericity.svg){#fig-cls-sphericity width=100%}
:::

::: {when-format="pdf"}
![(a) The null distribution of the uncorrected \( F \) statistic
with autoregressive errors, against the nominal \( F(5,115) \) density and the
moment-matched \( F(\epsilon d_1,\epsilon d_2) \) of @prp-cls-box, with the nominal
\( 5\% \) point marked. (b) Simulated level of four tests as \( \rho \) grows.](sphericity.pdf){width=100%}
:::

## The modern treatment

The corrections above patch a fixed analysis; the mixed-model approach fits the
covariance instead. Write the data as a linear mixed model (@def-mix-model) with a
fixed effect for each occasion and a within-subject covariance \( \R \) of chosen
shape, estimate its parameters by restricted maximum likelihood (@def-mix-reml),
and test the fixed effects by a Wald statistic with estimated degrees of
freedom (@prp-mix-inference). Compound symmetric \( \R \) reproduces the
uncorrected univariate \( F \) (@exr-cls-rm-mixed-equivalence), unstructured
\( \R \) the multivariate test up to the degrees-of-freedom approximation, and
anything between is a model the classical scheme has no room for. Unbalanced data
and missing visits are handled too, whereas @eq-cls-rm-ss and \( T^2 \) need a
complete rectangle. The price is that the covariance must be chosen, and the
choice affects the standard errors:
[Section 33.4](04-covariance-models.html).

::: {.remark}
[Between-subjects factors]

With a between-subjects factor the within-subject analysis is the same, applied to
the pooled within-group covariance, and the univariate tests are exact iff the
pooled \( \boldsymbol{\Psi} \) is spherical *and* the groups share a common
\( \bSigma \). Mauchly's (1940) test of sphericity itself, which most software
reports, is sensitive to nonnormality and, like any preliminary test, distorts
the level of what follows.
:::

## Exercises

### A. Check your understanding

::: {#exr-cls-sphericity-two}
[A1]

Show that every \( 2\times2 \) covariance matrix satisfies sphericity, so that with
\( m=2 \) the repeated-measures \( F \) test is exact whatever \( \bSigma \) is. Identify
the test in that case.
:::

::: {.solution}
With \( m=2 \), \( \boldsymbol{\Psi}=\bU\T\bSigma\bU \) is \( 1\times1 \), hence a multiple
of \( \I_1 \). By @thm-cls-sphericity(b) the test is exact. It is the paired \( t \)
test: the single contrast is \( (y_{i1}-y_{i2})/\sqrt2 \), and \( F \) is the square of
the paired \( t \) statistic on \( N-1 \) degrees of freedom.
:::

::: {#exr-cls-epsilon-ar1}
[A2]

For \( m=3 \) and an autoregressive covariance \( \sigma_{jk}=\rho^{\lvert j-k\rvert} \)
with \( 0<\rho<1 \), compute \( \boldsymbol{\Psi} \) for the basis
\( \bu_1=(1,-1,0)\T/\sqrt2 \), \( \bu_2=(1,1,-2)\T/\sqrt6 \) and show that sphericity
fails. Which pairwise difference has the largest variance?
:::

### B. Practice

::: {#exr-cls-cs-implies-sphericity}
[B1]

Show directly that compound symmetry
\( \bSigma=\sigma^2\{(1-\rho)\I+\rho\bone\bone\T\} \) gives
\( \boldsymbol{\Psi}=\sigma^2(1-\rho)\I_{m-1} \), and hence that the split-plot analysis
of @thm-cls-splitplot is the repeated-measures \( F \) test with \( \epsilon=1 \). Which
\( \lambda \) of @lem-cls-hf-form does it correspond to?
:::

::: {.solution}
\( \bU\T\bone=\bzero \), so
\( \boldsymbol{\Psi}=\bU\T\bSigma\bU=\sigma^2(1-\rho)\bU\T\bU=\sigma^2(1-\rho)\I_{m-1} \).
Hence \( \epsilon=1 \) and \( \lambda=\sigma^2(1-\rho) \), which is the subplot stratum
eigenvalue \( \lambda_s \) of @lem-cls-strata. In @lem-cls-hf-form(iv) take
\( \mathbf{a}=(\sigma^2\rho/2)\bone \).
:::

::: {#exr-cls-rm-mixed-equivalence}
[B2]

Consider the mixed model \( y_{ij}=\mu_j+u_i+e_{ij} \) with \( u_i \) and \( e_{ij} \)
independent normal, of variances \( \sigma_u^2 \) and \( \sigma_e^2 \). Show that
restricted maximum likelihood gives \( \hat{\sigma}_e^2=\text{MS}_E \) and
\( \hat{\sigma}_u^2=(\text{MS}_S-\text{MS}_E)/m \) when this is nonnegative, where
\( \text{MS}_S \) is the between-subject mean square, and that the Wald statistic for
\( H:\mu_1=\dots=\mu_m \), divided by \( m-1 \), is the uncorrected \( F \)
of @eq-cls-rm-ss.
:::

### C. Going deeper

::: {#exr-cls-power-comparison}
[C1]

For \( m=4 \), \( N=20 \) and \( \bSigma \) compound symmetric, compare by simulation
the power of the uncorrected univariate test, the Greenhouse–Geisser test and the
multivariate test of @prp-cls-multivariate against a linear trend in the occasion
means. Repeat with an autoregressive \( \bSigma \) of the same total variance.
Which would you choose if you knew the alternative were a linear trend, and why
is a single contrast better still?
:::
