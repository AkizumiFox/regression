# Spatial effects

Many regression data sets carry a location: the district a flat is in, the county a case
was reported from, the coordinates of a soil sample. Location is rarely of interest in
itself; it stands in for everything unmeasured that varies smoothly over space — climate,
soil, infrastructure, local labour markets — and leaving it out puts all of that into the
error and biases the coefficient of any spatially patterned covariate. Putting it in needs a
design block and a penalty, like any smooth covariate. The only question is what "smooth
over space" should mean, and there are two answers, one per kind of location data.

## The geoadditive model

::: {#def-add-spatial}
[Geoadditive model]

Let \( s_i \) be the location of observation \( i \). The **geoadditive model** is the
additive model @eq-add-model with one extra term,
\[
y_i=\x_{(i)}\T\bbeta+\sum_{j=1}^{q}f_j(z_{ij})+f_{\mathrm{geo}}(s_i)+\varepsilon_i ,
\qquad \sum_{i}f_{\mathrm{geo}}(s_i)=0 .
\]{#eq-add-geoadditive}

Two cases are distinguished by the kind of location.

::: {.enumerate options="label=(\alph*)"}
1. *(Areal data.)* \( s_i\in\{1,\dots,D\} \) names one of \( D \) regions. The design is
   the **incidence matrix** \( \Z_{\mathrm{geo}} \) with \( (\Z_{\mathrm{geo}})_{is}=1 \)
   if \( s_i=s \) and \( 0 \) otherwise, the coefficient \( \gamma_s \) is the effect of
   region \( s \), and the penalty is built from a **neighbourhood structure**: a
   symmetric relation \( s\sim t \) saying which regions share a border.

2. *(Point-referenced data.)* \( s_i\in\Real^2 \) is a coordinate pair. The design is
   \( (\Z_{\mathrm{geo}})_{il}=r(\norm{s_i-\kappa_l}) \) for a set of **knots**
   \( \kappa_1,\dots,\kappa_L \) and a radial function \( r \), and the penalty is
   \( \bP_{\mathrm{geo}}=\{r(\norm{\kappa_l-\kappa_m})\}_{lm} \).
:::

:::

In (a) the map is a graph and smoothness means that neighbouring regions do not differ by
much; in (b) the map is a plane and smoothness means what it means for a surface. Both are
pairs \( (\Z_{\mathrm{geo}},\bP_{\mathrm{geo}}) \), so
[Sections 44.1](01-additive-models.html) and [44.2](02-backfitting.html) apply without a
word of change.

## Markov random fields on a map

For areal data, write \( n_s \) for the number of neighbours of region \( s \) and define
the **neighbourhood penalty**
\[
\bgamma\T\bP\bgamma=\sum_{s\sim t}(\gamma_s-\gamma_t)^2,
\qquad
K_{st}=\begin{cases}
n_s, & t=s,\\
-1, & t\sim s,\\
0, & \text{otherwise,}
\end{cases}
\]{#eq-add-mrf}

where the sum runs over unordered neighbour pairs: the first-difference penalty of
@def-smo-pspline with "next knot" replaced by "neighbouring region", the two coinciding
exactly on a chain \( 1\sim2\sim\dots\sim D \).

::: {#exm-add-mrf-small}
[Four regions in a row]

Let four districts lie in a row, \( 1\sim2\sim3\sim4 \):
\[
\bP=\begin{pmatrix}
 1 & -1 & 0 & 0\\
-1 &  2 & -1 & 0\\
 0 & -1 &  2 & -1\\
 0 & 0 & -1 & 1
\end{pmatrix},
\]
with \( \bgamma\T\bP\bgamma=(\gamma_1-\gamma_2)^2+(\gamma_2-\gamma_3)^2+(\gamma_3-\gamma_4)^2 \),
rank \( 3 \) and null space \( \spn\{\bone\} \). If district \( 4 \) is instead an island
with no neighbours, its row and column vanish, the rank drops to \( 2 \) and the null space
is spanned by \( (1,1,1,0)\T \) and \( (0,0,0,1)\T \): the island is not penalized at all
and is estimated from its own data alone.
:::

::: {#prp-add-mrf}
[The neighbourhood penalty]

Let the regions \( 1,\dots,D \) carry a symmetric neighbourhood relation with no region its
own neighbour, and let \( \bP \) be as in @eq-add-mrf.

::: {.enumerate options="label=(\alph*)"}
1. \( \bP \) is symmetric and nonnegative definite, and
   \( \bgamma\T\bP\bgamma=\sum_{s\sim t}(\gamma_s-\gamma_t)^2 \).

2. \( \Null(\bP) \) is spanned by the indicator vectors of the connected components
   of the neighbourhood graph. In particular \( \rank(\bP)=D-c \) where \( c \) is
   the number of components, and \( \Null(\bP)=\spn\{\bone\} \) iff the graph is
   connected.

3. Treating \( \exp\{-\bgamma\T\bP\bgamma/(2\tau^2)\} \) as a density in
   \( \bgamma \) on \( \Null(\bP)\perpc \), the conditional distribution of
   \( \gamma_s \) given the others, *for a region with \( n_s\ge1 \)*, is normal with mean
   \( n_s^{-1}\sum_{t\sim s}\gamma_t \), the average of its neighbours, and variance
   \( \tau^2/n_s \). The prior is therefore a **Gaussian Markov random field**: each region
   is shrunk towards its neighbours, and a region with many neighbours is shrunk harder. A
   region with \( n_s=0 \) contributes nothing to the quadratic form: its coefficient is
   unpenalized, its conditional is flat, and it is estimated from its own data
   alone (@exm-add-mrf-small, @exr-add-island).
:::

:::

::: {.proof}
(a) Symmetry is immediate. Each unordered neighbour pair is counted twice in the double
sum over ordered pairs, so
\[
\sum_{s\sim t}(\gamma_s-\gamma_t)^2
=\tfrac12\sum_{s}\sum_{t\sim s}\bigl(\gamma_s^2-2\gamma_s\gamma_t+\gamma_t^2\bigr).
\]
The first inner sum gives \( \sum_sn_s\gamma_s^2 \), and so does the third, because
\( \sum_s\sum_{t\sim s}\gamma_t^2=\sum_tn_t\gamma_t^2 \) by symmetry of the relation.
Hence the right-hand side is
\( \sum_sn_s\gamma_s^2-\sum_s\sum_{t\sim s}\gamma_s\gamma_t=\bgamma\T\bP\bgamma \).
Nonnegativity follows because the left-hand side is a sum of squares.

(b) By (a), \( \bgamma\T\bP\bgamma=0 \) iff \( \gamma_s=\gamma_t \) across every
edge, iff \( \gamma \) is constant on each connected component. Since \( \bP \) is
nonnegative definite, that is equivalent to \( \bP\bgamma=\bzero \) by
@thm-mat-pd-characterizations, so the null space is the span of the component indicators,
which are independent because the components are disjoint.

(c) If \( n_s=0 \) then row \( s \) of \( \bP \) is zero and \( \gamma_s \) does not
appear in the quadratic form, which is the last sentence. For \( n_s\ge1 \), collect the
terms of \( \bgamma\T\bP\bgamma \) that involve \( \gamma_s \):
\[
n_s\gamma_s^2-2\gamma_s\sum_{t\sim s}\gamma_t+\text{const}
=n_s\Bigl(\gamma_s-\frac1{n_s}\sum_{t\sim s}\gamma_t\Bigr)^2+\text{const}' .
\]
Dividing by \( 2\tau^2 \) and exponentiating gives the normal density with the stated mean
and variance \( \tau^2/n_s \).
:::

Part (c) is the reason for the name. The penalty does not shrink a region towards zero —
that would be the ridge penalty \( \bP=\I \) — but towards its neighbourhood, so a
district with few observations borrows strength from those around it while the overall level
of the map stays free. That level is the one dimension the centring constraint of
@prp-add-identifiability(c) removes on a connected map; it belongs to the intercept.

::: {.remark}
[Variants]

Three variants fit the same framework. Weighting each pair by \( w_{st} \) (border length,
inverse centroid distance) replaces \( n_s \) by \( \sum_{t\sim s}w_{st} \) and \( -1 \)
by \( -w_{st} \). Adding \( \epsilon\I \) makes the penalty positive definite and the
prior proper, at the price of shrinking the map towards zero as well as towards itself. And
a *second*, unstructured term with \( \bP=\I \) on the same regions gives the
convolution model of Besag, York and Mollié (1991), separating smooth spatial variation from
region-specific noise.
:::

## Kriging bases for coordinates

For points rather than regions, smoothness comes from the theory of stationary random
fields, which is to say from [Part VII](../ch31-general-gauss-markov/index.html). Let the
spatial effect be a realization of a mean-zero process \( \{U(s)\} \) with
\( \Cov\{U(s),U(t)\}=\tau^2 r(\norm{s-t}) \): a covariance model in the sense of
@def-cls-covariance, on a continuous index set.

::: {#prp-add-kriging}
[Low-rank kriging is a penalized term]

Let \( \mathbf{R}=\{r(\norm{\kappa_l-\kappa_m})\}_{lm} \) be the correlation matrix of the
field at \( L \) knots, assumed positive definite, and let
\( \mathbf{u}=\{U(\kappa_1),\dots,U(\kappa_L)\}\T\sim\Normal_L(\bzero,\tau^2\mathbf{R}) \).

::: {.enumerate options="label=(\alph*)"}
1. The best linear predictor of \( U(s) \) from \( \mathbf{u} \) is
   \( \hat U(s)=\mathbf{r}(s)\T\mathbf{R}^{-1}\mathbf{u} \) with
   \( \mathbf{r}(s)=\{r(\norm{s-\kappa_1}),\dots,r(\norm{s-\kappa_L})\}\T \). Writing
   \( \bgamma=\mathbf{R}^{-1}\mathbf{u} \), the fitted spatial effect at the observations is
   \( \Z_{\mathrm{geo}}\bgamma \) with \( (\Z_{\mathrm{geo}})_{il}=r(\norm{s_i-\kappa_l}) \),
   and \( \bgamma\sim\Normal_L(\bzero,\tau^{2}\mathbf{R}^{-1}) \).

2. Consequently the log-density of \( \bgamma \) contributes
   \( \bgamma\T\mathbf{R}\bgamma/(2\tau^2) \), so the penalized criterion
   @eq-add-criterion with \( \bP_{\mathrm{geo}}=\mathbf{R} \) and
   \( \lambda_{\mathrm{geo}}=\sigma^2/\tau^2 \) has as its minimizer the posterior mode of
   \( \bgamma \) — the kriging predictor restricted to the span of the \( L \) knots.

3. Taking every distinct location as a knot recovers exact kriging, that is, generalized
   least squares with the covariance \( \tau^2\mathbf{R}+\sigma^2\I \) of
   [Chapter 31](../ch31-general-gauss-markov/index.html); taking \( L\ll n \) knots is the
   low-rank approximation that makes the term cheap.
:::

:::

::: {.proof}
(a) By @thm-cor-blup (or @thm-mix-blup with \( \X \) absent) the best linear predictor is
\( \Cov\{U(s),\mathbf{u}\}\Cov(\mathbf{u})^{-1}\mathbf{u}
=\tau^2\mathbf{r}(s)\T(\tau^2\mathbf{R})^{-1}\mathbf{u} \), the stated formula; evaluating
at \( s=s_i \) stacks into \( \Z_{\mathrm{geo}}\mathbf{R}^{-1}\mathbf{u} \), and the
distribution of \( \bgamma=\mathbf{R}^{-1}\mathbf{u} \) follows from @thm-rv-linear.
(b) \( -\log \) density of \( \Normal_L(\bzero,\tau^2\mathbf{R}^{-1}) \) is
\( \bgamma\T\mathbf{R}\bgamma/(2\tau^2) \) up to a constant. Adding it to
\( \norm{\y-\dots}^2/(2\sigma^2) \) and multiplying by \( 2\sigma^2 \) gives
@eq-add-criterion with the stated \( \lambda \).
(c) With \( L=n \) distinct locations, \( \Z_{\mathrm{geo}}=\mathbf{R} \) and the term
contributes \( \mathbf{R}\bgamma=\mathbf{u} \) with \( \mathbf{u} \) penalized by
\( \mathbf{u}\T\mathbf{R}^{-1}\mathbf{u} \); minimizing over \( \mathbf{u} \) is the mixed
model of @def-mix-model with \( \Z=\I \), \( \G=\tau^2\mathbf{R} \) and
\( \R=\sigma^2\I \), whose marginal covariance is
\( \tau^2\mathbf{R}+\sigma^2\I \) (@eq-mix-marginal).
:::

::: {.idea}
A spatial term is a random effect with a structured covariance — not an analogy: by
@prp-add-kriging the penalty *is* the inverse prior covariance, and
[Section 44.5](05-structured-additive.html) makes that general. The Markov random field
@eq-add-mrf says the same with a sparse precision matrix in place of a dense covariance,
which is why it is cheap on large maps.
:::

If \( \mathbf{R} \) and \( \Z_{\mathrm{geo}} \) are to be sparse, the correlation
function must have compact support or decay fast. The example below uses the Wendland
function \( r(h)=(1-h/\rho)_{+}^{4}(1+4h/\rho) \), positive definite in two dimensions and
zero beyond the range \( \rho \).

## A map

::: {#exm-add-map}
[Smoothing a district map]

A simulated country is a \( 12\times12 \) lattice of
\( 144 \) districts with between \( 1 \) and
\( 8 \) observations each, \( 604 \) in all. The
response has mean \( 2+0.7x+\text{surface}(s) \), with \( x \) an ordinary covariate and
the surface a smooth function of the district centroid that the model does not know; the
errors are normal with standard deviation \( 0.9 \).

[Figure 44.4.1](04-spatial-effects.html#fig-add-maps) compares three fits. Unpenalized
district effects use \( 143.0 \) degrees of freedom and give the noisy map
of panel (b), with root mean squared error \( 0.579 \) against the true
surface. The neighbourhood penalty @eq-add-mrf, with \( \lambda \) from the restricted
likelihood, uses \( 45.5 \) and brings that down to
\( 0.170 \). A kriging term on \( 21 \) knots uses
\( 17.7 \) and gives \( 0.228 \): smoother still, and
slightly too smooth, because the true surface has features the coarse knot grid cannot
represent. The coefficient of \( x \) comes out at \( 0.581 \) with a
standard error of \( 0.042 \) — against a true \( 0.7 \), about
\( 2.8 \) standard errors low on this draw — and the residual standard deviation at
\( 0.953 \).
:::

::: {when-format="html"}
![**Figure 44.4.1.** A simulated district map, all panels on one colour scale. (a) the surface used
to generate the data. (b) unpenalized district effects. (c) the Markov random field
@eq-add-mrf. (d) a kriging term on 21 knots.](spatial_maps.svg){#fig-add-maps width=100%}
:::

::: {when-format="pdf"}
![A simulated district map, all panels on one colour scale. (a) the surface used
to generate the data. (b) unpenalized district effects. (c) the Markov random field
@eq-add-mrf. (d) a kriging term on 21 knots.](spatial_maps.pdf){width=100%}
:::

The map and its penalty are built from the neighbourhood list alone:

```{.python .run #cell-spatial-lattice}
import numpy as np

def lattice(side):
    """Rook adjacency on a side x side grid of districts, and the district centroids."""
    coords = np.array([[c + 0.5, r + 0.5] for r in range(side) for c in range(side)])
    neighbours = []
    for s in range(side * side):
        r, c = divmod(s, side)
        nb = [(r + dr) * side + (c + dc) for dr, dc in ((-1, 0), (1, 0), (0, -1), (0, 1))
              if 0 <= r + dr < side and 0 <= c + dc < side]
        neighbours.append(sorted(nb))
    return coords, neighbours

def neighbourhood_penalty(neighbours):
    """K = D - W: the number of neighbours on the diagonal, -1 for each neighbour pair."""
    d = len(neighbours)
    K = np.zeros((d, d))
    for s, nbrs in enumerate(neighbours):
        K[s, s] = len(nbrs)
        for t in nbrs:
            K[s, t] = -1.0
    return K
```

::: {.warning}
[Spatial confounding]

A spatial term competes with every covariate that is itself spatially patterned: they are
concurve in the sense of @def-add-concurvity-index. In @exm-add-map the covariate \( x \)
was generated independently of location, so nothing was lost. When a covariate varies
smoothly over the map — income across districts, exposure across counties — the spatial term
can absorb a large part of its effect, and the fitted coefficient moves, sometimes
dramatically, when the term is added. Neither fit is "the" answer: without the spatial term
all spatial variation is attributed to the covariate, with it the smooth part goes to
location. Which is right is a question about unmeasured confounders, and the fit has no
opinion; [Chapter 25](../ch25-causal-interpretation/index.html) has the vocabulary. Report
both, and say which unmeasured causes the spatial term stands in for.
:::

## Exercises

### A. Check your understanding

::: {#exr-add-island}
[A1]

A map has two islands, each a connected group of districts, with no neighbour pairs between
them. What is \( \rank(\bP) \), what does the penalty leave unidentified, and what
does the centring constraint of @prp-add-identifiability(c) fix?
:::

### B. Practice

::: {#exr-add-mrf-chain-equals-difference}
[B1]

Show that on a chain \( 1\sim2\sim\dots\sim D \), the penalty @eq-add-mrf equals the
first-difference penalty \( \bgamma\T\mathbf{D}_1\T\mathbf{D}_1\bgamma \) of
@def-smo-pspline, and that the second-difference penalty is *not* of the form
@eq-add-mrf for any neighbourhood structure with nonnegative off-diagonal weights.
:::

::: {#exr-add-mrf-conditional}
[B2]

Using @prp-add-mrf(c), show that if district \( s \) has \( n_s \) neighbours and its own
\( m_s \) observations with mean \( \bar y_s \), and all other effects are held fixed, the
penalized criterion is minimized over \( \gamma_s \) at a weighted average of \( \bar y_s \)
and the neighbour mean, with weights \( m_s \) and \( \lambda n_s \). What happens when
\( m_s=0 \)?
:::

::: {.solution}
Up to constants the terms of @eq-add-criterion involving \( \gamma_s \) are
\( m_s(\bar y_s-\gamma_s)^2+\lambda n_s(\gamma_s-\bar\gamma_{(s)})^2 \) with
\( \bar\gamma_{(s)}=n_s^{-1}\sum_{t\sim s}\gamma_t \), so
\( \hat\gamma_s=(m_s\bar y_s+\lambda n_s\bar\gamma_{(s)})/(m_s+\lambda n_s) \). With
\( m_s=0 \) the district is set to its neighbour average: the map interpolates over the
gap.
:::

### C. Going deeper

::: {#exr-add-spatial-confounding}
[C1]

Let \( \bv \) be the vector of a covariate that is constant within each district, so that
\( \bv\in\C(\Z_{\mathrm{geo}}) \) for the incidence design. Show that in the *unpenalized*
geoadditive model the coefficient of \( \bv \) is not identified at all, and that with a
penalty it is identified but its estimate depends on \( \lambda_{\mathrm{geo}} \), tending
to the no-spatial-term estimate as \( \lambda_{\mathrm{geo}}\to\infty \). What does this say
about reporting such a coefficient?
:::
