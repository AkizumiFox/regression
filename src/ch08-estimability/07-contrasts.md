# Contrasts

In the one-way model the individual effects \( \alpha_k \) are not estimable, but their differences are. The
general form of such a difference is a contrast. Contrasts are the natural targets of inference about a factor.
They are estimable, they do not depend on the side condition or coding, and they are what a scientific question
about a factor usually asks: does treatment beat control, is there a trend across doses, do two groups of levels
differ? This section characterizes them, shows how an orthogonal set of contrasts splits the factor's sum of
squares, and constructs the Helmert and polynomial families.

## Definition and estimability

Throughout, the one-way model is \( \E(y_{kj})=\mu+\alpha_k \) for observation \( j \) at level \( k \), with \( n_k\ge1 \)
observations at level \( k=1,\dots,g \), level means \( \mu_k=\mu+\alpha_k \) and sample means \( \bar{y}_k \).

::: {#def-est-contrast}
[Contrast]

In a factor model, a **contrast** in the effects of a factor is a linear function \( \sum_kc_k\alpha_k \) whose coefficients
satisfy \( \sum_kc_k=0 \), not all zero. In the one-way model it equals \( \sum_kc_k\mu_k \), and we call either form a contrast in
the level means. Two contrasts are **orthogonal** (for the given design) if \( \sum_kc_kd_k/n_k=0 \).
:::

The two forms agree because \( \sum_kc_k\mu_k=\mu\sum_kc_k+\sum_kc_k\alpha_k \). This is where the condition
\( \sum_kc_k=0 \) comes from: it makes the unidentified \( \mu \) drop out. Simple differences \( \alpha_k-\alpha_l \),
comparisons of one level with the average of others such as \( \alpha_1-\tfrac12(\alpha_2+\alpha_3) \), and trend
coefficients are all contrasts.

::: {#thm-est-oneway-contrasts}
[Estimable functions in the one-way model]

In the one-way model with every \( n_k\ge1 \), the function \( \lambda_0\mu+\sum_k\lambda_k\alpha_k \) is estimable iff
\( \lambda_0=\sum_k\lambda_k \). In particular, a linear function of the effects alone, \( \sum_kc_k\alpha_k \), is estimable
iff it is a contrast (or zero). For a contrast, the least squares estimator and its variance under
\( \Cov(\Y)=\sigma^2\I \) are
\[
\hat{\psi}=\sum_kc_k\bar{y}_k,\qquad \Var(\hat{\psi})=\sigma^2\sum_k\frac{c_k^2}{n_k}.
\]{#eq-est-contrast-var}

:::

::: {.proof}
The null space of \( \X=[\bone,\Z] \) is spanned by \( \bv=(1,-1,\dots,-1)\T \) (@exm-proj-oneway-rank). By
@thm-est-characterization(c), the function is estimable iff \( \lambda_0-\sum_k\lambda_k=0 \). For an estimable function,
\( \lambda_0\mu+\sum_k\lambda_k\alpha_k=\sum_k\lambda_k(\mu+\alpha_k) \). Every least squares solution has
\( \hat{\mu}+\hat{\alpha}_k=\bar{y}_k \), since the fitted values are the level means (@exm-proj-oneway-M). So the estimate is
\( \sum_k\lambda_k\bar{y}_k \). The level means are uncorrelated with variances \( \sigma^2/n_k \), by @thm-rv-linear, which
gives the variance.
:::

The estimator in @eq-est-contrast-var is the same whatever coding or side condition was used to fit the model.
@prp-est-ls-estimator guarantees this, and it is the reason contrasts, rather than coefficients, should be
reported.

## Sums of squares and orthogonal contrasts

Each contrast defines a direction in the observation space, and its estimate is the component of \( \y \) along that
direction. With \( \bD=\diag(n_1,\dots,n_g) \) and \( \Z \) the indicator matrix, put
\[
\bu_{\mathbf{c}}=\Z\bD^{-1}\mathbf{c} ,
\]{#eq-est-contrast-vector}

the vector whose entries equal \( c_k/n_k \) for the observations at level \( k \).

::: {#prp-est-contrast-ss}
[Contrasts as projections]

Let \( \mathbf{c} \) and \( \mathbf{d} \) be contrasts, and write \( \hat{\psi}_{\mathbf{c}}=\mathbf{c}\T\hat{\bmu} \). Then:

::: {.enumerate options="label=(\alph*)"}
1. \( \hat{\psi}_{\mathbf{c}}=\bu_{\mathbf{c}}\T\y \), and \( \bu_{\mathbf{c}} \) is a nonzero vector in \( \C(\Z)\cap\bone\perpc \);

2. the **sum of squares for the contrast**,
    \[
    \operatorname{SS}(\mathbf{c})=\frac{\hat{\psi}_{\mathbf{c}}^2}{\sum_kc_k^2/n_k},
    \]
    equals \( \norm{\bP\y}^2 \), where \( \bP \) is the orthogonal projection onto the line spanned by \( \bu_{\mathbf{c}} \);

3. \( \Cov(\hat{\psi}_{\mathbf{c}},\hat{\psi}_{\mathbf{d}})=\sigma^2\sum_kc_kd_k/n_k=\sigma^2\bu_{\mathbf{c}}\T\bu_{\mathbf{d}} \) when \( \Cov(\Y)=\sigma^2\I \), so
           orthogonal contrasts have uncorrelated estimates and orthogonal vectors \( \bu \);

4. if \( \mathbf{c}_1,\dots,\mathbf{c}_{g-1} \) are pairwise orthogonal contrasts, then
    \[
    \sum_{j=1}^{g-1}\operatorname{SS}(\mathbf{c}_j)=\sum_{k=1}^gn_k(\bar{y}_k-\bar{y})^2 ,
    \]
    the sum of squares between levels.
:::

:::

::: {.proof}
(a) \( \hat{\bmu}=\bD^{-1}\Z\T\y \) (@eq-est-indicator-facts), so \( \mathbf{c}\T\hat{\bmu}=(\Z\bD^{-1}\mathbf{c})\T\y \). Clearly
\( \bu_{\mathbf{c}}\in\C(\Z) \), and \( \bone\T\bu_{\mathbf{c}}=\bone\T\Z\bD^{-1}\mathbf{c}=(n_1,\dots,n_g)\bD^{-1}\mathbf{c}=\sum_kc_k=0 \). It is nonzero because
\( \Z\bD^{-1} \) has full column rank and \( \mathbf{c}\ne\bzero \).
(b) By @exm-proj-line, \( \norm{\bP\y}^2=(\bu_{\mathbf{c}}\T\y)^2/\norm{\bu_{\mathbf{c}}}^2 \), and
\( \norm{\bu_{\mathbf{c}}}^2=\mathbf{c}\T\bD^{-1}\Z\T\Z\bD^{-1}\mathbf{c}=\mathbf{c}\T\bD^{-1}\mathbf{c}=\sum_kc_k^2/n_k \).
(c) @thm-rv-linear gives \( \Cov(\bu_{\mathbf{c}}\T\Y,\bu_{\mathbf{d}}\T\Y)=\sigma^2\bu_{\mathbf{c}}\T\bu_{\mathbf{d}} \), and the same computation as in (b)
gives \( \bu_{\mathbf{c}}\T\bu_{\mathbf{d}}=\mathbf{c}\T\bD^{-1}\mathbf{d} \).
(d) The space \( \C(\Z)\cap\bone\perpc \) has dimension \( g-1 \) (@thm-proj-nested(c), with \( \C(\bone)\subseteq\C(\Z) \)). By (a) and (c)
the vectors \( \bu_{\mathbf{c}_1},\dots,\bu_{\mathbf{c}_{g-1}} \) are \( g-1 \) nonzero, mutually orthogonal vectors in it, so they form an
orthogonal basis. By @thm-proj-sum the sum of the projections onto the lines they span is the projection onto that
space, which is \( \M-n^{-1}\bone\bone\T \) (@thm-proj-nested(b)). Applying it to \( \y \) and using Pythagoras gives
\( \sum_j\operatorname{SS}(\mathbf{c}_j)=\norm{(\M-n^{-1}\bone\bone\T)\y}^2 \). The vector \( (\M-n^{-1}\bone\bone\T)\y \) has entry
\( \bar{y}_k-\bar{y} \) at every observation of level \( k \).
:::

So a complete set of orthogonal contrasts decomposes the between-level sum of squares into \( g-1 \) one-dimensional
pieces. Under normality the pieces are independent, each is \( \sigma^2 \) times a noncentral \( \chi^2 \) with one degree
of freedom, and each is independent of the residual sum of squares (@thm-qf-orthogonal-projections). Each contrast can
therefore be tested on its own. [Chapter 9](../ch09-sums-of-squares/index.html) places this decomposition in the general
theory of orthogonal sums of squares, and Chapter 15 uses it for the one-way analysis of variance. In a *balanced* design,
with all \( n_k \) equal, orthogonality reduces to the ordinary \( \sum_kc_kd_k=0 \). In an unbalanced design, contrasts
that look orthogonal on paper need not be orthogonal for the data.

## Helmert contrasts

The Helmert matrix of @exm-mat-helmert has rows proportional to
\[
\mathbf{h}_k=(\underbrace{1,\dots,1}_{k-1},\,-(k-1),\,0,\dots,0),\qquad k=2,\dots,g .
\]
The contrast \( \mathbf{h}_k\T\bmu \) compares level \( k \) with the average of the levels before it, since it is \( (k-1) \) times
\( \tfrac1{k-1}\sum_{l<k}\mu_l-\mu_k \). The \( \mathbf{h}_k \) are mutually orthogonal in the ordinary sense, so in a balanced design they
are orthogonal contrasts. They suit factors whose levels are added in a meaningful order, such as successive refinements of a
treatment, where each new level is compared with everything before it.

## Polynomial contrasts

When the levels of a factor are ordered and equally spaced, such as doses \( 0,10,20,30 \) or the seven points of a party
scale, the natural questions are about the *shape* of the sequence of means. Is there a linear trend? Curvature? Give the levels scores \( s_k=k \). Apply Gram–Schmidt (@prp-proj-gram-schmidt)
to the vectors \( \bone,\mathbf{s},\mathbf{s}^{2},\dots,\mathbf{s}^{g-1} \) in \( \Real^g \), where \( \mathbf{s}^{j} \) has entries
\( s_k^j \). This gives orthogonal vectors \( \mathbf{p}_0,\dots,\mathbf{p}_{g-1} \), and \( \mathbf{p}_j \) is a polynomial of degree \( j \) in the score.
For \( j\ge1 \), \( \mathbf{p}_j\perp\bone \), so \( \mathbf{p}_j \) is a contrast. Rescaled to integers, the first three for \( g=7 \) are
\[
\text{linear }(-3,-2,-1,0,1,2,3),\quad
\text{quadratic }(5,0,-3,-4,-3,0,5),\quad
\text{cubic }(-1,1,1,0,-1,-1,1).
\]
By construction \( \mathbf{p}_j \) is orthogonal to every polynomial of degree less than \( j \). So if the level means follow a polynomial of degree
\( d \) in the score, every contrast of degree above \( d \) is zero. The linear contrast measures the trend, the quadratic the
curvature beyond it, and so on. This is the reason for their use.

With unequal group sizes the classical integer contrasts are no longer orthogonal for the design. The fix follows
@prp-est-contrast-ss(c). Orthogonality is needed in the inner product \( \sum_kc_kd_k/n_k \), which suggests setting
\( \mathbf{c}_j=\bD\mathbf{p}_j \), where now the \( \mathbf{p}_j \) are orthogonal polynomials in the weighted inner product \( \sum_kn_kp_kq_k \).
Then \( \sum_kc_{jk}=\sum_kn_kp_{jk}=0 \), and \( \sum_kc_{jk}c_{lk}/n_k=\sum_kn_kp_{jk}p_{lk}=0 \) for \( j\ne l \). In observation
space, \( \bu_{\mathbf{c}_j}=\Z\mathbf{p}_j \) is the polynomial of degree \( j \) in each observation's score, orthogonalized over the observations.
The weighted linear contrast gives exactly the regression sum of squares of \( \y \) on the score (@exr-est-linear-trend).

::: {#exm-est-party-trend}
[Trend in ideology across party identification]

Return to the seven party-identification groups of @exm-est-party-coding. The between-group sum of squares is
\( 772.87 \). The classical linear contrast has estimate \( \hat{\psi}=\sum_k(k-4)\bar{y}_k=11.175 \),
variance \( 0.1652\,\sigma^2 \) by @eq-est-contrast-var, and sum of squares
\( 755.85 \), which is \( 97.8 \) percent of the between-group sum of squares. The
quadratic and cubic contrasts have sums of squares \( 10.15 \) and \( 2.58 \). Ideology rises
almost linearly across the seven-point party scale, with some upward curvature
([Figure 8.7.1](#fig-est-trend)).

The group sizes range from \( 37 \) to \( 200 \), so the six classical contrasts are not orthogonal for this design. Their
sums of squares add to \( 784.02 \), which exceeds the between-group total. Part of the variation is counted twice.
The weighted orthogonal polynomial contrasts decompose the total exactly: \( 746.57 \) for the linear
trend and \( 26.30 \) for the five higher-degree contrasts together.
:::

::: {when-format="html"}
![**Figure 8.7.1.** Mean ideological self-placement by party identification (point area proportional to group size), with the
weighted least squares linear and quadratic trends in the score. The linear contrast accounts for almost all of the variation
between groups.](pid_trend.svg){#fig-est-trend width=62%}
:::

::: {when-format="pdf"}
![Mean ideological self-placement by party identification (point area proportional to group size), with the
weighted least squares linear and quadratic trends in the score. The linear contrast accounts for almost all of the variation
between groups.](pid_trend.pdf){width=62%}
:::

```{.python .run #cell-contrasts-data}
import numpy as np
import statsmodels.api as sm

anes = sm.datasets.anes96.load_pandas().data
y = anes["selfLR"].to_numpy()
level = anes["PID"].to_numpy().astype(int)       # 0 = strong Dem. ... 6 = strong Rep.
g = 7
n_k = np.bincount(level).astype(float)
means = np.bincount(level, weights=y) / n_k
```

```{.python .run #cell-contrasts-polynomial}
def orthogonal_polynomials(scores, weights):
    """Columns 1..g-1: polynomials in the scores, orthogonal in sum_k w_k u_k v_k."""
    V = np.vander(scores, len(scores), increasing=True)      # 1, s, s^2, ...
    W = np.sqrt(weights)[:, None]
    Q, _ = np.linalg.qr(W * V)                               # weighted Gram-Schmidt
    P = Q / W                                                # back to the original scale
    P = P / np.abs(P).max(axis=0)                            # cosmetic rescaling
    return P[:, 1:]                                          # drop the constant column

def contrast_ss(c, means, n_k):
    """Estimate c'mu-hat and its sum of squares (c'mu-hat)^2 / sum(c_k^2 / n_k)."""
    psi = c @ means
    return psi, psi ** 2 / np.sum(c ** 2 / n_k)

scores = np.arange(g, dtype=float)
C_unw = orthogonal_polynomials(scores, np.ones(g))           # the classical table (equal n)
C_wtd = n_k[:, None] * orthogonal_polynomials(scores, n_k)   # c_k = n_k p(k)
ss_between = np.sum(n_k * (means - y.mean()) ** 2)
for name, C in [("unweighted", C_unw), ("weighted", C_wtd)]:
    ss = [contrast_ss(C[:, j], means, n_k)[1] for j in range(g - 1)]
    print(f"{name:10s} SS by degree: {np.round(ss, 2)}  total {sum(ss):.2f}",
          f"(between groups {ss_between:.2f})")
```

## Coding matrices are not contrast matrices

[Section 8.5](05-factor-coding.html) coded a factor by an intercept and the columns \( \Z\mathbf{C} \). When the columns of \( \mathbf{C} \) are
contrasts, as for Helmert or polynomial coding, it is tempting to read the coefficient of column \( j \) as the estimate of the
contrast \( \mathbf{c}_j\T\bmu \). That reading is wrong in general. By @prp-est-coding(b) the coefficients are the rows of
\( \mathbf{K}^{-1} \), with \( \mathbf{K}=[\bone,\mathbf{C}] \), applied to \( \hat{\bmu} \). The rows of \( \mathbf{K}^{-1} \) and the columns of \( \mathbf{K} \) are related
simply only when the columns are orthogonal.

::: {#prp-est-coding-dual}
[Orthogonal codings]

If the columns of \( \mathbf{K}=[\bone_g,\mathbf{c}_1,\dots,\mathbf{c}_{g-1}] \) are mutually orthogonal, then the coefficients of the coded model are
\[
\hat{\gamma}_0=\frac1g\sum_k\bar{y}_k,\qquad \hat{\gamma}_j=\frac{\mathbf{c}_j\T\hat{\bmu}}{\mathbf{c}_j\T\mathbf{c}_j},\quad j=1,\dots,g-1 ,
\]
whatever the group sizes.
:::

::: {.proof}
With orthogonal columns, \( \mathbf{K}\T\mathbf{K}=\diag(g,\mathbf{c}_1\T\mathbf{c}_1,\dots,\mathbf{c}_{g-1}\T\mathbf{c}_{g-1}) \), so
\( \mathbf{K}^{-1}=(\mathbf{K}\T\mathbf{K})^{-1}\mathbf{K}\T \). Apply @prp-est-coding(b).
:::

The orthogonality required here is the plain one, \( \mathbf{c}_j\T\mathbf{c}_l=0 \), not the design orthogonality of
@def-est-contrast, because \( \hat{\bgamma}=\mathbf{K}^{-1}\hat{\bmu} \) is pure algebra and does not involve the \( n_k \). Helmert coding, in the usual software
convention, uses the columns \( \mathbf{c}_j=-\mathbf{h}_{j+1} \), which compare level \( j+1 \) with the levels before it, and
\( \mathbf{c}_j\T\mathbf{c}_j=j(j+1) \). In the party data the first three Helmert coefficients are \( 0.1408 \),
\( 0.0698 \) and \( 0.1271 \), while the contrasts themselves are
\( 0.2817 \), \( 0.4187 \) and \( 1.5253 \), larger by the factors \( 2 \),
\( 6 \) and \( 12 \). For a non-orthogonal coding there is no such shortcut. Reference coding has
\( \mathbf{c}_j=\mathbf{e}_{j+1} \), yet its coefficients are \( \mu_{j+1}-\mu_1 \), not \( \mu_{j+1} \).

```{.python .run #cell-contrasts-helmert}
def helmert(g):
    """Coding matrix whose column k compares level k+1 with the mean of levels 1..k."""
    C = np.zeros((g, g - 1))
    for k in range(1, g):
        C[:k, k - 1] = -1.0
        C[k, k - 1] = k
    return C

H = helmert(g)
K = np.column_stack([np.ones(g), H])
gamma = np.linalg.solve(K, means)                  # Helmert-coded coefficients
psi = H.T @ means                                  # the contrasts c_k' mu-hat
print("coefficients:", np.round(gamma[1:], 4))
print("contrasts:   ", np.round(psi, 4))
print("ratio:       ", np.round(psi / gamma[1:], 1))   # = c_k' c_k = k(k+1)
```

::: {.warning}
A coefficient from a contrast-coded regression estimates a contrast only up to a scale factor, and only when the coding
columns are orthogonal. To test or report a contrast, compute \( \mathbf{c}\T\hat{\bmu} \) and its variance @eq-est-contrast-var directly,
or read it from the coded fit using \( \mathbf{K}^{-1} \). Do not read it off the coefficient table.
:::

## Contrasts in the two-way additive model

In the additive model \( \mu+\alpha_i+\beta_j \) of [Section 8.6](06-several-factors.html), a **row contrast** is
\( \sum_ic_i\alpha_i \) with \( \sum_ic_i=0 \). It is a combination of the differences \( \alpha_i-\alpha_1 \), so in a connected
design it is estimable (@thm-est-connected). In a balanced design, with \( m\ge1 \) observations in every one of the \( ab \)
cells, the least squares fit is \( \bar{y}_{i\cdot\cdot}+\bar{y}_{\cdot j\cdot}-\bar{y}_{\cdot\cdot\cdot} \). This is shown in
@exm-proj-two-factor for \( m=1 \), and the averaging argument there is unchanged for \( m>1 \). Hence
\[
\widehat{\textstyle\sum_ic_i\alpha_i}=\sum_ic_i\bigl(\bar{y}_{i\cdot\cdot}+\bar{y}_{\cdot j\cdot}-\bar{y}_{\cdot\cdot\cdot}\bigr)=\sum_ic_i\bar{y}_{i\cdot\cdot},
\qquad \Var=\frac{\sigma^2}{bm}\sum_ic_i^2 ,
\]
since the terms not depending on \( i \) are multiplied by \( \sum_ic_i=0 \). The row contrast is estimated from the row means alone,
as if factor \( B \) were absent. In an unbalanced design this is false. The row means are then contaminated by the column effects
(@exr-est-unbalanced-rows), and the least squares estimate has to adjust for them. That adjustment is the subject of
Chapter 17. Interaction contrasts, and their estimability when cells are empty, were treated in @prp-est-interaction.

## Exercises

### A. Check your understanding

::: {#exr-est-poly4}
[A1]

Find the linear, quadratic and cubic orthogonal polynomial contrasts for \( g=4 \) equally spaced levels, scaled to integers.
Check that they are mutually orthogonal and that the quadratic contrast vanishes when \( \mu_k=a+bk \).
:::

::: {#exr-est-helmert-balance}
[A2]

Show that the Helmert contrasts \( \mathbf{h}_2,\dots,\mathbf{h}_g \) are orthogonal for the design iff \( n_1=n_2=\dots=n_{g-1} \), whatever
\( n_g \) is. *Hint:* compute \( \sum_kh_{jk}h_{lk}/n_k \) for \( j<l \).
:::

### B. Practice

::: {#exr-est-max-contrast}
[B1]

Show that for every contrast \( \mathbf{c} \), \( \operatorname{SS}(\mathbf{c})\le\sum_kn_k(\bar{y}_k-\bar{y})^2 \), with equality when
\( c_k=n_k(\bar{y}_k-\bar{y}) \). This fact is the basis of Scheffé's method in [Chapter 13](../ch13-multiplicity/index.html) (@thm-mc-scheffe).
:::

::: {.solution}
By @prp-est-contrast-ss(b), \( \operatorname{SS}(\mathbf{c})=\norm{\bP\y}^2 \) with \( \bP \) projecting onto a line inside
\( \mathcal S=\C(\Z)\cap\bone\perpc \). Let \( \bP_{\mathcal S}=\M-n^{-1}\bone\bone\T \). Because the line lies in \( \mathcal S \),
\( \bP\y=\bP\bP_{\mathcal S}\y \) (@thm-proj-nested(a)), and so
\( \norm{\bP\y}^2\le\norm{\bP_{\mathcal S}\y}^2=\sum_kn_k(\bar{y}_k-\bar{y})^2 \). Equality holds when the line contains
\( \bP_{\mathcal S}\y \), whose entries are \( \bar{y}_k-\bar{y} \) at level \( k \). That vector is \( \bu_{\mathbf{c}} \) for
\( c_k=n_k(\bar{y}_k-\bar{y}) \), and this \( \mathbf{c} \) is a contrast because \( \sum_kn_k(\bar{y}_k-\bar{y})=0 \).
:::

::: {#exr-est-unbalanced-rows}
[B2]

In the additive two-way model with \( n_{ij} \) observations in cell \( (i,j) \), show that
\( \E(\bar{y}_{i\cdot\cdot})=\mu+\alpha_i+\sum_jn_{ij}\beta_j/n_{i\cdot} \). Deduce that \( \sum_ic_i\bar{y}_{i\cdot\cdot} \) is unbiased for
every row contrast and every \( \bbeta \) iff the proportions \( n_{ij}/n_{i\cdot} \) do not depend on \( i \), that is, iff the cell counts are
*proportional*, \( n_{ij}=n_{i\cdot}n_{\cdot j}/n \).
:::

### C. Going deeper

::: {#exr-est-linear-trend}
[C1]

Let \( \x \) be the \( n \)-vector giving each observation's score \( s_k \). Show that the sum of squares of the weighted linear contrast
\( \mathbf{c}_1=\bD\mathbf{p}_1 \) equals the regression sum of squares \( \norm{(\M_1-n^{-1}\bone\bone\T)\y}^2 \) of the simple regression of \( \y \)
on \( \x \) with intercept, where \( \M_1 \) projects onto \( \C([\bone,\x]) \). Conclude that the remaining weighted contrasts together
measure the lack of fit of the straight line to the level means.
:::

::: {.solution}
In the weighted inner product, \( \mathbf{p}_1 \) is \( \mathbf{s} \) minus its
weighted mean, \( \mathbf{s}-\bar{s}_w\bone \) with \( \bar{s}_w=\sum_kn_ks_k/n \), up to scale. Then \( \bu_{\mathbf{c}_1}=\Z\bD^{-1}\bD\mathbf{p}_1=\Z\mathbf{p}_1 \), whose
entries are \( s_k-\bar{s}_w \). That is \( \x-\bar{x}\bone \), the centred score vector. By @prp-est-contrast-ss(b),
\( \operatorname{SS}(\mathbf{c}_1) \) is the squared length of the projection of \( \y \) onto the line spanned by \( \x-\bar{x}\bone \). By
@exm-proj-simple-nested this line is \( \C([\bone,\x])\cap\bone\perpc \), and the projection onto it is \( \M_1-n^{-1}\bone\bone\T \). The
remaining weighted contrasts span the orthogonal complement of that line in \( \C(\Z)\cap\bone\perpc \), which is
\( \C(\Z)\cap\C([\bone,\x])\perpc \). Their total sum of squares is \( \norm{(\M-\M_1)\y}^2 \), the reduction in residual sum of squares
from the straight-line model to the model with a separate mean for each level.
:::
