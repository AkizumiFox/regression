# Balanced two-way layouts via orthogonal projections

With \( m\ge2 \) observations in every cell, the replicates supply an estimate of \( \sigma^2 \) that is valid whatever the cell means are, and interaction can be tested without any assumption about its form. This section gives the complete analysis of the balanced two-way layout: the orthogonal decomposition, exact tests, contrasts, and simultaneous inference. @exr-proj-two-factor-interaction asked for the five projections, and @thm-tw-additive(d) has already met four of them.

## The decomposition

The model is the cell-means model
\[
\E(y_{ijk})=\mu_{ij},\qquad \Cov(\Y)=\sigma^2\I,\qquad i\le a,\ j\le b,\ k\le m,
\]{#eq-tw-cell-means}

written in the effects form of @def-tw-interaction as \( \mu_{ij}=\mu+\alpha_i+\beta_j+\gamma_{ij} \), with equal-weight effects. Its mean space \( \C(\X) \) is the set of vectors that are constant within cells, of dimension \( ab \) (@prp-est-interaction). With the ordering and notation of [Section 16.1](01-additive.html), define
\[
\begin{aligned}
\bP_0&=\bar{\mathbf{J}}_a\otimes\bar{\mathbf{J}}_b\otimes\bar{\mathbf{J}}_m, &
\bP_A&=\mathbf{C}_a\otimes\bar{\mathbf{J}}_b\otimes\bar{\mathbf{J}}_m, &
\bP_B&=\bar{\mathbf{J}}_a\otimes\mathbf{C}_b\otimes\bar{\mathbf{J}}_m,\\
\bP_{AB}&=\mathbf{C}_a\otimes\mathbf{C}_b\otimes\bar{\mathbf{J}}_m, &
\bP_E&=\I_a\otimes\I_b\otimes\mathbf{C}_m .
\end{aligned}
\]{#eq-tw-five-projections}

::: {#thm-tw-balanced}
[The balanced two-way layout]

Assume @eq-tw-cell-means with \( a,b\ge2 \) and \( m\ge2 \), and let \( \alpha_i \), \( \beta_j \), \( \gamma_{ij} \) be the equal-weight effects of @def-tw-interaction.

::: {.enumerate options="label=(\alph*)"}
1. The five matrices @eq-tw-five-projections are mutually orthogonal projections that sum to \( \I \), with ranks \( 1 \), \( a-1 \), \( b-1 \), \( (a-1)(b-1) \) and \( ab(m-1) \). The projection onto \( \C(\X) \) is \( \M=\bP_0+\bP_A+\bP_B+\bP_{AB} \), and that onto the additive model space is \( \bP_0+\bP_A+\bP_B \). So \( \C(\bP_{AB}) \) is the orthogonal complement of the additive model within the cell-means model.

2. The projections act by
   \[
\begin{aligned}
(\bP_A\y)_{ijk}&=\bar{y}_{i\cdot\cdot}-\bar{y}_{\cdot\cdot\cdot},\\
(\bP_{AB}\y)_{ijk}&=\bar{y}_{ij\cdot}-\bar{y}_{i\cdot\cdot}-\bar{y}_{\cdot j\cdot}+\bar{y}_{\cdot\cdot\cdot},\\
(\bP_E\y)_{ijk}&=y_{ijk}-\bar{y}_{ij\cdot},
\end{aligned}
\]
   and symmetrically for \( \bP_B \). The sums of squares are
   \[
\begin{aligned}
\text{SS}_A&=bm\sum_i(\bar{y}_{i\cdot\cdot}-\bar{y}_{\cdot\cdot\cdot})^2,\\
\text{SS}_B&=am\sum_j(\bar{y}_{\cdot j\cdot}-\bar{y}_{\cdot\cdot\cdot})^2,\\
\text{SS}_{AB}&=m\sum_{i,j}(\bar{y}_{ij\cdot}-\bar{y}_{i\cdot\cdot}-\bar{y}_{\cdot j\cdot}+\bar{y}_{\cdot\cdot\cdot})^2,\\
\text{SSE}&=\sum_{i,j,k}(y_{ijk}-\bar{y}_{ij\cdot})^2 .
\end{aligned}
\]

3. The mean squares have expectations
   \[
\begin{aligned}
\E\,\text{MS}_A&=\sigma^2+\frac{bm}{a-1}\sum_i\alpha_i^2,\\
\E\,\text{MS}_B&=\sigma^2+\frac{am}{b-1}\sum_j\beta_j^2,\\
\E\,\text{MS}_{AB}&=\sigma^2+\frac{m}{(a-1)(b-1)}\sum_{i,j}\gamma_{ij}^2,\\
\E\,\text{MSE}&=\sigma^2 .
\end{aligned}
\]

4. Under normal errors the four sums of squares are independent, \( \text{SSE}/\sigma^2\sim\chi^2(ab(m-1)) \), and
   \[
\begin{aligned}
F_A&=\frac{\text{MS}_A}{\text{MSE}}\sim F\Bigl(a-1,\,ab(m-1),\,\frac{bm}{\sigma^2}\sum_i\alpha_i^2\Bigr),\\
F_{AB}&=\frac{\text{MS}_{AB}}{\text{MSE}}\sim F\Bigl((a-1)(b-1),\,ab(m-1),\,\frac m{\sigma^2}\sum_{i,j}\gamma_{ij}^2\Bigr),
\end{aligned}
\]
   and likewise for \( F_B \). So \( F_{AB} \) is an exact test of additivity, \( F_A \) of \( \bar{\mu}_{1\cdot}=\dots=\bar{\mu}_{a\cdot} \), and \( F_B \) of \( \bar{\mu}_{\cdot1}=\dots=\bar{\mu}_{\cdot b} \).

5. \( \text{SS}_A \) is the sum of squares of the hypothesis \( \alpha_1=\dots=\alpha_a=0 \) in the cell-means model. It is also the sequential sum of squares of \( A \) whether \( A \) is entered before or after \( B \), and it equals the Type II and Type III sums of squares of \( A \). The same holds for \( B \), and \( \text{SS}_{AB} \) is the sequential, Type II and Type III sum of squares of the interaction.
:::

:::

::: {.proof}
(a) Orthogonality, ranks and the sum \( \I \) follow from @lem-tw-kron, as in the proof of @thm-tw-additive(d). The sum of the first four is \( (\bar{\mathbf{J}}_a+\mathbf{C}_a)\otimes(\bar{\mathbf{J}}_b+\mathbf{C}_b)\otimes\bar{\mathbf{J}}_m=\I_a\otimes\I_b\otimes\bar{\mathbf{J}}_m \). This matrix replaces each observation by its cell mean, so it is the projection onto the vectors that are constant within cells, which is \( \C(\X) \). The additive projection is @thm-tw-additive(a).

(b) Expand each Kronecker product with \( \mathbf{C}=\I-\bar{\mathbf{J}} \) and read off the means from the table in [Section 16.1](01-additive.html). For example,
\[
\begin{aligned}
\bP_{AB}={}&\I_a\otimes\I_b\otimes\bar{\mathbf{J}}_m-\I_a\otimes\bar{\mathbf{J}}_b\otimes\bar{\mathbf{J}}_m\\
&-\bar{\mathbf{J}}_a\otimes\I_b\otimes\bar{\mathbf{J}}_m+\bar{\mathbf{J}}_a\otimes\bar{\mathbf{J}}_b\otimes\bar{\mathbf{J}}_m .
\end{aligned}
\]
Each sum of squares is a sum over the \( n \) positions of a quantity that is constant over \( m \), \( bm \) or \( am \) of them.

(c) Apply (b) to the mean vector \( \bmu \), whose entries are \( \mu_{ij} \). Then \( \bP_A\bmu \) has entries \( \bar{\mu}_{i\cdot}-\bar{\mu}_{\cdot\cdot}=\alpha_i \), so \( \norm{\bP_A\bmu}^2=bm\sum_i\alpha_i^2 \). Likewise \( \norm{\bP_{AB}\bmu}^2=m\sum_{i,j}\gamma_{ij}^2 \), and \( \bP_E\bmu=\bzero \). The four subspaces \( \C(\bP_0),\dots,\C(\bP_{AB}) \) form an analysis of variance decomposition of \( \C(\X) \), so @thm-ss-expected-mean-squares gives the expectations.

(d) This is @cor-ss-ms-distributions with the noncentralities of (c). The noncentrality of \( F_A \) vanishes iff every \( \alpha_i=0 \), that is, iff the row means \( \bar{\mu}_{i\cdot} \) are equal. The noncentrality of \( F_{AB} \) vanishes iff every \( \gamma_{ij}=0 \), which is additivity by @prp-tw-interaction-equiv(vi).

(e) Each \( \alpha_i \) is a linear function of the cell means, so the hypothesis is testable. The mean vectors it allows are the \( \bmu\in\C(\X) \) with \( \bP_A\bmu=\bzero \). Since \( \bmu=\bP_0\bmu+\bP_A\bmu+\bP_B\bmu+\bP_{AB}\bmu \), these form \( \mathcal S_0=\C(\bP_0)\dirsum\C(\bP_B)\dirsum\C(\bP_{AB}) \). By @thm-glh-general-f(b) the sum of squares of the hypothesis is \( \norm{(\M-\bP_{\mathcal S_0})\y}^2=\norm{\bP_A\y}^2 \). The rest is Chapter 9: equal counts are proportional (@thm-ss-proportional), so @thm-ss-orthogonal-design applies, and @exr-ss-balanced-types shows that the sequential sums of squares of \( A \) in either order and its Type II and Type III sums of squares all equal \( \norm{\bP_A\y}^2 \). The interaction enters last in every order, so all its sums of squares are \( \norm{(\M-\bP_0-\bP_A-\bP_B)\y}^2=\norm{\bP_{AB}\y}^2 \).
:::

The theorem is summarized in the analysis of variance table:

| Source | df | Sum of squares | \( \E(\text{MS}) \) |
|---|---|---|---|
| \( A \) | \( a-1 \) | \( \norm{\bP_A\y}^2 \) | \( \sigma^2+bm\sum_i\alpha_i^2/(a-1) \) |
| \( B \) | \( b-1 \) | \( \norm{\bP_B\y}^2 \) | \( \sigma^2+am\sum_j\beta_j^2/(b-1) \) |
| \( AB \) | \( (a-1)(b-1) \) | \( \norm{\bP_{AB}\y}^2 \) | \( \sigma^2+m\sum_{i,j}\gamma_{ij}^2/((a-1)(b-1)) \) |
| error | \( ab(m-1) \) | \( \norm{\bP_E\y}^2 \) | \( \sigma^2 \) |
| corrected total | \( abm-1 \) | \( \norm{(\I-\bP_0)\y}^2 \) | |

The test for \( A \) is valid with or without interaction, because \( \C(\bP_{AB})\perp\C(\bP_A) \). What the interaction changes is the *meaning* of \( H_A \): equality of equal-weight row means is one hypothesis among many when the simple effects differ ([Section 16.2](02-interaction.html)). Part (e) is special to balance. With unequal cell counts the sequential, Type II and Type III sums of squares of \( A \) test hypotheses with different weights (@thm-ss-two-way-hypotheses), the subject of [Chapter 17](../ch17-unbalanced-data/index.html).

## Contrasts and their sums of squares

The \( F \) tests say whether a term matters. Contrasts say how. A **main-effect contrast** in \( A \) is \( \psi=\sum_ic_i\bar{\mu}_{i\cdot}=\sum_ic_i\alpha_i \) with \( \sum_ic_i=0 \). An **interaction contrast** is \( \sum_{i,j}d_{ij}\mu_{ij}=\sum_{i,j}d_{ij}\gamma_{ij} \), with \( (d_{ij}) \) having zero row and column sums (@prp-tw-interaction-equiv). Their least squares estimators replace means by sample means. Since \( \bar{y}_{i\cdot\cdot} \) averages \( bm \) independent observations and \( \bar{y}_{ij\cdot} \) averages \( m \),
\[
\Var\Bigl(\sum_ic_i\bar{y}_{i\cdot\cdot}\Bigr)=\frac{\sigma^2}{bm}\sum_ic_i^2,\qquad
\Var\Bigl(\sum_{i,j}d_{ij}\bar{y}_{ij\cdot}\Bigr)=\frac{\sigma^2}{m}\sum_{i,j}d_{ij}^2 .
\]
The sums of squares of the corresponding single-degree-of-freedom hypotheses are the squared estimates divided by these variances without the \( \sigma^2 \) (@thm-ss-single-df), as for a one-way contrast (@prp-aov-contrast-ss):
\[
\text{SS}(\mathbf{c})=\frac{\bigl(\sum_ic_i\bar{y}_{i\cdot\cdot}\bigr)^2}{\sum_ic_i^2/(bm)},\qquad
\text{SS}(\mathbf{d})=\frac{\bigl(\sum_{i,j}d_{ij}\bar{y}_{ij\cdot}\bigr)^2}{\sum_{i,j}d_{ij}^2/m}.
\]{#eq-tw-contrast-ss}

For main effects, a complete set of \( a-1 \) mutually orthogonal contrasts (\( \sum_ic_ic'_i=0 \)) splits \( \text{SS}_A \) into \( a-1 \) independent pieces. This is the one-way result @thm-aov-orthogonal-contrasts applied to the row means, which are one-way means with \( bm \) observations each. For interaction, orthogonal contrasts are most easily built as products.

::: {#prp-tw-product-contrasts}
[Product contrasts]

Let \( \mathbf{c}^{(1)},\dots,\mathbf{c}^{(a-1)} \) be mutually orthogonal contrast vectors in \( \Real^a \), and \( \mathbf{e}^{(1)},\dots,\mathbf{e}^{(b-1)} \) mutually orthogonal contrast vectors in \( \Real^b \). Then the \( (a-1)(b-1) \) tables \( \mathbf{d}^{(rs)}=\mathbf{c}^{(r)}\mathbf{e}^{(s)\top} \), with entries \( c^{(r)}_ie^{(s)}_j \), are coefficient tables of mutually orthogonal interaction contrasts, and
\[
\sum_{r=1}^{a-1}\sum_{s=1}^{b-1}\text{SS}\bigl(\mathbf{d}^{(rs)}\bigr)=\text{SS}_{AB}.
\]
Under normality the \( (a-1)(b-1) \) sums of squares are independent of each other and of \( \text{SSE} \).
:::

::: {.proof}
The row sums of \( \mathbf{d}^{(rs)} \) are \( c^{(r)}_i\sum_je^{(s)}_j=0 \) and the column sums vanish likewise. Let \( \bv_{rs}=\mathbf{c}^{(r)}\otimes\mathbf{e}^{(s)}\otimes\bone_m\in\Real^n \). Since \( \mathbf{C}_a\mathbf{c}^{(r)}=\mathbf{c}^{(r)} \) and \( \mathbf{C}_b\mathbf{e}^{(s)}=\mathbf{e}^{(s)} \), the mixed-product rule gives \( \bP_{AB}\bv_{rs}=\bv_{rs} \), so \( \bv_{rs}\in\C(\bP_{AB}) \). Also \( \bv_{rs}\T\bv_{r's'}=m\,(\mathbf{c}^{(r)\top}\mathbf{c}^{(r')})(\mathbf{e}^{(s)\top}\mathbf{e}^{(s')}) \), which is zero unless \( r=r' \) and \( s=s' \). So the \( \bv_{rs} \) are \( (a-1)(b-1) \) nonzero orthogonal vectors in a space of that dimension, an orthogonal basis of \( \C(\bP_{AB}) \). Hence \( \norm{\bP_{AB}\y}^2=\sum_{r,s}(\bv_{rs}\T\y)^2/\norm{\bv_{rs}}^2 \). Finally \( \bv_{rs}\T\y=m\sum_{i,j}d^{(rs)}_{ij}\bar{y}_{ij\cdot} \) and \( \norm{\bv_{rs}}^2=m\sum_{i,j}(d^{(rs)}_{ij})^2 \), so each term is \( \text{SS}(\mathbf{d}^{(rs)}) \) of @eq-tw-contrast-ss. Independence follows from @thm-qf-orthogonal-projections applied to the projections onto the lines \( \spn(\bv_{rs}) \) and \( \bP_E \).
:::

Product contrasts suit factors with structure. If \( \mathbf{c} \) compares a control with the treatments and \( \mathbf{e} \) is a linear trend in a quantitative factor, then \( \mathbf{c}\mathbf{e}\T \) asks whether the trend differs between them. The cell-means model of a balanced layout with quantitative levels \( w_i \) and \( z_j \) is the same as the polynomial model with all terms \( w^rz^s \), \( r<a \), \( s<b \) (@exr-tw-polynomial-surface), and the product contrasts of orthogonal polynomials give its sequential sums of squares.

## Simultaneous inference

The tests and contrasts above come in families, and the methods of [Chapter 13](../ch13-multiplicity/index.html) apply directly.

- **Tukey for the levels of one factor.** The row means \( \bar{y}_{i\cdot\cdot} \) are independent, each with variance \( \sigma^2/(bm) \). They are functions of \( (\I-\bP_E)\Y \), hence independent of \( \text{MSE} \), with or without interaction (@exr-tw-tukey-with-interaction). So @thm-mc-tukey applies with \( g=a \), multiplier \( q_\alpha(a,ab(m-1)) \) and standard error \( s/\sqrt{bm} \).
- **Scheffé for interaction contrasts.** Taking \( \mathcal V=\C(\bP_{AB}) \) in @thm-mc-scheffe gives intervals for all interaction contrasts at once, with multiplier \( \sqrt{(a-1)(b-1)F_\alpha((a-1)(b-1),ab(m-1))} \). By @cor-mc-scheffe-f, some interaction contrast is Scheffé-significant iff the \( F_{AB} \) test rejects. This is the honest way to report a contrast suggested by the interaction plot. The most extreme contrast always accounts for the whole interaction sum of squares (@exr-tw-best-contrast); the question is whether it is interpretable.
- **Bonferroni or Holm for planned families** (@thm-mc-bonferroni, @thm-mc-holm), such as the comparison of two levels of \( A \) separately at each level of \( B \).

::: {#exm-tw-bread-anova}
[Analysis of the bread trial]

For the data of @exm-tw-bread, the decomposition of @thm-tw-balanced gives

| Source | df | Sum of squares | Mean square | \( F \) | \( p \) |
|---|---|---|---|---|---|
| flour | 2 | 1627629 | 813814 | 255.74 | \( 6.6\times 10^{-17} \) |
| temperature | 3 | 140499 | 46833 | 14.72 | \( 1.2\times 10^{-5} \) |
| interaction | 6 | 98947 | 16491 | 5.18 | \( 1.5\times 10^{-3} \) |
| error | 24 | 76373 | 3182 | | |

so \( s=56.4 \) ml, and the interaction is clearly present. The script confirms that statsmodels reports the same sums of squares for Types I, II and III, as part (e) says.

*Where is the interaction?* Take the flour contrasts \( (2,-1,-1) \) (white against the other two) and \( (0,1,-1) \) (wholemeal against rye), and the orthogonal polynomial contrasts \( (-3,-1,1,3) \), \( (1,-1,-1,1) \), \( (-1,3,-3,1) \) for the equally spaced temperatures. The six products split \( \text{SS}_{AB} \) as in @prp-tw-product-contrasts:

| | linear | quadratic | cubic |
|---|---|---|---|
| white vs others | 52152 (\( F=16.39 \)) | 3134 (\( F=0.98 \)) | 6011 (\( F=1.89 \)) |
| wholemeal vs rye | 8118 (\( F=2.55 \)) | 18872 (\( F=5.93 \)) | 10660 (\( F=3.35 \)) |

The six add up to \( 98947 \). The largest piece, a fraction \( 0.527 \) of the interaction, says that white flour responds to temperature more steeply than the others. The second is the wholemeal peak at 32 °C. Does the white-by-linear component survive the most generous correction? Its estimate is \( 1444.3 \). The temperature contrast has \( \sum_je_jt_j=-72-28+32+108=40 \), and the flour contrast is twice "white minus the average of the others", so the estimate divided by \( 80 \) is the difference between white flour's slope and the average slope of the other two: \( 18.05 \) ml per °C. The Scheffé multiplier for interaction contrasts is \( 3.879 \), and the simultaneous interval for this slope difference is \( [0.75,\ 35.35] \) ml per °C. It excludes zero even though Scheffé's correction allows for every contrast the plot might have suggested. The temperature main effect is almost entirely linear, with sequential sums of squares \( 134098 \), \( 6373 \) and \( 28 \) for the linear, quadratic and cubic components.

*Comparing the flours.* The equal-weight flour means are \( 2074.8 \), \( 1727.1 \) and \( 1565.2 \) ml. Tukey's multiplier \( q_{0.05}(3,24)=3.532 \) with standard error \( s/\sqrt{12} \) gives half-width \( 57.5 \) ml. White minus wholemeal is \( 347.7 \), with interval \( [290.2,\ 405.3] \). White minus rye blend is \( 509.7 \) \( [452.2,\ 567.2] \), and wholemeal minus rye blend is \( 161.9 \) \( [104.4,\ 219.4] \). Because of the interaction, these average over the four temperatures the trial happened to use. The simple effects of white against wholemeal, with Bonferroni intervals for the four temperatures (multiplier \( t_{24,\,0.05/8}=2.700 \), half-width \( 124.4 \)), are \( 250.7 \) \( [126.3,\ 375.0] \), \( 376.7 \) \( [252.3,\ 501.0] \), \( 302.3 \) \( [178.0,\ 426.7] \) and \( 461.3 \) \( [337.0,\ 585.7] \). All four lie above zero, so with simultaneous confidence \( 95\% \) white flour is better at every temperature, and hence (@prp-tw-weighted-main(c)) under any weighting of the temperatures.
:::

```{.python .run #cell-bread-factorial-decompose}
import itertools

import numpy as np
from scipy import stats

flours = ["white", "wholemeal", "rye blend"]
temps = [24, 28, 32, 36]
true_means = np.array([[1880, 2010, 2130, 2190],
                       [1660, 1760, 1810, 1750],
                       [1500, 1540, 1580, 1590]], dtype=float)
a, b, m = 3, 4, 3
data_rng = np.random.default_rng(16_03)
volume = np.round(true_means[:, :, None] + data_rng.normal(scale=55, size=(a, b, m)))

def Jbar(k):
    return np.full((k, k), 1.0 / k)


def Cen(k):
    return np.eye(k) - Jbar(k)


y = volume.ravel()                                    # replicate fastest, then temperature, then flour
P = {"flour": np.kron(np.kron(Cen(a), Jbar(b)), Jbar(m)),
     "temperature": np.kron(np.kron(Jbar(a), Cen(b)), Jbar(m)),
     "interaction": np.kron(np.kron(Cen(a), Cen(b)), Jbar(m)),
     "error": np.kron(np.kron(np.eye(a), np.eye(b)), Cen(m))}
ss = {k: y @ Pk @ y for k, Pk in P.items()}
df = {k: int(round(np.trace(Pk))) for k, Pk in P.items()}
mse = ss["error"] / df["error"]
for k in ["flour", "temperature", "interaction", "error"]:
    F = ss[k] / df[k] / mse
    tail = f"F = {F:7.2f}  p = {stats.f.sf(F, df[k], df['error']):.2g}" if k != "error" else ""
    print(f"{k:12s} df {df[k]:2d}  SS {ss[k]:10.0f}  MS {ss[k] / df[k]:9.0f}  {tail}")
```

```{.python .run #cell-bread-factorial-contrasts}
cell = volume.mean(axis=2)
flour_c = {"white vs others": np.array([2.0, -1, -1]), "wholemeal vs rye": np.array([0.0, 1, -1])}
temp_c = {"linear": np.array([-3.0, -1, 1, 3]), "quadratic": np.array([1.0, -1, -1, 1]),
          "cubic": np.array([-1.0, 3, -3, 1])}
rows = []
for (fn, c), (tn, d) in itertools.product(flour_c.items(), temp_c.items()):
    D = np.outer(c, d)                                # an interaction contrast in the cell means
    est = (D * cell).sum()
    ss_c = est ** 2 / ((D ** 2).sum() / m)            # its one-degree-of-freedom sum of squares
    rows.append((fn, tn, est, ss_c, ss_c / mse))
    print(f"{fn:17s} x {tn:9s}  estimate {est:8.1f}  SS {ss_c:8.0f}  F {ss_c / mse:6.2f}")
print(f"sum of the six: {sum(r[3] for r in rows):.0f}   interaction SS: {ss['interaction']:.0f}")
```

```{.python .run #cell-bread-factorial-intervals}
fmean, tmean, grand = cell.mean(1), cell.mean(0), cell.mean()
s = np.sqrt(mse)
nu = df["error"]
q = stats.studentized_range.ppf(0.95, a, nu)
half_tukey = q * s / np.sqrt(b * m)                   # flour means average b*m loaves
for i, k in itertools.combinations(range(a), 2):
    d = fmean[i] - fmean[k]
    print(f"{flours[i]:9s} - {flours[k]:9s}: {d:7.1f}  Tukey [{d - half_tukey:7.1f}, {d + half_tukey:7.1f}]")
t_bonf = stats.t.ppf(1 - 0.05 / (2 * b), nu)          # white - wholemeal at each temperature
half_simple = t_bonf * s * np.sqrt(2 / m)
for j in range(b):
    d = cell[0, j] - cell[1, j]
    print(f"white - wholemeal at {temps[j]} C: {d:6.1f}  Bonferroni [{d - half_simple:6.1f}, {d + half_simple:6.1f}]")
c_scheffe = np.sqrt(df["interaction"] * stats.f.ppf(0.95, df["interaction"], nu))
```

::: {.warning}
It is tempting to drop a nonsignificant interaction and "pool" its sum of squares into the error. With \( ab(m-1) \) error degrees of freedom already available this gains little. If the interaction is real, the pooled mean square overestimates \( \sigma^2 \), and if the decision to pool depends on the data, the nominal levels no longer hold ([Section 16.4](04-higher-way.html)).
:::

## Exercises

### A. Check your understanding

::: {#exr-tw-cells-split}
[A1]

Let \( \text{SS}_{\text{cells}}=m\sum_{i,j}(\bar{y}_{ij\cdot}-\bar{y}_{\cdot\cdot\cdot})^2 \), the between-cells sum of squares of the one-way layout with \( ab \) groups. Show that \( \text{SS}_{\text{cells}}=\text{SS}_A+\text{SS}_B+\text{SS}_{AB} \), and identify each term with a projection.
:::

### B. Practice

::: {#exr-tw-single-cell}
[B1]

Suppose all cell means are equal except \( \mu_{11} \), which is larger by \( \delta \). Find \( \alpha_i \), \( \beta_j \) and \( \gamma_{ij} \), and show that the noncentralities of \( F_A \), \( F_B \) and \( F_{AB} \) are
\[
\frac{m\delta^2}{\sigma^2}\cdot\frac{a-1}{ab},\qquad
\frac{m\delta^2}{\sigma^2}\cdot\frac{b-1}{ab},\qquad
\frac{m\delta^2}{\sigma^2}\cdot\frac{(a-1)(b-1)}{ab}.
\]
Which test is most likely to detect a single deviant cell when \( a \) and \( b \) are large?
:::

::: {.solution}
The row means are \( \bar{\mu}_{1\cdot}=\mu+\delta/b \) and \( \bar{\mu}_{i\cdot}=\mu \) otherwise, and the grand mean is \( \mu+\delta/(ab) \). So \( \alpha_1=\delta(a-1)/(ab) \) and \( \alpha_i=-\delta/(ab) \) for \( i\ge2 \), and
\[
\sum_i\alpha_i^2=\frac{\delta^2}{a^2b^2}\bigl((a-1)^2+(a-1)\bigr)=\frac{\delta^2(a-1)}{ab^2},
\]
so \( bm\sum_i\alpha_i^2/\sigma^2=m\delta^2(a-1)/(ab\sigma^2) \). By symmetry the noncentrality for \( B \) is \( m\delta^2(b-1)/(ab\sigma^2) \). For the interaction, \( \gamma_{ij}=\delta(g_i-1/a)(h_j-1/b) \), where \( g_1=h_1=1 \) and \( g_i=h_j=0 \) for \( i,j\ge2 \), so \( \sum\gamma_{ij}^2=\delta^2\cdot\frac{a-1}{a}\cdot\frac{b-1}{b} \), and the noncentrality is \( m\delta^2(a-1)(b-1)/(ab\sigma^2) \). For large \( a,b \) the interaction noncentrality approaches \( m\delta^2/\sigma^2 \), while the main-effect noncentralities approach \( m\delta^2/(b\sigma^2) \) and \( m\delta^2/(a\sigma^2) \). Almost all of a single deviant cell shows up as interaction. It is, however, spread over \( (a-1)(b-1) \) degrees of freedom, and the Scheffé-type contrast aimed at that cell is the one to look at.
:::

::: {#exr-tw-tukey-with-interaction}
[B2]

In the cell-means model @eq-tw-cell-means with normal errors, show that \( \bar{y}_{1\cdot\cdot},\dots,\bar{y}_{a\cdot\cdot} \) are independent, each with variance \( \sigma^2/(bm) \), and jointly independent of \( \text{SSE} \), whatever the \( \mu_{ij} \). Deduce Tukey intervals for all \( \bar{\mu}_{i\cdot}-\bar{\mu}_{k\cdot} \) with exact simultaneous coverage.
:::

::: {.solution}
Different rows use disjoint sets of the independent observations, so the row means are independent, and each averages \( bm \) observations. The vector of row means is a linear function of \( (\I_a\otimes\bar{\mathbf{J}}_b\otimes\bar{\mathbf{J}}_m)\Y \), whose range lies in \( \C(\bP_0+\bP_A)\perp\C(\bP_E) \). By @thm-qf-orthogonal-projections, \( (\bP_0+\bP_A)\Y \) and \( \bP_E\Y \) are independent, and \( \text{SSE}=\norm{\bP_E\Y}^2 \). The row means have expectations \( \bar{\mu}_{i\cdot} \). These are exactly the conditions used in the proof of @thm-mc-tukey, with \( g=a \), common variance \( \sigma^2/(bm) \) and \( \nu=ab(m-1) \). So the intervals \( \bar{y}_{i\cdot\cdot}-\bar{y}_{k\cdot\cdot}\pm q_\alpha(a,ab(m-1))\,s/\sqrt{bm} \) have simultaneous coverage exactly \( 1-\alpha \).
:::

::: {#exr-tw-polynomial-surface}
[B3]

Let the levels of \( A \) and \( B \) be distinct numbers \( w_1,\dots,w_a \) and \( z_1,\dots,z_b \). Show that the vectors with entries \( w_i^rz_j^s \) (repeated over replicates), \( 0\le r<a \), \( 0\le s<b \), span the cell-means space \( \C(\X) \). Show that the additive model corresponds to the terms with \( r=0 \) or \( s=0 \). What does the coefficient of \( wz \) measure?
:::

### C. Going deeper

::: {#exr-tw-best-contrast}
[C1]

Let \( \hat{\gamma}_{ij}=\bar{y}_{ij\cdot}-\bar{y}_{i\cdot\cdot}-\bar{y}_{\cdot j\cdot}+\bar{y}_{\cdot\cdot\cdot} \). Show that among all interaction contrasts, \( \text{SS}(\mathbf{d}) \) of @eq-tw-contrast-ss is largest for \( d_{ij}\propto\hat{\gamma}_{ij} \), and that the maximum is \( \text{SS}_{AB} \). Relate this to @cor-mc-scheffe-f(a). Why does it mean that the \( F \) statistic of a contrast chosen after inspecting the interaction plot cannot be compared with \( F(1,ab(m-1)) \)?
:::

::: {.solution}
Since \( (d_{ij}) \) has zero margins, \( \sum d_{ij}\bar{y}_{ij\cdot}=\sum d_{ij}\hat{\gamma}_{ij} \) (@prp-tw-interaction-equiv applied to the table of cell means). By the Cauchy–Schwarz inequality, \( \bigl(\sum d_{ij}\hat{\gamma}_{ij}\bigr)^2\le\sum d_{ij}^2\sum\hat{\gamma}_{ij}^2 \), with equality iff \( \mathbf{d}\propto\hat{\bgamma} \). Hence \( \text{SS}(\mathbf{d})\le m\sum\hat{\gamma}_{ij}^2=\text{SS}_{AB} \), with equality at \( \mathbf{d}=\hat{\bgamma} \), which is itself a zero-margin table. In the notation of @cor-mc-scheffe-f(a) this says that the largest single-contrast \( F \) ratio is \( (a-1)(b-1)F_{AB} \). A contrast chosen for its large value is, in the extreme, this maximizing contrast. Its null distribution is that of \( (a-1)(b-1) \) times an \( F((a-1)(b-1),ab(m-1)) \) variable, not \( F(1,ab(m-1)) \).
:::

::: {#exr-tw-simple-effects-ss}
[C2]

For each column \( j \), let \( \text{SS}_{A\mid j}=m\sum_i(\bar{y}_{ij\cdot}-\bar{y}_{\cdot j\cdot})^2 \) be the sum of squares for \( A \) computed within column \( j \) alone. Show that \( \sum_j\text{SS}_{A\mid j}=\text{SS}_A+\text{SS}_{AB} \), and find the projection whose squared length is the left side. What hypothesis does the corresponding \( F \) test with \( b(a-1) \) numerator degrees of freedom address?
:::
