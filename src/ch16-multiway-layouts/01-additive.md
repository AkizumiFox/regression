# Two-way additive models

Two factors are observed together. Factor \( A \) has levels \( i=1,\dots,a \) and factor \( B \) has levels \( j=1,\dots,b \), with \( a,b\ge2 \). Every one of the \( ab \) combinations, or **cells**, is observed the same number \( m\ge1 \) of times. We write \( y_{ijk} \) for the \( k \)th observation in cell \( (i,j) \), so that \( n=abm \). @exm-proj-two-factor found the projections of the additive model for \( m=1 \) by a direct argument. This section redoes it with the tools used throughout the chapter, and adds replication, estimable functions, tests, and a test of additivity with one observation per cell.

## Averaging and centring matrices

Stack the observations with the replicate index \( k \) running fastest, then \( j \), then \( i \):
\[
\y=(y_{111},\dots,y_{11m},\ y_{121},\dots,y_{12m},\ \dots,\ y_{ab1},\dots,y_{abm})\T .
\]
With this order a matrix of the form \( \mathbf{E}_a\otimes\mathbf{E}_b\otimes\mathbf{E}_m \) acts on each index separately, as in @exm-mat-two-way-kronecker. For \( k\ge1 \) write
\[
\bar{\mathbf{J}}_k=\frac1k\bone_k\bone_k\T,\qquad \mathbf{C}_k=\I_k-\bar{\mathbf{J}}_k .
\]
The **averaging matrix** \( \bar{\mathbf{J}}_k \) replaces each entry of a \( k \)-vector by the mean of the entries. The **centring matrix** \( \mathbf{C}_k \) subtracts that mean. With the usual dot notation for means, the following products replace \( y_{ijk} \) by
\[
\begin{aligned}
(\I_a\otimes\I_b\otimes\bar{\mathbf{J}}_m)\y&:\ \bar{y}_{ij\cdot}, &
(\I_a\otimes\bar{\mathbf{J}}_b\otimes\bar{\mathbf{J}}_m)\y&:\ \bar{y}_{i\cdot\cdot},\\
(\bar{\mathbf{J}}_a\otimes\I_b\otimes\bar{\mathbf{J}}_m)\y&:\ \bar{y}_{\cdot j\cdot}, &
(\bar{\mathbf{J}}_a\otimes\bar{\mathbf{J}}_b\otimes\bar{\mathbf{J}}_m)\y&:\ \bar{y}_{\cdot\cdot\cdot},
\end{aligned}
\]
the cell, row, column and grand means. Each operator averages over the indices where it has a \( \bar{\mathbf{J}} \) and leaves the others alone. Every projection in this chapter is a Kronecker product of factors \( \I \), \( \bar{\mathbf{J}} \) and \( \mathbf{C} \), and the following lemma is all the algebra they need.

::: {#lem-tw-kron}
[Kronecker products of averaging and centring matrices]

Let \( l_1,\dots,l_t\ge1 \).

::: {.enumerate options="label=(\alph*)"}
1. \( \bar{\mathbf{J}}_k \) and \( \mathbf{C}_k \) are symmetric and idempotent, with ranks \( 1 \) and \( k-1 \), and \( \bar{\mathbf{J}}_k\mathbf{C}_k=\bzero \), \( \bar{\mathbf{J}}_k+\mathbf{C}_k=\I_k \).

2. If each \( \mathbf{E}_s \) is one of \( \I_{l_s} \), \( \bar{\mathbf{J}}_{l_s} \), \( \mathbf{C}_{l_s} \), then \( \mathbf{E}_1\otimes\dots\otimes\mathbf{E}_t \) is an orthogonal projection whose rank is the product of the ranks of the \( \mathbf{E}_s \).

3. If each \( \mathbf{E}_s \) and each \( \mathbf{F}_s \) is one of \( \bar{\mathbf{J}}_{l_s} \), \( \mathbf{C}_{l_s} \), and \( \mathbf{E}_s\ne\mathbf{F}_s \) for at least one \( s \), then \( (\mathbf{E}_1\otimes\dots\otimes\mathbf{E}_t)(\mathbf{F}_1\otimes\dots\otimes\mathbf{F}_t)=\bzero \).

4. The \( 2^t \) products \( \mathbf{E}_1\otimes\dots\otimes\mathbf{E}_t \) with each \( \mathbf{E}_s\in\{\bar{\mathbf{J}}_{l_s},\mathbf{C}_{l_s}\} \) sum to \( \I \).
:::

:::

::: {.proof}
(a) \( \bar{\mathbf{J}}_k \) is symmetric and \( \bar{\mathbf{J}}_k^2=k^{-2}\bone(\bone\T\bone)\bone\T=\bar{\mathbf{J}}_k \). Its trace is \( 1 \). Then \( \mathbf{C}_k=\I-\bar{\mathbf{J}}_k \) is the complementary projection, with trace \( k-1 \), and \( \bar{\mathbf{J}}_k\mathbf{C}_k=\bar{\mathbf{J}}_k-\bar{\mathbf{J}}_k^2=\bzero \). Ranks equal traces by @prp-proj-trace-rank. (b) By @prp-mat-kronecker(e) a Kronecker product of symmetric idempotent matrices is symmetric and idempotent, and by @prp-mat-kronecker(c) its rank is the product of the ranks. (c) By the mixed-product rule the product equals \( \mathbf{E}_1\mathbf{F}_1\otimes\dots\otimes\mathbf{E}_t\mathbf{F}_t \), and the factor with \( \mathbf{E}_s\ne\mathbf{F}_s \) is \( \bar{\mathbf{J}}\mathbf{C}=\bzero \) or \( \mathbf{C}\bar{\mathbf{J}}=\bzero \). (d) Expand \( \I=(\bar{\mathbf{J}}_{l_1}+\mathbf{C}_{l_1})\otimes\dots\otimes(\bar{\mathbf{J}}_{l_t}+\mathbf{C}_{l_t}) \) by bilinearity of the Kronecker product.
:::

The rank \( k-1 \) of \( \mathbf{C}_k \) is the source of every "levels minus one" in the analysis of variance, and part (d) is the source of every analysis of variance table in this chapter. Part (b) allows \( k=1 \): then \( \bar{\mathbf{J}}_1=1 \) and \( \mathbf{C}_1=0 \).

```{.python .run #cell-kron-projections-kron}
import itertools

import numpy as np

def Jbar(k):
    return np.full((k, k), 1.0 / k)          # averaging: replaces a k-vector by its mean


def Cen(k):
    return np.eye(k) - Jbar(k)               # centring: subtracts the mean


def kron_all(mats):
    out = np.ones((1, 1))
    for A in mats:
        out = np.kron(out, A)
    return out


def factorial_projections(levels, m):
    """P_S for every subset S of the factors (a 0/1 tuple), and the error projection."""
    P = {}
    for S in itertools.product([0, 1], repeat=len(levels)):
        P[S] = kron_all([Cen(l) if s else Jbar(l) for l, s in zip(levels, S)] + [Jbar(m)])
    P["error"] = kron_all([np.eye(l) for l in levels] + [Cen(m)])
    return P


def rank(A):
    return int(round(np.trace(A)))           # rank = trace for a projection
```

## The additive model

The **additive two-way model** is
\[
\E(y_{ijk})=\mu+\alpha_i+\beta_j,\qquad \Cov(\Y)=\sigma^2\I ,
\]{#eq-tw-additive}

with model matrix \( \X=[\bone_n,\ \Z_A,\ \Z_B] \), where \( \Z_A=\I_a\otimes\bone_b\otimes\bone_m \) and \( \Z_B=\bone_a\otimes\I_b\otimes\bone_m \) are the indicator matrices of the two factors. It says that the difference between two levels of \( A \) is the same in every column and, symmetrically, that the difference between two levels of \( B \) is the same in every row. The model has \( 1+a+b \) parameters and, by @thm-est-connected(a), rank \( a+b-1 \). Define
\[
\begin{aligned}
\bP_0&=\bar{\mathbf{J}}_a\otimes\bar{\mathbf{J}}_b\otimes\bar{\mathbf{J}}_m,\\
\bP_A&=\mathbf{C}_a\otimes\bar{\mathbf{J}}_b\otimes\bar{\mathbf{J}}_m,\\
\bP_B&=\bar{\mathbf{J}}_a\otimes\mathbf{C}_b\otimes\bar{\mathbf{J}}_m ,
\end{aligned}
\]{#eq-tw-additive-projections}

and \( \R=\I-\bP_0-\bP_A-\bP_B \). By the table of means above, \( \bP_A\y \) has entry \( \bar{y}_{i\cdot\cdot}-\bar{y}_{\cdot\cdot\cdot} \) in every position of row \( i \), and \( \bP_B\y \) has entry \( \bar{y}_{\cdot j\cdot}-\bar{y}_{\cdot\cdot\cdot} \) in every position of column \( j \).

## Additivity with one observation per cell

Many two-factor data sets have one observation per cell, such as methods each run once on several test problems. By @exm-proj-two-factor the residual space of the additive model then has dimension \( (a-1)(b-1) \), and the residuals are the interaction terms \( y_{ij}-\bar{y}_{i\cdot}-\bar{y}_{\cdot j}+\bar{y}_{\cdot\cdot} \). The residual sum of squares then serves two purposes. It estimates \( \sigma^2 \) if the model is additive, and it is the only place where nonadditivity could show. There is no model-free estimate of \( \sigma^2 \), and no test of additivity against a general alternative.

A test is still possible against a *specific* kind of departure. Suppose the cell means are additive on some other scale, \( \mu_{ij}=f(\nu+a_i+b_j) \), for a smooth increasing function \( f \) and small effects \( a_i \), \( b_j \) with \( \sum a_i=\sum b_j=0 \). A second-order Taylor expansion about \( \nu \) gives
\[
\mu_{ij}\approx f(\nu)+f'(\nu)(a_i+b_j)+\tfrac12f''(\nu)\bigl(a_i^2+b_j^2\bigr)+f''(\nu)\,a_ib_j .
\]
The terms in \( a_i^2 \) and \( b_j^2 \) are additive. Only the product term is not. To first order the row and column effects of \( \mu_{ij} \) are \( \alpha_i\approx f'(\nu)a_i \) and \( \beta_j\approx f'(\nu)b_j \), so the interaction is approximately
\[
\mu_{ij}-\bar{\mu}_{i\cdot}-\bar{\mu}_{\cdot j}+\bar{\mu}_{\cdot\cdot}\approx\theta\,\alpha_i\beta_j,\qquad \theta=\frac{f''(\nu)}{f'(\nu)^2}.
\]{#eq-tw-tukey-alternative}

Tukey (1949) proposed to look for interaction of exactly this shape. The row and column effects are estimated by \( r_i=\bar{y}_{i\cdot}-\bar{y}_{\cdot\cdot} \) and \( c_j=\bar{y}_{\cdot j}-\bar{y}_{\cdot\cdot} \) (with \( m=1 \) we drop the third subscript). The residuals are then regressed on the products \( r_ic_j \). The table \( (r_ic_j) \) has zero row and column sums, so it lies in the residual space. The regression through the origin of the residual table on it has slope and sum of squares
\[
\hat{\theta}=\frac{\sum_{i,j}r_ic_jy_{ij}}{\sum_ir_i^2\sum_jc_j^2},\qquad
\text{SS}_N=\frac{\bigl(\sum_{i,j}r_ic_jy_{ij}\bigr)^2}{\sum_ir_i^2\sum_jc_j^2}.
\]{#eq-tw-tukey-ss}

The regressor is computed from the same data, so @prp-cor-augmented-test does not apply. But a regressor that is a function of the fitted values alone leaves the \( F \) test exact (@thm-cor-fitted-regressors), and part (f) of the next theorem shows that Tukey's regressor is of this kind.

## Estimation and tests

The following theorem collects what the additive model delivers in a balanced layout, with and without replication.

::: {#thm-tw-additive}
[The additive two-way model]

Assume @eq-tw-additive with \( m\ge1 \) observations in every cell.

::: {.enumerate options="label=(\alph*)"}
1. \( \bP_0,\bP_A,\bP_B \) and \( \R \) are mutually orthogonal projections summing to \( \I \), with ranks \( 1 \), \( a-1 \), \( b-1 \) and \( n-a-b+1 \). The projection onto \( \C(\X) \) is \( \M=\bP_0+\bP_A+\bP_B \), and the fitted values are \( \hat{y}_{ijk}=\bar{y}_{i\cdot\cdot}+\bar{y}_{\cdot j\cdot}-\bar{y}_{\cdot\cdot\cdot} \).

2. \( \lambda_0\mu+\sum_i\lambda_i\alpha_i+\sum_j\kappa_j\beta_j \) is estimable iff \( \lambda_0=\sum_i\lambda_i=\sum_j\kappa_j \). The BLUE of a contrast \( \sum_ic_i\alpha_i \) is \( \sum_ic_i\bar{y}_{i\cdot\cdot} \), with variance \( \sigma^2\sum_ic_i^2/(bm) \). These are the estimate and variance of the one-way layout in \( A \) that ignores \( B \). The estimate of a contrast in the \( \alpha_i \) is uncorrelated with the estimate of any contrast in the \( \beta_j \).

3. Under normal errors, \( \text{SS}_A=\norm{\bP_A\y}^2=bm\sum_i(\bar{y}_{i\cdot\cdot}-\bar{y}_{\cdot\cdot\cdot})^2 \), \( \text{SS}_B=\norm{\bP_B\y}^2=am\sum_j(\bar{y}_{\cdot j\cdot}-\bar{y}_{\cdot\cdot\cdot})^2 \) and \( \text{SSE}=\norm{\R\y}^2 \) are independent, \( \text{SSE}/\sigma^2\sim\chi^2(n-a-b+1) \), and
   \[
F_A=\frac{\text{SS}_A/(a-1)}{\text{SSE}/(n-a-b+1)}\sim F\bigl(a-1,\ n-a-b+1,\ \gamma_A\bigr),\qquad
\gamma_A=\frac{bm}{\sigma^2}\sum_i(\alpha_i-\bar{\alpha})^2 ,
\]
   with \( \bar{\alpha}=a^{-1}\sum_i\alpha_i \). So \( F_A \) is an exact test of \( \alpha_1=\dots=\alpha_a \). Its numerator is that of the one-way \( F \) test in \( A \) that ignores \( B \) (@thm-aov-oneway-f); only the error term differs. The same holds for \( B \), with \( \gamma_B=(am/\sigma^2)\sum_j(\beta_j-\bar{\beta})^2 \).

4. If \( m\ge2 \), then \( \R=\bP_{AB}+\bP_E \) with \( \bP_{AB}=\mathbf{C}_a\otimes\mathbf{C}_b\otimes\bar{\mathbf{J}}_m \) and \( \bP_E=\I_a\otimes\I_b\otimes\mathbf{C}_m \), of ranks \( (a-1)(b-1) \) and \( ab(m-1) \). The squared lengths \( \norm{\bP_{AB}\y}^2 \) and \( \norm{\bP_E\y}^2 \) are the lack-of-fit and pure-error sums of squares of the additive model (@def-cor-pure-error). If the cell means \( \mu_{ij} \) are arbitrary, the lack-of-fit \( F \) test of additivity has noncentrality
   \[
\gamma_{AB}=\frac m{\sigma^2}\sum_{i,j}\bigl(\mu_{ij}-\bar{\mu}_{i\cdot}-\bar{\mu}_{\cdot j}+\bar{\mu}_{\cdot\cdot}\bigr)^2 .
\]

5. If \( m=1 \), then \( \R=\mathbf{C}_a\otimes\mathbf{C}_b \), and the residual space consists of the \( a\times b \) tables whose rows and columns all sum to zero. There is no pure error.

6. *(Tukey's test.)* If \( m=1 \), the errors are normal and \( (a-1)(b-1)\ge2 \), then with \( \text{SS}_N \) as in @eq-tw-tukey-ss,
   \[
F_N=\frac{\text{SS}_N}{(\text{SSE}-\text{SS}_N)/\bigl((a-1)(b-1)-1\bigr)}\sim F\bigl(1,\ (a-1)(b-1)-1\bigr).
\]
:::

:::

::: {.proof}
(a) By @lem-tw-kron(c) the three projections \( \bP_0,\bP_A,\bP_B \) are mutually orthogonal, and their ranks are \( 1 \), \( a-1 \), \( b-1 \) by @lem-tw-kron(b). Moreover \( \bP_0+\bP_A=\I_a\otimes\bar{\mathbf{J}}_b\otimes\bar{\mathbf{J}}_m \) replaces each observation by its row mean. It is the projection onto \( \C(\Z_A) \), the vectors that are constant on rows (@exm-proj-oneway-M). Likewise \( \bP_0+\bP_B \) projects onto \( \C(\Z_B) \). Hence
\[
\C(\X)=\C(\Z_A)+\C(\Z_B)=\bigl(\C(\bP_0)\dirsum\C(\bP_A)\bigr)+\bigl(\C(\bP_0)\dirsum\C(\bP_B)\bigr)=\C(\bP_0)\dirsum\C(\bP_A)\dirsum\C(\bP_B),
\]
a sum of mutually orthogonal subspaces. By @thm-proj-sum its projection is \( \bP_0+\bP_A+\bP_B \). Then \( \R=\I-\M \) is the complementary projection, orthogonal to the other three, with rank \( n-(a+b-1) \). The fitted values are \( \M\y=(\bP_0+\bP_A)\y+(\bP_0+\bP_B)\y-\bP_0\y \).

(b) The design is complete, hence connected, and the estimability condition is @exr-est-disconnected with a single component. For the estimator, \( \E(\bar{y}_{i\cdot\cdot})=\mu+\alpha_i+\bar{\beta} \), so \( \sum_ic_i\bar{y}_{i\cdot\cdot} \) is unbiased for \( \sum_ic_i\alpha_i \) when \( \sum_ic_i=0 \). It equals \( \boldsymbol{\uprho}\T\y \) with \( \boldsymbol{\uprho}=(bm)^{-1}\mathbf{c}\otimes\bone_b\otimes\bone_m \), and \( \boldsymbol{\uprho}=\bP_A\boldsymbol{\uprho}\in\C(\X) \) because \( \mathbf{c}\perp\bone_a \). So \( (\I-\M)\boldsymbol{\uprho}=\bzero \), and by @eq-opt-gm-split its variance equals that of the BLUE. By @thm-opt-gauss-markov(c) it is the BLUE. Its variance is \( \sigma^2\norm{\boldsymbol{\uprho}}^2=\sigma^2\sum_ic_i^2/(bm) \). The one-way layout in \( A \) has \( bm \) observations per level, and @thm-est-oneway-contrasts gives the same estimator and variance. A contrast in the \( \beta_j \) has coefficient vector \( \boldsymbol{\uprho}_B\in\C(\bP_B) \), and \( \Cov(\boldsymbol{\uprho}\T\Y,\boldsymbol{\uprho}_B\T\Y)=\sigma^2\boldsymbol{\uprho}\T\boldsymbol{\uprho}_B=0 \).

(c) The spaces \( \C(\bP_0),\C(\bP_A),\C(\bP_B) \) form an analysis of variance decomposition of \( \C(\X) \) (@def-ss-anova-decomposition), so @cor-ss-ms-distributions applies. The mean vector has entries \( \mu+\alpha_i+\beta_j \). Its row means minus its grand mean are \( \alpha_i-\bar{\alpha} \), so \( \norm{\bP_A\bmu}^2=bm\sum_i(\alpha_i-\bar{\alpha})^2 \). This is zero iff all \( \alpha_i \) are equal. The formula for \( \text{SS}_A \) is the same computation applied to \( \y \).

(d) By @lem-tw-kron(d) with \( t=3 \), \( \I \) is the sum of the eight products \( \mathbf{E}_a\otimes\mathbf{E}_b\otimes\mathbf{E}_m \). The four with \( \mathbf{E}_m=\bar{\mathbf{J}}_m \) are \( \bP_0,\bP_A,\bP_B,\bP_{AB} \). The four with \( \mathbf{E}_m=\mathbf{C}_m \) add up to \( \I_a\otimes\I_b\otimes\mathbf{C}_m=\bP_E \). So \( \R=\bP_{AB}+\bP_E \), with ranks from @lem-tw-kron(b). The matrix \( \I-\bP_E=\I_a\otimes\I_b\otimes\bar{\mathbf{J}}_m \) replaces each observation by its cell mean. It is the projection \( \M_Z \) onto the cell-indicator space, the \( \Z \) of @def-cor-pure-error. Hence \( \bP_E=\I-\M_Z \) and \( \bP_{AB}=\M_Z-\M \). By @thm-cor-lack-of-fit the noncentrality is \( \norm{(\I-\M)\boldsymbol{\uptheta}}^2/\sigma^2 \), where \( \boldsymbol{\uptheta} \) has entries \( \mu_{ij} \). Here \( (\I-\M)\boldsymbol{\uptheta}=\bP_{AB}\boldsymbol{\uptheta} \), because \( \bP_E\boldsymbol{\uptheta}=\bzero \). Its entries are the displayed interaction terms, each repeated \( m \) times.

(e) With \( m=1 \), \( \bar{\mathbf{J}}_1=1 \) and \( \mathbf{C}_1=0 \), so \( \bP_E=\bzero \) and \( \R=\bP_{AB}=\mathbf{C}_a\otimes\mathbf{C}_b \). A vector \( \bv \), read as an \( a\times b \) table \( \mathbf{V} \), satisfies \( (\mathbf{C}_a\otimes\mathbf{C}_b)\bv=\bv \) iff \( \mathbf{C}_a\mathbf{V}\mathbf{C}_b=\mathbf{V} \), which holds iff \( \mathbf{V}\bone_b=\bzero \) and \( \bone_a\T\mathbf{V}=\bzero\T \).

(f) Write \( \hat{y}_{ij}=\bar{y}_{\cdot\cdot}+r_i+c_j \), let \( \mathbf{r}=(r_i) \) and \( \mathbf{c}=(c_j) \), and put \( \mathbf{U}=\mathbf{r}\otimes\mathbf{c} \), the table \( (r_ic_j) \). Expanding,
\[
\hat{y}_{ij}^2=\bigl(\bar{y}_{\cdot\cdot}^2+r_i^2+2\bar{y}_{\cdot\cdot}r_i\bigr)+\bigl(c_j^2+2\bar{y}_{\cdot\cdot}c_j\bigr)+2r_ic_j .
\]
The first bracket depends on \( i \) alone and the second on \( j \) alone, so both lie in \( \C(\X) \). Since \( \bone_a\T\mathbf{r}=\bone_b\T\mathbf{c}=0 \), \( \R\mathbf{U}=\mathbf{C}_a\mathbf{r}\otimes\mathbf{C}_b\mathbf{c}=\mathbf{U} \) by (e). Hence the column \( \W=(\hat{y}_{ij}^2) \), a function of \( \M\Y \), satisfies \( (\I-\M)\W=2\mathbf{U} \). By @thm-ss-adding the extra sum of squares from adding \( \W \) to the additive model is \( (\mathbf{U}\T\y)^2/\mathbf{U}\T\mathbf{U} \), and \( \mathbf{U}\T\mathbf{U}=\sum_ir_i^2\sum_jc_j^2 \), so it is \( \text{SS}_N \), and \( \text{SSE}-\text{SS}_N \) is the residual sum of squares of the augmented model. The rank condition of @thm-cor-fitted-regressors holds because \( \mathbf{U}\ne\bzero \) with probability one: \( \mathbf{U}=\bzero \) only if \( \bP_A\Y=\bzero \) or \( \bP_B\Y=\bzero \), and each event has probability zero, because \( \norm{\bP_A\Y}^2/\sigma^2 \) has a continuous (noncentral \( \chi^2 \)) distribution on \( a-1\ge1 \) degrees of freedom, and likewise for \( B \). So \( \rank[\X,\W]=a+b \), which is less than \( n=ab \) because \( (a-1)(b-1)\ge2 \). @thm-cor-fitted-regressors with \( q=1 \) and \( n-r-q=(a-1)(b-1)-1 \) gives the distribution of \( F_N \).
:::

Part (b) is the practical content of balance: rows are compared exactly as if the columns did not exist. The columns still matter, because they remove their own variation from the error.

The listing builds the five projections of part (d) for \( a=3 \), \( b=4 \), \( m=2 \) and checks their ranks, orthogonality and sum. The script also checks each identification in the proof against projections computed from indicator columns.

```{.python .run #cell-kron-projections-twoway}
a, b, m = 3, 4, 2
P = factorial_projections([a, b], m)
names = {(0, 0): "mean", (1, 0): "A", (0, 1): "B", (1, 1): "AB", "error": "error"}
for key, name in names.items():
    print(f"{name:6s} rank {rank(P[key])}")
keys = list(names)
worst = max(np.abs(P[s] @ P[t]).max() for s, t in itertools.combinations(keys, 2))
print("largest entry of any product P_s P_t, s != t:", f"{worst:.1e}")
print("sum is the identity:", np.allclose(sum(P[k] for k in keys), np.eye(a * b * m)))
```

## Tukey's test in practice

Tukey's test (@thm-tw-additive(f)) is exact against every additive model, whatever the values of \( \mu \), \( \alpha_i \), \( \beta_j \). It has power against interactions roughly proportional to \( \alpha_i\beta_j \), the kind produced by measuring on the wrong scale, and little against anything else. The estimated slope says which scale. If \( f \) is the inverse of a power transformation, so that \( y^p \) is additive, then \( f(t)=t^{1/p} \), and @eq-tw-tukey-alternative gives \( \theta=(1-p)/f(\nu) \) (@exr-tw-tukey-power). Replacing \( f(\nu) \) by the grand mean suggests the power
\[
\hat{p}=1-\hat{\theta}\,\bar{y}_{\cdot\cdot},
\]{#eq-tw-tukey-power}

with \( \hat{p}=0 \) read as the logarithm. The rule is rough, but it turns a vague complaint into a specific suggestion.

::: {#exm-tw-solvers}
[Solvers on benchmark problems]

Six solvers were each run once on eight benchmark problems, and the running time was recorded in seconds. The data are synthetic, generated by the script `tukey_runtime.py` so that the logarithm of time is additive in solver and problem, with a little noise:

| | P1 | P2 | P3 | P4 | P5 | P6 | P7 | P8 |
|---|---|---|---|---|---|---|---|---|
| S1 | 1.7 | 3.4 | 5.6 | 14.1 | 28.1 | 50.2 | 82.9 | 249.4 |
| S2 | 2.1 | 5.7 | 9.0 | 14.7 | 36.9 | 76.0 | 116.7 | 292.0 |
| S3 | 1.3 | 3.1 | 5.1 | 10.5 | 18.8 | 37.6 | 82.2 | 148.3 |
| S4 | 3.3 | 9.1 | 12.7 | 19.1 | 52.8 | 110.0 | 183.1 | 371.8 |
| S5 | 2.0 | 4.7 | 7.2 | 13.5 | 31.4 | 65.2 | 111.1 | 289.3 |
| S6 | 1.5 | 3.4 | 5.4 | 10.4 | 26.0 | 43.2 | 106.3 | 172.0 |

The additive model in seconds leaves a residual sum of squares of \( 29396 \) on \( 35 \) degrees of freedom. Tukey's single direction takes \( \text{SS}_N=26126 \) of it, a fraction \( 0.889 \), and \( F_N=271.7 \) on \( 1 \) and \( 34 \) degrees of freedom (\( p=9\times 10^{-18} \)). The interaction has exactly the shape of @eq-tw-tukey-alternative: fast problems separate the solvers by fractions of a second, slow problems by minutes. The slope is \( \hat{\theta}=0.01563 \) and the grand mean is \( 61.46 \), so @eq-tw-tukey-power suggests the power \( 0.04 \), essentially the logarithm. [Figure 16.1.1](#fig-tw-tukey) shows the diagnostic plot. On the log scale, \( F_N=0.27 \) (\( p=0.60 \)), and the residuals show no trend against \( r_ic_j/\bar{y}_{\cdot\cdot} \).

The scale changes the conclusion about the solvers. In seconds, \( F_A=3.88 \) (\( p=0.007 \)), because the nonadditivity inflates the error term. In log seconds, \( F_A=77.7 \). The residual standard deviation is \( s=0.105 \), so differences between solvers are measured in proportions of about ten per cent. By @exr-mc-tukey-two-way the solver means of log time are independent with common variance \( \sigma^2/8 \) and independent of \( s^2 \). Tukey's method therefore gives simultaneous intervals for all \( 15 \) ratios of solver speeds, with multiplier \( q_{0.05}(6,35)=4.261 \). Solver 4, the slowest, takes \( 2.50 \) times as long as solver 3, the fastest, with simultaneous \( 95\% \) interval \( [2.14,\ 2.94] \).
:::

::: {when-format="html"}
![**Figure 16.1.1.** Tukey's diagnostic plot for the solver data: residuals of the additive fit against the comparison values \( r_ic_j/\bar{y}_{\cdot\cdot} \), with the fitted line of slope \( \hat{\theta}\bar{y}_{\cdot\cdot}=1-\hat{p} \). (a) In seconds the residuals lie along a line of slope near one, which points to the logarithm. (b) In log seconds there is no trend.](tukey_diagnostic.svg){#fig-tw-tukey width=100%}
:::

::: {when-format="pdf"}
![Tukey's diagnostic plot for the solver data: residuals of the additive fit against the comparison values \( r_ic_j/\bar{y}_{\cdot\cdot} \), with the fitted line of slope \( \hat{\theta}\bar{y}_{\cdot\cdot}=1-\hat{p} \). (a) In seconds the residuals lie along a line of slope near one, which points to the logarithm. (b) In log seconds there is no trend.](tukey_diagnostic.pdf){width=100%}
:::

```{.python .run #cell-tukey-runtime-tukey}
import numpy as np
from scipy import stats

gen_rng = np.random.default_rng(16_01)
a, b = 6, 8
solver_log = np.array([0.0, 0.25, -0.35, 0.55, 0.1, -0.2])            # log speed factors
problem_log = np.log([2.0, 4.5, 7.0, 12.0, 25.0, 60.0, 110.0, 240.0])  # typical seconds
noise = gen_rng.normal(scale=0.12, size=(a, b))
times = np.round(np.exp(solver_log[:, None] + problem_log[None, :] + noise), 1)

def additive_fit(Y):
    """Row effects, column effects and residuals of the additive fit to an a x b table."""
    g = Y.mean()
    r = Y.mean(axis=1) - g                    # ybar_i. - ybar..
    c = Y.mean(axis=0) - g                    # ybar_.j - ybar..
    resid = Y - g - r[:, None] - c[None, :]
    return g, r, c, resid


def tukey_test(Y):
    """Tukey's one-degree-of-freedom test for nonadditivity in an a x b table."""
    a, b = Y.shape
    g, r, c, resid = additive_fit(Y)
    u = np.outer(r, c)                        # lies in the interaction space
    ss_n = (u * Y).sum() ** 2 / (u * u).sum()
    sse = (resid ** 2).sum()
    df = (a - 1) * (b - 1) - 1
    F = ss_n / ((sse - ss_n) / df)
    theta = (u * Y).sum() / (u * u).sum()     # coefficient of a_i b_j
    return F, stats.f.sf(F, 1, df), theta, g, sse, ss_n


for label, Y in [("seconds", times), ("log seconds", np.log(times))]:
    F, p, theta, g, sse, ss_n = tukey_test(Y)
    print(f"{label:12s} F = {F:9.2f}  p = {p:.2g}  theta = {theta:.4f}"
          f"  suggested power 1 - theta*ybar = {1 - theta * g:.2f}")
```

The listing below checks @thm-tw-additive(f) by simulation. It draws \( 20000 \) tables from an additive normal model with \( a=6 \) and \( b=8 \). The \( 5\% \) test rejects in a fraction \( 0.0508 \) of them, and a Kolmogorov–Smirnov comparison of the simulated statistics with \( F(1,34) \) gives \( p=0.98 \).

```{.python .run #cell-tukey-runtime-null}
rng = np.random.default_rng(20260919)
reps = 20_000
Fs = np.empty(reps)
for t in range(reps):
    Y = 3.0 + np.arange(a)[:, None] * 0.5 + np.arange(b)[None, :] + rng.normal(size=(a, b))
    Fs[t] = tukey_test(Y)[0]
size = np.mean(Fs > stats.f.ppf(0.95, 1, (a - 1) * (b - 1) - 1))
ks = stats.kstest(Fs, stats.f(1, (a - 1) * (b - 1) - 1).cdf)
print(f"simulated size of the 5% test: {size:.4f}   Kolmogorov-Smirnov p-value: {ks.pvalue:.2f}")
```

::: {.warning}
A nonsignificant Tukey test does not establish additivity. It rules out only interaction of the form \( \theta\alpha_i\beta_j \). A single aberrant cell or row can leave \( F_N \) small and still distort every conclusion, so look at the residual table and an interaction plot ([Section 16.2](02-interaction.html)) as well.
:::

## Exercises

### A. Check your understanding

::: {#exr-tw-two-by-two}
[A1]

For \( a=b=2 \) and \( m=1 \), write out the \( 4\times4 \) matrices \( \bP_0 \), \( \bP_A \), \( \bP_B \) and \( \R \) of @eq-tw-additive-projections. Show that \( \R\y=\tfrac14(y_{11}-y_{12}-y_{21}+y_{22})(1,-1,-1,1)\T \), and explain why Tukey's test is not available here.
:::

::: {#exr-tw-row-contrast-ss}
[A2]

In the additive model with \( m \) observations per cell, show that the sum of squares for testing a single contrast \( \sum_ic_i\alpha_i=0 \) is \( \bigl(\sum_ic_i\bar{y}_{i\cdot\cdot}\bigr)^2/\bigl(\sum_ic_i^2/(bm)\bigr) \), the same as in the one-way layout in \( A \) (@prp-aov-contrast-ss). What differs between the two \( F \) tests of this contrast, the one here and the one-way test that ignores \( B \)?
:::

### B. Practice

::: {#exr-tw-additive-power}
[B1]

In the additive model, suppose \( \alpha_1=\delta \) and \( \alpha_2=\dots=\alpha_a=0 \). Find \( \gamma_A \). With \( a=4 \), \( b=5 \) and \( \delta=\sigma \), compare the power of the \( 5\% \) test for \( m=1 \) and \( m=2 \) using @thm-glh-power, and explain why doubling \( m \) more than doubles the error degrees of freedom.
:::

::: {.solution}
Here \( \bar{\alpha}=\delta/a \), and \( \sum_i(\alpha_i-\bar{\alpha})^2=\delta^2(1-1/a)^2+(a-1)\delta^2/a^2=\delta^2(a-1)/a \). So \( \gamma_A=bm\delta^2(a-1)/(a\sigma^2) \). With \( a=4 \), \( b=5 \), \( \delta=\sigma \) this is \( 3.75m \). The error degrees of freedom are \( n-a-b+1=20m-8 \), which is \( 12 \) for \( m=1 \) and \( 32 \) for \( m=2 \). The power is \( \Pr\{F(3,20m-8,3.75m)>F_{0.05}(3,20m-8)\} \), about \( 0.25 \) for \( m=1 \) and \( 0.56 \) for \( m=2 \) (computed by the script `exercise_checks.py`). The noncentrality doubles, and the error degrees of freedom more than double because the \( a+b-1 \) degrees of freedom spent on the mean structure are a fixed cost. Both effects raise the power.
:::

::: {#exr-tw-tukey-invariance}
[B2]

Show that Tukey's statistic \( F_N \) is unchanged if \( \y \) is replaced by \( \lambda\y+\mathbf{t} \), where \( \lambda\ne0 \) and \( \mathbf{t} \) is any additive table \( t_{ij}=u_i+v_j \). Conclude that \( F_N \) has the same distribution under every additive model with the same \( \sigma^2 \).
:::

::: {#exr-tw-tukey-power}
[B3]

Suppose \( y^p \) is additive for some \( p\ne0 \), so that \( \mu_{ij}=f(\nu+a_i+b_j) \) with \( f(t)=t^{1/p} \). Verify that \( \theta=f''(\nu)/f'(\nu)^2=(1-p)/f(\nu) \), and derive @eq-tw-tukey-power. What does \( \hat{\theta}\bar{y}_{\cdot\cdot}\approx1 \) suggest, and why is \( p=0 \) the right reading?
:::

::: {.solution}
\( f'(t)=p^{-1}t^{1/p-1} \) and \( f''(t)=p^{-1}(p^{-1}-1)t^{1/p-2} \). So
\[
\frac{f''(\nu)}{f'(\nu)^2}=\frac{p^{-1}(p^{-1}-1)\nu^{1/p-2}}{p^{-2}\nu^{2/p-2}}=p\Bigl(\frac1p-1\Bigr)\nu^{-1/p}=\frac{1-p}{f(\nu)} .
\]
Solving \( \theta=(1-p)/f(\nu) \) for \( p \), and estimating \( \theta \) by \( \hat{\theta} \) and \( f(\nu) \), the typical size of a cell mean, by \( \bar{y}_{\cdot\cdot} \), gives \( \hat{p}=1-\hat{\theta}\bar{y}_{\cdot\cdot} \). A value \( \hat{\theta}\bar{y}_{\cdot\cdot}\approx1 \) gives \( \hat{p}\approx0 \). For \( f(t)=e^t \), which makes \( \log y \) additive, \( f''/f'^2=e^{-\nu}=1/f(\nu) \), which is the limit \( p\to0 \) of the power formula. This matches the convention that \( (y^p-1)/p\to\log y \).
:::

::: {#exr-tw-tukey-attenuation}
[B4]

With \( m=1 \) and normal errors, suppose the cell means are \( \mu_{ij}=\mu+\alpha_i+\beta_j+\theta\alpha_i\beta_j \) with \( \sum_i\alpha_i=\sum_j\beta_j=0 \), exactly the shape of @eq-tw-tukey-alternative. Show that \( r_i=\alpha_i+\bar{\varepsilon}_{i\cdot}-\bar{\varepsilon}_{\cdot\cdot} \), and that, given \( \bP_A\Y \) and \( \bP_B\Y \),
\[
\E\bigl(\hat{\theta}\mid\bP_A\Y,\bP_B\Y\bigr)=\theta\,\frac{\sum_ir_i\alpha_i}{\sum_ir_i^2}\cdot\frac{\sum_jc_j\beta_j}{\sum_jc_j^2}.
\]
Explain why \( \hat{\theta} \) is nearly unbiased when the row and column effects are large compared with \( \sigma \), and why it tends to be pulled towards zero when they are not. What does this imply for the power rule @eq-tw-tukey-power?
:::

::: {.solution}
The interaction terms \( \theta\alpha_i\beta_j \) have zero row and column means, so \( \bar{\mu}_{i\cdot}-\bar{\mu}_{\cdot\cdot}=\alpha_i \) and \( r_i=\alpha_i+\bar{\varepsilon}_{i\cdot}-\bar{\varepsilon}_{\cdot\cdot} \); likewise for \( c_j \). With \( \mathbf{U}=\mathbf{r}\otimes\mathbf{c} \) as in the proof of @thm-tw-additive(f), \( \hat{\theta}=\mathbf{U}\T\R\Y/\mathbf{U}\T\mathbf{U} \). By @thm-qf-orthogonal-projections, \( \R\Y \) is independent of \( (\bP_A\Y,\bP_B\Y) \), whatever the mean, so given those two vectors \( \mathbf{U} \) is fixed and \( \E(\R\Y)=\R\bmu \), the table \( (\theta\alpha_i\beta_j) \). Hence the conditional mean is \( \theta\sum_{i,j}r_ic_j\alpha_i\beta_j/(\sum_ir_i^2\sum_jc_j^2) \), which factorizes as displayed. If the effects dominate the noise, \( \mathbf{r}\approx\boldsymbol{\upalpha} \) and both factors are near one. If not, \( \sum_ir_i^2 \) carries the extra noise \( \sigma^2(a-1)/b \) on average while \( \sum_ir_i\alpha_i \) does not: the attenuation of a regression on a noisy regressor. The power rule then suggests a \( \hat{p} \) too close to one, too weak a transformation.
:::

### C. Going deeper

::: {#exr-tw-tukey-replicated}
[C1]

Extend Tukey's test to \( m\ge2 \) observations per cell. Replace \( y_{ij} \) by the cell mean \( \bar{y}_{ij\cdot} \) in @eq-tw-tukey-ss and multiply the result by \( m \), and compare it with the pure-error mean square. Show that under the additive model the resulting statistic has the \( F(1,ab(m-1)) \) distribution. Show also that \( \text{SS}_{AB}-\text{SS}_N \) gives a test of the remaining nonadditivity.
:::
