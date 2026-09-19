# Latin squares and incomplete blocks

A **Latin square** removes two crossed sources of nuisance variation at once. A **balanced incomplete block design**
allows blocks smaller than the number of treatments, at a computable price in precision. The Latin square is orthogonal;
the incomplete block design is not, and its treatment comparisons are adjusted for blocks exactly as they were adjusted for a
covariate in [Section 18.1](01-ancova.html).

## Latin squares

::: {#def-dsn-latin-square}
[Latin square design]

An \( a\times a \) **Latin square** is an arrangement of \( a \) symbols in an \( a\times a \) array such that each symbol occurs exactly once
in each row and exactly once in each column. A **Latin square design** compares \( a \) treatments on \( a^2 \) units classified by
two crossed blocking factors with \( a \) levels each (rows and columns). The treatment applied to the unit in row \( i \) and column
\( j \) is the symbol in cell \( (i,j) \) of a Latin square chosen at random.
:::

The cyclic square \( k(i,j)=i+j \) modulo \( a \) exists for every \( a \); the usual randomization permutes its rows, columns and
symbols at random, which makes every unit equally likely to receive every treatment. With \( k(i,j) \) the treatment in cell
\( (i,j) \), the model is additive in all three classifications:
\[
y_{ij}=\mu+\rho_i+\kappa_j+\tau_{k(i,j)}+\varepsilon_{ij},\qquad \varepsilon_{ij}\ \iid\ \Normal(0,\sigma^2).
\]{#eq-dsn-latin-model}

::: {#thm-dsn-latin-square}
[Analysis of the Latin square]

Assume @eq-dsn-latin-model with \( a\ge3 \). Write \( \bar{y}_{i\cdot} \), \( \bar{y}_{\cdot j} \) and \( \bar{y}^{(k)} \) for the means of row \( i \),
column \( j \) and treatment \( k \), each an average of \( a \) observations, and \( \bar{y} \) for the grand mean.

::: {.enumerate options="label=(\alph*)"}
1. The centred row, column and treatment spaces are mutually orthogonal. So \( \I \) is the sum of five mutually orthogonal
   projections: onto the mean, rows, columns, treatments (ranks \( 1,a-1,a-1,a-1 \)) and error (rank \( (a-1)(a-2) \)). The
   corresponding sums of squares are
   \[
   \begin{aligned}
   \text{SS}_R&=a\sum_i(\bar{y}_{i\cdot}-\bar{y})^2,\qquad \text{SS}_C=a\sum_j(\bar{y}_{\cdot j}-\bar{y})^2,\\
   \text{SS}_T&=a\sum_k(\bar{y}^{(k)}-\bar{y})^2 ,
   \end{aligned}
   \]
   and \( \text{SSE}=\sum_{i,j}\bigl(y_{ij}-\bar{y}_{i\cdot}-\bar{y}_{\cdot j}-\bar{y}^{(k(i,j))}+2\bar{y}\bigr)^2 \).

2. Every contrast \( \sum_kc_k\tau_k \) is estimable, with least squares estimate \( \sum_kc_k\bar{y}^{(k)} \), free of the row and column
   effects, and variance \( \sigma^2\sum_kc_k^2/a \).

3. The treatment \( F \) statistic satisfies
   \[
   F_T=\frac{\text{SS}_T/(a-1)}{\text{SSE}/((a-1)(a-2))}\sim F\bigl(a-1,(a-1)(a-2),\gamma\bigr),
   \]
   with
   \[
   \gamma=\frac{a}{\sigma^2}\sum_k(\tau_k-\bar{\tau})^2 .
   \]
:::

:::

::: {.proof}
(a) Take any two of the three factors, say rows and treatments. Each treatment occurs exactly once in each row, so the two-way
table of counts has every entry \( 1 \), with all row and column totals \( a \), and it satisfies the proportional frequencies condition @eq-ss-proportional. By @thm-ss-proportional the two centred spaces are orthogonal. The same argument applies to rows and columns
(each cell once) and to columns and treatments. By @thm-proj-sum the projection onto the model space is
\( \bP_0+(\bP_R-\bP_0)+(\bP_C-\bP_0)+(\bP_T-\bP_0) \), of rank \( 1+3(a-1) \). The error projection is \( \I \) minus this, with rank
\( a^2-3a+2=(a-1)(a-2) \). Applying the projections to \( \y \) gives the stated sums of squares and the residual.

(b) Treatment \( k \) occurs once in each row and once in each column, so averaging @eq-dsn-latin-model over its \( a \) cells gives
\( \bar{y}^{(k)}=\mu+\bar{\rho}+\bar{\kappa}+\tau_k+\bar{\varepsilon}^{(k)} \). For a contrast the constant cancels. Estimability follows as in @thm-dsn-rcbd(b), and by (a) and @thm-ss-orthogonal-design the least squares estimate is the one computed from the treatment means
alone. (c) is @cor-ss-ms-distributions.
:::

The price is the error: two degrees of freedom for \( a=3 \) and six for \( a=4 \), so small squares are often replicated. The
error consists of pieces of the row-by-column, row-by-treatment and column-by-treatment interactions, which cannot be told apart (@exr-dsn-latin-interaction). An interaction, such as a treatment working better in some rows, distorts the treatment comparisons:
the design is only as good as its additivity assumption.

::: {#exm-dsn-typists}
[Keyboard layouts, typists and sessions]

In a synthetic experiment, five typists each type a standard passage in five sessions, comparing five keyboard layouts.
Typists differ greatly, and speed rises with practice. A randomized \( 5\times5 \) Latin square gives each layout once to each
typist and once in each session:
\[
\begin{array}{c|ccccc}
 & \text{S1} & \text{S2} & \text{S3} & \text{S4} & \text{S5}\\\hline
\text{typist 1} & 2 & 3 & 1 & 4 & 5\\
\text{typist 2} & 3 & 5 & 4 & 2 & 1\\
\text{typist 3} & 5 & 1 & 2 & 3 & 4\\
\text{typist 4} & 4 & 2 & 5 & 1 & 3\\
\text{typist 5} & 1 & 4 & 3 & 5 & 2
\end{array}
\]
The analysis of variance of the speeds (words per minute) is:

| source | df | sum of squares | mean square | \( F \) |
|---|---|---|---|---|
| typists | 4 | \( 689.2 \) | \( 172.3 \) | \( 82.7 \) |
| sessions | 4 | \( 89.8 \) | \( 22.44 \) | \( 10.8 \) |
| layouts | 4 | \( 73.6 \) | \( 18.40 \) | \( 8.83 \) |
| error | 12 | \( 24.99 \) | \( 2.083 \) | |
| corrected total | 24 | \( 877.6 \) | | |

The layouts differ (\( p=0.0015 \)). Their means are \( 47.06 \), \( 50.10 \), \( 45.92 \),
\( 50.30 \) and \( 48.94 \) words per minute, and a difference has standard error
\( \sqrt{2\times2.083/5}=0.913 \).

Had the sessions been ignored, with layouts randomized within typists, the error variance would have been about
\( (\text{MS}_C+(a-1)\text{MSE})/a \) (@exr-dsn-latin-efficiency). Relative to the Latin square this is
\[
\frac{22.44+4\times2.083}{5\times2.083}=2.95 .
\]
Ignoring the typists instead would have cost a factor of \( 17.35 \), and no blocking at all a factor of
\( 16.25 \).
:::

```{.python .run #cell-latin-square-anova}
import numpy as np
from scipy import stats

rng = np.random.default_rng(62)
a = 5
cyclic = (np.arange(a)[:, None] + np.arange(a)[None, :]) % a    # a standard square
square = cyclic[rng.permutation(a)][:, rng.permutation(a)]      # permute rows and columns
square = rng.permutation(a)[square]                             # and relabel the layouts
typist = rng.normal(0, 8, a)                          # typists differ a lot
session = np.array([-3.0, -1.0, 0.0, 1.5, 2.5])       # practice effect over sessions
layout = np.array([0.0, 2.0, -1.0, 3.0, 1.0])
y = np.round(55 + typist[:, None] + session[None, :] + layout[square]
             + rng.normal(0, 1.5, (a, a)), 1)

grand = y.mean()
row_m, col_m = y.mean(axis=1), y.mean(axis=0)
trt_m = np.array([y[square == k].mean() for k in range(a)])
ss_row = a * np.sum((row_m - grand) ** 2)
ss_col = a * np.sum((col_m - grand) ** 2)
ss_trt = a * np.sum((trt_m - grand) ** 2)
ss_tot = np.sum((y - grand) ** 2)
ss_err = ss_tot - ss_row - ss_col - ss_trt
df_err = (a - 1) * (a - 2)
ms_err = ss_err / df_err
F_trt = ss_trt / (a - 1) / ms_err
print(f"SS rows {ss_row:.1f}, columns {ss_col:.1f}, layouts {ss_trt:.1f}, error {ss_err:.2f} on {df_err} df")
print(f"F for layouts = {F_trt:.2f}, p = {stats.f.sf(F_trt, a - 1, df_err):.2g}")
print("layout means", trt_m.round(2), " se of a difference", np.sqrt(2 * ms_err / a).round(3))
```

Because the sessions are ordered in time, a layout may affect speed on the next one. Such **carryover** effects are not in @eq-dsn-latin-model; *crossover designs*, such as Latin squares balanced for the preceding treatment, are built to estimate them.

## Balanced incomplete block designs

Often a block cannot hold every treatment: a taster can judge only a few samples, a batch suffices for only a few runs. With
blocks of size \( k<t \), treatment comparisons are no longer free of block effects, and the fairest arrangement treats every pair
alike.

::: {#def-dsn-bibd}
[Balanced incomplete block design]

A **balanced incomplete block design** (BIBD) with parameters \( (t,b,r,k,\lambda) \) arranges \( t \) treatments in \( b \) blocks of size
\( k<t \), so that each block contains \( k \) distinct treatments, each treatment appears in \( r \) blocks, and each pair of distinct
treatments appears together in exactly \( \lambda \) blocks. Its **incidence matrix** is the \( t\times b \) matrix \( \mathbf{N} \) with
\( n_{ij}=1 \) if treatment \( i \) is in block \( j \) and \( 0 \) otherwise. Treatments are assigned to the units within each block at random,
and the labels of treatments and blocks are randomized as well.
:::

::: {#lem-dsn-bibd}
[Identities for balanced incomplete block designs]

For a BIBD:

::: {.enumerate options="label=(\alph*)"}
1. \( bk=tr \) and \( \lambda(t-1)=r(k-1) \);

2. \( \mathbf{N}\mathbf{N}\T=(r-\lambda)\I_t+\lambda\mathbf{J}_t \), where \( \mathbf{J}_t=\bone\bone\T \);

3. \( r>\lambda \), and \( b\ge t \) (**Fisher's inequality**).
:::

:::

::: {.proof}
(a) Both sides of the first identity count the units. For the second, fix a treatment \( i \). The \( r \) blocks containing it hold
\( r(k-1) \) other units, and each of the other \( t-1 \) treatments occurs in exactly \( \lambda \) of these blocks.

(b) The \( (i,i') \) entry of \( \mathbf{N}\mathbf{N}\T \) counts the blocks containing both \( i \) and \( i' \). It is \( r \) on the diagonal and \( \lambda \) off it.

(c) \( k<t \) gives \( k-1<t-1 \), so \( \lambda=r(k-1)/(t-1)<r \). The matrix in (b) has eigenvalue \( r-\lambda+\lambda t \) on \( \bone \) and
\( r-\lambda \) on \( \bone\perpc \), all positive. So \( \mathbf{N}\mathbf{N}\T \) is nonsingular, and \( t=\rank(\mathbf{N}\mathbf{N}\T)\le\rank(\mathbf{N})\le b \).
:::

The model is the additive block model, observed only on the occupied cells:
\[
y_{ij}=\mu+\tau_i+\beta_j+\varepsilon_{ij}\qquad\text{for the } bk \text{ pairs with } n_{ij}=1 .
\]{#eq-dsn-bibd-model}

Let \( T_i \) be the total of the observations on treatment \( i \), \( B_j \) the total of block \( j \), and define the **adjusted
treatment totals**
\[
Q_i=T_i-\frac1k\sum_jn_{ij}B_j ,
\]{#eq-dsn-bibd-q}

treatment \( i \)'s total minus the average level of the blocks it happened to be in. Then \( \sum_iQ_i=0 \).

::: {#thm-dsn-bibd}
[Intra-block analysis of a balanced incomplete block design]

Assume @eq-dsn-bibd-model for a BIBD, with \( \Cov(\Y)=\sigma^2\I \).

::: {.enumerate options="label=(\alph*)"}
1. \( \sum_ic_i\tau_i \) is estimable iff \( \sum_ic_i=0 \). Its least squares estimate and variance are
   \[
   \sum_ic_i\hat{\tau}_i=\frac{k}{\lambda t}\sum_ic_iQ_i,\qquad
   \Var\Bigl(\sum_ic_i\hat{\tau}_i\Bigr)=\frac{k}{\lambda t}\,\sigma^2\sum_ic_i^2 .
   \]{#eq-dsn-bibd-estimate}

2. The sum of squares for treatments adjusted for blocks is \( \text{SS}_{T\mid B}=\frac{k}{\lambda t}\sum_iQ_i^2 \) on \( t-1 \) degrees of
   freedom. The residual sum of squares is
   \( \text{SSE}=\text{SS}_{\text{tot}}-\text{SS}_B-\text{SS}_{T\mid B} \), where \( \text{SS}_B=\sum_jB_j^2/k-(\sum_jB_j)^2/(bk) \) is the
   unadjusted block sum of squares, on \( bk-b-t+1 \) degrees of freedom.

3. Under normality,
   \[
   F=\frac{\text{SS}_{T\mid B}/(t-1)}{\text{SSE}/(bk-b-t+1)}\sim F(t-1,\,bk-b-t+1,\,\gamma),
   \]
   with
   \[
   \gamma=\frac{\lambda t}{k\sigma^2}\sum_i(\tau_i-\bar{\tau})^2 .
   \]

4. A randomized complete block design with \( r \) blocks estimates each contrast with variance \( \sigma^2\sum_ic_i^2/r \). The ratio of that
   variance to @eq-dsn-bibd-estimate is the **efficiency factor**
   \[
   E=\frac{\lambda t}{rk}=\frac{t(k-1)}{(t-1)k}<1 .
   \]
:::

:::

::: {.proof}
Let \( \X_B \) be the \( n\times b \) indicator matrix of blocks (with the intercept in its span) and \( \Z \) the \( n\times t \) indicator matrix of
treatments, with \( n=bk \). Let \( \M_B \) project onto \( \C(\X_B) \), which replaces every observation by its block mean, and put
\( \tilde{\Z}=(\I-\M_B)\Z \). By @lem-ss-residualized, \( \C([\X_B,\Z])=\C(\X_B)\dirsum\C(\tilde{\Z}) \), orthogonally. This is the structure
of @thm-dsn-ancova with the treatment indicators in the role of covariates, except that \( \tilde{\Z} \) is not of full rank: \( \Z\bone=\bone \)
lies in \( \C(\X_B) \), so \( \tilde{\Z}\bone=\bzero \).

The matrix \( \mathbf{C}=\tilde{\Z}\T\tilde{\Z}=\Z\T\Z-\Z\T\M_B\Z \) is the *information matrix* of the design. Here \( \Z\T\Z=r\I \), and
\( \Z\T\M_B\Z=\mathbf{N}\mathbf{N}\T/k \), because \( \Z\T\X_B=\mathbf{N} \) and \( \X_B\T\X_B=k\I \). By @lem-dsn-bibd(b) and (a),
\[
\mathbf{C}=r\I-\frac{(r-\lambda)\I+\lambda\mathbf{J}}{k}=\frac{r(k-1)+\lambda}{k}\I-\frac{\lambda}{k}\mathbf{J}=\frac{\lambda t}{k}\Bigl(\I-\frac1t\mathbf{J}\Bigr).
\]
Also \( \tilde{\Z}\T\y=\Z\T\y-\Z\T\M_B\y=\mathbf{Q} \), since \( \Z\T\y \) is the vector of treatment totals and \( \Z\T\M_B\y \) has entries
\( \sum_jn_{ij}B_j/k \).

(a) A vector \( (\bzero,\mathbf{c}) \), zero on the block and intercept parameters, is estimable iff \( \mathbf{c}\T=\mathbf{a}\T\Z \) for some \( \mathbf{a} \) with
\( \mathbf{a}\T\X_B=\bzero \) (@thm-est-estimable-identifiable). Such an \( \mathbf{a} \) satisfies \( \mathbf{a}=(\I-\M_B)\mathbf{a} \), so
\( \mathbf{c}\in\C(\tilde{\Z}\T)=\C(\mathbf{C}) \), which is the space of contrasts because \( \mathbf{C} \) is a positive multiple of the centring projection.
Conversely, for a contrast \( \mathbf{c}=\mathbf{C}\mathbf{h} \), the vector \( \mathbf{a}=\tilde{\Z}\mathbf{h} \) works. Its least squares estimate is
\( \mathbf{a}\T\y=\mathbf{h}\T\mathbf{Q} \). With the generalized inverse \( \mathbf{C}\ginv=\frac{k}{\lambda t}\I \), which satisfies \( \mathbf{C}\mathbf{C}\ginv\mathbf{C}=\mathbf{C} \)
because \( \I-\mathbf{J}/t \) is idempotent, we have \( \mathbf{C}\mathbf{C}\ginv\mathbf{c}=\mathbf{c} \) for a contrast, so we can take \( \mathbf{h}=\mathbf{C}\ginv\mathbf{c} \), and the
estimate is \( \frac{k}{\lambda t}\mathbf{c}\T\mathbf{Q} \). Its variance is
\( \sigma^2\norm{\tilde{\Z}\mathbf{h}}^2=\sigma^2\mathbf{h}\T\mathbf{C}\mathbf{h}=\sigma^2\mathbf{c}\T\mathbf{C}\ginv\mathbf{c}=\frac{k}{\lambda t}\sigma^2\norm{\mathbf{c}}^2 \).

(b) The extra sum of squares for treatments after blocks is the squared length of the projection onto \( \C(\tilde{\Z}) \),
\( \y\T\tilde{\Z}\mathbf{C}\ginv\tilde{\Z}\T\y=\mathbf{Q}\T\mathbf{C}\ginv\mathbf{Q}=\frac{k}{\lambda t}\sum_iQ_i^2 \) (@lem-ss-residualized). The block sum of squares is
\( \norm{(\M_B-\bP_0)\y}^2 \), and the block, adjusted treatment and residual projections are mutually orthogonal and sum to \( \I-\bP_0 \) (@thm-proj-sum). The model has rank
\( b+t-1 \).

(c) is @thm-glh-f-test. The noncentrality is the same quadratic form evaluated at the mean, where
\( \tilde{\Z}\T\bmu=\mathbf{C}\boldsymbol{\uptau} \), giving \( \boldsymbol{\uptau}\T\mathbf{C}\mathbf{C}\ginv\mathbf{C}\boldsymbol{\uptau}/\sigma^2=\frac{\lambda t}{k}\sum_i(\tau_i-\bar{\tau})^2/\sigma^2 \).

(d) The ratio is \( (1/r)/(k/(\lambda t)) \), and \( \lambda=r(k-1)/(t-1) \) by @lem-dsn-bibd(a). Since \( k<t \), \( (k-1)/k<(t-1)/t \).
:::

The adjusted total \( Q_i \) compares treatment \( i \) with the blocks it appeared in, and
\( \E Q_i=\frac{\lambda t}{k}(\tau_i-\bar{\tau}) \) with no \( \beta \)'s left (@exr-dsn-bibd-eq): this is the **intra-block analysis**.
The efficiency factor is the fraction of information that survives the elimination of blocks: about one half for \( k=2 \), close to
one for \( k=t-1 \). Smaller blocks are often more homogeneous, and the smaller \( \sigma^2 \) can more than make up for \( E \). Every pair
sharing a block makes the design connected (@thm-est-connected), and balance gives every \( \tau_i-\tau_{i'} \) the same variance
\( 2k\sigma^2/(\lambda t) \). Unbalanced connected designs are analysed the same way with a different information matrix \( \mathbf{C} \).

::: {#exm-dsn-tasting}
[Seven recipes, fourteen tasters]

In a synthetic tasting experiment, seven recipes of a sports drink are rated on a 0–10 scale by \( 14 \) tasters, each of whom
can judge only three. The blocks are the lines of two Fano planes, the cyclic sets \( \{0,1,3\}+s \) and \( \{0,1,5\}+s \) modulo 7 for
\( s=0,\dots,6 \) (recipes numbered \( 0,\dots,6 \) in this description and \( 1,\dots,7 \) below). This is a BIBD with \( t=7 \), \( b=14 \), \( k=3 \),
\( r=6 \) and \( \lambda=2 \). There are \( 42 \) ratings and \( 22 \) error degrees of freedom.

Tasters use the scale very differently: the unadjusted block sum of squares is \( 69.72 \) out of a corrected total of
\( 80.71 \). The raw recipe means, \( 6.15 \), \( 6.10 \), \( 5.40 \), \( 6.60 \),
\( 6.30 \), \( 6.85 \) and \( 6.02 \), therefore mix recipe quality with the generosity of the tasters each recipe met.
The adjusted totals are \( Q=(-1.30,\ 0.90,\ -3.43,\ 3.13,\ 0.43,\ -1.50,\ 1.77) \), and the adjusted means
\( \bar{y}+\frac{k}{\lambda t}Q_i \) are \( 5.92 \), \( 6.40 \), \( 5.47 \), \( 6.87 \), \( 6.30 \),
\( 5.88 \) and \( 6.58 \). Recipe 6 had the highest raw mean. Its adjusted mean is
\[
6.202+\tfrac{3}{14}\times(-1.50)=5.88 ,
\]
among the lowest: recipe 6 happened to meet generous tasters, and the adjustment removes \( 0.97 \) from its mean.

The adjusted recipe sum of squares is \( 6.36 \), and the residual sum of squares is \( 4.64 \), so
\[
F=\frac{6.36/6}{4.64/22}=5.03
\]
on 6 and 22 degrees of freedom (\( p=0.0022 \)). A difference of two recipes has standard error
\( \sqrt{2\times3\times0.2108/14}=0.301 \). A complete block design with six tasters tasting all seven recipes, with the same
\( \sigma^2 \), would give \( 0.265 \); the variance ratio is the efficiency factor \( 7\times2/(6\times3)=0.7778 \). A
taster judging seven drinks would probably not rate as consistently, though.
:::

```{.python .run #cell-bibd-intrablock}
import numpy as np
from scipy import stats

t, k = 7, 3
blocks = [sorted((d + s) % t for d in base) for base in ([0, 1, 3], [0, 1, 5]) for s in range(t)]
b = len(blocks)                                       # 14 tasters
N = np.zeros((t, b), dtype=int)                       # incidence: recipe i tasted by taster j
for j, blk in enumerate(blocks):
    N[blk, j] = 1
r = N.sum(axis=1)[0]
lam = (N @ N.T)[0, 1]

rng = np.random.default_rng(77)
tau = np.array([0.0, 0.4, -0.3, 0.8, 0.1, -0.6, 0.5])   # recipe effects
taster = rng.normal(0, 1.2, b)                          # tasters use the scale differently
rows = [(i, j) for j, blk in enumerate(blocks) for i in blk]
recipe = np.array([i for i, j in rows])
block = np.array([j for i, j in rows])
y = np.round(6 + tau[recipe] + taster[block] + rng.normal(0, 0.5, len(rows)), 1)
n = len(y)

T_tot = np.bincount(recipe, weights=y, minlength=t)          # recipe totals
B_tot = np.bincount(block, weights=y, minlength=b)           # taster totals
Q = T_tot - N @ B_tot / k                                    # adjusted recipe totals
tau_hat = k * Q / (lam * t)                                  # intra-block estimates (sum to zero)
ss_trt_adj = k / (lam * t) * np.sum(Q ** 2)
ss_blk = np.sum(B_tot ** 2) / k - y.sum() ** 2 / n           # blocks, ignoring recipes
ss_tot = np.sum((y - y.mean()) ** 2)
ss_err = ss_tot - ss_blk - ss_trt_adj
df_err = n - b - t + 1
ms_err = ss_err / df_err
F = ss_trt_adj / (t - 1) / ms_err
se_diff = np.sqrt(2 * k * ms_err / (lam * t))
eff = lam * t / (r * k)                                      # efficiency factor
print("raw recipe means     ", (T_tot / r).round(2))
print("adjusted recipe means", (y.mean() + tau_hat).round(2))
print(f"adjusted recipe SS {ss_trt_adj:.3f}, error SS {ss_err:.3f} on {df_err} df,"
      f" F = {F:.2f}, p = {stats.f.sf(F, t - 1, df_err):.2g}")
print(f"se of a difference {se_diff:.3f}; efficiency factor {eff:.4f}")
```

::: {.remark}
[Recovering inter-block information]

The block totals, with \( \E B_j=k\mu+\sum_in_{ij}\tau_i+k\beta_j \), also carry treatment information. With random block effects
they give an independent *inter-block* estimate of every contrast, whose variance grows with \( \sigma_\beta^2 \). Combining the two
estimates by inverse-variance weights, as Yates proposed, is generalized least squares in a mixed model (Chapter 32). With blocks as
different as these tasters, the inter-block estimate adds little.
:::

## Exercises

### A. Check your understanding

::: {#exr-dsn-latin-count}
[A1]

Show that the cyclic array \( k(i,j)=i+j \bmod a \) is a Latin square, and that permuting its rows, its columns, or its symbols gives
a Latin square. For \( a=3 \), how many distinct Latin squares are there, and does the randomization described above reach all of them?
:::

::: {#exr-dsn-bibd-params}
[A2]

Check the identities of @lem-dsn-bibd for the design of @exm-dsn-tasting. Is there a BIBD with \( t=7 \), \( k=3 \) and \( b=7 \)? With
\( t=6 \), \( k=3 \), \( b=6 \)?
:::

::: {#exr-dsn-bibd-hand}
[A3]

From the numbers in @exm-dsn-tasting, compute the estimate of \( \tau_4-\tau_3 \) and its standard error, and decide whether recipe 4 is
better than recipe 3 at level \( 0.05 \), considered alone.
:::

### B. Practice

::: {#exr-dsn-latin-efficiency}
[B1]

In a Latin square, suppose the columns had not been used as a blocking factor, so that the treatments were randomized within rows only.
Following @prp-dsn-crd-efficiency, show that the expected error mean square of that design is \( \sigma^2+\sum_j(\kappa_j-\bar{\kappa})^2/(a-1) \),
and that \( (\text{MS}_C+(a-1)\text{MSE})/a \) is unbiased for it in the Latin square.
:::

::: {.solution}
Randomizing within rows, the units of row \( i \) are assigned to the treatments by a random permutation. Within a row, the column effects
\( \kappa_1,\dots,\kappa_a \) play the role of the unit values, with within-row variance \( S_\kappa^2=\sum_j(\kappa_j-\bar{\kappa})^2/(a-1) \),
the same for every row. By @prp-dsn-randomization(b), applied to the unit values \( \kappa_j+\varepsilon_{ij} \) and averaged over the errors, the
expected error mean square of the row-blocked design is \( \sigma^2+S_\kappa^2 \). In the Latin square,
\( \E\text{MS}_C=\sigma^2+aS_\kappa^2 \) and \( \E\text{MSE}=\sigma^2 \) (@thm-ss-expected-mean-squares), so
\( \E\{\text{MS}_C+(a-1)\text{MSE}\}/a=\sigma^2+S_\kappa^2 \).
:::

::: {#exr-dsn-latin-interaction}
[B2]

In a \( 4\times4 \) Latin square, the row-by-column interaction space of the \( 16 \) cells has dimension \( 9 \). Show that the treatment space
(centred) lies inside it, and that the error space is its orthogonal complement within it, of dimension \( 6 \). What does this imply about
the effect of a row-by-column interaction on the treatment comparisons?
:::

::: {.solution}
The row-by-column interaction space is the orthogonal complement, within \( \Real^{16} \), of the span of the mean, row and column
indicators; it has dimension \( 16-1-3-3=9 \). By @thm-dsn-latin-square(a) the centred treatment space is orthogonal to the mean, row and column
spaces, so it lies in this complement, with dimension \( 3 \). The error space is the orthogonal complement of the mean, row, column and
treatment spaces, so it is the rest of the interaction space, of dimension \( 6 \). A row-by-column interaction in the mean vector is a vector of
the \( 9 \)-dimensional space, and its component in the treatment space is indistinguishable from a treatment effect. Unless the interaction happens to
be orthogonal to the treatment space, it biases the treatment comparisons.
:::

::: {#exr-dsn-bibd-eq}
[B3]

Show that \( \E Q_i=\frac{\lambda t}{k}(\tau_i-\bar{\tau}) \) in a BIBD, directly from @eq-dsn-bibd-model and @lem-dsn-bibd.
:::

::: {.solution}
\( \E T_i=r\mu+r\tau_i+\sum_jn_{ij}\beta_j \) and \( \E B_j=k\mu+\sum_{i'}n_{i'j}\tau_{i'}+k\beta_j \). So
\[
\E Q_i=r\tau_i-\frac1k\sum_j\sum_{i'}n_{ij}n_{i'j}\tau_{i'}=r\tau_i-\frac1k\bigl(r\tau_i+\lambda\textstyle\sum_{i'\ne i}\tau_{i'}\bigr),
\]
the \( \mu \)'s and \( \beta \)'s cancelling. Writing \( \sum_{i'\ne i}\tau_{i'}=t\bar{\tau}-\tau_i \) gives
\( \E Q_i=\tau_i\bigl(r-(r-\lambda)/k\bigr)-\lambda t\bar{\tau}/k=\frac{\lambda t}{k}(\tau_i-\bar{\tau}) \), using \( r(k-1)+\lambda=\lambda t \).
:::

::: {#exr-dsn-graeco}
[B4]

A **Graeco-Latin square** superimposes two Latin squares of the same order so that every ordered pair of symbols occurs exactly once.
Write down one of order 3, state the model with rows, columns and two treatment factors, and give the analysis of variance table with its
degrees of freedom. How many error degrees of freedom remain for \( a=3 \), \( 4 \) and \( 5 \)?
:::

### C. Going deeper

::: {#exr-dsn-bibd-connected}
[C1]

For a general incomplete block design with incidence matrix \( \mathbf{N} \), replication \( r_i \) and block sizes \( k_j \), show that the
information matrix is \( \mathbf{C}=\diag(r_i)-\mathbf{N}\diag(1/k_j)\mathbf{N}\T \), that \( \mathbf{C}\bone=\bzero \), and that every treatment contrast is estimable
iff \( \rank(\mathbf{C})=t-1 \). Relate this to the connectedness of @thm-est-connected.
:::
