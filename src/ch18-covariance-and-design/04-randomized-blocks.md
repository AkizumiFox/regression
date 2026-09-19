# Randomized complete blocks

A covariate removes nuisance variation by regression. **Blocking** removes it by design: similar units (plots in one
part of a field, animals from one litter, strips cut from one board) are grouped into blocks, and treatments are compared
*within* blocks. This section analyses the randomized complete block design, explains how randomization justifies its
model, and measures what blocking gains.

## The design and its model

::: {#def-dsn-rcbd}
[Randomized complete block design]

There are \( t \) treatments and \( b \) blocks of \( t \) units each. In every block, the \( t \) treatments are assigned to the
\( t \) units by a random permutation, chosen uniformly and independently for different blocks. The response of the unit in block
\( j \) that received treatment \( i \) is written \( y_{ij} \).
:::

The usual model is the additive two-way model with one observation per cell,
\[
y_{ij}=\mu+\tau_i+\beta_j+\varepsilon_{ij},\qquad \varepsilon_{ij}\ \iid\ \Normal(0,\sigma^2),
\]{#eq-dsn-rcbd-model}

with treatment effects \( \tau_i \) and block effects \( \beta_j \). Its geometry is @exm-proj-two-factor and its analysis
of variance is that of @thm-tw-additive; we restate it in the form the design needs.

::: {#thm-dsn-rcbd}
[Analysis of the randomized complete block design]

Assume @eq-dsn-rcbd-model with \( t,b\ge2 \). Let \( \bar{y}_{i\cdot} \), \( \bar{y}_{\cdot j} \) and \( \bar{y}_{\cdot\cdot} \) be the treatment,
block and grand means.

::: {.enumerate options="label=(\alph*)"}
1. The sums of squares
   \[
   \begin{aligned}
   \text{SS}_T&=b\sum_i(\bar{y}_{i\cdot}-\bar{y}_{\cdot\cdot})^2,\qquad
   \text{SS}_B=t\sum_j(\bar{y}_{\cdot j}-\bar{y}_{\cdot\cdot})^2,\\
   \text{SSE}&=\sum_{i,j}(y_{ij}-\bar{y}_{i\cdot}-\bar{y}_{\cdot j}+\bar{y}_{\cdot\cdot})^2
   \end{aligned}
   \]
   are squared lengths of projections onto mutually orthogonal subspaces of dimensions \( t-1 \), \( b-1 \) and
   \( (t-1)(b-1) \), and with the correction for the mean they add up to \( \sum_{i,j}y_{ij}^2 \).

2. A function \( \sum_ic_i\tau_i \) is estimable iff \( \sum_ic_i=0 \). Its least squares estimate is \( \sum_ic_i\bar{y}_{i\cdot} \), and
   \[
   \sum_ic_i\bar{y}_{i\cdot}=\sum_ic_i\tau_i+\sum_ic_i\bar{\varepsilon}_{i\cdot}
   \]
   identically, whatever the block effects. Its variance is \( \sigma^2\sum_ic_i^2/b \).

3. The treatment \( F \) statistic satisfies
   \[
   F_T=\frac{\text{SS}_T/(t-1)}{\text{SSE}/((t-1)(b-1))}\sim F\bigl(t-1,(t-1)(b-1),\gamma\bigr),
   \]
   with
   \[
   \gamma=\frac{b}{\sigma^2}\sum_i(\tau_i-\bar{\tau})^2 .
   \]

4. \( \E(\text{SS}_B/(b-1))=\sigma^2+t\sum_j(\beta_j-\bar{\beta})^2/(b-1) \), and \( \E(\text{SSE}/((t-1)(b-1)))=\sigma^2 \).
:::

:::

::: {.proof}
Model @eq-dsn-rcbd-model is the additive two-way model of @thm-tw-additive with factor \( A \) the treatments, factor \( B \) the
blocks, \( a=t \) and \( m=1 \), so that \( n-a-b+1=(t-1)(b-1) \). Parts (a), (b) and (c) are @thm-tw-additive(a), (b) and (c), with (e)
for the form of the residual space. The only new statement is the identity in (b): averaging @eq-dsn-rcbd-model over \( j \) gives
\( \bar{y}_{i\cdot}=\mu+\tau_i+\bar{\beta}+\bar{\varepsilon}_{i\cdot} \), and \( \sum_ic_i(\mu+\bar{\beta})=0 \) for a contrast. Part (d) is
@thm-ss-expected-mean-squares applied to the decomposition in (a); the block projection sends the mean vector to the vector with
entries \( \beta_j-\bar{\beta} \), whose squared length is \( t\sum_j(\beta_j-\bar{\beta})^2 \).
:::

Part (b) is the reason to block. However large the block effects, they cancel exactly from every treatment comparison
and from the error mean square, which estimates the variation *within* blocks. The between-block variation moves into its own
line, \( \text{SS}_B \). A test of blocks is rarely of interest, since blocks were chosen because they differ; a large
\( \text{MS}_B \) measures what blocking achieved (@prp-dsn-crd-efficiency below).

Any block-by-treatment interaction is absorbed into the error, and without replication within cells it cannot be separated
from it. Tukey's one-degree-of-freedom test (@thm-tw-additive(f)) checks for its most common form. If treatment differences
genuinely vary between blocks, the comparisons refer to the average over the blocks used.

## Randomization as the justification of the model

The units of an experiment are not a sample from a normal population; they are the boards or patients at hand. What the
experimenter controls is the randomization, and it alone supports much of the analysis. Give each unit \( u \) a fixed value
\( \xi_u \), and assume **unit–treatment additivity**: unit \( u \) responds to treatment \( i \) with \( \xi_u+\tau_i \). The only
probability below is that of the randomization.

::: {#prp-dsn-randomization}
[Randomization moments of the block analysis]

Let the units of block \( j \) have values \( \xi_{j1},\dots,\xi_{jt} \), with mean \( \bar{\xi}_j \) and
\( W_j=\sum_s(\xi_{js}-\bar{\xi}_j)^2 \), and put
\[
S_W^2=\frac1{b(t-1)}\sum_jW_j ,
\]
the average within-block variance of the unit values. Under unit–treatment additivity and the randomization of @def-dsn-rcbd:

::: {.enumerate options="label=(\alph*)"}
1. every contrast estimate \( \sum_ic_i\bar{y}_{i\cdot} \) has randomization mean \( \sum_ic_i\tau_i \) and randomization variance
   \( S_W^2\sum_ic_i^2/b \);

2. \( \E(\text{SSE}/((t-1)(b-1)))=S_W^2 \);

3. \( \E(\text{SS}_T/(t-1))=S_W^2+b\sum_i(\tau_i-\bar{\tau})^2/(t-1) \).
:::

So the error mean square is unbiased for the variance that the model attributes to \( \sigma^2 \), and when the treatments
have no effect the treatment and error mean squares have the same expectation.
:::

::: {.proof}
Let \( \pi_j(i) \) be the unit of block \( j \) that receives treatment \( i \), and \( w_{ij}=\xi_{j\pi_j(i)}-\bar{\xi}_j \). Then
\( y_{ij}=\bar{\xi}_j+\tau_i+w_{ij} \). For fixed \( j \), the vector \( (w_{1j},\dots,w_{tj}) \) is a uniformly random permutation of the
deviations \( \xi_{js}-\bar{\xi}_j \), which sum to zero. Hence \( \E w_{ij}=0 \), \( \E w_{ij}^2=W_j/t \), and, since
\( \sum_kw_{kj}=0 \), \( \E w_{ij}w_{kj}=-W_j/\{t(t-1)\} \) for \( i\ne k \). Different blocks are independent.

(a) \( \sum_ic_i\bar{y}_{i\cdot}=\sum_ic_i\tau_i+b^{-1}\sum_j\sum_ic_iw_{ij} \), because \( \sum_ic_i=0 \) removes the \( \bar{\xi}_j \). The
mean follows. For the variance,
\[
\Var\Bigl(\sum_ic_iw_{ij}\Bigr)=\frac{W_j}{t}\sum_ic_i^2-\frac{W_j}{t(t-1)}\sum_{i\ne k}c_ic_k=\frac{W_j}{t-1}\sum_ic_i^2 ,
\]
since \( \sum_{i\ne k}c_ic_k=(\sum_ic_i)^2-\sum_ic_i^2=-\sum_ic_i^2 \). Summing over the independent blocks and dividing by
\( b^2 \) gives \( S_W^2\sum_ic_i^2/b \).

(b) The block means are \( \bar{y}_{\cdot j}=\bar{\xi}_j+\bar{\tau} \), and the treatment means are
\( \bar{y}_{i\cdot}=\bar{\xi}+\tau_i+\bar{w}_{i\cdot} \), with \( \bar{\xi}=b^{-1}\sum_j\bar{\xi}_j \). So the residual is
\( w_{ij}-\bar{w}_{i\cdot} \), and
\[
\text{SSE}=\sum_{i,j}w_{ij}^2-b\sum_i\bar{w}_{i\cdot}^2 .
\]
The first sum is \( \sum_jW_j \) for every randomization. By independence across blocks,
\( \E\bar{w}_{i\cdot}^2=b^{-2}\sum_j\E w_{ij}^2=b^{-2}\sum_jW_j/t \), so \( \E\bigl(b\sum_i\bar{w}_{i\cdot}^2\bigr)=b^{-1}\sum_jW_j \). Hence
\( \E(\text{SSE})=(1-1/b)\sum_jW_j=(t-1)(b-1)S_W^2 \).

(c) \( \text{SS}_T=b\sum_i(\tau_i-\bar{\tau}+\bar{w}_{i\cdot})^2 \), using \( \sum_i\bar{w}_{i\cdot}=0 \). The cross terms have mean zero, and
\( \E\bigl(b\sum_i\bar{w}_{i\cdot}^2\bigr)=(t-1)S_W^2 \) by the computation in (b).
:::

Independent measurement errors of variance \( \sigma_e^2 \) add \( \sigma_e^2 \) to every expectation. The proposition does not
make \( F_T \) exactly \( F \)-distributed, but it shows that the estimates and their estimated variances are unbiased because
of the randomization, whatever made the units differ. A systematic assignment could line treatments up with the \( \xi \)'s, and
then neither statement would hold.

**Randomization tests.**  Under the *sharp* null hypothesis that every unit would respond the same under every treatment,
the responses are fixed and only the labels are random. The proportion of the \( (t!)^b \) equally likely randomizations giving an
\( F_T \) at least as large as the observed one is an exact \( p \)-value, estimated by sampling randomizations when there are too many.
[Chapter 23](../ch23-resampling-inference/index.html) develops permutation tests for linear models.

## An example

::: {#exm-dsn-adhesives}
[Five wood adhesives on six boards]

In a synthetic experiment, five adhesives are compared for the shear strength of glued joints. Timber varies from board
to board, so each of six boards is cut into five strips, assigned to the adhesives by a random permutation. The strengths
(MPa) are:

| adhesive | board 1 | board 2 | board 3 | board 4 | board 5 | board 6 |
|---|---|---|---|---|---|---|
| 1 | 7.89 | 8.42 | 8.55 | 8.30 | 9.26 | 9.15 |
| 2 | 8.68 | 8.69 | 8.73 | 8.62 | 9.72 | 9.51 |
| 3 | 9.13 | 9.28 | 9.81 | 9.43 | 10.30 | 9.41 |
| 4 | 8.48 | 8.60 | 8.46 | 8.95 | 9.42 | 9.02 |
| 5 | 8.90 | 9.17 | 9.79 | 9.57 | 10.24 | 9.09 |

The analysis of variance is:

| source | df | sum of squares | mean square | \( F \) |
|---|---|---|---|---|
| adhesives | 4 | \( 4.106 \) | \( 1.0266 \) | \( 14.88 \) |
| boards | 5 | \( 4.068 \) | \( 0.8136 \) | \( 11.79 \) |
| error | 20 | \( 1.380 \) | \( 0.0690 \) | |
| corrected total | 29 | \( 9.555 \) | | |

The adhesives differ clearly (\( p=8.6\times 10^{-6} \)). The adhesive means are \( 8.595 \), \( 8.992 \),
\( 9.560 \), \( 8.822 \) and \( 9.460 \), and every difference has standard error
\( \sqrt{2\times0.0690/6}=0.152 \). Tukey's honestly significant difference at level \( 0.05 \),
\( q_{0.05}(5,20)\sqrt{0.0690/6} \), is \( 0.454 \) (@thm-mc-tukey applies, with the error mean square of the block
analysis, as in @exr-mc-tukey-two-way). It separates adhesives 3 and 5 from adhesives 1, 2 and 4: the significant pairs are 1–3, 1–5, 2–3, 2–5, 3–4 and 4–5.
Tukey's test for nonadditivity gives \( F=0.07 \) (\( p=0.80 \)), so there is no sign that the adhesive differences
depend on the strength of the board.

Of \( 20000 \) randomizations sampled from the \( 2.99\times 10^{12} \) possible plans, \( 1 \) gave an adhesive \( F \)
at least as large as \( 14.88 \): the randomization \( p \)-value is \( (1+1)/(20000+1)=0.00010 \), in line with the normal
theory.
:::

```{.python .run #cell-rcbd-anova}
import numpy as np
from scipy import stats

rng = np.random.default_rng(606)
t, b = 5, 6                                           # adhesives, boards
plan = np.array([rng.permutation(t) for _ in range(b)])   # which strip of each board gets which adhesive
tau = np.array([0.0, 0.40, 0.65, 0.15, 0.80])         # adhesive effects (unknown in practice)
board = rng.normal(0, 0.75, b)                        # board effects
y = np.round(9 + tau[:, None] + board[None, :] + rng.normal(0, 0.3, (t, b)), 2)

grand = y.mean()
trt_mean, blk_mean = y.mean(axis=1), y.mean(axis=0)
resid = y - trt_mean[:, None] - blk_mean[None, :] + grand
ss_trt = b * np.sum((trt_mean - grand) ** 2)
ss_blk = t * np.sum((blk_mean - grand) ** 2)
ss_err = np.sum(resid ** 2)
df_err = (t - 1) * (b - 1)
ms_trt, ms_blk, ms_err = ss_trt / (t - 1), ss_blk / (b - 1), ss_err / df_err
F_trt = ms_trt / ms_err
print(f"treatments SS {ss_trt:.2f}, blocks SS {ss_blk:.2f}, error SS {ss_err:.2f} on {df_err} df")
print(f"F for adhesives = {F_trt:.2f}, p = {stats.f.sf(F_trt, t - 1, df_err):.2g}")
print("adhesive means", trt_mean.round(2), " se of a difference", np.sqrt(2 * ms_err / b).round(3))
```

```{.python .run #cell-rcbd-randomization}
def f_stat(yy):
    tm, bm, gm = yy.mean(axis=1), yy.mean(axis=0), yy.mean()
    sst = b * np.sum((tm - gm) ** 2)
    sse = np.sum((yy - tm[:, None] - bm[None, :] + gm) ** 2)
    return (sst / (t - 1)) / (sse / df_err)

rng_perm = np.random.default_rng(7)
draws = 20_000
F_perm = np.empty(draws)
for r in range(draws):
    # re-randomize: permute the adhesive labels independently within every board
    perm = np.array([rng_perm.permutation(t) for _ in range(b)]).T
    F_perm[r] = f_stat(np.take_along_axis(y, perm, axis=0))
p_rand = (1 + np.sum(F_perm >= F_trt)) / (draws + 1)
print(f"randomization p-value {p_rand:.4f}  (F-test p-value {stats.f.sf(F_trt, t - 1, df_err):.2g})")
```

## How much did blocking gain?

Had the adhesives been assigned to the thirty strips completely at random, the completely randomized design (CRD) would
have had the board-to-board variation in its error. The block experiment itself can estimate how large that error would have
been, using a fact about sampling from a finite population.

::: {#lem-dsn-srs}
[Sample variance without replacement]

Let \( u_1,\dots,u_N \) be fixed numbers with mean \( \bar{u} \) and \( S^2=\sum_a(u_a-\bar{u})^2/(N-1) \). If \( \mathcal S \) is a simple
random sample of \( m\ge2 \) of the indices, drawn without replacement, then
\[
\E\sum_{a\in\mathcal S}(u_a-\bar{u}_{\mathcal S})^2=(m-1)S^2 ,
\]
where \( \bar{u}_{\mathcal S} \) is the sample mean.
:::

::: {.proof}
For any \( m \) numbers, \( \sum_a(u_a-\bar{u}_{\mathcal S})^2=\frac1{2m}\sum_{a,c}(u_a-u_c)^2 \), the sum being over ordered pairs of sample
members. The \( m(m-1) \) ordered pairs with \( a\ne c \) are each a uniformly random ordered pair of distinct population
members. For such a pair, \( \E(u_a-u_c)^2=\frac{1}{N(N-1)}\sum_{a\ne c}(u_a-u_c)^2=\frac{2N\sum_a(u_a-\bar{u})^2}{N(N-1)}=2S^2 \).
So the expectation is \( \frac1{2m}\,m(m-1)\,2S^2=(m-1)S^2 \).
:::

::: {#prp-dsn-crd-efficiency}
[What a completely randomized design would have given]

Let the \( N=bt \) units follow @eq-dsn-rcbd-model, with fixed block effects \( \beta_j \), and write
\( S_\beta^2=t\sum_j(\beta_j-\bar{\beta})^2/(N-1) \) for the variance of the block effects over the \( N \) units.

::: {.enumerate options="label=(\alph*)"}
1. If instead the treatments were assigned to the same units completely at random, \( b \) units each, the within-treatment
   mean square of the resulting one-way analysis would have expectation \( \sigma^2+S_\beta^2 \), over the errors and the
   randomization.

2. In the block design,
   \[
   \hat{\sigma}^2_{\text{CR}}=\frac{(b-1)\text{MS}_B+b(t-1)\text{MSE}}{bt-1}
   \]{#eq-dsn-crd-variance}

   is unbiased for \( \sigma^2+S_\beta^2 \).
:::

:::

::: {.proof}
(a) The \( b \) units given treatment \( i \) form a simple random sample of the \( N \) units. Their responses are
\( \mu+\tau_i+\beta_{j(u)}+\varepsilon_u \), so the within-treatment sum of squares of treatment \( i \) is that of the numbers
\( \beta_{j(u)}+\varepsilon_u \) over the sample. Given the assignment, the errors contribute \( (b-1)\sigma^2 \) in expectation,
and the cross terms have mean zero. By @lem-dsn-srs the block effects contribute \( (b-1)S_\beta^2 \) in expectation over the
assignment. Summing over the \( t \) treatments and dividing by \( t(b-1) \) gives \( \sigma^2+S_\beta^2 \).

(b) By @thm-dsn-rcbd(d), \( \E[(b-1)\text{MS}_B]=(b-1)\sigma^2+t\sum_j(\beta_j-\bar{\beta})^2 \) and
\( \E[b(t-1)\text{MSE}]=b(t-1)\sigma^2 \). The sum has expectation \( (bt-1)\sigma^2+(N-1)S_\beta^2 \).
:::

The **relative efficiency** \( \hat{\sigma}^2_{\text{CR}}/\text{MSE} \) estimates how many units a CRD would need per unit of the block
design for the same precision. Some authors multiply it by Fisher's factor \( (\nu_B+1)(\nu_C+3)/\{(\nu_B+3)(\nu_C+1)\} \), with
\( \nu_B=(t-1)(b-1) \) and \( \nu_C=t(b-1) \) the two error degrees of freedom, to allow for estimating \( \sigma^2 \) with fewer degrees
of freedom; this is a convention, not a theorem.

For the adhesives,
\[
\hat{\sigma}^2_{\text{CR}}=\frac{5\times0.8136+24\times0.0690}{29}=0.1974 ,
\]
and the relative efficiency is \( 0.1974/0.0690=2.86 \) (\( 2.81 \) with Fisher's factor). Blocking
was worth almost three times as many strips; without it the adhesive \( F \) would have been about \( 1.0266/0.1974=5.20 \).

```{.python .run #cell-rcbd-efficiency}
s2_crd = ((b - 1) * ms_blk + b * (t - 1) * ms_err) / (b * t - 1)   # estimated CRD error variance
re = s2_crd / ms_err
nu_b, nu_c = df_err, b * t - t                        # error df: blocked, completely randomized
re_fisher = re * (nu_b + 1) * (nu_c + 3) / ((nu_b + 3) * (nu_c + 1))
F_crd = ms_trt / s2_crd                               # what the treatment F would roughly have been
print(f"estimated CRD error variance {s2_crd:.2f} against {ms_err:.2f}:"
      f" relative efficiency {re:.2f} ({re_fisher:.2f} with Fisher's df factor)")
```

::: {.remark}
[Fixed or random blocks]

Treating the boards as a random sample of boards changes nothing about treatment comparisons in a complete block design,
since block effects cancel from every contrast (@thm-dsn-rcbd(b)). It matters for incomplete blocks
([Section 18.5](05-latin-squares-incomplete-blocks.html)) and for predicting a new block; Chapter 32 treats such mixed models.
:::

## Exercises

### A. Check your understanding

::: {#exr-dsn-rcbd-df}
[A1]

An experiment has \( t=4 \) treatments in \( b=8 \) blocks. Give the degrees of freedom of the block design's error, and of the
error of a completely randomized design with the same 32 units. How large must the relative efficiency (without Fisher's factor)
be for blocking to be worthwhile?
:::

::: {#exr-dsn-rcbd-hand}
[A2]

Verify from the table in @exm-dsn-adhesives that the treatment mean of adhesive 3 is \( 9.560 \), and compute the
difference between adhesives 3 and 5 with its standard error. Is it significant by Tukey's criterion?
:::

### B. Practice

::: {#exr-dsn-block-order}
[B1]

The strips cut from a board differ systematically: strips from the outer edge are stronger. How would you modify the design and
analysis to remove this source of variation as well? Which design of the next section does your answer lead to?
:::

::: {#exr-dsn-randomization-var}
[B2]

In @prp-dsn-randomization, show that \( \Cov\bigl(\sum_ic_i\bar{y}_{i\cdot},\sum_id_i\bar{y}_{i\cdot}\bigr)=S_W^2\sum_ic_id_i/b \) for two
contrasts. Conclude that orthogonal contrasts are uncorrelated over the randomization, as they are under the normal model.
:::

::: {.solution}
As in the proof, \( \Cov(\sum_ic_iw_{ij},\sum_kd_kw_{kj})=\frac{W_j}{t}\sum_ic_id_i-\frac{W_j}{t(t-1)}\sum_{i\ne k}c_id_k \), and
\( \sum_{i\ne k}c_id_k=(\sum_ic_i)(\sum_kd_k)-\sum_ic_id_i=-\sum_ic_id_i \). So the covariance is \( W_j\sum_ic_id_i/(t-1) \).
Summing over blocks and dividing by \( b^2 \) gives \( S_W^2\sum_ic_id_i/b \), which is zero for orthogonal contrasts.
:::

::: {#exr-dsn-crd-efficiency-proof}
[B3]

Suppose the block effects are random, \( \beta_j\ \iid\ \Normal(0,\sigma_\beta^2) \), independent of the errors. Show that
\( \E(\hat{\sigma}^2_{\text{CR}})=\sigma^2+\sigma_\beta^2\,t(b-1)/(bt-1) \), and interpret the factor \( t(b-1)/(bt-1) \).
:::

### C. Going deeper

::: {#exr-dsn-random-blocks}
[C1]

In @eq-dsn-rcbd-model let the block effects be random, \( \beta_j\ \iid\ \Normal(0,\sigma_\beta^2) \), independent of the errors.
Show that the treatment means \( \bar{y}_{i\cdot} \) are equicorrelated, with variance \( (\sigma^2+\sigma_\beta^2)/b \) and covariance
\( \sigma_\beta^2/b \), while every contrast estimate keeps the distribution given by @thm-dsn-rcbd(b). Show that
\( \E(\text{MS}_B)=\sigma^2+t\sigma_\beta^2 \), and that \( \bigl(\text{MS}_B+(t-1)\text{MSE}\bigr)/(tb) \) is unbiased for
\( \Var(\bar{y}_{i\cdot}) \). Why is \( s/\sqrt b \) the right scale for comparing treatments but the wrong standard error for a
single treatment mean?
:::

::: {.solution}
Averaging over blocks, \( \bar{y}_{i\cdot}=\mu+\tau_i+\bar{\beta}+\bar{\varepsilon}_{i\cdot} \), where \( \bar{\beta}\sim\Normal(0,\sigma_\beta^2/b) \) is
common to all treatments and the \( \bar{\varepsilon}_{i\cdot} \) are independent \( \Normal(0,\sigma^2/b) \) and independent of
\( \bar{\beta} \). Hence \( \Var(\bar{y}_{i\cdot})=(\sigma^2+\sigma_\beta^2)/b \) and \( \Cov(\bar{y}_{i\cdot},\bar{y}_{k\cdot})=\Var(\bar{\beta})=\sigma_\beta^2/b \)
for \( i\ne k \). In a contrast \( \bar{\beta} \) is multiplied by \( \sum_ic_i=0 \), so \( \sum_ic_i\bar{y}_{i\cdot}=\sum_ic_i\tau_i+\sum_ic_i\bar{\varepsilon}_{i\cdot} \)
exactly as with fixed blocks. The block means are \( \bar{y}_{\cdot j}=\mu+\bar{\tau}+\beta_j+\bar{\varepsilon}_{\cdot j} \), and the
\( \beta_j+\bar{\varepsilon}_{\cdot j} \) are independent with variance \( \sigma_\beta^2+\sigma^2/t \), so
\( \E(\text{SS}_B)=t(b-1)(\sigma_\beta^2+\sigma^2/t) \) and \( \E(\text{MS}_B)=\sigma^2+t\sigma_\beta^2 \). The error sum of squares is unchanged,
because the residual projection annihilates every vector that is constant within blocks, so \( \E(\text{MSE})=\sigma^2 \). Then
\[
\E\Bigl(\frac{\text{MS}_B+(t-1)\text{MSE}}{tb}\Bigr)=\frac{t\sigma^2+t\sigma_\beta^2}{tb}=\frac{\sigma^2+\sigma_\beta^2}b .
\]

A comparison of treatments uses only within-block information, in which \( \bar{\beta} \) cancels, so \( s/\sqrt b \) is its scale and the
Tukey intervals of @exr-mc-tukey-two-way remain exact. A single treatment mean is a statement about the population of blocks, and its
uncertainty includes the sampling of the \( b \) blocks: \( s^2/b \) misses the term \( \sigma_\beta^2/b \). The estimate above mixes two
mean squares with different degrees of freedom, so an interval built on it is only approximately \( t \); the proper treatment is the
subject of Chapter 32.
:::
