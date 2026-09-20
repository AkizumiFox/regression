# Identifiability

[Section 6.4](../ch06-projections/04-least-squares.html) ended with a table of four coefficient
vectors that disagreed in every entry, some even in sign, and still produced
the same fitted values. Nothing had gone wrong with the computation. The model
matrix had linearly dependent columns, and so the least squares problem had a whole affine
subspace of solutions. This chapter asks what can be learned about \( \bbeta \) in that
situation. We begin with a question that comes before estimation and does not involve data
at all. Suppose we knew the distribution of \( \Y \) exactly, as if we had infinitely many
replications. Would we then know \( \bbeta \)? If not, which features of \( \bbeta \) would we know?

## What the distribution of the data can reveal

Throughout, \( \Y=\X\bbeta+\be \) is the linear model of @def-lm-linear-model: \( \X \) is a fixed
\( n\times p \) matrix of rank \( r \), and \( \E(\be)=\bzero \). We do not assume \( r=p \). The
coefficient vector enters the model only through the product \( \X\bbeta \). Suppose the
distribution of the error vector \( \be \) is allowed to depend on other parameters, such as
\( \sigma^2 \), but not on \( \bbeta \). Then the distribution of \( \Y \) is the distribution of
\( \be \) shifted by the vector \( \X\bbeta \). Two coefficient vectors with the same product
\( \X\bbeta \) therefore give the same distribution of \( \Y \), with the other parameters held fixed.
Conversely, the distribution of \( \Y \) determines its mean \( \E(\Y)=\X\bbeta \). So the
most we can learn about \( \bbeta \) from the distribution is the value of \( \X\bbeta \). The
mean vector
\[
\bmu=\X\bbeta\in\C(\X)
\]
is the object the data speak about. A value of \( \bbeta \) is one way of writing \( \bmu \) in
coordinates.

This observation reduces the question to linear algebra. The map
\( \bbeta\mapsto\X\bbeta \) from \( \Real^p \) to \( \Real^n \) is linear. Two coefficient vectors have the same
image iff their difference lies in the null space \( \Null(\X) \). So
\[
\{\bb\in\Real^p:\X\bb=\X\bbeta\}=\bbeta+\Null(\X),
\]{#eq-est-fibre}

an affine subspace of dimension \( p-r \) by rank–nullity (@eq-proj-rank-nullity). All
vectors in @eq-est-fibre are observationally equivalent: no experiment with this model matrix
can distinguish between them.

::: {#def-est-identifiable}
[Identifiability]

In the linear model \( \E(\Y)=\X\bbeta \), \( \bbeta\in\Real^p \), the parameter \( \bbeta \) is
**identifiable** if \( \X\bbeta_1=\X\bbeta_2 \) implies \( \bbeta_1=\bbeta_2 \). A function
\( g \) defined on \( \Real^p \), with values in any set, is **identifiable** if
\( \X\bbeta_1=\X\bbeta_2 \) implies \( g(\bbeta_1)=g(\bbeta_2) \).
:::

The definition is phrased through the mean, not the whole distribution. By the argument
above the two versions agree whenever the error distribution does not involve \( \bbeta \),
which is the situation throughout this book until the generalized linear models of [Part VIII](../ch34-exponential-families-glm/index.html).
There the mean is \( h(\X\bbeta) \) for a known one-to-one function \( h \), and the same
definition applies unchanged. It says that \( g(\bbeta) \) is a well-defined function of the
point \( \bmu \), whichever coordinates were used to describe that point.

::: {#thm-est-identifiable-mean}
[Identifiable means “a function of the mean”]

A function \( g \) on \( \Real^p \) is identifiable iff there is a function \( h \) on \( \C(\X) \) with
\[
g(\bbeta)=h(\X\bbeta)\qquad\text{for every }\bbeta\in\Real^p .
\]
Equivalently, \( g \) is constant on each affine subspace \( \bbeta+\Null(\X) \).
:::

::: {.proof}
If \( g=h\circ\X \) and \( \X\bbeta_1=\X\bbeta_2 \), then \( g(\bbeta_1)=h(\X\bbeta_1)=h(\X\bbeta_2)=g(\bbeta_2) \).
Conversely, let \( g \) be identifiable. For \( \bmu\in\C(\X) \) choose any \( \bbeta \) with
\( \X\bbeta=\bmu \) and put \( h(\bmu)=g(\bbeta) \). The value does not depend on the choice,
because any two choices have the same image under \( \X \). By construction
\( g(\bbeta)=h(\X\bbeta) \). The last sentence restates the definition using @eq-est-fibre.
:::

The theorem makes identifiability a statement about information. A function
that can be computed from \( \bmu \) is one that could, in principle, be estimated
consistently, provided that \( \bmu \) itself can be estimated. A function that cannot
be computed from \( \bmu \) is one about which the model is silent. Its value can be anything
along the affine subspace, and more data from the same design cannot help.

## When the parameter is identifiable

::: {#prp-est-identifiable-rank}
The following are equivalent: (a) \( \bbeta \) is identifiable; (b) \( \Null(\X)=\{\bzero\} \);
(c) \( \rank(\X)=p \); (d) \( \X\T\X \) is nonsingular. In that case the function \( h \) of
@thm-est-identifiable-mean for \( g(\bbeta)=\bbeta \) is \( h(\bmu)=(\X\T\X)^{-1}\X\T\bmu \).
:::

::: {.proof}
By @eq-est-fibre, \( \bbeta \) is identifiable iff every affine subspace \( \bbeta+\Null(\X) \) is a
single point, which is (b). (b)\( \Leftrightarrow \)(c) is rank–nullity, and
(c)\( \Leftrightarrow \)(d) is @cor-proj-gram-colspace(c). If \( \X\T\X \) is nonsingular and
\( \bmu=\X\bbeta \), then \( (\X\T\X)^{-1}\X\T\bmu=(\X\T\X)^{-1}\X\T\X\bbeta=\bbeta \).
:::

When \( \rank(\X)=r<p \) the parameter is not identifiable, but parts of it may be. The
simplest parts are the individual coordinates.

::: {#prp-est-coordinate}
[Identifiability of a single coefficient]

Write \( \x_1,\dots,\x_p \) for the columns of \( \X \). For each \( j \) the following are equivalent:

::: {.enumerate options="label=(\alph*)"}
1. \( \beta_j \) is identifiable;

2. every \( \bv\in\Null(\X) \) has \( v_j=0 \);

3. \( \x_j \) is not a linear combination of the other columns of \( \X \);

4. deleting column \( j \) from \( \X \) lowers the rank.
:::

:::

::: {.proof}
(a)\( \Leftrightarrow \)(b): by @eq-est-fibre, the \( j \)th coordinates of the vectors observationally
equivalent to \( \bbeta \) are \( \beta_j+v_j \), \( \bv\in\Null(\X) \). They all equal \( \beta_j \) iff
\( v_j=0 \) for every null vector. (b)\( \Leftrightarrow \)(c): a relation
\( \x_j=\sum_{k\ne j}a_k\x_k \) is the same as a null vector with \( j \)th entry \( -1 \), and any
null vector with \( v_j\neq0 \) can be rescaled to have \( v_j=-1 \). (c)\( \Leftrightarrow \)(d): the
column space loses a dimension when \( \x_j \) is removed iff \( \x_j \) is not in the span of the
remaining columns.
:::

Part (c) gives a practical test. A coefficient is lost exactly when its regressor
is *aliased*, that is, when it can be rebuilt from the others. In a model that is only partly
identifiable, some coefficients can still be perfectly well identified. If \( \X=[\X_1,\x_j] \) and
\( \x_j\notin\C(\X_1) \), then \( \beta_j \) is identifiable however badly \( \X_1 \) is
parameterized. A typical case is a numeric covariate added to an overparameterized factor
model: the covariate's coefficient is identifiable as long as the covariate is not constant
within every level of the factor ([Section 8.6](06-several-factors.html)).

## Where rank deficiency comes from

Linear dependence among the columns of \( \X \) arises in three quite different ways, and it
helps to distinguish them.

*By construction.* A factor with \( g \) levels entered as an intercept together with all \( g \)
indicator columns gives \( g+1 \) columns spanning a \( g \)-dimensional space
(@exm-proj-oneway-rank). Such *overparameterized* models are written down deliberately,
because they treat the levels symmetrically. The dependence is harmless once it is understood,
and most of this chapter is about understanding it.

*By the logic of the variables.* Some sets of variables satisfy an identity that holds
for every observation. Examples are shares of a whole that add to one, a total entered
together with all its parts, and, the best-known case, age, period and birth cohort.

*By the design.* Two regressors that are logically unrelated may be exactly
collinear in the sample at hand. A factor level may never be observed. The design may
have fewer distinct rows than parameters. In each case another data set with the same
model could identify the parameter, but this one cannot.

Identifiability is therefore a property of the model matrix, and through it of the design,
not of the population being studied. The next example shows the second kind of
dependence.

::: {#exm-est-apc}
[Age, period and cohort]

A survey is fielded every five years from 2000 to 2020, and its respondents are grouped in
eleven five-year age bands from 20 to 70. For each of the \( 55 \) combinations we
have one average response. A demographer wants to separate three linear trends: in age
(people change as they grow older), in period (everyone changes with the times), and in birth
cohort (generations differ). With age, period and cohort centred at \( 45 \), \( 2010 \)
and \( 1965 \), the model is
\[
\E(y)=\beta_0+\beta_A\,\mathrm{age}+\beta_P\,\mathrm{period}+\beta_C\,\mathrm{cohort}.
\]
But cohort is period minus age, for every respondent. The four columns satisfy
\( \mathrm{age}-\mathrm{period}+\mathrm{cohort}=\bzero \), so \( \bv=(0,1,-1,1)\T \) spans \( \Null(\X) \)
and \( \rank(\X)=3 \). By @prp-est-coordinate, none of the three slopes is identifiable.

The listing fits the model three times, each time dropping one of the three columns. All
three fits have residual sum of squares \( 157.00 \), because they have the same
column space. They tell three different stories:

| Fit | age slope | period slope | cohort slope |
|---|---|---|---|
| no cohort term | 0.102 | 0.092 | 0.000 |
| no period term | 0.194 | 0.000 | 0.092 |
| no age term | 0.000 | 0.194 | -0.102 |

One analyst reports a steady ageing effect and no generational change. Another reports a
strong ageing effect and a positive cohort trend. A third reports no ageing at all, a strong
secular trend, and younger generations scoring lower. The data were in fact simulated with
slopes \( 0.30 \), \( -0.10 \) and \( 0.20 \). That vector lies on the parallel line \( \bbeta+t\bv \) of
equally valid truths, and the data cannot prefer it. Only the identifiable
combinations are recovered: the fitted \( 0.194 \) and \( 0.092 \) estimate the true \( 0.20 \) and \( 0.10 \).

What the three fits share are the linear functions orthogonal to \( \bv=(0,1,-1,1)\T \):
\[
\beta_A+\beta_P=0.194,\qquad \beta_P+\beta_C=0.092,\qquad \beta_A-\beta_C=0.102,
\]
in every row of the table, and \( \beta_0 \). Any statement about the separate trends
requires information from outside these data, for instance a constraint justified by subject
knowledge, and the conclusions then rest on that constraint (@exr-est-apc-factors).
[Figure 8.1.1](#fig-est-apc) shows how each fit divides the same fitted surface among the three
terms.
:::

::: {when-format="html"}
![**Figure 8.1.1.** Three equally good fits of the age–period–cohort model. Each panel
shows the contribution of one term to the fitted mean under the three fits of @exm-est-apc.
The sum of the three contributions is the same for every respondent, so the
data cannot choose between the fits.](apc_stories.svg){#fig-est-apc width=100%}
:::

::: {when-format="pdf"}
![Three equally good fits of the age–period–cohort model. Each panel
shows the contribution of one term to the fitted mean under the three fits of @exm-est-apc.
The sum of the three contributions is the same for every respondent, so the
data cannot choose between the fits.](apc_stories.pdf){width=100%}
:::

```{.python .run #cell-apc-design}
import numpy as np

ages = np.arange(20, 75, 5)                  # 11 age bands: 20, 25, ..., 70
periods = np.arange(2000, 2025, 5)           # 5 survey waves: 2000, ..., 2020
A, P = np.meshgrid(ages, periods, indexing="ij")
age, period = A.ravel() - 45.0, P.ravel() - 2010.0   # centred at 45 and 2010
cohort = period - age                        # birth year minus 1965, exactly
X = np.column_stack([np.ones(age.size), age, period, cohort])
print("n =", X.shape[0], " p =", X.shape[1], " rank =", np.linalg.matrix_rank(X))
print("X @ (0, 1, -1, 1) =", np.abs(X @ [0, 1, -1, 1]).max())
```

```{.python .run #cell-apc-stories}
rng = np.random.default_rng(8)
truth = np.array([50.0, 0.30, -0.10, 0.20])
y = X @ truth + rng.normal(scale=1.5, size=age.size)

def fit_without(j):
    """Least squares after deleting column j (one of age, period, cohort)."""
    keep = [k for k in range(4) if k != j]
    b = np.zeros(4)
    b[keep] = np.linalg.lstsq(X[:, keep], y, rcond=None)[0]
    return b

stories = {"no cohort term": fit_without(3),
           "no period term": fit_without(2),
           "no age term": fit_without(1)}
for name, b in stories.items():
    sse = np.sum((y - X @ b) ** 2)
    print(f"{name:15s} b = {np.round(b, 3)}  SSE = {sse:.3f}",
          f" age+period = {b[1] + b[2]:.3f}  period+cohort = {b[2] + b[3]:.3f}")
```

## Linear and nonlinear identifiable functions

In the example the identifiable linear functions were exactly those whose coefficient
vectors are orthogonal to the null space. That is true in general.

::: {#prp-est-linear-identifiable}
A linear function \( \blambda\T\bbeta \) is identifiable iff \( \blambda\perp\Null(\X) \), that is, iff
\( \blambda\in\C(\X\T) \). The identifiable linear functions form a subspace of dimension
\( r=\rank(\X) \).
:::

::: {.proof}
By @thm-est-identifiable-mean, \( \blambda\T\bbeta \) is identifiable iff
\( \blambda\T(\bbeta+\bv)=\blambda\T\bbeta \) for all \( \bbeta \) and all \( \bv\in\Null(\X) \), that is,
iff \( \blambda\T\bv=0 \) for all \( \bv\in\Null(\X) \). By @cor-proj-gram-colspace(d),
\( \Null(\X)\perpc=\C(\X\T) \), which has dimension \( \rank(\X\T)=r \).
:::

Identifiable functions need not be linear. If \( \blambda_1\T\bbeta \) and \( \blambda_2\T\bbeta \) are
identifiable, so are their ratio (where defined), their product, and \( \max\{\blambda_1\T\bbeta,0\} \).
Every function of \( \bmu \) is identifiable, for example \( \norm{\bmu}^2 \) or the largest mean. For
instance, in a quadratic regression \( \E(y)=\beta_0+\beta_1x+\beta_2x^2 \) with a full-rank design,
the location \( -\beta_1/(2\beta_2) \) of the turning point is identifiable. It is not linear, and
[Section 12.5](../ch12-intervals-and-bands/05-calibration.html) has to estimate it with care. The next section shows that among
the *linear* functions, identifiability coincides with a much more concrete property,
the existence of a linear unbiased estimator.

::: {.idea}
The data determine the mean vector \( \bmu=\X\bbeta \) and nothing else about \( \bbeta \). A function
of \( \bbeta \) is identifiable exactly when it is a function of \( \bmu \). A *linear* function
\( \blambda\T\bbeta \) is identifiable exactly when \( \blambda \) lies in the row space \( \C(\X\T) \),
the orthogonal complement of the invisible directions \( \Null(\X) \).
:::

## Exercises

### A. Check your understanding

::: {#exr-est-apc-basis}
[A1]

In @exm-est-apc, show directly that \( \beta_0 \), \( \beta_A+\beta_P \) and \( \beta_P+\beta_C \) are
identifiable, that every identifiable linear function is a combination of these three, and that
\( \beta_A-\beta_C \) is one of them. Is \( \beta_A+\beta_C \) identifiable?
:::

::: {#exr-est-shares}
[A2]

A household's spending is split into food, housing and other, with shares \( s_1+s_2+s_3=1 \).
A regression of savings on an intercept and all three shares is proposed. Find \( \Null(\X) \) (assuming
the shares vary enough that there is no other relation among the columns), decide which
coefficients are identifiable, and describe two standard ways to make the model identifiable.
:::

### B. Practice

::: {#exr-est-ancova-identifiable}
[B1]

Let \( \Z=[\bz_1,\dots,\bz_g] \) be the indicator matrix of a factor and \( \x \) a numeric covariate,
and consider \( \X=[\bone,\Z,\x] \). Show that the slope of \( \x \) is identifiable iff \( \x\notin\C(\Z) \),
that is, iff \( \x \) is not constant within every level. Show that no intercept-type parameter
\( \beta_0 \) or \( \alpha_k \) is identifiable.
:::

::: {.solution}
The columns \( \bone,\bz_1,\dots,\bz_g \) satisfy
\( \bone-\sum_k\bz_k=\bzero \), so they are dependent, and \( \C([\bone,\Z])=\C(\Z) \). By
@prp-est-coordinate(c), the slope is identifiable iff \( \x \) is not a combination of
\( \bone,\bz_1,\dots,\bz_g \), that is, iff \( \x\notin\C(\Z) \). A vector lies in \( \C(\Z) \) iff it is constant
within each level. For the other coefficients, \( \bv=(1,-1,\dots,-1,0)\T \) is a null vector with
nonzero entries in positions \( 0,1,\dots,g \), so by @prp-est-coordinate(b) none of
\( \beta_0,\alpha_1,\dots,\alpha_g \) is identifiable.
:::

::: {#exr-est-identifiable-cosets}
[B2]

Let \( \bP \) be the orthogonal projection onto \( \C(\X\T) \). Show that \( g \) is identifiable iff
\( g(\bbeta)=g(\bP\bbeta) \) for every \( \bbeta \). Deduce that the identifiable functions are exactly
the functions of \( \bP\bbeta \), and that \( \bP\bbeta \) is the unique coefficient vector in \( \C(\X\T) \)
representing the mean \( \X\bbeta \).
:::

::: {#exr-est-new-row}
[B3]

A design with model matrix \( \X \) of rank \( r<p \) is augmented by one more observation at a
point \( \x_0 \). Show that the rank increases iff \( \x_0\notin\C(\X\T) \). Show that, in the
age–period–cohort design, no additional observation can make the slopes identifiable, and explain
why in words.
:::

### C. Going deeper

::: {#exr-est-apc-factors}
[C1]

Replace the three linear terms of @exm-est-apc by three factors: an indicator for each age
band, each period and each cohort (with an intercept). Show that the null space of the resulting
model matrix has dimension \( 3 \) plus exactly one extra dimension. Identify the extra null vector: it
assigns the linear sequences \( a-\bar{a} \), \( -(t-\bar{t}) \) and \( c-\bar{c} \) to the age, period and
cohort effects. Deduce that every second difference of age effects, such as
\( \alpha_{a+1}-2\alpha_a+\alpha_{a-1} \), is identifiable, while no linear trend in any single one of
the three factors is.
:::

::: {.solution}
Each factor's indicators sum to \( \bone \), which gives the usual three
null vectors (intercept minus all indicators of one factor). The identity
\( \mathrm{age}-\mathrm{period}+\mathrm{cohort}=\text{const} \) holds for every observation. Writing the
age of an observation as \( \sum_a a\,z^A_a \), and similarly for period and cohort, it becomes a relation
\( \sum_a a\,\bz^A_a-\sum_t t\,\bz^P_t+\sum_c c\,\bz^C_c-\kappa\bone=\bzero \) among the columns. Subtracting
suitable multiples of the three standard relations centres the three sequences, and this gives the
stated vector. It is not a combination of the other three, because those are constant on each factor's
effects and this one is linear and nonconstant. To see that there are no more, check that the
matrix has rank \( 1+(A-1)+(P-1)+(C-1)-1 \). This rank check can be done directly or numerically for the
design of @exm-est-apc. A function \( \sum_a d_a\alpha_a \) is identifiable iff \( \mathbf{d} \) is orthogonal to
all four null vectors, that is, iff \( \sum_a d_a=0 \) and \( \sum_a d_a\,a=0 \). A second difference
satisfies both. A linear trend \( \sum_a(a-\bar{a})\alpha_a \) fails the second condition.
:::
