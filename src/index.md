# Preface

Regression is taught under at least five different descriptions, and each one is
correct. In one, regression is least squares together with the distribution
theory of quadratic forms. In another, it is orthogonal projection onto a
subspace, and everything else follows as a corollary. A third starts from
likelihood and treats the normal linear model as the Gaussian member of the
generalized linear family. A fourth defines regression as estimating a
conditional distribution, with linearity one structural assumption among many.
A fifth treats it as a data-analytic procedure whose assumptions must be
checked, one at a time, against the data.

Most textbooks pick one of these descriptions and stay with it. This book uses
all of them, in order, as rungs of a single ladder. Parts I–III set up the
classical linear model and its exact theory. After that, *each part relaxes
exactly one assumption*, and asks what survives, what must be rebuilt, and what
new tools the relaxation requires.

| Part | What changes | What we study |
|---|---|---|
| I | — | The apparatus: matrices, random vectors, the normal distribution, quadratic forms |
| II | — | \( \Y=\X\bbeta+\be \), \( \E(\be)=\bzero \), \( \Cov(\be)=\sigma^2\I \), as a projection problem |
| III | adds \( \be\sim\Normal(\bzero,\sigma^2\I) \) | Exact finite-sample inference |
| IV | the model matrix is rank deficient by construction | Analysis of variance as subspace comparison, multi-way layouts, unbalanced data, designed experiments |
| V | audits every assumption | What breaks, how to see it, what inference survives |
| VI | drops unbiasedness | Shrinkage, selection, regularization, \( p>n \) |
| VII | drops \( \Cov(\be)=\sigma^2\I \) | Generalized least squares, general Gauss–Markov, mixed models |
| VIII | drops normality of \( \Y \) | Exponential families and generalized linear models |
| IX | drops the linear predictor | Polynomials, splines, smoothing, additive models |
| X | drops modelling only the mean | Quantile regression and distributional models |

## Choices this book makes

**Geometry before distribution theory.**  Least squares is introduced
as a nearest-point problem ([Chapter 6](ch06-projections/index.html)) before any sampling theory.
Fitted values, residuals and sums of squares are projections. Coefficients are
coordinates, and coordinates depend on the basis you choose. This is the only
framing in which rank-deficient models need no special treatment.

**Estimability is an estimation idea.**  It appears in Part II, not as a
preliminary to the analysis of variance. Putting it there is what allows [Part IV](ch15-anova-subspaces/index.html)
to treat ANOVA as ordinary subspace comparison.

**Computation is part of the theory.**  Conditioning, QR and the
singular value decomposition come immediately after least squares
([Chapter 10](ch10-computation/index.html)), because numerical stability and statistical
stability are closely linked.

**The Bayesian treatment is a thread, not a quarantined chapter.**
Boxed “Bayesian thread” sections run through Parts II, III, VII and VIII —
[Section 7.6](ch07-optimality/06-bayes-conjugate.html),
[Section 12.6](ch12-intervals-and-bands/06-bayes.html),
[Section 32.6](ch32-linear-mixed-models/06-bayes.html) and
[Sections 39.3–39.4](ch39-glms-in-practice-bayes/03-bayesian-linear-revisited.html) — each
developing the Bayesian counterpart of the result just proved, and the thread
runs on into Part IX, where a roughness penalty turns out to be a normal prior
([Section 43.6](ch43-smoothing/06-mixed-model-and-bayes.html)).

**Departures have a theory.**  [Part V](ch19-theory-of-departures/index.html) opens by deriving what omitted
variables, wrong covariance, non-normality, outliers and collinearity actually
do to the estimator. Diagnostics and remedies come after that.

**Material most regression texts omit.**  This book adds chapters on
bootstrap and permutation inference ([Chapter 23](ch23-resampling-inference/index.html)),
errors in variables ([Chapter 24](ch24-errors-in-variables/index.html)), the causal
interpretation of coefficients ([Chapter 25](ch25-causal-interpretation/index.html)),
high-dimensional regression ([Chapter 28](ch28-high-dimensional/index.html)), missing data
([Chapter 41](ch41-missing-data/index.html)) and structured additive and distributional
regression ([Chapter 44](ch44-additive-models/index.html), [Chapter 45](ch45-quantile-gamlss/index.html)).

## Code

Every figure and every number obtained by computation is produced by a Python
script in the book's `code/` directory, and the scripts check their own claims
with assertions. The code in the text is taken from those scripts. On this
website each code block is a live cell: press **Run** (or Shift+Enter) to execute
it in your browser, edit it, and run it again. Nothing needs to be installed;
the first run downloads a Python environment with numpy, scipy, pandas and
statsmodels. Cells on a page share their variables, and running a cell first
runs the cells above it that have not run yet.

## Sources

The coverage of this book was planned against six standard texts:
Rencher and Schaalje (2008), Seber and Lee (2003),
Christensen (2020), Sen and Srivastava (1990),
Agresti (2015) and Fahrmeir et al. (2021).
Their chapters were mapped onto a single prerequisite graph, and the order of
this book respects that graph. Each chapter's *Notes and sources* section
names the treatments it drew on and the primary literature. The exposition,
proofs, examples and exercises are original to this book.

## How to Use This Book

### Parts and chapters

Parts I–VI (Chapters 1–30) are the classical linear model, complete in itself:
estimation, exact inference, designed experiments, diagnostics, and the move from
unbiased estimation to regularization. Parts [VII](ch31-general-gauss-markov/index.html)–[X](ch45-quantile-gamlss/index.html) (Chapters [31](ch31-general-gauss-markov/index.html)–[45](ch45-quantile-gamlss/index.html)) climb the rest
of the ladder: correlated errors, non-normal responses, nonlinear predictors, and
models for entire conditional distributions. All forty-five chapters are here;
each one links to the results it uses, so any of them can be entered directly.

### Reading paths

- **A first course on linear models.** Chapters [5](ch05-model-and-least-squares/index.html)–[9](ch09-sums-of-squares/index.html), [11](ch11-general-linear-hypothesis/index.html)–[12](ch12-intervals-and-bands/index.html), [15](ch15-anova-subspaces/index.html), [19](ch19-theory-of-departures/index.html)–[20](ch20-residuals-leverage-influence/index.html) and [29](ch29-model-selection/index.html).
  Consult Part I as needed.
- **The theory spine.** Part I, then Chapters [6](ch06-projections/index.html), [7](ch07-optimality/index.html), [8](ch08-estimability/index.html), [9](ch09-sums-of-squares/index.html), [11](ch11-general-linear-hypothesis/index.html) and [31](ch31-general-gauss-markov/index.html).
  This is the shortest route to a full structural understanding of the linear model.
- **For implementers.** Chapters [5](ch05-model-and-least-squares/index.html), [6](ch06-projections/index.html), [10](ch10-computation/index.html) and [20](ch20-residuals-leverage-influence/index.html), [Part VI](ch26-collinearity/index.html) and [Chapter 34](ch34-exponential-families-glm/index.html).
- **A course on generalized models.** Chapters [5](ch05-model-and-least-squares/index.html), [6](ch06-projections/index.html) and [11](ch11-general-linear-hypothesis/index.html), then Parts [VII](ch31-general-gauss-markov/index.html)–[X](ch45-quantile-gamlss/index.html).
- **Smoothing and flexible regression.** Chapters [6](ch06-projections/index.html), [27](ch27-shrinkage/index.html), [29](ch29-model-selection/index.html) and [30](ch30-regularization-boosting/index.html), then [Part IX](ch42-polynomials-piecewise/index.html);
  add [Chapter 32](ch32-linear-mixed-models/index.html) for the mixed-model reading of a penalty and
  [Chapter 34](ch34-exponential-families-glm/index.html) for generalized additive models.
- **Beyond the conditional mean.** Chapters [12](ch12-intervals-and-bands/index.html), [21](ch21-nonnormality-heteroscedasticity-serial/index.html) and [22](ch22-transformations/index.html), then [Chapter 45](ch45-quantile-gamlss/index.html).

### Prerequisites

Each chapter opens with a list of the earlier chapters it uses, taken from the
prerequisite graph the book was planned against. Part I is written as
reference material. Read it quickly the first time and come back to it when a
later chapter sends you there.

### Exercises

Exercises come at the end of each section in three groups: **A. Check your
understanding** (routine: a definition or computation has been understood),
**B. Practice** (core: the result is used later in the book) and **C. Going deeper**
(hard: may need an idea not in the text). Many have a worked solution, folded
under the exercise; try the exercise before opening it.

## Notation

| Symbol | Meaning |
|---|---|
| \( \Y \), \( \y \) | response vector (random, observed), \( n\times 1 \) |
| \( \X \) | model matrix, \( n\times p \); columns \( \x_1,\dots,\x_p \) |
| \( \mathbf{X} \) | a *random* vector of regressors, in the population sections [5.1](ch05-model-and-least-squares/01-what-a-model-claims.html), [6.11](ch06-projections/11-population.html) and [14.2](ch14-correlation-lack-of-fit-prediction/02-best-linear-prediction.html)–14.5 only |
| \( \bbeta \), \( \hbeta \) | coefficient vector and a least squares estimate |
| \( \be \) | error vector, \( \be=\Y-\X\bbeta \) |
| \( \hY \), \( \he \) | fitted values \( \X\hbeta \) and residuals \( \Y-\hY \) |
| \( \sigma^2 \) | error variance |
| \( \Real^n \) | \( n \)-dimensional real space with \( \inner{\bu}{\bv}=\bu\T\bv \) |
| \( \C(\A) \), \( \Null(\A) \) | column space and null space of \( \A \) |
| \( \mathcal{S}\perpc \) | orthogonal complement of a subspace \( \mathcal S \) |
| \( \A\T \), \( \A\ginv \) | transpose; a generalized inverse (\( \A\A\ginv\A=\A \)) |
| \( \G \) | a generalized inverse of \( \X\T\X \) (Chapters [6](ch06-projections/index.html)–[30](ch30-regularization-boosting/index.html)); from [Chapter 31](ch31-general-gauss-markov/index.html) on, \( \Cov(\bu) \) in a mixed model |
| \( \M \) | orthogonal projection onto \( \C(\X) \) |
| \( \Mo \) | orthogonal projection onto the column space of a reduced model |
| \( \rank \), \( \tr \), \( \diag \) | rank, trace, diagonal matrix |
| \( \bone \), \( \bzero \), \( \I \) | vector of ones, zero vector, identity |
| \( \norm{\cdot} \) | Euclidean norm unless stated otherwise |
| \( \E \), \( \Var \), \( \Cov \) | expectation, variance, covariance (matrix) |
| \( \Normal_n(\bmu,\bSigma) \) | \( n \)-variate normal distribution |
| \( \chi^2(r,\gamma) \), \( F(r,s,\gamma) \) | noncentral distributions, noncentrality \( \gamma \) |
| \( \V \) | error covariance: a known shape in \( \Cov(\Y)=\sigma^2\V \) ([Chapter 31](ch31-general-gauss-markov/index.html)), or \( \V=\Z\G\Z\T+\R \) in a mixed model (Chapters [32](ch32-linear-mixed-models/index.html)–[33](ch33-clustered-longitudinal-splitplot/index.html)) |
| \( \Z \), \( \bu \), \( \G \), \( \R \) | random-effect design matrix, random effects, \( \Cov(\bu) \) and \( \Cov(\be) \) ([Part VII](ch31-general-gauss-markov/index.html)) |
| \( \theta \), \( b(\theta) \), \( \phi \), \( V(\mu) \) | natural parameter, cumulant function, dispersion and variance function of an exponential dispersion family ([Part VIII](ch34-exponential-families-glm/index.html)) |
| \( \eta \), \( g \), \( \W \), \( D \) | linear predictor \( \x_{(i)}\T\bbeta \), link function, working weights and deviance of a generalized linear model ([Part VIII](ch34-exponential-families-glm/index.html)) |
| \( \R_i \), \( \R \) | working correlation matrix of a cluster ([Chapter 40](ch40-glmm-gee/index.html)) and missingness indicator array ([Chapter 41](ch41-missing-data/index.html)): local uses of \( \R \) — as are the residual projection of Chapter 16 and the constraint matrix of Chapter 17 — each announced where it starts |
| \( \B \), \( \bgamma \) | matrix of basis functions, \( B_{ij}=B_j(x_i) \), and its coefficient vector ([Part IX](ch42-polynomials-piecewise/index.html)) |
| \( \kappa_j \), \( t_j \), \( K \) | knots and their number: \( \kappa_1<\dots<\kappa_K \) are the interior knots in [Chapter 42](ch42-polynomials-piecewise/index.html); [Chapter 43](ch43-smoothing/index.html) writes \( t_1<\dots<t_K \) for the interior knots and keeps \( \kappa \) for the extended B-spline sequence; [Chapter 44](ch44-additive-models/index.html) writes \( \kappa_l \) for the knots of a kriging term |
| \( \bP \), \( \lambda \), \( \bS_\lambda \) | penalty matrix, smoothing parameter and smoother matrix of a penalized fit (Chapters [43](ch43-smoothing/index.html)–[45](ch45-quantile-gamlss/index.html)); \( \tr(\bS_\lambda) \) is its effective degrees of freedom |
| \( \rho_{\tau} \), \( \tau \) | the check loss \( \rho_{\tau}(u)=u\{\tau-\mathbf{1}\{u<0\}\} \) and the quantile level ([Chapter 45](ch45-quantile-gamlss/index.html)) |

Vectors are columns. Matrices and vectors are set in bold upright type and
scalars in italic. A bare “projection” always means an orthogonal
projection with respect to the Euclidean inner product. [Section 6.9](ch06-projections/09-inner-products.html)
covers the other kinds.
