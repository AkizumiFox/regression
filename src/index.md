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
| IV | drops \( \rank(\X)=p \) | Estimability, analysis of variance, designed structure |
| V | audits every assumption | What breaks, how to see it, what inference survives |
| VI | drops unbiasedness | Shrinkage, selection, regularization, \( p>n \) |
| VII | drops \( \Cov(\be)=\sigma^2\I \) | Generalized least squares, general Gauss–Markov, mixed models |
| VIII | drops normality of \( \Y \) | Exponential families and generalized linear models |
| IX | drops the linear predictor | Polynomials, splines, smoothing, additive models |
| X | drops modeling only the mean | Quantile regression and distributional models |

## Choices this book makes

**Geometry before distribution theory.**  Least squares is introduced
as a nearest-point problem ([Chapter 6](ch06-projections/index.html)) before any sampling theory.
Fitted values, residuals and sums of squares are projections. Coefficients are
coordinates, and coordinates depend on the basis you choose. This is the only
framing in which rank-deficient models need no special treatment.

**Estimability is an estimation idea.**  It appears in Part II, not as a
preliminary to the analysis of variance. Putting it there is what allows Part IV
to treat ANOVA as ordinary subspace comparison.

**Computation is part of the theory.**  Conditioning, QR and the
singular value decomposition come immediately after least squares
([Chapter 10](ch10-computation/index.html)), because numerical stability and statistical
stability are closely linked.

**The Bayesian treatment is a thread, not a quarantined chapter.**
Boxed “Bayesian thread” sections close Parts II, III, VII and VIII. Each
develops the Bayesian counterpart of the result just proved.

**Departures have a theory.**  Part V opens by deriving what omitted
variables, wrong covariance, non-normality, outliers and collinearity actually
do to the estimator. Diagnostics and remedies come after that.

**Material most regression texts omit.**  This book adds chapters on
bootstrap and permutation inference (Chapter 23),
errors in variables (Chapter 24), the causal
interpretation of coefficients (Chapter 25),
high-dimensional regression (Chapter 28) and missing data
(Chapter 41).

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
unbiased estimation to regularization. Parts VII–X (Chapters 31–45) climb the rest
of the ladder: correlated errors, non-normal responses, nonlinear predictors, and
models for entire conditional distributions. Chapters are published here as they
are written; the ones not yet on the site are named in plain text.

### Reading paths

- **A first course on linear models.** Chapters [5](ch05-model-and-least-squares/index.html)–[9](ch09-sums-of-squares/index.html), 11–12, 15, 19–20 and 29.
  Consult Part I as needed.
- **The theory spine.** Part I, then Chapters [6](ch06-projections/index.html), [7](ch07-optimality/index.html), [8](ch08-estimability/index.html), [9](ch09-sums-of-squares/index.html), 11 and 31.
  This is the shortest route to a full structural understanding of the linear model.
- **For implementers.** Chapters [5](ch05-model-and-least-squares/index.html), [6](ch06-projections/index.html), [10](ch10-computation/index.html) and 20, Part VI and Chapter 34.
- **A course on generalized models.** Chapters [5](ch05-model-and-least-squares/index.html), [6](ch06-projections/index.html) and 11, then Parts VII–X.

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
| \( \bbeta \), \( \hbeta \) | coefficient vector and a least squares estimate |
| \( \be \) | error vector, \( \be=\Y-\X\bbeta \) |
| \( \hY \), \( \he \) | fitted values \( \X\hbeta \) and residuals \( \Y-\hY \) |
| \( \sigma^2 \) | error variance |
| \( \Real^n \) | \( n \)-dimensional real space with \( \inner{\bu}{\bv}=\bu\T\bv \) |
| \( \C(\A) \), \( \Null(\A) \) | column space and null space of \( \A \) |
| \( \mathcal{S}\perpc \) | orthogonal complement of a subspace \( \mathcal S \) |
| \( \A\T \), \( \A\ginv \) | transpose; a generalized inverse (\( \A\A\ginv\A=\A \)) |
| \( \M \) | orthogonal projection onto \( \C(\X) \) |
| \( \Mo \) | orthogonal projection onto the column space of a reduced model |
| \( \rank \), \( \tr \), \( \diag \) | rank, trace, diagonal matrix |
| \( \bone \), \( \bzero \), \( \I \) | vector of ones, zero vector, identity |
| \( \norm{\cdot} \) | Euclidean norm unless stated otherwise |
| \( \E \), \( \Var \), \( \Cov \) | expectation, variance, covariance (matrix) |
| \( \Normal_n(\bmu,\bSigma) \) | \( n \)-variate normal distribution |
| \( \chi^2(r,\gamma) \), \( F(r,s,\gamma) \) | noncentral distributions, noncentrality \( \gamma \) |

Vectors are columns. Matrices and vectors are set in bold upright type and
scalars in italic. A bare “projection” always means an orthogonal
projection with respect to the Euclidean inner product. [Section 6.9](ch06-projections/09-inner-products.html)
covers the other kinds.
