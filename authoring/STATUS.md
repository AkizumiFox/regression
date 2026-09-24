# Book status

The plan is `blueprint/book.yaml` (45 chapters in 10 parts); each chapter's brief is in
`authoring/blueprints/`. Chapter order is checked against the prerequisite DAG by
`blueprint/check_order.py`.

| Ch | Dir | Topic | State |
|---|---|---|---|
| 1 | ch01-matrix-algebra | Matrix algebra for linear models | deployed |
| 2 | ch02-random-vectors | Random vectors | deployed |
| 3 | ch03-multivariate-normal | The multivariate normal distribution | deployed |
| 4 | ch04-quadratic-forms | Quadratic forms and their distributions | deployed |
| 5 | ch05-model-and-least-squares | The model, least squares and the normal equations | deployed |
| 6 | ch06-projections | Projections: the geometry of least squares | deployed (pilot) |
| 7 | ch07-optimality | Optimality: Gauss–Markov, maximum likelihood and minimum variance | deployed |
| 8 | ch08-estimability | Estimability, identifiability and indicator variables | deployed |
| 9 | ch09-sums-of-squares | Sums of squares and orthogonal decomposition | deployed |
| 10 | ch10-computation | Computation: QR, SVD and conditioning | deployed |
| 11 | ch11-general-linear-hypothesis | Testing models: the general linear hypothesis | deployed |
| 12 | ch12-intervals-and-bands | Confidence intervals, regions and prediction bands | deployed |
| 13 | ch13-multiplicity | Multiplicity: Scheffé, Tukey, Bonferroni and the false discovery rate | deployed |
| 14 | ch14-correlation-lack-of-fit-prediction | Correlation, lack of fit and prediction theory | deployed |
| 15 | ch15-anova-subspaces | Analysis of variance as subspace comparison | deployed |
| 16 | ch16-multiway-layouts | Multi-way layouts and interaction | deployed |
| 17 | ch17-unbalanced-data | Unbalanced data, cell means and empty cells | deployed |
| 18 | ch18-covariance-and-design | Analysis of covariance and designed experiments | deployed |
| 19 | ch19-theory-of-departures | What goes wrong: a theory of departures | deployed |
| 20 | ch20-residuals-leverage-influence | Residuals, leverage and influence | deployed |
| 21 | ch21-nonnormality-heteroscedasticity-serial | Non-normality, heteroscedasticity and serial correlation | deployed |
| 22 | ch22-transformations | Transformations | deployed |
| 23 | ch23-resampling-inference | Bootstrap and permutation inference | deployed |
| 24 | ch24-errors-in-variables | Errors in variables | deployed |
| 25 | ch25-causal-interpretation | Causal interpretation of coefficients | deployed |
| 26 | ch26-collinearity | Collinearity: diagnosis and consequences | deployed |
| 27 | ch27-shrinkage | Shrinkage: ridge, principal components and the lasso | deployed |
| 28 | ch28-high-dimensional | High-dimensional regression: p > n | deployed |
| 29 | ch29-model-selection | Model selection and prediction | deployed |
| 30 | ch30-regularization-boosting | Regularization and boosting | deployed |
| 31 | ch31-general-gauss-markov | The general Gauss–Markov model | deployed |
| 32 | ch32-linear-mixed-models | Random effects and linear mixed models | deployed |
| 33 | ch33-clustered-longitudinal-splitplot | Clustered, longitudinal and split-plot data | deployed |
| 34 | ch34-exponential-families-glm | Exponential families and the generalized linear model | deployed |
| 35 | ch35-binary-responses | Binary responses | deployed |
| 36 | ch36-multinomial-ordinal | Multinomial and ordinal responses | deployed |
| 37 | ch37-counts | Counts | deployed |
| 38 | ch38-quasi-likelihood | Quasi-likelihood and overdispersion | deployed |
| 39 | ch39-glms-in-practice-bayes | GLMs in practice; the Bayesian thread | deployed |
| 40 | ch40-glmm-gee | GLMMs and generalized estimating equations | deployed |
| 41 | ch41-missing-data | Missing data | deployed |
| 42 | ch42-polynomials-piecewise | Polynomials, piecewise fits and the road to splines | deployed |
| 43 | ch43-smoothing | Smoothing: kernels, splines and penalties | deployed |
| 44 | ch44-additive-models | Additive, geoadditive and structured additive models | deployed |
| 45 | ch45-quantile-gamlss | Quantile regression and GAMLSS | deployed |

All forty-five chapters were deployed on 2026-09-21 to https://regression.akizumifox.com
(site repo AkizumiFox/regression, GitHub Pages from `main`, DNS through Cloudflare). The site
carries 350 chapter pages and the combined PDF is 2,002 pages, with 3,195 labels and
8,744 references.
Republish with `./build.py deploy --push`.

## The book is complete

All forty-five chapters of `blueprint/book.yaml` are written. There are no unwritten
chapters left to point at, and no chapter promises one. The forward references that were
carried as plain text while Parts VII–X were being written have all been turned into links
or `@label` citations; `./build.py check` fails on any that were missed, and
`grep -rnE 'Chapter 4[2-5]|Chapters 4[2-5]' src` should show only links.

## Keeping it true

- Every computed number comes from `code/chNN/`. After changing a script, run
  `MPLBACKEND=Agg .venv/bin/python tools/check_numbers.py`; use `--accept` only when a value
  is meant to change.
- `./build.py check` must pass with no failures and no spelling warnings: new words go into
  `spelling.txt` (sorted case-insensitively, base forms only — possessives are stripped by
  the checker).
- New notation is recorded in `authoring/NOTATION.md` before it is used. Part IX fixes the
  basis matrix, knots, the penalty matrix, the smoothing parameter and the smoother matrix;
  Part X the check loss and the quantile level.
- `blueprint/coverage/chNN.md` records, source topic by source topic, where the book covers
  it; `blueprint/make_blueprints.py` regenerates the briefs from `book.yaml`.

## Known imperfections

A read-only audit of the finished book (six dimensions: notation, duplication,
cross-references, narrative, conventions, rigour) found 71 items. All 17 major ones are fixed, and the minor ones are either fixed or
recorded below. This is what is left, written down rather than hidden:

- **Exercise solutions.** Every section page now has all three groups (112 exercises were
  written for the 38 pages that were short, 105 of them with worked solutions). Across the
  book about half of the older exercises still carry no solution; the book only promises
  them for the core ones.
- **References.** Normalized: 117 entries in 36 chapters were corrected, journal names
  brought to one form, and works cited in several chapters made identical. 28 first names
  could not be confirmed and keep their initials rather than be guessed at.
- **Advisory warnings.** `./build.py check` prints 26 "names a result without a reference"
  warnings. The heuristic in `build/manifest.py` now drops the three unactionable classes
  (a title that reads as a common noun, a result proved later, a result on the same page),
  which took it from 186 to 26; what is left is a phrase used as a name where a citation
  would be an improvement rather than a necessity.
- **The dependency graph.** Fixed. It counts only backward citations now, so `graph.html`
  shows the transitive reduction again (60 of 456 dependencies drawn, 396 implied by a
  longer chain) and the "chapters cite each other in a circle" warning is gone. A forward
  pointer in prose -- chapter 1 naming a chapter 26 label, of which 120 chapter pairs have
  one -- is a signpost, not a dependency, and used to put an arrow back into the reading
  order; a cycle has no unique transitive reduction, so the page fell back to all 576 edges.
  Forward pointers in prose stay allowed. Only a proof may not point forward, and none does
  -- four did, and `./build.py check` now fails if one ever does again.
- **The forward-citation gate.** `tools/check_forward_deps.py` is kind-aware and wired into
  `./build.py check` beside `check_optional.py`. A proof, a proof idea or a claim may not
  cite a result the book proves later; an example, an exercise, a solution, a remark, a
  warning or a statement may, and 80 do. A citation whose kind was not recorded counts as
  load-bearing, so the gate stops rather than guesses. The four that were repaired: a proof
  of @prp-mat-orthogonal(c) that used the multiplicativity of the determinant from Section
  1.6 (the claim now lives only where it is proved, as @prp-mat-det(f)), and three Key Idea
  boxes -- in 19.2, 32.3 and 42.4 -- whose closing sentence named a later definition; each
  now links the later section instead of citing its label. `--signposts` lists the 80.
- **The preface's reading paths.** The six hand-written chapter lists are gone. Every one
  named too few sections: the theory spine was 4 short of its own proofs, the generalized
  models course 43. The preface points at the generated pages instead, and a new profile,
  `linear-model-theory` (48 sections), carries the theory spine, which none of the seven
  other profiles aimed at. Chapters 17, 22, 23, 26 and 36 to 40 are still on no profile's path.
- **Source hygiene, left deliberately.** 26 section files write `\mathbf{Z}`, `\mathbf{A}`
  and the like where a macro exists. Normalizing them is *not* wanted: the output is
  identical, and the macros carry book-wide meanings (`\Z` is the random-effect design
  matrix from Part VII on), so rewriting a generic bold Z in Chapter 2 as `\Z` would
  manufacture false hits for exactly the greps a notation audit runs.
