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
| 34–45 | | see `blueprint/book.yaml` | not started |

Deployed 2026-09-20 to https://regression.akizumifox.com (site repo AkizumiFox/regression, GitHub
Pages from `main`, DNS through Cloudflare). Republish with `./build.py deploy --push`.
Chapters 1–33 (Parts I–VII) are deployed; the combined PDF is 1,485 pages.

## Open forward promises

Chapters 1–33 point forward to unwritten chapters in prose ("Chapter 35 treats binary
responses"). Those references are plain text until the target chapter is published; when a
chapter is added, search `src/` for "Chapter N" and link it. The forward references to
Chapters 31–33 were linked when Part VII was written.
