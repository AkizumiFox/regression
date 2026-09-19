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
| 19–45 | | see `blueprint/book.yaml` | not started |

Deployed 2026-09-19 to https://regression.akizumifox.com (site repo AkizumiFox/regression, GitHub
Pages from `main`, DNS through Cloudflare). Republish with `./build.py deploy --push`.
Chapters 1–18 (Parts I–IV) are deployed.

## Open forward promises

Chapters 1–18 point forward to unwritten chapters in prose ("Chapter 32 treats random and mixed effects"). Those references are plain text until the target chapter is published;
when a chapter is added, search `src/` for "Chapter N" and link it.
