# Book status

The plan is `blueprint/book.yaml` (45 chapters in 10 parts); each chapter's brief is in
`authoring/blueprints/`. Chapter order is checked against the prerequisite DAG by
`blueprint/check_order.py`.

| Ch | Dir | Topic | State |
|---|---|---|---|
| 1 | ch01-matrix-algebra | Matrix algebra for linear models | written |
| 2 | ch02-random-vectors | Random vectors | written |
| 3 | ch03-multivariate-normal | The multivariate normal distribution | written |
| 4 | ch04-quadratic-forms | Quadratic forms and their distributions | written |
| 5 | ch05-model-and-least-squares | The model, least squares and the normal equations | not started |
| 6 | ch06-projections | Projections: the geometry of least squares | written (pilot) |
| 7–45 | | see `blueprint/book.yaml` | not started |

Not deployed yet: the deploy target (AkizumiFox/regression, regression.akizumifox.com) is
configured but has never been pushed.

## Open forward promises

Chapters 1–4 and 6 point forward to unwritten chapters in prose ("Chapter 11 uses this space
for the F statistic"). Those references are plain text until the target chapter is published;
when a chapter is added, search `src/` for "Chapter N" and link it.
