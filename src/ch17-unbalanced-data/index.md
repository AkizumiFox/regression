# Unbalanced Data, Cell Means and Empty Cells

Balanced layouts give orthogonal decompositions, and then every sum of squares answers one
clear question. Real experiments lose plots and patients, and observational data are unbalanced
from the start. With unequal cell counts the same data produce several different "sums of squares
for \( A \)", each testing its own hypothesis. This chapter starts from the parameterization in which
nothing is ambiguous, a mean for every cell, and writes every hypothesis as linear constraints on
the cell means. That shows exactly what the Type I, II and III sums of squares test, why Type III
needs contrast codings, how the additive model is a constrained cell-means model, and which questions
survive when some cells are empty.

**What you need.** Projections and nested models from [Chapter 6](../ch06-projections/index.html) (@thm-proj-nested, @thm-proj-constraint-space). Estimability, factor codings and connectedness from
[Chapter 8](../ch08-estimability/index.html) (@thm-est-characterization, @def-est-coding,
@thm-est-connected, @prp-est-interaction). The definitions of the Type I, II and III sums of squares and
the first cell-mean results in [Chapter 9](../ch09-sums-of-squares/index.html) (@def-ss-sequential,
@def-ss-partial, @thm-ss-two-way-hypotheses, @thm-ss-hypothesis). Testable hypotheses and restricted
least squares from [Chapter 11](../ch11-general-linear-hypothesis/index.html) (@thm-glh-general-f, @prp-glh-restricted-ls). The one-way layout of [Chapter 15](../ch15-anova-subspaces/index.html)
and the balanced two-way layout of [Chapter 16](../ch16-multiway-layouts/index.html). Four earlier
facts are each used once: in [Section 17.1](01-cell-means.html), the Gauss–Markov theorem
(@thm-opt-gauss-markov) and the mean of a quadratic form (@thm-rv-quadform-mean); in
[Section 17.2](02-types-of-sums-of-squares.html), the algebra of Kronecker products
(@prp-mat-kronecker); and in [Section 17.3](03-constrained-models.html), the moments of least
squares (@thm-lm-moments).

## Roadmap

- [The cell-means model](01-cell-means.html)
- [Types of sums of squares and what they test](02-types-of-sums-of-squares.html)
- [Constrained models](03-constrained-models.html)
- [Empty cells and estimability](04-empty-cells.html)
- [Summary and notes](05-summary-and-notes.html)
