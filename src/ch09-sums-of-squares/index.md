# Sums of Squares and Orthogonal Decomposition

An analysis of variance table looks like bookkeeping. Behind every row there is a subspace,
and behind the whole table an orthogonal decomposition of the space in which the data live:
sums of squares are squared lengths of projections, and degrees of freedom are dimensions.
Sequential and partial sums of squares, and the Type I, II and III sums of squares printed by
software, differ only in which subspaces they project onto, and so in which hypotheses they
address. When the subspaces are orthogonal, as in balanced designs and with orthogonal
polynomials, the distinctions disappear. The coefficient of determination is a ratio of two of
these squared lengths, and we mark what it does not measure. Finally, the expectation
of each sum of squares shows what an \( F \) ratio compares, preparing the tests of Part III.

**What you need.** [Chapter 6](../ch06-projections/index.html) throughout, especially
projections and their ranks (@thm-proj-sum, @prp-proj-trace-rank), nested models
(@thm-proj-nested, @thm-proj-constraint-space), the Frisch–Waugh–Lovell theorem (@thm-proj-fwl)
and \( R^2 \) as a squared cosine (@thm-proj-r2-cosine). The linear model and
the moments of least squares estimates from [Chapter 5](../ch05-model-and-least-squares/index.html)
(@def-lm-linear-model, @thm-lm-moments). The mean of a quadratic form (@thm-rv-quadform-mean)
and the distribution theory of [Chapter 4](../ch04-quadratic-forms/index.html)
(@thm-qf-orthogonal-projections, @thm-qf-cochran-algebra). Estimable functions and
contrasts from [Chapter 8](../ch08-estimability/index.html) appear in [Sections 9.2](02-extra-sums.html) to [9.4](04-orthogonal-designs.html).

## Roadmap

- [Decomposing the total sum of squares](01-total-decomposition.html)
- [Extra sums of squares](02-extra-sums.html)
- [Sequential and partial sums of squares](03-sequential-partial.html)
- [Orthogonal designs](04-orthogonal-designs.html)
- [The coefficient of determination and its limits](05-r-squared.html)
- [Expected mean squares](06-expected-mean-squares.html)
- [Summary and notes](07-summary-and-notes.html)
