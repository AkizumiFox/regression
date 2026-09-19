# The Multivariate Normal Distribution

A random vector is normal when every linear combination of its entries
is normal. From that definition come the moment generating function, the density, and
closure under all linear maps, singular ones included. Two further properties make the
normal distribution the working assumption of linear-model theory: zero covariance
between jointly normal blocks means independence, and conditioning one block on another
gives a linear regression with constant error variance. Partial correlation, the
correlation that remains after linear adjustment, connects these facts to least
squares.

**What you need.** [Chapter 1](../ch01-matrix-algebra/index.html) (spectral theorem, square roots, generalized and
partitioned inverses), [Chapter 2](../ch02-random-vectors/index.html) (mean vectors and covariance matrices).
Some examples refer forward to least squares (Chapter 5 and [Chapter 6](../ch06-projections/index.html)).

## Roadmap

- [Definition through linear combinations](01-definition.html)
- [Moment generating function and density](02-mgf-density.html)
- [Linear transformations and marginals](03-linear.html)
- [Independence and zero covariance](04-independence.html)
- [Conditional distributions](05-conditional.html)
- [Partial correlation](06-partial-correlation.html)
- [Summary and notes](07-summary-and-notes.html)
