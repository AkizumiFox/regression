# Missing Data

Every data set in this book so far has arrived complete. Real ones do not. A
respondent skips the question about income, a patient stops coming to the clinic,
a sensor fails for three weeks, a laboratory loses a sample. The holes are rarely
scattered at random, and the software's silent default — drop every record with a
hole in it — is a decision, not a neutral act.

This chapter treats the holes as data. The pattern of what is missing has a
distribution of its own, and the subject turns on what that distribution may depend
on. If it depends on nothing, almost any analysis survives. If it depends only on
what was recorded, likelihood and Bayesian inference survive and most shortcuts do
not. If it depends on what was *not* recorded, nothing in the observed data can tell
you so. Those three cases organize what follows: when dropping incomplete records is
safe, why single imputation lies about precision, what EM is and what it does not
give you, how multiple imputation restores the uncertainty imputation destroys, and
how to report a conclusion resting on an untestable assumption.

**What you need.** [Chapter 5](../ch05-model-and-least-squares/index.html) and
[Chapter 6](../ch06-projections/index.html) (the linear model and least squares),
[Section 7.3](../ch07-optimality/03-maximum-likelihood.html) (maximum likelihood
@thm-opt-mle, the score and information @prp-opt-information),
[Section 7.6](../ch07-optimality/06-bayes-conjugate.html) and
[Section 12.6](../ch12-intervals-and-bands/06-bayes.html) (the conjugate
Bayesian linear model @thm-opt-bayes-conjugate and credible sets @def-ci-credible),
[Chapter 11](../ch11-general-linear-hypothesis/index.html) (the likelihood ratio
test @thm-glh-lrt). [Section
18.6](../ch18-covariance-and-design/06-missing-observations.html) (@thm-dsn-missing, filling in a lost plot of a designed experiment) is the special
case this chapter generalizes. [Chapter 23](../ch23-resampling-inference/index.html)
supplies the bootstrap (@def-bs-bootstrap), which is one of the two
honest routes to a standard error after EM;
[Chapter 25](../ch25-causal-interpretation/index.html) the habit of asking what a
coefficient is meant to answer; and
[Chapter 34](../ch34-exponential-families-glm/index.html) the generalized linear
model (@def-glm-model, @thm-glm-score), since nothing here is special to a normal
response. [Section 41.4](04-em.html) works out the EM algorithm for a normal response, and for that
it needs the conditional distribution of one normal block given another
(@thm-mvn-conditional, [Chapter 3](../ch03-multivariate-normal/index.html)). Only three earlier
results are cited inside a proof here: that one, @thm-opt-mle in
[Section 41.3](03-likelihood.html), and @thm-dsn-missing in [Section 41.4](04-em.html). The rest of
the list is a place to go back to — the linear model and its likelihood, the tests and credible sets
a completed data set would be analysed with, and the bootstrap and the causal question of
Chapters 23 and 25.

## Roadmap

- [Missingness mechanisms](01-mechanisms.html)
- [Complete cases and quick repairs](02-complete-case.html)
- [Likelihood under missingness at random](03-likelihood.html)
- [The EM algorithm](04-em.html)
- [Multiple imputation](05-multiple-imputation.html)
- [Sensitivity analysis](06-sensitivity.html)
- [Summary and notes](07-summary-and-notes.html)
