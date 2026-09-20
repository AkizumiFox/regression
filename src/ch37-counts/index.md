# Counts

A count has no upper bound, cannot be negative, and is nearly always more variable
when it is large than when it is small. The linear model fits none of these
descriptions and no transformation repairs all of them at once. The generalized
linear model does better: put the linear predictor on the log scale, so that the
mean is multiplicative, and let the distribution supply a variance that grows with
the mean. That is the Poisson log-linear model, the starting point of this
chapter.

It is only a starting point. Real counts are almost always more variable than the
Poisson allows, often have far more zeros, and are sometimes not a sample of
independent events at all but a contingency table in disguise. Each complaint has
a model attached to it: exposure and rates, the exact connection between Poisson
counts and multinomial proportions, the negative binomial, and the two standard
ways to model an excess of zeros.

**What you need.** [Chapter 34](../ch34-exponential-families-glm/index.html) (the
generalized linear model: @def-glm-model, the likelihood equations @thm-glm-score,
iteratively reweighted least squares @thm-glm-irls, the deviance @def-glm-deviance
and @thm-glm-deviance, and the three test statistics @prp-glm-three-tests).
[Chapter 35](../ch35-binary-responses/index.html) supplies the separation theorem
@thm-bin-separation, whose Poisson counterpart is @prp-cnt-existence, and
[Chapter 36](../ch36-multinomial-ordinal/index.html) the baseline-category logit
model @def-mlt-baseline, which [Section 37.3](03-poisson-multinomial.html) shows is a log-linear model in
disguise. From earlier parts: the two-way layout and its interaction (@def-tw-interaction
and @prp-tw-interaction-equiv), variance-stabilizing
transformations (@thm-tr-variance-stabilizing), the sandwich covariance
estimator (@thm-het-sandwich), the likelihood ratio test on the boundary of the
parameter space (@thm-mix-boundary), and model selection by AIC and BIC (@def-sel-aic-bic).


## Roadmap

- [Poisson regression](01-poisson-regression.html)
- [Offsets and rates](02-offsets-and-rates.html)
- [Poisson and multinomial](03-poisson-multinomial.html)
- [Negative binomial models](04-negative-binomial.html)
- [Zero inflation](05-zero-inflation.html)
- [Summary and notes](06-summary-and-notes.html)
