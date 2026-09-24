# Multinomial and Ordinal Responses

A response with two categories has one logit. A response with \( c \) categories has
\( c-1 \) of them, and the modelling question is which \( c-1 \) comparisons to make. If the
categories are names — a brand, a party, a mode of travel — the natural comparisons are
with a fixed baseline, and the result is the *baseline-category logit model*, a
generalized linear model whose response is a vector of counts. If the categories are
ordered, that choice wastes information: a model ignoring the order needs \( (c-1)p \)
coefficients to say what a single slope may say. Ordinal models get back down to \( p \) by
cutting the scale and giving every cut the same coefficients. Where the cut is made is
the whole story. Cut *at* a level and you get the cumulative, or proportional odds,
model; cut *conditionally on having got that far* and you get the sequential, or
continuation-ratio, model. These answer different questions, and the data can tell them
apart.

**What you need.** [Chapter 34](../ch34-exponential-families-glm/index.html) (exponential
families, the likelihood equations, Fisher scoring, the deviance and the three test
statistics) and [Chapter 35](../ch35-binary-responses/index.html) (logistic regression, the
alternative links, separation, goodness of fit). From earlier parts:
[Chapter 11](../ch11-general-linear-hypothesis/index.html) for the likelihood ratio
test @thm-glh-lrt and the three flavours of test statistic @prp-glh-trinity,
[Chapter 21](../ch21-nonnormality-heteroscedasticity-serial/index.html) for weighted least
squares and the sandwich estimator (@thm-het-sandwich), and
[Chapter 29](../ch29-model-selection/index.html) for the information criteria (@def-sel-aic-bic).
Those three are named for continuity with the normal-theory chapters and are pointed at in the
text; the proofs here cite only Chapters 34 and 35. Nothing here needs normality of anything.

## Roadmap

- [Baseline-category logits](01-baseline-category-logits.html)
- [Cumulative logit models](02-cumulative-logit-models.html)
- [Proportional odds and its checks](03-proportional-odds.html)
- [Sequential models](04-sequential-models.html)
- [Summary and notes](05-summary-and-notes.html)
