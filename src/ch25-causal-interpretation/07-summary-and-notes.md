# Summary and Notes

## Summary

::: {.idea}

1. \( \E\bigl(Y\mid\operatorname{do}(X=x)\bigr) \), defined by a structural causal model (@def-cau-scm), differs from
   \( \E(Y\mid X=x) \) under confounding, at every sample size (@prp-cau-confounded).

2. Potential outcomes (@def-cau-potential-outcomes): the naive comparison is the effect on the treated plus a selection
   bias (@prp-cau-selection-bias).

3. Under complete randomization the difference in means is unbiased, with Neyman's variance; the estimate
   \( s_1^2/n_1+s_0^2/n_0 \) is conservative and equals the HC2 sandwich (@thm-cau-randomized).

4. An intervention deletes one factor of the Markov factorization (@prp-cau-markov); linear effects are sums over paths (@prp-cau-path-tracing);
   d-separation reads conditional independence off the graph (@prp-cau-three-structures, @thm-cau-d-separation).

5. Confounding is an open back-door path, and in linear models its bias is a sum over common ancestors (@prp-cau-confounding-bias).

6. An adjustment set identifies the effect; in linear models the coefficient of \( X \) is then the total effect (@thm-cau-backdoor).
   Causes of \( Y \) improve precision, pure causes of \( X \) worsen it (@prp-cau-adjustment-precision).

7. Mediators, common effects and pre-treatment colliders are bad controls (@thm-cau-bad-controls).

8. In experiments, analysis of covariance can lose to no adjustment, and the interacted regression is asymptotically
   optimal (@lem-cau-fixed-coefficients, @thm-cau-adjustment).

:::

## Notes and sources

**Structural models and graphs.** Wright (1921, 1934) introduced path analysis (@prp-cau-path-tracing), and Haavelmo (1943)
read structural equations as autonomous mechanisms that interventions change one at a time. The intervention as an
explicit operation and the truncated factorization (@prp-cau-markov) are from Pearl (1995), and in parallel the manipulation
theorem of Spirtes, Glymour and Scheines (2000; first edition 1993); Pearl (2009) is the systematic account. d-separation is due to Pearl
(1988), with proofs by Geiger, Verma and Pearl (1990) and Lauritzen, Dawid, Larsen and Leimer (1990); Spirtes, Glymour
and Scheines (2000) develop faithfulness.

**Potential outcomes.** @def-cau-potential-outcomes and @eq-cau-neyman-variance are Neyman's (1923), in a paper on field
trials translated in 1990. Rubin (1974) extended the framework to observational studies, Holland (1986) named the
fundamental problem of causal inference, and Imbens and Rubin (2015) set out the finite-population view. Finite-population central limit theorems go back to Hájek (1960); Li and Ding (2017) give forms
suited to causal inference. Samii and Aronow (2012) show that Neyman's estimate is the HC2 sandwich (@exr-cau-hc2).

**Adjustment and bad controls.** The back-door and front-door criteria are from Pearl (1995); Rosenbaum and Rubin (1983)
introduced strong ignorability and the propensity score; Shpitser, VanderWeele and Robins (2010) give a necessary and
sufficient adjustment criterion; and Henckel, Perković and Maathuis (2022) find the most efficient adjustment set in linear
models. @exr-cau-stratified is Angrist (1998)'s observation, and the "Table 2 fallacy" is named by Westreich and Greenland
(2013). Angrist and Pischke (2009) popularized the term bad control. Berkson (1946) described selection on a common effect,
Greenland, Pearl and Robins (1999) showed that pre-treatment colliders can bias adjustment, Ding and Miratrix (2015)
quantify M-bias and butterfly bias, Elwert and Winship (2014) survey collider bias, and Pearl (2010) analyses bias
amplification.

**Adjustment in experiments.** Covariance adjustment of experiments is Fisher's (1935). The design-based critique is
Freedman (2008a, 2008b), and the interacted estimator and its sandwich standard errors are Lin (2013). The route through the
exact variance of fixed-coefficient estimators (@lem-cau-fixed-coefficients), which reduces the comparison to completing a
square, is this book's.

## References

- Angrist, Joshua D. (1998). Estimating the Labor Market Impact of Voluntary Military Service Using Social Security Data on Military Applicants. *Econometrica* 66(2), 249–288.
- Angrist, Joshua D. and Pischke, Jörn-Steffen (2009). *Mostly Harmless Econometrics: An Empiricist's Companion*. Princeton, NJ: Princeton University Press.
- Berkson, Joseph (1946). Limitations of the Application of Fourfold Table Analysis to Hospital Data. *Biometrics Bulletin* 2(3), 47–53.
- Billingsley, Patrick (1995). *Probability and Measure*. 3rd edition. New York: Wiley.
- Ding, Peng and Miratrix, Luke W. (2015). To Adjust or Not to Adjust? Sensitivity Analysis of M-Bias and Butterfly-Bias. *Journal of Causal Inference* 3(1), 41–57.
- Elwert, Felix and Winship, Christopher (2014). Endogenous Selection Bias: The Problem of Conditioning on a Collider Variable. *Annual Review of Sociology* 40, 31–53.
- Fisher, Ronald A. (1935). *The Design of Experiments*. Edinburgh: Oliver and Boyd.
- Freedman, David A. (2008a). On Regression Adjustments to Experimental Data. *Advances in Applied Mathematics* 40(2), 180–193.
- Freedman, David A. (2008b). On Regression Adjustments in Experiments with Several Treatments. *The Annals of Applied Statistics* 2(1), 176–196.
- Geiger, Dan, Verma, Thomas and Pearl, Judea (1990). Identifying Independence in Bayesian Networks. *Networks* 20(5), 507–534.
- Greenland, Sander, Pearl, Judea and Robins, James M. (1999). Causal Diagrams for Epidemiologic Research. *Epidemiology* 10(1), 37–48.
- Haavelmo, Trygve (1943). The Statistical Implications of a System of Simultaneous Equations. *Econometrica* 11(1), 1–12.
- Hájek, Jaroslav (1960). Limiting Distributions in Simple Random Sampling from a Finite Population. *Publications of the Mathematical Institute of the Hungarian Academy of Sciences* 5, 361–374.
- Henckel, Leonard, Perković, Emilija and Maathuis, Marloes H. (2022). Graphical Criteria for Efficient Total Effect Estimation via Adjustment in Causal Linear Models. *Journal of the Royal Statistical Society, Series B* 84(2), 579–599.
- Holland, Paul W. (1986). Statistics and Causal Inference. *Journal of the American Statistical Association* 81(396), 945–960.
- Imbens, Guido W. and Rubin, Donald B. (2015). *Causal Inference for Statistics, Social, and Biomedical Sciences: An Introduction*. Cambridge: Cambridge University Press.
- Lauritzen, Steffen L., Dawid, A. Philip, Larsen, Birgitte N. and Leimer, Hanns-Georg (1990). Independence Properties of Directed Markov Fields. *Networks* 20(5), 491–505.
- Li, Xinran and Ding, Peng (2017). General Forms of Finite Population Central Limit Theorems with Applications to Causal Inference. *Journal of the American Statistical Association* 112(520), 1759–1769.
- Lin, Winston (2013). Agnostic Notes on Regression Adjustments to Experimental Data: Reexamining Freedman's Critique. *The Annals of Applied Statistics* 7(1), 295–318.
- Neyman, Jerzy (1923). On the Application of Probability Theory to Agricultural Experiments. Essay on Principles. Section 9. Translated from the Polish and edited by D. M. Dabrowska and T. P. Speed (1990), *Statistical Science* 5(4), 465–472.
- Pearl, Judea (1988). *Probabilistic Reasoning in Intelligent Systems: Networks of Plausible Inference*. San Mateo, CA: Morgan Kaufmann.
- Pearl, Judea (1995). Causal Diagrams for Empirical Research. *Biometrika* 82(4), 669–688.
- Pearl, Judea (2009). *Causality: Models, Reasoning, and Inference*. 2nd edition. Cambridge: Cambridge University Press.
- Pearl, Judea (2010). On a Class of Bias-Amplifying Variables that Endanger Effect Estimates. In *Proceedings of the Twenty-Sixth Conference on Uncertainty in Artificial Intelligence*. Corvallis, OR: AUAI Press.
- Rosenbaum, Paul R. and Rubin, Donald B. (1983). The Central Role of the Propensity Score in Observational Studies for Causal Effects. *Biometrika* 70(1), 41–55.
- Rubin, Donald B. (1974). Estimating Causal Effects of Treatments in Randomized and Nonrandomized Studies. *Journal of Educational Psychology* 66(5), 688–701.
- Samii, Cyrus and Aronow, Peter M. (2012). On Equivalencies between Design-Based and Regression-Based Variance Estimators for Randomized Experiments. *Statistics & Probability Letters* 82(2), 365–370.
- Shpitser, Ilya, VanderWeele, Tyler J. and Robins, James M. (2010). On the Validity of Covariate Adjustment for Estimating Causal Effects. In *Proceedings of the Twenty-Sixth Conference on Uncertainty in Artificial Intelligence*. Corvallis, OR: AUAI Press.
- Spirtes, Peter, Glymour, Clark and Scheines, Richard (2000). *Causation, Prediction, and Search*. 2nd edition. Cambridge, MA: MIT Press.
- Westreich, Daniel and Greenland, Sander (2013). The Table 2 Fallacy: Presenting and Interpreting Confounder and Modifier Coefficients. *American Journal of Epidemiology* 177(4), 292–298.
- Wright, Sewall (1921). Correlation and Causation. *Journal of Agricultural Research* 20(7), 557–585.
- Wright, Sewall (1934). The Method of Path Coefficients. *The Annals of Mathematical Statistics* 5(3), 161–215.
