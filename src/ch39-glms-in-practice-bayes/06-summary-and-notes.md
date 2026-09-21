# Summary and Notes

## Summary

::: {.idea}

1. A response family is chosen from the support and from the mean–variance relationship;
   with replicates, the variance-function power is estimated by regressing log group
   variance on log group mean (@prp-prc-model-building(a)).

2. Additivity belongs to the link, not to the data: a mean additive on one link scale is
   additive on another only when the two differ by an affine
   map (@prp-prc-model-building(b)).

3. Likelihood criteria compare models for the same random variable. A model fitted to
   \( t(\y) \) needs the correction \( -2\sum_i\log t'(y_i) \), itself of order
   \( n \), before its AIC can stand beside one for
   \( \y \) (@prp-prc-model-building(c), @exm-prc-strikes).

4. Diagnostics should be aimed. Pearson residuals have approximate variance
   \( 1-h_{ii} \); a simulated envelope has pointwise level \( 2/(B+1) \); the link is
   tested by adding \( \hat\eta^2 \); a covariate's scale by partial-residual and
   added-variable plots; the variance function by the height and trend of
   \( |r^{P}| \) (@prp-prc-diagnostics). In @exm-prc-tvnews only the last is wrong.

5. The log-likelihood is quadratic in \( \bbeta \) for every design exactly when the
   family is normal (@prp-prc-quadratic): in the Bayesian thread every definition
   survives, every closed form does not.

6. The Bayesian generalized linear model is @def-prc-posterior. A proper prior with a
   bounded likelihood gives a proper posterior; a flat prior on a logistic model gives one
   exactly when the data are in overlap, that is exactly when the maximum likelihood
   estimate exists (@prp-prc-propriety).

7. Under local smoothness, a growing Hessian and a concentration condition, the Laplace
   approximation to the marginal likelihood has relative error
   \( O(n^{-1}) \) (@thm-prc-laplace); its logarithm is the BIC with an \( O(1) \)
   remainder carrying the whole influence of the prior (@cor-prc-bic).

8. Metropolis–Hastings and Gibbs leave the posterior invariant by detailed
   balance (@prp-prc-mh-invariance, @def-prc-mcmc); Albert–Chib augmentation makes every
   probit full conditional normal or truncated normal (@prp-prc-albert-chib). Convergence
   is judged by \( \hat R \), mixing by effective sample size (@exm-prc-anes-mcmc).

9. Posterior predictive checks compare the data with replicates from @eq-prc-ppd; using
   the data twice makes them conservative, so an extreme \( p \)-value is evidence and a
   middling one is not (@def-prc-ppc, @exm-prc-ppc).

10. A normal prior is ridge, a Laplace prior the lasso at its mode, and a normal prior on
    a logit is bimodal on the probability scale exactly when
    \( \sigma^2>2 \) (@prp-prc-priors). A proper prior makes the mode exist under
    separation, where maximum likelihood has none (@exm-prc-separated-bayes).

:::

## Notes and sources

**Model building.**  The three-part structure — random component, systematic component,
link — is Nelder and Wedderburn's (1972); McCullagh and Nelder (1989) remains the standard
account of how the pieces are chosen, their chapter 8 for the log link being preferred to
the canonical reciprocal in gamma models and their chapter 11 for mean–variance
diagnostics. Fahrmeir, Kneib, Lang and Marx (2021, §5.3) develop the log-normal and gamma
alternatives for nonnegative responses. The Jacobian correction
of @prp-prc-model-building(c) is elementary, widely known, and as widely omitted from
software output.

**Diagnostics.**  Pregibon (1981) transferred the linear-model diagnostics of
[Chapter 20](../ch20-residuals-leverage-influence/index.html) to generalized linear models
by the substitution \( \X\mapsto\hat{\W}^{1/2}\X \), and Pregibon (1980) introduced the
goodness-of-link test of @prp-prc-diagnostics(d). The half-normal plot with a simulated
envelope is Atkinson's (1981, 1985), and the exchangeability argument behind
@eq-prc-envelope is the one that justifies Monte Carlo tests generally. Cook and Weisberg
(1983) read a trend in absolute residuals as evidence about the variance function.

**The Bayesian thread.**  Lindley and Smith (1972) is the origin of the hierarchical
linear model. Diaconis and Ylvisaker (1979) characterized the conjugate priors of an
exponential family; the generalized linear form @eq-prc-dy, with a pseudo-response and
per-observation weights, is Bedrick, Christensen and Johnson (1996) and Chen and Ibrahim
(2003). Rencher and Schaalje (2008, chapter 11) and Agresti (2015, chapter 10) are the
source treatments closest to this chapter, and Gelman, Carlin, Stern, Dunson, Vehtari and
Rubin (2013) the standard modern reference for its last two sections.

**Laplace and the BIC.**  Laplace's method is two centuries old; its use for posterior
moments and marginal likelihoods is Tierney and Kadane's (1986). Schwarz (1978) derived
the BIC, and Kass and Raftery (1995) made the link between Bayes factors, the Laplace
approximation and the BIC standard; their section 4 discusses the regularity conditions
that @thm-prc-laplace assumes, and the set used here — a shrinking ball, bounded scaled
derivatives, a concentration bound — is convenient rather than weakest. The deviance
information criterion is Spiegelhalter, Best, Carlin and van der Linde (2002); the
Bernstein–von Mises theorem in the form quoted is van der Vaart (1998, chapter 10).

**MCMC.**  The algorithm of @eq-prc-mh is Metropolis, Rosenbluth, Rosenbluth, Teller and
Teller (1953) as generalized by Hastings (1970); Tierney (1994) gave the rigorous account
for continuous state spaces and Robert and Casella (2004) the textbook treatment. The
Gibbs sampler is Geman and Geman (1984), brought into mainstream statistics by Gelfand and
Smith (1990). Optimal scaling is Roberts, Gelman and Gilks (1997), and the iteratively
reweighted proposal follows Gamerman (1997). Data augmentation for the probit model is
Albert and Chib (1993); for the logit, Holmes and Held (2006) use a scale mixture and
Polson, Scott and Windle (2013) the Pólya–gamma identity now in general use. @eq-prc-rhat
is Gelman and Rubin (1992) and the truncation rule for the effective sample size is
Geyer's (1992). Posterior predictive checking is Rubin (1984), Meng (1994) and Gelman,
Meng and Stern (1996), its \( p \)-values analysed by Bayarri and Berger (2000).
Hamiltonian Monte Carlo, which is what modern software actually runs, is omitted because
its justification needs machinery this book does not develop.

**Priors as penalties.**  The ridge–normal correspondence is as old as ridge
regression (@thm-shr-ridge); the lasso–Laplace one was noted by Tibshirani (1996) and
turned into a sampler by Park and Casella (2008). Gelman, Jakulin, Pittau and Su (2008)
proposed the Cauchy default prior for logistic regression with the rescaling to standard
deviation \( 1/2 \), and gave the separation argument that @prp-prc-priors(d) formalizes.
That Firth's (1993) penalty is the Jeffreys (1946) prior was observed by Firth himself;
Heinze and Schemper (2002) made it the standard remedy. The horseshoe is Carvalho, Polson
and Scott (2010). The propriety argument of @prp-prc-propriety(b) is the Bayesian face
of @thm-bin-separation.

## References

- Agresti, Alan (2015). *Foundations of Linear and Generalized Linear Models*. Hoboken, NJ: Wiley.
- Albert, James H. and Chib, Siddhartha (1993). Bayesian Analysis of Binary and Polychotomous Response Data. *Journal of the American Statistical Association* 88(422), 669–679.
- Atkinson, Anthony C. (1981). Two Graphical Displays for Outlying and Influential Observations in Regression. *Biometrika* 68(1), 13–20.
- Atkinson, Anthony C. (1985). *Plots, Transformations, and Regression: An Introduction to Graphical Methods of Diagnostic Regression Analysis*. Oxford: Clarendon Press.
- Bedrick, Edward J., Christensen, Ronald and Johnson, Wesley (1996). A New Perspective on Priors for Generalized Linear Models. *Journal of the American Statistical Association* 91(436), 1450–1460.
- Bayarri, M. J. and Berger, James O. (2000). P Values for Composite Null Models. *Journal of the American Statistical Association* 95(452), 1127–1142.
- Carvalho, Carlos M., Polson, Nicholas G. and Scott, James G. (2010). The Horseshoe Estimator for Sparse Signals. *Biometrika* 97(2), 465–480.
- Chen, Ming-Hui and Ibrahim, Joseph G. (2003). Conjugate Priors for Generalized Linear Models. *Statistica Sinica* 13(2), 461–476.
- Cook, R. Dennis and Weisberg, Sanford (1983). Diagnostics for Heteroscedasticity in Regression. *Biometrika* 70(1), 1–10.
- Diaconis, Persi and Ylvisaker, Donald (1979). Conjugate Priors for Exponential Families. *The Annals of Statistics* 7(2), 269–281.
- Fahrmeir, Ludwig, Kneib, Thomas, Lang, Stefan and Marx, Brian D. (2021). *Regression: Models, Methods and Applications*. 2nd edition. Berlin: Springer.
- Firth, David (1993). Bias Reduction of Maximum Likelihood Estimates. *Biometrika* 80(1), 27–38.
- Gamerman, Dani (1997). Sampling from the Posterior Distribution in Generalized Linear Mixed Models. *Statistics and Computing* 7(1), 57–68.
- Gelfand, Alan E. and Smith, Adrian F. M. (1990). Sampling-Based Approaches to Calculating Marginal Densities. *Journal of the American Statistical Association* 85(410), 398–409.
- Gelman, Andrew, Carlin, John B., Stern, Hal S., Dunson, David B., Vehtari, Aki and Rubin, Donald B. (2013). *Bayesian Data Analysis*. 3rd edition. Boca Raton, FL: CRC Press.
- Gelman, Andrew, Jakulin, Aleks, Pittau, Maria Grazia and Su, Yu-Sung (2008). A Weakly Informative Default Prior Distribution for Logistic and Other Regression Models. *The Annals of Applied Statistics* 2(4), 1360–1383.
- Gelman, Andrew, Meng, Xiao-Li and Stern, Hal (1996). Posterior Predictive Assessment of Model Fitness via Realized Discrepancies. *Statistica Sinica* 6(4), 733–760.
- Gelman, Andrew and Rubin, Donald B. (1992). Inference from Iterative Simulation Using Multiple Sequences. *Statistical Science* 7(4), 457–472.
- Geyer, Charles J. (1992). Practical Markov Chain Monte Carlo. *Statistical Science* 7(4), 473–483.
- Geman, Stuart and Geman, Donald (1984). Stochastic Relaxation, Gibbs Distributions, and the Bayesian Restoration of Images. *IEEE Transactions on Pattern Analysis and Machine Intelligence* 6(6), 721–741.
- Hastings, W. Keith (1970). Monte Carlo Sampling Methods Using Markov Chains and Their Applications. *Biometrika* 57(1), 97–109.
- Heinze, Georg and Schemper, Michael (2002). A Solution to the Problem of Separation in Logistic Regression. *Statistics in Medicine* 21(16), 2409–2419.
- Holmes, Chris C. and Held, Leonhard (2006). Bayesian Auxiliary Variable Models for Binary and Multinomial Regression. *Bayesian Analysis* 1(1), 145–168.
- Jeffreys, Harold (1946). An Invariant Form for the Prior Probability in Estimation Problems. *Proceedings of the Royal Society of London. Series A* 186(1007), 453–461.
- Kass, Robert E. and Raftery, Adrian E. (1995). Bayes Factors. *Journal of the American Statistical Association* 90(430), 773–795.
- Lindley, Dennis V. and Smith, Adrian F. M. (1972). Bayes Estimates for the Linear Model. *Journal of the Royal Statistical Society, Series B* 34(1), 1–41.
- McCullagh, Peter and Nelder, John A. (1989). *Generalized Linear Models*. 2nd edition. London: Chapman and Hall.
- Meng, Xiao-Li (1994). Posterior Predictive \( p \)-Values. *The Annals of Statistics* 22(3), 1142–1160.
- Metropolis, Nicholas, Rosenbluth, Arianna W., Rosenbluth, Marshall N., Teller, Augusta H. and Teller, Edward (1953). Equation of State Calculations by Fast Computing Machines. *The Journal of Chemical Physics* 21(6), 1087–1092.
- Nelder, John A. and Wedderburn, Robert W. M. (1972). Generalized Linear Models. *Journal of the Royal Statistical Society, Series A* 135(3), 370–384.
- Park, Trevor and Casella, George (2008). The Bayesian Lasso. *Journal of the American Statistical Association* 103(482), 681–686.
- Polson, Nicholas G., Scott, James G. and Windle, Jesse (2013). Bayesian Inference for Logistic Models Using Pólya–Gamma Latent Variables. *Journal of the American Statistical Association* 108(504), 1339–1349.
- Pregibon, Daryl (1980). Goodness of Link Tests for Generalized Linear Models. *Journal of the Royal Statistical Society, Series C* 29(1), 15–24.
- Pregibon, Daryl (1981). Logistic Regression Diagnostics. *The Annals of Statistics* 9(4), 705–724.
- Rencher, Alvin C. and Schaalje, G. Bruce (2008). *Linear Models in Statistics*. 2nd edition. Hoboken, NJ: Wiley.
- Robert, Christian P. and Casella, George (2004). *Monte Carlo Statistical Methods*. 2nd edition. New York: Springer.
- Roberts, Gareth O., Gelman, Andrew and Gilks, Walter R. (1997). Weak Convergence and Optimal Scaling of Random Walk Metropolis Algorithms. *The Annals of Applied Probability* 7(1), 110–120.
- Rubin, Donald B. (1984). Bayesianly Justifiable and Relevant Frequency Calculations for the Applied Statistician. *The Annals of Statistics* 12(4), 1151–1172.
- Schwarz, Gideon (1978). Estimating the Dimension of a Model. *The Annals of Statistics* 6(2), 461–464.
- Spiegelhalter, David J., Best, Nicola G., Carlin, Bradley P. and van der Linde, Angelika (2002). Bayesian Measures of Model Complexity and Fit. *Journal of the Royal Statistical Society, Series B* 64(4), 583–639.
- Tibshirani, Robert (1996). Regression Shrinkage and Selection via the Lasso. *Journal of the Royal Statistical Society, Series B* 58(1), 267–288.
- Tierney, Luke (1994). Markov Chains for Exploring Posterior Distributions. *The Annals of Statistics* 22(4), 1701–1728.
- Tierney, Luke and Kadane, Joseph B. (1986). Accurate Approximations for Posterior Moments and Marginal Densities. *Journal of the American Statistical Association* 81(393), 82–86.
- van der Vaart, Aad W. (1998). *Asymptotic Statistics*. Cambridge: Cambridge University Press.
