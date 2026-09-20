# Summary and Notes

## Summary

::: {.idea}

1. An exponential dispersion family has density
   \( \exp[w\{y\theta-b(\theta)\}/\phi+c(y,\phi,w)] \) with natural parameter \( \theta \),
   cumulant function \( b \), dispersion \( \phi \) and known weight \( w \)
   (@def-glm-edf). The natural parameter space is convex and \( b \) is strictly convex on it
   (@lem-glm-convex).

2. All moments come from \( b \): \( \mu=b'(\theta) \), \( \Var(Y)=\phi b''(\theta)/w \), and the
   \( r \)th cumulant is \( (\phi/w)^{r-1}b^{(r)}(\theta) \) (@prp-glm-moments). Reparameterizing by
   the mean gives the variance function \( V(\mu) \) (@def-glm-variance), which determines the
   family: \( V=1,\ \mu,\ \mu^2,\ \mu^3,\ \mu(1-\mu),\ \mu(1+\mu/k) \).

3. A generalized linear model is a random component from such a family, a linear predictor
   \( \X\bbeta \), and a link \( g(\mu_i)=\eta_i \) (@def-glm-model); the canonical link
   \( g=(b')^{-1} \) makes \( \theta_i=\eta_i \) (@def-glm-link). The link and the variance
   function are separate decisions, and neither transforms the data.

4. The score is \( \phi^{-1}\X\T\W\bD^{-1}(\y-\bmu) \) with working weights
   \( \W_{ii}=w_ih'(\eta_i)^2/V(\mu_i) \); the expected information is \( \phi^{-1}\X\T\W\X \), and
   under the canonical link the observed information equals it (@thm-glm-score). The
   likelihood equations do not involve \( \phi \).

5. Under the canonical link the fit reproduces the weighted column totals of \( \X \) exactly
   (@cor-glm-marginals), and the log-likelihood is strictly concave, so the maximum is unique
   when it exists (@thm-glm-concave). Existence can fail, as in separation.

6. Fisher scoring is weighted least squares on the working response
   \( z_i=\eta_i+(y_i-\mu_i)g'(\mu_i) \) with weights \( \W \), iterated (@thm-glm-irls). It is
   Newton's method exactly when the link is canonical, which is why canonical-link fits
   converge quadratically and others linearly.

7. Under stated regularity conditions — notably \( \max_iq_{ni}\to0 \), the leverage-like
   condition @eq-glm-leverage — the estimate exists with probability tending to one, is
   consistent, and \( \A_n^{1/2}(\hbeta-\bbeta^0) \) is asymptotically standard normal
   (@thm-glm-asymptotics for the canonical link, @prp-glm-noncanonical for the rest). Wald,
   score and likelihood ratio statistics are asymptotically equivalent and \( \chi^2(q) \)
   (@prp-glm-three-tests), but the Wald statistic is not invariant and can fall as the
   evidence grows (the Hauck–Donner effect).

8. The deviance is twice the gap to the saturated log-likelihood
   (@def-glm-deviance); every contribution is nonnegative
   (@prp-glm-deviance-positive). Differences of deviances between nested models are
   \( \chi^2 \) in the limit; the deviance itself is a goodness-of-fit statistic only for
   grouped data with large groups, and for ungrouped binary data it is a function of
   \( \hbeta \) alone (@thm-glm-deviance).

9. The dispersion is estimated by \( X^2/(n-p) \) (@eq-glm-phi-hat), which leaves \( \hbeta \)
   unchanged. Pearson, deviance and Anscombe residuals (@def-glm-residuals) and the weighted
   hat matrix carry the diagnostics of
   [Chapter 20](../ch20-residuals-leverage-influence/index.html) across, and AIC applies as
   in [Chapter 29](../ch29-model-selection/index.html).

:::

## Notes and sources

**The synthesis.** Nelder and Wedderburn (1972) created the subject. Every piece existed
before — probit analysis, logistic regression, log-linear models for contingency tables, the
gamma model for survival times, weighted least squares — and the paper's contribution was to
see that the exponential family structure made all of them one algorithm. That is what let a
single program, GLIM, fit them all, and what makes the code
of [Section 34.4](04-irls.html) family-agnostic. McCullagh and Nelder (1989) remains the
definitive treatment; its chapter 2 is the source for the deviance and residual definitions
used here, chapter 12 for the diagnostics.

**Exponential dispersion families.** The class of @def-glm-edf, with \( \phi \) a genuine
second parameter rather than a nuisance folded into \( \theta \), is due to
Jørgensen (1987), developed in Jørgensen (1997); that paper also classifies the power
variance functions \( V(\mu)=\mu^{\,\xi} \), the Tweedie families of @exr-glm-tweedie. The
analytic facts behind @prp-glm-moments are in Widder (1941); the statistical apparatus built
on them, including the convexity of the natural parameter space and the mean-value mapping,
is Barndorff-Nielsen (1978) and Brown (1986); the existence criterion quoted
after @exm-glm-no-maximum is Barndorff-Nielsen's.

**Fitting.** The scoring method is older than the subject: it appears in Fisher's appendix
to Bliss (1935) on probit analysis, as a device for fitting a dose–response curve by hand.
Nelder and Wedderburn (1972) identified the general form and recognized it as weighted least
squares. Green (1984) is the standard reference for iteratively reweighted least squares as
a general estimation device and for the convergence issues
of [Section 34.4](04-irls.html). Existence and uniqueness link by link are settled by
Wedderburn (1976), the source of the gamma/identity warning
in [Section 34.3](03-likelihood-equations.html); that only the variance function enters the
estimating equations, so one may drop the density, is Wedderburn (1974), which Chapter 38
follows.

**Asymptotics.** The complete theory — consistency, asymptotic normality, and existence
with probability tending to one — is Fahrmeir and Kaufmann (1985), under conditions
considerably weaker than (G1)–(G4) and covering non-canonical links and stochastic
regressors; @prp-glm-noncanonical quotes what this book needs from it. The proof
of @thm-glm-asymptotics takes the shortcut the canonical link allows: the log-likelihood is
concave and its Hessian nonrandom, so the argument reduces to the central limit theorem
of @lem-het-clt for the score, plus the convexity lemma of Pollard (1991), which makes
concavity do the work usually done by uniform convergence arguments. The condition
\( \max_iq_{ni}\to0 \) is the generalized-linear-model form of the bounded-leverage condition
used throughout the book, and \( \sum_iq_{ni}=p \) is @prp-proj-leverage(b) in a new setting.

**The three tests.** Cox and Hinkley (1974, chapter 9) is the standard general account of
their asymptotic equivalence. The failure of @exm-glm-hauck-donner was identified by Hauck
and Donner (1977) and named after them; it is one case of a general phenomenon, since the
Wald statistic uses a quadratic approximation to the log-likelihood centred at \( \hbeta \),
worst exactly where the estimate is most extreme. None of this arises in the normal linear
model, where the log-likelihood is exactly quadratic in \( \bbeta \) for fixed \( \sigma^2 \),
and @prp-glh-trinity records the exact relationships that result.

**The deviance and its limits.** That the residual deviance of a binary logistic fit is a
function of \( \hbeta \) alone, hence useless as a goodness-of-fit statistic, is in McCullagh
and Nelder (1989, section 4.4). The regime in which the deviance *is* a goodness-of-fit
statistic — groups fixed, group sizes growing — is sometimes called grouped or
small-dispersion asymptotics, against the \( n\to\infty \) asymptotics
of @thm-glm-asymptotics. The two give different answers for the same statistic, the source
of most of the confusion in practice, and
[Chapter 35](../ch35-binary-responses/index.html) returns to it with @prp-bin-fit. Anscombe
residuals are from Anscombe (1953); the extension of the diagnostics of
[Chapter 20](../ch20-residuals-leverage-influence/index.html), including the weighted hat
matrix and one-step deletion approximations, is Pregibon (1981).

**Textbook treatments.** Agresti (2015, chapter 4) covers the same ground at a similar
level and is the closest companion; its chapters 5–7 correspond to
[Chapter 35](../ch35-binary-responses/index.html),
[Chapter 36](../ch36-multinomial-ordinal/index.html) and
[Chapter 37](../ch37-counts/index.html). Fahrmeir, Kneib, Lang and Marx (2021, section 5.4)
handles the weight \( w \) and grouped data carefully, and is the source of the convention
used here that \( y \) may be a group average. Dobson and Barnett (2018) is a gentler
introduction; Rencher and Schaalje (2008, chapter 18) treats the subject briefly, as an
extension of linear-model theory.

**What comes next.** Chapter 38 keeps the score equations and drops the density, which
handles data more variable than the family allows. Chapter 39 takes up model building in
practice and the Bayesian treatment, whose priors are the counterpart of the penalized
likelihoods of [Chapter 27](../ch27-shrinkage/index.html) and
[Chapter 30](../ch30-regularization-boosting/index.html). Chapter 40 joins the subject to the
mixed models of [Chapter 32](../ch32-linear-mixed-models/index.html), where the likelihood is
no longer available in closed form.

## References

- Agresti, Alan (2015). *Foundations of Linear and Generalized Linear Models*. Hoboken, NJ: Wiley.
- Anscombe, Francis J. (1953). Contribution to the Discussion of H. Hotelling's Paper. *Journal of the Royal Statistical Society, Series B* 15(2), 229–230.
- Barndorff-Nielsen, Ole (1978). *Information and Exponential Families in Statistical Theory*. Chichester: Wiley.
- Bliss, Chester I. (1935). The Calculation of the Dosage–Mortality Curve. *Annals of Applied Biology* 22(1), 134–167.
- Brown, Lawrence D. (1986). *Fundamentals of Statistical Exponential Families, with Applications in Statistical Decision Theory*. Hayward, CA: Institute of Mathematical Statistics.
- Cox, David R. and Hinkley, David V. (1974). *Theoretical Statistics*. London: Chapman and Hall.
- Dobson, Annette J. and Barnett, Adrian G. (2018). *An Introduction to Generalized Linear Models*. 4th edition. Boca Raton, FL: Chapman and Hall/CRC.
- Fahrmeir, Ludwig and Kaufmann, Heinz (1985). Consistency and Asymptotic Normality of the Maximum Likelihood Estimator in Generalized Linear Models. *The Annals of Statistics* 13(1), 342–368.
- Fahrmeir, Ludwig, Kneib, Thomas, Lang, Stefan and Marx, Brian D. (2021). *Regression: Models, Methods and Applications*. 2nd edition. Berlin: Springer.
- Green, Peter J. (1984). Iteratively Reweighted Least Squares for Maximum Likelihood Estimation, and Some Robust and Resistant Alternatives. *Journal of the Royal Statistical Society, Series B* 46(2), 149–192.
- Hauck, Walter W. and Donner, Allan (1977). Wald's Test as Applied to Hypotheses in Logit Analysis. *Journal of the American Statistical Association* 72(360), 851–853.
- Jørgensen, Bent (1987). Exponential Dispersion Models. *Journal of the Royal Statistical Society, Series B* 49(2), 127–162.
- Jørgensen, Bent (1997). *The Theory of Dispersion Models*. London: Chapman and Hall.
- McCullagh, Peter and Nelder, John A. (1989). *Generalized Linear Models*. 2nd edition. London: Chapman and Hall.
- Nelder, John A. and Wedderburn, Robert W. M. (1972). Generalized Linear Models. *Journal of the Royal Statistical Society, Series A* 135(3), 370–384.
- Pollard, David (1991). Asymptotics for Least Absolute Deviation Regression Estimators. *Econometric Theory* 7(2), 186–199.
- Pregibon, Daryl (1981). Logistic Regression Diagnostics. *The Annals of Statistics* 9(4), 705–724.
- Rencher, Alvin C. and Schaalje, G. Bruce (2008). *Linear Models in Statistics*. 2nd edition. Hoboken, NJ: Wiley.
- Rockafellar, R. Tyrrell (1970). *Convex Analysis*. Princeton, NJ: Princeton University Press.
- Wedderburn, Robert W. M. (1974). Quasi-Likelihood Functions, Generalized Linear Models, and the Gauss–Newton Method. *Biometrika* 61(3), 439–447.
- Wedderburn, Robert W. M. (1976). On the Existence and Uniqueness of the Maximum Likelihood Estimates for Certain Generalized Linear Models. *Biometrika* 63(1), 27–32.
- Widder, David V. (1941). *The Laplace Transform*. Princeton, NJ: Princeton University Press.
