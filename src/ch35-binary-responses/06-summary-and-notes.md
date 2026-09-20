# Summary and Notes

## Summary

::: {.idea}

1. A binary response is a binomial exponential dispersion family with variance function
   \( \pi(1-\pi) \), dispersion one and natural parameter the log odds; the canonical
   link is the logit (@def-bin-logistic).

2. A logistic coefficient is a log odds ratio, constant over the other covariates;
   \( |\beta_j|/4 \) bounds its effect on the probability scale, and the average marginal
   effect converts it to a change in probability (@prp-bin-interpretation).

3. Grouping observations with identical covariates changes the log-likelihood by a
   constant, so the estimate, its covariance and every likelihood ratio are unchanged;
   the deviance is not (@prp-bin-grouping).

4. Fisher scoring is iteratively reweighted least squares with weights
   \( m_i\pi_i(1-\pi_i) \); observed and expected information coincide, and the
   log-likelihood is strictly concave when \( \rank(\X)=p \) (@prp-bin-irls).

5. Every binary regression is a threshold model: logistic, normal and smallest
   extreme-value latent errors give the logit, probit and complementary log–log links,
   and a proportional hazards process observed at a fixed time gives the last of them
   and no other (@def-bin-latent and @prp-bin-links).

6. The latent scale is not identified, so probit and logit coefficients are comparable
   only after rescaling, by \( 1.8138 \) or \( 1.5958 \)
   depending on the convention (@prp-bin-latent-scale).

7. The maximum likelihood estimate exists if and only if the data are in overlap; under
   separation some coefficients diverge and the Wald test collapses while the likelihood
   ratio does not (@def-bin-separation and @thm-bin-separation).

8. Firth's penalty \( \tfrac12\log\det(\X\T\W\X) \) removes the leading bias and gives
   finite estimates; a ridge or Cauchy penalty does the same more crudely
   (@prp-bin-firth).

9. Grouped Pearson and deviance statistics have \( \chi^2(N-p) \) limits as the group
   sizes grow. The ungrouped deviance is exactly
   \( -2\sum_i\{\hat\pi_i\log\hat\pi_i+(1-\hat\pi_i)\log(1-\hat\pi_i)\} \), a function of
   the fit alone with no \( \chi^2 \) limit; differences of deviances remain valid
   (@prp-bin-fit).

10. In-sample calibration is automatic, so calibration must be judged out of sample, and
    the area under the ROC curve is a rank statistic, blind to calibration
    (@prp-bin-calibration and @prp-bin-auc).

11. Retrospective sampling on \( Y \) alone shifts the intercept by
    \( \log(\rho_1/\rho_0) \) and leaves the slopes alone, and the logit is the only link
    with that property (@thm-bin-casecontrol).

12. Matched sets are handled by conditioning; fitting the stratum intercepts instead is
    inconsistent, and for matched pairs gives exactly twice the conditional estimate
    (@prp-bin-conditional).

:::

## Notes and sources

**Coverage.** The topics follow Agresti (2015, chapter 5) and Fahrmeir, Kneib, Lang
and Marx (2021, sections 5.1–5.2), with proofs, examples and exercises of this book's
own. Collett (2003) and Hosmer, Lemeshow and Sturdivant (2013) are the standard
applied accounts, and the source of the events-per-variable guidance alluded to in
[Section 35.1](01-logistic-regression.html); McCullagh and Nelder (1989, chapter 4) is
the compact theoretical reference.

**Origins.** The logistic function entered statistics through bioassay: probit
analysis came first, with Bliss (1934) and Finney (1971), and Berkson (1944) argued
for the logistic alternative against considerable resistance. Cox (1958a) gave the
modern regression treatment of binary sequences and Cox (1970) the first book on
binary data. The complementary log–log link is the binary face of the proportional
hazards model of Cox (1972), with the discrete-time version of
@exr-bin-cloglog-interval due to Prentice and Gloeckler (1978); the link family of
@exr-bin-link-family is Aranda-Ordaz's (1981). The "divide by four" rule is from Gelman
and Hill (2007).

**Separation.** Silvapulle (1981) proved existence of the estimate under a convexity
condition equivalent to overlap; Albert and Anderson (1984) introduced the terminology
of complete and quasi-complete separation and the linear-programming diagnostic, and
Santner and Duffy (1986) corrected a gap in their argument. The proof of
@thm-bin-separation given here is the direct one: strict concavity leaves only
coercivity to establish, and coercivity is exactly overlap. The counting formula of
@exr-bin-separation-high-dimensional is Cover's (1965). Firth (1993) introduced the
Jeffreys-prior penalty as a bias-reduction device; Heinze and Schemper (2002)
recognized it as a remedy for separation, and Kosmidis and Firth (2021) proved the
resulting estimate always finite. The Cauchy prior is from Gelman, Jakulin, Pittau and
Su (2008); exact conditional inference is due to Cox (1970) and Hirji, Mehta and Patel
(1987); and the Wald pathology of @exr-bin-hauck-donner is Hauck and Donner's (1977).

**Goodness of fit.** That the ungrouped deviance depends on the data only through the
fitted values is standard and appears in McCullagh and Nelder (1989, section 4.4); the
exact identity @eq-bin-entropy-deviance and the simulation of @exm-bin-deviance-sim
are this book's way of making the point. The Hosmer–Lemeshow statistic is from Hosmer
and Lemeshow (1980), its instability discussed at length in Hosmer, Lemeshow and
Sturdivant (2013); a smoothed alternative is le Cessie and van Houwelingen (1991).
Calibration in the sense of @prp-bin-calibration goes back to Cox (1958b). The ROC
curve entered medical statistics through Hanley and McNeil (1982), and Pepe, Janes,
Longton, Leisenring and Newcomb (2004) documented how insensitive its area is to
genuinely informative covariates.

**Case–control studies.** The retrospective likelihood argument that justifies using
the prospective fit is Prentice and Pyke (1979); Breslow and Day (1980) is the
classical monograph and Breslow (1996) a readable retrospective on the design. The
conditional likelihood of @prp-bin-conditional comes from the tradition of Andersen
(1970), itself an instance of the incidental-parameter problem of Neyman and Scott
(1948) met in @exm-opt-neyman-scott. Noncollapsibility is analysed by Greenland,
Robins and Pearl (1999), who separate it carefully from confounding.

**What comes next.** Chapter 38 takes up quasi-likelihood and overdispersion, which
matter for grouped binary data whose variance exceeds \( m\pi(1-\pi) \); Chapter 39
treats fitting in practice and the Bayesian thread, including the priors touched on in
[Section 35.3](03-separation.html); Chapter 40 joins generalized linear
models to the random effects of
[Chapter 32](../ch32-linear-mixed-models/index.html); and
Chapter 41 covers missing data. Responses with more than two categories are the subject
of [Chapter 36](../ch36-multinomial-ordinal/index.html). The 1996 American National
Election Studies survey used throughout is public domain; it illustrates method only,
and nothing here should be read as a claim about the election.

## References

- Agresti, Alan (2015). *Foundations of Linear and Generalized Linear Models*. Hoboken, NJ: Wiley.
- Albert, Adelin and Anderson, John A. (1984). On the Existence of Maximum Likelihood Estimates in Logistic Regression Models. *Biometrika* 71(1), 1–10.
- Andersen, Erling B. (1970). Asymptotic Properties of Conditional Maximum-Likelihood Estimators. *Journal of the Royal Statistical Society, Series B* 32(2), 283–301.
- Aranda-Ordaz, Francisco J. (1981). On Two Families of Transformations to Additivity for Binary Response Data. *Biometrika* 68(2), 357–363.
- Berkson, Joseph (1944). Application of the Logistic Function to Bio-Assay. *Journal of the American Statistical Association* 39(227), 357–365.
- Bliss, Chester I. (1934). The Method of Probits. *Science* 79(2037), 38–39.
- Breslow, Norman E. (1996). Statistics in Epidemiology: The Case-Control Study. *Journal of the American Statistical Association* 91(433), 14–28.
- Breslow, Norman E. and Day, Nicholas E. (1980). *Statistical Methods in Cancer Research, Volume I: The Analysis of Case-Control Studies*. Lyon: International Agency for Research on Cancer.
- Collett, David (2003). *Modelling Binary Data*. 2nd edition. Boca Raton: Chapman and Hall/CRC.
- Cover, Thomas M. (1965). Geometrical and Statistical Properties of Systems of Linear Inequalities with Applications in Pattern Recognition. *IEEE Transactions on Electronic Computers* EC-14(3), 326–334.
- Cox, David R. (1958a). The Regression Analysis of Binary Sequences. *Journal of the Royal Statistical Society, Series B* 20(2), 215–242.
- Cox, David R. (1958b). Two Further Applications of a Model for Binary Regression. *Biometrika* 45(3/4), 562–565.
- Cox, David R. (1970). *The Analysis of Binary Data*. London: Methuen.
- Cox, David R. (1972). Regression Models and Life-Tables. *Journal of the Royal Statistical Society, Series B* 34(2), 187–220.
- Fahrmeir, Ludwig, Kneib, Thomas, Lang, Stefan and Marx, Brian D. (2021). *Regression: Models, Methods and Applications*. 2nd edition. Berlin: Springer.
- Finney, David J. (1971). *Probit Analysis*. 3rd edition. Cambridge: Cambridge University Press.
- Firth, David (1993). Bias Reduction of Maximum Likelihood Estimates. *Biometrika* 80(1), 27–38.
- Gelman, Andrew and Hill, Jennifer (2007). *Data Analysis Using Regression and Multilevel/Hierarchical Models*. Cambridge: Cambridge University Press.
- Gelman, Andrew, Jakulin, Aleks, Pittau, Maria Grazia and Su, Yu-Sung (2008). A Weakly Informative Default Prior Distribution for Logistic and Other Regression Models. *The Annals of Applied Statistics* 2(4), 1360–1383.
- Greenland, Sander, Robins, James M. and Pearl, Judea (1999). Confounding and Collapsibility in Causal Inference. *Statistical Science* 14(1), 29–46.
- Hanley, James A. and McNeil, Barbara J. (1982). The Meaning and Use of the Area under a Receiver Operating Characteristic (ROC) Curve. *Radiology* 143(1), 29–36.
- Hauck, Walter W. and Donner, Allan (1977). Wald's Test as Applied to Hypotheses in Logit Analysis. *Journal of the American Statistical Association* 72(360), 851–853.
- Heinze, Georg and Schemper, Michael (2002). A Solution to the Problem of Separation in Logistic Regression. *Statistics in Medicine* 21(16), 2409–2419.
- Hirji, Karim F., Mehta, Cyrus R. and Patel, Nitin R. (1987). Computing Distributions for Exact Logistic Regression. *Journal of the American Statistical Association* 82(400), 1110–1117.
- Hosmer, David W. and Lemeshow, Stanley (1980). Goodness of Fit Tests for the Multiple Logistic Regression Model. *Communications in Statistics — Theory and Methods* 9(10), 1043–1069.
- Hosmer, David W., Lemeshow, Stanley and Sturdivant, Rodney X. (2013). *Applied Logistic Regression*. 3rd edition. Hoboken, NJ: Wiley.
- Kosmidis, Ioannis and Firth, David (2021). Jeffreys-Prior Penalty, Finiteness and Shrinkage in Binomial-Response Generalized Linear Models. *Biometrika* 108(1), 71–82.
- le Cessie, Saskia and van Houwelingen, Johannes C. (1991). A Goodness-of-Fit Test for Binary Regression Models Based on Smoothing Methods. *Biometrics* 47(4), 1267–1282.
- McCullagh, Peter and Nelder, John A. (1989). *Generalized Linear Models*. 2nd edition. London: Chapman and Hall.
- Neyman, Jerzy and Scott, Elizabeth L. (1948). Consistent Estimates Based on Partially Consistent Observations. *Econometrica* 16(1), 1–32.
- Pepe, Margaret S., Janes, Holly, Longton, Gary, Leisenring, Wendy and Newcomb, Polly (2004). Limitations of the Odds Ratio in Gauging the Performance of a Diagnostic, Prognostic, or Screening Marker. *American Journal of Epidemiology* 159(9), 882–890.
- Prentice, Ross L. and Gloeckler, Lynn A. (1978). Regression Analysis of Grouped Survival Data with Application to Breast Cancer Data. *Biometrics* 34(1), 57–67.
- Prentice, Ross L. and Pyke, Ronald (1979). Logistic Disease Incidence Models and Case-Control Studies. *Biometrika* 66(3), 403–411.
- Santner, Thomas J. and Duffy, Diane E. (1986). A Note on A. Albert and J. A. Anderson's Conditions for the Existence of Maximum Likelihood Estimates in Logistic Regression Models. *Biometrika* 73(3), 755–758.
- Silvapulle, Mervyn J. (1981). On the Existence of Maximum Likelihood Estimators for the Binomial Response Models. *Journal of the Royal Statistical Society, Series B* 43(3), 310–313.
