# Summary and Notes

## Summary

::: {.idea}

1. If \( \Z\sim\Normal_r(\bmu,\I) \) then \( \norm{\Z}^2\sim\chi^2(r,\gamma) \) with
           \( \gamma=\norm{\bmu}^2 \), mean \( r+\gamma \) and variance \( 2r+4\gamma \). The law is a
           Poisson\( (\gamma/2) \) mixture of \( \chi^2(r+2k) \) laws, and its tails increase strictly
           with \( \gamma \) (@thm-qf-ncchisq and @prp-qf-ncchisq-monotone).

2. \( F(r,s,\gamma) \) and \( t(s,\delta) \) put a noncentral numerator over an independent
           central denominator. The power of an \( F \) test increases strictly with \( \gamma \), and
           \( T^2\sim F(1,s,\delta^2) \) (@thm-qf-f-power).

3. \( \E(\Y\T\A\Y)=\tr(\A\bSigma)+\bmu\T\A\bmu \) for any distribution. Under normality
           \[
\Var(\Y\T\A\Y)=2\tr\{(\A\bSigma)^2\}+4\bmu\T\A\bSigma\A\bmu,\qquad
          \Cov(\mathbf{K}\Y,\Y\T\A\Y)=2\mathbf{K}\bSigma\A\bmu
\]
           (@thm-qf-mean-var).

4. A normal quadratic form is a constant plus independent terms
           \( \lambda_j(W_j+\nu_j)^2 \), and for singular \( \bSigma \) possibly normal terms. The
           weights \( \lambda_j \) are the nonzero eigenvalues of \( \A\bSigma \) (@thm-qf-canonical).

5. For \( \Y\sim\Normal_n(\bmu,\I) \), \( \Y\T\A\Y \) is chi-squared iff \( \A \) is idempotent. Then
           the degrees of freedom are \( \rank\A=\tr\A \) and the noncentrality is \( \bmu\T\A\bmu \)
           (@thm-qf-chisq-identity). For positive definite \( \bSigma \) the condition is that
           \( \A\bSigma \) be idempotent. For projections,
           \( \norm{\bP\Y}^2/\sigma^2\sim\chi^2(\rank\bP,\norm{\bP\bmu}^2/\sigma^2) \) (@thm-qf-chisq).

6. \( \B\Y \) and \( \Y\T\A\Y \) are independent if \( \B\bSigma\A=\mathbf{0} \), and \( \Y\T\A\Y \) and
           \( \Y\T\B\Y \) are independent if \( \A\bSigma\B=\mathbf{0} \), because a quadratic form is a function
           of \( \A\Y \) (@thm-qf-indep-linear and @thm-qf-indep-quadratic).

7. Projections of a spherical normal vector onto mutually orthogonal subspaces are
           independent. This gives the independence of \( \bar{Y} \) and \( S^2 \) and the exact
           \( F(r-r_0,n-r,\gamma) \) law of the nested-model \( F \) statistic
           (@thm-qf-orthogonal-projections, @cor-qf-sample-variance and @thm-qf-nested-f).

8. Cochran: if symmetric matrices add to an idempotent matrix and their ranks add up,
           the corresponding sums of squares are independent chi-squared variables
           (@thm-qf-cochran-algebra and @thm-qf-cochran).

:::

## Notes and sources

**Source texts.**  Rencher and Schaalje (2008, chapter 5) cover the material of
[Section 4.2](02-moments.html), [Section 4.3](03-chisq.html), [Section 4.4](04-independence.html) and [Section 4.5](05-cochran.html) for positive definite
covariance matrices. They derive the mgf of a quadratic form by integrating against the
normal density and obtain its variance by differentiating the logarithm of the mgf; the converse
of the chi-squared theorem is referred to the literature. They use the noncentrality
\( \bmu\T\bmu/2 \). Seber and Lee (2003, section 2.4)
reduce the general case to a spherical normal vector by an orthogonal diagonalization. They
prove necessity of idempotence for central forms by factorizing the moment generating
function, and give the eigenvalue criterion for singular covariance matrices.
Christensen (2020, section 1.3) states the three conditions of
@thm-qf-chisq(c) for arbitrary nonnegative definite \( \bSigma \), and obtains independence
by viewing quadratic forms as functions of linear ones, the approach followed in
[Section 4.4](04-independence.html). The canonical representation (@thm-qf-canonical), the proof of
necessity through the coefficients @eq-qf-cumulant-match, the monotonicity proof through
the Poisson mixture, and the treatment of the Cochran algebra are this book's own
arrangements of standard material.

**Noncentral distributions.**  The noncentral chi-squared distribution first appeared in
Fisher (1928). Patnaik (1949) studied its approximation by a scaled
central chi-squared variable, and Satterthwaite (1946) used the same two-moment
idea for linear combinations of mean squares, which reappears in the mixed models of
[Chapter 32](../ch32-linear-mixed-models/index.html) (@prp-mix-inference).
Johnson et al. (1995, chapters 29–31) collect the properties of the noncentral
chi-squared, \( F \) and \( t \) families. The monotonicity of power in the noncentrality and in the
degrees of freedom is treated systematically by Ghosh (1973).

**General quadratic forms.**  The exact distribution of a quadratic form that is not
chi-squared can be computed by numerical inversion of its characteristic function
(Gil-Pelaez 1951; Imhof 1961). Mathai and Provost (1992) give a
book-length treatment of quadratic forms in normal and nonnormal variables. The mean square
successive difference of @exm-qf-drift goes back to von Neumann (1941),
who derived the exact distribution of its ratio to \( S^2 \).

**Independence and Cochran's theorem.**  Cochran (1934) proved the
decomposition theorem that bears his name, in the context of the analysis of covariance.
The necessity of \( \A\bSigma\B=\mathbf{0} \) for independence of two forms is attributed to
Craig (1943). Its full proof for arbitrary symmetric matrices and arbitrary means
turned out to be subtle, and Driscoll and Gundberg (1986) trace the flawed and correct
arguments. Searle (1971) gives general proofs of the independence and Cochran
theorems in the style most linear-models texts follow. For nonnegative definite \( \A \) and \( \B \),
the condition \( \bSigma\A\bSigma\B\bSigma=\mathbf{0} \) of Christensen (2020, Theorem 1.3.8)
is in fact equivalent to \( \A\bSigma\B=\mathbf{0} \) (@exr-qf-singular-independence), so that
theorem adds nothing to @thm-qf-indep-quadratic even for singular \( \bSigma \), whereas his
Theorem 1.3.9, for general symmetric \( \A \) and \( \B \), does. That the independence of \( \bar{Y} \) and
\( S^2 \) characterizes the normal distribution was shown by Lukacs (1942).

## References

- Christensen, Ronald (2020). *Plane Answers to Complex Questions: The Theory of Linear Models*. 5th edition. Cham: Springer.
- Cochran, William G. (1934). The Distribution of Quadratic Forms in a Normal System, with Applications to the Analysis of Covariance. *Mathematical Proceedings of the Cambridge Philosophical Society* 30(2), 178–191.
- Craig, Allen T. (1943). Note on the Independence of Certain Quadratic Forms. *The Annals of Mathematical Statistics* 14(2), 195–197.
- Driscoll, Michael F. and Gundberg, William R. (1986). A History of the Development of Craig's Theorem. *The American Statistician* 40(1), 65–70.
- Fisher, Ronald A. (1928). The General Sampling Distribution of the Multiple Correlation Coefficient. *Proceedings of the Royal Society of London. Series A* 121(788), 654–673.
- Ghosh, Bhaskar K. (1973). Some Monotonicity Theorems for \( \chi^2 \), F and t Distributions with Applications. *Journal of the Royal Statistical Society, Series B* 35(3), 480–492.
- Gil-Pelaez, J. (1951). Note on the Inversion Theorem. *Biometrika* 38(3/4), 481–482.
- Imhof, Jean-Pierre (1961). Computing the Distribution of Quadratic Forms in Normal Variables. *Biometrika* 48(3/4), 419–426.
- Johnson, Norman L., Kotz, Samuel and Balakrishnan, Narayanaswamy (1995). *Continuous Univariate Distributions*. 2nd edition. New York: Wiley.
- Lukacs, Eugene (1942). A Characterization of the Normal Distribution. *The Annals of Mathematical Statistics* 13(1), 91–93.
- Mathai, Arak M. and Provost, Serge B. (1992). *Quadratic Forms in Random Variables: Theory and Applications*. New York: Marcel Dekker.
- Patnaik, P. B. (1949). The Non-Central \( \chi^2 \)- and F-Distributions and Their Applications. *Biometrika* 36(1/2), 202–232.
- Rencher, Alvin C. and Schaalje, G. Bruce (2008). *Linear Models in Statistics*. 2nd edition. Hoboken, NJ: Wiley.
- Satterthwaite, Franklin E. (1946). An Approximate Distribution of Estimates of Variance Components. *Biometrics Bulletin* 2(6), 110–114.
- Searle, Shayle R. (1971). *Linear Models*. New York: Wiley.
- Seber, George A. F. and Lee, Alan J. (2003). *Linear Regression Analysis*. 2nd edition. Hoboken, NJ: Wiley.
- von Neumann, John (1941). Distribution of the Ratio of the Mean Square Successive Difference to the Variance. *The Annals of Mathematical Statistics* 12(4), 367–395.
