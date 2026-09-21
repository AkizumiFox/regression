# Summary and Notes

## Summary

::: {.idea}

1. A polynomial regression is a linear model in a transformed
   basis (@def-ply-polynomial): estimation, the \( F \) tests, the bands, lack of
   fit and the diagnostics of Parts II and III apply unchanged, and only the
   reading of an individual coefficient is lost.

2. Raising the degree buys bias and sells variance unevenly: the standard error
   of the fit grows fastest at the edges, leverage concentrates there, and a
   single observation moves the whole curve (@prp-ply-instability). The monomial
   basis is also one of the worst conditioned in numerical analysis, its Gram
   matrix converging to a Hilbert matrix; centring and rescaling to
   \( [-1,1] \) recover most of that for free.

3. The orthogonal polynomials of a design exist, are unique, and satisfy a
   three-term recurrence, because multiplication by the argument is
   self-adjoint and raises degree by one (@def-ply-orthogonal,
   @lem-ply-shift and @thm-ply-recurrence). In the orthonormal basis
   \( \X\T\X=\I \): coefficients do not change with the degree, sequential sums of
   squares are order-free, and \( \kappa(\X)=1 \) (@cor-ply-sequential).

4. Changing the basis fixes the arithmetic and nothing else: the oscillation,
   the boundary variance and the extrapolation belong to \( \mathcal P_d \), not
   to any list of functions spanning it.

5. Cutting the range at \( K \) knots and fitting degree-\( d \) pieces with
   \( C^{d-1} \) joins gives a space of dimension \( d+1+K \) with the truncated
   power functions as a basis (@def-ply-piecewise, @thm-ply-spline-space). The
   continuity constraints are independent linear restrictions, so each level of
   smoothness is a testable hypothesis.

6. A regression spline is ordinary least squares in that
   basis (@def-ply-regression-spline). Knots buy flexibility one degree of
   freedom at a time and spend it locally; their placement controls where the
   variance goes, and quantile knots are the defensible
   default (@prp-ply-knots).

7. The truncated power basis is for proving things, not for computing with: its
   condition number grows geometrically with \( K \), while a locally supported
   basis for the same space, on a design that occupies every support, does not.

8. A fitted curve is still a conditional mean, still says nothing about
   intervention, and still must not be extrapolated: of the bases compared on
   the carbon dioxide series, the one that fitted best inside the data was among
   the worst outside it ([Section 42.4](04-regression-splines.html)).

:::

## Notes and sources

**Polynomial regression and ill-conditioning.**  That the Gram matrix of the
monomials is essentially a Hilbert matrix is an old observation; Forsythe (1957)
is the paper that drew from it the recommendation to fit in an orthogonal basis,
and Seber and Lee (2003, section 7.1.1) cite him for the formula. The
conditioning facts are from numerical analysis: Todd (1954) for the Hilbert
matrix, and Higham (2002), whose chapters on Vandermonde systems and on test
matrices cover the two cases. Runge (1901) gave the interpolation example and
Trefethen (2019) is the modern account, including Chebyshev nodes and the
Lebesgue constants of @exr-ply-lebesgue; a least squares polynomial of moderate
degree on a large sample is not an interpolant and misbehaves more mildly, but
by the same mechanism.

**Orthogonal polynomials.**  The general theory is classical, and Szegő (1939)
remains the reference for the continuous case. The discrete case — orthogonality
with respect to a design rather than a measure — was brought into statistical
computing by Forsythe (1957), and the recurrence is sometimes named for him
there. The integer contrasts tabulated for equally spaced levels are older, and
are the ones [Section 8.7](../ch08-estimability/07-contrasts.html) uses.

**Splines.**  The name and the theory begin with Schoenberg (1946), and the rank
condition in its sharp form with Schoenberg and Whitney (1953). Least squares
fitting with fixed knots grew in the 1970s: Poirier (1973) and Wold (1974) are
early accounts aimed at data analysis, Smith (1979) a readable introduction.
Wold's advice — few knots, several observations between them, at most one
extremum per interval — is still sound, and Harrell (2015, chapter 2) gives the
modern version. The truncated power basis and its B-spline remedy are treated by
de Boor (2001), and in statistical terms by Eilers and Marx (1996) and Ruppert,
Wand and Carroll (2003, chapter 3).

**Segmented regression.**  A knot at an unknown location is a nonlinear problem
and is not treated here: see Hudson (1966) for the fitting, Feder (1975) for the
asymptotics when the change point is identified, and Seber and Wild (1989,
chapter 9) for a survey. Cobb (1978) shows what a conditional analysis of the
Nile series of @exm-ply-nile-ladder requires.

**Where the thread goes.**  Two loose ends are tied in
[Chapter 43](../ch43-smoothing/index.html): a basis that is stable however many
knots there are, and a continuous rather than discrete flexibility parameter.
Fahrmeir, Kneib, Lang and Marx (2021, chapter 8), the source for the rest of
this part, treats them together, as does Hastie and Tibshirani (1990) from the
additive-model side.

**Data.**  The Mauna Loa carbon dioxide series and the Nile flows are
public-domain data sets distributed with statsmodels, already used in
@exm-lm-co2 and @exm-ss-nile-polynomials; the skewed design of
@exm-ply-knot-placement is simulated with a fixed seed.

## References

- Cobb, George W. (1978). The Problem of the Nile: Conditional Solution to a Changepoint Problem. *Biometrika* 65(2), 243–251.
- de Boor, Carl (2001). *A Practical Guide to Splines*. Revised edition. New York: Springer.
- Eilers, Paul H. C. and Marx, Brian D. (1996). Flexible Smoothing with B-splines and Penalties. *Statistical Science* 11(2), 89–121.
- Fahrmeir, Ludwig, Kneib, Thomas, Lang, Stefan and Marx, Brian D. (2021). *Regression: Models, Methods and Applications*. 2nd edition. Berlin: Springer.
- Feder, Paul I. (1975). On Asymptotic Distribution Theory in Segmented Regression Problems: Identified Case. *The Annals of Statistics* 3(1), 49–83.
- Forsythe, George E. (1957). Generation and Use of Orthogonal Polynomials for Data-Fitting with a Digital Computer. *Journal of the Society for Industrial and Applied Mathematics* 5(2), 74–88.
- Harrell, Frank E. (2015). *Regression Modeling Strategies*. 2nd edition. Cham: Springer.
- Hastie, Trevor and Tibshirani, Robert (1990). *Generalized Additive Models*. London: Chapman and Hall.
- Higham, Nicholas J. (2002). *Accuracy and Stability of Numerical Algorithms*. 2nd edition. Philadelphia: SIAM.
- Hudson, Derek J. (1966). Fitting Segmented Curves Whose Join Points Have to Be Estimated. *Journal of the American Statistical Association* 61(316), 1097–1129.
- Poirier, Dale J. (1973). Piecewise Regression Using Cubic Splines. *Journal of the American Statistical Association* 68(343), 515–524.
- Runge, Carl (1901). Über empirische Funktionen und die Interpolation zwischen äquidistanten Ordinaten. *Zeitschrift für Mathematik und Physik* 46, 224–243.
- Ruppert, David, Wand, Matthew P. and Carroll, Raymond J. (2003). *Semiparametric Regression*. Cambridge: Cambridge University Press.
- Schoenberg, Isaac J. (1946). Contributions to the Problem of Approximation of Equidistant Data by Analytic Functions. *Quarterly of Applied Mathematics* 4, 45–99 and 112–141.
- Schoenberg, Isaac J. and Whitney, Anne (1953). On Pólya Frequency Functions III: The Positivity of Translation Determinants with an Application to the Interpolation Problem by Spline Curves. *Transactions of the American Mathematical Society* 74(2), 246–259.
- Seber, George A. F. and Lee, Alan J. (2003). *Linear Regression Analysis*. 2nd edition. Hoboken, NJ: Wiley.
- Seber, George A. F. and Wild, Chris J. (1989). *Nonlinear Regression*. New York: Wiley.
- Smith, Patricia L. (1979). Splines as a Useful and Convenient Statistical Tool. *The American Statistician* 33(2), 57–62.
- Szegő, Gábor (1939). *Orthogonal Polynomials*. American Mathematical Society Colloquium Publications 23. New York: American Mathematical Society.
- Todd, John (1954). The Condition of the Finite Segments of the Hilbert Matrix. In *Contributions to the Solution of Systems of Linear Equations and the Determination of Eigenvalues*, National Bureau of Standards Applied Mathematics Series 39, 109–116.
- Trefethen, Lloyd N. (2019). *Approximation Theory and Approximation Practice*. Extended edition. Philadelphia: SIAM.
- Wold, Svante (1974). Spline Functions in Data Analysis. *Technometrics* 16(1), 1–11.
