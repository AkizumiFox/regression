# Book-wide notation

Fixed for all chapters. Macros are in `latex/macros.tex` (they work in the PDF and on the web).
A chapter that needs a new symbol adds it here and there before using it.

| Type | Meaning |
|---|---|
| `\Y`, `\y` | response vector (random, observed), n × 1 |
| `\X`, `\x_j`, `\x_{(i)}` | model matrix n × p; its j-th column; its i-th row as a column |
| `\bbeta`, `\hbeta`, `\be` | coefficients, a least squares estimate, the error vector ε |
| `\hY`, `\he` | fitted values, residuals |
| `\A, \B, \G, \M, \Q, \R, \V, \bP, \bD, \bL, \bS, \bT, \bU, \bH` | bold matrices (`\bH` is the hat matrix, written `\bH` in Ch. 5 before Ch. 6 renames it `\M`; locally, `\bH_{\G}=\G\X\T\X` in Ch. 8 and a Householder reflection in Ch. 10) |
| `\bu, \bv, \bw, \bz, \bb, \br, \bmu, \bgamma, \blambda, \bzero, \bone` | bold vectors |
| `\bSigma`, `\bLambda`, `\mathbf{a}`, `\boldsymbol{\uptheta}` | Matrices and vectors are **bold upright**, Latin and Greek alike. Latin: `\mathbf{a}`. Lowercase Greek: `\boldsymbol{\uptheta}` (upright Greek from `upgreek`: `\upalpha`, `\upbeta`, `\upvarepsilon`, …); capital Greek is upright already: `\boldsymbol{\Sigma}`. Never `\bm`, and never `\boldsymbol{\theta}` (that is bold *italic*) |
| `\I` | identity. `\M` the orthogonal projection onto C(X), `\Mo` onto a reduced model |
| `\T` | transpose, rendered with `\top`. Never `^T`, `'` or `\mathsf{T}` |
| `\ginv` | generalized inverse, `\A\ginv` means A⁻ (any matrix with A A⁻ A = A) |
| `\C(\X)`, `\Null(\X)` | column space, null space |
| `\perpc`, `\dirsum` | orthogonal complement, direct sum |
| `\rank, \tr, \diag, \spn, \argmin, \vecop` | operators |
| `\E, \Var, \Cov` | expectation, variance, covariance (matrix). `\Cov(\bu,\bv)` is the cross-covariance |
| `\Normal_n(\bmu,\bSigma)`, `\iid` | normal distribution, "independent and identically distributed" |
| `\chi^2(r,\gamma)`, `F(r,s,\gamma)`, `t(r,\delta)` | noncentral distributions |
| `\norm{\cdot}`, `\inner{\cdot}{\cdot}` | Euclidean norm and inner product unless stated |
| `\Real` | the real numbers |

## Conventions

- **Noncentrality without a factor ½.** If Z ~ N_r(μ, I) then ZᵀZ ~ χ²(r, μᵀμ), with mean
  r + γ and variance 2r + 4γ. Some books (Rencher, Christensen) use γ/2; say so when citing them.
- **Nonnegative definite**, never "positive semidefinite" (except when naming other books' terms).
- **The multivariate normal is defined through linear combinations**, so singular covariances are
  included.
- Accents always take braces: `\hat{\y}`, `\bar{\x}`, `\tilde{\X}`.
- British spelling (centring, centred, behaviour, modelling), as in the chapters written so far.
- **Correlated errors (Part VII).** `\V` is the error covariance. In ch. 31 it is a *shape*,
  Cov(Y) = σ²V, with V known, known up to a few parameters, or nonnegative definite and
  singular; in ch. 32–33 it is the marginal covariance V = ZGZᵀ + R itself, with σ² absorbed.
  `\Z` is the random-effect design matrix, `\bu` the random effects, `\G` = Cov(u) and
  `\R` = Cov(ε). Variance components are `\boldsymbol{\uptheta}`.

- **Generalized linear models (Part VIII).** Fixed across chapters 34–37. The response
  distribution is an exponential dispersion family with *natural parameter* `\theta_i`
  (scalar; `\vartheta` is free for an unrelated parameter in an exercise), *cumulant
  function* `b(\theta)`, *dispersion* `\phi`, *prior weight* `w_i` and *variance function*
  `V(\mu)`, so `\mu_i=b'(\theta_i)` and `\Var(Y_i)=\phi b''(\theta_i)/w_i`. The *link* is
  `g`, its inverse `h=g^{-1}`, the *linear predictor* is `\eta_i=\x_{(i)}\T\bbeta` with
  vector `\boldsymbol{\upeta}=\X\bbeta`, and the mean vector is `\bmu`. In the fitting,
  `\W` is the diagonal matrix of *working weights* `w_ih'(\eta_i)^2/V(\mu_i)`,
  `\bD=\diag(d\mu_i/d\eta_i)` and `\bz` is the *working response*; `\X\T\W\X` is the
  expected information. `D` is the *deviance* and `D^{*}=D/\phi` the *scaled deviance*.
  Category probabilities are `\boldsymbol{\uppi}`; in chapter 36 `\boldsymbol{\uptheta}`
  is the vector of *cutpoints* of a cumulative-link model (not the variance components of
  Part VII) and `\boldsymbol{\uppsi}` collects all parameters of a fit. Chapter 37 writes
  `\kappa` for the negative binomial shape and `\alpha=1/\kappa` for its reciprocal.

- **Estimating equations, random effects and missingness (chapters 38–41).** The symbols of
  chapters 34–37 above are unchanged. Chapter 38 writes `Q(\mu;y)` for the *quasi-likelihood*
  of @def-ql-quasi and `\bU(\bbeta)` for an *estimating function*; the dispersion is the same
  `\phi`, estimated by `\hat{\phi}`, and the sandwich covariance keeps the bread/meat form of
  chapter 21. Chapter 40 keeps Part VII's `\Z`, `\bu\sim\Normal_q\{\bzero,\G(\boldsymbol{\uptheta})\}`
  and `\G` for the random effects of a GLMM, and writes `\R_i(\boldsymbol{\upalpha})` for the
  **working correlation** matrix of a GEE, with `\V_i=\phi\A_i^{1/2}\R_i\A_i^{1/2}` the working
  covariance and `\A_i=\diag\{V(\mu_{ij})\}` — a local use of `\R`, not the `\R=\Cov(\be)` of
  chapters 32–33, flagged in the chapter where it starts. In chapter 41 `\R` is the
  **missingness-indicator** array, with `R_{ij}=1` when the entry is recorded, `\mathbf{D}` the
  complete data, `\mathbf{D}_{\mathrm{obs}}` and `\mathbf{D}_{\mathrm{mis}}` its recorded and
  unrecorded parts, `\boldsymbol{\uptheta}` the parameter of the data model and
  `\boldsymbol{\uppsi}` that of the mechanism; `\bR_i` in section 40.4.6 is the same indicator
  in diagonal form.

- **Bases, penalties and smoothers (Part IX).** A *basis matrix* is `\B`, with
  \( B_{ij}=B_j(x_i) \) the \( j \)th basis function at the \( i \)th design point, and
  `\bgamma` its coefficients; unpenalized columns stay in `\X` with coefficients `\bbeta`.
  *Knots* are `\kappa`: chapter 42 writes \( \kappa_1<\dots<\kappa_K \) for the \( K \)
  interior knots of a piecewise polynomial, chapter 43 keeps `\kappa` for the *extended*
  sequence of a B-spline basis and writes \( t_1<\dots<t_K \) for the interior knots inside
  it; chapter 44 writes \( \kappa_l \) for the knots of a kriging term. `\kappa(\cdot)` with
  an argument is the condition number of @def-mat-condition-number, as in chapters 1 and 10.
  The *penalty matrix* is `\bP`, nonnegative definite, with *smoothing parameter*
  \( \lambda\ge0 \); the penalty is \( \lambda\bgamma\T\bP\bgamma \) and a difference matrix
  is `\bD_k`, so that \( \bP=\bD_k\T\bD_k \) for a P-spline. A fit that is linear in the
  response has *smoother matrix* `\bS_\lambda` (`\bS` when \( \lambda \) is fixed), and its
  *effective degrees of freedom* is \( \tr(\bS_\lambda) \), the trace of chapters 27 and 30.
  Chapter 43 writes \( K \) for a kernel function and \( h \) for a bandwidth in section 43.1
  only; from section 43.2 on \( K \) is the number of interior knots. Chapters 44 and 45
  index a term by \( j \): a design block `\Z_j` with penalty `\bP_j`, coefficients
  `\bgamma_j` and smoothing parameter \( \lambda_j \), so `\Z` keeps the meaning it has in
  Part VII, the design of the coefficients that are shrunk.

- **Beyond the mean (Part X).** The *quantile level* is \( \tau\in(0,1) \), the *check loss*
  is \( \rho_{\tau}(u)=u\{\tau-1\{u<0\}\} \) and \( Q_Y(\tau) \) is a quantile function;
  \( s(\tau)=1/f\{Q(\tau)\} \) is the sparsity function. Asymmetric squared loss gives
  *expectiles* at level \( \tau \). In a GAMLSS the \( k \)th parameter of the response
  distribution is \( \theta_k \) with its own predictor \( \boldsymbol{\upeta}_k \) — a local
  use of \( \theta \), not the natural parameter of @def-glm-model.
  The indicator of an event is written \( 1\{\cdot\} \) throughout.
