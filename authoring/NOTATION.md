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
| `\bSigma`, `\bLambda`, `\bm{...}` | bold Greek; `\bm` for anything else |
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
