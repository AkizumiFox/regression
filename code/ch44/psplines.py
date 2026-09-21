"""Chapter 44: the chapter's toolkit for structured additive models.

Every model term is a pair (design matrix, penalty matrix). This module builds the common
terms -- P-spline blocks, varying-coefficient blocks, Markov random fields, radial
(kriging-type) bases and random-effect blocks -- fits the resulting model by one penalized
solve or by penalized IRLS for a generalized response, and chooses the smoothing parameters
by the restricted-likelihood updates of Section 44.5. Nothing here is specific to one
example, so the other scripts of the chapter import it.
"""
# <<basis>>
import numpy as np
from scipy.interpolate import BSpline


def bspline_basis(x, n_basis=15, degree=3, lo=None, hi=None):
    """n_basis B-splines of the given degree on equally spaced knots covering the data."""
    lo = x.min() if lo is None else lo
    hi = x.max() if hi is None else hi
    inner = np.linspace(lo, hi, n_basis - degree + 1)
    h = inner[1] - inner[0]
    knots = np.r_[lo - h * np.arange(degree, 0, -1), inner, hi + h * np.arange(1, degree + 1)]
    return BSpline.design_matrix(np.clip(x, lo, hi), knots, degree, extrapolate=False).toarray()


def difference_penalty(d, order=2):
    """K = D'D for the order-th difference matrix D: the P-spline penalty."""
    D = np.diff(np.eye(d), n=order, axis=0)
    return D.T @ D


def spline_term(x, n_basis=15, order=2, weight=None):
    """A P-spline term: its design, its penalty, and the map from new covariate values to
    design rows. The constraint sum_i f(x_i) = 0 is absorbed into the basis, and
    weight = u turns the term into the varying-coefficient term f(x) u."""
    lo, hi = x.min(), x.max()
    B = bspline_basis(x, n_basis, lo=lo, hi=hi)
    Q, _ = np.linalg.qr(B.sum(axis=0)[:, None], mode="complete")
    U = Q[:, 1:]                                  # a basis of {gamma : sum_i f(x_i) = 0}
    Z = B if weight is None else B * np.asarray(weight)[:, None]

    def design(xnew, wnew=None):
        Bn = bspline_basis(np.asarray(xnew, float), n_basis, lo=lo, hi=hi)
        return (Bn if wnew is None else Bn * np.asarray(wnew)[:, None]) @ U

    return Z @ U, U.T @ difference_penalty(n_basis, order) @ U, design


def spline_block(x, n_basis=15, order=2, weight=None):
    """The (design, penalty) pair of spline_term, without the prediction map."""
    return spline_term(x, n_basis, order, weight)[:2]
# <</basis>>


# <<spatial>>
def centre(Z, K):
    """Absorb the identifiability constraint 1'Z gamma = 0 into a basis."""
    Q, _ = np.linalg.qr(Z.sum(axis=0)[:, None], mode="complete")
    U = Q[:, 1:]
    return Z @ U, U.T @ K @ U


def markov_block(region, neighbours, d):
    """Areal spatial term: an incidence design and the neighbourhood penalty K = N - B."""
    Z = np.zeros((len(region), d))
    Z[np.arange(len(region)), region] = 1.0
    K = np.zeros((d, d))
    for s, nbrs in enumerate(neighbours):
        K[s, s] = len(nbrs)
        for t in nbrs:
            K[s, t] = -1.0
    return centre(Z, K)


def radial_block(coords, knots, rho=None):
    """Kriging-type term: the radial basis r(||s - k||), penalized by the same function
    evaluated between the knots."""
    def r(h):
        h = np.maximum(h, 1e-12) / rho
        return np.where(h < 1, (1 - h) ** 4 * (1 + 4 * h), 0.0)   # Wendland correlation

    dist = lambda A, B: np.sqrt(((A[:, None, :] - B[None, :, :]) ** 2).sum(-1))
    if rho is None:
        rho = 1.5 * dist(knots, knots).max() / np.sqrt(len(knots))
    return centre(r(dist(coords, knots)), r(dist(knots, knots)))


def random_block(group, d):
    """A random-effect term: a 0/1 incidence design with the ridge penalty K = I."""
    Z = np.zeros((len(group), d))
    Z[np.arange(len(group)), group] = 1.0
    return centre(Z, np.eye(d))
# <</spatial>>


# <<fit>>
def assemble(blocks, X, lambdas):
    """Stack the parametric columns and the terms into one design and one block penalty."""
    Z = np.column_stack([X] + [Zj for Zj, _ in blocks])
    K = np.zeros((Z.shape[1], Z.shape[1]))
    start = X.shape[1]
    for (_, Kj), lam in zip(blocks, lambdas):
        d = Kj.shape[0]
        K[start:start + d, start:start + d] = lam * Kj
        start += d
    return Z, K


class Fit:
    """One penalized solve: coefficients, components, operators and degrees of freedom."""

    def __init__(self, y, blocks, X, lambdas, W=None, z=None, operators=True):
        Z, K = assemble(blocks, X, lambdas)
        w = np.ones(len(y)) if W is None else W
        Zw = Z * w[:, None]
        B = Z.T @ Zw                                           # the weighted cross-product
        A = B + K
        Ainv = np.linalg.inv(A)
        self.Z, self.K, self.A, self.w = Z, K, A, w
        self.gamma = Ainv @ (Zw.T @ (y if z is None else z))
        self.fitted = Z @ self.gamma
        self.edf_total = np.trace(Ainv @ B)                    # tr(S) without forming S
        self.parts, self.edf, self.comp = [], [], []
        start = X.shape[1]
        for Zj, _ in blocks:
            sl = slice(start, start + Zj.shape[1])
            self.parts.append(Zj @ self.gamma[sl])
            self.edf.append(np.trace(Ainv[sl] @ B[:, sl]))
            start += Zj.shape[1]
        if operators:                                          # the n x n smoothers
            C = Ainv @ Zw.T
            self.S = Z @ C
            start = X.shape[1]
            for Zj, _ in blocks:
                self.comp.append(Zj @ C[start:start + Zj.shape[1]])
                start += Zj.shape[1]
        self.beta = self.gamma[:X.shape[1]]
# <</fit>>


# <<reml>>
def null_dim(K, tol=1e-9):
    """The dimension of the penalty null space: the part of a term that is never penalized."""
    e = np.linalg.eigvalsh(K)
    return int(np.sum(e < tol * max(e.max(), 1.0)))


def reml_lambdas(y, blocks, X, lambdas=None, iters=80, tol=1e-8):
    """Restricted-likelihood updates for a normal response: refit, then set each variance
    component to tau_j^2 = gamma_j' K_j gamma_j / (df_j - m_j) and sigma^2 to the penalized
    residual sum of squares over n - df, and put lambda_j = sigma^2 / tau_j^2."""
    lam = np.ones(len(blocks)) if lambdas is None else np.asarray(lambdas, float)
    for _ in range(iters):
        fit = Fit(y, blocks, X, lam, operators=False)
        sigma2 = np.sum((y - fit.fitted) ** 2) / (len(y) - fit.edf_total)
        new = lam.copy()
        start = X.shape[1]
        for j, (_, Kj) in enumerate(blocks):
            d = Kj.shape[0]
            g = fit.gamma[start:start + d]
            tau2 = max(g @ Kj @ g, 1e-12) / max(fit.edf[j] - null_dim(Kj), 1e-6)
            new[j] = min(max(sigma2 / tau2, 1e-8), 1e10)
            start += d
        done = np.max(np.abs(np.log(new) - np.log(lam))) < tol
        lam = new
        if done:
            break
    return lam, Fit(y, blocks, X, lam), sigma2
# <</reml>>


# <<pirls>>
def penalized_irls(y, blocks, X, lambdas, family="poisson", offset=0.0, iters=60, tol=1e-10):
    """Penalized IRLS: at each step a penalized weighted least squares fit of the working
    response on the same design -- Chapter 34's IRLS with the penalty added."""
    eta = np.log(np.maximum(y, 0.5)) if family == "poisson" else np.zeros(len(y))
    for _ in range(iters):
        if family == "poisson":
            mu = np.exp(eta + offset)
            w = dmu = mu                            # w = (dmu/deta)^2 / V(mu) = mu
        else:                                       # binomial with the logit link
            mu = 1 / (1 + np.exp(-(eta + offset)))
            w = dmu = mu * (1 - mu)
        w = np.maximum(w, 1e-10)
        z = eta + (y - mu) / dmu                    # the working response
        fit = Fit(z, blocks, X, lambdas, W=w, z=z, operators=False)
        done = np.max(np.abs(fit.fitted - eta)) < tol
        eta = fit.fitted
        if done:
            break
    fit.eta, fit.mu = eta, mu
    return fit


def pearson_dispersion(y, fit, family="poisson"):
    """The Pearson estimate of the dispersion, with the effective degrees of freedom."""
    v = fit.mu if family == "poisson" else fit.mu * (1 - fit.mu)
    return np.sum((y - fit.mu) ** 2 / v) / (len(y) - fit.edf_total)


def poisson_deviance(y, mu):
    """The Poisson deviance of Chapter 34, with the 0 log 0 convention."""
    t = np.where(y > 0, y * np.log(np.maximum(y, 1e-300) / mu), 0.0)
    return 2 * np.sum(t - (y - mu))


def irls_reml(y, blocks, X, family="poisson", offset=0.0, outer=40, tol=1e-7,
              dispersion=False):
    """Penalized IRLS with the smoothing parameters updated at every step: each working
    problem gets the update of reml_lambdas, with the dispersion of Chapter 38 in place of
    the nominal phi = 1 when dispersion is True."""
    lam = np.ones(len(blocks))
    for _ in range(outer):
        fit = penalized_irls(y, blocks, X, lam, family=family, offset=offset)
        phi = pearson_dispersion(y, fit, family) if dispersion else 1.0
        new = lam.copy()
        start = X.shape[1]
        for j, (_, Kj) in enumerate(blocks):
            d = Kj.shape[0]
            g = fit.gamma[start:start + d]
            tau2 = max(g @ Kj @ g, 1e-12) / max(fit.edf[j] - null_dim(Kj), 1e-6)
            new[j] = min(max(phi / tau2, 1e-8), 1e10)
            start += d
        done = np.max(np.abs(np.log(new) - np.log(lam))) < tol
        lam = new
        if done:
            break
    fit = penalized_irls(y, blocks, X, lam, family=family, offset=offset)
    fit.phi = pearson_dispersion(y, fit, family)
    return lam, fit
# <</pirls>>


if __name__ == "__main__":
    # self-checks on the pieces the other scripts rely on
    rng = np.random.default_rng(1)
    x = rng.uniform(0, 1, 200)
    B = bspline_basis(x, 12)
    assert B.shape == (200, 12)
    assert np.allclose(B.sum(axis=1), 1.0)                      # partition of unity
    assert np.all(B >= 0)
    Z, K, design = spline_term(x, 12)
    assert Z.shape == (200, 11) and np.allclose(Z.sum(axis=0), 0.0, atol=1e-9)
    assert np.allclose(K, K.T) and null_dim(K) == 1
    assert np.allclose(design(x), Z)
    nb = [[1], [0, 2], [1]]
    Zs, Ks = markov_block(np.array([0, 1, 2, 2]), nb, 3)
    Kfull = np.array([[1.0, -1, 0], [-1, 2, -1], [0, -1, 1]])
    assert np.allclose(Kfull @ np.ones(3), 0.0)                 # constants are unpenalized
    assert np.min(np.linalg.eigvalsh(Ks)) > -1e-10
    print("psplines: basis, penalty and block constructions check out")
