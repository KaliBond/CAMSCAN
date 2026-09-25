"""Shared model helpers: per-society factor fits and the one-factor dynamic null."""
import numpy as np

from load import NODES

W8 = 2 ** np.arange(8)


def codes(X):
    """8-bit integer code for each sign configuration row."""
    return (X * W8).sum(1)


def eff_repertoire(c):
    _, n = np.unique(c, return_counts=True)
    p = n / n.sum()
    return float(np.exp(-(p * np.log(p)).sum()))


def switching(c):
    return float((c[1:] != c[:-1]).mean())


def ar1(x):
    x0, x1 = x[:-1] - x.mean(), x[1:] - x.mean()
    phi = float((x0 * x1).sum() / (x0 * x0).sum()) if (x0 * x0).sum() > 0 else 0.0
    phi = float(np.clip(phi, -0.99, 0.99))
    e = x1 - phi * x0
    return phi, float(e.std())


def factor_fit(D, k=1):
    """PCA on within-society centred D (unscaled). Returns mu, L (8xk), G (Txk, unit variance), share."""
    mu = D.mean(0)
    Z = D - mu
    U, s, Vt = np.linalg.svd(Z, full_matrices=False)
    share = s ** 2 / (s ** 2).sum()
    G = U[:, :k] * np.sqrt(len(D))            # unit-variance scores
    L = Vt[:k].T * s[:k] / np.sqrt(len(D))     # loadings in D units
    for j in range(k):                        # orient: mean loading positive
        if L[:, j].mean() < 0:
            L[:, j] *= -1; G[:, j] *= -1
    return mu, L, G, share


def fit_dynamic_1f(D):
    mu, L, G, share = factor_fit(D, 1)
    R = D - mu - G @ L.T
    gphi, gsd = ar1(G[:, 0])
    rpar = [ar1(R[:, i]) for i in range(8)]
    return dict(mu=mu, L=L[:, 0], gphi=gphi, gsd=gsd, rpar=rpar, share=share)


def simulate_1f(p, T, rng, burn=100):
    n = T + burn
    g = np.zeros(n); r = np.zeros((n, 8))
    eg = rng.standard_normal(n) * p["gsd"]
    er = rng.standard_normal((n, 8))
    for t in range(1, n):
        g[t] = p["gphi"] * g[t - 1] + eg[t]
        for i in range(8):
            phi, sd = p["rpar"][i]
            r[t, i] = phi * r[t - 1, i] + sd * er[t, i]
    g, r = g[burn:], r[burn:]
    g = (g - g.mean()) / (g.std() + 1e-12)  # match unit-variance scale of fitted G
    return p["mu"] + np.outer(g, p["L"]) + r
