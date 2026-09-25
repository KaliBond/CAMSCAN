"""Shared-parameter dynamic factor vs relational grammar, scored on held-out societies.

All models are fitted on training societies only and evaluated on societies they never saw.

Contemporaneous (leave-one-node-out) test
    Predict node i's K-S in a held-out society from the other 7 nodes at the same time.
    Every model is a covariance structure; the prediction is the Gaussian conditional mean
    Sigma[i,-i] Sigma[-i,-i]^-1 x[-i]:
      FA1 / FA2 / FA3 : shared k-factor structure (8k+8 parameters, node-specific uniquenesses)
      FULL            : unconstrained shared covariance = explicit node-to-node relations (36 parameters)
    Sign-level analogue: logistic prediction of node i's sign from
      MF   : mean-field (number of other nodes positive; 2 params per node)
      ISING: pairwise couplings to each of the 7 other signs (8 params per node)

Dynamic (one-step-ahead) test on Delta D
      AR   : node-specific shared AR(1) coefficients
      AR+G : AR plus node-specific response to the shared FA1 field
      AR+G2: AR plus response to two shared fields
      VAR  : full shared 8x8 cross-lagged matrix (node-to-node lagged relations)
      AR+G+VAR_resid : field model plus cross-lagged relations among field residuals

Parameter-sharing test (within held-out society, first 60% -> last 40%)
      own-FA1   : society's own 1-factor structure and node means
      shared-FA1: loadings/uniquenesses from other societies, node means from own first 60%
      shared-FA1+offsets: also node offsets shared; only the society's overall level is local
      own-FULL / shared-FULL: unconstrained covariance

Usage: python key_test.py [corpus|blind]
"""
import numpy as np
from sklearn.decomposition import FactorAnalysis
from sklearn.linear_model import LogisticRegression, Ridge
from sklearn.metrics import log_loss

from load import load, derive, NODES

import sys
MODE = sys.argv[1] if len(sys.argv) > 1 else "corpus"
disc, repl = load("discovery"), load("replication")
ALL = {**disc, **repl}
BL = load("blind") if MODE == "blind" else {}
D = {k: derive(a)["D"] for k, (y, a) in {**ALL, **BL}.items()}
rng = np.random.default_rng(1)


def centred(k):
    return D[k] - D[k].mean(0)


def cov_full(train):
    Z = np.vstack([centred(k) for k in train])
    return np.cov(Z.T)


def cov_fa(train, q):
    Z = np.vstack([centred(k) for k in train])
    fa = FactorAnalysis(n_components=q, random_state=0).fit(Z)
    return fa.components_.T @ fa.components_ + np.diag(fa.noise_variance_)


def cond_pred(Sig, Z):
    P = np.zeros_like(Z)
    for i in range(8):
        o = [j for j in range(8) if j != i]
        b = np.linalg.solve(Sig[np.ix_(o, o)], Sig[o, i])
        P[:, i] = Z[:, o] @ b
    return P


def r2(Y, P):
    return 1 - ((Y - P) ** 2).sum() / (Y ** 2).sum()


def signacc(k, P):
    mu = D[k].mean(0)
    return ((mu + P > 0) == (D[k] > 0)).mean()


def contemporaneous(train, test):
    models = {"FA1": cov_fa(train, 1), "FA2": cov_fa(train, 2), "FA3": cov_fa(train, 3), "FULL": cov_full(train)}
    res = {m: [] for m in models}; acc = {m: [] for m in models}
    for k in test:
        Z = centred(k)
        for m, Sig in models.items():
            P = cond_pred(Sig, Z); res[m].append(r2(Z, P)); acc[m].append(signacc(k, P))
    return res, acc


def sign_models(train, test):
    out = {"MF": [], "ISING": []}
    for i in range(8):
        o = [j for j in range(8) if j != i]
        Xtr = np.vstack([(D[k] > 0)[:, o] for k in train]).astype(float)
        ytr = np.concatenate([(D[k] > 0)[:, i] for k in train]).astype(int)
        mf = LogisticRegression(max_iter=2000).fit(Xtr.sum(1, keepdims=True), ytr)
        isg = LogisticRegression(max_iter=2000).fit(Xtr, ytr)
        for k in test:
            X = (D[k] > 0)[:, o].astype(float); y = (D[k] > 0)[:, i].astype(int)
            out["MF"].append(log_loss(y, mf.predict_proba(X.sum(1, keepdims=True))[:, 1], labels=[0, 1]) * len(y))
            out["ISING"].append(log_loss(y, isg.predict_proba(X)[:, 1], labels=[0, 1]) * len(y))
    return {m: float(np.sum(v)) for m, v in out.items()}


def fa_scores(train, q):
    Z = np.vstack([centred(k) for k in train])
    return FactorAnalysis(n_components=q, random_state=0).fit(Z)


def design(k, fa1, fa2, kind):
    Z = centred(k); x, y = Z[:-1], np.diff(Z, axis=0)
    if kind == "AR":
        F = [x]
    elif kind == "AR+G":
        F = [x, fa1.transform(x)]
    elif kind == "AR+G2":
        F = [x, fa2.transform(x)]
    elif kind == "VAR":
        F = [x]
    elif kind == "AR+G+VARres":
        g = fa1.transform(x); F = [x, g, x - g @ fa1.components_]
    return np.hstack(F), x, y


def dynamic(train, test):
    fa1, fa2 = fa_scores(train, 1), fa_scores(train, 2)
    res = {}
    for kind in ("AR", "AR+G", "AR+G2", "VAR", "AR+G+VARres"):
        sse = sst = 0.0
        for i in range(8):
            def cols(F, x):
                if kind in ("VAR",):
                    return F
                if kind == "AR+G+VARres":  # own lag, field, other nodes' field residuals
                    o = [j for j in range(8) if j != i]
                    return np.column_stack([x[:, i], F[:, 8:9], F[:, 9:][:, o]])
                return np.column_stack([x[:, i], F[:, 8:]])
            tr = [design(k, fa1, fa2, kind) for k in train]
            Xtr = np.vstack([cols(F, x) for F, x, y in tr]); ytr = np.concatenate([y[:, i] for F, x, y in tr])
            m = Ridge(alpha=1.0).fit(Xtr, ytr)
            for k in test:
                F, x, y = design(k, fa1, fa2, kind)
                p = m.predict(cols(F, x)); sse += ((y[:, i] - p) ** 2).sum(); sst += (y[:, i] ** 2).sum()
        res[kind] = 1 - sse / sst
    return res


def sharing(train, k, frac=0.6):
    """Time split within held-out society k; returns MSE and sign accuracy on the last 40%."""
    Dk = D[k]; T = len(Dk); c = int(T * frac)
    a, b = Dk[:c], Dk[c:]
    mu_own = a.mean(0)
    Zown = a - mu_own
    out = {}
    # shared offsets: node deviation from society mean, averaged over training societies
    delta = np.mean([D[j].mean(0) - D[j].mean() for j in train], 0)
    mu_off = a.mean() + delta
    Sig_own_fa = FactorAnalysis(1, random_state=0).fit(Zown)
    Sig_own_fa = Sig_own_fa.components_.T @ Sig_own_fa.components_ + np.diag(Sig_own_fa.noise_variance_)
    specs = {
        "own-FA1": (mu_own, Sig_own_fa),
        "own-FULL": (mu_own, np.cov(Zown.T) + 1e-6 * np.eye(8)),
        "shared-FA1": (mu_own, cov_fa(train, 1)),
        "shared-FULL": (mu_own, cov_full(train)),
        "shared-FA1+offsets": (mu_off, cov_fa(train, 1)),
    }
    for m, (mu, Sig) in specs.items():
        P = mu + cond_pred(Sig, b - mu)
        out[m] = (((b - P) ** 2).mean(), ((P > 0) == (b > 0)).mean())
    return out


def report(tag, train_sets):
    print(f"\n==== {tag} ====")
    agg_r, agg_a, wins = {}, {}, {}
    sgn = {"MF": 0.0, "ISING": 0.0}
    dyn = []
    for train, test in train_sets:
        r, a = contemporaneous(train, test)
        for m in r:
            agg_r.setdefault(m, []).extend(r[m]); agg_a.setdefault(m, []).extend(a[m])
        s = sign_models(train, test)
        for m in s: sgn[m] += s[m]
        dyn.append(dynamic(train, test))
    n = len(agg_r["FA1"])
    print(f"Contemporaneous leave-one-node-out, {n} held-out societies (mean R2 | mean sign acc | #societies FULL beats model):")
    for m in agg_r:
        beats = sum(f > x for f, x in zip(agg_r["FULL"], agg_r[m]))
        print(f"  {m:5s} R2={np.mean(agg_r[m]):.3f} | acc={np.mean(agg_a[m]):.3f} | FULL better in {beats}/{n}")
    print(f"Sign-level held-out log-loss (total nats): mean-field {sgn['MF']:.1f}  pairwise Ising {sgn['ISING']:.1f}  "
          f"(Ising gain {sgn['MF'] - sgn['ISING']:+.1f})")
    print("One-step Delta D out-of-sample R2 (pooled over held-out societies):")
    for m in dyn[0]:
        print(f"  {m:12s} {np.mean([d[m] for d in dyn]):+.4f}")
    return agg_r


def sharing_table(pairs, title):
    print(f"\n==== {title} ====")
    tab = {}
    for train, k in pairs:
        for m, (mse, acc) in sharing(train, k).items():
            tab.setdefault(m, []).append((mse, acc))
    base = np.array([x[0] for x in tab["own-FA1"]])
    for m, v in tab.items():
        v = np.array(v)
        print(f"  {m:20s} median MSE ratio vs own-FA1 {np.median(v[:, 0] / base):.3f}  mean sign acc {v[:, 1].mean():.3f}  "
              f"better than own-FA1 in {(v[:, 0] < base).sum()}/{len(base)}")


if MODE == "blind":
    nb = list(BL)
    report("LOSO within blind (12)", [([j for j in nb if j != k], [k]) for k in nb])
    report("Train 31 corpus societies -> test blind 12", [(list(ALL), nb)])
    sharing_table([([j for j in nb if j != k], k) for k in nb], "Parameter sharing within blind (shared = other blind societies)")
    sharing_table([(list(ALL), k) for k in nb], "Parameter sharing: shared parameters from the 31 corpus societies")
    sys.exit()

names_d, names_r = list(disc), list(repl)
# 1) leave-one-society-out within discovery
report("LOSO within discovery (11)", [([j for j in names_d if j != k], [k]) for k in names_d])
# 2) train on discovery, test on 20 unseen replication societies
report("Train discovery -> test replication (20 unseen)", [(names_d, names_r)])
# 3) leave-one-society-out across all 31
allr = report("LOSO across all 31", [([j for j in ALL if j != k], [k]) for k in ALL])

# 4) parameter sharing within held-out society
sharing_table([([j for j in ALL if j != k], k) for k in ALL],
              "Parameter sharing: fit on first 60% (own) or on other societies (shared); test last 40%")
