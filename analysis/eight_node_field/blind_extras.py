"""Blind-panel supplements: weighted mean-field vs Ising, node-order agreement with the
31-society corpus, residual vs scorer-noise variance (blind ENV), per-society consensus fit.

Usage: BLIND=... BLIND_ENV=... python blind_extras.py
"""
import os
import numpy as np, pandas as pd
from scipy.stats import spearmanr, rankdata
from sklearn.decomposition import FactorAnalysis
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import log_loss

from load import load, derive, NODES
from common import factor_fit

ENV = os.environ.get("BLIND_ENV", "/root/.claude/uploads/7fc0cf52-48f9-5e67-80f0-06feacaf7d85/d2fd1448-JUNO_Blind_motif_ENV.csv")
BL = load("blind"); CORP = {**load("discovery"), **load("replication")}
DD = {k: derive(a)["D"] for k, (y, a) in {**BL, **CORP}.items()}
nb = list(BL)

# weighted mean-field vs Ising, LOSO within blind
tot = {"MF": 0.0, "MFw": 0.0, "ISING": 0.0}; per = []
for k in nb:
    train = [j for j in nb if j != k]
    fa = FactorAnalysis(1, random_state=0).fit(np.vstack([DD[j] - DD[j].mean(0) for j in train]))
    w = fa.components_[0] / fa.noise_variance_; w *= np.sign(w.mean())
    g = {}
    for i in range(8):
        o = [j for j in range(8) if j != i]
        Xtr = np.vstack([(DD[j] > 0)[:, o] for j in train]).astype(float)
        ytr = np.concatenate([(DD[j] > 0)[:, i] for j in train]).astype(int)
        X = (DD[k] > 0)[:, o].astype(float); y = (DD[k] > 0)[:, i].astype(int)
        for m, f in (("MF", lambda A: A.sum(1, keepdims=True)), ("MFw", lambda A: (A @ w[o])[:, None]), ("ISING", lambda A: A)):
            ll = log_loss(y, LogisticRegression(max_iter=2000).fit(f(Xtr), ytr).predict_proba(f(X))[:, 1], labels=[0, 1]) * len(y)
            tot[m] += ll; g[m] = g.get(m, 0) + ll
    per.append(g["MFw"] - g["ISING"])
print(f"[signs] blind LOSO log-loss: MF {tot['MF']:.0f}, weighted MF {tot['MFw']:.0f}, Ising {tot['ISING']:.0f}; "
      f"Ising beats weighted-MF in {sum(p > 0 for p in per)}/{len(per)}")

# node order agreement
def ranks(S, f):
    return np.mean([rankdata(f(k)) for k in S], 0)
thr = lambda k: (lambda mu, L: -mu / L[:, 0])(*factor_fit(DD[k], 1)[:2])
lam = lambda k: factor_fit(DD[k], 1)[1][:, 0]
fp = lambda k: (DD[k] > 0).mean(0)
for name, f in (("threshold", thr), ("loading", lam), ("fraction positive", fp)):
    print(f"[order] {name:17s} blind vs corpus-31 mean-rank Spearman {spearmanr(ranks(nb, f), ranks(CORP, f))[0]:+.2f}")

# scorer noise vs residual
env = pd.read_csv(ENV)
print("[noise] residual variance / scorer-mean error variance ((SD_K^2+SD_S^2)/5), median over nodes")
r1, r2 = [], []
for k in nb:
    e = env[env.Society_ID == k]
    ev = ((e.SD_K ** 2 + e.SD_S ** 2) / 5).groupby(e.Node).mean().reindex(NODES).values
    D = DD[k]; out = []
    for q, store in ((1, r1), (2, r2)):
        mu, L, G, _ = factor_fit(D, q); rv = (D - mu - G @ L.T).var(0)
        store.append(np.median(rv / ev)); out.append(f"k={q} {np.median(rv / ev):.1f}")
    print(f"[noise]   {k}  noise var {ev.mean():.3f}; " + "; ".join(out))
print(f"[noise] median across societies: k=1 {np.median(r1):.1f}, k=2 {np.median(r2):.1f}")
# is scorer disagreement higher when the field is changing fast or near zero crossings?
dis = []
for k in nb:
    e = env[env.Society_ID == k].pivot(index="T_Offset", columns="Node", values="SD_K").reindex(columns=NODES).values
    D = DD[k]
    near = np.abs(D) <= 0.5
    dis.append((e[near].mean(), e[~near].mean()))
dis = np.array(dis)
print(f"[noise] mean SD_K near threshold (|K-S|<=0.5) {dis[:,0].mean():.3f} vs away {dis[:,1].mean():.3f}; higher near threshold in {(dis[:,0]>dis[:,1]).sum()}/12")

# per-society consensus fit
R = {k: rankdata(thr(k)) for k in nb}
print("[fit] threshold order vs leave-one-out blind consensus: " +
      ", ".join(f"{k[-1]}:{spearmanr(R[k], np.mean([R[j] for j in R if j != k], 0))[0]:+.2f}" for k in nb))
