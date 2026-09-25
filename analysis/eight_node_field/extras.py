"""Supplementary checks.

(1) Weighted mean-field vs pairwise Ising for signs (held-out, LOSO across all societies):
    is the Ising gain just heterogeneous node weights (which a factor model also implies)?
(2) Node-order replication: discovery vs replication mean ranks.
(3) What drives PC1 share (P11)? Ties, exact repeated rows, series length.
(4) Measurement floor: scorer disagreement (ENV files, 5 scorers) vs variance of the
    residual left after removing 1 and 2 common factors.
"""
import glob, os
import numpy as np, pandas as pd
from scipy.stats import spearmanr, rankdata
from sklearn.decomposition import FactorAnalysis
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import log_loss

from load import load, derive, NODES, WM
from common import factor_fit, codes, eff_repertoire

disc, repl = load("discovery"), load("replication")
ALL = {**disc, **repl}
DD = {k: derive(a)["D"] for k, (y, a) in ALL.items()}

# (1) -------------------------------------------------------------
tot = {"MF": 0.0, "MFw": 0.0, "ISING": 0.0}; per = []
for k in ALL:
    train = [j for j in ALL if j != k]
    Z = np.vstack([DD[j] - DD[j].mean(0) for j in train])
    fa = FactorAnalysis(1, random_state=0).fit(Z)
    w = fa.components_[0] / fa.noise_variance_          # factor-score weights
    w = w * np.sign(w.mean())
    gain = {}
    for i in range(8):
        o = [j for j in range(8) if j != i]
        Xtr = np.vstack([(DD[j] > 0)[:, o] for j in train]).astype(float)
        ytr = np.concatenate([(DD[j] > 0)[:, i] for j in train]).astype(int)
        X = (DD[k] > 0)[:, o].astype(float); y = (DD[k] > 0)[:, i].astype(int)
        for m, f in (("MF", lambda A: A.sum(1, keepdims=True)), ("MFw", lambda A: (A @ w[o])[:, None]), ("ISING", lambda A: A)):
            mdl = LogisticRegression(max_iter=2000).fit(f(Xtr), ytr)
            ll = log_loss(y, mdl.predict_proba(f(X))[:, 1], labels=[0, 1]) * len(y)
            tot[m] += ll; gain[m] = gain.get(m, 0) + ll
    per.append(gain["MFw"] - gain["ISING"])
print(f"[signs] held-out log-loss: MF {tot['MF']:.0f}, weighted MF {tot['MFw']:.0f}, Ising {tot['ISING']:.0f}; "
      f"Ising beats weighted-MF in {sum(p > 0 for p in per)}/{len(per)} societies")

# (2) -------------------------------------------------------------
def thr_rank(S):
    R = []
    for k in S:
        D = DD[k]; mu, L, G, _ = factor_fit(D, 1)
        R.append(rankdata(-mu / L[:, 0]))
    return np.mean(R, 0)
rd, rr = thr_rank(disc), thr_rank(repl)
pos = lambda S: np.mean([rankdata((DD[k] > 0).mean(0)) for k in S], 0)
print(f"[order] threshold mean-rank agreement discovery vs replication: Spearman {spearmanr(rd, rr)[0]:.2f}; "
      f"fraction-positive rank agreement {spearmanr(pos(disc), pos(repl))[0]:.2f}")
print("[order] fraction positive (disc/repl): " + ", ".join(
    f"{n}:{np.median([(DD[k]>0).mean(0)[i] for k in disc]):.2f}/{np.median([(DD[k]>0).mean(0)[i] for k in repl]):.2f}" for i, n in enumerate(NODES)))

# (3) -------------------------------------------------------------
rows = []
for k, (y, a) in ALL.items():
    d = derive(a); D = d["D"]
    rows.append(dict(soc=k, set="disc" if k in disc else "repl", T=len(D), share=factor_fit(D, 1)[3][0],
                     ties=(D == 0).mean(), repeat=(np.abs(np.diff(a.reshape(len(a), -1), axis=0)).sum(1) == 0).mean(),
                     eff=eff_repertoire(codes(d["X"]))))
df = pd.DataFrame(rows)
print("[P11 drivers] Spearman with PC1 share: " + ", ".join(f"{c} {spearmanr(df.share, df[c])[0]:+.2f}" for c in ("ties", "repeat", "T", "eff")))
print("[P11 drivers] Spearman with eff repertoire: " + ", ".join(f"{c} {spearmanr(df.eff, df[c])[0]:+.2f}" for c in ("ties", "repeat", "T")))
import statsmodels.formula.api as smf
df["lr"] = np.log(df.eff)
fit = smf.ols("lr ~ share + repeat + np.log(T)", data=df).fit()
print(f"[P11 drivers] OLS log(eff) ~ share + repeat + log T (n={len(df)}): " +
      ", ".join(f"{p} {fit.params[p]:+.2f} (p={fit.pvalues[p]:.3f})" for p in fit.params.index[1:]))
print(df.sort_values("share", ascending=False).round(3).to_string(index=False))

# (4) -------------------------------------------------------------
print("[noise] residual variance after k factors vs scorer-mean error variance (K_sd^2+S_sd^2)/5")
ratios = {1: [], 2: []}
for k in disc:
    f = glob.glob(f"{WM}/data/v2.3/canonical/{k}_ENV_*_cleaned.csv")
    if not f: continue
    env = pd.read_csv(f[0]); yrs = ALL[k][0]
    env = env[env.Year.isin(yrs)]
    ev = ((env.K_sd ** 2 + env.S_sd ** 2) / 5).groupby(env.Node).mean().reindex(NODES).values
    D = DD[k]
    line = []
    for q in (1, 2):
        mu, L, G, _ = factor_fit(D, q)
        rv = (D - mu - G @ L.T).var(0)
        ratios[q].append(np.median(rv / ev)); line.append(f"k={q}: median resid/noise {np.median(rv/ev):.1f}")
    print(f"[noise]   {k:10s} noise var {ev.mean():.3f}; " + "; ".join(line))
print(f"[noise] across societies median ratio k=1 {np.median(ratios[1]):.1f}, k=2 {np.median(ratios[2]):.1f}")
