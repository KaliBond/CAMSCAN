"""Robustness checks for the proposition battery.

(a) Tie coding: K==S coded as positive instead of non-positive.
(b) Series-length confound: morphology (P10/P11) on a common 1950-2025 window.
(c) Field-preserving surrogate: keep the observed common field(s) G, loadings and
    each node's residual series, but circularly shift residuals independently per node.
    This destroys only cross-node alignment of what the field does not explain.

Usage: python robustness.py [discovery|replication|all]
"""
import sys
import numpy as np
from scipy.stats import spearmanr

from load import load, derive, NODES
from common import codes, eff_repertoire, switching, factor_fit

which = sys.argv[1] if len(sys.argv) > 1 else "discovery"
rng = np.random.default_rng(7)
data = {**load("discovery"), **load("replication")} if which == "all" else load(which)
S = {k: derive(a) for k, (y, a) in data.items()}
names = list(S)

# (a) tie coding
for rule, f in (("K==S -> non-positive", lambda D: D > 0), ("K==S -> positive", lambda D: D >= 0)):
    er = [eff_repertoire(codes(f(S[k]["D"]).astype(int))) for k in names]
    sw = [switching(codes(f(S[k]["D"]).astype(int))) for k in names]
    nrep = [len(np.unique(codes(f(S[k]["D"]).astype(int)))) for k in names]
    print(f"[ties] {rule:22s} mean repertoire {np.mean(nrep):.1f}, eff {np.mean(er):.2f}, switching {np.mean(sw):.3f}, rho(er,sw)={spearmanr(er, sw)[0]:.2f}")

# (b) common window
print("[window] 1950-2025 common window")
er, sw, sh, keep = [], [], [], []
for k, (yrs, a) in data.items():
    m = (yrs >= 1950) & (yrs <= 2025)
    if m.sum() < 70: continue
    d = derive(a[m]); c = codes(d["X"])
    er.append(eff_repertoire(c)); sw.append(switching(c)); sh.append(factor_fit(d["D"], 1)[3][0]); keep.append(k)
er, sw, sh = map(np.array, (er, sw, sh))
print(f"[window] n={len(keep)} eff rep {er.min():.2f}-{er.max():.2f}, switching {sw.min():.3f}-{sw.max():.3f}")
print(f"[window] rho(er,sw)={spearmanr(er, sw)[0]:.2f}  rho(PC1share,er)={spearmanr(sh, er)[0]:.2f}  rho(PC1share,sw)={spearmanr(sh, sw)[0]:.2f}")
for k, a, b, c in sorted(zip(keep, sh, er, sw), key=lambda r: -r[1]):
    print(f"[window]   {k:13s} PC1={a:.2f} eff={b:5.2f} sw={c:.3f}")

# (c) field-preserving surrogates
print("[surrogate] observed repertoire / Hamming vs field-preserving surrogate (k factors kept)")
for q in (1, 2):
    lo, hlo, tot = 0, 0, []
    for k in names:
        D = S[k]["D"]; T = len(D)
        mu, L, G, _ = factor_fit(D, q)
        base = mu + G @ L.T; R = D - base
        obs = len(np.unique(codes(S[k]["X"])))
        hobs = np.abs(np.diff(S[k]["X"], axis=0)).sum(1).mean()
        sims, hs = [], []
        for _ in range(500):
            Rs = np.column_stack([np.roll(R[:, i], rng.integers(1, T)) for i in range(8)])
            Xs = (base + Rs > 1e-9).astype(int)
            sims.append(len(np.unique(codes(Xs)))); hs.append(np.abs(np.diff(Xs, axis=0)).sum(1).mean())
        p = np.mean(np.array(sims) <= obs); ph = np.mean(np.array(hs) <= hobs)
        lo += p < 0.05; hlo += ph < 0.05
        tot.append((k, obs, np.mean(sims), p, hobs, np.mean(hs)))
    print(f"[surrogate] k={q}: repertoire below 5% tail in {lo}/{len(names)}; Hamming below 5% tail in {hlo}/{len(names)}")
    print("[surrogate]   " + "; ".join(f"{k} {o}/{m:.1f}" for k, o, m, p, h, hm in tot))
