"""Lagged directed edges after removing 3 factors, controlling own persistence (robust to constant residuals)."""
import sys
import numpy as np
from scipy.stats import binomtest
from load import load, derive, NODES
from common import factor_fit

MODE = sys.argv[1] if len(sys.argv) > 1 else "blind"
data = load("blind") if MODE == "blind" else {**load("discovery"), **load("replication")}
sgn = []
for k, (y, a) in data.items():
    D = derive(a)["D"]; mu, L, G, _ = factor_fit(D, 3); R = D - mu - G @ L.T
    S = np.full((8, 8), np.nan)
    for i in range(8):
        for j in range(8):
            if i == j: continue
            Z = np.column_stack([np.ones(len(R) - 1), R[:-1, i], R[:-1, j]])
            beta, *_ = np.linalg.lstsq(Z, R[1:, i], rcond=None)
            if abs(beta[2]) > 1e-10: S[i, j] = np.sign(beta[2])
    sgn.append(S)
sgn = np.array(sgn)
rows = []
for i in range(8):
    for j in range(8):
        if i == j: continue
        v = sgn[:, i, j]; v = v[~np.isnan(v)]; npos = int((v > 0).sum())
        rows.append((f"{NODES[j]}->{NODES[i]}", npos, len(v), binomtest(npos, len(v), .5).pvalue))
p = np.array([r[3] for r in rows]); o = np.argsort(p); n = len(p)
q = np.empty(n); q[o] = np.minimum(np.minimum.accumulate((p[o] * n / np.arange(1, n + 1))[::-1])[::-1], 1)
best = np.argsort(p)[:5]
print(f"[F12] {MODE}: top edges " + "; ".join(f"{rows[b][0]} {rows[b][1]}/{rows[b][2]} q={q[b]:.3f}" for b in best))
for lab in ("Helm->Shield", "Shield->Helm"):
    b = [x[0] for x in rows].index(lab); print(f"[F12] {MODE}: {lab} positive {rows[b][1]}/{rows[b][2]}, q={q[b]:.3f}")
print(f"[F12] {MODE}: edges with q<0.05: {int((q < .05).sum())}/56")
