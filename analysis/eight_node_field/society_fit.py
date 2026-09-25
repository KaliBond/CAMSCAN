"""Per-society agreement of node threshold ordering with the leave-one-out consensus of all others."""
import numpy as np
from scipy.stats import spearmanr, rankdata
from load import load, derive
from common import factor_fit

ALL = {**load("discovery"), **load("replication")}
R = {}
for k, (y, a) in ALL.items():
    D = derive(a)["D"]; mu, L, G, _ = factor_fit(D, 1)
    R[k] = rankdata(-mu / L[:, 0])
rows = sorted(((k, spearmanr(R[k], np.mean([R[j] for j in R if j != k], 0))[0]) for k in R), key=lambda r: -r[1])
for k, r in rows:
    print(f"{k:13s} rho={r:+.2f}")
print(f"median {np.median([r for _, r in rows]):+.2f}; negative in {sum(r < 0 for _, r in rows)}/{len(rows)}")
