"""How well identified is the shared threshold order (blind panel)? Compares the user's order,
the recovery-crossing order and the chain-optimal order on pooled split states."""
from itertools import permutations
import numpy as np
from scipy.stats import kendalltau
from load import load, derive, NODES

B = load("blind"); W8 = 2 ** np.arange(8)
X = np.vstack([(derive(a)["D"] > 0).astype(int) for y, a in B.values()])
s = X.sum(1); X = X[(s > 0) & (s < 8)]


def dist(order):
    xp = X[:, order]; cs = np.cumsum(xp, 1); tot = xp.sum(1)[:, None]
    err = (np.arange(1, 9) - cs) + (tot - cs)
    return np.minimum(err.min(1), tot[:, 0]).mean()


ords = list(permutations(range(8)))
d = np.array([dist(list(o)) for o in ords])
idx = lambda names: [NODES.index(n) for n in names]
user = idx(["Shield", "Lore", "Archive", "Craft", "Helm", "Flow", "Stewards", "Hands"])
recov = idx(["Archive", "Shield", "Lore", "Helm", "Craft", "Flow", "Stewards", "Hands"])
best = ords[int(d.argmin())]
print(f"split states n={len(X)}; best distance {d.min():.3f}; random mean {d.mean():.3f}")
for lab, o in (("user order", user), ("recovery order", recov), ("chain optimum", list(best))):
    print(f"{lab:15s} dist {dist(o):.3f}  percentile {np.mean(d < dist(o)):.4f}  " + " > ".join(NODES[i] for i in o))
near = np.flatnonzero(d <= d.min() + 0.02)
pos = np.array([np.argsort(ords[i]) for i in near])  # position of each node
print(f"orders within 0.02 of optimum: {len(near)}; node position range among them:")
for i, n in enumerate(NODES):
    print(f"  {n:9s} positions {pos[:, i].min()}-{pos[:, i].max()} (mean {pos[:, i].mean():.1f})")
rank = lambda o: np.argsort(o)
print(f"Kendall tau user vs recovery order {kendalltau(rank(user), rank(recov))[0]:.2f}; user vs optimum {kendalltau(rank(user), rank(list(best)))[0]:.2f}")
