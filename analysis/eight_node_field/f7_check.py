"""F7/F8 re-check with NaN-safe handling: societies whose half has < MIN split states are excluded
(an empty half makes every order tie, so its order is undefined)."""
import sys
from itertools import permutations
import numpy as np
from scipy.stats import kendalltau, wilcoxon
from load import load, derive

MODE = sys.argv[1] if len(sys.argv) > 1 else "blind"
MIN = 10
data = load("blind") if MODE == "blind" else {**load("discovery"), **load("replication")}
X = {k: (derive(a)["D"] > 0).astype(int) for k, (y, a) in data.items()}
ORD = np.array(list(permutations(range(8))))
BITS = (np.arange(256)[:, None] >> np.arange(8)) & 1
DIST = np.empty((len(ORD), 256), np.float32)
tot = BITS.sum(1)[None, :, None]; m = np.arange(1, 9)[None, None, :]
for c in range(0, len(ORD), 4032):
    cs = np.cumsum(BITS[:, ORD[c:c + 4032]].transpose(1, 0, 2), axis=2)
    DIST[c:c + 4032] = np.minimum(((m - cs) + (tot - cs)).min(2), tot[..., 0])
W8 = 2 ** np.arange(8)


def split(x):
    s = x.sum(1); return x[(s > 0) & (s < 8)]


def ranks(x):
    H = np.bincount((split(x) * W8).sum(1), minlength=256).astype(np.float32)
    md = DIST @ H; b = np.flatnonzero(md == md.min())
    return np.argsort(ORD[b], axis=1).mean(0)


def nested(rank, xs):
    idx = np.argsort(rank, kind="stable"); p = np.zeros_like(xs)
    for r, n in enumerate(xs.sum(1)): p[r, idx[:n]] = 1
    return p


rng = np.random.default_rng(3)
H = {}
for k, x in X.items():
    h = len(x) // 2
    if len(split(x[:h])) >= MIN and len(split(x[h:])) >= MIN:
        H[k] = (ranks(x[:h]), ranks(x[h:]))
ks = list(H)
same = np.mean([kendalltau(*H[k])[0] for k in ks])
cross = np.mean([kendalltau(H[a][0], H[b][1])[0] for a in ks for b in ks if a != b])
perm = [np.mean([kendalltau(H[ks[i]][0], H[ks[j]][1])[0] for i, j in enumerate(rng.permutation(len(ks)))]) for _ in range(2000)]
print(f"[F7] {MODE}: {len(ks)}/{len(X)} societies with >= {MIN} split states per half; tau same {same:.3f} vs different {cross:.3f}; "
      f"p={(1 + np.sum(np.array(perm) >= same)) / 2001:.4f}")

eo, es, ho, hs = [], [], [], []
for k in ks:
    x = X[k]; h = len(x) // 2
    sh = ranks(np.vstack([X[j] for j in X if j != k]))
    for own, te in ((H[k][0], x[h:]), (H[k][1], x[:h])):
        t = split(te); po, ps = nested(own, t), nested(sh, t)
        eo.append((po == t).all(1).mean()); es.append((ps == t).all(1).mean())
        ho.append(np.abs(po - t).sum(1).mean()); hs.append(np.abs(ps - t).sum(1).mean())
print(f"[F8] {MODE}: exact nested match own-half {np.mean(eo):.3f} vs shared {np.mean(es):.3f} (p={wilcoxon(eo, es).pvalue:.3f}, own better in {sum(a > b for a, b in zip(eo, es))}/{len(eo)}); "
      f"Hamming {np.mean(ho):.3f} vs {np.mean(hs):.3f} (p={wilcoxon(ho, hs).pvalue:.3f})")
