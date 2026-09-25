"""Independent validation of the round-2 findings (conditional grammar, nonlinear grammar,
sign-state grammar, node-local memory, threshold chains, society-specific orderings,
recovery hierarchy, recovery-window couplings, innovation covariance, lagged edges).

All predictive comparisons are leave-one-society-out (LOSO): models are fitted on the other
societies and scored on the held-out one. Latent factors are each society's own PCA on
centred K-S; to stay rotation-invariant across societies, pooled models use each node's
k-factor fitted value F^k_i (and residual R^k_i), plus the oriented first field G1.

Usage: python validate_round2.py [blind|corpus]
"""
import sys
from itertools import permutations
import numpy as np
from scipy.stats import wilcoxon, binomtest, kendalltau
from sklearn.linear_model import Ridge, LogisticRegression
from sklearn.ensemble import HistGradientBoostingRegressor
from sklearn.metrics import log_loss

from load import load, derive, NODES
from common import factor_fit

MODE = sys.argv[1] if len(sys.argv) > 1 else "blind"
data = load("blind") if MODE == "blind" else {**load("discovery"), **load("replication")}
names = list(data)
rng = np.random.default_rng(11)
D = {k: derive(a)["D"] for k, (y, a) in data.items()}
X = {k: (D[k] > 0).astype(int) for k in names}
M = {k: derive(a)["M"] for k, (y, a) in data.items()}


def recon(k, q):
    mu, L, G, _ = factor_fit(D[k], q)
    return mu + G @ L.T, G[:, 0]


REC = {(k, q): recon(k, q) for k in names for q in (1, 3)}


def p_w(a, b):
    """two-sided Wilcoxon on paired per-society values (a vs b)."""
    d = np.asarray(a) - np.asarray(b)
    return wilcoxon(d).pvalue if np.any(d != 0) else 1.0


def pr(tag, s):
    print(f"[{tag}] {s}", flush=True)


print(f"==== {MODE}: {len(names)} societies ====")

# ------------------------------------------------------------------ F1 conditional grammar
def f1_rows(k, q, kind):
    F, g1 = REC[(k, q)]; Dk = D[k]; R = Dk - F
    rows = {i: [] for i in range(8)}; ys = {i: [] for i in range(8)}
    for i in range(8):
        o = [j for j in range(8) if j != i]
        base = np.column_stack([Dk[:-1, i], F[:-1, i]])
        if kind == "base":
            Z = base
        elif kind == "linear":
            Z = np.column_stack([base, R[:-1][:, o]])
        elif kind == "field-regime":
            up = (g1[:-1] > 0)[:, None]
            Z = np.column_stack([base, R[:-1][:, o] * up, R[:-1][:, o] * ~up])
        elif kind == "sign-config":
            s = X[k][:-1][:, o]
            Z = np.column_stack([base, R[:-1][:, o] * s, R[:-1][:, o] * (1 - s)])
        rows[i] = Z; ys[i] = Dk[1:, i]
    return rows, ys


def loso_rmse(builder, model=lambda: Ridge(alpha=1e-3)):
    per = {}
    B = {k: builder(k) for k in names}
    for k in names:
        tr = [j for j in names if j != k]; se = n = 0.0
        for i in range(8):
            m = model().fit(np.vstack([B[j][0][i] for j in tr]), np.concatenate([B[j][1][i] for j in tr]))
            e = B[k][1][i] - m.predict(B[k][0][i]); se += (e ** 2).sum(); n += len(e)
        per[k] = np.sqrt(se / n)
    return per


for q in (1, 3):
    res = {kind: loso_rmse(lambda k, kind=kind: f1_rows(k, q, kind)) for kind in ("base", "linear", "field-regime", "sign-config")}
    b = np.array([res["base"][k] for k in names])
    line = []
    for kind, v in res.items():
        v = np.array([v[k] for k in names])
        line.append(f"{kind} {v.mean():.5f}" + ("" if kind == "base" else f" (better in {(v < b).sum()}/{len(b)}, p={p_w(v, b):.3f})"))
    pr("F1", f"k={q} RMSE D(t+1): " + "; ".join(line))

# ------------------------------------------------------------------ F2 nonlinear grammar
def f2_rows(k, others):
    F, g1 = REC[(k, 3)]; Dk = D[k]; rows, ys = {}, {}
    for i in range(8):
        cols = [F[:-1, i], Dk[:-1, i], g1[:-1]]
        if others:
            cols += [Dk[:-1, j] for j in range(8) if j != i]
        rows[i] = np.column_stack(cols); ys[i] = Dk[1:, i]
    return rows, ys


gbm = lambda: HistGradientBoostingRegressor(max_iter=200, learning_rate=0.05, max_leaf_nodes=15, random_state=0)
a = loso_rmse(lambda k: f2_rows(k, False), gbm); b = loso_rmse(lambda k: f2_rows(k, True), gbm)
va, vb = np.array([a[k] for k in names]), np.array([b[k] for k in names])
pr("F2", f"GBM RMSE: 3-factor+own {va.mean():.4f} vs +all other nodes {vb.mean():.4f}; others worse in {(vb > va).sum()}/{len(va)}, p={p_w(vb, va):.4f}")

# ------------------------------------------------------------------ F3 sign-state grammar
def loso_logit(builder, C=1.0):
    B = {k: builder(k) for k in names}
    out = {}
    for k in names:
        tr = [j for j in names if j != k]; P = np.zeros_like(B[k][1], dtype=float)
        for i in range(8):
            m = LogisticRegression(C=C, max_iter=3000).fit(np.vstack([B[j][0][i] for j in tr]), np.concatenate([B[j][1][:, i] for j in tr]))
            P[:, i] = m.predict_proba(B[k][0][i])[:, 1]
        out[k] = (B[k][1], np.clip(P, 1e-6, 1 - 1e-6))
    return out


def f3_rows(k, others):
    F, g1 = REC[(k, 3)]; Xk = X[k]; rows = {}
    for i in range(8):
        cols = [F[:-1, i], g1[:-1], Xk[:-1, i]]
        if others:
            cols += [Xk[:-1, j] for j in range(8) if j != i]
        rows[i] = np.column_stack(cols)
    return rows, Xk[1:]


def summarise(out):
    ll = {k: log_loss(y.ravel(), p.ravel(), labels=[0, 1]) for k, (y, p) in out.items()}
    acc = np.mean([((p > .5) == y).mean() for y, p in out.values()])
    ex = np.mean([((p > .5) == y).all(1).mean() for y, p in out.values()])
    return ll, acc, ex


la, acca, exa = summarise(loso_logit(lambda k: f3_rows(k, False)))
lb, accb, exb = summarise(loso_logit(lambda k: f3_rows(k, True)))
va, vb = np.array([la[k] for k in names]), np.array([lb[k] for k in names])
pr("F3", f"X(t+1): 3f+own sign logloss {va.mean():.4f} acc {acca:.4f} exact {exa:.3f} | +7 other signs logloss {vb.mean():.4f} "
         f"acc {accb:.4f} exact {exb:.3f}; worse in {(vb > va).sum()}/{len(va)}, p={p_w(vb, va):.4f}")

# ------------------------------------------------------------------ F4 node-local memory (current sign)
def g1_minus(k, i):
    o = [j for j in range(8) if j != i]
    Z = D[k][:, o] - D[k][:, o].mean(0)
    U, s, Vt = np.linalg.svd(Z, full_matrices=False)
    g = U[:, 0] * np.sqrt(len(Z))
    return g * np.sign(Vt[0].mean())


GM = {(k, i): g1_minus(k, i) for k in names for i in range(8)}


def f4_rows(k, kind):
    rows = {}
    for i in range(8):
        g = GM[(k, i)]
        cols = [g[1:]]
        if "lag" in kind: cols.append(X[k][:-1, i])
        if "dir" in kind: cols.append((np.diff(g) > 0).astype(float))
        rows[i] = np.column_stack(cols)
    return rows, X[k][1:]


res4 = {kind: summarise(loso_logit(lambda k, kind=kind: f4_rows(k, kind)))[0] for kind in ("field", "field+dir", "field+lag", "field+lag+dir")}
base = np.array([res4["field"][k] for k in names])
for kind, v in res4.items():
    v = np.array([v[k] for k in names])
    pr("F4", f"X(t) from leave-node-out field: {kind:14s} logloss {v.mean():.4f}" +
       ("" if kind == "field" else f" (better than field-only in {(v < base).sum()}/{len(v)}, p={p_w(v, base):.2e})"))
v1, v2 = np.array([res4["field+lag"][k] for k in names]), np.array([res4["field+lag+dir"][k] for k in names])
pr("F4", f"direction added to field+lag: better in {(v2 < v1).sum()}/{len(v1)}, p={p_w(v2, v1):.3f}")

# ------------------------------------------------------------------ F5-F8 threshold chains
ORD = np.array(list(permutations(range(8))))            # 40320 x 8, easiest -> hardest
codes = np.arange(256)
BITS = (codes[:, None] >> np.arange(8)) & 1             # 256 x 8
# distance of each config to nearest nested chain for each order: 40320 x 256 (computed in chunks)
DIST = np.empty((len(ORD), 256), dtype=np.int8)
tot = BITS.sum(1).astype(np.int8)[None, :, None]
m = np.arange(1, 9, dtype=np.int8)[None, None, :]
for c0 in range(0, len(ORD), 4032):
    cs = np.cumsum(BITS.astype(np.int8)[:, ORD[c0:c0 + 4032]].transpose(1, 0, 2), axis=2, dtype=np.int8)
    err = (m - cs) + (tot - cs)                         # chain whose first m nodes are positive
    DIST[c0:c0 + 4032] = np.minimum(err.min(2), tot[..., 0])  # m=0 chain: all negative
DIST = DIST.astype(np.float32)
del cs, err
W8 = 2 ** np.arange(8)


def hist(xs, split_only=False):
    c = (xs * W8).sum(1)
    if split_only:
        s = xs.sum(1); c = c[(s > 0) & (s < 8)]
    return np.bincount(c, minlength=256).astype(float)


def mean_dist(H):  # H: 256 histogram -> mean distance for every order
    return DIST @ H / H.sum()


def best_order(H):
    H = H.copy(); H[0] = H[255] = 0          # learn from split states only
    md = mean_dist(H); best = np.flatnonzero(md == md.min())
    return best, md


# F5: learned shared order on held-out societies vs random orders
pct, learned_d, rand_d, orders = [], [], [], []
for k in names:
    Htr = sum(hist(X[j]) for j in names if j != k)
    b, _ = best_order(Htr)
    mdk = mean_dist(hist(X[k]))
    learned_d.append(mdk[b].mean()); rand_d.append(mdk.mean())
    pct.append(np.mean(mdk <= mdk[b].mean())); orders.append(b)
null = np.array([np.mean([mean_dist(hist(X[k]))[rng.integers(len(ORD))] for k in names]) for _ in range(2000)])
p5 = (1 + (null <= np.mean(learned_d)).sum()) / 2001
Hall = sum(hist(X[k]) for k in names)
ball, mdall = best_order(Hall)
ranks = np.argsort(ORD[ball], axis=1).mean(0)   # mean position of each node
shared = [NODES[i] for i in np.argsort(ranks)]
pr("F5", f"held-out mean chain distance: learned order {np.mean(learned_d):.3f} vs random {np.mean(rand_d):.3f}; "
         f"mean held-out percentile {np.mean(pct):.3f}; permutation p={p5:.4f}")
pr("F5", f"pooled best order ({len(ball)} tied), easiest->hardest: {' > '.join(shared)}")

# F6: split configurations only
on, within1, dl, dr = [], [], [], []
for k, b in zip(names, orders):
    H = hist(X[k], split_only=True)
    if H.sum() == 0: continue
    dd = DIST[b].mean(0)                      # average over tied learned orders
    on.append((dd == 0) @ H / H.sum() if len(b) == 1 else (DIST[b] == 0).mean(0) @ H / H.sum())
    within1.append((DIST[b] <= 1).mean(0) @ H / H.sum()); dl.append(dd @ H / H.sum()); dr.append(mean_dist(H).mean())
pr("F6", f"split states: exactly on shared chain {np.mean(on):.3f}; within 1 flip {np.mean(within1):.3f}; "
         f"mean distance {np.mean(dl):.3f} vs random orders {np.mean(dr):.3f}")

# F7: split-half order stability
def order_ranks(H):
    b, _ = best_order(H)
    return np.argsort(ORD[b], axis=1).mean(0)   # mean position of each node, averaged over ties


halves = {}
for k in names:
    T = len(X[k]); h = T // 2
    halves[k] = (order_ranks(hist(X[k][:h])), order_ranks(hist(X[k][h:])))
same = np.mean([kendalltau(*halves[k])[0] for k in names])
cross = np.mean([kendalltau(halves[a][0], halves[b][1])[0] for a in names for b in names if a != b])
perm = []
for _ in range(2000):
    p = rng.permutation(len(names))
    perm.append(np.mean([kendalltau(halves[names[i]][0], halves[names[p[i]]][1])[0] for i in range(len(names))]))
p7 = (1 + (np.array(perm) >= same).sum()) / 2001
pr("F7", f"Kendall tau between halves: same society {same:.3f}, different societies {cross:.3f}; permutation p={p7:.4f}")

# F8: own-half order vs shared order on temporally held-out split states
def nested_pred(order_rank, xs):
    """Given count of positives, predict the nested configuration from an order (lower rank = easier)."""
    s = xs.sum(1); pred = np.zeros_like(xs)
    idx = np.argsort(order_rank, kind="stable")
    for r, n in enumerate(s):
        pred[r, idx[:n]] = 1
    return pred


ex_own, ex_sh, hd_own, hd_sh = [], [], [], []
for k in names:
    T = len(X[k]); h = T // 2
    sh_rank = order_ranks(sum(hist(X[j]) for j in names if j != k))
    for tr_part, te_part in ((slice(0, h), slice(h, T)), (slice(h, T), slice(0, h))):
        xt = X[k][te_part]; s = xt.sum(1); xt = xt[(s > 0) & (s < 8)]
        if len(xt) == 0: continue
        own = order_ranks(hist(X[k][tr_part]))
        po, ps = nested_pred(own, xt), nested_pred(sh_rank, xt)
        ex_own.append((po == xt).all(1).mean()); ex_sh.append((ps == xt).all(1).mean())
        hd_own.append(np.abs(po - xt).sum(1).mean()); hd_sh.append(np.abs(ps - xt).sum(1).mean())
pr("F8", f"held-out split states, exact nested match: own-half order {np.mean(ex_own):.3f} vs shared {np.mean(ex_sh):.3f} "
         f"(p={p_w(ex_own, ex_sh):.3f}); Hamming {np.mean(hd_own):.3f} vs {np.mean(hd_sh):.3f} (p={p_w(hd_own, hd_sh):.3f})")

# ------------------------------------------------------------------ F9 recovery hierarchy
ev, fail = [], []
for k in names:
    Mk, Xk = M[k], X[k]; T = len(Mk)
    for t in range(1, T):
        if Mk[t - 1] <= 0 < Mk[t]:
            ev.append(Xk[t])
            e = t
            while e + 1 < T and Mk[e + 1] > 0: e += 1
            fail.append(~Xk[t:e + 1].any(0))
ev, fail = np.array(ev), np.array(fail)
o = np.argsort(-ev.mean(0))
pr("F9", f"{len(ev)} up-crossings; already positive: " + ", ".join(f"{NODES[i]} {ev.mean(0)[i]:.3f}" for i in o))
pr("F9", "never positive during the spell: " + ", ".join(f"{NODES[i]} {fail.mean(0)[i]:.3f}" for i in np.argsort(-fail.mean(0))))

# ------------------------------------------------------------------ F10 recovery/decline-window couplings
def windows(k, up=True):
    Mk = M[k]; w = np.zeros(len(Mk), bool)
    for t in range(1, len(Mk)):
        if (Mk[t - 1] <= 0 < Mk[t]) if up else (Mk[t - 1] > 0 >= Mk[t]):
            w[t:t + 3] = True
    return w


for up, lab in ((True, "recovery"), (False, "decline")):
    for q in (1, 3):
        def build(k, grammar, q=q):
            F, _ = REC[(k, q)]; Dk = D[k]; R = Dk - F; w = windows(k, up)[:-1, None]
            rows, ys = {}, {}
            for i in range(8):
                o = [j for j in range(8) if j != i]
                Z = [Dk[:-1, i], F[:-1, i], w[:, 0]]
                if grammar: Z += [R[:-1][:, o] * w]
                rows[i] = np.column_stack(Z); ys[i] = (Dk[1:, i], w[:, 0])
            return rows, ys
        errs = {}
        for g in (False, True):
            B = {k: build(k, g) for k in names}; per = {}
            for k in names:
                tr = [j for j in names if j != k]; se = n = 0.0
                for i in range(8):
                    mdl = Ridge(alpha=1e-3).fit(np.vstack([B[j][0][i] for j in tr]), np.concatenate([B[j][1][i][0] for j in tr]))
                    y, w = B[k][1][i]
                    if w.sum() == 0: continue
                    e = (y - mdl.predict(B[k][0][i]))[w]
                    se += ((e / D[k][:, i].std()) ** 2).sum(); n += len(e)
                if n: per[k] = np.sqrt(se / n)
            errs[g] = per
        ks = [k for k in errs[False] if k in errs[True]]
        a, b = np.array([errs[False][k] for k in ks]), np.array([errs[True][k] for k in ks])
        pr("F10", f"{lab:8s} windows k={q}: standardized RMSE no-grammar {a.mean():.4f} vs window-grammar {b.mean():.4f}; "
                  f"better in {(b < a).sum()}/{len(a)}, p={p_w(b, a):.3f}")

# ------------------------------------------------------------------ F11 innovation covariance
for q in (1, 3, 5):
    E = {}
    for k in names:
        if q == 5:
            mu, L, G, _ = factor_fit(D[k], 5); F = mu + G @ L.T
        else:
            F = REC[(k, q)][0]
        Dk = D[k]; e = np.zeros((len(Dk) - 1, 8))
        for i in range(8):
            Z = np.column_stack([np.ones(len(Dk) - 1), Dk[:-1, i], F[:-1, i]])
            beta = np.linalg.lstsq(Z, Dk[1:, i], rcond=None)[0]
            e[:, i] = Dk[1:, i] - Z @ beta
        E[k] = (e - e.mean(0)) / e.std(0)
    gain = []
    for k in names:
        C = np.corrcoef(np.vstack([E[j] for j in names if j != k]).T)
        sign, logdet = np.linalg.slogdet(C)
        Ci = np.linalg.inv(C); z = E[k]
        ll_full = -0.5 * (logdet + np.einsum("ti,ij,tj->t", z, Ci, z)).mean()
        ll_diag = -0.5 * (z ** 2).sum(1).mean()
        gain.append(ll_full - ll_diag)
    gain = np.array(gain)
    pr("F11", f"k={q}: held-out log-lik gain of shared innovation correlation vs independent (nats/step) {gain.mean():+.3f}; "
              f"positive in {(gain > 0).sum()}/{len(gain)}")

# ------------------------------------------------------------------ F12 lagged directed edges
S = np.zeros((len(names), 8, 8)); Tt = np.zeros_like(S)
for a_, k in enumerate(names):
    mu, L, G, _ = factor_fit(D[k], 3); R = D[k] - mu - G @ L.T
    for i in range(8):
        for j in range(8):
            if i == j: continue
            Z = np.column_stack([np.ones(len(R) - 1), R[:-1, i], R[:-1, j]])
            y = R[1:, i]
            beta, *_ = np.linalg.lstsq(Z, y, rcond=None)
            res = y - Z @ beta
            se = np.sqrt(res @ res / (len(y) - 3) * np.linalg.inv(Z.T @ Z)[2, 2])
            S[a_, i, j] = np.sign(beta[2]); Tt[a_, i, j] = beta[2] / se
edges = []
for i in range(8):
    for j in range(8):
        if i == j: continue
        npos = int((S[:, i, j] > 0).sum())
        p = binomtest(npos, len(names), 0.5).pvalue
        edges.append((f"{NODES[j]}->{NODES[i]}", npos, p))
ps = np.array([e[2] for e in edges]); o = np.argsort(ps); n = len(ps)
q = np.minimum.accumulate((ps[o] * n / np.arange(1, n + 1))[::-1])[::-1]
qv = np.empty(n); qv[o] = np.minimum(q, 1)
top = sorted(zip(edges, qv), key=lambda r: r[0][2])[:6]
pr("F12", "most consistent lagged edges (source->target, #societies positive, sign-test p, BH q): " +
   "; ".join(f"{e[0]} {e[1]}/{len(names)} p={e[2]:.4f} q={qq:.3f}" for e, qq in top))
for lab in ("Helm->Shield", "Shield->Helm"):
    e = [x for x in zip(edges, qv) if x[0][0] == lab][0]
    pr("F12", f"{lab}: positive in {e[0][1]}/{len(names)}, p={e[0][2]:.4f}, q={e[1]:.3f}")
pr("F12", f"edges with q<0.05: {int((qv < 0.05).sum())}/56")
