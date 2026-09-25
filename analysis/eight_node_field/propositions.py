"""Test propositions 1-22 on a CAMS ensemble panel.

Usage: python propositions.py [discovery|replication|all] [n_sims]
"""
import json, sys
from itertools import combinations
import numpy as np
from scipy.stats import spearmanr, rankdata, chi2
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score

from load import load, derive, NODES
from common import codes, eff_repertoire, switching, factor_fit, fit_dynamic_1f, simulate_1f

which = sys.argv[1] if len(sys.argv) > 1 else "discovery"
NS = int(sys.argv[2]) if len(sys.argv) > 2 else 200
rng = np.random.default_rng(20260925)

if which == "all":
    data = {**load("discovery"), **load("replication")}
else:
    data = load(which)
names = list(data)
S = {k: derive(a) for k, (y, a) in data.items()}
out = {"set": which, "societies": names, "n": len(names)}


def pr(tag, msg):
    print(f"[{tag}] {msg}")


def kendall_w(R):
    """R: m societies x n items ranks."""
    m, n = R.shape
    Rs = R.sum(0)
    return 12 * ((Rs - Rs.mean()) ** 2).sum() / (m ** 2 * (n ** 3 - n))


def w_perm(V, nperm=10000):
    R = np.array([rankdata(v) for v in V])
    w = kendall_w(R)
    null = np.array([kendall_w(np.array([rng.permutation(r) for r in R])) for _ in range(nperm)])
    return w, (1 + (null >= w).sum()) / (nperm + 1), R.mean(0)


# ---------- fit per-society models ----------
FIT = {k: fit_dynamic_1f(S[k]["D"]) for k in names}
SIMS = {k: [simulate_1f(FIT[k], len(S[k]["D"]), rng) for _ in range(NS)] for k in names}

# ---------- P1 repertoire vs independence / persistence-preserving null ----------
rows = []
for k in names:
    X = S[k]["X"]; T = len(X); obs = len(np.unique(codes(X)))
    ind = [len(np.unique(codes(np.array([rng.permutation(X[:, i]) for i in range(8)]).T))) for _ in range(NS)]
    circ = [len(np.unique(codes(np.array([np.roll(X[:, i], rng.integers(T)) for i in range(8)]).T))) for _ in range(NS)]
    f1 = [len(np.unique(codes((d > 0).astype(int)))) for d in SIMS[k]]
    rows.append((k, T, obs, np.mean(ind), np.mean(circ), np.mean(f1), np.mean(np.array(f1) <= obs)))
out["P1_P2_repertoire"] = rows
pr("P1", "society T observed | indep | circular-shift | 1-factor | P(1f<=obs)")
for r in rows:
    pr("P1", f"{r[0]:13s} {r[1]:4d} {r[2]:4d} | {r[3]:6.1f} | {r[4]:6.1f} | {r[5]:6.1f} | {r[6]:.3f}")
pr("P1", f"all below circular null: {all(r[2] < r[4] for r in rows)}")
pr("P2", f"mean PC1 share of K-S var = {np.mean([FIT[k]['share'][0] for k in names]):.3f}; "
         f"societies below 5% tail of 1f null = {sum(r[6] < 0.05 for r in rows)}/{len(rows)}")

# ---------- P3 low-rank sign reconstruction ----------
acc = {1: [], 2: [], 3: []}; exact = {1: [], 2: [], 3: []}
for k in names:
    D = S[k]["D"]; X = S[k]["X"]
    for q in (1, 2, 3):
        mu, L, G, _ = factor_fit(D, q)
        Xh = (mu + G @ L.T > 0).astype(int)
        acc[q].append((Xh == X).mean()); exact[q].append((Xh == X).all(1).mean())
out["P3"] = {q: (float(np.mean(acc[q])), float(np.mean(exact[q]))) for q in acc}
for q in acc:
    pr("P3", f"rank {q}: sign acc {np.mean(acc[q]):.3f}, exact config {np.mean(exact[q]):.3f}")

# ---------- P4/P5 node identity: Kendall W ----------
negfrac, lam, thr = [], [], []
for k in names:
    X = S[k]["X"]; split = (X.sum(1) > 0) & (X.sum(1) < 8)
    negfrac.append(1 - X[split].mean(0) if split.any() else np.full(8, .5))
    lam.append(FIT[k]["L"])
    G1 = factor_fit(S[k]["D"], 1)[2][:, 0]
    # threshold = level of the common field at which node crosses zero (in G SD units)
    thr.append(-FIT[k]["mu"] / np.where(np.abs(FIT[k]["L"]) > 1e-9, FIT[k]["L"], 1e-9))
for tag, V in (("neg-in-split", negfrac), ("loading", lam), ("threshold", thr)):
    w, p, mr = w_perm(np.array(V))
    order = [NODES[i] for i in np.argsort(-mr)]
    out[f"P4_{tag}"] = dict(W=w, p=p, order_high_to_low=order, mean_rank=dict(zip(NODES, mr.round(2).tolist())))
    pr("P4", f"{tag:12s} W={w:.3f} p={p:.4f}  high->low: {order}")
pr("P5", "median threshold (G SD) by node: " + ", ".join(f"{n}:{v:+.2f}" for n, v in zip(NODES, np.median(np.clip(thr, -10, 10), 0))))
pr("P5", "median loading by node: " + ", ".join(f"{n}:{v:.2f}" for n, v in zip(NODES, np.median(lam, 0))))
pr("P5", "median fraction positive by node: " + ", ".join(f"{n}:{v:.2f}" for n, v in zip(NODES, np.median([S[k]['X'].mean(0) for k in names], 0))))

# ---------- P6 pairwise sign correlations ----------
allpos, means, negpairs = True, [], []
for k in names:
    X = S[k]["X"].astype(float)
    sd = X.std(0)
    ok = sd > 0
    C = np.corrcoef(X[:, ok].T)
    iu = np.triu_indices(ok.sum(), 1)
    means.append(np.nanmean(C[iu]))
    nneg = int((C[iu] <= 0).sum()); negpairs.append((k, nneg, int(8 - ok.sum())))
out["P6"] = dict(mean_phi=float(np.mean(means)), nonpositive_pairs=negpairs)
pr("P6", f"mean pairwise sign corr {np.mean(means):.3f}; (society, #pairs<=0, #constant nodes): {negpairs}")

# ---------- P7/P8 cross-society motif overlap ----------
reps = {k: set(np.unique(codes(S[k]["X"]))) for k in names}
pc = np.array([bin(c).count("1") for c in range(256)])
by_k = {j: np.where(pc == j)[0] for j in range(9)}


def overlap_stat(R):
    return sum(len(R[a] & R[b]) for a, b in combinations(R, 2))


def jacc(R):
    return np.mean([len(R[a] & R[b]) / len(R[a] | R[b]) for a, b in combinations(R, 2)])


def split_only(R):
    return {k: {c for c in v if 0 < pc[c] < 8} for k, v in R.items()}


obs_ov = overlap_stat(split_only(reps))
null_ov = []
for _ in range(1000):
    R = {}
    for k, v in split_only(reps).items():
        cnt = np.bincount([pc[c] for c in v], minlength=9)
        R[k] = set(np.concatenate([rng.choice(by_k[j], cnt[j], replace=False) for j in range(9) if cnt[j]]))
    null_ov.append(overlap_stat(R))
null_ov = np.array(null_ov)
f1_ov, f1_j = [], []
for s in range(NS):
    R = split_only({k: set(np.unique(codes((SIMS[k][s] > 0).astype(int)))) for k in names})
    f1_ov.append(overlap_stat(R)); f1_j.append(jacc(R))
f1_ov = np.array(f1_ov)
out["P7"] = dict(obs=obs_ov, null_mean=float(null_ov.mean()), p=float((1 + (null_ov >= obs_ov).sum()) / 1001))
out["P8"] = dict(jaccard=float(jacc(split_only(reps))), f1_overlap_mean=float(f1_ov.mean()),
                 f1_overlap_p=float((1 + (f1_ov >= obs_ov).sum()) / (NS + 1)), f1_jaccard=float(np.mean(f1_j)))
pr("P7", f"split-config pairwise shared count obs={obs_ov}, stratified random null mean={null_ov.mean():.1f}, p={out['P7']['p']:.4f}")
pr("P8", f"mean Jaccard (split repertoires)={out['P8']['jaccard']:.3f}; 1-factor null overlap mean={f1_ov.mean():.1f} "
         f"(obs {obs_ov}), p(1f>=obs)={out['P8']['f1_overlap_p']:.3f}; 1f Jaccard={np.mean(f1_j):.3f}")
# per-config recurrence vs one-factor null
cnt_obs = np.zeros(256)
for v in split_only(reps).values():
    for c in v: cnt_obs[c] += 1
cnt_null = np.zeros((NS, 256))
for s in range(NS):
    for k in names:
        for c in set(np.unique(codes((SIMS[k][s] > 0).astype(int)))):
            if 0 < pc[c] < 8: cnt_null[s, c] += 1
pv = (1 + (cnt_null >= cnt_obs).sum(0)) / (NS + 1)
tested = np.where(cnt_obs >= 3)[0]
bonf = [(int(c), int(cnt_obs[c]), float(pv[c])) for c in tested if pv[c] * len(tested) < 0.05]
fmt = lambda c: "".join("+" if (c >> i) & 1 else "-" for i in range(8))
top = sorted(range(256), key=lambda c: -cnt_obs[c])[:6]
pr("P8", f"configs in >=5 societies: {int((cnt_obs >= 5).sum())}; top: " +
   "; ".join(f"{fmt(c)} n={int(cnt_obs[c])} null={cnt_null[:, c].mean():.1f}" for c in top))
pr("P8", f"configs surviving Bonferroni vs 1f null: {[(fmt(c), n, p) for c, n, p in bonf]}  (node order {NODES})")
out["P8"]["bonferroni_survivors"] = [(fmt(c), n, p) for c, n, p in bonf]

# ---------- P10/P11 morphology ----------
er = np.array([eff_repertoire(codes(S[k]["X"])) for k in names])
sw = np.array([switching(codes(S[k]["X"])) for k in names])
sh = np.array([FIT[k]["share"][0] for k in names])
out["P10"] = dict(eff_range=[float(er.min()), float(er.max())], sw_range=[float(sw.min()), float(sw.max())],
                  rho=float(spearmanr(er, sw)[0]))
out["P11"] = dict(rho_share_er=float(spearmanr(sh, er)[0]), rho_share_sw=float(spearmanr(sh, sw)[0]))
pr("P10", f"eff repertoire {er.min():.2f}-{er.max():.2f}; switching {sw.min():.3f}-{sw.max():.3f}; rho={out['P10']['rho']:.2f}")
pr("P11", f"rho(PC1 share, eff rep)={out['P11']['rho_share_er']:.2f}; rho(PC1 share, switching)={out['P11']['rho_share_sw']:.2f}")
# length confound
Ts = np.array([len(S[k]["X"]) for k in names])
pr("P11", f"length check: rho(T, eff rep)={spearmanr(Ts, er)[0]:.2f}")

# ---------- P12 path continuity ----------
ham, hsh, hf1 = [], [], []
for k in names:
    X = S[k]["X"]
    ham.append(np.abs(np.diff(X, axis=0)).sum(1).mean())
    hsh.append(np.mean([np.abs(np.diff(X[rng.permutation(len(X))], axis=0)).sum(1).mean() for _ in range(200)]))
    hf1.append(np.mean([np.abs(np.diff((d > 0).astype(int), axis=0)).sum(1).mean() for d in SIMS[k]]))
out["P12"] = dict(obs=float(np.mean(ham)), shuffled=float(np.mean(hsh)), f1=float(np.mean(hf1)),
                  all_below=bool(all(a < b for a, b in zip(ham, hsh))))
pr("P12", f"adjacent Hamming obs {np.mean(ham):.2f}, shuffled {np.mean(hsh):.2f} (range {min(hsh):.2f}-{max(hsh):.2f}), 1-factor {np.mean(hf1):.2f}")

# ---------- P13 AC changes the carrier ----------
pos_el = pos_ch = neg_el = neg_ch = 0
for k in names:
    D, sg = S[k]["D"], S[k]["sig"]
    for t in range(len(D)):
        p = D[t] > 0; n = D[t] < 0
        if p.sum() >= 2:
            pos_el += 1
            pos_ch += np.argmax(np.where(p, D[t], -np.inf)) != np.argmax(np.where(p, sg[t], -np.inf))
        if n.sum() >= 2:
            neg_el += 1
            neg_ch += np.argmin(np.where(n, D[t], np.inf)) != np.argmin(np.where(n, sg[t], np.inf))
out["P13"] = dict(pos=(int(pos_ch), pos_el), neg=(int(neg_ch), neg_el))
pr("P13", f"largest + carrier changes {pos_ch}/{pos_el}={pos_ch/pos_el:.3f}; largest - carrier {neg_ch}/{neg_el}={neg_ch/max(neg_el,1):.3f}")

# ---------- P14 AC lead/lag ----------
lags = range(-3, 4)
best = []; pooled = {n: {l: [] for l in lags} for n in NODES}
for k in names:
    dAC = np.diff(S[k]["AC"], axis=0); dD = np.diff(S[k]["D"], axis=0)
    for i, n in enumerate(NODES):
        rs = {}
        for l in lags:  # l>0: AC change at t precedes D change at t+l
            a = dAC[: len(dAC) - l, i] if l >= 0 else dAC[-l:, i]
            b = dD[l:, i] if l >= 0 else dD[: len(dD) + l, i]
            rs[l] = np.corrcoef(a, b)[0, 1] if a.std() > 0 and b.std() > 0 else np.nan
            pooled[n][l].append(rs[l])
        best.append(max(rs, key=lambda l: -1 if np.isnan(rs[l]) else abs(rs[l])))
best = np.array(best)
out["P14"] = dict(n=len(best), lag0=int((best == 0).sum()), pos=int((best > 0).sum()), neg=int((best < 0).sum()),
                  pooled={n: {l: float(np.nanmean(v)) for l, v in d.items()} for n, d in pooled.items()})
pr("P14", f"argmax|r| lag: 0 in {(best==0).sum()}/{len(best)}, AC-leads {(best>0).sum()}, D-leads {(best<0).sum()}")
pr("P14", "pooled r(lag0) / r(+1): " + ", ".join(f"{n}:{np.nanmean(pooled[n][0]):.2f}/{np.nanmean(pooled[n][1]):.2f}" for n in NODES))

# ---------- P15-P17 recovery ----------
def recoveries(D, M, AC=None):
    X = (D > 0).astype(int); T = len(M); res = []
    sdM = M.std()
    for t in range(1, T):
        if M[t - 1] <= 0 < M[t]:
            split = X[t].sum() < 8
            e = t
            while e + 1 < T and M[e + 1] > 0: e += 1
            full = X[t:e + 1].all(1).any()
            prev = [u for u in range(t - 1) if M[u] > 0 and M[u + 1] <= 0]
            if prev:
                u = prev[-1]
                h = int(np.abs(X[t] - X[u]).sum())
                res.append((split, full, h, abs(M[t] - M[u]) / sdM))
            else:
                res.append((split, full, None, None))
    return res


R = [r for k in names for r in recoveries(S[k]["D"], S[k]["M"])]
nsplit = sum(r[0] for r in R); split_nofull = sum(r[0] and not r[1] for r in R)
ided = [r for r in R if r[2] is not None]; diff = sum(r[2] > 0 for r in ided)
close = [r for r in ided if r[3] <= 0.1]
out["P15_17"] = dict(n=len(R), split=nsplit, split_no_full=split_nofull, identified=len(ided), different=diff,
                     median_h=float(np.median([r[2] for r in ided])) if ided else None,
                     close=len(close), close_diff=sum(r[2] > 0 for r in close))
pr("P15", f"M<=0 -> M>0 transitions {len(R)}; split {nsplit} ({nsplit/len(R):.3f}); of split, never reach all-8+ in spell {split_nofull} ({split_nofull/max(nsplit,1):.3f})")
pr("P16", f"identified prior regimes {len(ided)}; different config {diff} ({diff/max(len(ided),1):.3f}); median Hamming {out['P15_17']['median_h']}; |dM|<=0.1SD: {out['P15_17']['close_diff']}/{len(close)}")
sim_split, sim_diff = [], []
for s in range(NS):
    rr = [r for k in names for r in recoveries(SIMS[k][s], (S[k]["AC"] * SIMS[k][s]).mean(1))]
    if rr:
        sim_split.append(np.mean([r[0] for r in rr]))
        ii = [r for r in rr if r[2] is not None]
        if ii: sim_diff.append(np.mean([r[2] > 0 for r in ii]))
out["P17"] = dict(sim_split=float(np.mean(sim_split)), sim_diff=float(np.mean(sim_diff)))
pr("P17", f"1-factor sim (observed AC): split-recovery rate {np.mean(sim_split):.3f}, config mismatch {np.mean(sim_diff):.3f}")

# ---------- P18 many-to-one ----------
tot = dif = 0
for k in names:
    M, c = S[k]["M"], codes(S[k]["X"]); sd = M.std()
    for a in range(len(M)):
        for b in range(a + 3, len(M)):
            if abs(M[a] - M[b]) < 0.01 * sd:
                tot += 1; dif += c[a] != c[b]
out["P18"] = dict(pairs=tot, different=dif)
pr("P18", f"near-identical-M pairs {tot}; different configs {dif/max(tot,1):.3f}")

# ---------- P19 cancellation ----------
kap = np.concatenate([S[k]["kappa"] for k in names])
Iz = np.concatenate([(S[k]["I"] - S[k]["I"].mean()) / S[k]["I"].std() for k in names])
hi = kap >= 0.8
q75 = np.quantile(Iz, 0.75)
out["P19"] = dict(n_k08=int(hi.sum()), frac_in_topI=float((Iz[hi] >= q75).mean()) if hi.any() else None, N=len(kap))
pr("P19", f"obs with kappa>=0.8: {hi.sum()}/{len(kap)}; share of these in top intensity quartile: {out['P19']['frac_in_topI']}")

# ---------- P20 rising vs falling LOSO ----------
Xs, ys, gs = [], [], []
for k in names:
    M = S[k]["M"]; z = (M - M.mean()) / M.std()
    Xs.append(np.column_stack([z[1:], S[k]["X"][1:]])); ys.append((M[1:] > M[:-1]).astype(int)); gs += [k] * (len(M) - 1)
Xa, ya, ga = np.vstack(Xs), np.concatenate(ys), np.array(gs)
auc_m, auc_ms = [], []
for k in names:
    tr, te = ga != k, ga == k
    if len(set(ya[te])) < 2: continue
    for cols, store in (([0], auc_m), (list(range(9)), auc_ms)):
        m = LogisticRegression(max_iter=1000).fit(Xa[tr][:, cols], ya[tr])
        store.append(roc_auc_score(ya[te], m.predict_proba(Xa[te][:, cols])[:, 1]))
out["P20"] = dict(auc_M=float(np.mean(auc_m)), auc_M_signs=float(np.mean(auc_ms)))
pr("P20", f"LOSO AUC rising-vs-falling: M only {np.mean(auc_m):.3f}; M + 8 signs {np.mean(auc_ms):.3f}")

# ---------- P21 residual pair structure ----------
for q in (1, 2):
    sg = []
    for k in names:
        D = S[k]["D"]; mu, L, G, _ = factor_fit(D, q)
        Rr = D - mu - G @ L.T
        C = np.corrcoef(Rr.T); sg.append(np.sign(C[np.triu_indices(8, 1)]))
    sg = np.array(sg)
    cons = np.abs(sg.sum(0)) == len(names)
    pairs = [f"{NODES[i]}-{NODES[j]}({'+' if sg[0, n] > 0 else '-'})" for n, (i, j) in enumerate(zip(*np.triu_indices(8, 1))) if cons[n]]
    # binomial expectation of fully consistent pairs by chance
    exp = 28 * 2 * 0.5 ** len(names)
    out[f"P21_k{q}"] = pairs
    pr("P21", f"after removing {q} factor(s): pairs with same residual sign in all {len(names)} societies: {pairs} (chance ~{exp:.3f})")

json.dump(out, open(f"results_{which}.json", "w"), indent=1, default=lambda o: o.item() if hasattr(o, "item") else str(o))
