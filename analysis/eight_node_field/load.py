"""Load CAMS ensemble-mean (ENS) panels into an annual society x year x node array.

Data source: KaliBond/wintermute (cloned to $WINTERMUTE, default /home/user/kalibond/wintermute).
Discovery set  = v2.3 canonical ENS files with annual resolution.
Replication set = data/nations/*_ENS.csv societies NOT in the discovery set, annual resolution.
"""
import glob, os, re
import numpy as np, pandas as pd

WM = os.environ.get("WINTERMUTE", "/home/user/kalibond/wintermute")
NODES = ["Helm", "Shield", "Lore", "Stewards", "Craft", "Hands", "Archive", "Flow"]
DIMS = ["Coherence", "Capacity", "Stress", "Abstraction"]
MIN_YEARS = 50


def _read(path):
    df = pd.read_csv(path, dtype=str)
    df.columns = [c.strip() for c in df.columns]
    if "Society" not in df.columns:
        for alt in ("Nation", "Country", "society"):
            if alt in df.columns:
                df = df.rename(columns={alt: "Society"})
    df["Node"] = df["Node"].astype(str).str.strip()
    df = df[df["Node"].isin(NODES)].copy()
    df["Year"] = pd.to_numeric(df["Year"], errors="coerce")
    for d in DIMS:
        df[d] = pd.to_numeric(df[d], errors="coerce")
    df = df.dropna(subset=["Year"] + DIMS)
    df["Year"] = df["Year"].astype(int)
    df = df.groupby(["Year", "Node"], as_index=False)[DIMS].mean()
    return df


def _annual_block(df):
    """Longest run of consecutive years with all 8 nodes present."""
    full = df.groupby("Year")["Node"].nunique()
    yrs = np.array(sorted(full[full == 8].index))
    if len(yrs) == 0:
        return None
    runs, start = [], 0
    for k in range(1, len(yrs) + 1):
        if k == len(yrs) or yrs[k] != yrs[k - 1] + 1:
            runs.append(yrs[start:k]); start = k
    best = max(runs, key=len)
    if len(best) < MIN_YEARS:
        return None
    sub = df[df["Year"].isin(best)]
    arr = np.stack([sub.pivot(index="Year", columns="Node", values=d).loc[best, NODES].values for d in DIMS], -1)
    return best, arr  # arr: T x 8 x 4 (C,K,S,A)


def load(which="discovery"):
    canon = sorted(glob.glob(f"{WM}/data/v2.3/canonical/*_ENS_*_cleaned.csv"))
    canon = [p for p in canon if "LatimVetus" not in p]
    canon_names = {os.path.basename(p).split("_ENS_")[0] for p in canon}
    if which == "discovery":
        files = {os.path.basename(p).split("_ENS_")[0]: p for p in canon}
    else:
        files = {}
        for p in sorted(glob.glob(f"{WM}/data/nations/*_ENS.csv")):
            name = os.path.basename(p).replace("_ENS.csv", "")
            if name not in canon_names:
                files[name] = p
    out = {}
    for name, p in files.items():
        try:
            blk = _annual_block(_read(p))
        except Exception as e:  # malformed file
            print(f"skip {name}: {e}"); continue
        if blk is not None:
            out[name] = blk
    return out


def derive(arr):
    C, K, S, A = (arr[..., j] for j in range(4))
    D = K - S                      # node deficit/surplus
    AC = A * C / 100.0             # expression / amplitude term
    sig = AC * D                   # v2.4 node field sigma_i
    M = sig.mean(1)                # aggregate field
    I = np.abs(sig).mean(1)        # intensity
    kappa = np.where(I > 0, 1 - np.abs(M) / np.where(I > 0, I, 1), 0.0)  # cancellation
    X = (D > 0).astype(int)        # sign configuration (positive vs non-positive)
    return dict(D=D, AC=AC, sig=sig, M=M, I=I, kappa=kappa, X=X)


if __name__ == "__main__":
    for w in ("discovery", "replication"):
        d = load(w)
        print(w, len(d))
        for k, (yrs, a) in d.items():
            print(f"  {k:14s} {yrs[0]}-{yrs[-1]} T={len(yrs)}")
