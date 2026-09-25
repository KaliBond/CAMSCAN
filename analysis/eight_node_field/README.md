# Eight-node field vs relational grammar: testing the 24 propositions

Tested 2026-09-25 against the CAMS ensemble-mean (five-scorer) panels in
[`KaliBond/wintermute`](https://github.com/KaliBond/wintermute).

## Data and definitions

| Set | Source | Societies (annual, ≥50 consecutive complete years) |
|---|---|---|
| **Discovery (11)** | `data/v2.3/canonical/*_ENS_*` | Argentina, Australia, Chile, China, France, Germany, Iran, Norway, Poland, Russia, UK |
| **Replication (20, unseen)** | `data/nations/*_ENS.csv` not in the canonical set | Brazil, Cambodia, Colombia, Greece, Hong Kong, Indonesia, Israel, Italy, Japan, Laos, Latvia, Lebanon, Mongolia, New Zealand, Philippines, Saudi Arabia, Singapore, South Africa, Ukraine, Venezuela |

Canada, Sweden and Thailand (5–10 year steps) and LatimVetus (20-year steps) are excluded because they aren't annual.
The original analysis used 12 societies, and which 12 isn't recorded, so the numbers below are a
re-derivation, not a reproduction.

* `D = K − S` for each node. A node is positive when `D > 0`. Ties (`K = S`) count as non-positive (a sensitivity check is in `robustness.py`).
* `AC = A·C/100`, and the node field is `σ = AC·D` (`cams_framework_v2_4.py`).
* `M = mean σ`, `I = mean |σ|`, `κ = 1 − |M|/I`.
* One-factor dynamic null: PCA on centred D, AR(1) common field, and independent AR(1) node residuals (200 simulations).

Scripts: `load.py`, `common.py`, `propositions.py`, `robustness.py`, `key_test.py`, `extras.py`, `society_fit.py`.
Raw outputs: `out_*.txt`, `results_*.json`. Everything runs in under 10 minutes.

## Verdicts

✅ supported · 🟡 supported with an important qualification · ❌ not supported in these data

| # | Proposition (short) | Discovery (11) | Replication (20) | Verdict |
|---|---|---|---|---|
| 1 | Repertoire far below independence and persistence nulls | 11/11 below the circular-shift null | 20/20 | ✅ |
| 2 | One-factor null reproduces the restriction (only 1/12 below the 5% tail) | PC1 = 81% of variance, but **7/11** below the 5% tail | PC1 = 85%, **10/20** below | ❌ The restriction is stronger than one Gaussian field produces in about half the societies. A field-preserving surrogate (observed G, residuals shifted per node) still leaves 6/11 and 10/20 below. |
| 3 | Rank-1/2/3 sign accuracy about 0.89/0.92/0.94, exact configuration about 0.52 | 0.880/0.899/0.917, exact 0.50 | 0.886/0.905/0.924, exact 0.54 | ✅ |
| 4 | Node rankings concordant across societies | W = 0.42 (neg-in-split), 0.51 (loading), 0.50 (threshold), all p = 1e-4 | W = 0.14, 0.30, 0.22, p ≤ 0.006 | ✅ (weaker out of sample) |
| 5 | Hands/Stewards/Flow/Helm harder to turn positive; Archive easiest | Hands > Helm > Flow > Stewards > Shield > Craft > Lore > Archive | Hands > Craft > Flow > Stewards > Helm > Shield > Archive > Lore | 🟡 **Hands** is the only node positive < 50% of the time in both sets (32%/32%). Lore/Archive sit at the easy end in both. The middle of the order moves (Shield, Craft). Discovery vs replication rank agreement ρ = 0.63. |
| 6 | Every pair positively sign-correlated in every society | mean φ = 0.50; exceptions in Argentina (3 pairs) and UK (1) | mean φ = 0.56; exceptions in Laos (1), Latvia (4), Ukraine (7) | 🟡 No antagonistic blocs, but "every pair, every society" is too strong |
| 7 | Exact-motif overlap above a stratified random null | 355 vs 181, p = 0.001 | 404 vs 262, p = 0.001 | ✅ |
| 8 | Motifs aren't the conserved object; no universal coalition survives correction | Jaccard 0.149 (one-factor null 0.162); 0 configurations survive Bonferroni | Jaccard 0.072 (null 0.092); 0 survive | ✅ Observed repertoires are, if anything, *more* idiosyncratic than the one-factor null |
| 9 | Conservation hierarchy: identities > coalitions > morphologies | — | — | ✅ Refined: **loadings are conserved better than thresholds** (see parameter sharing below) |
| 10 | Wide spread in repertoire and switching; strongly associated | eff 4.6–21.8, switching 0.19–0.60, ρ = 0.61 | eff 2.5–21.6, switching 0.06–0.61, ρ = 0.76 | ✅ (ρ < 0.91) |
| 11 | Synchronisation (PC1 share) is the main axis of morphology | ρ = −0.27 / −0.29; common 1950–2025 window: +0.12 / +0.26 | ρ = −0.87 / −0.66; window −0.74 / −0.75 | 🟡 Absent within the canonical 11. Across all 31, it survives controls for length and repeated rows (β = −3.6, p = 0.001), but it is carried by the PC1 > 0.9 societies, several of which are 70–89% repeated rows (see caveats) |
| 12 | Paths are continuous; the one-factor model reproduces most of this | adjacent Hamming 0.74 vs shuffled 2.64; one-factor 1.17 | 0.72 vs 3.05; one-factor 1.17 | 🟡 Continuity is real, but the data are *smoother* than the one-factor model. Repeated rows contribute. |
| 13 | AC changes the carrier, more on the negative side | + carrier 23.5%, − carrier 41.8% | 20.3%, 56.3% | ✅ (larger than originally reported) |
| 14 | AC does not lead K−S | largest correlation at lag 0 in 85/88 | 157/160 | ✅ Strong |
| 15 | Aggregate recovery ≠ structural recovery | 71 recoveries, 94.4% split; 67% never reach all-8-positive within the spell | 93 recoveries, 95.7%; 65% | ✅ |
| 16 | Recovery lands in a different configuration | 94.2%, median Hamming 3; 16/18 when \|ΔM\| ≤ 0.1 SD | 77.6%, median 2; 11/22 | 🟡 Holds in discovery, weaker in replication |
| 17 | The one-factor model reproduces the recovery mismatch | simulated split 0.90, mismatch 0.92 | 0.92, 0.84 | ✅ |
| 18 | M is many-to-one (~58%) | 34% | 8% (inflated by identical repeated years) | 🟡 Qualitatively true; the magnitude doesn't replicate |
| 19 | Cancellation occurs, but isn't a high-intensity phenomenon | 66 observations with κ ≥ 0.8; 3% of them in the top intensity quartile | 120; 0.8% | ✅ Stronger than stated: cancellation is essentially a **low-intensity** state |
| 20 | Configuration doesn't help predict rising vs falling | leave-one-society-out AUC: M 0.568, M + signs 0.558 | 0.565 vs 0.560 | ✅ |
| 21 | No universal pairwise residual channels | after 1 factor: Helm–Flow(−) and Shield–Craft(−) consistent in 11/11; none after 2 factors; Lore–Archive not reproduced | none after 1 or 2 factors | ✅ |

## The decisive test: shared field vs relational grammar (`key_test.py`)

Every model was fitted on training societies and scored on societies it never saw.

**Contemporaneous (predict each node from the other seven).** Each model is a covariance structure,
and the prediction is its Gaussian conditional mean.

| Held-out design | FA1 (16 params) | FA2 | FA3 | FULL relations (36 params) | FULL beats FA1 in |
|---|---|---|---|---|---|
| Leave one society out, discovery (11) | R² 0.727 | 0.736 | 0.745 | 0.752 | 9/11 |
| Train discovery → test replication (20) | 0.768 | 0.791 | 0.800 | 0.806 | 17/20 |
| Leave one society out, all 31 | 0.751 | 0.767 | 0.775 | 0.782 | 20/31 |

Sign level: pairwise Ising couplings beat a loading-weighted mean-field model on held-out log-loss in **23/31** societies.

**Dynamic (predict the one-step change ΔD).** Out-of-sample R² is only about 6%, and almost all of it
comes from each node's own mean reversion. Adding the common field adds at most 0.7 points (and nothing
on the replication set). A full cross-lagged node→node VAR **never** beats the field model, and adding
cross-lagged relations among field residuals adds nothing.

**Parameter sharing (last 40% of each held-out society; median test MSE relative to the society's own FA1).**

| Parameters | MSE ratio | Sign accuracy |
|---|---|---|
| Own loadings and means (FA1) | 1.00 | 0.840 |
| **Loadings from other societies**, own means | 1.02 (better in 15/31) | 0.838 |
| Shared loadings **and shared node offsets** | 1.15 | 0.828 |
| Own full covariance | 0.83 | 0.840 |
| Shared full covariance | 0.91 | 0.848 |

### What this establishes

1. **Loadings are conserved and thresholds are partly local.** Importing the node loadings from other
   societies costs almost nothing. Importing the node offsets (thresholds) costs about 15% in MSE.
2. **There is reproducible structure beyond one field.** It appears in held-out societies, both
   continuous and sign-level, and part of it is shared across societies.
3. **That structure is contemporaneous, not causal.** No lagged node-to-node relation predicts out of sample.
4. **Eight nodes cannot tell "a few fields" from "a relational grammar".** A three-factor model on eight
   variables has only 7 fewer free parameters than the unconstrained covariance, and it captures most of the
   gain. Covariance data alone therefore cannot separate three conserved fields from a pairwise grammar.
   Separating them needs interventions or exogenous shocks that hit one node, not more nulls.

A model consistent with everything here is:

  D(i,s,t) = μ(i,s) + λ(i)·G(s,t) + [low-rank or relational term, partly shared] + ε(i,s,t)

The loading λ is **shared across societies** (unlike the original box, where it varied by society),
the threshold μ is society-specific, and there are no lagged cross-node terms.

## Blind panel (JUNO_Blind_motif, 12 societies × 61 steps)

This is the panel the propositions were originally derived on. Society identities are withheld, and no
attempt was made to recover them. The data aren't committed; set `BLIND` / `BLIND_ENV` to the CSV paths.
Outputs: `out_blind.txt`, `out_robust_blind.txt`, `out_key_test_blind.txt`, `out_blind_extras.txt`.

**Reproduction.** The battery reproduces the original figures almost exactly: PC1 82.9%; rank-1/2/3
sign accuracy 0.887/0.920/0.936; W = 0.33–0.38 (p = 1e-4); effective repertoire 4.96–20.40; switching
0.283–0.700; ρ = 0.91, −0.83, −0.74; Hamming 1.07; +/− carrier changes 101 and 156; lag 0 in 88/96;
44 recoveries, 40 split, 19/40 terminate early; 34/36 and 14/14 different configurations; 58%
many-to-one; 28 cases of κ ≥ 0.8; AUC 0.628 vs 0.621; Lore–Archive(+) residual in 12/12.
Every pair is positively sign-correlated in all 12 societies.

This panel has almost no repeated rows (except Society_K, 20%) and few ties (3%). Every series is 61
steps, so length can't confound P11. That is why P2, P6, P11 and P18 hold here but weaken on the
wintermute corpora. The differences between panels come from the data, not the method.

- **P2:** 3/12 below the 5% tail of the Gaussian one-factor null (original 1/12). With the
  field-preserving surrogate, 1/12. ✅ on this panel.

**Decisive test on the blind panel**

| Design | FA1 | FA2 | FA3 | FULL | Ising vs mean-field (log-loss gain) |
|---|---|---|---|---|---|
| Leave one society out, blind | R² 0.753 | 0.777 | **0.787** | 0.783 | −24 (Ising worse) |
| Train on 31 corpus societies → blind | 0.751 | 0.758 | **0.786** | **0.786** | −45 (Ising worse) |

- **Contemporaneous:** the structure beyond one field is fully captured by 2–3 shared factors.
  Unconstrained pairwise relations add nothing. Pairwise sign couplings *lose* to a mean-field model out
  of sample (and beat a loading-weighted one in only 7/12).
- **Dynamic:** the common field adds 1.3 points of out-of-sample ΔD R² over own-node AR (0.088 → 0.101).
  Cross-lagged relations among field residuals add nothing (0.1005 → 0.1013).
- **Parameter sharing** (median MSE relative to the society's own FA1):

  | Shared parameters | from other blind societies | from the 31 corpus societies |
  |---|---|---|
  | loadings (FA1) | **0.94** (better in 8/12) | **0.99** (7/12) |
  | loadings + node offsets | 1.12 | 1.18 |
  | full covariance | **0.86** (10/12) | **0.85** (10/12) |
  | own full covariance | 1.13 (overfits) | — |

  Loadings are conserved, and they transfer *across corpora*. Thresholds are society-specific. The shared
  second-order structure (the 2–3 factors) transfers too.
- **Node order across corpora:** blind vs corpus-31 mean ranks agree at ρ = 0.82 (threshold), 0.83
  (loading) and 0.83 (fraction positive). Some blind societies may also appear in the corpus, so this isn't
  a fully independent replication.
- **Measurement floor:** residual variance is 3.4× the scorer-mean error variance after one factor, but
  only **2.0×** after two (≤ 1.3× in 5/12 societies). Anything beyond two fields is close to what five
  scorers can resolve. Scorer disagreement isn't higher near the K = S threshold (7/12, 0.656 vs 0.669).

**Verdict on the central question (blind panel).** The data support *conserved node response functions
(shared loadings) to two or three shared latent fields, with society-specific thresholds*. They do not
support an additional relational grammar. Held out, explicit node-to-node relations (same-time or lagged,
continuous or sign-level) never beat the low-rank field model, and the residual they would have to explain
sits at about twice scorer noise.

## Caveats the data raise

* **Repeated rows.** In many series a large share of year-to-year transitions repeat all 32 scores exactly:
  Italy 80%, Israel 89%, Lebanon 88%, Indonesia 87%, Hong Kong 69%, and France, Norway, Poland and Australia
  34–36%. This looks like carry-forward in the scoring runs. It inflates persistence, PC1 share and
  many-to-one counts, and restricts repertoires. P2, P11, P12, P16 and P18 are the most exposed.
* **Ties.** 9% of discovery K−S values (up to 30%; up to 59% in replication) are exact ties from the
  integer scoring grid. A Gaussian null can't reproduce this, which is part of why P2 fails.
* **Halo versus society.** Residuals after the field are 6–7× the five-scorer error variance (median), so
  they aren't scorer noise. But inter-scorer SD can't detect a *shared* prior. If all LLM scorers carry the
  same "era mood", that would itself look like a strong common field. The robust node order (Hands hardest,
  Lore and Archive easiest) could likewise be a universal scoring prior rather than a universal social
  structure. The discriminating test is to regress G and the node offsets on independent external series
  (conflict, fiscal, demographic, V-Dem-type indices).
* **No bloc structure.** Each society's agreement with the consensus node order (`out_society_fit.txt`)
  shows no geopolitical grouping. China (ρ = 0.81) sits beside Germany, the UK and Japan. Russia and Iran
  sit mid-range. The four weakest fits are Argentina, Saudi Arabia, Hong Kong and Venezuela. The shared
  functional basis crosses every political system in the panel.
