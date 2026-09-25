# Eight-node field: propositions, redrafted after testing

Evidence base (see `README.md` and the `out_*.txt` files):
**B** = JUNO blind motif panel (12 societies × 61 steps; the panel the propositions were derived on);
**C** = canonical v2.3 ensembles (11 annual societies);
**R** = other `data/nations` ensembles (20 annual societies, unseen by any model fitted on B or C).

Status key: **Stands** · **Stands, revised** (numbers or wording changed) · **Provisional** (panel-dependent)

| # | Status | Main change |
|---|---|---|
| 1 | Stands | Holds in 43/43 societies |
| 2 | Stands, revised | 3/12 below the one-factor null on B (not 1/12); about half of C and R |
| 3 | Stands | Replicates in C and R |
| 4 | Stands | Weaker but significant in R |
| 5 | Stands, revised | Hands and Archive/Lore robust; Shield/Craft positions don't replicate |
| 6 | Stands, revised | "Every pair, every society" holds on B only (fails in 5/31 in C+R) |
| 7 | Stands | p = 0.001 on all three panels |
| 8 | Stands, revised | The null matches overlap in Jaccard, not in raw counts |
| 9 | Stands, refined | Loadings are the most conserved layer |
| 10 | Stands | ρ weaker in C/R (0.61, 0.76) |
| 11 | Provisional | Absent in C; confounded by repeated rows in R |
| 12 | Stands, revised | One-factor model reproduces about ¾ of the continuity, not "most" |
| 13 | Stands, revised | Eligibility-dependent numbers; negative side more affected in every panel |
| 14 | Stands | Lag 0 in 330/344 society × node cases |
| 15 | Stands | Reproduced exactly; replicates |
| 16 | Stands, revised | Weaker in R |
| 17 | Stands, revised | Simulated rates 0.83–0.92, slightly below observed |
| 18 | Provisional | 58% on B, 34% on C, 8% on R |
| 19 | Stands, strengthened | Cancellation is a low-intensity phenomenon |
| 20 | Stands | Replicates in C and R |
| 21 | Stands, strengthened | The surviving pairs differ between panels |
| 22 | Stands, sharpened | Loading conservation now shown by out-of-sample transfer |
| 23 | Revised | λ shared across societies; 2–3 fields; no lagged cross-node terms |
| 24 | Revised | The decisive test has been run |
| 25–29 | New | Decisive test, measurement floor, data conditions, halo caveat, no bloc structure |

---

## Redrafted conclusions

**1. The eight nodes are strongly coupled, not independent indicators.**
Repertoires are far smaller than under node-wise independence, and smaller than under a null that keeps
each node's own persistence (circular shifts). This holds in every society tested (12/12 B, 11/11 C,
20/20 R). This establishes a constrained joint state space.

**2. A restricted repertoire is not by itself evidence of an eight-node relational grammar.**
A one-factor dynamic model fitted to each society explains 82.9% of the variance in K−S (B; 81% C,
85% R). On B, only 3/12 societies fall below the 5% tail of the Gaussian one-factor repertoire null. Only
1/12 falls below a stricter surrogate that keeps the observed common field and shifts each node's
residual independently. On the larger corpora the restriction is stronger than one field produces in
about half the societies (surrogate: 6/11 C, 10/20 R). Part of that excess comes from integer-grid ties
(K = S) and carried-forward rows. Coupling is established; its mechanism is not.

**3. The minimal description is a common field plus reusable node-specific response functions.**
A rank-1 reconstruction predicts 88.7% of node signs and 51.9% of complete configurations exactly. Ranks
2 and 3 raise sign accuracy to 92.0% and 93.6% (B). The same pattern replicates in C (88.0/89.9/91.7%)
and R (88.6/90.5/92.4%). The system is strongly low-dimensional despite eight distinguishable functions.

**4. Node identity matters reproducibly.**
Node rankings are concordant across societies (node-label permutation, p = 10⁻⁴):
- tendency to stay negative in split configurations: W = 0.33;
- loading strength: W = 0.38;
- sign-crossing threshold: W = 0.38.

Concordance is stronger in C (W = 0.42–0.51) and weaker but still significant in the unseen R
(W = 0.14–0.30, p ≤ 0.006). Node rankings from B and from the 31-society corpus agree at ρ ≈ 0.82–0.83.
The node labels are not interchangeable.

**5. Recurrent roles are relative thresholds on a shared field, not coalitions.**
The most robust feature is that **Hands** is the hardest node to turn positive in every panel. It is the
only node positive less than half the time in both C and R (32%), and 50% in B. **Archive and Lore** sit
at the easy-to-positive end in every panel. Helm, Stewards and Flow usually sit on the harder side. The
positions of Shield and Craft shift between panels, so they shouldn't be assigned fixed roles. The
recurring order matters, not a four-versus-four partition.

**6. There are no antagonistic blocs.**
The mean pairwise sign correlation is 0.57 (B; 0.50 C, 0.56 R). On B every one of the 28 pairs is
positive in every society. Across the corpora, 5/31 societies have a few zero or negative pairs,
concentrated in short or turbulent series (Ukraine 7 pairs, Latvia 4, Argentina 3). The general picture
is nodes co-moving with a common field and crossing their thresholds at different points.

**7. Exact configurations recur across societies more than arbitrary states would.**
Against a null that keeps each society's repertoire size and its count of configurations with 1…7
positive nodes, cross-society overlap is well above chance: p = 0.001 on B (179 vs 109), C and R.

**8. Exact motifs are not the conserved object.**
Mean pairwise Jaccard overlap of split-state repertoires is 0.095 (B; 0.149 C, 0.072 R). The one-factor
null gives comparable or slightly *higher* Jaccard overlap (0.111 B). Raw overlap counts are higher still
under the null, because it generates larger repertoires. No individual configuration survives Bonferroni
correction against the one-factor null in any panel. There is no evidence for a universal coalition.

**9. Hierarchy of conservation.**
From most to least conserved:
1. **node loadings** (transferable across societies and corpora; see 25);
2. node order and relative thresholds;
3. absolute thresholds (society-specific);
4. exact coalitions;
5. whole-society repertoires and trajectories.

This fits "conserved local response functions + variable global morphology", not "one universal network".

**10. Societies differ widely in repertoire geometry.**
Effective repertoire runs from 4.96 to 20.40 states and switching probability from 0.283 to 0.700.
The two are strongly associated: ρ = 0.91 on B, 0.61 C, 0.76 R. Societies that explore more states move
between them faster.

**11. (Provisional) Synchronisation strength as an axis of morphology.**
On B, PC1 share correlates with effective repertoire at ρ = −0.83 and with switching at ρ = −0.74, with
series length held equal. The relation replicates in R (ρ = −0.87) and survives pooled controls for
length and repeated rows (β = −3.6, p = 0.001, n = 31). However, it is **absent within C** (ρ = −0.27;
+0.12 on a common 1950–2025 window). In R it is carried by series with PC1 > 0.9, which are mostly
carry-forward-heavy. Treat it as a candidate definition of global morphology, not an established one.

**12. Transitions are path-continuous, and mostly, but not entirely, explained by the field.**
About 1.07 nodes change sign per step, against 3.26 under random reordering (B; 0.74 vs 2.64–3.05 in
C/R). Every society is far below its shuffled null. The one-factor AR model gives 1.52. It accounts for
about three-quarters of the reduction, and the observed paths are smoother than the model. Continuity
does not imply pairwise causation.

**13. AC acts mainly as an expression/amplitude modifier, especially on the negative side.**
Weighting by AC changes the largest positive carrier in 18.5% of eligible observations and the largest
negative carrier in 41.5% (B, with eligibility defined as ≥ 2 same-sign nodes; the original analysis
reported 16.4% and 33.1% under its own eligibility rule). The same holds in C (23.5% / 41.8%) and R
(20.3% / 56.3%). AC matters much more for which node dominates the negative side.

**14. AC does not temporally lead K−S.**
The largest |ΔAC–Δ(K−S)| correlation is contemporaneous in 88/96 society × node cases (B), 85/88 (C)
and 157/160 (R). Pooled by node, every node peaks at lag 0 (r ≈ 0.54–0.70 on B), and one-step-ahead
correlations are near zero. "Scale/intensity coupled to the current state" is the safe reading.

**15. Aggregate and structural recovery are different observables.**
There are 44 transitions from M ≤ 0 to M > 0. In 40/44 (90.9%), aggregate positivity returns while at
least one node stays non-positive, and 19/40 of those spells end before all eight nodes are positive.
This replicates at scale: C has 71 recoveries (94.4% split, 67% ending early) and R has 93 (95.7%, 65%).

**16. Recovery usually does not recreate the previous configuration.**
On B, 34/36 recoveries with an identifiable prior positive regime return in a different configuration
(median difference 3 nodes), including 14/14 where M differs by ≤ 0.1 SD. C matches (94.2%; 16/18).
The effect is weaker in R (77.6%; 11/22), whose series often repeat identical rows.

**17. This non-equivalence does not require a special hysteretic mechanism.**
A one-factor dynamic simulation driven by the observed AC gives a split-recovery rate of 0.83–0.92 and a
configuration-mismatch rate of 0.84–0.92 across panels. These are close to, though slightly below, the
observed rates. The claim "aggregate recovery ≠ configuration recovery" is sound. Thresholded
low-dimensional dynamics is a sufficient mechanism.

**18. (Provisional in magnitude) M is many-to-one.**
Among same-society observations at least 3 steps apart with |ΔM| < 0.01 SD, 58% have different
configurations on B. The figure is 34% in C and only 8% in R, where near-identical M values often come
from carried-forward identical years. The qualitative claim holds, but the rate depends on data quality.

**19. M, I, κ and configuration carry non-redundant information, and cancellation is a low-intensity
state.**
There are 28 observations with κ ≥ 0.8 on B (66 in C, 120 in R). Only 3.6% of them are in the top
intensity quartile (3.0% C, 0.8% R). Strong signed cancellation happens mostly when the field is weak,
not when it is intense.

**20. Structural path dependence, not classical hysteresis.**
In leave-one-society-out prediction of rising vs falling, M alone gives AUC 0.628, and adding eight node
signs gives 0.621 (B; C 0.568 vs 0.558; R 0.565 vs 0.560). There is no society-general rising/falling
configuration signature.

**21. No universal pairwise channels.**
After one factor is removed, a few residual pairs keep the same sign across a whole panel, but they
differ between panels:
- B: Lore–Archive (+) and Craft–Archive (−);
- C: Helm–Flow (−) and Shield–Craft (−);
- R: none.

After two factors are removed, none remains in any panel. This is additional low-rank structure, not
invariant node-pair causation.

**22. The defensible "eight conserved functions" claim.**
Eight repeatedly identifiable response channels exist. Their loadings are conserved well enough to
transfer between societies and between independently scored corpora. Their thresholds are partly
society-specific. Societies differ in coupling, repertoire breadth, persistence and trajectory. Eight
independently interacting causal degrees of freedom are not demonstrated, and on the cleanest panel
they are actively disfavoured (see 25).

**23. Revised model.**

  D(i,s,t) = μ(i,s) + Σₖ λ(i,k)·G(k,s,t) + ε(i,s,t),  with k = 1…2 (at most 3)

- The G(k,s,t) are society-level common fields: one dominant field and one or two secondary ones.
- The loadings **λ(i,k) are shared across societies**. Sharing costs nothing out of sample.
- The thresholds **μ(i,s) are society-specific** but ordered: Hands high, Archive/Lore low.
- There are **no lagged node-to-node terms**.
- ε is close to scorer noise once two fields are removed (see 26).

**24. Central conclusion.**

> The same eight functional channels, with conserved response loadings, recur across societies within a
> field of two or three dimensions. Societies differ in the thresholds, coupling strength, repertoire and
> trajectories through which those channels are expressed. The data favour *conserved response
> functions* over *conserved interacting relations*.

## New propositions from the decisive test

**25. Shared-field vs relational-grammar test (held-out societies).**
- **Same-time structure, B:** predicting each node from the other seven, R² is 0.753 (1 factor),
  0.777 (2), **0.787** (3) and 0.783 (unconstrained pairwise relations). Trained on the 31 corpus
  societies and tested on B, three factors and full relations tie at 0.786.
- **Sign level, B:** pairwise (Ising) couplings do *worse* than mean-field on held-out log-loss (−24 and
  −45 nats).
- **Across the corpora:** full relations beat one factor in 20/31 societies (R² 0.751 → 0.782), but three
  factors recover most of that (0.775). Ising does beat weighted mean-field there (23/31), possibly
  because of grid ties.
- **Dynamics:** the field improves one-step prediction over own-node AR (0.088 → 0.101 on B).
  Cross-lagged node → node terms add nothing in any panel.
- **Parameter sharing:** loadings borrowed from other societies predict as well as or better than a
  society's own (median error ratio 0.94 within B, 0.99 imported from the corpus, 1.02 within the
  corpus). Borrowing node thresholds costs 12–18%. A shared full covariance beats a society's own
  one-factor model (0.85–0.91).

**26. Measurement floor.**
On B, residual variance is 3.4× the five-scorer error variance after one field, but only **2.0×** after
two (≤ 1.3× in 5/12 societies). Anything a relational grammar would have to explain beyond two fields
is within a factor of about 2 of what the ensemble can resolve. Scorer disagreement is not higher near
the K = S threshold.

**27. Data conditions for replication.**
The propositions reproduce fully on B, which has 3% ties, almost no repeated rows and equal-length
series. They weaken in proportion to data defects on the corpora: up to 30–59% ties, and 34–89% of
year-to-year transitions repeating all 32 scores exactly (Italy, Israel, Lebanon, Indonesia, Hong Kong;
about a third in France, Norway, Poland and Australia). Propositions 2, 6, 11, 16 and 18 are the most
sensitive. Future panels should record whether a row was scored or carried forward.

**28. Halo caveat (untested).**
Inter-scorer disagreement cannot detect a prior shared by all five LLM scorers. A shared "era mood"
would look like a common field, and a shared belief about institutions (for example, "labour is always
stressed; archives are always capable") would look like a conserved node order. The discriminating test
is to regress G and μ on independent external series: conflict, fiscal, demographic and
institutional indices.

**29. The shared basis crosses political systems.**
Agreement between each society's threshold order and the consensus of the others shows no geopolitical
grouping. China (ρ = 0.81) sits beside Germany, the UK and Japan, and Russia and Iran sit mid-range. The
weakest fits (Argentina, Saudi Arabia, Hong Kong, Venezuela) belong to no single bloc. The eight-channel
response basis behaves as a common structure of societies, not a property of any one political family.
