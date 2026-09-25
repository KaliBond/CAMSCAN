# Synthesis: a functional coordinate system, not a causal circuit

This is the final interpretation after two rounds of testing. The revised wording follows the reviewer's
edits, and the round-2 validation (`validate_round2.py`, `f7_check.py`, `f12_check.py`) adds several
corrections. **B** = JUNO blind panel (12 societies). **W** = wintermute corpus (31 societies, independently
scored).

"Conserved" and "portable" mean recurrence and stability across societies and panels. They do not mean
mathematical invariance.

## Round-2 validation: what replicated

| # | Finding | B | W | Status |
|---|---|---|---|---|
| 1 | Conditional grammar (regime- or sign-switched cross-node coefficients) doesn't beat a 3-factor baseline | worse (p = 0.027, 0.009) | neutral or worse (sign-switching p = 0.009) | **Replicates** |
| 2 | Nonlinear (gradient-boosted) grammar doesn't help | other nodes worse in 10/12, p = 0.016 | neutral, 16/31, p = 0.43 | **Replicates** (neutral in W) |
| 3 | Other nodes' current signs add nothing to X(t+1) once field and own sign are known | worse in 10/12, p = 0.002 | **slightly better** in 21/31, p = 0.012 (log-loss 0.2778 → 0.2730) | **B only.** In W, other signs add a small amount |
| 4 | Node-local memory dominates | own lag improves all 12, p = 4.9e-4 (0.430 → 0.307) | 30/31, p = 2e-9 (0.507 → 0.259) | **Replicates strongly** |
| 4b | The field's direction adds nothing | **adds** in 12/12 (0.307 → 0.300) | **adds** in 25/31 (0.259 → 0.241) | **Not supported.** Direction adds a small, consistent amount on top of memory |
| 5 | A shared order learned elsewhere fits held-out societies better than random orders | p = 0.0005 | p = 0.0005 | **Replicates** |
| 5b | Order: Shield > Lore > Archive > Craft > Helm > Flow > Stewards > Hands | exact chain optimum on B's split states (5 orders within 0.02) | optimum is Lore > Archive > Craft > Flow > Stewards > Helm > **Shield** > Hands | **Partly.** Hands last and Lore/Archive near the start are robust; **Shield swings from first to seventh** |
| 6 | Among split states, ~14% exactly on the chain and ~71% within one flip | 13.3%, 71.2%; distance 1.27 vs 1.62 random | 37.5%, 70.8%; 1.01 vs 1.55 | **Replicates** (within one flip ≈ 71% in both) |
| 7 | Society-specific order is persistent (split-half) | τ 0.425 vs 0.126, p = 0.0015 (0.420 vs 0.211, p = 0.019 with ≥ 10 split states per half) | τ 0.399 vs 0.138, p = 0.0005 | **Replicates** |
| 8 | Own-society order predicts held-out split states better than the shared order | 27.9% vs 11.8%, p = 0.019; with ≥ 10 split states per half, 28.3% vs 13.9%, p = 0.099 | 37.9% vs 37.6%, p = 0.75 | **Fragile.** Significant only on B's full sample; absent in W |
| 9 | Recovery follows the threshold hierarchy | Archive 82% … Hands 39% already positive | Helm, Lore, Archive, Shield ≈ 63–67% … Stewards 49%, **Hands 24%** | **Replicates** at the hard end (Stewards, Hands); the easy end reorders |
| 10 | No recovery-specific couplings | window grammar worse (p = 0.003, 0.012) | worse (p < 0.001) | **Replicates**; the couplings actively hurt prediction |
| 11 | Residual innovations share a correlation structure | +3.0 nats/step, 12/12, at 1, 3 and 5 factors | +3.5, 30/31 | **Replicates** |
| 12 | Helm→Shield / Shield→Helm as candidate lagged motifs | 8/12 and 3/12 positive; 0/56 edges at q < 0.05 | 14/31 each; 0/56 | **Motif not reproduced.** "No edge survives FDR" replicates |

Four findings need to change from the round-2 text:
- **4b:** direction of the field is not useless. Node-local memory dominates, but direction adds a small, consistent amount.
- **8:** the advantage of a society's own order over the shared order is not robust.
- **3:** holds on B but not in W, where the other signs carry a little information (as the pairwise sign-coupling test found earlier).
- **5b:** the exact shared order is panel-dependent at the Shield position.

## Revised conclusions (incorporating the reviewer's edits)

**What was ruled out.** Several attractive but unsupported pictures have been ruled out *for these
datasets and model families*:
- a universal pairwise grammar;
- relations that switch on in particular regimes or during recovery;
- a transportable 8×8 lagged causal graph.

Unrestricted and conditional relational models lose to the factor architecture out of sample. No directed
lagged edge survives false-discovery correction on either panel. This gives little support for treating
CAMS as a fixed causal graph that can be estimated once and transported across societies.

**The central result.** The eight functions behave as a shared response basis, non-exchangeable
coordinates in a common low-dimensional field, not as a fixed relational circuit. The same functional
basis appears in substantially different dynamical morphologies.

**The threshold chain.** Among split configurations (all-positive and all-negative states excluded), about
71% lie within one node flip of the shared threshold chain in both panels. Between 13% (B) and 38% (W) lie
exactly on it. The broad structure is robust: Hands at the hard end, Stewards next to it, Lore and Archive
towards the easy end. Fine positions are not, and Shield in particular moves from easiest (B) to
near-hardest (W).

**Synchronisation.** Synchronisation is a major empirical axis of global morphology *on the blind panel*.
There, societies dominated more strongly by the first common factor have smaller effective repertoires
(ρ ≈ −0.83) and lower switching rates (ρ ≈ −0.74). This is consistent with stronger synchronisation
constraining the accessible state space. The association replicates in the 20 replication societies but is
absent within the 11 canonical ones, so it should be reported as an association, and one that depends on
the panel.

**Morphology.** Repertoires and switching rates differ widely. Split-half tests show that at least part of
the society-specific morphology, namely the threshold order, persists through time rather than behaving
like arbitrary sampling variation (B p ≈ 0.002–0.02; W p = 0.0005).

**Hierarchy from portable structure to contingent morphology.**
1. Functional identity and response roles
2. Relative threshold and loading structure
3. Society-specific threshold order
4. Configuration repertoire and transition geometry
5. Exact coalitions and pairwise edges

- Levels 1–2 show the strongest cross-society portability (loadings transfer out of sample; the broad order recurs).
- Levels 3–4 characterise persistent society-specific morphology.
- Level 5 shows little evidence of universal portability.

**Structured deformation, stated at the strength the evidence allows.** Variation around the shared basis
is not arbitrary: each society's threshold order is stable across the two halves of its history in both
panels. Whether that society-specific order *predicts* later configurations better than the shared order is
unresolved. It does on the blind panel's full sample (27.9% vs 11.8% exact, p = 0.019), but not with a
minimum-data filter (p = 0.099) and not in the corpus (37.9% vs 37.6%). The supported claim is:

> shared broad response hierarchy + persistent society-specific deformation

Proven predictive superiority of the deformation is not yet part of it.

**Common state vs functional composition.** The latent dimensions capture much of the common state. The
signed configuration and node histories retain information about functional composition. Near-identical
aggregate values M correspond to different sign configurations: 58% of temporally separated pairs on B,
34% and 8% on the corpus sets, where carried-forward rows inflate the matches. CAMS may therefore be useful
not because eight variables are needed to forecast the common trajectory, but because they show how a
low-dimensional macrostate is realised functionally.

> Low dimensionality of dynamics does not imply low dimensionality of functional description.

**History.** History is retained mainly through the persistence of individual functional states, which is
the largest effect in every test, and through persistent society-specific morphology. The field's recent
direction adds a small further amount. There is no transportable network of pairwise causal edges.

**Interventions.** Any intervention interpretation that treats an individual node in isolation would omit
two empirically important features of the observed dynamics: the common latent field and node-specific
persistence. These are observational forecasting results and establish nothing about intervention effects.

## Final formulation

> CAMS is not supported as a universal eight-node causal circuit. It is better supported as an
> eight-channel functional coordinate system embedded in a low-dimensional dynamical field. Node identities,
> relative response roles and broad threshold structure recur across societies. Exact threshold ordering,
> repertoire geometry and temporal traversal acquire persistent society-specific form. Aggregate recovery is
> therefore many-to-one with respect to functional configuration. Historical state is retained chiefly
> through node-specific persistence and society-specific morphology, not through a transportable network of
> pairwise causal edges.

**Epistemic boundary.** The present analysis does not identify what generates the common field. Political,
institutional, economic, demographic, ecological and other historical processes may all contribute to it.
The statistical result is only that much of the eight-node covariance can be represented by a small number
of shared dimensions. The tests establish statistical organisation and cross-society portability of parts of
the response basis. They do not establish:
- that the latent dimensions are fundamental social causes;
- that the eight channels are mathematically scale-invariant;
- that intervening on one channel would have the effects implied by the observational dynamics.

Two measurement caveats also remain:
- all scorers are LLMs that may share priors, so part of the common field or of the node order could reflect a scoring prior;
- after two factors, the residual structure is within about 2× the scorer disagreement.
