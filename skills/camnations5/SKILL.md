---
name: camnations5
description: Five-scorer ensemble CAMS assessment pipeline for a single nation across multiple years. Runs cams-scorer five times independently under an identical prompt, averages raw Coherence/Capacity/Stress/Abstraction at the node-year level, then runs cams-calc once on the mean scores. Outputs a long-format ensemble-mean CSV plus a companion envelope CSV capturing inter-scorer variance. Use this skill — even proactively — whenever the user says "CAMNATIONS5", "camnations5", "five-scorer CAMS", "ensemble CAMS", or asks for a scorer-robust, uncertainty-aware, or publication-ready version of a CAMNATIONS run. Does not replace camnations (v1); coexists with it.
---

# CAMNATIONS5 — Five-Scorer Ensemble CAMS Pipeline

Ensemble variant of the CAMNATIONS pipeline. Run five independent scoring passes under an identical prompt, average the raw scores at each node-year cell, then compute Node Value and Bond Strength once from the means. Emit a consensus dataset plus the scorer-disagreement envelope.

## Input Parameters

| Parameter | Type | Required | Default | Notes |
|-----------|------|----------|---------|-------|
| nation | string | yes | — | Society/nation to assess |
| start_year | integer | yes | — | First year to assess |
| end_year | integer | no | 2026 | Last year. Cap to 2026 if > 2026 |
| interval | string | no | yearly | Use `yearly` unless the user explicitly asks for `6_month` intervals |
| n_scorers | integer | no | 5 | Independent passes. Range: 3–7; default is always exactly 5 |

### Validation rules
- `n_scorers < 3` → fail with explicit error
- `n_scorers > 7` → cap to 7
- `start_year > end_year` → fail with explicit error
- `end_year > 2026` → cap to 2026
- If `end_year` is omitted, set it to 2026
- If `interval` is omitted, use `yearly`
- Only use `6_month` intervals when explicitly requested by the user

## Time coverage rules
- Default coverage is **every calendar year** from the mentioned `start_year` through 2026, inclusive.
- If the user gives a start date/year but no end date/year, score `start_year, start_year + 1, ... , 2026`; do not stop early.
- Do **not** hop across milestone years and later fill gaps. No decade sampling, representative-year sampling, or interpolation is allowed.
- Every subsequent year in the range must be scored explicitly with all 8 nodes.
- Use `6_month` intervals only when the user explicitly asks for six-month intervals; otherwise remain yearly.
- For `6_month` intervals, score every half-year point from the mentioned start date through the specified end date, or through 2026 if no end date is specified. Do not skip half-year points and do not backfill inferred values.
- The `n_scorers` default is always **5** unless the user explicitly supplies another valid value.

## Output contract
Emit exactly two CSV code blocks in this order:

1) **Block 1 — Ensemble Mean**
```csv
Society,Year,Node,Coherence,Capacity,Stress,Abstraction,Node Value,Bond Strength
```
- Mean C/K/S/A across scorers, rounded to **1 decimal place**
- Run `cams-calc` (or `cams-rater2`) exactly once on the mean scores to compute Node Value and Bond Strength
- Sort rows by `Society`, then `Year`, then canonical `Node` order: `Helm`, `Shield`, `Lore`, `Stewards`, `Craft`, `Hands`, `Archive`, `Flow`

2) **Block 2 — Envelope**
```csv
Society,Year,Node,C_sd,K_sd,S_sd,A_sd,V_range,V_min,V_max
```
- Sample SD per dimension across passes (`ddof=1`)
- Per-scorer Node Value: `V_i = C_i + K_i - S_i + 0.5*A_i`
- Envelope: `V_min`, `V_max`, `V_range = V_max - V_min`
- Sort rows by `Society`, then `Year`, then canonical `Node` order: `Helm`, `Shield`, `Lore`, `Stewards`, `Craft`, `Hands`, `Archive`, `Flow`

## Core workflow

1. Run `cams-scorer` (CAMS Scoring Protocol v2.1) **sequentially** for each pass (`pass_1 ... pass_n`) using identical prompt settings.
2. For each `(Society, Year, Node)` cell, compute arithmetic means for C/K/S/A and round to 1 decimal.
3. Compute SD envelope stats and V-range diagnostics from pass-level scores.
4. If any pass has missing rows, use available rows cell-by-cell.
5. If any cell has fewer than 3 contributing passes, fail the run and report failing cells.
6. Run `cams-calc` once over ensemble means.
7. Sort both final output blocks by `Society`, `Year`, and canonical `Node` order (`Helm`, `Shield`, `Lore`, `Stewards`, `Craft`, `Hands`, `Archive`, `Flow`).
8. Output only the two CSV blocks; no interleaved narration unless requested.

## Scorer isolation rules
- No cross-pass peeking
- No pre-scoring historical narrative
- No pass-to-pass correction/smoothing
- Same scorer prompt for all passes
- Sequential execution only (not parallel)
- No pass-level commentary between scoring runs

## Statistical conventions
- Integer scores inside each pass (1–10)
- Fractional values allowed only in post-ensemble means
- Node Value linearity implies mean(V) == V(mean scores)
- Bond Strength is non-linear; canonical output computes it from mean scores

## Recommended use
- Use CAMNATIONS5 for publication-facing, comparative, or uncertainty-aware work
- Keep CAMNATIONS v1 for rapid exploratory workflows
