---
name: cams-calc
description: "CALCULATES NODE VALUE AND BOND STRENGTH AND"
---

# CAMS-CALC

Alias skill for CAMS node-level derived metrics.

Use this skill when a workflow needs downstream calculation from raw CAMS scores.

- Input columns: `Society,Year,Node,Coherence,Capacity,Stress,Abstraction`
- Output columns: `Society,Year,Node,Coherence,Capacity,Stress,Abstraction,Node Value,Bond Strength`
- Equations and strict output behavior are defined in `skills/cams-rater2/SKILL.md`.
