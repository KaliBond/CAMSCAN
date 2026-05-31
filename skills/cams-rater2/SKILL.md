---
name: cams-rater2
version: 2.0.0
spec: CAMS-CAN v1.0-FINAL
description: >
  Compute Node Value and Bond Strength from any raw CAMS score dataset and
  deliver the result as a single CSV code block. Trigger on: "rate", "rater 2",
  "cams-rater2", "node value", "bond strength", "run rater".
---

# CAMS RATER 2

Takes raw CAMS scores. Calculates Node Value and Bond Strength. Delivers one CSV code block. Nothing else.

---

## Input

Required columns: `Society, Year, Node, Coherence, Capacity, Stress, Abstraction`

Nodes per year (8): `Helm, Shield, Lore, Stewards, Craft, Hands, Archive, Flow`

Source: uploaded CSV, inline CSV block, or named MARKER file.

---

## Equations

### Node Value — Equation 1

```text
Node Value = Coherence + Capacity − Stress + 0.5 × Abstraction
```

### Bond Strength — Equations 2 & 4

Build the 8×8 pairwise matrix for each year:

```text
B_ij = [0.6 × C_i × C_j + 0.4 × A_i × A_j] × exp(−1.0 × (S_i + S_j) / 20)
```

Bond Strength per node = mean of its 7 pairwise B_ij values (j ≠ i).

λ = 1.0. Fixed. Do not change.

---

## Implementation

```python
import pandas as pd
import numpy as np

NODE_ORDER = ['Helm','Shield','Lore','Stewards','Craft','Hands','Archive','Flow']

df = pd.read_csv(input_path)
for col in ['Coherence','Capacity','Stress','Abstraction']:
    df[col] = df[col].astype(float)

df['Node Value'] = (
    df['Coherence'] + df['Capacity'] - df['Stress'] + 0.5 * df['Abstraction']
).round(1)

rows = []
for (society, year), group in df.groupby(['Society', 'Year'], sort=True):
    group = group.set_index('Node').reindex(NODE_ORDER).reset_index()
    if group[['Coherence', 'Capacity', 'Stress', 'Abstraction']].isna().any().any():
        continue

    C = group['Coherence'].values.astype(float)
    A = group['Abstraction'].values.astype(float)
    S = group['Stress'].values.astype(float)
    n = len(group)
    B = np.zeros((n, n))

    for i in range(n):
        for j in range(i + 1, n):
            linear = 0.6 * C[i] * C[j] + 0.4 * A[i] * A[j]
            decay = np.exp(-1.0 * (S[i] + S[j]) / 20.0)
            B[i, j] = linear * decay
            B[j, i] = B[i, j]

    group['Bond Strength'] = np.round(np.sum(B, axis=1) / (n - 1), 3)
    rows.append(group)

out = pd.concat(rows)[
    ['Society','Year','Node','Coherence','Capacity','Stress','Abstraction','Node Value','Bond Strength']
]
out['Node'] = pd.Categorical(out['Node'], categories=NODE_ORDER, ordered=True)
out = out.sort_values(['Society', 'Year', 'Node']).reset_index(drop=True)
```

---

## Output

**One CSV code block. Nothing else.**

```csv
Society,Year,Node,Coherence,Capacity,Stress,Abstraction,Node Value,Bond Strength
```

- All years. All 8 nodes per year. No truncation.
- Sort output rows by `Society`, then `Year`, then canonical `Node` order: `Helm`, `Shield`, `Lore`, `Stewards`, `Craft`, `Hands`, `Archive`, `Flow`.
- Node Value: 1 decimal place.
- Bond Strength: 3 decimal places.
- No files. No commentary. No `present_files`.

---

## Constraints

- Never alter raw scores
- Never skip years in the source data
- If a year has fewer than 8 nodes, skip its bond calculation
