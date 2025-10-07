import demes
import numpy as np
import json
import sys
from pathlib import Path

demes_yaml = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("src/demes.yaml")
out_json   = Path(sys.argv[2]) if len(sys.argv) > 2 else Path("src/migration_matrices.json")

model = demes.load(str(demes_yaml))

matrices, times = model.migration_matrices()

def adjust_diagonal(m):
    m = np.array(m, dtype=float)
    row_sums_off = np.sum(m, axis=1) - np.diag(m)
    np.fill_diagonal(m, 1.0 - row_sums_off)
    return np.round(m, 5).tolist()


pairs = sorted(zip(times, matrices), key=lambda x: x[0])
payload = { str(int(t)): adjust_diagonal(m) for t, m in pairs }

with open(out_json, "w") as f:
    json.dump(payload, f, indent=4)

print(f"Wrote {out_json} with {len(payload)} matrices:", list(payload.keys()))
