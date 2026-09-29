#!/bin/bash
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT"
03_analysis/16_rsfc_ale.sh 04_results/rsfc_ale/rsfc_main_v4.txt 04_results/rsfc_ale/rsfc_plus_sensitivity_v4.txt
.venv/bin/python 03_analysis/33_rsfc_nimare_jackknife.py 04_results/rsfc_ale/rsfc_main_v4.txt > 04_results/rsfc_ale/nimare_v4.log 2>&1
.venv/bin/python - <<'PY'
import pandas as pd, glob
d = "04_results/rsfc_ale/nimare/rsfc_main_v4/"
t = pd.read_csv(glob.glob(d + "jackknife_*counts_tail-positive*.csv")[0], index_col=0).set_index("id")
rows = []
for c in t.columns:
    v = t[c].sort_values(ascending=False); tot = v.sum()
    rows.append(dict(cluster=c, n_contrib_gt5pct=int((v / tot > 0.05).sum()) if tot > 0 else 0, max_share_pct=round(100 * v.iloc[0] / tot, 1) if tot > 0 else 0, top_experiment=v.index[0].split("-")[0]))
pd.DataFrame(rows).to_csv(d + "jackknife_summary.csv", index=False); print(pd.DataFrame(rows).to_string(index=False))
PY
echo V4_DONE
