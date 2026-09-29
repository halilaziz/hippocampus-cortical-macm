#!/bin/bash
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
# Full ALE-MACM batch, 1000 Monte Carlo iterations each
cd "$ROOT"
export NCORES=12
for db in pooled neurosynth nq_indep; do
  for seed in BIL_HIPP L_HIPP R_HIPP BIL_aHIPP BIL_pHIPP L_aHIPP L_pHIPP R_aHIPP R_pHIPP; do
    if [ -f 04_results/macm/$db/$seed/info.json ]; then continue; fi
    echo "=== $db $seed $(date)"
    .venv/bin/python 03_analysis/04_macm_ale.py $db $seed 1000 2>&1 | grep -E "^\{|Error|Traceback"
  done
done
echo "BATCH DONE $(date)"
