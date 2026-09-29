"""Collect key numbers from all result files into one JSON used by the report/manuscript builders."""
import json, glob, os, pandas as pd, numpy as np
R = {}
R["prisma"] = json.load(open("01_literature/prisma_counts.json"))
m = []
for f in glob.glob("04_results/macm/*/*/info.json"):
    i = json.load(open(f)); m.append(i)
R["macm_info"] = sorted(m, key=lambda x: (x["db"], x["seed"]))
s = []
for f in glob.glob("04_results/specificity/*/*/info.json"):
    s.append(json.load(open(f)))
R["spec_info"] = sorted(s, key=lambda x: (x["db"], x["seed"]))
for k, f in [("macm_summary", "04_results/tables/macm_summary.csv"), ("spec_summary", "04_results/tables/specificity_summary.csv"),
             ("replication", "04_results/tables/replication.csv"), ("distal", "04_results/tables/network_profiles_distal12mm.csv")]:
    if os.path.exists(f): R[k] = pd.read_csv(f).to_dict(orient="records")
for db in ["pooled", "neurosynth", "nq_indep"]:
    f = f"04_results/tables/gradient_{db}_trend.csv"
    if os.path.exists(f): R[f"gradient_trend_{db}"] = pd.read_csv(f).to_dict(orient="records")
R["contrasts"] = [json.load(open(f)) for f in glob.glob("04_results/contrasts/*/*/info.json")]
for f in glob.glob("04_results/tables/contrasts/*/*_summary.json"):
    R.setdefault("contrast_summaries", []).append(json.load(open(f)))
if os.path.exists("04_results/tables/gingerale_validation.csv"):
    R["gingerale_validation"] = pd.read_csv("04_results/tables/gingerale_validation.csv").to_dict(orient="records")
for f in glob.glob("04_results/rsfc_ale/nimare/*/info.json"):
    R.setdefault("rsfc", []).append(json.load(open(f)))
json.dump(R, open("04_results/results_summary.json", "w"), indent=1, default=float)
print({k: (len(v) if isinstance(v, list) else "dict") for k, v in R.items()})
