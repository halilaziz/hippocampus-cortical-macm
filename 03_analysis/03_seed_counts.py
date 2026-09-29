import pandas as pd, json
from nimare.dataset import Dataset
res = {}
dss = {n: Dataset.load(f"02_data/{n}_dataset.pkl.gz") for n in ["neurosynth", "neuroquery"]}
for n, ds in dss.items():
    res[n] = {"n_studies": len(ds.ids), "n_foci": len(ds.coordinates)}
    for s in ["L_HIPP","R_HIPP","BIL_HIPP","L_aHIPP","L_pHIPP","R_aHIPP","R_pHIPP","BIL_aHIPP","BIL_pHIPP"]:
        ids = ds.get_studies_by_mask(f"04_results/seeds/{s}.nii.gz")
        res[n][s] = len(ids)
    print(n, res[n])
# overlap by PMID
p = {n: set(i.split("-")[0] for i in ds.ids) for n, ds in dss.items()}
print("PMID overlap:", len(p["neurosynth"] & p["neuroquery"]), "NS only", len(p["neurosynth"]-p["neuroquery"]), "NQ only", len(p["neuroquery"]-p["neurosynth"]))
print(list(dss["neurosynth"].ids[:3]), list(dss["neuroquery"].ids[:3]))
print(dss["neurosynth"].coordinates.head()); print(dss["neurosynth"].metadata.columns.tolist())
