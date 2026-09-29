"""Download Neurosynth (v7) and NeuroQuery databases as NiMARE Datasets (MNI152 2mm)."""
import os
from nimare.extract import fetch_neurosynth, fetch_neuroquery

out = "02_data"
for name, fetch in [("neurosynth", fetch_neurosynth), ("neuroquery", fetch_neuroquery)]:
    kw = dict(data_dir=out, overwrite=False, return_type="dataset")
    kw.update(dict(source="abstract", vocab="terms") if name == "neurosynth" else dict(source="combined", vocab="neuroquery6308", type="tfidf"))
    if name == "neurosynth":
        kw["version"] = "7"
    ds = fetch(**kw)[0]
    ds.save(os.path.join(out, f"{name}_dataset.pkl.gz"))
    print(name, "studies:", len(ds.ids), "foci:", len(ds.coordinates),
          "annotations:", ds.annotations.shape, "space:", ds.space)
