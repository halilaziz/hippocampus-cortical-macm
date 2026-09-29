"""Network profiles restricted to cortex >12 mm from the hippocampal seed (removes kernel spill-over into adjacent
parahippocampal/ventral temporal cortex), for MACM and specificity maps."""
import sys, glob, os, numpy as np, nibabel as nib, pandas as pd
from scipy import ndimage
sys.path.insert(0, "03_analysis")
from importlib import import_module
nt = import_module("08_networks_tables")
H = np.asarray(nib.load("04_results/seeds/BIL_HIPP.nii.gz").dataobj) > 0
dist = ndimage.distance_transform_edt(~H, sampling=2.0)
far = dist > 12
rows = []
for kind in ["macm", "specificity"]:
    for f in sorted(glob.glob(f"04_results/tables/{kind}/*/*_thr.nii.gz")):
        db = f.split("/")[-2]; seed = os.path.basename(f).replace("_thr.nii.gz", "")
        m = (np.asarray(nib.load(f).dataobj) > 0) & far
        prof, ncort = nt.network_profile(m)
        rows.append(dict(analysis=kind, db=db, seed=seed, n_cortical_vox_distal=ncort, **{r.network: round(r.pct_of_cortical_map, 1) for r in prof.itertuples()}))
df = pd.DataFrame(rows); df.to_csv("04_results/tables/network_profiles_distal12mm.csv", index=False)
print(df[df.db.isin(["pooled", "nq_indep", "neurosynth"])].to_string(index=False))
