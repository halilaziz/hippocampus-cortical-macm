"""Label native-GingerALE rsFC ALE clusters (cluster-FWE p<.05) and test overlap with MACM specificity map."""
import sys, os, json, numpy as np, nibabel as nib, pandas as pd
import os
RSV = "v4" if os.path.exists("04_results/rsfc_ale/gingerale/rsfc_main_v4_p001_C05_1k_Z.nii") else "v3"
from nilearn import image
sys.path.insert(0, "03_analysis")
from importlib import import_module
nt = import_module("08_networks_tables")
ref = nib.load("04_results/macm/neurosynth/BIL_HIPP/z.nii.gz")
spec = np.asarray(nib.load("04_results/tables/specificity/pooled/BIL_HIPP_thr.nii.gz").dataobj) > 0
macm = np.asarray(nib.load("04_results/tables/specificity/pooled/BIL_HIPP_thr.nii.gz").dataobj) > 0
rows = []
NAMES = [f"rsfc_main_{RSV}", f"rsfc_plus_sensitivity_{RSV}"]
for name in NAMES:
    gdir = "04_results/rsfc_ale/gingerale"
    z = image.resample_to_img(f"{gdir}/{name}_p001_C05_1k_Z.nii", ref, interpolation="nearest", force_resample=True, copy_header=True)
    zz = np.squeeze(np.asarray(z.dataobj)).astype(np.float32)
    f = f"{gdir}/{name}_sig_resampled.nii.gz"; nib.save(nib.Nifti1Image(zz, ref.affine), f)
    ale = image.resample_to_img(f"{gdir}/{name}_ALE.nii", ref, interpolation="nearest", force_resample=True, copy_header=True)
    af = f"{gdir}/{name}_ALE_resampled.nii.gz"; nib.save(nib.Nifti1Image(np.squeeze(np.asarray(ale.dataobj)).astype(np.float32), ref.affine), af)
    s, prof = nt.process(f, "04_results/tables/rsfc_ale", name, zcrit=1e-6, min_vox=1, alefile=af)
    m = zz > 0
    n_exp = open(f"04_results/rsfc_ale/{name}.txt").read().count("Subjects=")
    n_sub = sum(int(l.split("=")[1]) for l in open(f"04_results/rsfc_ale/{name}.txt") if "Subjects=" in l)
    n_foci = sum(1 for l in open(f"04_results/rsfc_ale/{name}.txt") if l.strip() and not l.startswith("//"))
    s.pop("name", None); rows.append(dict(name=name, n_experiments=n_exp, n_subjects=n_sub, n_foci=n_foci, **s, pct_in_macm_specific=round(100 * (m & spec).sum() / max(m.sum(), 1), 1)))
    print(rows[-1]); print(pd.read_csv(f"04_results/tables/rsfc_ale/{name}_clusters.csv").to_string(index=False))
pd.DataFrame(rows).to_csv("04_results/tables/rsfc_ale_summary.csv", index=False)
os.makedirs("04_results/rsfc_ale/nimare/main", exist_ok=True)
json.dump(rows[0], open("04_results/rsfc_ale/nimare/main/info.json", "w"))
