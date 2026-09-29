"""NiMARE ALE of rsFC coordinate set (true N per experiment), cluster-FWE, plus Jackknife contribution analysis and
leave-one-experiment-out robustness (Eickhoff et al., 2016)."""
import sys, os, json, numpy as np, pandas as pd, nibabel as nib
from nimare.io import convert_sleuth_to_dataset
from nimare.meta.cbma.ale import ALE
from nimare.correct import FWECorrector
from nimare.diagnostics import Jackknife
f = sys.argv[1]; name = os.path.basename(f).replace(".txt", ""); out = f"04_results/rsfc_ale/nimare/{name}"; os.makedirs(out, exist_ok=True)
ds = convert_sleuth_to_dataset(f, target="mni152_2mm")
res = ALE(null_method="approximate").fit(ds)
cres = FWECorrector(method="montecarlo", voxel_thresh=0.001, n_iters=1000, n_cores=12, vfwe_only=False).transform(res)
cres.save_maps(output_dir=out)
jk = Jackknife(target_image="z_desc-size_level-cluster_corr-FWE_method-montecarlo", voxel_thresh=None, n_cores=12)
cr = jk.transform(cres)
tabs = cr.tables if hasattr(cr, "tables") else {}
for k, t in tabs.items():
    if t is not None and len(t): t.to_csv(f"{out}/jackknife_{k}.csv", index=True); print(k); print(t.round(3).to_string()[:3000])
json.dump(dict(name=name, n_experiments=len(ds.ids), n_foci=len(ds.coordinates)), open(f"{out}/info.json", "w"))
