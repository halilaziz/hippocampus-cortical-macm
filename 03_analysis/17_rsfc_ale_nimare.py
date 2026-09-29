"""NiMARE ALE of the published rsFC foci (same settings as GingerALE, true per-experiment N), plus
leave-one-experiment-out jackknife and per-cluster contribution analysis (Eickhoff et al., 2016)."""
import sys, os, json, numpy as np, pandas as pd, nibabel as nib
from nimare.io import convert_sleuth_to_dataset
from nimare.meta.cbma.ale import ALE
from nimare.correct import FWECorrector
from nimare.diagnostics import Jackknife
sys.path.insert(0, "03_analysis")
from importlib import import_module
nt = import_module("08_networks_tables")
f = sys.argv[1]; name = os.path.basename(f).replace(".txt", "")
out = f"04_results/rsfc_ale/nimare/{name}"; os.makedirs(out, exist_ok=True)
ds = convert_sleuth_to_dataset(f, target="mni152_2mm")
res = ALE(null_method="approximate").fit(ds)
cres = FWECorrector(method="montecarlo", voxel_thresh=0.001, n_iters=1000, n_cores=12, vfwe_only=False).transform(res)
cres.save_maps(output_dir=out)
n_sub = ds.metadata["sample_sizes"].apply(lambda x: x[0] if isinstance(x, (list, tuple)) else x).sum() if "sample_sizes" in ds.metadata else np.nan
info = dict(name=name, n_experiments=len(ds.ids), n_foci=len(ds.coordinates), n_subjects=float(n_sub))
json.dump(info, open(f"{out}/info.json", "w")); print(info)
cl = f"{out}/logp_desc-size_level-cluster_corr-FWE_method-montecarlo.nii.gz"
L = nib.load(cl); lp = np.asarray(L.dataobj)
sig = np.where(lp > 1.3, np.asarray(nib.load(f"{out}/z.nii.gz").dataobj), 0).astype(np.float32)
nib.save(nib.Nifti1Image(sig, L.affine), f"{out}/sig_clusterFWE.nii.gz")
s, prof = nt.process(f"{out}/sig_clusterFWE.nii.gz", "04_results/tables/rsfc_ale", name, zcrit=1e-6, min_vox=1, statfile=f"{out}/z.nii.gz", alefile=f"{out}/stat.nii.gz")
print(s); print(prof.round(1).to_string(index=False))
try:
    jk = Jackknife(target_image="z_desc-size_level-cluster_corr-FWE_method-montecarlo", voxel_thresh=None)
    cr, _ = jk.transform(cres)
    t = cr.tables[list(cr.tables)[0]] if hasattr(cr, "tables") else cr[0]
    t.to_csv(f"04_results/tables/rsfc_ale/{name}_jackknife_contributions.csv", index=False); print(t.head())
except Exception as e:
    print("jackknife failed:", e)
