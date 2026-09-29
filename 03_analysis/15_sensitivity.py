"""Sensitivity analyses:
 (1) amygdala-proximity: anterior hippocampus specificity after excluding every study with a focus in the amygdala (HO thr50);
 (2) seed-definition: 'core' hippocampus (Harvard-Oxford probability >= 75%) instead of maxprob-thr50;
 both via MKDA chi2 vs foci-matched reference, FDR q<.01."""
import sys, os, json, numpy as np, nibabel as nib
sys.path.insert(0, "03_analysis")
from lib import *
from importlib import import_module
spec = import_module("06_specificity")
from nimare.meta.cbma.mkda import MKDAChi2
from nimare.meta.kernel import MKDAKernel
from nimare.correct import FDRCorrector
from nilearn import datasets, image
db = sys.argv[1]
ds = pooled_dataset() if db == "pooled" else load(db)
ho = datasets.fetch_atlas_harvard_oxford("sub-maxprob-thr50-2mm"); lab = list(ho.labels)
ref = nib.load(f"{SEEDDIR}/BIL_HIPP.nii.gz")
H = np.asarray(image.resample_to_img(ho.maps, ref, interpolation="nearest", force_resample=True, copy_header=True).dataobj)
amy = np.isin(H, [lab.index("Left Amygdala"), lab.index("Right Amygdala")]).astype(np.uint8)
nib.save(nib.Nifti1Image(amy, ref.affine), f"{SEEDDIR}/BIL_AMY.nii.gz")
prob = datasets.fetch_atlas_harvard_oxford("sub-prob-2mm"); pl = list(prob.labels)
P = image.resample_to_img(prob.maps, ref, interpolation="nearest", force_resample=True, copy_header=True).get_fdata()
def vol(name): return P[..., pl.index(name) - (1 if pl[0].lower().startswith("background") else 0)]
core = ((vol("Left Hippocampus") >= 75) | (vol("Right Hippocampus") >= 75)).astype(np.uint8)
nib.save(nib.Nifti1Image(core, ref.affine), f"{SEEDDIR}/BIL_HIPP_core75.nii.gz")
print("core voxels", int(core.sum()), "amygdala voxels", int(amy.sum()))
amy_ids = set(ds.get_studies_by_mask(f"{SEEDDIR}/BIL_AMY.nii.gz"))
jobs = {"BIL_aHIPP_noAmyg": sorted(set(seed_ids(ds, "BIL_aHIPP")) - amy_ids),
        "BIL_pHIPP_noAmyg": sorted(set(seed_ids(ds, "BIL_pHIPP")) - amy_ids),
        "BIL_HIPP_core75": list(ds.get_studies_by_mask(f"{SEEDDIR}/BIL_HIPP_core75.nii.gz"))}
for name, ids in jobs.items():
    out = f"04_results/specificity/{db}/{name}"
    if os.path.exists(f"{out}/info.json"): continue
    os.makedirs(out, exist_ok=True)
    refids, m1, m2 = spec.matched_reference(ds, ids, 42)
    res = MKDAChi2(kernel_transformer=MKDAKernel(r=10)).fit(ds.slice(ids), ds.slice(refids))
    FDRCorrector(method="indep", alpha=0.01).transform(res).save_maps(output_dir=out)
    json.dump(dict(db=db, seed=name, n_seed=len(ids), n_ref=len(refids), mean_foci_seed=round(m1, 1), mean_foci_ref=round(m2, 1), seconds=0), open(f"{out}/info.json", "w"))
    print(name, len(ids))
