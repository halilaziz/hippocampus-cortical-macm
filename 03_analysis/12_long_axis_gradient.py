"""Long-axis gradient: bilateral hippocampus split into 6-mm coronal bins; for each bin, specific co-activation
(MKDA chi2 vs foci-matched reference) and Yeo-7 network profile of the positive association z-map."""
import sys, os, json, numpy as np, nibabel as nib, pandas as pd
sys.path.insert(0, "03_analysis")
from lib import *
from importlib import import_module
spec = import_module("06_specificity")
from nimare.meta.cbma.mkda import MKDAChi2
from nimare.meta.kernel import MKDAKernel
from nimare.correct import FDRCorrector
db = sys.argv[1] if len(sys.argv) > 1 else "pooled"
ds = pooled_dataset() if db == "pooled" else (nq_independent() if db == "nq_indep" else load(db))
m = nib.load(f"{SEEDDIR}/BIL_HIPP.nii.gz"); M = np.asarray(m.dataobj) > 0
ijk = np.argwhere(M); Y = nib.affines.apply_affine(m.affine, ijk)[:, 1]
edges = np.arange(np.floor(Y.min()/2)*2, Y.max()+6, 6)
out = f"04_results/gradient/{db}"; os.makedirs(out, exist_ok=True)
rows = []
for lo, hi in zip(edges[:-1], edges[1:]):
    sel = ijk[(Y >= lo) & (Y < hi)]
    if len(sel) < 20: continue
    B = np.zeros(M.shape, np.uint8); B[tuple(sel.T)] = 1
    name = f"y{int(lo)}_{int(hi)}"; f = f"{out}/{name}_mask.nii.gz"; nib.save(nib.Nifti1Image(B, m.affine), f)
    d = f"{out}/{name}"
    if not os.path.exists(f"{d}/info.json"):
        os.makedirs(d, exist_ok=True)
        ids = list(ds.get_studies_by_mask(f)); ref, m1, m2 = spec.matched_reference(ds, ids, 42)
        res = MKDAChi2(kernel_transformer=MKDAKernel(r=10)).fit(ds.slice(ids), ds.slice(ref))
        FDRCorrector(method="indep", alpha=0.01).transform(res).save_maps(output_dir=d)
        json.dump(dict(bin=name, y_lo=float(lo), y_hi=float(hi), n_vox=int(len(sel)), n_seed=len(ids), mean_foci=m1), open(f"{d}/info.json", "w"))
    print(open(f"{d}/info.json").read())
