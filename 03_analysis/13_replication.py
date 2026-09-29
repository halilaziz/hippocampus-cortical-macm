"""Replication between Neurosynth (discovery) and the non-overlapping NeuroQuery subset (independent replication):
Dice overlap of thresholded maps, Pearson r of unthresholded z-maps (in-brain), and conjunction (minimum statistic) maps."""
import glob, os, json, numpy as np, nibabel as nib, pandas as pd
from nilearn import datasets, image
ref = nib.load(glob.glob("04_results/macm/*/*/z.nii.gz")[0])
brain = np.asarray(image.resample_to_img(datasets.load_mni152_brain_mask(), ref, interpolation="nearest", force_resample=True, copy_header=True).dataobj) > 0
rows = []; os.makedirs("04_results/tables/conjunction", exist_ok=True)
for kind, unthr, tdir in [("macm", "04_results/macm/{db}/{s}/z.nii.gz", "04_results/tables/macm/{db}/{s}_thr.nii.gz"),
                          ("specificity", "04_results/specificity/{db}/{s}/z_desc-association.nii.gz", "04_results/tables/specificity/{db}/{s}_thr.nii.gz")]:
    for s in ["BIL_HIPP", "L_HIPP", "R_HIPP", "BIL_aHIPP", "BIL_pHIPP", "L_aHIPP", "L_pHIPP", "R_aHIPP", "R_pHIPP"]:
        f = {db: tdir.format(db=db, s=s) for db in ["neurosynth", "nq_indep"]}
        u = {db: unthr.format(db=db, s=s) for db in ["neurosynth", "nq_indep"]}
        if not all(os.path.exists(x) for x in list(f.values()) + list(u.values())): continue
        A, B = [np.asarray(nib.load(f[db]).dataobj) for db in ["neurosynth", "nq_indep"]]
        ua, ub = [np.asarray(nib.load(u[db]).dataobj)[brain] for db in ["neurosynth", "nq_indep"]]
        a, b = A > 0, B > 0
        dice = 2 * (a & b).sum() / (a.sum() + b.sum())
        conj = np.where(a & b, np.minimum(A, B), 0).astype(np.float32)
        nib.save(nib.Nifti1Image(conj, ref.affine), f"04_results/tables/conjunction/{kind}_{s}_conj.nii.gz")
        rows.append(dict(analysis=kind, seed=s, vox_neurosynth=int(a.sum()), vox_nq_indep=int(b.sum()), vox_conjunction=int((a & b).sum()),
                         dice=round(dice, 3), pct_ns_replicated=round(100 * (a & b).sum() / a.sum(), 1), pct_nq_in_ns=round(100 * (a & b).sum() / b.sum(), 1), r_unthresholded=round(np.corrcoef(ua, ub)[0, 1], 3)))
df = pd.DataFrame(rows); df.to_csv("04_results/tables/replication.csv", index=False); print(df.to_string(index=False))
