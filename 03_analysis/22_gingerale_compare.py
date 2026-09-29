"""Voxelwise equivalence of native GingerALE 3.0.2 ALE maps vs NiMARE ALE maps (same foci, same kernel N=20, non-additive)."""
import glob, os, numpy as np, nibabel as nib, pandas as pd
from nilearn import image
rows = []
for g in sorted(glob.glob("04_results/gingerale/out/*_ALE.nii")):
    b = os.path.basename(g).replace("_ALE.nii", "")
    if b.startswith("test_") or "_p001_" in b: continue
    db, seed = b.split("_", 1)
    n = f"04_results/macm/{db}/{seed}/stat.nii.gz"
    if not os.path.exists(n): continue
    N = nib.load(n); G = image.resample_to_img(g, N, interpolation="nearest", force_resample=True, copy_header=True)
    a, c = np.asarray(N.dataobj), np.squeeze(np.asarray(G.dataobj))
    m = (a > 0) & (c > 0)
    r = np.corrcoef(a[m], c[m])[0, 1]
    ia, ic = np.unravel_index(np.argmax(a), a.shape), np.unravel_index(np.argmax(c), c.shape)
    rows.append(dict(db=db, seed=seed, n_vox=int(m.sum()), pearson_r=round(r, 5), max_ALE_nimare=round(float(a.max()), 4), max_ALE_gingerale=round(float(c.max()), 4),
                     median_abs_rel_diff_pct=round(100 * float(np.median(np.abs(a[m] - c[m]) / np.maximum(c[m], 1e-9))), 3),
                     peak_distance_mm=round(float(np.linalg.norm((np.array(ia) - np.array(ic)) * 2)), 1)))
df = pd.DataFrame(rows); df.to_csv("04_results/tables/gingerale_validation.csv", index=False); print(df.to_string(index=False))
