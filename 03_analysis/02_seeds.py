"""Build hippocampal seed masks (MNI152 2mm) from Harvard-Oxford (maxprob thr50) and split long axis at y=-21 (uncal apex; Poppenk et al., 2013)."""
import numpy as np, nibabel as nib, os
from nilearn import datasets, image
out = "04_results/seeds"; os.makedirs(out, exist_ok=True)
ho = datasets.fetch_atlas_harvard_oxford("sub-maxprob-thr50-2mm")
img = image.load_img(ho.maps); lab = list(ho.labels)
data = np.asarray(img.dataobj)
print(lab)
aff = img.affine
ijk = np.indices(data.shape).reshape(3, -1).T
xyz = nib.affines.apply_affine(aff, ijk).reshape(data.shape + (3,))
Y = xyz[..., 1]
seeds = {}
for hemi in ["Left", "Right"]:
    idx = lab.index(f"{hemi} Hippocampus")
    m = data == idx
    h = hemi[0]
    seeds[f"{h}_HIPP"] = m
    seeds[f"{h}_aHIPP"] = m & (Y >= -21)
    seeds[f"{h}_pHIPP"] = m & (Y < -21)
seeds["BIL_HIPP"] = seeds["L_HIPP"] | seeds["R_HIPP"]
seeds["BIL_aHIPP"] = seeds["L_aHIPP"] | seeds["R_aHIPP"]
seeds["BIL_pHIPP"] = seeds["L_pHIPP"] | seeds["R_pHIPP"]
for k, m in seeds.items():
    nib.save(nib.Nifti1Image(m.astype(np.uint8), aff), f"{out}/{k}.nii.gz")
    print(k, int(m.sum()), "voxels", int(m.sum()) * 8, "mm3")
