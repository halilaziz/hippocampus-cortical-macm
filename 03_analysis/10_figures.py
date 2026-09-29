"""Publication figures: surface renderings of thresholded maps, network profiles."""
import sys, os, glob, json, numpy as np, pandas as pd, nibabel as nib, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from nilearn import plotting, datasets, surface
FIG = "04_results/figures"; os.makedirs(FIG, exist_ok=True)
fs = datasets.fetch_surf_fsaverage("fsaverage5")
def surf4(img, title, fname, cmap="hot", vmax=None, symmetric=False, threshold=1e-3):
    fig, axes = plt.subplots(1, 4, figsize=(14, 3.4), subplot_kw={"projection": "3d"})
    views = [("left", "lateral"), ("left", "medial"), ("right", "medial"), ("right", "lateral")]
    data = {h: surface.vol_to_surf(img, fs[f"pial_{h}"], radius=3, interpolation="linear", kind="line") for h in ["left", "right"]}
    if vmax is None: vmax = max(np.nanmax(np.abs(d)) for d in data.values())
    for ax, (h, v) in zip(axes, views):
        plotting.plot_surf_stat_map(fs[f"infl_{h}"], data[h], hemi=h, view=v, bg_map=fs[f"sulc_{h}"], axes=ax,
                                    cmap=cmap, vmax=vmax, threshold=threshold, colorbar=(ax is axes[-1]), symmetric_cbar=symmetric)
    fig.suptitle(title, fontsize=13); plt.savefig(f"{FIG}/{fname}", dpi=200, bbox_inches="tight"); plt.close(fig)
def slices(img, title, fname, cmap="hot", cut=(-54, -40, -26, -12, 2, 16, 30, 44)):
    d = plotting.plot_stat_map(img, display_mode="z", cut_coords=cut, cmap=cmap, title=title, threshold=1e-3, colorbar=True, black_bg=False, draw_cross=False)
    d.savefig(f"{FIG}/{fname}", dpi=200); d.close()
if __name__ == "__main__":
    which = sys.argv[1:] or ["macm"]
    if "macm" in which:
        for f in sorted(glob.glob("04_results/tables/macm/*/*_thr.nii.gz")):
            db = f.split("/")[-2]; seed = os.path.basename(f).replace("_thr.nii.gz", "")
            surf4(f, f"ALE-MACM {seed} ({db}; voxel-FWE p<.05)", f"macm_{db}_{seed}_surf.png")
            slices(f, f"{seed} {db}", f"macm_{db}_{seed}_axial.png")
    if "spec" in which:
        for f in sorted(glob.glob("04_results/tables/specificity/*/*_thr.nii.gz")):
            db = f.split("/")[-2]; seed = os.path.basename(f).replace("_thr.nii.gz", "")
            surf4(f, f"Specific co-activation {seed} ({db}; FDR q<.01)", f"spec_{db}_{seed}_surf.png", cmap="YlOrRd")
    if "contrast" in which:
        for f in sorted(glob.glob("04_results/tables/contrasts/*/*_diff.nii.gz")):
            name = os.path.basename(f).replace("_diff.nii.gz", "")
            surf4(f, f"{name} (red = first > second)", f"contrast_{name}_surf.png", cmap="RdBu_r", symmetric=True)
            slices(f, name, f"contrast_{name}_axial.png", cmap="RdBu_r")
