"""Anatomical labelling, cluster tables, cortical network profiles (Yeo-7 via Schaefer-400/7Net) for thresholded maps."""
import sys, os, glob, json, numpy as np, pandas as pd, nibabel as nib
from nilearn import datasets, image
from nilearn.reporting import get_clusters_table
from scipy import ndimage
sys.path.insert(0, "03_analysis")
VFWE = "z_level-voxel_corr-FWE_method-montecarlo.nii.gz"
Z_CRIT = 1.6449  # one-sided p<.05 FWE
MIN_VOX = 10     # 80 mm3 minimum cluster extent for tables
ref = nib.load(glob.glob("04_results/macm/*/*/z.nii.gz")[0])
sch = datasets.fetch_atlas_schaefer_2018(n_rois=400, yeo_networks=7, resolution_mm=2)
schimg = image.resample_to_img(sch.maps, ref, interpolation="nearest", force_resample=True, copy_header=True)
S = np.asarray(schimg.dataobj).astype(int)
labels = [l.decode() if isinstance(l, bytes) else l for l in sch.labels]
if labels[0].lower().startswith("background"): labels = labels[1:]
NETS = ["Vis", "SomMot", "DorsAttn", "SalVentAttn", "Limbic", "Cont", "Default"]
NETNAMES = {"Vis":"Visual","SomMot":"Somatomotor","DorsAttn":"Dorsal attention","SalVentAttn":"Salience/ventral attention","Limbic":"Limbic","Cont":"Frontoparietal control","Default":"Default mode"}
parcel_net = {i+1: l.split("_")[2] for i, l in enumerate(labels)}
parcel_hemi = {i+1: l.split("_")[1] for i, l in enumerate(labels)}
hoc = datasets.fetch_atlas_harvard_oxford("cort-maxprob-thr25-2mm"); hos = datasets.fetch_atlas_harvard_oxford("sub-maxprob-thr25-2mm")
HOC = np.asarray(image.resample_to_img(hoc.maps, ref, interpolation="nearest", force_resample=True, copy_header=True).dataobj).astype(int)
HOS = np.asarray(image.resample_to_img(hos.maps, ref, interpolation="nearest", force_resample=True, copy_header=True).dataobj).astype(int)
aff = ref.affine; inv = np.linalg.inv(aff)
BAD = ("White Matter", "Cerebral Cortex", "Ventricle", "Brain-Stem")
def _lookup(i, j, k):
    c, s_ = HOC[i, j, k], HOS[i, j, k]
    if c > 0: return hoc.labels[c]
    if s_ > 0 and not any(b in hos.labels[s_] for b in BAD): return hos.labels[s_]
    return ""
def label_xyz(x, y, z):
    i, j, k = np.round(inv @ [x, y, z, 1])[:3].astype(int)
    lab = _lookup(i, j, k); near = ""
    if not lab:   # nearest grey-matter label within 3 voxels (6 mm)
        best = None
        for di in range(-3, 4):
            for dj in range(-3, 4):
                for dk in range(-3, 4):
                    d = di*di + dj*dj + dk*dk
                    if d > 9: continue
                    l = _lookup(i+di, j+dj, k+dk)
                    if l and (best is None or d < best[0]): best = (d, l)
        if best: lab, near = best[1], " (nearest)"
    hemi = "L" if x < 0 else "R"
    if not lab and y < -30 and z < -15: lab = "Cerebellum"
    if lab and not lab.startswith(("Left", "Right")): lab = f"{hemi} {lab}"
    lab = lab.replace("Left ", "L ").replace("Right ", "R ")
    p = S[i, j, k]
    return lab + near, (NETNAMES[parcel_net[p]] if p > 0 else "")
def network_profile(mask):
    cort = S > 0
    rows = []
    tot = (mask & cort).sum()
    for n in NETS:
        netvox = np.isin(S, [k for k, v in parcel_net.items() if v == n])
        ov = (mask & netvox).sum()
        rows.append(dict(network=NETNAMES[n], vox_in_map=int(ov), pct_of_cortical_map=100*ov/max(tot,1), pct_of_network_covered=100*ov/netvox.sum()))
    return pd.DataFrame(rows), int(tot)
def parcel_profile(mask, zmap):
    rows = []
    for p in range(1, len(labels)+1):
        pv = S == p; ov = (mask & pv).sum()
        if ov: rows.append(dict(parcel=labels[p-1], network=NETNAMES[parcel_net[p]], pct_parcel=100*ov/pv.sum(), mean_z=float(zmap[mask & pv].mean())))
    return pd.DataFrame(rows).sort_values("pct_parcel", ascending=False)
def process(sigfile, outdir, name, zcrit=Z_CRIT, min_vox=MIN_VOX, statfile=None, alefile=None):
    """sigfile: map defining significance (voxel-FWE z); statfile: map used for peaks (ALE z)."""
    os.makedirs(outdir, exist_ok=True)
    img = nib.load(sigfile); zs = np.asarray(img.dataobj).astype(float)
    m = zs > zcrit
    lab, n = ndimage.label(m, structure=np.ones((3,3,3)))
    sizes = ndimage.sum(m, lab, range(1, n+1))
    keep = np.isin(lab, [i+1 for i, s in enumerate(sizes) if s >= min_vox])
    z = np.asarray(nib.load(statfile).dataobj).astype(float) if statfile else zs
    z2 = np.where(keep, z, 0); nib.save(nib.Nifti1Image(z2.astype(np.float32), img.affine), f"{outdir}/{name}_thr.nii.gz")
    tab = get_clusters_table(nib.Nifti1Image(z2, img.affine), stat_threshold=1e-6, cluster_threshold=min_vox, min_distance=12, two_sided=False)
    if len(tab):
        lbl = [label_xyz(r.X, r.Y, r.Z) for r in tab.itertuples()]
        tab["Region (Harvard-Oxford)"] = [l[0] for l in lbl]; tab["Yeo-7 network"] = [l[1] for l in lbl]
        if alefile:
            A = np.asarray(nib.load(alefile).dataobj)
            tab["ALE"] = [float(A[tuple(np.round(inv @ [r.X, r.Y, r.Z, 1])[:3].astype(int))]) for r in tab.itertuples()]
        tab = tab.rename(columns={"Peak Stat": "Z-score"})
    tab.to_csv(f"{outdir}/{name}_clusters.csv", index=False)
    prof, ncort = network_profile(keep); prof.to_csv(f"{outdir}/{name}_networks.csv", index=False)
    parcel_profile(keep, z).to_csv(f"{outdir}/{name}_parcels.csv", index=False)
    summ = dict(name=name, n_sig_vox=int(keep.sum()), vol_mm3=int(keep.sum()*8), n_cortical_vox=ncort,
                n_clusters=int(sum(s >= min_vox for s in sizes)))
    json.dump(summ, open(f"{outdir}/{name}_summary.json", "w"))
    return summ, prof
if __name__ == "__main__" and len(sys.argv) == 1:
    rows = []
    for d in sorted(glob.glob("04_results/macm/*/*/info.json")):
        d = os.path.dirname(d); db, seed = d.split("/")[-2:]
        s, prof = process(f"{d}/{VFWE}", f"04_results/tables/macm/{db}", seed, statfile=f"{d}/z.nii.gz", alefile=f"{d}/stat.nii.gz")
        info = json.load(open(f"{d}/info.json"))
        rows.append({**info, **s, **{f"pct_{r.network}": round(r.pct_of_cortical_map, 1) for r in prof.itertuples()}})
    pd.DataFrame(rows).to_csv("04_results/tables/macm_summary.csv", index=False)
    print(pd.DataFrame(rows)[["db","seed","n_experiments","n_foci","n_sig_vox","n_clusters","n_cortical_vox"] + [c for c in rows[0] if c.startswith("pct_")]].to_string(index=False))

def run_specificity(q=0.01):
    rows = []
    for d in sorted(glob.glob("04_results/specificity/*/*/info.json")):
        d = os.path.dirname(d); db, seed = d.split("/")[-2:]
        z = nib.load(f"{d}/z_desc-association.nii.gz"); zz = np.asarray(z.dataobj)
        p = np.asarray(nib.load(f"{d}/p_desc-association_level-voxel_corr-FDR_method-indep.nii.gz").dataobj)
        sig = np.where((p < q) & (zz > 0), zz, 0).astype(np.float32)
        tmp = f"{d}/sig_pos_q{q}.nii.gz"; nib.save(nib.Nifti1Image(sig, z.affine), tmp)
        s, prof = process(tmp, f"04_results/tables/specificity/{db}", seed, zcrit=1e-6, statfile=f"{d}/z_desc-association.nii.gz")
        info = json.load(open(f"{d}/info.json"))
        rows.append({**info, **s, **{f"pct_{r.network}": round(r.pct_of_cortical_map, 1) for r in prof.itertuples()}})
    df = pd.DataFrame(rows); df.to_csv("04_results/tables/specificity_summary.csv", index=False); print(df.drop(columns=["name"]).to_string(index=False))

def run_contrasts(p_unc=0.001, min_vox=25):
    """ALE subtraction: voxels with p<.001 (permutation) and >=200 mm3, masked by the minuend's thresholded MACM (GingerALE convention)."""
    from scipy.stats import norm
    zc = norm.isf(p_unc)
    for d in sorted(glob.glob("04_results/contrasts/*/*/info.json")):
        info = json.load(open(d)); d = os.path.dirname(d); db = d.split("/")[-2]
        z = nib.load(f"{d}/z_desc-group1MinusGroup2.nii.gz"); zz = np.asarray(z.dataobj).astype(float)
        mA = np.asarray(nib.load(f"04_results/tables/macm/{db}/{info['a']}_thr.nii.gz").dataobj) > 0
        mB = np.asarray(nib.load(f"04_results/tables/macm/{db}/{info['b']}_thr.nii.gz").dataobj) > 0
        out = f"04_results/tables/contrasts/{db}"; os.makedirs(out, exist_ok=True)
        diff = np.zeros_like(zz)
        for name, arr, mask, sign in [(f"{info['a']}_gt_{info['b']}", zz, mA, 1), (f"{info['b']}_gt_{info['a']}", -zz, mB, -1)]:
            sig = np.where((arr > zc) & mask, arr, 0).astype(np.float32)
            tmp = f"{d}/{name}_sig.nii.gz"; nib.save(nib.Nifti1Image(sig, z.affine), tmp)
            s, prof = process(tmp, out, name, zcrit=1e-6, min_vox=min_vox)
            thr = np.asarray(nib.load(f"{out}/{name}_thr.nii.gz").dataobj)
            diff += sign * thr
            print(name, s); print(prof.round(1).to_string(index=False))
        nib.save(nib.Nifti1Image(diff.astype(np.float32), z.affine), f"{out}/{info['a']}_vs_{info['b']}_diff.nii.gz")

if __name__ == "__main__" and len(sys.argv) > 1:
    if "spec" in sys.argv: run_specificity()
    if "contrast" in sys.argv: run_contrasts()
