"""Meta-analytic connectivity gradients of the human hippocampus.

For every hippocampal voxel, the K experiments whose nearest hippocampal focus lies closest to the voxel are selected
(spatial filter; Chase et al., 2015). Their cortical co-activation profile is summarised over the Schaefer-400 parcels,
using only foci > 12 mm from the hippocampus (no kernel/proximity spill-over). The expected profile given each
experiment's number of foci (estimated from all non-hippocampal experiments of the same database) is subtracted,
giving a hippocampus-specific co-activation profile per voxel. Voxel-by-voxel affinity (normalised angle of cosine
similarity, row-wise top-10% sparsification) is decomposed with diffusion-map embedding (Margulies et al., 2016).

Usage: python 26_meta_gradient.py <db> <K> [split]   db in {pooled, neurosynth, nq_indep}; split in {half1, half2}
"""
import sys, os, json, time, numpy as np, pandas as pd, nibabel as nib
from scipy import ndimage
from scipy.spatial import cKDTree
from scipy.sparse.linalg import eigsh
sys.path.insert(0, "03_analysis")
from lib import *
from nilearn import datasets, image

db = sys.argv[1]; K = int(sys.argv[2]); split = sys.argv[3] if len(sys.argv) > 3 else None
NPERM = int(sys.argv[4]) if len(sys.argv) > 4 else 0
if split == "none": split = None
OUT = f"04_results/meta_gradient/{db}_K{K}" + (f"_{split}" if split else ""); os.makedirs(OUT, exist_ok=True)
t0 = time.time()
ds = pooled_dataset() if db == "pooled" else (nq_independent() if db == "nq_indep" else load(db))
C = ds.coordinates[["id", "x", "y", "z"]].copy()

# --- geometry
H_img = nib.load(f"{SEEDDIR}/BIL_HIPP.nii.gz"); H = np.asarray(H_img.dataobj) > 0; aff = H_img.affine; inv = np.linalg.inv(aff)
sch = datasets.fetch_atlas_schaefer_2018(n_rois=400, yeo_networks=7, resolution_mm=2)
S = np.asarray(image.resample_to_img(sch.maps, H_img, interpolation="nearest", force_resample=True, copy_header=True).dataobj).astype(int)
dist_to_H = ndimage.distance_transform_edt(~H, sampling=2.0)
ijk = np.round(nib.affines.apply_affine(inv, C[["x", "y", "z"]].values)).astype(int)
ok = np.all((ijk >= 0) & (ijk < np.array(H.shape)), axis=1); C, ijk = C[ok].reset_index(drop=True), ijk[ok]
C["dH"] = dist_to_H[tuple(ijk.T)]

# --- study x parcel binary activation (foci > 12 mm from hippocampus; focus counts toward parcels within 6 mm)
offs = np.array([(a, b, c) for a in range(-3, 4) for b in range(-3, 4) for c in range(-3, 4) if a * a + b * b + c * c <= 9])
studies = np.array(sorted(C.id.unique())); sidx = {s: i for i, s in enumerate(studies)}
far = C.dH.values > 12
fi = ijk[far]; fs = C.id.values[far]
rows, cols = [], []
for o in offs:
    q = fi + o; q = np.clip(q, 0, np.array(S.shape) - 1); lab = S[tuple(q.T)]
    m = lab > 0; rows.append(np.fromiter((sidx[s] for s in fs[m]), int)); cols.append(lab[m] - 1)
rows = np.concatenate(rows); cols = np.concatenate(cols)
A = np.zeros((len(studies), 400), np.float32); A[rows, cols] = 1
nfoci = C.groupby("id").size().reindex(studies).values

# --- hippocampal experiments: nearest-focus distance of each study to each hippocampal voxel
Hvox = np.argwhere(H); Hxyz = nib.affines.apply_affine(aff, Hvox)
near = C[C.dH <= 4]                                     # foci in or within 4 mm of the hippocampal mask
hip_ids = np.array(sorted(near.id.unique()))
if split:
    rng = np.random.default_rng(11); perm = rng.permutation(len(hip_ids)); half = perm[: len(perm) // 2] if split == "half1" else perm[len(perm) // 2:]
    hip_ids = np.sort(hip_ids[half]); near = near[near.id.isin(hip_ids)]
_hipall = set(C[C.dH <= 4].id); nonhip = np.array([s for s in studies if s not in _hipall])
# expected profile given number of foci (20 quantile bins of non-hippocampal studies)
nh_idx = np.array([sidx[s] for s in nonhip]); edges = np.unique(np.quantile(nfoci[nh_idx], np.linspace(0, 1, 21)))
binof = lambda n: np.clip(np.searchsorted(edges, n, side="right") - 1, 0, len(edges) - 2)
Ebin = np.vstack([A[nh_idx[binof(nfoci[nh_idx]) == b]].mean(0) for b in range(len(edges) - 1)])
tree = cKDTree(near[["x", "y", "z"]].values); near_ids = near.id.values
hid = {s: i for i, s in enumerate(hip_ids)}
D = np.full((len(Hvox), len(hip_ids)), np.inf, np.float32)
pairs = tree.query_ball_point(Hxyz, r=40)
for v, lst in enumerate(pairs):
    if not lst: continue
    d = np.linalg.norm(near[["x", "y", "z"]].values[lst] - Hxyz[v], axis=1)
    ids = near_ids[lst]
    dfv = pd.DataFrame({"i": [hid[s] for s in ids], "d": d}).groupby("i").d.min()
    D[v, dfv.index.values] = dfv.values
order = np.argsort(D, axis=1)[:, :K]
hip_sidx = np.array([sidx[s] for s in hip_ids])
prof = np.zeros((len(Hvox), 400), np.float32); obs = np.zeros_like(prof)
for v in range(len(Hvox)):
    si = hip_sidx[order[v]]
    obs[v] = A[si].mean(0); prof[v] = obs[v] - Ebin[binof(nfoci[si])].mean(0)
radius = np.take_along_axis(D, order[:, -1:], 1).ravel()

# --- diffusion map embedding
def diffusion_map(X, sparsity=0.9, n_comp=10, alpha=0.5):
    Xn = X / np.linalg.norm(X, axis=1, keepdims=True)
    cs = np.clip(Xn @ Xn.T, -1, 1); Aff = 1 - np.arccos(cs) / np.pi
    thr = np.quantile(Aff, sparsity, axis=1, keepdims=True); Aff = np.where(Aff >= thr, Aff, 0); Aff = (Aff + Aff.T) / 2
    d = Aff.sum(1); L = Aff / np.outer(d ** alpha, d ** alpha)
    d2 = L.sum(1); M = L / np.sqrt(np.outer(d2, d2))
    vals, vecs = np.linalg.eigh(M); idx = np.argsort(vals)[::-1]; vals, vecs = vals[idx], vecs[:, idx]
    psi = vecs / vecs[:, [0]]
    lam = vals[1:n_comp + 1]; emb = psi[:, 1:n_comp + 1] * (lam / (1 - lam))
    return emb, lam
emb, lam = diffusion_map(prof)
if NPERM:
    # long-axis coordinate = first principal axis of voxel coordinates (within hemisphere, anterior positive)
    LA = np.zeros(len(Hvox))
    for hm in [Hxyz[:, 0] < 0, Hxyz[:, 0] >= 0]:
        X = Hxyz[hm] - Hxyz[hm].mean(0); u = np.linalg.svd(X, full_matrices=False)[2][0]; u = u if u[1] > 0 else -u; LA[hm] = X @ u
    rng = np.random.default_rng(123); null = []
    for it in range(NPERM):
        perm = rng.permutation(len(hip_sidx)); ps = hip_sidx[perm]            # study location <-> cortical profile link broken
        P = np.zeros_like(prof)
        for v in range(len(Hvox)):
            si = ps[order[v]]; P[v] = A[si].mean(0) - Ebin[binof(nfoci[si])].mean(0)
        e, l = diffusion_map(P)
        null.append((l[0], l[0] / l.sum(), max(abs(np.corrcoef(e[:, c], LA)[0, 1]) for c in range(3))))
    null = np.array(null); e0 = emb
    obs_r = abs(np.corrcoef(e0[:, 0], LA)[0, 1])
    res = dict(n_perm=NPERM, obs_lambda1=float(lam[0]), obs_varexp1=float(lam[0] / lam.sum()), obs_abs_r_G1_longaxis=float(obs_r),
               p_lambda1=float((1 + (null[:, 0] >= lam[0]).sum()) / (NPERM + 1)), p_varexp1=float((1 + (null[:, 1] >= lam[0] / lam.sum()).sum()) / (NPERM + 1)),
               p_r_longaxis_maxG1to3=float((1 + (null[:, 2] >= obs_r).sum()) / (NPERM + 1)),
               null_lambda1_95=float(np.quantile(null[:, 0], .95)), null_varexp1_95=float(np.quantile(null[:, 1], .95)), null_r_95=float(np.quantile(null[:, 2], .95)))
    json.dump(res, open(f"{OUT}/permutation.json", "w"), indent=1); print(json.dumps(res, indent=1)); sys.exit()
# orient G1 so that it correlates positively with MNI y (posterior -> anterior)
for c in range(3):
    if np.corrcoef(emb[:, c], Hxyz[:, 1])[0, 1] < 0 and c == 0: emb[:, c] *= -1
np.savez_compressed(f"{OUT}/gradient_data.npz", emb=emb, lam=lam, prof=prof, obs=obs, Hvox=Hvox, Hxyz=Hxyz, radius=radius,
                    hip_ids=hip_ids, order=order, studies=studies, A=A, nfoci=nfoci, Ebin=Ebin, edges=edges)
for c in range(3):
    g = np.zeros(H.shape, np.float32); g[tuple(Hvox.T)] = emb[:, c]; nib.save(nib.Nifti1Image(g, aff), f"{OUT}/G{c+1}.nii.gz")
info = dict(db=db, K=K, split=split, n_hip_experiments=int(len(hip_ids)), n_voxels=int(len(Hvox)),
            median_filter_radius_mm=float(np.median(radius)), lambdas=[float(x) for x in lam[:5]],
            var_explained=[float(x) for x in (lam / lam.sum())[:5]],
            r_G1_y=float(np.corrcoef(emb[:, 0], Hxyz[:, 1])[0, 1]), r_G1_z=float(np.corrcoef(emb[:, 0], Hxyz[:, 2])[0, 1]),
            r_G1_absx=float(np.corrcoef(emb[:, 0], np.abs(Hxyz[:, 0]))[0, 1]),
            r_G2_y=float(np.corrcoef(emb[:, 1], Hxyz[:, 1])[0, 1]), r_G2_z=float(np.corrcoef(emb[:, 1], Hxyz[:, 2])[0, 1]),
            r_G2_absx=float(np.corrcoef(emb[:, 1], np.abs(Hxyz[:, 0]))[0, 1]), seconds=round(time.time() - t0))
json.dump(info, open(f"{OUT}/info.json", "w"), indent=1); print(json.dumps(info, indent=1))
