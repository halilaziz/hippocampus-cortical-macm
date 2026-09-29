"""Characterise meta-analytic hippocampal gradients: long axis, cortical projection, gradient-vs-parcel tests,
replication across databases/filters/halves, and functional decoding along G1."""
import sys, os, json, glob, numpy as np, pandas as pd, nibabel as nib
from scipy.stats import spearmanr, pearsonr
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score
from statsmodels.stats.multitest import multipletests
sys.path.insert(0, "03_analysis")
from importlib import import_module
nt = import_module("08_networks_tables")
OUT = "04_results/meta_gradient"; TAB = "04_results/tables/meta_gradient"; os.makedirs(TAB, exist_ok=True)
Z = {os.path.basename(os.path.dirname(f)): np.load(f, allow_pickle=True) for f in glob.glob(f"{OUT}/*/gradient_data.npz")}
P = Z["pooled_K150"]; emb, lam, prof, Hxyz = P["emb"], P["lam"], P["prof"], P["Hxyz"]
res = {}

# 1. long axis coordinate (first principal axis within each hemisphere; anterior positive)
LA = np.zeros(len(Hxyz))
for hm in [Hxyz[:, 0] < 0, Hxyz[:, 0] >= 0]:
    X = Hxyz[hm] - Hxyz[hm].mean(0); u = np.linalg.svd(X, full_matrices=False)[2][0]; u = u if u[1] > 0 else -u; LA[hm] = X @ u
    res.setdefault("long_axis_vectors", []).append([round(float(x), 3) for x in u])
# residual transverse axes: medial-lateral (|x|) and dorsal-ventral after removing long axis
res["r_G1_longaxis"] = float(pearsonr(emb[:, 0], LA)[0])
res["rho_G1_longaxis"] = float(spearmanr(emb[:, 0], LA)[0])
for c in range(1, 3):
    res[f"r_G{c+1}_longaxis"] = float(pearsonr(emb[:, c], LA)[0]); res[f"r_G{c+1}_absx"] = float(pearsonr(emb[:, c], np.abs(Hxyz[:, 0]))[0])
    res[f"r_G{c+1}_hemi"] = float(pearsonr(emb[:, c], np.sign(Hxyz[:, 0]))[0])
res["r_G1_left_vs_right_longaxis"] = [float(pearsonr(emb[Hxyz[:, 0] < 0, 0], LA[Hxyz[:, 0] < 0])[0]), float(pearsonr(emb[Hxyz[:, 0] >= 0, 0], LA[Hxyz[:, 0] >= 0])[0])]
res["var_explained"] = [float(x) for x in (lam / lam.sum())[:5]]

# 2. cortical projection: correlation of each parcel's specific co-activation with G1 across hippocampal voxels
r_parcel = np.array([pearsonr(emb[:, 0], prof[:, p])[0] if prof[:, p].std() > 0 else 0 for p in range(400)])
lab = [l.decode() if isinstance(l, bytes) else l for l in nt.sch.labels]
if lab[0].lower().startswith("background"): lab = lab[1:]
net = [nt.NETNAMES[l.split("_")[2]] for l in lab]
cp = pd.DataFrame({"parcel": lab, "network": net, "r_with_G1": r_parcel, "mean_specific_coact": prof.mean(0)})
cp.to_csv(f"{TAB}/G1_cortical_projection_parcels.csv", index=False)
netsum = cp.groupby("network").agg(mean_r=("r_with_G1", "mean"), pct_parcels_r_pos=("r_with_G1", lambda x: 100 * (x > 0.3).mean()),
                                   pct_parcels_r_neg=("r_with_G1", lambda x: 100 * (x < -0.3).mean()), mean_specific=("mean_specific_coact", "mean")).round(3)
netsum.to_csv(f"{TAB}/G1_cortical_projection_networks.csv"); print(netsum.to_string())
# volumetric parcel map of r for plotting
img = np.zeros(nt.S.shape, np.float32)
for p in range(400): img[nt.S == p + 1] = r_parcel[p]
nib.save(nib.Nifti1Image(img, nt.ref.affine), f"{OUT}/G1_cortical_projection_r.nii.gz")
print("top anterior-end parcels:", cp.sort_values("r_with_G1", ascending=False).head(10)[["parcel", "r_with_G1"]].values.tolist())
print("top posterior-end parcels:", cp.sort_values("r_with_G1").head(10)[["parcel", "r_with_G1"]].values.tolist())

# network profile along G1 (deciles)
dec = pd.qcut(emb[:, 0], 10, labels=False)
rows = []
for d in range(10):
    m = dec == d
    for n in sorted(set(net)):
        idx = [i for i, x in enumerate(net) if x == n]
        rows.append(dict(G1_decile=d + 1, G1_mean=float(emb[m, 0].mean()), longaxis_mm=float(LA[m].mean()), network=n, specific_coact=float(100 * prof[m][:, idx].mean())))
pd.DataFrame(rows).to_csv(f"{TAB}/G1_decile_network_profile.csv", index=False)

# 3. gradient vs parcels
Xn = prof / np.linalg.norm(prof, axis=1, keepdims=True)
sil = {k: float(silhouette_score(Xn, KMeans(k, n_init=20, random_state=0).fit_predict(Xn), metric="cosine")) for k in range(2, 8)}
res["silhouette_by_k"] = sil
try:
    import diptest
    dip, pdip = diptest.diptest(emb[:, 0]); res["dip_G1"] = float(dip); res["p_dip_G1"] = float(pdip)
except Exception as e:
    res["dip_error"] = str(e)
# cross-validated prediction of voxel profiles from G1: step function (k bins, k-means on G1) vs smooth polynomial (same #params)
from sklearn.model_selection import GroupKFold
groups = np.digitize(LA, np.quantile(LA, np.linspace(0, 1, 11)[1:-1]))   # spatial blocks along the long axis
def cv_r2(kind, k):
    ss_res = ss_tot = 0.0
    for tr, te in GroupKFold(n_splits=5).split(prof, groups=groups):
        g_tr, g_te = emb[tr, 0], emb[te, 0]
        if kind == "step":
            km = KMeans(k, n_init=20, random_state=0).fit(g_tr[:, None]); lt, le = km.labels_, km.predict(g_te[:, None])
            means = np.vstack([prof[tr][lt == j].mean(0) for j in range(k)]); pred = means[le]
        else:
            V = lambda g: np.vander(g, k, increasing=True); B = np.linalg.lstsq(V(g_tr), prof[tr], rcond=None)[0]; pred = V(g_te) @ B
        ss_res += ((prof[te] - pred) ** 2).sum(); ss_tot += ((prof[te] - prof[tr].mean(0)) ** 2).sum()
    return 1 - ss_res / ss_tot
res["cv_r2"] = {f"{kind}_k{k}": round(float(cv_r2(kind, k)), 4) for k in [2, 3, 4, 5] for kind in ["step", "smooth"]}
print("silhouette", sil); print("cv_r2", res["cv_r2"])

# 4. replication of G1 across analyses (voxelwise r)
rep = []
for k2, Z2 in Z.items():
    if k2 == "pooled_K150": continue
    e2 = Z2["emb"]; r = max((pearsonr(emb[:, 0], s * e2[:, c])[0], c) for c in range(3) for s in [1, -1])
    rep.append(dict(analysis=k2, n_experiments=int(len(Z2["hip_ids"])), matched_component=f"G{r[1]+1}", r_with_pooled_G1=round(float(r[0]), 3),
                    r_G1_longaxis=round(float(abs(pearsonr(e2[:, 0], LA)[0])), 3)))
if "neurosynth_K150" in Z and "nq_indep_K100" in Z:
    a, b = Z["neurosynth_K150"]["emb"][:, 0], Z["nq_indep_K100"]["emb"][:, 0]
    res["r_G1_neurosynth_vs_nqindep"] = float(abs(pearsonr(a, b)[0]))
if "pooled_K150_half1" in Z and "pooled_K150_half2" in Z:
    res["r_G1_half1_vs_half2"] = float(abs(pearsonr(Z["pooled_K150_half1"]["emb"][:, 0], Z["pooled_K150_half2"]["emb"][:, 0])[0]))
pd.DataFrame(rep).to_csv(f"{TAB}/G1_replication.csv", index=False); print(pd.DataFrame(rep).to_string(index=False))

# 5. functional decoding along G1 (Neurosynth studies with hippocampal foci; study G1 = mean G1 at its hippocampal foci)
from lib import load
ns = load("neurosynth"); ann = ns.annotations.set_index("id")
terms = [c for c in ann.columns if c.startswith("terms_abstract_tfidf__")]
Xt = (ann[terms] > 0.001); Xt.columns = [c.split("__")[1] for c in terms]
Hv = {tuple(v): i for i, v in enumerate(P["Hvox"])}
inv = np.linalg.inv(nt.ref.affine); c = ns.coordinates
ijk = np.round(nib.affines.apply_affine(inv, c[["x", "y", "z"]].values)).astype(int)
g_of = [emb[Hv[tuple(v)], 0] if tuple(v) in Hv else np.nan for v in ijk]
c = c.assign(G1=g_of).dropna(subset=["G1"]); sg = c.groupby("id").G1.mean()
Xs = Xt.loc[sg.index]; Xs = Xs.loc[:, Xs.sum() >= 30]
drows = []
for t in Xs.columns:
    y = Xs[t].values.astype(float); r, p = pearsonr(sg.values, y); drows.append((t, int(y.sum()), r, p, sg[y == 1].mean(), sg[y == 0].mean()))
dd = pd.DataFrame(drows, columns=["term", "n_studies", "r_pointbiserial", "p", "mean_G1_with", "mean_G1_without"])
dd["q_fdr"] = multipletests(dd.p, method="fdr_bh")[1]; dd = dd.sort_values("r_pointbiserial")
dd.to_csv(f"{TAB}/G1_term_decoding.csv", index=False)
res["n_decoding_studies"] = int(len(sg)); res["n_terms_q05"] = int((dd.q_fdr < 0.05).sum())
print("posterior-end terms:", dd[dd.q_fdr < 0.05].head(15)[["term", "r_pointbiserial"]].values.tolist())
print("anterior-end terms:", dd[dd.q_fdr < 0.05].tail(15)[["term", "r_pointbiserial"]].values.tolist())
json.dump(res, open(f"{TAB}/gradient_summary.json", "w"), indent=1); print(json.dumps({k: v for k, v in res.items() if k != "long_axis_vectors"}, indent=1))
