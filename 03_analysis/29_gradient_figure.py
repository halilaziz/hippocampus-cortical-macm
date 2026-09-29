"""Figure 5: meta-analytic hippocampal connectivity gradient."""
import sys, json, numpy as np, pandas as pd, nibabel as nib, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from nilearn import plotting, datasets, surface
sys.path.insert(0, "03_analysis")
exec(open("03_analysis/28_gradient_content_null.py").read().split("NPERM =")[0])   # helpers: load, longaxis, stats, net, NETS
plt.rcParams.update({"font.family": "Arial", "font.size": 9})
NETCOL = {"Visual": "#781C86", "Somatomotor": "#4682B4", "Dorsal attention": "#00760E", "Salience/ventral attention": "#C43AFA",
          "Limbic": "#b8a44a", "Frontoparietal control": "#E69422", "Default mode": "#CD3E4E"}
Zp, hp = load("pooled_K150"); H = Zp["Hxyz"]; emb = Zp["emb"]; LA = longaxis(H); g1 = emb[:, 0]
proj, _ = stats(Zp["prof"], LA)
fig = plt.figure(figsize=(13, 11)); gs = fig.add_gridspec(3, 4, height_ratios=[1, 1, 1], hspace=0.45, wspace=0.35)
# A: G1 on hippocampus (sagittal projection y-z), left and right
for k, (hm, name) in enumerate([(H[:, 0] < 0, "Left"), (H[:, 0] >= 0, "Right")]):
    ax = fig.add_subplot(gs[0, k]); o = np.argsort(np.abs(H[hm, 0]))
    sc = ax.scatter(H[hm, 1][o], H[hm, 2][o], c=g1[hm][o], cmap="Spectral_r", s=22, marker="s", edgecolors="none")
    ax.set_aspect("equal"); ax.set_xlabel("MNI y (mm)"); ax.set_ylabel("MNI z (mm)"); ax.set_title(f"{'A  ' if k == 0 else ''}G1, {name} hippocampus (sagittal view)", loc="left", fontweight="bold" if k == 0 else None, fontsize=10)
    ax.spines[["top", "right"]].set_visible(False)
plt.colorbar(sc, ax=ax, fraction=0.04, label="G1 (a.u.)")
# B: G1 vs long axis
ax = fig.add_subplot(gs[0, 2:])
ax.scatter(LA, g1, c=np.where(H[:, 0] < 0, "#4c72b0", "#dd8452"), s=6, alpha=0.6)
ax.set_xlabel("Long-axis position (mm; posterior → anterior)"); ax.set_ylabel("G1")
perm = json.load(open("04_results/meta_gradient/pooled_K150/permutation.json"))
ax.set_title("B  G1 follows the long axis", loc="left", fontweight="bold", fontsize=10)
ax.text(0.02, 0.95, f"r = {np.corrcoef(LA, g1)[0,1]:.2f} (geometric null 95th pct: {perm['null_r_95']:.2f})\nG1 variance {100*perm['obs_varexp1']:.1f}% vs null 95th pct {100*perm['null_varexp1_95']:.1f}% (p = {perm['p_varexp1']:.3f})",
        transform=ax.transAxes, va="top", fontsize=8)
ax.spines[["top", "right"]].set_visible(False)
# C: cortical projection on surface
fs = datasets.fetch_surf_fsaverage("fsaverage5")
img = "04_results/meta_gradient/G1_cortical_projection_r.nii.gz"
data = {h: surface.vol_to_surf(img, fs[f"pial_{h}"], radius=2, interpolation="nearest_most_frequent", kind="line") for h in ["left", "right"]}
for c, (h, v) in enumerate([("left", "lateral"), ("left", "medial"), ("right", "medial"), ("right", "lateral")]):
    ax = fig.add_subplot(gs[1, c], projection="3d")
    plotting.plot_surf_stat_map(fs[f"infl_{h}"], data[h], hemi=h, view=v, bg_map=fs[f"sulc_{h}"], axes=ax, cmap="RdBu_r", vmax=0.7, threshold=0.1, colorbar=(c == 3), symmetric_cbar=True, cbar_tick_format="%.1f")
    if c == 0: ax.set_title("C  Cortical projection of G1 (red: anterior end, blue: posterior end)", loc="left", fontweight="bold", fontsize=10, x=0)
# D: network–long-axis association with null interval
nd = pd.read_csv("04_results/tables/meta_gradient/network_longaxis_permutation.csv").sort_values("r_longaxis_vs_specific_coact")
ax = fig.add_subplot(gs[2, :2]); y = np.arange(len(nd))
ax.barh(y, nd.r_longaxis_vs_specific_coact, color=[NETCOL[n] for n in nd.network])
for i, r in enumerate(nd.itertuples()): ax.plot([r.null_2p5, r.null_97p5], [i, i], color="k", lw=3, alpha=0.25, solid_capstyle="butt")
ax.set_yticks(y); ax.set_yticklabels(nd.network); ax.axvline(0, c="k", lw=0.6)
ax.set_xlabel("r (long-axis position × specific co-activation)\ngrey: 95% permutation null"); ax.set_title("D  Network content along the long axis", loc="left", fontweight="bold", fontsize=10)
ax.spines[["top", "right"]].set_visible(False)
# E: replication of projection pattern
Zn, _ = load("neurosynth_K150"); Zq, _ = load("nq_indep_K100")
pn, _ = stats(Zn["prof"], longaxis(Zn["Hxyz"])); pq, _ = stats(Zq["prof"], longaxis(Zq["Hxyz"]))
rr = json.load(open("04_results/tables/meta_gradient/projection_replication_null.json"))
ax = fig.add_subplot(gs[2, 2]); ax.scatter(pn, pq, c=[NETCOL[n] for n in net], s=8)
ax.set_xlabel("Neurosynth projection (r)"); ax.set_ylabel("Independent NeuroQuery (r)")
ax.set_title("E  Independent replication", loc="left", fontweight="bold", fontsize=10)
ax.text(0.97, 0.03, f"r = {rr['r_projection_neurosynth_vs_nqindep']:.2f}\nnull 99th pct = {rr['null_pairwise_projection_r_99']:.2f}", transform=ax.transAxes, va="bottom", ha="right", fontsize=8, bbox=dict(fc="white", ec="none", alpha=0.8))
ax.spines[["top", "right"]].set_visible(False)
# F: term decoding
td = pd.read_csv("04_results/tables/meta_gradient/G1_term_decoding.csv")
top = pd.concat([td.sort_values("r_pointbiserial").head(6), td.sort_values("r_pointbiserial").tail(6)])
ax = fig.add_subplot(gs[2, 3]); yy = np.arange(len(top))
ax.barh(yy, top.r_pointbiserial, color=["#4575b4" if v < 0 else "#d73027" for v in top.r_pointbiserial])
ax.set_yticks(yy); ax.set_yticklabels([f"{t.replace('autobiographical memory', 'autobiogr. memory')}{'*' if q < .05 else ''}" for t, q in zip(top.term, top.q_fdr)], fontsize=7.5); ax.axvline(0, c="k", lw=0.6)
ax.set_xlabel("r with study G1"); ax.set_title("F  Terms along G1", loc="left", fontweight="bold", fontsize=10)
ax.spines[["top", "right"]].set_visible(False)
plt.savefig("04_results/figures/pub/Fig5_meta_gradient.png", dpi=300, bbox_inches="tight"); print("saved")
