"""Composite publication figures."""
import os, json, glob, numpy as np, pandas as pd, nibabel as nib, matplotlib
import os
RSV = "v4" if os.path.exists("04_results/rsfc_ale/gingerale/rsfc_main_v4_p001_C05_1k_Z.nii") else "v3"
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
from nilearn import plotting, datasets, surface
FIG = "04_results/figures/pub"; os.makedirs(FIG, exist_ok=True)
fs = datasets.fetch_surf_fsaverage("fsaverage5")
NETCOL = {"Visual": "#781C86", "Somatomotor": "#4682B4", "Dorsal attention": "#00760E", "Salience/ventral attention": "#C43AFA",
          "Limbic": "#b8a44a", "Frontoparietal control": "#E69422", "Default mode": "#CD3E4E"}
NETS = list(NETCOL)
plt.rcParams.update({"font.family": "Arial", "font.size": 9})
def surf_row(fig, gs_row, img, cmap, vmax, title, symmetric=False, ylabel=None, threshold=1e-3):
    views = [("left", "lateral"), ("left", "medial"), ("right", "medial"), ("right", "lateral")]
    data = {h: surface.vol_to_surf(img, fs[f"pial_{h}"], radius=3, interpolation="linear", kind="line") for h in ["left", "right"]}
    axes = []
    for c, (h, v) in enumerate(views):
        ax = fig.add_subplot(gs_row[c], projection="3d")
        plotting.plot_surf_stat_map(fs[f"infl_{h}"], data[h], hemi=h, view=v, bg_map=fs[f"sulc_{h}"], axes=ax, cmap=cmap,
                                    vmax=vmax, threshold=threshold, colorbar=(c == 3), symmetric_cbar=symmetric)
        axes.append(ax)
    axes[0].set_title(title, loc="left", fontsize=10, fontweight="bold", x=0.0)
    return axes
def fig_workflow():
    fig, ax = plt.subplots(figsize=(10, 5.2)); ax.axis("off")
    P = json.load(open("01_literature/prisma_counts.json"))
    V2 = json.load(open("01_literature/prisma_counts_v2.json")) if os.path.exists("01_literature/prisma_counts_v2.json") else None
    n_inc = (V2["included"] + V2["included_other_methods"]) if V2 else P["included"]["total_included"]
    boxes = [
        (0.02, 0.62, 0.30, 0.33, f"Systematic literature search\nPubMed + Europe PMC (24 Sep 2026)\n{P['identification']['records_identified_total']:,} records → {P['screening']['records_screened']:,} screened\ndual screening; {n_inc} included\nafter full-text assessment", "#e8eef7"),
        (0.35, 0.62, 0.30, 0.33, "Coordinate databases\nNeurosynth v7: 14,371 studies\nNeuroQuery: 13,459 studies\nPooled (PMID-deduplicated): 19,148 studies,\n660,152 foci", "#e8f4ea"),
        (0.68, 0.62, 0.30, 0.33, "Seeds (MNI152, Harvard-Oxford)\nL/R/bilateral hippocampus\nanterior/posterior (uncal apex y = −21)\n6 coronal long-axis bins (6 mm)", "#fbeee6"),
        (0.02, 0.12, 0.205, 0.38, "ALE-MACM\n(GingerALE algorithm)\nconsistent co-activation\nvoxel-FWE p<.05,\n1,000 permutations", "#f3f3f3"),
        (0.26, 0.12, 0.205, 0.38, "Specific co-activation\nMKDA χ² vs\nfoci-matched reference\nFDR q<.01", "#f3f3f3"),
        (0.50, 0.12, 0.205, 0.38, "Contrasts & gradient\nALE subtraction\n(10,000 permutations)\nlong-axis network\nprofiles", "#f3f3f3"),
        (0.74, 0.12, 0.24, 0.38, "Validation\nindependent replication\n(non-overlapping NeuroQuery)\nGingerALE 3.0.2 equivalence\nsensitivity analyses", "#f3f3f3")]
    for x, y, w, h, t, c in boxes:
        ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.01", fc=c, ec="#555", lw=1))
        ax.text(x + w / 2, y + h / 2, t, ha="center", va="center", fontsize=8.5)
    for x in [0.17, 0.5, 0.83]:
        ax.annotate("", xy=(0.5, 0.52), xytext=(x, 0.61), arrowprops=dict(arrowstyle="->", color="#555"))
    plt.savefig(f"{FIG}/Fig1_workflow.png", dpi=300, bbox_inches="tight"); plt.close()
def netbars(ax, files, labels, title):
    w = 0.8 / len(files); x = np.arange(len(NETS))
    for i, (f, l) in enumerate(zip(files, labels)):
        d = pd.read_csv(f).set_index("network").reindex(NETS)
        ax.bar(x + i * w - 0.4 + w / 2, d.pct_of_cortical_map, w, label=l, color=[NETCOL[n] for n in NETS], alpha=1 - 0.45 * i, edgecolor="k", lw=0.4, hatch=["", "///", "..."][i])
    ax.set_xticks(x); ax.set_xticklabels([n.replace("/", "/\n").replace(" ", "\n", 1) for n in NETS], fontsize=7.5)
    ax.set_ylabel("% of significant cortical voxels"); ax.set_title(title, loc="left", fontweight="bold", fontsize=10)
    ax.spines[["top", "right"]].set_visible(False); ax.legend(frameon=False, fontsize=8)
def fig_main(db="pooled"):
    fig = plt.figure(figsize=(12, 9.5)); gs = fig.add_gridspec(3, 4, height_ratios=[1, 1, 1.1], hspace=0.25)
    surf_row(fig, [gs[0, i] for i in range(4)], f"04_results/tables/macm/{db}/BIL_HIPP_thr.nii.gz", "hot", 20, "A  Consistent co-activation (ALE-MACM, voxel-FWE p<.05)")
    surf_row(fig, [gs[1, i] for i in range(4)], f"04_results/tables/specificity/{db}/BIL_HIPP_thr.nii.gz", "YlOrRd", 10, "B  Specific co-activation (vs foci-matched reference, FDR q<.01)")
    ax = fig.add_subplot(gs[2, :])
    netbars(ax, [f"04_results/tables/macm/{db}/BIL_HIPP_networks.csv", f"04_results/tables/specificity/{db}/BIL_HIPP_networks.csv"],
            ["Consistent (ALE-MACM)", "Specific (MKDA χ², matched)"], "C  Yeo-7 network composition of hippocampal cortical connectivity")
    plt.savefig(f"{FIG}/Fig2_bilateral_{db}.png", dpi=300, bbox_inches="tight"); plt.close()
def fig_axis(db="pooled"):
    fig = plt.figure(figsize=(12, 12)); gs = fig.add_gridspec(4, 4, height_ratios=[1, 1, 1, 1.2], hspace=0.3)
    surf_row(fig, [gs[0, i] for i in range(4)], f"04_results/tables/specificity/{db}/BIL_aHIPP_thr.nii.gz", "YlOrRd", 10, "A  Anterior hippocampus: specific co-activation")
    surf_row(fig, [gs[1, i] for i in range(4)], f"04_results/tables/specificity/{db}/BIL_pHIPP_thr.nii.gz", "YlGnBu", 10, "B  Posterior hippocampus: specific co-activation")
    f = f"04_results/tables/contrasts/{db}/BIL_aHIPP_vs_BIL_pHIPP_diff.nii.gz"
    if os.path.exists(f): surf_row(fig, [gs[2, i] for i in range(4)], f, "RdBu_r", 4, "C  ALE subtraction (Z): anterior > posterior (red) / posterior > anterior (blue)", symmetric=True, threshold=1.0)
    g = pd.read_csv(f"04_results/tables/gradient_{db}.csv")
    for k, (col, lab) in enumerate([("pct_network_sig", "% network voxels significant"), ("delta_p", "ΔP activation vs matched ref. (%)")]):
        ax = fig.add_subplot(gs[3, 2 * k:2 * k + 2])
        for n, gg in g.groupby("network"):
            gg = gg.sort_values("y_center"); ax.plot(gg.y_center, gg[col], "-o", color=NETCOL[n], label=n, lw=2, ms=4)
        ax.axvline(-21, ls="--", c="grey", lw=1); ax.axhline(0, c="k", lw=0.4); ax.set_xlabel("Seed position along long axis, MNI y (mm)  (posterior → anterior)")
        ax.set_ylabel(lab); ax.spines[["top", "right"]].set_visible(False)
        ax.set_title("D  Long-axis connectivity gradient" if k == 0 else "", loc="left", fontweight="bold", fontsize=10)
        if k == 1: ax.legend(fontsize=7, frameon=False, ncol=2)
    plt.savefig(f"{FIG}/Fig3_long_axis_{db}.png", dpi=300, bbox_inches="tight"); plt.close()
def fig_replication():
    r = pd.read_csv("04_results/tables/replication.csv")
    fig = plt.figure(figsize=(12, 7)); gs = fig.add_gridspec(2, 4, height_ratios=[1, 1.1], hspace=0.3)
    surf_row(fig, [gs[0, i] for i in range(4)], "04_results/tables/conjunction/specificity_BIL_HIPP_conj.nii.gz", "YlOrRd", 8, "A  Conjunction: Neurosynth ∩ independent NeuroQuery (specific co-activation)")
    for k, kind in enumerate(["macm", "specificity"]):
        ax = fig.add_subplot(gs[1, 2 * k:2 * k + 2]); d = r[r.analysis == kind]
        x = np.arange(len(d)); w = 0.27
        ax.bar(x - w, d.dice, w, label="Dice (thresholded)", color="#4c72b0"); ax.bar(x, d.pct_nq_in_ns / 100, w, label="Replication voxels in discovery map", color="#55a868"); ax.bar(x + w, d.r_unthresholded, w, label="r (unthresholded)", color="#dd8452")
        ax.set_xticks(x); ax.set_xticklabels(d.seed, rotation=45, fontsize=7.5); ax.set_ylim(0, 1.25); ax.set_yticks([0, .2, .4, .6, .8, 1]); ax.spines[["top", "right"]].set_visible(False)
        ax.set_title(f"{'B' if k == 0 else 'C'}  Replication: {'ALE-MACM' if kind == 'macm' else 'specific co-activation'}", loc="left", fontweight="bold", fontsize=10)
        ax.legend(frameon=False, fontsize=7, ncol=3, loc="upper center")
    plt.savefig(f"{FIG}/Fig4_replication.png", dpi=300, bbox_inches="tight"); plt.close()
def fig_rsfc():
    """Figure 6: GingerALE meta-analysis of published hippocampal-seed rsFC studies vs MACM specific co-activation."""
    import nibabel as nib
    fig = plt.figure(figsize=(12, 7.5)); gs = fig.add_gridspec(2, 4, hspace=0.25)
    surf_row(fig, [gs[0, i] for i in range(4)], f"04_results/rsfc_ale/gingerale/rsfc_main_{RSV}_sig_resampled.nii.gz", "hot", 0.045,
             "A  Resting-state FC meta-analysis of published hippocampal-seed studies (GingerALE, cluster-FWE p<.05)")
    a = np.asarray(nib.load(f"04_results/rsfc_ale/gingerale/rsfc_main_{RSV}_sig_resampled.nii.gz").dataobj) > 0
    b = np.asarray(nib.load("04_results/tables/specificity/pooled/BIL_HIPP_thr.nii.gz").dataobj) > 0
    ov = np.zeros(a.shape, np.float32); ov[b] = 1; ov[a & b] = 3; ov[a & ~b] = 2
    img = nib.Nifti1Image(ov, nib.load("04_results/tables/specificity/pooled/BIL_HIPP_thr.nii.gz").affine)
    from matplotlib.colors import ListedColormap
    from matplotlib.patches import Patch
    cmap = ListedColormap(["#fdd49e", "#3182bd", "#de2d26"])
    dat = {h: surface.vol_to_surf(img, fs[f"pial_{h}"], radius=2, interpolation="nearest_most_frequent", kind="line") for h in ["left", "right"]}
    for c, (h, v) in enumerate([("left", "lateral"), ("left", "medial"), ("right", "medial"), ("right", "lateral")]):
        ax = fig.add_subplot(gs[1, c], projection="3d")
        roi = np.round(dat[h]).astype(float); roi[roi < 1] = np.nan
        plotting.plot_surf_roi(fs[f"infl_{h}"], roi, hemi=h, view=v, bg_map=fs[f"sulc_{h}"], axes=ax, cmap=cmap, vmin=1, vmax=3, colorbar=False)
        if c == 0: ax.set_title("B  Overlap with MACM specific co-activation", loc="left", fontweight="bold", fontsize=10, x=0)
    fig.legend(handles=[Patch(color="#fdd49e", label="MACM specific only"), Patch(color="#3182bd", label="rsFC meta-analysis only"), Patch(color="#de2d26", label="Both")],
               loc="lower center", ncol=3, frameon=False, fontsize=9, bbox_to_anchor=(0.5, 0.05))
    plt.savefig(f"{FIG}/Fig6_rsfc_meta.png", dpi=300, bbox_inches="tight"); plt.close()

if __name__ == "__main__":
    import sys
    for w in (sys.argv[1:] or ["workflow", "main", "axis", "replication"]):
        try:
            {"workflow": fig_workflow, "main": fig_main, "axis": fig_axis, "replication": fig_replication, "rsfc": fig_rsfc}[w](); print(w, "ok")
        except Exception as e:
            print(w, "FAILED", repr(e))
