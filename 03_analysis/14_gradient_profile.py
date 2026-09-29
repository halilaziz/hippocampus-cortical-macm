"""Summarise long-axis bins: fraction of significant (FDR q<.01, positive) cortical voxels per Yeo-7 network,
and mean association z per network; plot the anterior-posterior connectivity gradient."""
import sys, glob, os, json, numpy as np, nibabel as nib, pandas as pd, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from scipy.stats import spearmanr
sys.path.insert(0, "03_analysis")
from importlib import import_module
nt = import_module("08_networks_tables")
db = sys.argv[1] if len(sys.argv) > 1 else "neurosynth"
rows = []
for d in sorted(glob.glob(f"04_results/gradient/{db}/*/info.json")):
    info = json.load(open(d)); d = os.path.dirname(d)
    z = np.asarray(nib.load(f"{d}/z_desc-association.nii.gz").dataobj)
    p = np.asarray(nib.load(f"{d}/p_desc-association_level-voxel_corr-FDR_method-indep.nii.gz").dataobj)
    sig = (p < 0.01) & (z > 0)
    dP = np.asarray(nib.load(f"{d}/prob_desc-group1.nii.gz").dataobj) - np.asarray(nib.load(f"{d}/prob_desc-group2.nii.gz").dataobj)
    cort = nt.S > 0
    for n in nt.NETS:
        net = np.isin(nt.S, [k for k, v in nt.parcel_net.items() if v == n])
        rows.append(dict(**info, y_center=(info["y_lo"] + info["y_hi"]) / 2, network=nt.NETNAMES[n],
                         pct_of_sig_cortex=100 * (sig & net).sum() / max((sig & cort).sum(), 1),
                         pct_network_sig=100 * (sig & net).sum() / net.sum(), mean_z=float(z[net].mean()), delta_p=100*float(dP[net].mean())))
df = pd.DataFrame(rows); df.to_csv(f"04_results/tables/gradient_{db}.csv", index=False)
stats = []
for n, g in df.groupby("network"):
    rho, pv = spearmanr(g.y_center, g.mean_z); rho2, pv2 = spearmanr(g.y_center, g.delta_p)
    stats.append(dict(network=n, rho_y_vs_meanz=round(rho, 3), p=round(pv, 4), rho_y_vs_deltaP=round(rho2, 3), p_deltaP=round(pv2, 4)))
st = pd.DataFrame(stats); st.to_csv(f"04_results/tables/gradient_{db}_trend.csv", index=False)
print(df.pivot_table(index="network", columns="y_center", values="pct_network_sig").round(1).to_string())
print(df.pivot_table(index="network", columns="y_center", values="mean_z").round(2).to_string()); print(df.pivot_table(index="network", columns="y_center", values="delta_p").round(2).to_string()); print(st.to_string(index=False))
print(df.drop_duplicates("bin")[["bin", "n_vox", "n_seed", "mean_foci"]].to_string(index=False))
cols = {"Visual": "#781C86", "Somatomotor": "#4682B4", "Dorsal attention": "#00760E", "Salience/ventral attention": "#C43AFA",
        "Limbic": "#b8a44a", "Frontoparietal control": "#E69422", "Default mode": "#CD3E4E"}
fig, ax = plt.subplots(1, 2, figsize=(12, 4.2))
for n, g in df.groupby("network"):
    g = g.sort_values("y_center")
    ax[0].plot(g.y_center, g.pct_network_sig, "-o", color=cols[n], label=n, lw=2)
    ax[1].plot(g.y_center, g.delta_p, "-o", color=cols[n], label=n, lw=2)
for a, t in zip(ax, ["% of network voxels with specific co-activation (FDR q<.01)", "Effect size: ΔP(activation) seed vs matched reference (%)"]):
    a.set_xlabel("Hippocampal seed position, MNI y (mm)  [posterior ← → anterior]"); a.set_title(t, fontsize=10); a.axvline(-21, ls="--", c="grey", lw=1)
    a.spines[["top", "right"]].set_visible(False)
ax[1].axhline(0, c="k", lw=0.5); ax[1].legend(fontsize=8, frameon=False, loc="best")
plt.tight_layout(); plt.savefig(f"04_results/figures/gradient_{db}.png", dpi=200); print("saved")
