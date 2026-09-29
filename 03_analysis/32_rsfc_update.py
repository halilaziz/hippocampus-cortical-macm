"""Merge coordinates extracted during full-text assessment into the rsFC GingerALE foci files (MNI; TAL converted with
the Lancaster icbm2tal inverse, SPM/'other' variant), then write updated main and sensitivity sets."""
import os, glob, numpy as np, pandas as pd
L = "01_literature"; OUT = "04_results/rsfc_ale"
# Lancaster et al. (2007): inverse of icbm_other2tal gives Talairach -> MNI
TAL2ICBM_OTHER = np.linalg.inv(np.array([[0.9357, 0.0029, -0.0072, -1.0423], [-0.0065, 0.9396, -0.0726, -1.3940], [0.0103, 0.0752, 0.8967, 3.6475], [0, 0, 0, 1]]))
new = [pd.read_csv(f) for f in [f"{L}/second_screen/fulltext/rsfc_new_coords.csv", f"{L}/second_screen/fulltext/round2/rsfc_new_coords_round2.csv"] if os.path.exists(f)]
new = pd.concat(new, ignore_index=True) if new else pd.DataFrame()
new["inclusion"] = "main"
r3f = f"{L}/second_screen/fulltext/coords_round3/rsfc_new_coords_round3.csv"
if os.path.exists(r3f):
    r3 = pd.read_csv(r3f); r3.loc[r3.study_id == "Blessing2016", "N"] = 50   # conservative: conjunction of two samples of 50
    new = pd.concat([new, r3], ignore_index=True)
pdff = f"{L}/second_screen/fulltext/pdf_round/rsfc_new_coords_pdf.csv"
if os.path.exists(pdff):
    pdfc = pd.read_csv(pdff)
    # Vincent 2006 Table 2 lists parietal-restricted ROI centres, not whole-brain peaks: excluded, consistent with Vincent 2008
    new = pd.concat([new, pdfc[pdfc.study_id != "Vincent2006"]], ignore_index=True)
for f in sorted(glob.glob(f"{L}/second_screen/fulltext/pdf_round2/b*/rsfc_coords_b*.csv")):
    d = pd.read_csv(f)
    if len(d): new = pd.concat([new, d], ignore_index=True)
def block(study, seed, n, xyz):
    s = f"// {study}: {seed}\n// Subjects={int(n)}\n" + "".join(f"{x:.2f}\t{y:.2f}\t{z:.2f}\n" for x, y, z in xyz) + "\n"
    return s
base_ids = {l.split(":")[0][3:].strip() for l in open(f"{OUT}/rsfc_plus_sensitivity.txt") if l.startswith("// ") and ":" in l}
add = ""; add_sens = ""
new["exp_key"] = new["same_sample_id"].where(new.get("same_sample_id").notna(), new["study_id"]) if "same_sample_id" in new else new["study_id"]
new["exp_key"] = new["exp_key"].astype(str).str.replace(r"_(RW|Jena|NKIRS|S1|MGH31|Stanford36|NKI153|SALD262|Pitt7T15|Auburn7T23|Jena100|YA15|OA15)$", "", regex=True).where(~new["exp_key"].astype(str).str.startswith("Hermiller"), new["exp_key"])
for sid, g in new.groupby("exp_key"):
    if sid in base_ids: print(sid, "already in base set -> skipped (duplicate experiment)"); continue
    xyz = g[["x", "y", "z"]].astype(float).values
    if str(g.space.iloc[0]).upper().startswith("TAL"):
        xyz = (TAL2ICBM_OTHER @ np.c_[xyz, np.ones(len(xyz))].T).T[:, :3]
    ok = (np.abs(xyz[:, 0]) <= 80) & (xyz[:, 1] >= -120) & (xyz[:, 1] <= 90) & (xyz[:, 2] >= -60) & (xyz[:, 2] <= 90)
    b = block(sid, str(g.seed.iloc[0])[:60], g.N.iloc[0], xyz[ok])
    if (g.inclusion == "main").all(): add += b
    else: add_sens += b
    print(sid, g.inclusion.iloc[0], int(g.N.iloc[0]), int(ok.sum()), "foci", "TAL->MNI" if str(g.space.iloc[0]).upper().startswith("TAL") else "")
for base, name, extra in [("rsfc_main.txt", "rsfc_main_v4.txt", add), ("rsfc_plus_sensitivity.txt", "rsfc_plus_sensitivity_v4.txt", add + add_sens)]:
    txt = open(f"{OUT}/{base}").read().rstrip("\n") + "\n\n" + extra
    open(f"{OUT}/{name}", "w").write(txt)
    print(name, txt.count("Subjects="), "experiments", sum(1 for l in txt.splitlines() if l.strip() and not l.startswith("//")), "foci")
