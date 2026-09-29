"""Export seed-selected experiments as GingerALE (Sleuth-style) foci text files, MNI space."""
import sys, os
sys.path.insert(0, "03_analysis")
from lib import *
out = "04_results/gingerale/foci"; os.makedirs(out, exist_ok=True)
seeds = ["BIL_HIPP","L_HIPP","R_HIPP","BIL_aHIPP","BIL_pHIPP","L_aHIPP","L_pHIPP","R_aHIPP","R_pHIPP"]
def write(ds, ids, path):
    c = ds.coordinates[ds.coordinates["id"].isin(ids)]
    meta = ds.metadata.set_index("id")
    with open(path, "w") as f:
        f.write("// Reference=MNI\n")
        for i, g in c.groupby("id", sort=True):
            a = str(meta.loc[i, "authors"]).split(",")[0] if i in meta.index else ""
            y = meta.loc[i, "year"] if i in meta.index else ""
            f.write(f"// {a}, {y}: PMID {i}\n// Subjects={SAMPLE_SIZE}\n")
            for x, yy, z in g[["x","y","z"]].values:
                f.write(f"{x:.1f}\t{yy:.1f}\t{z:.1f}\n")
            f.write("\n")
    return c["id"].nunique(), len(c)
if __name__ == "__main__":
    for db in ["pooled", "neurosynth"]:
        ds = pooled_dataset() if db == "pooled" else load(db)
        sets = {s: set(seed_ids(ds, s)) for s in seeds}
        for s in seeds:
            print(db, s, write(ds, sets[s], f"{out}/{db}_{s}.txt"))
        # exclusive sets for contrasts
        for a, b in [("BIL_aHIPP","BIL_pHIPP"), ("L_HIPP","R_HIPP")]:
            A, B = sets[a] - sets[b], sets[b] - sets[a]
            print(db, a, "excl", write(ds, A, f"{out}/{db}_{a}_excl.txt"), b, "excl", write(ds, B, f"{out}/{db}_{b}_excl.txt"))
