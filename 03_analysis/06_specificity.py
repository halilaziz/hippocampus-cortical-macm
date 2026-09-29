"""Specific co-activation (MKDA chi-square; Wager et al., 2007): seed-activating studies vs. a foci-count-MATCHED
reference set of non-seed studies from the same database. Matching removes the confound that studies reporting
hippocampal foci report more foci overall (mean 56.5 vs 32.1 in Neurosynth), which otherwise inflates
association everywhere. Positive association = activation more likely when the hippocampus is co-activated."""
import sys, os, json, time, numpy as np
sys.path.insert(0, "03_analysis")
from lib import *
from nimare.meta.cbma.mkda import MKDAChi2
from nimare.meta.kernel import MKDAKernel
from nimare.correct import FDRCorrector

def matched_reference(ds, ids, rng):
    n = ds.coordinates.groupby("id").size()
    seed = n[n.index.isin(ids)]
    pool = n[~n.index.isin(ids)].sample(frac=1, random_state=rng).sort_values(kind="stable")
    pool_ids, pool_n = list(pool.index), pool.values.astype(float)
    avail = np.ones(len(pool_ids), bool); chosen = []
    for sid, k in seed.sample(frac=1, random_state=rng).items():
        cand = np.where(avail)[0]
        j = cand[np.argmin(np.abs(pool_n[cand] - k))]
        avail[j] = False; chosen.append(pool_ids[j])
    return chosen, float(seed.mean()), float(n[chosen].mean())

if __name__ == "__main__":
    db = sys.argv[1]; seeds = sys.argv[2].split(",")
    ds = pooled_dataset() if db == "pooled" else (load(db) if db != "nq_indep" else nq_independent())
    for seed in seeds:
        out = f"04_results/specificity/{db}/{seed}"; os.makedirs(out, exist_ok=True)
        if os.path.exists(f"{out}/info.json"): continue
        t = time.time()
        ids = seed_ids(ds, seed)
        ref, m1, m2 = matched_reference(ds, ids, 42)
        est = MKDAChi2(kernel_transformer=MKDAKernel(r=10))
        res = est.fit(ds.slice(ids), ds.slice(ref))
        FDRCorrector(method="indep", alpha=0.01).transform(res).save_maps(output_dir=out)
        json.dump(dict(db=db, seed=seed, n_seed=len(ids), n_ref=len(ref), mean_foci_seed=round(m1,1), mean_foci_ref=round(m2,1),
                       seconds=round(time.time()-t)), open(f"{out}/info.json", "w"))
        open(f"{out}/reference_ids.txt", "w").write("\n".join(ref))
        print(open(f"{out}/info.json").read())
