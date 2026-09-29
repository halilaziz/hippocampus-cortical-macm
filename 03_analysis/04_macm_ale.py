"""ALE-MACM (GingerALE-equivalent: Turkeltaub 2012 non-additive MA, Eickhoff 2012 analytic null,
cluster-level FWE p<.05 with voxel-level cluster-forming p<.001, Monte Carlo permutations)."""
import sys, os, time, json
sys.path.insert(0, "03_analysis")
from lib import *
from nimare.meta.cbma.ale import ALE
from nimare.meta.kernel import ALEKernel
from nimare.correct import FWECorrector

db, seed, n_iters = sys.argv[1], sys.argv[2], int(sys.argv[3])
cores = int(os.environ.get("NCORES", 12))
ds = {"neurosynth": load, "neuroquery": load}.get(db, None)
ds = ds(db) if ds else (pooled_dataset() if db == "pooled" else nq_independent())
ids = seed_ids(ds, seed)
sub = ds.slice(ids)
out = f"04_results/macm/{db}/{seed}"; os.makedirs(out, exist_ok=True)
t = time.time()
est = ALE(kernel_transformer=ALEKernel(sample_size=SAMPLE_SIZE), null_method="approximate")
res = est.fit(sub)
cres = FWECorrector(method="montecarlo", voxel_thresh=0.001, n_iters=n_iters, n_cores=cores, vfwe_only=False).transform(res)
cres.save_maps(output_dir=out)
info = dict(db=db, seed=seed, n_experiments=len(ids), n_foci=int(len(sub.coordinates)), n_iters=n_iters,
            sample_size=SAMPLE_SIZE, seconds=round(time.time() - t, 1))
json.dump(info, open(f"{out}/info.json", "w"), indent=1)
sub.coordinates.to_csv(f"{out}/foci.csv", index=False)
print(info)
