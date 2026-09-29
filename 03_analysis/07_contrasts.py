"""ALE subtraction contrasts (Eickhoff et al., 2011 permutation approach, as in GingerALE):
anterior- vs posterior-exclusive and left- vs right-exclusive experiment sets.
Studies activating both seeds are excluded so groups are independent."""
import sys, os, json, time
sys.path.insert(0, "03_analysis")
from lib import *
from nimare.meta.cbma.ale import ALESubtraction
from nimare.meta.kernel import ALEKernel
db, a, b, n_iters = sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4])
ds = pooled_dataset() if db == "pooled" else load(db)
A, B = set(seed_ids(ds, a)), set(seed_ids(ds, b))
Ae, Be = sorted(A - B), sorted(B - A)
out = f"04_results/contrasts/{db}/{a}_vs_{b}"; os.makedirs(out, exist_ok=True)
t = time.time()
est = ALESubtraction(kernel_transformer=ALEKernel(sample_size=SAMPLE_SIZE), n_iters=n_iters, n_cores=int(os.environ.get("NCORES", 12)))
res = est.fit(ds.slice(Ae), ds.slice(Be))
res.save_maps(output_dir=out)
json.dump(dict(db=db, a=a, b=b, n_a=len(Ae), n_b=len(Be), n_shared_excluded=len(A & B), n_iters=n_iters, seconds=round(time.time()-t)), open(f"{out}/info.json", "w"), indent=1)
print(open(f"{out}/info.json").read())
