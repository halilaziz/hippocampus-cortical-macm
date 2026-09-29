"""Content-based null tests for the meta-analytic gradient.

Geometry alone makes the leading eigenvector of any smooth field on an elongated structure follow its long axis, so the
long-axis alignment of G1 is not by itself evidence. Here we test the *content* of the gradient against a null in which
cortical profiles are shuffled across hippocampal experiments (same spatial filters, same geometry):
 (1) network-level association between long-axis position and specific co-activation (per Yeo-7 network);
 (2) replicability of the G1 cortical projection pattern across independent databases, compared with the similarity
     expected between null projections.
"""
import sys, json, numpy as np, pandas as pd
from scipy.stats import pearsonr
sys.path.insert(0, "03_analysis")
from importlib import import_module
g26 = open("03_analysis/26_meta_gradient.py").read()
ns_ = {}; exec(g26[g26.index("def diffusion_map"):g26.index("emb, lam = diffusion_map(prof)")], {"np": np}, ns_); diffusion_map = ns_["diffusion_map"]
nt = import_module("08_networks_tables")
lab = [l.decode() if isinstance(l, bytes) else l for l in nt.sch.labels]
if lab[0].lower().startswith("background"): lab = lab[1:]
net = np.array([nt.NETNAMES[l.split("_")[2]] for l in lab]); NETS = sorted(set(net))

def load(name):
    Z = np.load(f"04_results/meta_gradient/{name}/gradient_data.npz", allow_pickle=True)
    studies = list(Z["studies"]); sidx = {s: i for i, s in enumerate(studies)}
    hip_sidx = np.array([sidx[s] for s in Z["hip_ids"]])
    return Z, hip_sidx
def longaxis(Hxyz):
    LA = np.zeros(len(Hxyz))
    for hm in [Hxyz[:, 0] < 0, Hxyz[:, 0] >= 0]:
        X = Hxyz[hm] - Hxyz[hm].mean(0); u = np.linalg.svd(X, full_matrices=False)[2][0]; u = u if u[1] > 0 else -u; LA[hm] = X @ u
    return LA
def profiles(Z, ps):
    A, Ebin, edges, nf, order = Z["A"], Z["Ebin"], Z["edges"], Z["nfoci"], Z["order"]
    binof = lambda n: np.clip(np.searchsorted(edges, n, side="right") - 1, 0, len(edges) - 2)
    P = np.zeros((len(order), 400), np.float32)
    for v in range(len(order)):
        si = ps[order[v]]; P[v] = A[si].mean(0) - Ebin[binof(nf[si])].mean(0)
    return P
def stats(P, LA):
    e, l = diffusion_map(P)
    g = e[:, 0] * np.sign(np.corrcoef(e[:, 0], LA)[0, 1])
    proj = np.array([np.corrcoef(g, P[:, p])[0, 1] if P[:, p].std() > 0 else 0 for p in range(400)])
    netr = {n: np.corrcoef(LA, P[:, net == n].mean(1))[0, 1] for n in NETS}
    return proj, netr

NPERM = int(sys.argv[1]) if len(sys.argv) > 1 else 200
Zp, hp = load("pooled_K150"); LA = longaxis(Zp["Hxyz"])
obs_proj, obs_netr = stats(Zp["prof"], LA)
rng = np.random.default_rng(99); null_proj, null_netr = [], []
for it in range(NPERM):
    pr, nr = stats(profiles(Zp, hp[rng.permutation(len(hp))]), LA); null_proj.append(pr); null_netr.append(nr)
null_proj = np.array(null_proj); NN = pd.DataFrame(null_netr)
rows = []
for n in NETS:
    o = obs_netr[n]; nd = NN[n].values
    rows.append(dict(network=n, r_longaxis_vs_specific_coact=round(o, 3), null_2p5=round(np.quantile(nd, .025), 3), null_97p5=round(np.quantile(nd, .975), 3),
                     p_two_sided=round((1 + (np.abs(nd) >= abs(o)).sum()) / (NPERM + 1), 4)))
net_df = pd.DataFrame(rows); net_df.to_csv("04_results/tables/meta_gradient/network_longaxis_permutation.csv", index=False); print(net_df.to_string(index=False))
# replication of projection pattern across independent databases vs null similarity
Zn, _ = load("neurosynth_K150"); Zq, _ = load("nq_indep_K100")
pn, _ = stats(Zn["prof"], longaxis(Zn["Hxyz"])); pq, _ = stats(Zq["prof"], longaxis(Zq["Hxyz"]))
r_obs = pearsonr(pn, pq)[0]
iu = np.triu_indices(NPERM, 1); null_pair = np.corrcoef(null_proj)[iu]
res = dict(n_perm=NPERM, r_projection_neurosynth_vs_nqindep=float(r_obs), r_projection_pooled_vs_nqindep=float(pearsonr(obs_proj, pq)[0]),
           null_pairwise_projection_r_mean=float(null_pair.mean()), null_pairwise_projection_r_95=float(np.quantile(null_pair, .95)),
           null_pairwise_projection_r_99=float(np.quantile(null_pair, .99)),
           p_replication=float((1 + (null_pair >= r_obs).sum()) / (len(null_pair) + 1)),
           r_obs_projection_vs_null_mean_projection=float(pearsonr(obs_proj, null_proj.mean(0))[0]))
json.dump(res, open("04_results/tables/meta_gradient/projection_replication_null.json", "w"), indent=1); print(json.dumps(res, indent=1))
