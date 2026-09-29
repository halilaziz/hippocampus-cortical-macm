"""Functional characterisation: term enrichment (Neurosynth v7 abstract tf-idf > 0.001) of seed-activating
studies vs. the rest of Neurosynth (chi-square, FDR). Also anterior-exclusive vs posterior-exclusive comparison."""
import sys, os, numpy as np, pandas as pd
sys.path.insert(0, "03_analysis")
from lib import *
from scipy.stats import chi2_contingency
from statsmodels.stats.multitest import multipletests
ds = load("neurosynth")
ann = ds.annotations.set_index("id")
terms = [c for c in ann.columns if c.startswith("terms_abstract_tfidf__")]
X = (ann[terms] > 0.001)
X.columns = [c.split("__")[1] for c in terms]
X = X.loc[:, X.sum() >= 50]   # terms used in >=50 studies
out = "04_results/decoding"; os.makedirs(out, exist_ok=True)
def enrich(A, B, label):
    a = X.loc[A]; b = X.loc[B]
    rows = []
    for t in X.columns:
        k1, n1, k2, n2 = a[t].sum(), len(a), b[t].sum(), len(b)
        if k1 + k2 < 10: continue
        tab = np.array([[k1, n1-k1], [k2, n2-k2]])
        chi, p, _, _ = chi2_contingency(tab, correction=True)
        rows.append((t, k1, k1/n1, k2/n2, (k1/n1)/(max(k2,0.5)/n2), chi, p))
    df = pd.DataFrame(rows, columns=["term","n_seed","p_term_seed","p_term_ref","ratio","chi2","p"])
    df["q_fdr"] = multipletests(df.p, method="fdr_bh")[1]
    df = df.sort_values("chi2", ascending=False)
    df.to_csv(f"{out}/{label}.csv", index=False)
    top = df[(df.q_fdr < 0.001) & (df.ratio > 1)].head(25)
    print(label, "enriched terms q<.001:", ((df.q_fdr<0.001)&(df.ratio>1)).sum()); print(top[["term","n_seed","ratio","q_fdr"]].to_string(index=False))
sets = {s: set(seed_ids(ds, s)) for s in ["BIL_HIPP","BIL_aHIPP","BIL_pHIPP","L_HIPP","R_HIPP"]}
allids = set(ds.ids)
enrich(sorted(sets["BIL_HIPP"]), sorted(allids - sets["BIL_HIPP"]), "BIL_HIPP_vs_rest")
A = sets["BIL_aHIPP"] - sets["BIL_pHIPP"]; P = sets["BIL_pHIPP"] - sets["BIL_aHIPP"]
enrich(sorted(A), sorted(P), "aHIPPexcl_vs_pHIPPexcl")
enrich(sorted(P), sorted(A), "pHIPPexcl_vs_aHIPPexcl")
L = sets["L_HIPP"] - sets["R_HIPP"]; R = sets["R_HIPP"] - sets["L_HIPP"]
enrich(sorted(L), sorted(R), "Lexcl_vs_Rexcl")
enrich(sorted(R), sorted(L), "Rexcl_vs_Lexcl")
