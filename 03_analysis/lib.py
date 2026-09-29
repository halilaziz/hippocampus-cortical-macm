"""Shared helpers for hippocampal MACM analyses."""
import os, numpy as np, pandas as pd
from nimare.dataset import Dataset

SEEDDIR = "04_results/seeds"
SAMPLE_SIZE = 20  # fixed N for the Eickhoff (2009) kernel; Neurosynth/NeuroQuery lack per-study N

def load(name):
    return Dataset.load(f"02_data/{name}_dataset.pkl.gz")

def pooled_dataset():
    """Union of Neurosynth + NeuroQuery, deduplicated by PMID (Neurosynth coordinates kept for overlap)."""
    path = "02_data/pooled_dataset.pkl.gz"
    if os.path.exists(path):
        return Dataset.load(path)
    ns, nq = load("neurosynth"), load("neuroquery")
    ns_pmids = set(ns.coordinates.study_id.astype(str))
    nq_sub = nq.slice([i for i in nq.ids if i.split("-")[0] not in ns_pmids])
    ns = ns.copy()
    for d, src in [(ns, "neurosynth"), (nq_sub, "neuroquery")]:
        d.annotations = d.annotations[["id", "study_id", "contrast_id"]]
        d.metadata = d.metadata.assign(source=src)
    ds = ns.merge(nq_sub)
    ds.save(path)
    return ds

def nq_independent():
    ns, nq = load("neurosynth"), load("neuroquery")
    ns_pmids = set(ns.coordinates.study_id.astype(str))
    return nq.slice([i for i in nq.ids if i.split("-")[0] not in ns_pmids])

def seed_ids(ds, seed):
    return list(ds.get_studies_by_mask(f"{SEEDDIR}/{seed}.nii.gz"))
