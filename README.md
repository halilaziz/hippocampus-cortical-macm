# Specific cortical partners and long-axis gradient of the hippocampus

Code, derived data and statistical maps for:

> Velioglu HA. *Specific cortical partners and long-axis gradient of the hippocampus.* (manuscript submitted).

The project combines a PRISMA-guided systematic review of hippocampal functional-connectivity studies with three meta-analytic approaches:

1. **Meta-analytic connectivity modeling (MACM)** of hippocampal seeds in Neurosynth v7 and NeuroQuery (2,598 experiments in the pooled bilateral analysis), using activation likelihood estimation (ALE; the GingerALE algorithm as implemented in NiMARE) and voxel-wise equivalence tests against native GingerALE 3.0.2.
2. **Specific co-activation** (MKDA χ²) against reference experiments matched for the number of reported foci, and a **voxel-wise meta-analytic connectivity gradient** (diffusion-map embedding) tested against a geometry-preserving permutation null.
3. **Coordinate-based meta-analysis (GingerALE)** of published hippocampal-seed resting-state studies in healthy adults (20 experiments, 785 participants).

## Repository layout

| Path | Content |
|---|---|
| `03_analysis/` | Analysis pipeline (numbered in execution order) |
| `04_results/seeds/` | Hippocampal seed masks (MNI152, 2 mm) |
| `04_results/gingerale/foci/` | Seed-selected experiments in GingerALE (Sleuth) text format, MNI |
| `04_results/rsfc_ale/` | Resting-state coordinate sets (`rsfc_main_v4.txt`, `rsfc_plus_sensitivity_v4.txt`) and GingerALE/NiMARE outputs |
| `04_results/macm/`, `specificity/`, `contrasts/`, `meta_gradient/` | Unthresholded and FWE/FDR-corrected statistical maps |
| `04_results/tables/` | Cluster tables, network profiles, thresholded maps, replication and gradient statistics |
| `04_results/figures/pub/` | Figures |
| `01_literature/` | Search strings and counts, screening protocol, first/second-screener decisions, author adjudication, full-text decisions, key-study table |
| `01_literature/rsfc_coords/`, `01_literature/second_screen/fulltext/**/rsfc_*coords*.csv` | Extracted resting-state peak coordinates with source table for every focus |

Full texts, PDFs and abstracts are **not** redistributed (copyright). The GingerALE whole-brain mask is extracted automatically from a local GingerALE installation by the shell scripts.

## Reproducing the analyses

```bash
python -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
python 03_analysis/01_fetch_databases.py      # downloads Neurosynth and NeuroQuery (NiMARE)
python 03_analysis/02_seeds.py
bash   03_analysis/run_macm_batch.sh          # ALE-MACM, 1,000 Monte Carlo permutations per analysis
python 03_analysis/06_specificity.py pooled BIL_HIPP,L_HIPP,R_HIPP,BIL_aHIPP,BIL_pHIPP
python 03_analysis/07_contrasts.py pooled BIL_aHIPP BIL_pHIPP 10000
python 03_analysis/26_meta_gradient.py pooled 150
python 03_analysis/28_gradient_content_null.py 1000
bash   03_analysis/run_v4.sh                  # resting-state GingerALE meta-analysis + jackknife
```

Native GingerALE runs require GingerALE 3.0.2 (BrainMap) and Java; set `GINGERALE_JAR` if it is not at the default macOS location. Scripts expect to be run from the repository root. The literature scripts (`30_`, `31_`) document the screening workflow; they depend on raw records that are not redistributed.

## Software

Python 3.12, NiMARE 0.21, nilearn 0.13 and other packages listed in `requirements.txt`; GingerALE 3.0.2.

## License

Code: MIT License (see `LICENSE`). Derived data, maps and tables: CC BY 4.0. Neurosynth and NeuroQuery data are subject to their own licenses.

## Citation

Please cite the article (see `CITATION.cff`) and the underlying resources: Neurosynth (Yarkoni et al. 2011), NeuroQuery (Dockès et al. 2020), NiMARE (Salo et al. 2023) and GingerALE (Eickhoff et al. 2009, 2012; Turkeltaub et al. 2012).
