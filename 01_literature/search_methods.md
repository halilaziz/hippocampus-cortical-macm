# Search methods: HIPPCON literature search

*Draft text for the Methods section. Every number below is generated from `prisma_counts.json` and `search_log.csv`.*

## Information sources and search strategy

We searched PubMed/MEDLINE on 24 September 2026 through the NCBI E-utilities API (`esearch`, `efetch`, `esummary`). All requests were capped at 3 per second, used `retmax=10000`, and were limited to publication dates up to 2026/09/24.

Ten search strings were run, grouped into four blocks:
- **Block A: coordinate-based and meta-analytic connectivity.** Covers MACM, ALE, co-activation, BrainMap, Neurosynth and connectivity-based parcellation of the hippocampus and medial temporal lobe (A1-A3).
- **Block B: resting-state or intrinsic FC with hippocampal or hippocampal-subregion seeds.** Subregions include anterior/posterior, head/body/tail, and CA/subiculum/dentate gyrus. Block B includes a MeSH-based query and a sensitivity query for alternative wording such as "spontaneous activity", "effective connectivity", "functional coupling" and "cortico-hippocampal" (B1-B4).
- **Block C: task-based FC.** Covers PPI/gPPI, beta-series, background connectivity and task-state coupling (C1-C2).
- **Block D: the pre-specified title/abstract query** (D1).

To capture preprints, we ran one extra query through the Europe PMC REST API, restricted to preprint sources (`SRC:PPR`) (E1).

Search terms were combined with the Boolean operators shown in Table S1. Titles and abstracts were searched with the `[tiab]` field tag and wildcard truncation (`hippocamp*`). The `humans[mh]` filter was applied where stated.

After screening showed that several landmark studies were missed by the connectivity terminology, we added the sensitivity query B4. Known landmark studies that were still not retrieved were then located with targeted author/title PubMed look-ups and reference knowledge ("identification via other methods" in PRISMA 2020; n = 19). Each of these records was verified against its PubMed entry.

## Records and de-duplication

The PubMed queries returned 6107 hits, corresponding to 3341 unique PMIDs. Europe PMC returned 278 preprint records.

Duplicates were removed in three steps:
1. By PMID.
2. By DOI, then by normalised title. When a preprint and its journal version were both indexed in PubMed (n = 7 pairs), the peer-reviewed version was kept.
3. Europe PMC preprints already represented in PubMed by DOI or title were removed (n = 101).

In total, 2874 duplicates were removed, leaving **3511 unique records** for screening. For each record we retrieved PMID, title, authors, year, journal, DOI and abstract from PubMed XML (`records_all.csv`).

## Eligibility criteria

We included studies that met all of the following:
- **(I1)** Human participants.
- **(I2)** Healthy adults (≥18 years; young, older or lifespan samples), or a healthy control group whose hippocampal connectivity was analysed and reported separately.
- **(I3)** Hippocampus used as a seed or ROI whose connectivity was analysed. This could be the whole hippocampus, left/right, a long-axis segment, head/body/tail, subfields, or a hippocampal gradient or parcel.
- **(I4)** Connectivity reported with the whole brain or broadly with the cerebral cortex. This included seed-to-voxel analyses, seed-to-many-parcel analyses, network assignment of hippocampal components, and connectivity-based parcellations or gradients defined from hippocampus-to-brain connectivity.
- **(I5)** fMRI (resting-state, naturalistic or task, BOLD), or a coordinate-based meta-analysis of functional neuroimaging (MACM, ALE co-activation, BrainMap, Neurosynth).

We excluded the following (codes are used in `screening.csv`):
- E1: animal studies.
- E2: non-fMRI modalities.
- E5: studies whose focus was a clinical, at-risk, genetic-risk or exposure population (patient vs control).
- E6: paediatric- or adolescent-only samples.
- E7: studies in which the hippocampus was not a connectivity seed.
- E8: reviews, commentaries and protocols.
- E9: connectivity restricted to a few ROIs, not whole-brain or cortical.
- E10: studies whose primary outcome was an intervention or state manipulation, with no baseline healthy map.
- E12: preprints.
- E13: duplicate reports.

## Screening procedure

Screening used two stages. **Be transparent about this in the paper.**

**Stage 1: keyword-rule exclusion.** Deterministic, conservative regular-expression rules (`scripts/04_screen_rules.py`, `05_screen_stage2.py`) excluded only records that clearly met an exclusion criterion from the title and abstract:
- animal-only studies;
- no fMRI or connectivity content;
- errata, editorials and letters;
- no mention of the hippocampus;
- reviews without meta-analytic connectivity analyses;
- a clinical, at-risk or exposure population named in the title;
- paediatric-only samples named in the title.

This stage excluded 2026 records. A first version of the rules also auto-excluded records with no hippocampal "seed" wording. A blinded audit of a stratified random sample of rule-based exclusions (n = 130) found a false-exclusion rate of about 1.7%. The audit also exposed a word-boundary bug: "lobes" matched the term "obes". We therefore fixed the regular expressions and sent every record without seed wording to full abstract reading. In the audited strata that remain rule-based (animal, non-fMRI, no-hippocampus, review, clinical title and paediatric title; n = 90), the only disagreement was the one caused by the fixed bug.

**Stage 2: full abstract reading.** The remaining 1485 records were screened by reading the title and abstract against a written protocol (`raw/manual/SCREENING_PROTOCOL.txt`). The screener recorded:
- the decision (include / maybe / exclude) and the reason code;
- the sample, N, seed definition and cortical targets named in the abstract;
- whether a healthy control group with a hippocampal seed was present.

**Disclosure required.** This abstract reading was done by an AI assistant (Claude, Anthropic) working through the protocol in parallel batches, not by human raters. Before publication, a human reviewer should double-screen the included and "maybe" records and a random sample of exclusions, and report inter-rater agreement (for example Cohen's κ).

A reviewer pass then adjudicated eligibility at abstract level. This removed preprints and duplicate reports and applied documented overrides, including overrides arising from verification of key studies. Every override is recorded, with its reason, in `screening.csv`.

## Results of the search (PRISMA 2020 counts)

- **Identified:** 6385 records (PubMed 6107 hits across 10 queries, plus Europe PMC 278), plus 19 through other methods.
- **Duplicates removed:** 2874.
- **Screened:** 3511.
- **Excluded at title/abstract stage:** 2026 by keyword rules and 1192 after abstract reading.
- **Assessed for eligibility:** 293 reports. Of these, 25 were excluded (E5: 3, E7: 1, E13: 4, E12: 17).
- **Pending full-text confirmation:** 149 reports, plus 3 from other methods.
- **Included:** 119 reports from the database search (meta 3, mixed 8, rsfc 65, task 43), plus 2 through other methods, for a total of **121**.

The pending reports need full-text retrieval before the PRISMA flow is final. Up to 273 reports could be eligible.

## Data extraction (planned)

For each included report we will extract:
- sample characteristics;
- acquisition details (field strength, resolution);
- seed definition (atlas, long-axis landmark, subfield segmentation);
- connectivity method (seed-based, ICA, gradient, PPI/gPPI, beta-series, MACM);
- analysis space;
- statistical thresholding;
- peak coordinates of cortical connectivity effects.

These data will feed coordinate-based meta-analyses (ALE/MKDA in NiMARE, or GingerALE). MACM will be replicated in both BrainMap and Neurosynth/NeuroQuery.

## Table S1: search strings (run 2026-09-24)

| ID | Block | Query | Hits |
|---|---|---|---|
| A1 | a) coordinate-based meta-analysis / MACM / ALE / co-activation / BrainMap / Neurosynth / CBP of hippocampus-MTL | `(hippocamp*[tiab] OR "medial temporal lobe"[tiab] OR parahippocamp*[tiab]) AND ("meta-analytic connectivity"[tiab] OR MACM[tiab] OR "activation likelihood estimation"[tiab] OR coactivation[tiab] OR co-activation[tiab] OR BrainMap[tiab] OR Neurosynth[tiab] OR "connectivity-based parcellation"[tiab] OR "coordinate-based meta-analysis"[tiab]) AND ("1900/01/01"[dp] : "2026/09/24"[dp])` | 520 |
| A2 | a) hippocampal parcellation / gradients / long-axis organisation from connectivity | `hippocamp*[tiab] AND (parcellat*[tiab] OR gradient*[tiab] OR connectopic*[tiab] OR "long axis"[tiab] OR "longitudinal axis"[tiab]) AND ("functional connectivity"[tiab] OR coactivation[tiab] OR co-activation[tiab] OR "meta-analytic"[tiab]) AND (fMRI[tiab] OR "resting state"[tiab] OR "resting-state"[tiab] OR "magnetic resonance"[tiab] OR "meta-analytic"[tiab]) AND ("1900/01/01"[dp] : "2026/09/24"[dp])` | 119 |
| A3 | a) meta-analyses of hippocampal/MTL connectivity (publication type or tiab) | `(hippocamp*[tiab] OR "medial temporal lobe"[tiab]) AND ("meta-analysis"[pt] OR "meta-analysis"[tiab] OR "meta-analytic"[tiab] OR "meta-analyses"[tiab]) AND ("functional connectivity"[tiab] OR coactivation[tiab] OR co-activation[tiab] OR "connectivity"[tiab]) AND (fMRI[tiab] OR neuroimaging[tiab] OR "functional magnetic resonance"[tiab] OR "resting state"[tiab] OR "resting-state"[tiab]) AND ("1900/01/01"[dp] : "2026/09/24"[dp])` | 81 |
| B1 | b) resting-state FC of hippocampus | `hippocamp*[tiab] AND ("resting state"[tiab] OR "resting-state"[tiab] OR "intrinsic connectivity"[tiab] OR "intrinsic functional connectivity"[tiab]) AND "functional connectivity"[tiab] AND humans[mh] AND ("1900/01/01"[dp] : "2026/09/24"[dp])` | 1303 |
| B2 | b) FC of hippocampal subregions (long axis / head-body-tail / subfields) | `hippocamp*[tiab] AND ("anterior hippocampus"[tiab] OR "posterior hippocampus"[tiab] OR "hippocampal head"[tiab] OR "hippocampal body"[tiab] OR "hippocampal tail"[tiab] OR "long axis"[tiab] OR "longitudinal axis"[tiab] OR subfield*[tiab] OR subiculum[tiab] OR "dentate gyrus"[tiab] OR CA1[tiab] OR CA3[tiab] OR "cornu ammonis"[tiab] OR subregion*[tiab]) AND ("functional connectivity"[tiab] OR "resting state"[tiab] OR "resting-state"[tiab] OR "intrinsic connectivity"[tiab]) AND (fMRI[tiab] OR "functional MRI"[tiab] OR "functional magnetic resonance"[tiab] OR "BOLD"[tiab] OR "resting state"[tiab] OR "resting-state"[tiab]) AND humans[mh] AND ("1900/01/01"[dp] : "2026/09/24"[dp])` | 359 |
| B3 | b) MeSH-based: hippocampus + rest + MRI | `Hippocampus[mh] AND Rest[mh] AND "Magnetic Resonance Imaging"[mh] AND humans[mh] AND ("1900/01/01"[dp] : "2026/09/24"[dp])` | 120 |
| B4 | b) sensitivity supplement: alternative connectivity terminology (spontaneous activity, effective connectivity, coupling, cortico-hippocampal) | `hippocamp*[tiab] AND ("effective connectivity"[tiab] OR "spontaneous activity"[tiab] OR "spontaneous BOLD"[tiab] OR "spontaneous fluctuations"[tiab] OR "low-frequency fluctuations"[tiab] OR "intrinsic activity"[tiab] OR "functional coupling"[tiab] OR "functional interaction"[tiab] OR "functional interactions"[tiab] OR "connectivity MRI"[tiab] OR fcMRI[tiab] OR "functional networks"[tiab] OR "functional network"[tiab] OR "network connectivity"[tiab] OR "cortico-hippocampal"[tiab] OR "hippocampal-cortical"[tiab] OR "hippocampal-neocortical"[tiab] OR "hippocampo-cortical"[tiab] OR "hippocampal connectivity"[tiab]) AND ("Magnetic Resonance Imaging"[mh] OR fMRI[tiab] OR "functional MRI"[tiab] OR "functional magnetic resonance"[tiab] OR BOLD[tiab]) AND humans[mh] AND ("1900/01/01"[dp] : "2026/09/24"[dp])` | 1130 |
| C1 | c) task-based FC: PPI / beta-series / background connectivity | `hippocamp*[tiab] AND ("psychophysiological interaction"[tiab] OR PPI[tiab] OR gPPI[tiab] OR "beta series"[tiab] OR "beta-series"[tiab] OR "background connectivity"[tiab] OR "task-based functional connectivity"[tiab] OR "task-related functional connectivity"[tiab] OR "task-dependent functional connectivity"[tiab] OR "task-evoked functional connectivity"[tiab]) AND (fMRI[tiab] OR "functional magnetic resonance"[tiab] OR "functional MRI"[tiab]) AND ("1900/01/01"[dp] : "2026/09/24"[dp])` | 125 |
| C2 | c) task-based hippocampal coupling during memory/navigation tasks | `hippocamp*[tiab] AND ("functional coupling"[tiab] OR "functional connectivity"[tiab]) AND (encoding[tiab] OR retrieval[tiab] OR "memory task"[tiab] OR "navigation"[tiab]) AND (fMRI[tiab] OR "functional magnetic resonance"[tiab]) AND humans[mh] NOT ("resting state"[tiab] OR "resting-state"[tiab]) AND ("1900/01/01"[dp] : "2026/09/24"[dp])` | 181 |
| D1 | precise tiab query (as specified in protocol) | `(hippocamp*[tiab]) AND ("functional connectivity"[tiab] OR coactivation OR co-activation OR "meta-analytic connectivity") AND (fMRI OR "resting state" OR "resting-state") AND humans[mh] AND ("1900/01/01"[dp] : "2026/09/24"[dp])` | 2169 |
| E1 | supplementary: preprints (bioRxiv/medRxiv/Research Square etc.) via Europe PMC | `(TITLE_ABS:hippocamp*) AND (TITLE_ABS:"functional connectivity" OR TITLE_ABS:coactivation OR TITLE_ABS:"co-activation" OR TITLE_ABS:"meta-analytic connectivity") AND (TITLE_ABS:fMRI OR TITLE_ABS:"resting state" OR TITLE_ABS:"resting-state") AND SRC:PPR AND FIRST_PDATE:[1900-01-01 TO 2026-09-24]` | 278 |
