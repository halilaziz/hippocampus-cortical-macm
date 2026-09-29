"""Updated PRISMA 2020 counts after dual screening and full-text assessment, plus a PRISMA flow diagram."""
import json, os, pandas as pd, numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
L = "01_literature"; SS = f"{L}/second_screen"; FT = f"{SS}/fulltext"
P = json.load(open(f"{L}/prisma_counts.json")); A = json.load(open(f"{SS}/agreement_summary.json"))
def dec(path, col_candidates=("decision", "ft_decision", "final_decision")):
    if not os.path.exists(path): return None
    d = pd.read_csv(path); c = [x for x in d.columns if x.lower() in col_candidates][0]
    d["_dec"] = d[c].astype(str).str.lower(); return d
ft1 = dec(f"{FT}/fulltext_decisions.csv"); ft2 = dec(f"{FT}/round2/fulltext_decisions_round2.csv"); inc = dec(f"{FT}/included_check/included_fulltext_decisions.csv")
# institutional-access PDFs (lab machine) override "unresolved" decisions
import glob as _g
_prs = [dec(f) for f in [f"{FT}/pdf_round/pdf_decisions.csv"] + sorted(_g.glob(f"{FT}/pdf_round2/b*/decisions_b*.csv"))]
_prs = [x for x in _prs if x is not None]
PR = pd.concat(_prs, ignore_index=True).drop_duplicates("record_id", keep="last") if _prs else None
if PR is not None:
    for d in (ft1, ft2, inc):
        if d is None: continue
        m = d.record_id.isin(PR.record_id)
        for rid in d.loc[m, "record_id"]:
            r = PR[PR.record_id == rid].iloc[0]
            d.loc[d.record_id == rid, "_dec"] = r["_dec"]
            if "code" in d and "code" in PR: d.loc[d.record_id == rid, "code"] = r["code"]
def cnt(d):
    if d is None: return None
    inc_ = d._dec.str.contains("include").sum(); unr = d._dec.str.contains("unresolved").sum(); exc = len(d) - inc_ - unr
    code = [c for c in d.columns if c.lower() in ("code", "exclusion_code", "ft_code")]
    codes = d[d._dec.str.startswith("exclude")][code[0]].value_counts().to_dict() if code else {}
    return dict(n=int(len(d)), include=int(inc_), exclude=int(exc), unresolved=int(unr), exclusion_codes={str(k): int(v) for k, v in codes.items()})
C = dict(pending_round1=cnt(ft1), reforwarded_round2=cnt(ft2), abstract_includes_verified=cnt(inc))
assessed = sum(c["n"] for c in C.values() if c); included = sum(c["include"] for c in C.values() if c)
unresolved = sum(c["unresolved"] for c in C.values() if c); excluded_ft = sum(c["exclude"] for c in C.values() if c)
V2 = dict(date="2026-09-25", identification=P["identification"], screened=P["screening"]["records_screened"],
          dual_screened=A["n_compared"], kappa_binary=A["all"]["kappa_binary"], kappa_3class=A["all"]["kappa_3class"],
          pct_agreement=A["all"]["pct_agreement_binary"], disagreements=A["n_disagreements"],
          excluded_title_abstract=P["screening"]["records_screened"] - assessed - (P["eligibility"]["excluded"]),
          excluded_preprint_duplicate_at_eligibility=P["eligibility"]["excluded"],
          fulltext_sought=assessed, fulltext_not_retrieved=unresolved, fulltext_assessed=assessed - unresolved,
          fulltext_excluded=excluded_ft, included=included, included_other_methods=P["included"]["other_methods"], stages=C)
json.dump(V2, open(f"{L}/prisma_counts_v2.json", "w"), indent=1); print(json.dumps(V2, indent=1))
# flow diagram
fig, ax = plt.subplots(figsize=(9, 10)); ax.axis("off"); ax.set_xlim(0, 10); ax.set_ylim(0, 12)
def box(x, y, w, h, t, fc="#eef3fa"):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.05", fc=fc, ec="#444")); ax.text(x + w / 2, y + h / 2, t, ha="center", va="center", fontsize=8.5)
def arr(x1, y1, x2, y2): ax.annotate("", xy=(x2, y2), xytext=(x1, y1), arrowprops=dict(arrowstyle="->", color="#444"))
I = P["identification"]
box(0.3, 10.4, 5.2, 1.3, f"Records identified\nPubMed (10 queries): {I['pubmed_hits_total_across_queries']:,} hits\nEurope PMC preprints: {I['europepmc_preprint_hits']}")
box(6.2, 10.4, 3.5, 1.3, f"Duplicates removed\nn = {I['duplicates_removed']:,}", "#f7f7f7")
box(0.3, 8.3, 5.2, 1.5, f"Records screened (title/abstract)\nn = {V2['screened']:,}\nDual screened: {V2['dual_screened']:,} (κ = {V2['kappa_binary']:.2f})")
box(6.2, 8.3, 3.5, 1.5, f"Excluded at title/abstract\nn = {V2['excluded_title_abstract']:,}\n(preprints/duplicates at\neligibility: {V2['excluded_preprint_duplicate_at_eligibility']})", "#f7f7f7")
box(0.3, 6.3, 5.2, 1.3, f"Reports sought for retrieval\nn = {V2['fulltext_sought']}")
box(6.2, 6.3, 3.5, 1.3, f"Reports not retrieved\n(no open access)\nn = {V2['fulltext_not_retrieved']}", "#f7f7f7")
box(0.3, 4.2, 5.2, 1.4, f"Reports assessed for eligibility\nn = {V2['fulltext_assessed']}")
codes = {}
for c in C.values():
    if c:
        for k, v in c["exclusion_codes"].items(): codes[k] = codes.get(k, 0) + v
box(6.2, 3.9, 3.5, 2.0, "Reports excluded\nn = %d\n" % V2["fulltext_excluded"] + "\n".join(f"{k}: {v}" for k, v in sorted(codes.items(), key=lambda x: -x[1])[:6]), "#f7f7f7")
box(0.3, 1.8, 5.2, 1.5, f"Studies included\nn = {V2['included']} (database search)\n+ {V2['included_other_methods']} via other methods", "#e6f2e6")
arr(2.9, 10.4, 2.9, 9.8); arr(5.5, 11.05, 6.2, 11.05); arr(2.9, 8.3, 2.9, 7.6); arr(5.5, 9.05, 6.2, 9.05); arr(2.9, 6.3, 2.9, 5.6); arr(5.5, 6.95, 6.2, 6.95); arr(2.9, 4.2, 2.9, 3.3); arr(5.5, 4.9, 6.2, 4.9)
plt.savefig("04_results/figures/pub/FigS1_PRISMA.png", dpi=300, bbox_inches="tight"); print("saved")
