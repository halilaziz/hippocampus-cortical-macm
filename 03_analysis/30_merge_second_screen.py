"""Compare blinded second screening with the first screening (agreement, Cohen's kappa), list disagreements for human
adjudication, and fold full-text decisions into updated PRISMA counts."""
import json, glob, os, pandas as pd, numpy as np
from sklearn.metrics import cohen_kappa_score, confusion_matrix
L = "01_literature"; SS = f"{L}/second_screen"
first = pd.read_csv(f"{L}/screening.csv", low_memory=False)
sec = pd.concat([pd.read_json(f, lines=True) for f in sorted(glob.glob(f"{SS}/batch_*_decisions.jsonl"))], ignore_index=True).drop_duplicates("record_id", keep="last")
m = first.merge(sec[["record_id", "decision", "code", "reason"]].rename(columns={"decision": "second_decision", "code": "second_code", "reason": "second_reason"}), on="record_id", how="inner")
m["first_decision"] = m.screen_decision.fillna("exclude").str.lower().replace({"included": "include", "excluded": "exclude"})
m["stratum"] = np.where(m.screen_method.str.contains("manual", na=False), "abstract-read", "rule-excluded sample")
bin_ = lambda x: np.where(x.isin(["include", "maybe"]), "forward", "exclude")
res = {"n_compared": int(len(m))}
for st, g in [("all", m)] + list(m.groupby("stratum")):
    a, b = bin_(g.first_decision), bin_(g.second_decision)
    r = dict(n=int(len(g)), pct_agreement_binary=round(100 * float((a == b).mean()), 1))
    if len(set(a) | set(b)) > 1: r["kappa_binary"] = round(float(cohen_kappa_score(a, b)), 3)
    if st != "rule-excluded sample":
        r["kappa_3class"] = round(float(cohen_kappa_score(g.first_decision, g.second_decision)), 3)
    r["first_forward_second_exclude"] = int(((a == "forward") & (b == "exclude")).sum())
    r["first_exclude_second_forward"] = int(((a == "exclude") & (b == "forward")).sum())
    res[st] = r
cm = pd.crosstab(m.first_decision, m.second_decision); cm.to_csv(f"{SS}/confusion_matrix.csv")
dis = m[bin_(m.first_decision) != bin_(m.second_decision)][["record_id", "pmid", "first_author", "year", "title", "stratum", "first_decision", "screen_code", "screen_reason",
                                                               "second_decision", "second_code", "second_reason", "final_status"]]
dis.to_csv(f"{SS}/disagreements_for_adjudication.csv", index=False)
res["n_disagreements"] = int(len(dis))
# full-text stage
ft_path = f"{SS}/fulltext/fulltext_decisions.csv"
if os.path.exists(ft_path):
    ft = pd.read_csv(ft_path)
    dcol = [c for c in ft.columns if c.lower() in ("decision", "ft_decision", "final_decision")][0]
    res["fulltext"] = {str(k): int(v) for k, v in ft[dcol].value_counts().items()}
    ccol = [c for c in ft.columns if c.lower() in ("code", "exclusion_code", "ft_code")]
    if ccol: res["fulltext_exclusion_codes"] = {str(k): int(v) for k, v in ft[ft[dcol].astype(str).str.startswith("exclude")][ccol[0]].value_counts().items()}
    P = json.load(open(f"{L}/prisma_counts.json"))
    inc_ft = int((ft[dcol] == "include").sum()); unresolved = int(ft[dcol].astype(str).str.contains("unresolved").sum())
    res["prisma_updated"] = dict(included_before=P["included"]["total_included"], included_from_fulltext=inc_ft,
                                 total_included_after=P["included"]["total_included"] + inc_ft, still_unresolved_no_access=unresolved)
json.dump(res, open(f"{SS}/agreement_summary.json", "w"), indent=1); print(json.dumps(res, indent=1)); print(cm)
