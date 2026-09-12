"""Bind the current scoped M4 reports without declaring the full paper done."""
import json
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import csv
import numpy as np
from src.data_io import ROOT,digest,save_json,verify_sources
from validation.endface_run import fingerprint


def main():
    inputs=verify_sources();prior=0
    for q in (1,2,3,4):
        d=json.loads((ROOT/f"reports/q{q}_export_audit.json").read_text())
        for p,h in d["artifacts_sha256"].items():assert digest(ROOT/p)==h,p;prior+=1
    geo=json.loads((ROOT/"reports/endface_audit.json").read_text())
    budget=json.loads((ROOT/"reports/unified_error_budget.json").read_text())
    assert geo["fingerprint"]==fingerprint() and geo["complete"] and budget["complete"]
    assert budget["script_sha256"]==digest(ROOT/"validation/unified_error_budget.py")
    for p,h in budget["evidence_sha256"].items():assert digest(ROOT/p)==h,p
    for p,h in geo["cache_artifacts_sha256"].items():assert digest(ROOT/p)==h,p
    assert len(geo["cases"])==17 and len(geo["cache_artifacts_sha256"])==58
    assert all(r["water_balance_relative"]<1e-6 for r in geo["cases"])
    with (ROOT/"tables/unified_error_budget.csv").open(encoding="utf-8-sig",newline="") as f:
        rows=list(csv.DictReader(f))
    assert len(rows)==len(budget["records"])
    for a,b in zip(rows,budget["records"]):
        for k,v in b.items():
            if k=="value":assert float(a[k])==v and np.isfinite(v)
            else:assert a[k]==v
    assert budget["confidence_interval"] is None and not budget["root_sum_of_squares_permitted"]
    assert len(budget["unquantified"])==7 and all(r["value"] is None for r in budget["unquantified"])
    paths=["reports/endface_assessment.md","reports/error_budget.md","reports/endface_audit.json",
        "reports/endface_kernel_checks.json","reports/endface_diagnostics.json","reports/q1_nonlinear_independent.json",
        "reports/unified_error_budget.json","tables/endface_comparison.csv","tables/unified_error_budget.csv",
        "validation/endface_suite.py","validation/m4_scope_audit.py","README.md","执行计划.md","reports/model_contract.md"]
    report={"status":"passed_requested_endface_assessment_and_unified_error_budget_not_full_paper_completion",
        "script_sha256":digest(__file__),"original_inputs_unchanged":len(inputs["files"]),"prior_audited_artifacts_unchanged":prior,
        "geometry_cached_runs":29,"matched_geometry_contrasts":17,"budget_rows_read_back":len(rows),
        "unquantified_model_categories":7,"primary_workbooks_changed":False,
        "completed_scope":["conditional end-face geometry audit with paired spatial/time comparisons",
            "Q1 original nonlinear independent weak-form sample check","unified traceable error ledger and scoped written assessment"],
        "not_completed_or_not_claimed":["experimental physical validation","full 2D absolute solution at submission precision",
            "universal error bound across end exposure/anisotropy/sorption assumptions","full manuscript figures references and final review"],
        "artifacts_sha256":{p:digest(ROOT/p) for p in paths}}
    save_json(ROOT/"reports/m4_validation_audit.json",report)
    print(json.dumps({k:v for k,v in report.items() if k!="artifacts_sha256"},ensure_ascii=False,indent=2))


if __name__=="__main__":main()
