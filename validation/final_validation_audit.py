#!/usr/bin/env python3
"""Read-only audit of the final validation bundle and its stated scope."""
from pathlib import Path
import csv
import json
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from src.data_io import ROOT, digest, save_json, verify_sources


def main():
    manifest = verify_sources()
    report_path = ROOT/"reports/final_result_validation.json"
    markdown_path = ROOT/"reports/final_result_validation.md"
    table_path = ROOT/"tables/final_result_validation_summary.csv"
    holdout_path = ROOT/"reports/input_holdout_validation.json"
    report = json.loads(report_path.read_text(encoding="utf-8"))
    holdout = json.loads(holdout_path.read_text(encoding="utf-8"))

    assert report["status"] == "conditionally_accepted_for_figure_generation_under_model_contract"
    assert report["script_sha256"] == digest(ROOT/"validation/final_result_validation.py")
    assert report["original_inputs_unchanged"] == len(manifest["files"]) == 8
    assert report["prior_evidence_integrity"]["status"].startswith("passed_")
    assert report["dimension_audit"]["status"].startswith("passed_")
    assert report["bounds_and_events"]["status"].startswith("passed_")
    assert report["conservation"]["status"].startswith("passed_")
    assert report["extreme_and_counterexample_tests"]["status"].startswith("passed_")
    assert holdout["status"].startswith("passed_")
    assert all(row["validation_role"].endswith("not_ground_truth") for row in report["source_templates"].values())
    assert all(row["status"].startswith("passed_") for row in report["workbooks"].values())
    assert all(row["status"].startswith("passed_") for row in report["independent_references"])

    numeric = sum(row["numeric_result_cells_checked"] for row in report["workbooks"].values())
    blanks = sum(row["outside_domain_blank_cells_checked"] for row in report["workbooks"].values())
    assert numeric == 8_882_410 and blanks == 23_165
    for name, row in report["workbooks"].items():
        assert digest(ROOT/row["path"]) == row["sha256"], name
    for name, expected in report["full_precision_artifacts_sha256"].items():
        assert digest(ROOT/"results/final"/name) == expected, name
    for relative, expected in report["artifacts_sha256"].items():
        assert digest(ROOT/relative) == expected, relative
    for relative, expected in holdout["artifacts_sha256"].items():
        assert digest(ROOT/relative) == expected, relative

    with table_path.open(encoding="utf-8-sig", newline="") as handle:
        rows = list(csv.DictReader(handle))
    expected_rows = 3*len(report["workbooks"])
    expected_rows += 2*sum(len(row["comparisons"]) for row in report["independent_references"])
    expected_rows += 2*len(report["conservation"]["rows"])
    assert len(rows) == expected_rows
    assert all(row["status"].startswith("passed_") for row in rows)
    prose = markdown_path.read_text(encoding="utf-8")
    for phrase in ("条件验收通过", "未当作标准答案", "不是置信区间", "可以开始图表生成"):
        assert phrase in prose

    output = {
        "status": "passed_final_validation_bundle_readback",
        "script_sha256": digest(__file__),
        "original_inputs_unchanged": len(manifest["files"]),
        "numeric_workbook_cells_reconciled": numeric,
        "outside_domain_blank_cells_reconciled": blanks,
        "summary_csv_rows_reconciled": len(rows),
        "verified_artifacts_sha256": {
            str(report_path.relative_to(ROOT)): digest(report_path),
            str(markdown_path.relative_to(ROOT)): digest(markdown_path),
            str(table_path.relative_to(ROOT)): digest(table_path),
            str(holdout_path.relative_to(ROOT)): digest(holdout_path),
            **{row["path"]: row["sha256"] for row in report["workbooks"].values()},
        },
        "scope": "bundle consistency and evidence readback; not new experimental validation",
    }
    save_json(ROOT/"reports/final_validation_audit.json", output)
    print(json.dumps(output, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
