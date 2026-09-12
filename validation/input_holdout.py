#!/usr/bin/env python3
"""Hold out measured input knots to audit interpolation, not the PDE prediction.

The problem supplies environmental and radius drivers but no measurements inside
the material.  These checks therefore quantify interpolation fidelity only and
must not be described as experimental validation of the heat/moisture fields.
"""
from pathlib import Path
import csv
import json
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np
from openpyxl import load_workbook
from scipy.interpolate import PchipInterpolator

from src.data_io import ROOT, digest, load_config, save_json, verify_sources


def read_numeric_sheet(path):
    workbook = load_workbook(path, read_only=True, data_only=True)
    rows = list(workbook.active.values)
    workbook.close()
    header = list(rows[0])
    data = np.asarray(rows[1:], dtype=float)
    return header, data


def metrics(observed, predicted):
    residual = np.asarray(predicted) - np.asarray(observed)
    return {
        "count": int(residual.size),
        "bias": float(np.mean(residual)),
        "mae": float(np.mean(np.abs(residual))),
        "rmse": float(np.sqrt(np.mean(residual**2))),
        "max_abs": float(np.max(np.abs(residual))),
        "p95_abs": float(np.quantile(np.abs(residual), 0.95)),
    }


def leave_one_out(times, values, method):
    predicted = np.empty(len(times)-2)
    for out, i in enumerate(range(1, len(times)-1)):
        keep = np.ones(len(times), dtype=bool)
        keep[i] = False
        if method == "linear":
            predicted[out] = np.interp(times[i], times[keep], values[keep])
        else:
            predicted[out] = PchipInterpolator(times[keep], values[keep], extrapolate=False)(times[i])
    return values[1:-1], predicted


def interleaved_holdout(times, values, method):
    # Non-adjacent interior points: every fifth record beginning at index 5.
    held = np.arange(5, len(times)-1, 5)
    keep = np.ones(len(times), dtype=bool)
    keep[held] = False
    if method == "linear":
        predicted = np.interp(times[held], times[keep], values[keep])
    else:
        predicted = PchipInterpolator(times[keep], values[keep], extrapolate=False)(times[held])
    return held, values[held], predicted


def audit_series(name, times, values, unit):
    value_range = float(np.ptp(values))
    result = {
        "name": name,
        "unit": unit,
        "samples": int(len(times)),
        "time_start_s": float(times[0]),
        "time_end_s": float(times[-1]),
        "time_strictly_increasing": bool(np.all(np.diff(times) > 0)),
        "time_step_s": {
            "minimum": float(np.min(np.diff(times))),
            "maximum": float(np.max(np.diff(times))),
            "uniform": bool(np.ptp(np.diff(times)) == 0),
        },
        "value_min": float(np.min(values)),
        "value_max": float(np.max(values)),
        "finite": bool(np.isfinite(values).all()),
        "leave_one_out": {},
        "interleaved_every_fifth": {},
    }
    for method in ("linear", "pchip"):
        observed, predicted = leave_one_out(times, values, method)
        item = metrics(observed, predicted)
        item["max_abs_fraction_of_observed_range"] = item["max_abs"]/value_range
        result["leave_one_out"][method] = item
        held, observed, predicted = interleaved_holdout(times, values, method)
        item = metrics(observed, predicted)
        item["held_indices"] = held.tolist()
        item["max_abs_fraction_of_observed_range"] = item["max_abs"]/value_range
        result["interleaved_every_fifth"][method] = item

    held = np.arange(5, len(times)-1, 5)
    persistence = values[held-1]
    baseline = metrics(values[held], persistence)
    baseline["definition"] = "previous_observation_persistence"
    result["simple_baseline"] = baseline
    linear_rmse = result["interleaved_every_fifth"]["linear"]["rmse"]
    result["linear_rmse_improvement_over_persistence_fraction"] = float(
        1-linear_rmse/baseline["rmse"]
    )
    return result


def main():
    manifest = verify_sources()
    config = load_config()
    source = Path(config["source_directory"])
    env_header, env = read_numeric_sheet(source/"附件1.xlsx")
    radius_header, radius = read_numeric_sheet(source/"附件2.xlsx")
    assert env_header == ["时间", "温度", "水分浓度"]
    assert radius_header == ["时间", "半径"]

    series = [
        audit_series("environment_temperature", env[:, 0], env[:, 1], "degC"),
        audit_series("environment_moisture_input", env[:, 0], env[:, 2], "kg/kg_air_basis_as_supplied"),
        audit_series("observed_radius", radius[:, 0], radius[:, 1], "cm"),
    ]
    assert all(row["finite"] and row["time_strictly_increasing"] for row in series)
    assert np.all(radius[:, 1] > 0) and np.all(np.diff(radius[:, 1]) <= 0)

    output = {
        "status": "passed_input_integrity_and_interpolation_holdout_checks",
        "script_sha256": digest(__file__),
        "original_inputs_unchanged": len(manifest["files"]),
        "series": series,
        "radius_monotone_nonincreasing": True,
        "boundary_contract": {
            "observed_environment_support_s": [float(env[0, 0]), float(env[-1, 0])],
            "post_observation_environment_rule": "mean_10800_14400_inclusive_then_hold",
            "post_observation_rule_has_measurement_holdout_validation": False,
            "observed_radius_support_s": [float(radius[0, 0]), float(radius[-1, 0])],
            "radius_main_interpolation": "linear_no_extrapolation",
            "air_to_material_map": "C_eq=beta*Y_is_a_boundary_closure_assumption_not_a_unit_identity",
        },
        "scope": {
            "what_is_validated": "fidelity of interpolation between measured input knots",
            "what_is_not_validated": "internal temperature/moisture predictions or post-4h environment extension",
            "why_no_supervised_k_fold": "the attachments contain drivers, not internal material response labels",
        },
    }
    report_path = ROOT/"reports/input_holdout_validation.json"
    save_json(report_path, output)

    table_path = ROOT/"tables/input_holdout_metrics.csv"
    table_path.parent.mkdir(exist_ok=True)
    with table_path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["series", "design", "method", "count", "bias", "mae", "rmse", "max_abs", "p95_abs", "unit"])
        for row in series:
            for design in ("leave_one_out", "interleaved_every_fifth"):
                for method, values in row[design].items():
                    writer.writerow([row["name"], design, method, values["count"], values["bias"],
                                     values["mae"], values["rmse"], values["max_abs"], values["p95_abs"], row["unit"]])
            values = row["simple_baseline"]
            writer.writerow([row["name"], "interleaved_every_fifth", "persistence", values["count"], values["bias"],
                             values["mae"], values["rmse"], values["max_abs"], values["p95_abs"], row["unit"]])
    output["artifacts_sha256"] = {
        str(table_path.relative_to(ROOT)): digest(table_path),
    }
    save_json(report_path, output)
    print(json.dumps({"status": output["status"], "series": [
        {"name": r["name"], "linear_loo_rmse": r["leave_one_out"]["linear"]["rmse"],
         "linear_interleaved_rmse": r["interleaved_every_fifth"]["linear"]["rmse"]}
        for r in series]}, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
