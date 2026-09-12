#!/usr/bin/env python3
"""Whole-process Q23 convergence; no claim from output spacing alone."""
import json
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np

from scripts.question3 import read_case, run_case, fingerprint
from src.data_io import ROOT, digest, save_json


def get_case(label, **kwargs):
    try:
        result, meta = read_case(label)
        for key, value in kwargs.items():
            mapped = {"parent":"q2_parent_label", "max_step":"max_step_s", "output_step":"output_step_s",
                      "scenario":"boundary_scenario"}.get(key, key)
            if meta.get(mapped) != value:
                raise ValueError("Parameters changed")
        return result, meta
    except (ValueError, FileNotFoundError):
        return run_case(label=label, **kwargs)


def difference(candidate, reference):
    # Compare only common regular output times, never interpolate across events.
    ts, ic, ir = np.intersect1d(candidate["time_s"], reference["time_s"], return_indices=True)
    out = {"common_output_count": len(ts), "fields": {}}
    for f in ("temperature", "moisture"):
        delta = np.abs(candidate[f][ic]-reference[f][ir])
        where = np.unravel_index(np.argmax(delta), delta.shape)
        out["fields"][f] = {"max_abs":float(delta[where]), "time_s":float(ts[where[0]]),
            "radius_cm":float(candidate["radius_m"][where[1]]*100),
            "four_decimal_rounding_differences":int(np.count_nonzero(np.round(candidate[f][ic],4)!=np.round(reference[f][ir],4)))}
    return out


def main():
    report = {"source_fingerprint":fingerprint(), "script_sha256":digest(__file__),
              "complete":False, "mesh_series":[], "temporal":[], "is_strict_error_bound":False}
    outputs = {}
    previous = None
    for n, parent in ((320,"q2_grid_320"),(640,"q2_grid_640"),(1280,"q2_final"),(2560,"q2_grid_2560")):
        label = "q3_final" if n == 1280 else f"q3_mesh{n}"
        c, m = get_case(label, parent=parent, max_step=30., rtol=1e-10,
                        method="BDF", output_step=1. if n >= 640 else 60., scenario="main")
        assert m["statistics"]["event_found"]
        row = {"label":label, "n":n, "event_time_s":m["statistics"]["event_time_s"]}
        if previous is not None:
            row["event_delta_s"] = row["event_time_s"]-previous[1]["statistics"]["event_time_s"]
            row["difference_from_previous"] = difference(c, previous[0])
        outputs[n] = (c,m)
        previous = (c,m)
        report["mesh_series"].append(row)
        save_json(ROOT/"reports/q3_convergence.json", report)
    baseline, bm = outputs[1280]
    cases = [("q3_cap60",60.,1e-10,"BDF"), ("q3_cap15",15.,1e-10,"BDF"),
             ("q3_radau",30.,1e-10,"Radau"), ("q3_tol6",60.,1e-6,"BDF"),
             ("q3_tol7",60.,1e-7,"BDF"), ("q3_tol8",60.,1e-8,"BDF")]
    for label, cap, tolerance, method in cases:
        c,m = get_case(label, parent="q2_final", max_step=cap, rtol=tolerance,
                       method=method, output_step=1. if label in ("q3_cap60","q3_cap15","q3_radau") else 60., scenario="main")
        row = {"label":label, "max_step_s":cap, "rtol":tolerance, "method":method,
               "event_delta_s":m["statistics"]["event_time_s"]-bm["statistics"]["event_time_s"],
               "difference_from_baseline":difference(c,baseline), "statistics":m["statistics"]}
        if label in ("q3_cap60", "q3_cap15"):
            assert abs(m["statistics"]["accepted_max_step_s"]-cap) < 1e-7
            assert m["statistics"]["accepted_steps"] != bm["statistics"]["accepted_steps"]
        report["temporal"].append(row)
        save_json(ROOT/"reports/q3_convergence.json", report)
    report["baseline_event_time_s"] = bm["statistics"]["event_time_s"]
    report["baseline_statistics"] = bm["statistics"]
    report["complete"] = True
    save_json(ROOT/"reports/q3_convergence.json", report)
    print(json.dumps({"mesh":report["mesh_series"],"temporal_event_deltas":[[r["label"],r["event_delta_s"]] for r in report["temporal"]]}),flush=True)


if __name__ == "__main__":
    main()
