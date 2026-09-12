#!/usr/bin/env python3
"""Q2 complete-export-grid refinement and binding temporal-cap comparisons."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np

from scripts.question2 import run_case, fingerprint
from src.data_io import ROOT, digest, save_json


def compare(candidate, reference):
    out = {}
    for field in ("temperature", "moisture"):
        delta = np.abs(candidate[field][1:]-reference[field][1:])
        where = np.unravel_index(np.argmax(delta), delta.shape)
        out[field] = {
            "max_abs": float(delta[where]), "time_s": float(candidate["time_s"][where[0]+1]),
            "radius_cm": float(candidate["radius_m"][where[1]]*100),
            "rounding_differences": int(np.count_nonzero(np.round(candidate[field][1:], 4) != np.round(reference[field][1:], 4))),
            "paper_rounding_differences": int(np.count_nonzero(
                np.round(candidate[field][np.ix_(np.arange(1800, 10801, 1800), [0, 5, 10, 15, 20])], 4)
                != np.round(reference[field][np.ix_(np.arange(1800, 10801, 1800), [0, 5, 10, 15, 20])], 4)))
        }
    return out


def main():
    report = {"source_fingerprint": fingerprint(), "script_sha256": digest(__file__),
              "scope": "Q2 first 3h; differences are not strict error bounds",
              "mesh_series": [], "temporal": [], "complete": False}
    baseline = baseline_meta = previous = None
    for n in (80, 160, 320, 640, 1280, 2560):
        label = "q2_final" if n == 1280 else f"q2_grid_{n}"
        cache, meta = run_case(n=n, rtol=1e-10, atol_temperature=1e-11, atol_moisture=1e-13, max_step=2., label=label)
        row = {"label": label, "n": n, "statistics": meta["statistics"]}
        if previous is not None:
            row["difference_from_previous"] = compare(cache, previous)
        if n == 1280:
            baseline, baseline_meta = cache, meta
        report["mesh_series"].append(row)
        previous = cache
        save_json(ROOT/"reports/q2_convergence.json", report)
    cases = [("q2_tol_1e6", {"rtol": 1e-6, "atol_temperature": 1e-8, "atol_moisture": 1e-10, "max_step": 60.}),
             ("q2_tol_1e7", {"rtol": 1e-7, "atol_temperature": 1e-9, "atol_moisture": 1e-11, "max_step": 60.}),
             ("q2_tol_1e8", {"rtol": 1e-8, "atol_temperature": 1e-10, "atol_moisture": 1e-12, "max_step": 60.}),
             ("q2_cap4", {"max_step": 4.}), ("q2_cap1", {"max_step": 1.}),
             ("q2_radau", {"method": "Radau"})]
    for label, options in cases:
        args = dict(n=1280, rtol=1e-10, atol_temperature=1e-11, atol_moisture=1e-13, max_step=2.)
        args.update(options)
        cache, meta = run_case(**args, label=label)
        report["temporal"].append({"label": label, "difference_from_baseline": compare(cache, baseline), "statistics": meta["statistics"]})
        save_json(ROOT/"reports/q2_convergence.json", report)
    report["baseline_statistics"] = baseline_meta["statistics"]
    report["complete"] = True
    save_json(ROOT/"reports/q2_convergence.json", report)
    print("Q2 mesh and temporal comparisons complete", flush=True)


if __name__ == "__main__":
    main()
