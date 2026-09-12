#!/usr/bin/env python3
"""Q2 single-factor parameter scenarios and explicitly controlled freezes."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from scripts.question2 import run_case, fingerprint
from src.data_io import ROOT, digest, save_json


def main():
    report = {"source_fingerprint": fingerprint(), "script_sha256": digest(__file__),
              "interpretation": "scenario contrasts, not confidence intervals or experimental causality",
              "cases": [], "complete": False}
    scenarios = [(f"q2_{key}_{multiplier:g}", {key: multiplier})
                 for key in ("h_multiplier", "hm_multiplier", "D_multiplier", "beta_multiplier")
                 for multiplier in (.9, 1.1)]
    scenarios += [("q2_freeze_D_temperature", {"freeze_D_temperature": True}),
                  ("q2_freeze_thermal_moisture", {"freeze_thermal_moisture": True})]
    for label, options in scenarios:
        cache, meta = run_case(n=1280, rtol=1e-10, atol_temperature=1e-11, atol_moisture=1e-13,
                               max_step=2., label=label, perturbation=options)
        report["cases"].append({"label": label, "options": options,
            "temperature_centre_3h": float(cache["temperature"][-1, 0]),
            "temperature_surface_3h": float(cache["temperature"][-1, -1]),
            "moisture_centre_3h": float(cache["moisture"][-1, 0]),
            "moisture_surface_3h": float(cache["moisture"][-1, -1]),
            "moisture_mean_3h": float(cache["moisture_mean"][-1]), "statistics": meta["statistics"]})
        save_json(ROOT/"reports/q2_sensitivity.json", report)
    report["complete"] = True
    save_json(ROOT/"reports/q2_sensitivity.json", report)


if __name__ == "__main__":
    main()
