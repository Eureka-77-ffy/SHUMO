#!/usr/bin/env python3
"""Propagate Q2 parameter/freeze scenarios and boundary alternatives to Q3."""
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))

from scripts.question3 import fingerprint
from src.data_io import ROOT, digest, save_json
from validation.q3_convergence import get_case


def main():
    report={"source_fingerprint":fingerprint(),"script_sha256":digest(__file__),
            "complete":False,"cases":[],"interpretation":"controlled_scenarios_not_statistical_confidence_intervals"}
    cases=[]
    for key in ("h_multiplier","hm_multiplier","D_multiplier","beta_multiplier"):
        for value in (.9,1.1):
            cases.append((f"q3_{key}_{value:g}",f"q2_{key}_{value:g}","main"))
    for key in ("freeze_D_temperature","freeze_thermal_moisture"):
        cases.append((f"q3_{key}",f"q2_{key}","main"))
    for scenario in ("last_observation","mean_9000_14400"):
        cases.append((f"q3_{scenario}","q2_final",scenario))
    for label,parent,scenario in cases:
        c,m=get_case(label,parent=parent,max_step=30.,rtol=1e-10,method="BDF",output_step=60.,scenario=scenario)
        report["cases"].append({"label":label,"parent":parent,"perturbation":m["perturbation"],
            "boundary_scenario":scenario,"event_time_s":m["statistics"]["event_time_s"],
            "moisture_mean_at_event":float(c["moisture_mean"][-1]),
            "moisture_surface_at_event":float(c["moisture"][-1,-1]),"statistics":m["statistics"]})
        save_json(ROOT/"reports/q3_sensitivity.json",report)
    report["complete"]=True
    save_json(ROOT/"reports/q3_sensitivity.json",report)


if __name__=="__main__":
    main()
