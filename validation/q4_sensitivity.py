#!/usr/bin/env python3
"""Q4 physical-parameter, radius-interpolation and environment scenarios."""
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from scripts.question4 import fingerprint
from src.data_io import ROOT,digest,save_json
from validation.q4_convergence import get_case


def main():
    report={"source_fingerprint":fingerprint(),"script_sha256":digest(__file__),"complete":False,"cases":[]}
    cases=[(f"q4_{key}_{v:g}",{"perturbation":{key:v}}) for key in
           ("h_multiplier","hm_multiplier","D_multiplier","beta_multiplier") for v in (.9,1.1)]
    cases += [("q4_pchip",{"geometry":"pchip"}),
              ("q4_last_observation",{"scenario":"last_observation"}),
              ("q4_mean_9000_14400",{"scenario":"mean_9000_14400"})]
    for label,options in cases:
        c,m=get_case(label,n=1280,max_step=30.,rtol=1e-10,**options)
        assert m["statistics"]["event_found"]
        report["cases"].append({"label":label,"options":options,"event_time_s":m["statistics"]["event_time_s"],
            "surface_moisture_at_event":float(c["moisture"][-1,-1]),"mean_moisture_at_event":float(c["moisture_mean"][-1]),
            "statistics":m["statistics"]})
        save_json(ROOT/"reports/q4_sensitivity.json",report)
    report["complete"]=True;save_json(ROOT/"reports/q4_sensitivity.json",report)


if __name__=="__main__":main()
