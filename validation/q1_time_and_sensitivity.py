#!/usr/bin/env python3
"""Q1 tolerance, binding step caps, time-integrator and face-flux crosschecks.

These comparisons are not an independent nonlinear spatial solver.
"""
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))

import numpy as np
from scripts.question1 import run_case
from src.data_io import ROOT, save_json

FIELDS=("temperature","moisture")
PAPER_TIMES=[100,300,600,900,1200,1500,1800]
PAPER_RADII=[0,5,10,15,20]


def compare(candidate,reference):
    out={}
    for field in FIELDS:
        delta=np.abs(candidate[field][1:]-reference[field][1:])
        where=np.unravel_index(np.argmax(delta),delta.shape)
        out[field]={"max_abs":float(delta[where]),"time_s":int(where[0]+1),
                    "radius_cm":round(float(candidate["radius_m"][where[1]]*100),5),
                    "rounding_differences":int(np.count_nonzero(np.round(candidate[field][1:],4)!=np.round(reference[field][1:],4))),
                    "paper_rounding_differences":int(np.count_nonzero(
                        np.round(candidate[field][np.ix_(PAPER_TIMES,PAPER_RADII)],4)!=
                        np.round(reference[field][np.ix_(PAPER_TIMES,PAPER_RADII)],4)))}
    return out


def main():
    _,reference,refmeta=run_case(1280,2.,rtol=1e-10,atol=1e-12,max_step=1.,label="q1_final")
    report={"reference":"q1_final","scope":"Q1 local numerical acceptance, not full independent nonlinear spatial validation",
            "comparisons":[],"step_cap_series":[],"sensitivities":[]}
    cases=[("tol_1e-6",{"rtol":1e-6,"atol":1e-8}),
           ("tol_1e-7",{"rtol":1e-7,"atol":1e-9}),
           ("tol_1e-8",{"rtol":1e-8,"atol":1e-10}),
           ("cap_60",{}),("cap_2",{"max_step":2.}),
           ("cap_0p5",{"max_step":.5}),
           ("radau",{"method":"Radau","max_step":1.}),
           ("kirchhoff",{"scheme":"kirchhoff","max_step":1.}),
           ("fine2560",{"n":2560,"max_step":1.})]
    step_results={1.:(reference,refmeta)}
    for label,options in cases:
        args=dict(n=1280,grading=2.,rtol=1e-10,atol=1e-12,max_step=60.,label=label,store_cells=False)
        args.update(options)
        _,cache,meta=run_case(**args)
        row={"label":label,"comparison_to_reference":compare(cache,reference),"fields":meta["fields"]}
        print(row,flush=True)
        report["comparisons"].append(row)
        if label in ("cap_2","cap_0p5"):
            step_results[options["max_step"]]=(cache,meta)
        save_json(ROOT/"reports/q1_time_and_sensitivity.json",report)
    for coarse,fine in [(2.,1.),(1.,.5)]:
        report["step_cap_series"].append({"coarse_s":coarse,"fine_s":fine,
            "coarse_fields":step_results[coarse][1]["fields"],"fine_fields":step_results[fine][1]["fields"],
            "difference":compare(step_results[coarse][0],step_results[fine][0])})
    for parameter in ("h_multiplier","hm_multiplier","D_multiplier","beta_multiplier"):
        for multiplier in (.9,1.1):
            label=f"sensitivity_{parameter}_{multiplier:g}"
            _,cache,meta=run_case(1280,2.,rtol=1e-10,atol=1e-12,max_step=1.,label=label,
                                  store_cells=False,perturbation={parameter:multiplier})
            row={"parameter":parameter,"multiplier":multiplier,
                 "scenario_not_statistical_confidence_interval":True,
                 "comparison_to_baseline":compare(cache,reference)}
            for field in FIELDS:
                row[field+"_centre_1800"]=float(cache[field][-1,0])
                row[field+"_surface_1800"]=float(cache[field][-1,-1])
                row[field+"_mean_1800"]=float(cache[field+"_mean"][-1])
            report["sensitivities"].append(row)
            save_json(ROOT/"reports/q1_time_and_sensitivity.json",report)
    report["complete"]=True
    save_json(ROOT/"reports/q1_time_and_sensitivity.json",report)


if __name__=="__main__":
    main()
