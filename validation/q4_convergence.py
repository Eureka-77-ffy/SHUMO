#!/usr/bin/env python3
"""Spatial and binding time refinement of Q4, including moving output masks."""
import json
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import numpy as np
from scripts.question4 import read_case,run_case,fingerprint
from src.data_io import ROOT,digest,save_json


def get_case(label,**options):
    try:
        c,m=read_case(label)
        aliases={"group":"property_set","extension":"radius_extension","scenario":"boundary_scenario",
                 "max_step":"max_step_s","limit":"requested_limit_s"}
        for k,v in options.items():
            if m.get(aliases.get(k,k))!=v:raise ValueError("Parameters changed")
        return c,m
    except (ValueError,FileNotFoundError):return run_case(label=label,**options)


def compare(c,cm,r,rm):
    ts,i,j=np.intersect1d(c["time_s"],r["time_s"],return_indices=True)
    result={"common_times":len(ts),"event_delta_s":cm["statistics"]["event_time_s"]-rm["statistics"]["event_time_s"],"fields":{}}
    for f in ("temperature","moisture"):
        a,b=c[f][i],r[f][j]
        assert np.array_equal(np.isnan(a),np.isnan(b))
        delta=np.abs(a-b)
        loc=np.unravel_index(np.nanargmax(delta),delta.shape)
        finite=np.isfinite(a)
        result["fields"][f]={"max_abs":float(delta[loc]),"time_s":float(ts[loc[0]]),
            "position": "surface" if loc[1]==21 else f"{loc[1]/10:g} cm",
            "four_decimal_rounding_differences":int(np.count_nonzero(np.round(a[finite],4)!=np.round(b[finite],4)))}
    return result


def main():
    report={"source_fingerprint":fingerprint(),"script_sha256":digest(__file__),"mesh_series":[],"temporal":[],"complete":False}
    prior=None
    for n in (160,320,640,1280,2560):
        label="q4_final" if n==1280 else f"q4_grid{n}"
        c,m=get_case(label,n=n,max_step=30.,rtol=1e-10,geometry="linear",group="q4")
        row={"label":label,"n":n,"event_time_s":m["statistics"]["event_time_s"]}
        if prior is not None:row["difference_from_previous"]=compare(c,m,*prior)
        prior=(c,m)
        report["mesh_series"].append(row)
        save_json(ROOT/"reports/q4_convergence.json",report)
    ref,rm=read_case("q4_final")
    for label,cap,rtol,method in (("q4_cap60",60.,1e-10,"BDF"),("q4_cap15",15.,1e-10,"BDF"),
            ("q4_radau",30.,1e-10,"Radau"),("q4_tol6",60.,1e-6,"BDF"),
            ("q4_tol7",60.,1e-7,"BDF"),("q4_tol8",60.,1e-8,"BDF")):
        c,m=get_case(label,n=1280,max_step=cap,rtol=rtol,method=method)
        if label in ("q4_cap60","q4_cap15"):
            assert abs(m["statistics"]["accepted_max_step_s"]-cap)<1e-7
            assert m["statistics"]["accepted_steps"]!=rm["statistics"]["accepted_steps"]
        report["temporal"].append({"label":label,"difference":compare(c,m,ref,rm),"statistics":m["statistics"]})
        save_json(ROOT/"reports/q4_convergence.json",report)
    report["complete"]=True
    save_json(ROOT/"reports/q4_convergence.json",report)
    print(json.dumps(report),flush=True)


if __name__=="__main__":main()
