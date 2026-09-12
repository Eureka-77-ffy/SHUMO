#!/usr/bin/env python3
"""Refine independent FE space and time, including the actual drying event."""
import json
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import numpy as np

from scripts.question3 import read_case, fingerprint
from src.data_io import ROOT, digest, save_json, verify_sources
from validation.q3_independent_fem import run, fingerprint as fem_fingerprint


def get_case(n, dt):
    label=f"q3_fem_n{n}_dt{dt:g}_initial0.9375"
    folder=ROOT/"results/cache/q3_independent"
    try:
        meta=json.loads((folder/(label+".json")).read_text())
        path=folder/(label+".npz")
        assert meta["source_fingerprint"]==fem_fingerprint() and meta["cache_sha256"]==digest(path)
        assert meta["initial_parent_sha256"]==digest(ROOT/"results/cache/q2_independent"/(meta["initial_parent_label"]+".npz"))
        with np.load(path) as raw:
            return {k:raw[k] for k in raw.files},meta
    except (FileNotFoundError,AssertionError):
        return run(n,dt)


def main():
    verify_sources()
    reference,rm=read_case("q3_final")
    report={"source_fingerprint":fingerprint(),"independent_fingerprint":fem_fingerprint(),
        "script_sha256":digest(__file__),"complete":False,"cases":[],"separate_refinements":[]}
    outputs={}
    for n,dt in ((640,30.),(1280,30.),(1280,15.),(1280,7.5),(640,7.5),(2560,7.5),(2560,3.75)):
        c,m=get_case(n,dt)
        outputs[(n,dt)]=(c,m)
        times,ii,jj=np.intersect1d(c["time_s"][:-1],reference["time_s"],return_indices=True)
        row={"label":m["label"],"elements":n,"max_step_s":dt,"event_time_s":m["event_time_s"],
            "event_delta_from_primary_s":m["event_time_s"]-rm["statistics"]["event_time_s"],
            "common_profile_times_s":times.tolist(),"sampled_field_max_difference":{
                f:float(np.max(np.abs(c[f][ii]-reference[f][jj]))) for f in ("temperature","moisture")}}
        report["cases"].append(row)
        save_json(ROOT/"reports/q3_independent_verification.json",report)
    for coarse,fine in (((1280,30.),(1280,15.)),((1280,15.),(1280,7.5)),
                        ((640,7.5),(1280,7.5)),((1280,7.5),(2560,7.5)),((2560,7.5),(2560,3.75))):
        a,am=outputs[coarse]; b,bm=outputs[fine]
        times,ii,jj=np.intersect1d(a["time_s"][:-1],b["time_s"][:-1],return_indices=True)
        report["separate_refinements"].append({"coarse":[*coarse],"fine":[*fine],
            "event_delta_s":bm["event_time_s"]-am["event_time_s"],
            "field_differences":{f:float(np.max(np.abs(a[f][ii]-b[f][jj]))) for f in ("temperature","moisture")}})
    report["complete"]=True
    report["scope"]="sampled_tail_profiles_and_global_event; independent_P1_BDF2; not_experimental_validation"
    save_json(ROOT/"reports/q3_independent_verification.json",report)
    print(json.dumps(report),flush=True)


if __name__=="__main__":
    main()
