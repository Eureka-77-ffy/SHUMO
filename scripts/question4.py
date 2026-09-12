#!/usr/bin/env python3
"""Q4 and controlled radius/property combinations, always from original t=0."""
import argparse
import json
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import numpy as np

from scripts.question3 import fingerprint as previous_fingerprint
from src.data_io import ROOT,digest,load_config,save_json,verify_sources
from src.events.integration import ScenarioEnvironment
from src.fvm import RadialMesh
from src.moving.model import MovingCoupled,Radius
from src.moving.integration import integrate


def fingerprint():
    fp=previous_fingerprint()
    fp.update({str(p.relative_to(ROOT)):digest(p) for p in sorted((ROOT/"src/moving").glob("*.py"))+[Path(__file__)]})
    return fp


def read_case(label):
    folder=ROOT/"results/cache/q4"; path=folder/(label+".npz")
    m=json.loads((folder/(label+".json")).read_text())
    if m["source_fingerprint"]!=fingerprint() or m["cache_sha256"]!=digest(path):
        raise ValueError("Stale Q4 cache: "+label)
    with np.load(path) as raw: return {k:raw[k] for k in raw.files},m


def run_case(label="q4_trial",n=320,group="q4",geometry="linear",extension="error",
             scenario="main",max_step=30.,rtol=1e-10,method="BDF",limit=259200.,perturbation=None):
    cfg=load_config(); manifest=verify_sources()
    expected={"properties":"q4","radius":"observed_shrinkage","initialization":"fresh","end":"global_moisture_event_or_radius_data_end"}
    if cfg["problems"]["q4"]!=expected: raise ValueError("Q4 timeline contract changed")
    radius=Radius(cfg,geometry,extension)
    model=MovingCoupled(RadialMesh(n,.02,2.),cfg,ScenarioEnvironment(cfg,scenario),radius,group,perturbation)
    c,stats=integrate(model,limit=limit,max_step=max_step,rtol=rtol,method=method)
    assert np.array_equal(c["profile_state"][0],np.tile([28.,2.55],n))
    folder=ROOT/"results/cache/q4"; folder.mkdir(parents=True,exist_ok=True)
    path=folder/(label+".npz"); np.savez_compressed(path,**c)
    m={"label":label,"n":n,"grading":2.,"property_set":group,"geometry":geometry,
        "radius_extension":extension,"boundary_scenario":scenario,"initialization":"fresh",
        "physical_time_origin_s":0.,"rtol":rtol,"atol_temperature":1e-11,"atol_moisture":1e-13,
        "max_step_s":max_step,"method":method,"requested_limit_s":limit,"output_step_s":60.,
        "perturbation":dict(perturbation or {}),"source_fingerprint":fingerprint(),
        "cache_sha256":digest(path),"source_inputs":manifest["files"],"statistics":stats,
        "status":"computed_pending_Q4_numerical_and_independent_validation"}
    save_json(folder/(label+".json"),m)
    print(json.dumps({"label":label,"statistics":stats}),flush=True)
    return c,m


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--label",default="q4_trial");p.add_argument("--n",type=int,default=320)
    p.add_argument("--group",choices=["q4","q23"],default="q4")
    p.add_argument("--geometry",choices=["linear","pchip","fixed"],default="linear")
    p.add_argument("--extension",choices=["error","hold_last"],default="error")
    p.add_argument("--scenario",choices=["main","last_observation","mean_9000_14400"],default="main")
    p.add_argument("--max-step",type=float,default=30.);p.add_argument("--rtol",type=float,default=1e-10)
    p.add_argument("--method",choices=["BDF","Radau"],default="BDF")
    p.add_argument("--limit",type=float,default=259200.)
    run_case(**vars(p.parse_args()))
