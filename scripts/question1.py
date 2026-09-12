#!/usr/bin/env python3
"""Run Q1 only; cache full precision and per-run metadata, not a final claim."""
import argparse
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))

import numpy as np
from src.data_io import ROOT,load_config,verify_sources,save_json,source_fingerprint,digest
from src.boundary import Environment
from src.properties import Diffusivity
from src.fvm import RadialMesh,ScalarFVM
from src.solver import integrate_scalar
from src.reconstruction import reconstruct


def run_case(n=160,grading=1.0,rtol=1e-9,atol=1e-11,max_step=60,method="BDF",scheme="harmonic",label=None,store_cells=True,perturbation=None):
    cfg=load_config(); manifest=verify_sources(); env=Environment(cfg)
    # Q1 always starts from the original uniform state, never another question's cache.
    if cfg["problems"]["q1"] != {"properties":"q1","radius":"fixed","initialization":"fresh","end_s":1800}:
        raise ValueError("Q1 timeline/property contract changed; review before running")
    p=cfg["property_sets"]["q1"]; boundary=cfg["boundary"]; initial_cfg=cfg["initial"]
    if (p["rho_kg_m3"]["slope_C"] != 0 or p["cp_J_kg_K"]["fraction_C_over_1_plus_C"] != 0
            or p["k_W_m_K"]["fraction_C_over_1_plus_C"] != 0 or p["D_m2_s"]["thermal_exponent_K"] != 0):
        raise ValueError("The Q1 scalar solver requires constant thermal properties and D(C)")
    perturbation={} if perturbation is None else dict(perturbation)
    if set(perturbation)-{"h_multiplier","hm_multiplier","D_multiplier","beta_multiplier"}:
        raise ValueError("Unknown Q1 sensitivity parameter")
    if any(value<=0 for value in perturbation.values()):
        raise ValueError("Sensitivity multipliers must be positive")
    label=label or f"n{n}_g{grading:g}_{method}_{scheme}_r{rtol:g}_h{max_step:g}"
    folder=ROOT/"results/cache/q1"; folder.mkdir(parents=True,exist_ok=True)
    npz_path=folder/(label+".npz")
    meta_path=folder/(label+".json")
    times=np.arange(0.,1801.)
    radii=np.arange(21,dtype=float)*.001
    mesh=RadialMesh(n,cfg["geometry"]["reference_radius_m"],grading)
    meta={"label":label,"problem":"q1","property_set":"q1","start_s":0,"end_s":1800,
          "n":n,"grading":grading,"internal_flux":scheme,"rtol":rtol,"atol":atol,
          "max_step_s":max_step,"method":method,"source_fingerprint":source_fingerprint(),
          "source_inputs":manifest["files"],"script_sha256":digest(__file__),
          "perturbation":perturbation,
          "status":"computed_pending_comparative_validation","fields":{}}
    cache={"time_s":times,"radius_m":radii,"cell_centres_m":mesh.centres,"faces_m":mesh.faces,"weights":mesh.weights}
    specs=[("temperature",Diffusivity(p["k_W_m_K"]["offset"]),
            boundary["heat_transfer_W_m2_K"]*perturbation.get("h_multiplier",1.),
            p["rho_kg_m3"]["offset"]*p["cp_J_kg_K"]["offset"],initial_cfg["temperature_C"]),
           ("moisture",Diffusivity(p["D_m2_s"]["prefactor"]*perturbation.get("D_multiplier",1.),p["D_m2_s"]["moisture_exponent"]),
            boundary["normalized_mass_transfer_m_s"]*perturbation.get("hm_multiplier",1.),1.,initial_cfg["moisture_dry_basis"])]
    for field,law,transfer,capacity,initial in specs:
        ambient_scale=perturbation.get("beta_multiplier",1.) if field=="moisture" else 1.
        model=ScalarFVM(mesh,law,transfer,lambda t,f=field,b=ambient_scale:b*env.value(t,f),capacity,scheme)
        result,stats=integrate_scalar(model,initial,1800,times,env.times,method,rtol,atol,max_step)
        sampled=reconstruct(model,times,result[:,:n],radii,uniform_initial=initial)
        cache[field]=sampled
        cache[field+"_mean"]=np.sum(result[:,:n]*mesh.weights,axis=1)
        cache[field+"_outward_integral"]=result[:,-1]
        if store_cells:
            cache[field+"_cells"]=result[:,:n]
        meta["fields"][field]=stats
        print(label,field,"done",round(stats["elapsed_s"],3),"s; balance",stats["independent_balance_max_abs"],flush=True)
    np.savez_compressed(npz_path,**cache)
    meta["cache_sha256"]=digest(npz_path)
    save_json(meta_path,meta)
    return label,cache,meta


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--n",type=int,default=160)
    parser.add_argument("--grading",type=float,default=1.)
    parser.add_argument("--rtol",type=float,default=1e-9)
    parser.add_argument("--atol",type=float,default=1e-11)
    parser.add_argument("--max-step",type=float,default=60.)
    parser.add_argument("--method",choices=["BDF","Radau"],default="BDF")
    parser.add_argument("--scheme",choices=["harmonic","kirchhoff"],default="harmonic")
    parser.add_argument("--label")
    args=parser.parse_args()
    run_case(args.n,args.grading,args.rtol,args.atol,args.max_step,args.method,args.scheme,args.label)
