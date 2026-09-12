#!/usr/bin/env python3
"""Check B/C counterfactual event precision without reusing the primary RHS.

These are imposed-geometry scenarios, not observed alternate experiments.
The independent material weak-form solver starts from its own uniform state.
"""
import json
import sys
import time
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))

import numpy as np
from scipy.optimize import brentq

from scripts.question4 import fingerprint, read_case
from src.data_io import ROOT,digest,load_config,save_json,verify_sources
from validation.q4_convergence import get_case,compare
from validation.q4_independent_fem import MaterialFEM,fingerprint as fem_fingerprint


def independent_fingerprint():
    return {**fem_fingerprint(),str(Path(__file__).relative_to(ROOT)):digest(__file__)}


def independent_run(letter,n):
    label=f"factorial_{letter}_fem_n{n}"
    folder=ROOT/"results/cache/q4_factorial_independent"
    path=folder/(label+".npz");mp=folder/(label+".json")
    try:
        meta=json.loads(mp.read_text())
        assert meta["source_fingerprint"]==independent_fingerprint() and meta["cache_sha256"]==digest(path)
        with np.load(path) as raw:return {k:raw[k] for k in raw.files},meta
    except (FileNotFoundError,AssertionError):pass
    started=time.perf_counter();cfg=load_config();model=MaterialFEM(n,cfg)
    if letter=="B":
        model.p=cfg["property_sets"]["q23"];limit=259200.
    elif letter=="C":
        model.radius_at=lambda t:.02
        limit=604800.
    else:raise ValueError("Only B/C need additional counterfactual tests")
    T=np.full(n+1,28.);C=np.full(n+1,2.55)
    t=0.;older=previous=None;steps=0
    sampled=np.unique(np.r_[0.,1.,10.,60.,600.,1800.,np.arange(3600.,limit+1,1800.)])
    cuts=np.unique(np.r_[sampled,model.original_air[:,0],model.radius_observations[:,0],limit])
    times,Ts,Cs=[0.],[T.copy()],[C.copy()]
    event_time=None
    for stop in cuts[1:]:
        if t>=14400.:model.air=np.array([[14400.,*model.platform],[limit,*model.platform]])
        while t<stop-1e-8:
            # Resolve the driven thermal transient independently; after 5 h,
            # a separate late cap controls the much slower moisture trajectory.
            cap=1.875 if t<18000. else 15.
            dt=min(cap,.04*(1.875/30)*max(t,.01),stop-t)
            dt=min(dt,.001) if previous is None else min(dt,1.8*previous)
            newT,newC,_=model.step(t+dt,dt,T,C,older,previous)
            steps+=1
            if newC.max()<=.15:
                def fun(h):
                    if h==0:return float(C.max()-.15)
                    return float(model.step(t+h,h,T,C,older,previous)[1].max()-.15)
                tau=brentq(fun,0.,dt,xtol=1e-6)
                event_time=t+tau
                T,C,_=model.step(event_time,tau,T,C,older,previous)
                t=event_time;break
            older,previous=(T,C),dt;T,C=newT,newC;t+=dt
        if event_time is not None:
            times.append(t);Ts.append(T.copy());Cs.append(C.copy());break
        t=float(stop)
        if stop in sampled:times.append(t);Ts.append(T.copy());Cs.append(C.copy())
        if stop==14400.:older=previous=None
        if stop%86400==0:print(f"Factorial {letter} independent n={n}: {t/3600:g} h, Cmax={C.max():.7f}",flush=True)
    if event_time is None:raise RuntimeError("Counterfactual event outside stated horizon")
    rr=np.arange(21)*.001
    times=np.asarray(times);radii=np.array([model.radius_at(t) for t in times])
    fields=[]
    for rows in (Ts,Cs):
        out=np.full((len(times),22),np.nan)
        for i,(r,row) in enumerate(zip(radii,rows)):
            valid=rr<=r;out[i,np.flatnonzero(valid)]=np.interp(rr[valid]/r,model.nodes/.02,row)
            out[i,-1]=row[-1]
        fields.append(out)
    cache={"time_s":times,"radius_current_m":radii,"temperature":fields[0],"moisture":fields[1]}
    folder.mkdir(parents=True,exist_ok=True);np.savez_compressed(path,**cache)
    meta={"label":label,"letter":letter,"n":n,"source_fingerprint":independent_fingerprint(),
        "cache_sha256":digest(path),"event_time_s":event_time,"steps":steps,
        "early_max_step_s":1.875,"late_max_step_s":15.,"initialization":"fresh",
        "elapsed_s":time.perf_counter()-started}
    save_json(mp,meta);print(meta,flush=True)
    return cache,meta


def main():
    verify_sources()
    report={"source_fingerprint":fingerprint(),"script_sha256":digest(__file__),
        "independent_fingerprint":independent_fingerprint(),"cases":[],"complete":False}
    for letter,group,geometry,limit in (("B","q23","linear",259200.),("C","q4","fixed",604800.)):
        primary,pm=read_case(f"factorial_{letter}")
        cap,cm=get_case(f"factorial_{letter}_cap15",n=1280,group=group,geometry=geometry,
                       limit=limit,max_step=15.,rtol=1e-10)
        row={"letter":letter,"primary_label":pm["label"],"cap15_label":cm["label"],
            "primary_time_comparison":compare(cap,cm,primary,pm),"independent":[]}
        for n in (1280,2560):
            c,m=independent_run(letter,n)
            times,i,j=np.intersect1d(c["time_s"][:-1],primary["time_s"],return_indices=True)
            fields={f:float(np.nanmax(np.abs(c[f][i]-primary[f][j]))) for f in ("temperature","moisture")}
            row["independent"].append({"label":m["label"],"n":n,"event_time_s":m["event_time_s"],
                "event_delta_primary_s":m["event_time_s"]-pm["statistics"]["event_time_s"],
                "sampled_field_max_differences":fields,"common_times":len(times)})
        row["independent_grid_event_change_s"]=row["independent"][1]["event_time_s"]-row["independent"][0]["event_time_s"]
        report["cases"].append(row);save_json(ROOT/"reports/q4_factorial_verification.json",report)
    report["complete"]=True
    save_json(ROOT/"reports/q4_factorial_verification.json",report)
    print(json.dumps(report),flush=True)


if __name__=="__main__":main()
