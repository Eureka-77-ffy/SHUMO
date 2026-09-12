#!/usr/bin/env python3
"""Independent moving FE space/time refinements and material-mass checks."""
import copy
import json
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import numpy as np
from scripts.question4 import fingerprint,read_case
from src.data_io import ROOT,digest,load_config,save_json,verify_sources
from validation.q2_independent_fem import IndependentFEM
from validation.q4_independent_fem import run,fingerprint as independent_fingerprint,MaterialFEM


def get_case(n,dt):
    label=f"q4_fem_n{n}_dt{dt:g}"
    folder=ROOT/"results/cache/q4_independent"
    try:
        m=json.loads((folder/(label+".json")).read_text());path=folder/(label+".npz")
        assert m["source_fingerprint"]==independent_fingerprint() and m["cache_sha256"]==digest(path)
        with np.load(path) as raw:return {k:raw[k] for k in raw.files},m
    except (FileNotFoundError,AssertionError):return run(n,dt)


def main():
    verify_sources();cfg=load_config();primary,pm=read_case("q4_final")
    m=MaterialFEM(20,cfg);fixed=IndependentFEM(20,cfg);fixed.p=cfg["property_sets"]["q4"]
    m.radius_at=lambda t:.02
    T=np.linspace(31.,33.,21);C=np.linspace(2.,1.,21)
    Ta,Ca,_=m.step(100.,2.,T,C);Tb,Cb,_=fixed.step(100.,2.,T,C)
    assert np.array_equal(Ta,Tb) and np.array_equal(Ca,Cb)
    m=MaterialFEM(20,cfg);m.p=copy.deepcopy(m.p);m.p["D_m2_s"]["prefactor"]=0.;m.h=m.hm=0.
    T=np.full(21,31.);C=1+.2*np.cos(np.pi*m.nodes/.02)
    newT,newC,_=m.step(259200.,1800.,T,C)
    assert np.max(np.abs(newC-C))<1e-12 and np.max(np.abs(newT-T))<1e-10
    report={"source_fingerprint":fingerprint(),"script_sha256":digest(__file__),
        "independent_fingerprint":independent_fingerprint(),"complete":False,
        "independent_fixed_radius_step_exact_match":True,
        "independent_pure_shrink_moisture_change":float(np.max(np.abs(newC-C))),
        "cases":[],"separate_refinements":[]}
    outputs={}
    cases=((320,15.),(640,15.),(1280,15.),(1280,7.5),(1280,3.75),(2560,3.75),(2560,1.875))
    for n,dt in cases:
        c,meta=get_case(n,dt);outputs[(n,dt)]=(c,meta)
        times,i,j=np.intersect1d(c["time_s"][:-1],primary["time_s"],return_indices=True)
        fields={}
        for f in ("temperature","moisture"):
            a,b=c[f][i],primary[f][j]
            assert np.array_equal(np.isnan(a),np.isnan(b))
            d=np.abs(a-b);loc=np.unravel_index(np.nanargmax(d),d.shape)
            fields[f]={"max_abs":float(d[loc]),"time_s":float(times[loc[0]]),
                "position":"surface" if loc[1]==21 else f"{loc[1]/10:g} cm"}
        report["cases"].append({"label":meta["label"],"n":n,"max_step_s":dt,
            "event_time_s":meta["event_time_s"],"event_delta_primary_s":meta["event_time_s"]-pm["statistics"]["event_time_s"],
            "sampled_comparisons":fields,"common_times_s":times.tolist()})
        save_json(ROOT/"reports/q4_independent_verification.json",report)
    for coarse,fine in (((320,15.),(640,15.)),((640,15.),(1280,15.)),
         ((1280,15.),(1280,7.5)),((1280,7.5),(1280,3.75)),((1280,3.75),(2560,3.75)),((2560,3.75),(2560,1.875))):
        a,am=outputs[coarse];b,bm=outputs[fine]
        times,i,j=np.intersect1d(a["time_s"][:-1],b["time_s"][:-1],return_indices=True)
        report["separate_refinements"].append({"coarse":list(coarse),"fine":list(fine),
            "event_delta_s":bm["event_time_s"]-am["event_time_s"],
            "field_differences":{f:float(np.nanmax(np.abs(a[f][i]-b[f][j]))) for f in ("temperature","moisture")}})
    report["complete"]=True
    save_json(ROOT/"reports/q4_independent_verification.json",report);print(json.dumps(report),flush=True)


if __name__=="__main__":main()
