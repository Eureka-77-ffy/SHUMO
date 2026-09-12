"""Complete Q1's pending independent nonlinear P1 water-field check."""
import json
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import numpy as np
from src.data_io import ROOT,load_config,digest,save_json,verify_sources
from validation.q2_independent_fem import IndependentFEM
from src.fvm import RadialMesh,ScalarFVM
from src.properties import Diffusivity
from src.boundary import Environment
from src.reconstruction import reconstruct


def fingerprint():
    return {str(p.relative_to(ROOT)):digest(p) for p in (Path(__file__),ROOT/"validation/q2_independent_fem.py",ROOT/"config/model_config.json")}


def run(n,dt):
    label=f"q1_fem_n{n}_dt{dt:g}";folder=ROOT/"results/cache/q1_independent";folder.mkdir(parents=True,exist_ok=True)
    path=folder/(label+".npz");mp=path.with_suffix(".json")
    if path.exists() and mp.exists():
        meta=json.loads(mp.read_text());assert meta["fingerprint"]==fingerprint() and meta["cache_sha256"]==digest(path)
        return dict(np.load(path)),meta
    cfg=load_config();m=IndependentFEM(n,cfg);m.p=cfg["property_sets"]["q1"]
    targets=np.array([0.,1.,10.,100.,300.,600.,900.,1200.,1500.,1800.])
    cuts=np.unique(np.r_[targets,m.air[m.air[:,0]<=1800,0]])
    T=np.full(n+1,28.);C=np.full(n+1,2.55);ts=[0.];Ts=[T.copy()];Cs=[C.copy()]
    t=0.;older=previous=None;steps=0
    for stop in cuts[1:]:
        while t<stop-1e-10:
            h=min(dt,.04*(dt/30)*max(t,.01),stop-t)
            if previous is not None:h=min(h,1.8*previous)
            newT,newC,_=m.step(t+h,h,T,C,older,previous)
            older,previous=(T,C),h;T,C=newT,newC;t+=h;steps+=1
        t=float(stop)
        if stop in targets:ts.append(t);Ts.append(T.copy());Cs.append(C.copy())
    rr=np.arange(21)*.001
    cache={"time_s":np.array(ts),"nodes_m":m.nodes,"temperature_nodes":np.array(Ts),"moisture_nodes":np.array(Cs),
        "temperature":np.array([np.interp(rr,m.nodes,row) for row in Ts]),"moisture":np.array([np.interp(rr,m.nodes,row) for row in Cs])}
    np.savez_compressed(path,**cache)
    meta={"label":label,"elements":n,"max_step_s":dt,"steps":steps,"fingerprint":fingerprint(),"cache_sha256":digest(path),"group":"q1","initialization":"fresh"}
    save_json(mp,meta);print(label,"completed",steps,"steps",flush=True);return cache,meta


def main():
    verify_sources();cfg=load_config();primary_path=ROOT/"results/cache/q1/q1_final.npz"
    pm=json.loads(primary_path.with_suffix(".json").read_text());assert digest(primary_path)==pm["cache_sha256"]
    primary=dict(np.load(primary_path));outputs={};cases=[]
    for n,dt in ((320,1.875),(640,1.875),(1280,1.875),(1280,.9375),(1280,.46875),(2560,.46875)):
        a,meta=run(n,dt);outputs[(n,dt)]=a;indices=a["time_s"].astype(int)
        comparisons={f:float(np.max(abs(a[f]-primary[f][indices]))) for f in ("temperature","moisture")}
        cases.append({**meta,"sampled_max_differences":comparisons})
    final=outputs[(2560,.46875)];ts=final["time_s"];mesh=RadialMesh(1280,.02,2);env=Environment(cfg)
    p=cfg["property_sets"]["q1"];dense={}
    for f,law,transfer,capacity,initial in (("temperature",Diffusivity(.36),25.,820*2600,28.),("moisture",Diffusivity(7e-9,.89),8e-7,1.,2.55)):
        model=ScalarFVM(mesh,law,transfer,lambda t,field=f:env.value(t,field),capacity)
        values=reconstruct(model,ts,primary[f+"_cells"][ts.astype(int)],final["nodes_m"],uniform_initial=initial)
        dense[f]=float(np.max(abs(values-final[f+"_nodes"])))
    refinements={}
    for key,a,b in (("space_640_1280_dt1p875",(640,1.875),(1280,1.875)),("time_1p875_0p9375",(1280,1.875),(1280,.9375)),
                     ("time_0p9375_0p46875",(1280,.9375),(1280,.46875)),("space_1280_2560_dt0p46875",(1280,.46875),(2560,.46875))):
        refinements[key]={f:float(np.max(abs(outputs[a][f]-outputs[b][f]))) for f in ("temperature","moisture")}
    assert max(dense.values())<5e-5 and max(cases[-1]["sampled_max_differences"].values())<5e-5
    report={"status":"passed_Q1_original_nonlinear_independent_sample_verification","fingerprint":fingerprint(),"primary_cache_sha256":digest(primary_path),
            "time_s":ts.tolist(),"cases":cases,"separate_refinements":refinements,"dense_profiles":{"radial_nodes":2561,"max_abs":dense},
            "scope":"Q1 0-1800s; independent nodal weak form and own BDF2; finite samples, not rigorous error bound","complete":True}
    save_json(ROOT/"reports/q1_nonlinear_independent.json",report);print(json.dumps({"dense":dense,"refinements":refinements}),flush=True)


if __name__=="__main__":main()
