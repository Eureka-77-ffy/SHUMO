"""Run and cache a conditional geometry audit separately from final answers."""
import argparse
import json
import sys
import time
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import numpy as np
from scipy.integrate import solve_ivp
from numpy.polynomial.legendre import leggauss
from src.data_io import ROOT,load_config,verify_sources,digest,save_json
from src.solver import InitializedBDF
from scripts.question4 import fingerprint as primary_fingerprint
from validation.endface_model import AxisymmetricFVM


def fingerprint():
    d=primary_fingerprint()
    for p in (ROOT/"validation/endface_model.py",Path(__file__)):d[str(p.relative_to(ROOT))]=digest(p)
    return d


def run(group,nr,nz,end_multiplier=1.,rtol=1e-8,max_step=120.,label=None,limit=None):
    cfg=load_config();verify_sources();started=time.perf_counter()
    label=label or f"{group}_r{nr}_z{nz}_e{end_multiplier:g}_tol{rtol:g}_dt{max_step:g}"
    folder=ROOT/"results/cache/endface";folder.mkdir(parents=True,exist_ok=True)
    path=folder/(label+".npz");mp=path.with_suffix(".json")
    if path.exists() and mp.exists():
        meta=json.loads(mp.read_text())
        assert meta["fingerprint"]==fingerprint(),"Stale endface cache"
        assert meta["cache_sha256"]==digest(path)
        assert [meta[k] for k in ("group","nr","nz","end_multiplier","rtol","max_step_s","requested_limit")]==[group,nr,nz,end_multiplier,rtol,max_step,limit]
        return dict(np.load(path)),meta
    model=AxisymmetricFVM(cfg,group,nr,nz,end_multiplier)
    until=limit if limit is not None else (1800. if group=="q1" else 259200.)
    targets=np.unique(np.r_[0.,1.,10.,100.,300.,600.,900.,1200.,1500.,np.arange(1800.,until+1,1800.),until])
    cuts=np.unique(np.r_[targets[targets<=until],model.environment.times[model.environment.times<=until]])
    y=np.tile([28.,2.55],model.n);times=[0.];states=[y.copy()];samples=[model.sample(0.,y)]
    steps=nfev=nlu=0;balance=waterint=0.;event_time=None
    gx,gw=leggauss(3)
    def event(t,y):return model.maximum(t,y)-.15
    event.terminal=True;event.direction=-1
    for a,b in zip(cuts[:-1],cuts[1:]):
        model.environment.right_at_join=a>=14400.
        out=solve_ivp(model.rhs,(a,b),y,method=InitializedBDF,jac=model.jacobian,
                      rtol=rtol,atol=np.tile([rtol*.01,rtol*.0001],model.n),
                      max_step=max_step,dense_output=True,events=None if group=="q1" else event)
        if not out.success:raise RuntimeError(out.message)
        steps+=len(out.t)-1;nfev+=out.nfev;nlu+=out.nlu
        if np.min(out.y[1::2])<=0:raise RuntimeError("Nonpositive accepted 2D moisture")
        # Independent integration over the accepted trajectory of BOTH side
        # and end-face fluxes; no claim of an independent physics model.
        for j,dt in enumerate(np.diff(out.t)):
            tq=(out.t[j]+out.t[j+1])/2+dt*gx/2
            rate=[model.evaluate(float(t),out.sol(t))[1] for t in tq]
            waterint+=dt/2*float(gw@rate)
        y=out.y[:,-1];stop=float(out.t[-1])
        balance=max(balance,abs(float(model.weights@y[1::2])-2.55+waterint))
        if stop in targets or stop<b:
            times.append(stop);states.append(y.copy());samples.append(model.sample(stop,y))
        if stop<b or (group!="q1" and len(out.t_events[0])):
            event_time=float(out.t_events[0][0]);break
        if b%21600==0 or (group=="q1" and b==until):
            print(f"{label}: {b/3600:.2f} h, max C={model.maximum(b,y):.7f}, elapsed={time.perf_counter()-started:.1f}s",flush=True)
    cache={"time_s":np.array(times),"state":np.array(states),"xi_centres":model.rc,"z_centres_m":model.zc}
    for key in samples[0]:cache[key]=np.array([s[key] for s in samples])
    np.savez_compressed(path,**cache)
    meta={"label":label,"fingerprint":fingerprint(),"cache_sha256":digest(path),"group":group,"nr":nr,"nz":nz,
          "end_multiplier":end_multiplier,"rtol":rtol,"max_step_s":max_step,"requested_limit":limit,
          "event_time_s":event_time,"end_s":float(times[-1]),"steps":steps,"nfev":nfev,"nlu":nlu,
          "water_balance_max_at_cuts":balance,"water_balance_relative":balance/2.55,
          "initialization":"fresh","conditional_geometry_not_independent_primary_validation":True,
          "elapsed_s":time.perf_counter()-started}
    save_json(mp,meta);print(json.dumps(meta,ensure_ascii=False),flush=True)
    return cache,meta


if __name__=="__main__":
    p=argparse.ArgumentParser();p.add_argument("--group",choices=["q1","q23","q4"],required=True)
    p.add_argument("--nr",type=int,default=32);p.add_argument("--nz",type=int,default=32)
    p.add_argument("--end-multiplier",type=float,default=1.);p.add_argument("--rtol",type=float,default=1e-8)
    p.add_argument("--max-step",type=float,default=120.);p.add_argument("--label");p.add_argument("--limit",type=float)
    run(**vars(p.parse_args()))
