#!/usr/bin/env python3
"""Independent moving P1 weak form: current mass, fixed reference stiffness.

On the material mesh r=R*xi, the time term is R^2 M(b) udot, the
stiffness is K(a), and Robin is R*h. No derivative of R^2*u is used.
Does not import primary moving geometry, RHS, surface or event code.
"""
import argparse
import json
import sys
import time
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import numpy as np
from openpyxl import load_workbook
from scipy.optimize import brentq

from src.data_io import ROOT,digest,load_config,save_json,verify_sources
from validation.q2_independent_fem import IndependentFEM


class MaterialFEM(IndependentFEM):
    def __init__(self,n,cfg):
        super().__init__(n,cfg)
        self.p=cfg["property_sets"]["q4"]
        self.R0=self.R
        wb=load_workbook(Path(cfg["source_directory"])/"附件2.xlsx",read_only=True,data_only=True)
        self.radius_observations=np.asarray(list(wb.active.values)[1:],float);wb.close()
        self.radius_observations[:,1]*=.01
        self.original_air=self.air.copy()
        self.platform=self.air[self.air[:,0]>=10800.,1:].mean(axis=0)

    def radius_at(self,t):
        if t<0 or t>259200.:raise ValueError("Independent radius outside observed interval")
        return float(np.interp(t,self.radius_observations[:,0],self.radius_observations[:,1]))

    def step(self,t,dt,oldT,oldC,older=None,previous_dt=None):
        self.R=self.radius_at(t)
        scale=(self.R/self.R0)**2
        ta=np.interp(t,self.air[:,0],self.air[:,1])
        ce=self.config["boundary"]["beta"]*np.interp(t,self.air[:,0],self.air[:,2])
        T,C=oldT.copy(),oldC.copy()
        effectiveT,effectiveC,effective_dt=oldT,oldC,dt
        if older is not None:
            ratio=dt/previous_dt;a0=(1+2*ratio)/(1+ratio)
            effectiveT=((1+ratio)*oldT-ratio**2/(1+ratio)*older[0])/a0
            effectiveC=((1+ratio)*oldC-ratio**2/(1+ratio)*older[1])/a0
            effective_dt=dt/a0
        for it in range(60):
            Cg,Tg=self.gauss_state(C),self.gauss_state(T)
            b,k,_=self.coefficients(Cg,Tg)
            # Geometry multiplies CURRENT mass times the material derivative;
            # multiplying old state by the old volume would create spurious dilution.
            newT=self.linear_step(self.mass_matrix(scale*b),self.stiffness(k),effectiveT,effective_dt,self.h,ta)
            D=self.coefficients(Cg,self.gauss_state(newT))[2]
            newC=self.linear_step(tuple(scale*m for m in self.mass),self.stiffness(D),effectiveC,effective_dt,self.hm,ce)
            if np.min(newC)<=0 or not np.isfinite(newT).all():raise RuntimeError("Nonphysical independent iterate")
            errorT,errorC=np.max(np.abs(newT-T)),np.max(np.abs(newC-C))
            T,C=newT,newC
            if errorT<2e-10 and errorC<2e-12:return T,C,it+1
        raise RuntimeError(f"Independent moving Picard failed at {t}, dt={dt}, {errorT}, {errorC}")


def fingerprint():
    paths=[Path(__file__),ROOT/"validation/q2_independent_fem.py",ROOT/"config/model_config.json"]
    return {str(p.relative_to(ROOT)):digest(p) for p in paths}


def run(n=640,max_step=15.):
    started=time.perf_counter();verify_sources();cfg=load_config();model=MaterialFEM(n,cfg)
    T=np.full(n+1,28.);C=np.full(n+1,2.55)
    t=0.;older=previous_dt=None;steps=iterations=0
    times,Ts,Cs=[0.],[T.copy()],[C.copy()]
    samples=np.unique(np.r_[0.,1.,10.,60.,100.,300.,600.,900.,1200.,1800.,np.arange(3600.,259201.,1800.)])
    cuts=np.unique(np.r_[model.original_air[:,0],model.radius_observations[:,0],samples])
    threshold=.15;event_time=None;actual_max=0.
    for stop in cuts[1:]:
        if t>=14400.:model.air=np.array([[14400.,*model.platform],[259200.,*model.platform]])
        while t<stop-1e-8:
            dt=min(max_step,.04*(max_step/30)*max(t,.01),stop-t)
            if previous_dt is None:dt=min(dt,.001)
            else:dt=min(dt,1.8*previous_dt)
            newT,newC,count=model.step(t+dt,dt,T,C,older,previous_dt)
            iterations+=count;steps+=1;actual_max=max(actual_max,dt)
            if newC.max()<=threshold:
                def residual(h):
                    if h==0:return float(C.max()-threshold)
                    return float(model.step(t+h,h,T,C,older,previous_dt)[1].max()-threshold)
                tau=brentq(residual,0.,dt,xtol=1e-6)
                event_time=t+tau
                T,C,_=model.step(event_time,tau,T,C,older,previous_dt)
                t=event_time
                break
            older,previous_dt=(T,C),dt;T,C=newT,newC;t+=dt
        if event_time is not None:
            times.append(t);Ts.append(T.copy());Cs.append(C.copy());break
        t=float(stop)
        if stop in samples:times.append(t);Ts.append(T.copy());Cs.append(C.copy())
        if stop==14400.:older=previous_dt=None
        if stop%43200==0:print(f"Q4 FEM n={n}, dt={max_step}, {t/3600:.0f} h, Cmax={C.max():.7f}",flush=True)
    if event_time is None:raise RuntimeError("Independent Q4 did not cross within observations")
    radii=np.arange(21)*.001;times=np.asarray(times)
    current=np.array([model.radius_at(t) for t in times])
    fields=[]
    for rows in (Ts,Cs):
        array=np.full((len(times),22),np.nan)
        for i,(row,R) in enumerate(zip(rows,current)):
            valid=radii<=R
            array[i,np.flatnonzero(valid)]=np.interp(radii[valid]/R,model.nodes/model.R0,row)
            array[i,-1]=row[-1]
        fields.append(array)
    label=f"q4_fem_n{n}_dt{max_step:g}"
    folder=ROOT/"results/cache/q4_independent";folder.mkdir(parents=True,exist_ok=True)
    path=folder/(label+".npz")
    cache={"time_s":times,"radius_current_m":current,"radius_fixed_m":radii,
        "nodes_xi":model.nodes/model.R0,"temperature_nodes":np.asarray(Ts),"moisture_nodes":np.asarray(Cs),
        "temperature":fields[0],"moisture":fields[1]}
    np.savez_compressed(path,**cache)
    meta={"label":label,"n":n,"max_step_s":max_step,"event_time_s":event_time,
        "event_radius_m":model.radius_at(event_time),"critical_maximum":float(C.max()),
        "maximum_xi":float(model.nodes[np.argmax(C)]/model.R0),"accepted_max_step_s":actual_max,
        "steps":steps,"picard_iterations":iterations,"initialization":"fresh","time_origin_s":0.,
        "source_fingerprint":fingerprint(),"cache_sha256":digest(path),"elapsed_s":time.perf_counter()-started,
        "method":"independent_material_P1_consistent_current_mass_variable_BDF2_and_substep_event"}
    save_json(folder/(label+".json"),meta);print(json.dumps(meta),flush=True)
    return cache,meta


if __name__=="__main__":
    p=argparse.ArgumentParser();p.add_argument("--n",type=int,default=640);p.add_argument("--max-step",type=float,default=15.)
    run(**vars(p.parse_args()))
