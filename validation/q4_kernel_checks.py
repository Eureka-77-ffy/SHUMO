#!/usr/bin/env python3
"""Check geometry, Jacobians, pure material shrinkage and manufactured PDEs."""
import copy
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import numpy as np
from scipy.integrate import solve_ivp

from scripts.question4 import fingerprint
from src.coupled.model import CoupledFVM
from src.data_io import ROOT,digest,load_config,save_json,verify_sources
from src.events.integration import ScenarioEnvironment
from src.fvm import RadialMesh
from src.moving.model import MovingCoupled,Radius
from src.moving.integration import reconstruct


class SmoothRadius:
    kind="manufactured"
    extension="error"
    times=np.array([0.,600.])
    def value(self,t): return .02*(1-.2*np.asarray(t)/600.)


class ManufacturedEnvironment:
    times=np.array([0.,600.])
    def __init__(self,cfg):
        self.radius=SmoothRadius()
        self.p=cfg["property_sets"]["q4"]
    def exact(self,t,xi):
        T=30+2*np.sin(t/600)+3*(1+.1*t/600)*xi**2
        C=1.6-.1*t/600+.2*xi**2
        return T,C
    def value(self,t,field):
        T,C=self.exact(np.asarray(t),1.)
        R=self.radius.value(t)
        if field=="temperature":
            k=.12+.20*C/(1+C)
            return T+k/R*(6*(1+.1*np.asarray(t)/600))/25.
        D=.00042*np.exp(-.30/C-3850/(T+273.15))
        return C+D/R*.4/8e-7


def manufactured(n):
    cfg=load_config(); env=ManufacturedEnvironment(cfg)
    model=MovingCoupled(RadialMesh(n,.02,2.),cfg,env,env.radius)
    xi=model.mesh.centres/.02
    T,C=env.exact(0.,xi)
    initial=np.empty(2*n);initial[::2]=T;initial[1::2]=C
    def forcing(t):
        T,C=env.exact(t,xi)
        b,k,_,kp=model.thermal(C)
        D,DC,DT=model.diffusion(C,T)
        B=.2; E=3*(1+.1*t/600)
        LT=(4*E*k+4*E*B*xi**2*kp)/env.radius.value(t)**2/b
        LC=(4*B*D+4*B*xi**2*(DC*B+DT*E))/env.radius.value(t)**2
        out=np.empty(2*n)
        out[::2]=2*np.cos(t/600)/600+.3/600*xi**2-LT
        out[1::2]=-.1/600-LC
        return out
    result=solve_ivp(lambda t,y:model.rhs(t,y)+forcing(t),(0.,600.),initial,
        method="Radau",rtol=1e-11,atol=1e-13,jac=model.jacobian,max_step=2.)
    assert result.success
    T,C=env.exact(600.,xi)
    return {"n":n,"temperature_max_error":float(np.max(np.abs(result.y[::2,-1]-T))),
            "moisture_max_error":float(np.max(np.abs(result.y[1::2,-1]-C)))}


def main():
    verify_sources();cfg=load_config();env=ScenarioEnvironment(cfg)
    report={"source_fingerprint":fingerprint(),"script_sha256":digest(__file__),"jacobians":[]}
    radius=Radius(cfg)
    assert np.allclose(radius.value(radius.times),radius.radii,rtol=0,atol=1e-17)
    rejected=False
    try: radius.value(259201.)
    except ValueError: rejected=True
    assert rejected
    report["radius_extrapolation_rejected"]=True
    for group in ("q23","q4"):
        mesh=RadialMesh(9,.02,2.)
        moving=MovingCoupled(mesh,cfg,env,Radius(cfg,"fixed"),group)
        fixed=CoupledFVM(mesh,cfg,env,group)
        y=np.empty(18);y[::2]=np.linspace(32,38,9);y[1::2]=np.linspace(2.5,.8,9)
        assert np.array_equal(moving.rhs(1000.,y),fixed.rhs(1000.,y))
        assert np.array_equal(moving.jacobian(1000.,y).toarray(),fixed.jacobian(1000.,y).toarray())
        moving.radius=radius
        for t in (100.,10000.,100000.):
            numeric=np.column_stack([(moving.rhs(t,y+1e-5*v)-moving.rhs(t,y-1e-5*v))/2e-5 for v in np.eye(18)])
            analytic=moving.jacobian(t,y).toarray()
            error=float(np.max(np.abs(numeric-analytic))/np.max(np.abs(analytic)))
            assert error<1e-7
            rhs=moving.rhs(t,y);ts,cs=moving.surface(t,y[-2],y[-1]);ta,ce=moving.ambient(t)
            water=np.sum(mesh.weights*rhs[1::2])+2*moving.physical_hm/radius.value(t)*(cs-ce)
            heat=np.sum(mesh.weights*moving.thermal(y[1::2])[0]*rhs[::2])+2*moving.physical_h/radius.value(t)*(ts-ta)
            assert abs(water)<1e-14 and abs(heat)<1e-7
            # A second calculation on the current physical mesh must agree;
            # its distances, face areas and volumes all change explicitly.
            physical=CoupledFVM(RadialMesh(9,float(radius.value(t)),2.),cfg,env,group)
            physical_delta=float(np.max(np.abs(physical.rhs(t,y)-rhs)))
            assert physical_delta<1e-10
            report["jacobians"].append({"group":group,"time_s":t,"relative_error":error,
                "physical_reference_operator_max_difference":physical_delta})
    report["fixed_radius_degenerates_exactly_to_same_property_fixed_model"]=True
    pure_cfg=copy.deepcopy(cfg)
    pure_cfg["property_sets"]["q4"]["D_m2_s"]["prefactor"]=0.
    pure_cfg["boundary"]["heat_transfer_W_m2_K"]=0.
    pure_cfg["boundary"]["normalized_mass_transfer_m_s"]=0.
    model=MovingCoupled(RadialMesh(40,.02,2.),pure_cfg,env,radius)
    xi=model.mesh.centres/.02
    state=np.empty(80);state[::2]=31.;state[1::2]=1+.2*np.cos(np.pi*xi)
    result=solve_ivp(model.rhs,(0.,259200.),state,method="Radau",jac=model.jacobian,
                      rtol=1e-10,atol=1e-12,max_step=1800.)
    assert result.success
    delta=float(np.max(np.abs(result.y-state[:,None])))
    assert delta<1e-11
    R=radius.value(result.t)
    dry_measure=(.02/R)**2*R**2/.02**2
    water_mass=np.sum(model.mesh.weights[:,None]*result.y[1::2],axis=0)*dry_measure
    assert np.max(np.abs(dry_measure-1))<1e-14 and np.ptp(water_mass)<1e-12
    report["pure_D0_hm0_material_shrinkage"]={"profile_max_change":delta,
        "relative_dry_mass_max_change":float(np.max(np.abs(dry_measure-1))),
        "normalized_water_mass_range":float(np.ptp(water_mass)),
        "incorrect_volume_concentration_compression_would_change_C":True}
    # Vectorized-time geometry and actual-domain blanks, not stale scalar R.
    normal=MovingCoupled(RadialMesh(40,.02,2.),cfg,env,radius)
    yy=np.tile([28.,2.55],40)
    tt=np.array([0.,1800.,21600.,259200.])
    Y=np.tile(yy[:,None],(1,len(tt)))
    batched=normal.rhs(tt,Y)
    separate=np.column_stack([normal.rhs(t,yy) for t in tt])
    assert np.allclose(batched,separate,rtol=1e-13,atol=1e-13)
    values,extra=reconstruct(normal,tt,Y)
    valid=np.arange(21)[None,:]*.001<=radius.value(tt)[:,None]
    assert np.array_equal(np.isfinite(values[1][:,:21]),valid)
    assert np.isfinite(values[1][:,-1]).all()
    report["batched_geometry_and_outside_domain_blank_mask"]=True
    report["manufactured_solution"]=[manufactured(n) for n in (40,80,160,320)]
    for a,b in zip(report["manufactured_solution"][:-1],report["manufactured_solution"][1:]):
        assert b["temperature_max_error"]<a["temperature_max_error"]*.4
        assert b["moisture_max_error"]<a["moisture_max_error"]*.4
    report["manufactured_problem"]="smooth R(t), variable appendix4 b/k/D, prescribed nonuniform T,C, analytical sources and physical Robin data"
    report["status"]="passed_Q4_material_coordinate_kernel_and_manufactured_checks"
    save_json(ROOT/"reports/q4_kernel_checks.json",report)
    print(report,flush=True)


if __name__=="__main__":main()
