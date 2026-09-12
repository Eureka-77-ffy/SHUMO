#!/usr/bin/env python3
"""Q4 material-coordinate, dense-profile and error-budget diagnostics."""
import json
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import numpy as np
from scipy.interpolate import PchipInterpolator

from scripts.question4 import read_case,fingerprint
from src.data_io import ROOT,digest,load_config,save_json,verify_sources
from src.events.integration import ScenarioEnvironment,global_maximum
from src.fvm import RadialMesh
from src.moving.model import MovingCoupled,Radius
from validation.q4_independent_fem import fingerprint as independent_fingerprint


def reference_nodes(model,times,states):
    cells=[states[:,::2],states[:,1::2]]
    surfaces=model.surface(times,cells[0][:,-1],cells[1][:,-1])
    x=model.mesh.centres/.02;x0,x1=x[:2]
    out=[]
    for cell,surface,initial in zip(cells,surfaces,(28.,2.55)):
        centre=(x1*x1*cell[:,0]-x0*x0*cell[:,1])/(x1*x1-x0*x0)
        values=np.column_stack((centre,cell,surface));values[times==0]=initial
        out.append(values)
    return np.r_[0.,x,1.],out


def main():
    verify_sources();cfg=load_config();cache,meta=read_case("q4_final")
    env=ScenarioEnvironment(cfg);rad=Radius(cfg)
    model=MovingCoupled(RadialMesh(meta["n"],.02,2.),cfg,env,rad)
    t=cache["time_s"];R=cache["radius_current_m"]
    assert np.array_equal(R,rad.value(t))
    valid=cache["radius_fixed_m"][None,:]<=R[:,None]
    for field in ("temperature","moisture"):
        assert np.array_equal(np.isfinite(cache[field][:,:21]),valid)
        assert np.isfinite(cache[field][:,-1]).all()
    assert np.array_equal(cache["profile_state"][0],np.tile([28.,2.55],1280))
    assert np.min(cache["profile_state"][:,1::2])>0
    assert np.max(cache["moisture_global_max"])<=2.55+1e-10
    assert np.max(cache["moisture_radial_increase_max"])<1e-9
    st=meta["statistics"];event=st["event_time_s"]
    assert event<259200. and st["post_event_g"][0]>0 and max(st["post_event_g"][2:])<0
    assert abs(st["post_event_g"][1])<1e-11
    pt=cache["profile_time_s"];py=cache["profile_state"]
    points,nodes=reference_nodes(model,pt,py)
    # Sample physical-domain masks at exact radius knots as well as all exports.
    physical_mask_checks=int(valid.size)
    initial_surface_duplicate=np.max(np.abs(cache["moisture"][0,20]-cache["moisture"][0,21]))
    final_surface_duplicate=np.max(np.abs(cache["moisture"][-1,12]-cache["moisture"][-1,21]))
    assert initial_surface_duplicate<1e-13 and final_surface_duplicate<1e-13
    # Independent P1 nodes versus primary material reconstruction on dense profiles.
    label="q4_fem_n2560_dt1.875"
    folder=ROOT/"results/cache/q4_independent"
    im=json.loads((folder/(label+".json")).read_text())
    assert im["source_fingerprint"]==independent_fingerprint()
    assert im["cache_sha256"]==digest(folder/(label+".npz"))
    with np.load(folder/(label+".npz")) as raw: fe={k:raw[k] for k in raw.files}
    common,i,j=np.intersect1d(pt[:-1],fe["time_s"][:-1],return_indices=True)
    dense={}
    for field,values in zip(("temperature","moisture"),nodes):
        sampled=PchipInterpolator(points,values[i],axis=1)(fe["nodes_xi"])
        near=fe["nodes_xi"]<points[1]
        x0,x1=points[1:3]
        sampled[:,near]=values[i,0,None]+((values[i,2]-values[i,1])/(x1*x1-x0*x0))[:,None]*fe["nodes_xi"][None,near]**2
        delta=np.abs(sampled-fe[field+"_nodes"][j])
        loc=np.unravel_index(np.argmax(delta),delta.shape)
        dense[field]={"max_abs":float(delta[loc]),"time_s":float(common[loc[0]]),"xi":float(fe["nodes_xi"][loc[1]])}
    assert all(r["max_abs"]<5e-5 for r in dense.values())
    # Independent comparison at their separate terminal times is explicitly
    # not a same-time field comparison; it checks terminal shape consistency.
    terminal={}
    for f,values in zip(("temperature","moisture"),nodes):
        sampled=PchipInterpolator(points,values[-1])(fe["nodes_xi"])
        terminal[f]=float(np.max(np.abs(sampled-fe[f+"_nodes"][-1])))
    # Full reconstructed output versus piecewise-linear in physical coordinates.
    # Profiles rather than all seconds are used, consistently with Q4 output frequency.
    linear_diffs={}
    for f,values in zip(("temperature","moisture"),nodes):
        diffs=[]
        for time,vals in zip(pt,values):
            idx=np.flatnonzero(t==time)
            if not len(idx):continue
            physical=cache["radius_fixed_m"];inside=physical<=float(rad.value(time))
            linear=np.interp(physical[inside]/rad.value(time),points,vals)
            diffs.extend(np.abs(cache[f][idx[0],:21][inside]-linear).tolist())
        linear_diffs[f]=max(diffs)
    assert max(linear_diffs.values())<5e-5
    T,C=cache["temperature"][-1,[0,-1]],cache["moisture"][-1,[0,-1]]
    b,k,_,_=model.thermal(C);D=model.diffusion(C,T)[0]
    ta,ce=model.ambient(event)
    air=env.value(pt,"temperature")
    air_delta=np.max(np.abs(nodes[0]-air[:,None]),axis=1)
    iso=None
    for i in np.flatnonzero(air_delta<.1):
        if np.all(air_delta[i:]<.1) and pt[-1]-pt[i]>=1800:iso=float(pt[i]);break
    interval=np.minimum(np.searchsorted(rad.times,pt,side="right")-1,len(rad.times)-2)
    slopes=np.diff(rad.radii)/np.diff(rad.times)
    radius_profile=rad.value(pt);rdot=slopes[interval]
    centreC=nodes[1][:,0];surfaceC=nodes[1][:,-1]
    centreD=model.diffusion(centreC,nodes[0][:,0])[0]
    surfaceD=model.diffusion(surfaceC,nodes[0][:,-1])[0]
    shrink_scale=np.full_like(radius_profile,np.inf)
    np.divide(radius_profile,np.abs(rdot),out=shrink_scale,where=rdot!=0)
    mechanisms={"source_fingerprint":fingerprint(),"script_sha256":digest(__file__),
        "model_conditional_not_experimental":True,"critical_time_s":event,"critical_time_h":event/3600,
        "radius_at_event_cm":100*R[-1],"radius_ratio_at_event":R[-1]/.02,
        "diffusion_geometric_multiplier_at_event":(.02/R[-1])**2,
        "surface_to_volume_ratio_multiplier_at_event":.02/R[-1],
        "dry_density_ratio_at_event":(.02/R[-1])**2,
        "centre_temperature_C":float(T[0]),"surface_temperature_C":float(T[1]),
        "centre_moisture":float(C[0]),"surface_moisture":float(C[1]),
        "mean_moisture":float(cache["moisture_mean"][-1]),
        "fraction_initial_water_removed":float(1-cache["moisture_mean"][-1]/2.55),
        "centre_D_m2_s":float(D[0]),"surface_D_m2_s":float(D[1]),
        "centre_to_surface_D_ratio":float(D[0]/D[1]),"alpha_over_D":(k/b/D).tolist(),
        "event_log_temperature_D_contribution":(3850*(1/301.15-1/(T+273.15))).tolist(),
        "event_log_moisture_D_contribution":(.30*(1/2.55-1/C)).tolist(),
        "internal_moisture_drop_fraction":float((C[0]-C[1])/(C[0]-ce)),
        "normalized_water_flux_m_s":float(model.physical_hm*(C[1]-ce)),
        "mean_moisture_rate_per_s":float(-2*model.physical_hm/R[-1]*(C[1]-ce)),
        "sampled_near_isothermal_since_s":iso,"temperature_PDE_retained":True,
        "maximizer_is_centre_for_t_ge_6h":bool(np.all(cache["moisture_argmax_xi"][t>=21600]==0)),
        "profile_scales":{"time_s":pt.tolist(),"R_m":radius_profile.tolist(),"Rdot_m_s":rdot.tolist(),
             "shrink_time_s":[float(v) if np.isfinite(v) else None for v in shrink_scale],
             "centre_R2_over_D_s":(radius_profile**2/centreD).tolist(),
             "surface_R2_over_D_s":(radius_profile**2/surfaceD).tolist(),"max_air_temperature_gap_C":air_delta.tolist()}}
    save_json(ROOT/"reports/q4_mechanisms.json",mechanisms)
    fine,fm=read_case("q4_grid2560");step,sm=read_case("q4_cap15")
    target=np.floor(event/60)*60
    i=np.flatnonzero(t==target)[0];j=np.flatnonzero(fine["time_s"]==target)[0]
    local=cache["moisture_global_max"][i]-fine["moisture_global_max"][j]
    spatial=event-fm["statistics"]["event_time_s"];slope=st["event_slope_per_s"]
    report={"source_fingerprint":fingerprint(),"script_sha256":digest(__file__),
        "status":"passed_Q4_dense_independent_profiles_geometry_and_error_diagnostics",
        "mask_entries_checked":physical_mask_checks,"fixed_1p2cm_matches_final_surface":True,
        "dense_independent":{"time_s":common.tolist(),"radial_nodes":len(fe["nodes_xi"]),"fields":dense},
        "terminal_shape_differences_at_separate_critical_times":terminal,
        "pchip_vs_linear_at_exportable_profiles":linear_diffs,
        "numerical_event_differences_s":{"1280_to_2560":-spatial,"step30_to_15":sm["statistics"]["event_time_s"]-event,
             "independent_fine_minus_primary":im["event_time_s"]-event},
        "conditioning":{"slope_per_s":slope,"comparison_time_s":target,"C_grid_difference":float(local),
             "propagated_time_difference_s":abs(float(local)/slope),"five_e_minus5_C_to_seconds":5e-5/abs(slope)},
        "spatial_error_estimate_s_assuming_second_order":abs(spatial)*4/3,
        "not_rigorous_error_bound":True,"not_statistical_confidence_interval":True,
        "near_agreement_may_include_cancelling_spatial_errors":True,
        "hours_rounding_comparison":[round(v/3600,4) for v in (event,fm["statistics"]["event_time_s"],sm["statistics"]["event_time_s"],im["event_time_s"])],
        "all_export_last_digits_proven_stable":False}
    assert len(set(report["hours_rounding_comparison"]))==1
    assert abs(abs(spatial)-abs(float(local)/slope))<.002
    save_json(ROOT/"reports/q4_diagnostics.json",report)
    print(json.dumps({"diagnostics":report,"mechanisms":{k:v for k,v in mechanisms.items() if k not in ("source_fingerprint","profile_scales")}},ensure_ascii=False),flush=True)


if __name__=="__main__":main()
