"""Geometry effect resolution, local scales and remaining reconstruction checks.

All comparisons are model-conditional. Paired increments can cancel common
radial errors; they are neither high-resolution 2D answers nor error bounds.
"""
import json
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import numpy as np

from src.data_io import ROOT,digest,load_config,save_json,verify_sources
from src.coupled.model import CoupledFVM
from src.events.integration import ScenarioEnvironment
from src.fvm import RadialMesh
from scripts.question2 import read_case as read_q2
from scripts.question3 import read_case as read_q3
from scripts.question4 import read_case as read_q4
from validation.endface_run import fingerprint
from validation.endface_summary import read
from validation.q4_diagnostics import reference_nodes


def main():
    verify_sources();cfg=load_config()
    geo=json.loads((ROOT/"reports/endface_audit.json").read_text())
    assert geo["fingerprint"]==fingerprint() and geo["complete"]
    evidence={"reports/endface_audit.json":digest(ROOT/"reports/endface_audit.json")}
    rows=geo["cases"]
    def pick(g,nr,nz,rt=1e-8,dt=120.):
        return next(r for r in rows if (r["group"],r["nr"],r["nz"],r["rtol"],r["max_step_s"])==(g,nr,nz,rt,dt))
    refinements=[]
    for g,nr,nz,rt,dt,nr2,nz2,rt2,dt2,kind in (
        ("q23",64,32,1e-8,120.,128,64,1e-8,120.,"both_space_axes"),
        ("q23",64,64,1e-8,120.,128,64,1e-8,120.,"radial"),
        ("q23",128,64,1e-8,120.,128,128,1e-8,120.,"axial"),
        ("q23",64,64,1e-8,120.,64,64,1e-9,60.,"time_and_tolerance"),
        ("q4",64,64,1e-8,120.,128,64,1e-8,120.,"radial"),
        ("q4",128,64,1e-8,120.,128,128,1e-8,120.,"axial"),
        ("q4",64,64,1e-8,120.,64,64,1e-10,30.,"time_and_tolerance"),
        ("q4",128,64,1e-8,120.,128,64,1e-10,30.,"time_and_tolerance"),
        ("q4",64,64,1e-10,30.,128,64,1e-10,30.,"radial_at_tight_time")):
        a,b=pick(g,nr,nz,rt,dt),pick(g,nr2,nz2,rt2,dt2)
        refinements.append({"group":g,"kind":kind,"coarse_label":a["exposed_label"],"fine_label":b["exposed_label"],
            "coarse_paired_effect_s":a["paired_event_difference_s"],"fine_paired_effect_s":b["paired_event_difference_s"],
            "change_in_paired_effect_s":b["paired_event_difference_s"]-a["paired_event_difference_s"],
            "change_in_absolute_2D_time_s":b["event_exposed_s"]-a["event_exposed_s"]})
    symmetry_checks=[]
    for label in sorted({r["insulated_label"] for r in rows}):
        a,m=read(label);y=a["state"].reshape(len(a["time_s"]),m["nz"],m["nr"],2)
        diff=float(np.max(abs(y-y[:,:1])))
        assert diff<1e-9
        symmetry_checks.append({"label":label,"axial_state_max_variation":diff})
    # Add the missing fixed-domain reconstruction audit without changing the
    # already frozen Q2/Q3 solver or any previous acceptance artifact.
    reconstruction={};primaries={}
    for name,loader,label in (("Q2_front",read_q2,"q2_final"),("Q23_tail",read_q3,"q3_final")):
        cache,meta=loader(label);primaries[name]=(cache,meta)
        model=CoupledFVM(RadialMesh(meta["n"],.02,2),cfg,ScenarioEnvironment(cfg),"q23")
        pt=cache["profile_time_s"];points,nodes=reference_nodes(model,pt,cache["profile_state"])
        ix=np.searchsorted(cache["time_s"],pt);assert np.array_equal(cache["time_s"][ix],pt)
        result={}
        for field,values in zip(("temperature","moisture"),nodes):
            linear=np.array([np.interp(cache["radius_m"]/.02,points,v) for v in values])
            diff=abs(linear-cache[field][ix]);loc=np.unravel_index(np.argmax(diff),diff.shape)
            result[field]={"max_abs":float(diff[loc]),"time_s":float(pt[loc[0]]),"radius_m":float(cache["radius_m"][loc[1]])}
        assert max(v["max_abs"] for v in result.values())<5e-5
        reconstruction[name]={"profile_count":len(pt),"fields":result,"scope":"saved internal profiles and 21 output radii; not every exported second",
                              "primary_cache_sha256":meta["cache_sha256"]}
    q4,q4m=read_q4("q4_final");primaries["Q4"]=(q4,q4m)
    scales=[];H=.125
    for g in ("q1","q23","q4"):
        p=cfg["property_sets"][g]
        states=[("initial",28.,2.55,.02,1800. if g=="q1" else 10800.)]
        if g!="q1":
            c,m=primaries["Q23_tail" if g=="q23" else "Q4"]
            for j,label in ((0,"event_centre"),(-1,"event_surface")):
                states.append((label,c["temperature"][-1,j],c["moisture"][-1,j],.02 if g=="q23" else .012,m["statistics"]["event_time_s"]))
        for name,T,C,R,t in states:
            rho=p["rho_kg_m3"]["offset"]+p["rho_kg_m3"]["slope_C"]*C
            cp=p["cp_J_kg_K"]["offset"]+p["cp_J_kg_K"]["fraction_C_over_1_plus_C"]*C/(1+C)
            k=p["k_W_m_K"]["offset"]+p["k_W_m_K"]["fraction_C_over_1_plus_C"]*C/(1+C)
            a=k/(rho*cp);d=p["D_m2_s"];D=d["prefactor"]*np.exp(-d["moisture_exponent"]/C-d["thermal_exponent_K"]/(T+273.15))
            scales.append({"group":g,"state":name,"temperature_C":float(T),"moisture":float(C),"R_m":R,
                "alpha_m2_s":float(a),"D_m2_s":float(D),"radial_heat_scale_h":float(R*R/a/3600),
                "axial_heat_scale_h":float(H*H/a/3600),"radial_water_scale_h":float(R*R/D/3600),
                "axial_water_scale_h":float(H*H/D/3600),"axial_to_radial_ratio":(H/R)**2,
                "diagnostic_time_s":float(t),"frozen_coefficient_sqrt_alpha_t_m":float(np.sqrt(a*t)),
                "frozen_coefficient_sqrt_D_t_m":float(np.sqrt(D*t)),
                "note":"local coefficient screen; neither cumulative propagation length nor drying time prediction"})
    # Explicitly show why subtracting a coarse exposed result from the fine
    # 1D baseline would confuse numerical grid bias with end-face effects.
    official={"q23":primaries["Q23_tail"][1]["statistics"]["event_time_s"],"q4":q4m["statistics"]["event_time_s"]}
    resolution=[]
    for row in (pick("q23",128,128),pick("q4",128,64,1e-10,30.)):
        resolution.append({"group":row["group"],"label":row["exposed_label"],"fine_1D_time_s":official[row["group"]],
            "insulated_coarse_minus_fine_1D_s":row["event_insulated_s"]-official[row["group"]],
            "wrong_unmatched_subtraction_s":row["event_exposed_s"]-official[row["group"]],
            "correct_matched_grid_effect_s":row["paired_event_difference_s"],
            "do_not_replace_final_1D_time_with_this_coarse_2D_time":True})
    report={"status":"passed_conditional_endface_resolution_and_reconstruction_diagnostics","complete":True,
        "fingerprint":fingerprint(),"script_sha256":digest(__file__),"evidence_sha256":evidence,
        "primary_reconstruction_helper_sha256":digest(ROOT/"validation/q4_diagnostics.py"),
        "paired_refinements":refinements,"insulated_axial_uniformity":symmetry_checks,
        "reconstruction":reconstruction,"local_scales":scales,"matched_grid_resolution":resolution,
        "endface_error_bound":None,"all_four_decimal_digits_justified_by_geometry":False}
    save_json(ROOT/"reports/endface_diagnostics.json",report)
    print(json.dumps({"refinements":refinements,"reconstruction":reconstruction,"resolution":resolution},ensure_ascii=False),flush=True)


if __name__=="__main__":main()
