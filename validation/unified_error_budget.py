"""Evidence-indexed error ledger; differences are not independent variances."""
import csv
import json
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from src.data_io import ROOT,digest,save_json,verify_sources
from validation.endface_run import fingerprint as geometry_fingerprint


def main():
    inputs=verify_sources();provenance={};checks=0
    for q in (1,2,3,4):
        p=ROOT/f"reports/q{q}_export_audit.json";audit=json.loads(p.read_text());provenance[str(p.relative_to(ROOT))]=digest(p)
        for rel,h in audit["artifacts_sha256"].items():assert digest(ROOT/rel)==h,rel;checks+=1
    def read(name):
        p=ROOT/f"reports/{name}.json";provenance[str(p.relative_to(ROOT))]=digest(p)
        return json.loads(p.read_text())
    q1=read("q1_time_and_sensitivity");q1d=read("q1_diagnostics");q1i=read("q1_nonlinear_independent")
    q2=read("q2_convergence");q2i=read("q2_independent_verification");q2s=read("q2_sensitivity")
    q3=read("q3_convergence");q3e=read("q3_error_budget");q3i=read("q3_independent_verification")
    q4=read("q4_convergence");q4d=read("q4_diagnostics");q4i=read("q4_independent_verification");q4s=read("q4_sensitivity")
    geo=read("endface_audit");kernel=read("endface_kernel_checks");geometry_diagnostics=read("endface_diagnostics")
    for r in (q1i,geo):assert r["complete"]
    for r in (geo,kernel,geometry_diagnostics):assert r["fingerprint"]==geometry_fingerprint()
    assert geo["script_sha256"]==digest(ROOT/"validation/endface_summary.py")
    assert kernel["script_sha256"]==digest(ROOT/"validation/endface_kernel_checks.py")
    assert kernel["heat_reference_script_sha256"]==digest(ROOT/"validation/endface_heat_reference.py")
    assert geometry_diagnostics["script_sha256"]==digest(ROOT/"validation/endface_diagnostics.py")
    for rel,h in geometry_diagnostics["evidence_sha256"].items():assert digest(ROOT/rel)==h,rel
    for rel,h in geo["cache_artifacts_sha256"].items():assert digest(ROOT/rel)==h,rel
    for rel,h in q1i["fingerprint"].items():assert digest(ROOT/rel)==h,rel
    for row in q1i["cases"]:
        assert digest(ROOT/f"results/cache/q1_independent/{row['label']}.npz")==row["cache_sha256"]
    assert digest(ROOT/"results/cache/q1/q1_final.npz")==q1i["primary_cache_sha256"]
    records=[]
    def add(q,category,quantity,value,unit,comparison,scope,source,kind="observed_difference_not_bound"):
        records.append({"question":q,"category":category,"quantity":quantity,"value":value,"unit":unit,
                        "comparison":comparison,"scope":scope,"source":source,"interpretation":kind})
    select=lambda rows,label:next(r for r in rows if r["label"]==label)
    scopes={"Q1":"0-1800 s; 21 radii, every second unless stated",
            "Q2_front":"0-10800 s; 21 radii, every second unless stated",
            "Q23_tail":"10800 s to common pre-event times; 21 radii, every second unless stated",
            "Q4":"common 60 s times plus diagnostic 0/1/10 s; valid fixed radii and true surface"}
    datasets=[("Q1",select(q1["comparisons"],"fine2560")["comparison_to_reference"],select(q1["comparisons"],"cap_0p5")["comparison_to_reference"],"q1_time_and_sensitivity"),
              ("Q2_front",q2["mesh_series"][-1]["difference_from_previous"],select(q2["temporal"],"q2_cap1")["difference_from_baseline"],"q2_convergence"),
              ("Q23_tail",q3["mesh_series"][-1]["difference_from_previous"]["fields"],select(q3["temporal"],"q3_cap15")["difference_from_baseline"]["fields"],"q3_convergence"),
              ("Q4",q4["mesh_series"][-1]["difference_from_previous"]["fields"],select(q4["temporal"],"q4_cap15")["difference"]["fields"],"q4_convergence")]
    for q,space,time,source in datasets:
        for f,unit in (("temperature","degC"),("moisture","kg_water_per_kg_dry")):
            add(q,"numerical_space",f,space[f]["max_abs"],unit,"1280 vs 2560 cells",scopes[q],source)
            add(q,"numerical_time",f,time[f]["max_abs"],unit,"halve effective max_step at fixed mesh",scopes[q],source)
            count=space[f].get("four_decimal_rounding_differences",space[f].get("rounding_differences"))
            add(q,"format_rounding",f,count,"cells","four-decimal values changed by mesh refinement",scopes[q],source,"format_stability_count_not_error_magnitude")
    for q,values,scope,source in (
        ("Q1",q1i["dense_profiles"]["max_abs"],"10 saved times; 2561-node dense profiles","q1_nonlinear_independent"),
        ("Q2_front",q2i["dense_paper_profile_comparison"]["max_abs"],"six paper times; 1281-node dense profiles","q2_independent_verification"),
        ("Q23_tail",q3i["cases"][-1]["sampled_field_max_difference"],"common saved time profiles; independent 2560-element trajectory","q3_independent_verification"),
        ("Q4",{f:v["max_abs"] for f,v in q4i["cases"][-1]["sampled_comparisons"].items()},"listed common times/physical radii including early boundary layer","q4_independent_verification")):
        for f,unit in (("temperature","degC"),("moisture","kg_water_per_kg_dry")):
            add(q,"independent_implementation",f,values[f],unit,"primary FVM vs independent nodal FEM/BDF2",scope,source)
    add("Q1","independent_analytic","temperature",q1d["heat_bessel_full_export_max_abs_C"],"degC","primary vs 512 Fourier-Bessel modes",scopes["Q1"],"q1_diagnostics")
    for f,unit in (("temperature","degC"),("moisture","kg_water_per_kg_dry")):
        add("Q1","output_reconstruction",f,q1d["fields"][f]["pchip_vs_piecewise_linear_export_max"],unit,"PCHIP vs piecewise linear",scopes["Q1"],"q1_diagnostics")
        add("Q4","output_reconstruction",f,q4d["pchip_vs_linear_at_exportable_profiles"][f],unit,"PCHIP vs piecewise linear","saved profiles, not every internal solver time","q4_diagnostics")
        for q,row in geometry_diagnostics["reconstruction"].items():
            add(q,"output_reconstruction",f,row["fields"][f]["max_abs"],unit,"PCHIP vs piecewise linear",row["scope"],"endface_diagnostics")
        # These are alternative face-flux implementations, not new physics.
        flux=select(q1["comparisons"],"kirchhoff")["comparison_to_reference"]
        add("Q1","spatial_flux_variant",f,flux[f]["max_abs"],unit,"harmonic vs Kirchhoff interior-face flux",scopes["Q1"],"q1_time_and_sensitivity")
    for q,field,stat,key,source in (
        ("Q1","moisture",q1d["fields"]["moisture"],"independent_balance_relative","q1_diagnostics"),
        ("Q1","temperature",q1d["fields"]["temperature"],"independent_balance_relative","q1_diagnostics"),
        ("Q2_front","moisture",q2["baseline_statistics"],"water_balance_relative","q2_convergence"),
        ("Q2_front","effective_heat",q2["baseline_statistics"],"heat_corrected_balance_relative","q2_convergence"),
        ("Q23_tail","moisture",q3["baseline_statistics"],"water_balance_relative","q3_convergence"),
        ("Q23_tail","effective_heat",q3["baseline_statistics"],"heat_corrected_balance_relative","q3_convergence")):
        add(q,"balance_residual",field,stat[key],"normalized","trajectory quadrature of model boundary flux/corrected heat",scopes[q],source,"balance_residual_not_solution_error")
    q4path=ROOT/"results/cache/q4/q4_final.json";q4meta=json.loads(q4path.read_text())
    provenance[str(q4path.relative_to(ROOT))]=digest(q4path)
    for rel,h in q4meta["source_fingerprint"].items():assert digest(ROOT/rel)==h,rel
    assert digest(q4path.with_suffix(".npz"))==q4meta["cache_sha256"]
    for field,key in (("moisture","water_balance_relative"),("effective_heat","heat_corrected_balance_relative")):
        add("Q4","balance_residual",field,q4meta["statistics"][key],"normalized","moving material balance with variable-capacity correction",scopes["Q4"],"results/cache/q4/q4_final.json","balance_residual_not_solution_error")
        add("Q23_full","balance_segment_triangle_envelope",field,q2["baseline_statistics"][key]+q3["baseline_statistics"][key],
            "normalized","front max residual plus tail max residual","triangle envelope over recorded segment checks",
            "q2_convergence + q3_convergence","bound_on_recorded_balance_residual_only_not_solution_error")
    for q,grid,step,ind,estimate,slope,source in (
        ("Q3",q3e["numerical_differences_s"]["primary_1280_2560"],q3e["numerical_differences_s"]["primary_step30_step15"],q3e["numerical_differences_s"]["independent_fine_minus_primary"],q3e["asymptotic_estimates"]["primary_1280_spatial_error_scale_s_if_second_order"],q3e["event_conditioning"]["slope_per_s"],"q3_error_budget"),
        ("Q4",abs(q4d["numerical_event_differences_s"]["1280_to_2560"]),q4d["numerical_event_differences_s"]["step30_to_15"],q4d["numerical_event_differences_s"]["independent_fine_minus_primary"],q4d["spatial_error_estimate_s_assuming_second_order"],q4d["conditioning"]["slope_per_s"],"q4_diagnostics")):
        for cat,val,label in (("numerical_space",grid,"1280 vs 2560"),("numerical_time",step,"max_step 30 vs 15 s"),("independent_implementation",ind,"fine independent vs primary"),("asymptotic_estimate",estimate,"assuming second order")):
            add(q,cat,"critical_time",val,"s",label,"global-max threshold event",source)
        add(q,"event_conditioning","slope",slope,"kg_per_kg_per_s","local slope near stable maximizer","threshold neighborhood",source)
        add(q,"event_conditioning","time_equivalent_of_C_5e-5",5e-5/abs(slope),"s","delta_t approx delta_C / abs(slope)","not an observed event error",source)
    # A near-zero threshold residual describes the root on the numerical
    # trajectory only. Record the actual strict-inequality checks separately.
    for q,stat,integer_key in (("Q3",q3["baseline_statistics"],"first_checked_integer_second"),
                               ("Q4",q4meta["statistics"],"first_strict_integer_second")):
        source="q3_convergence" if q=="Q3" else "results/cache/q4/q4_final.json"
        add(q,"root_residual","max_C_minus_threshold",stat["post_event_g"][1],"kg_water_per_kg_dry",
            "g on numerical event state","not discretization or model accuracy",source,"floating_point_root_residual_not_total_error")
        add(q,"discrete_execution_time","first_verified_strict_integer_second",stat[integer_key],"s",
            "ceil critical time then explicitly verify max_C < 0.15","selected numerical model",source,"verified_integer_time_not_confidence_bound")
    # Q1/Q2 also need model-scenario evidence in a four-question budget.
    for row in q1["sensitivities"]:
        label=f"{row['parameter']}={row['multiplier']}"
        for f,unit in (("temperature","degC"),("moisture","kg_water_per_kg_dry")):
            add("Q1","model_parameter_or_boundary",f,row["comparison_to_baseline"][f]["max_abs"],unit,label,
                scopes["Q1"],"q1_time_and_sensitivity","conditional_max_field_difference_not_confidence_interval")
    # Store absolute scenario observables, distinctly labelled; do not invent
    # a difference by subtracting a mismatched mean/centre quantity.
    for row in q2s["cases"]:
        for key,unit in (("temperature_centre_3h","degC"),("temperature_surface_3h","degC"),
                         ("moisture_centre_3h","kg_water_per_kg_dry"),("moisture_surface_3h","kg_water_per_kg_dry"),
                         ("moisture_mean_3h","kg_water_per_kg_dry")):
            add("Q2_front","model_parameter_or_boundary",key,row[key],unit,row["label"],"state at 10800 s", "q2_sensitivity",
                "conditional_absolute_scenario_value_not_error")
    for row in q3e["model_scenarios"]:
        add("Q3","model_parameter_or_boundary","critical_time",row["event_delta_s"],"s",row["label"],"one-factor/frozen-feedback scenarios","q3_error_budget","conditional_scenario_not_confidence_interval")
    for row in q4s["cases"]:
        add("Q4","model_parameter_or_boundary","critical_time",row["event_time_s"]-q4d["critical_time_s"] if "critical_time_s" in q4d else row["event_time_s"]-q4meta["statistics"]["event_time_s"],"s",row["label"],"one-factor/geometry interpolation/environment scenarios","q4_sensitivity","conditional_scenario_not_confidence_interval")
    for r in geo["cases"]:
        scope=f"matched nr={r['nr']},nz={r['nz']},rtol={r['rtol']}; common saved times"
        if r["paired_event_difference_s"] is not None:
            add("Q3" if r["group"]=="q23" else "Q4","model_endface","critical_time",r["paired_event_difference_s"],"s","exposed ends minus insulated ends",scope,"endface_audit","paired_scenario_difference_not_rigorous_error_bound")
        for f,unit in (("temperature","degC"),("moisture","kg_water_per_kg_dry")):
            add(r["group"],"model_endface",f,r["fields"][f]["midplane_max_abs"],unit,"midplane exposed vs insulated",scope,"endface_audit","conditional_geometry_difference")
        for snapshot in r["snapshots"]:
            add(r["group"],"model_endface","mean_moisture",snapshot["mean_C_difference"],"kg_water_per_kg_dry",
                "whole-cylinder average: exposed minus insulated",scope+f"; at {snapshot['time_s']} s","endface_audit","conditional_geometry_difference")
    for r in geometry_diagnostics["paired_refinements"]:
        add("Q3" if r["group"]=="q23" else "Q4","endface_effect_numerical_resolution","change_in_paired_time_effect",
            r["change_in_paired_effect_s"],"s",r["kind"],r["coarse_label"]+" -> "+r["fine_label"],"endface_diagnostics",
            "resolution_of_scenario_difference_not_bound_on_true_2D_time")
    balance_normalizations={"Q1_temperature":"divide temperature-equivalent heat residual by 22 degC",
        "all_moisture":"divide dry-mass-normalized average-C balance residual by initial C0=2.55",
        "Q23_effective_heat":"divide corrected H residual by initial effective b(C0)*22 K",
        "Q4_effective_heat":"same initial-b normalization; reference/material-coordinate H, not total mixture enthalpy",
        "Q23_tail":"tail integration residual since 10800 s, not a recomputed whole-run residual"}
    unquantified=[{"item":x,"value":None,"reason":y} for x,y in (
        ("real air-to-material sorption mapping","no measured sorption/desorption isotherm or internal moisture observations"),
        ("actual end contact/exposure and axial anisotropy","not specified by problem or attachments; equal-transfer scenario is not calibration"),
        ("latent heat, moisture-carried enthalpy, deformation work","not in selected effective energy model; residual cannot estimate missing physics"),
        ("true property-law/model discrepancy","no experimental reference or uncertainty distribution"),
        ("measurement noise in air histories and prescribed radius","no uncertainty model or repeat observations; interpolation variants do not calibrate measurement error"),
        ("Q4 strict mixture density and length/shrinkage relation","strict total density interpretation conflicts with prescribed radius under fixed length/dry mass"),
        ("continuous whole-space whole-time error bound","finite mesh/time comparisons and samples only"))]
    report={"status":"completed_unified_evidence_ledger_with_explicit_unquantified_model_errors","script_sha256":digest(__file__),
            "evidence_sha256":provenance,"prior_audited_artifacts_unchanged":checks,"original_inputs_unchanged":len(inputs["files"]),
            "records":records,"unquantified":unquantified,"confidence_interval":None,"root_sum_of_squares_permitted":False,
            "balance_normalizations":balance_normalizations,
            "all_four_decimal_digits_physically_validated":False,"complete":True}
    save_json(ROOT/"reports/unified_error_budget.json",report)
    with (ROOT/"tables/unified_error_budget.csv").open("w",encoding="utf-8-sig",newline="") as f:
        w=csv.DictWriter(f,fieldnames=list(records[0]));w.writeheader();w.writerows(records)
    print(json.dumps({"rows":len(records),"prior_artifacts_unchanged":checks,"unquantified_categories":len(unquantified)}),flush=True)


if __name__=="__main__":main()
