#!/usr/bin/env python3
"""Keep discretization, event conditioning and modeling scenarios distinct."""
import json
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import numpy as np

from scripts.question3 import read_case, fingerprint
from src.data_io import ROOT, digest, save_json, verify_sources


def main():
    verify_sources()
    baseline,bm=read_case("q3_final")
    fine,fm=read_case("q3_mesh2560")
    cap,cm=read_case("q3_cap15")
    convergence=json.loads((ROOT/"reports/q3_convergence.json").read_text())
    independent=json.loads((ROOT/"reports/q3_independent_verification.json").read_text())
    scenarios=json.loads((ROOT/"reports/q3_sensitivity.json").read_text())
    for report in (convergence,independent,scenarios):
        assert report["complete"] and report["source_fingerprint"]==fingerprint()
    t=bm["statistics"]["event_time_s"]
    slope=bm["statistics"]["event_slope_per_s"]
    target=60*np.floor(t/60)
    i=np.flatnonzero(baseline["time_s"]==target)[0]
    j=np.flatnonzero(fine["time_s"]==target)[0]
    local_C_delta=float(baseline["moisture_global_max"][i]-fine["moisture_global_max"][j])
    spatial=t-fm["statistics"]["event_time_s"]
    temporal=t-cm["statistics"]["event_time_s"]
    grid_deltas=[abs(r["event_delta_s"]) for r in convergence["mesh_series"][1:]]
    independent_final=independent["cases"][-1]
    error_budget={"source_fingerprint":fingerprint(),"script_sha256":digest(__file__),
        "critical_time_s":t,"critical_time_h":t/3600,
        "event_conditioning":{"slope_per_s":slope,"comparison_time_s":target,
            "local_maximum_C_grid_difference":local_C_delta,
            "linearized_time_difference_s":abs(local_C_delta/slope),
            "five_e_minus5_C_translates_to_s":5e-5/abs(slope),
            "assumptions":"stable maximizer and nonzero downward crossing slope"},
        "numerical_differences_s":{"primary_1280_2560":abs(spatial),
            "primary_step30_step15":abs(temporal),
            "BDF_Radau":abs(next(r["event_delta_s"] for r in convergence["temporal"] if r["label"]=="q3_radau")),
            "independent_fine_minus_primary":independent_final["event_delta_from_primary_s"],
            "independent_step7p5_step3p75":independent["separate_refinements"][-1]["event_delta_s"]},
        "asymptotic_estimates":{"primary_observed_orders":[float(np.log2(a/b)) for a,b in zip(grid_deltas[:-1],grid_deltas[1:])],
            "primary_1280_spatial_error_scale_s_if_second_order":abs(spatial)*4/3,
            "not_a_rigorous_error_bound":True,
            "independent_small_difference_can_include_cancelling_discretization_errors":True},
        "hours_rounding_comparison":{
            "primary1280":round(t/3600,4),"primary2560":round(fm["statistics"]["event_time_s"]/3600,4),
            "step15":round(cm["statistics"]["event_time_s"]/3600,4),
            "independent2560_step3p75":round(independent_final["event_time_s"]/3600,4)},
        "model_scenarios":[{"label":r["label"],"event_delta_s":r["event_time_s"]-t,
                           "relative_delta":r["event_time_s"]/t-1} for r in scenarios["cases"]],
        "unquantified_model_errors":["air-to-material sorption closure beyond beta scenarios", "neglected end faces", "latent heat and complete mixture energy balance", "empirical property uncertainty beyond multiplier scenarios"],
        "confidence_interval":None,"do_not_combine_as_root_sum_squares":True,
        "all_export_last_digits_stable":False}
    assert len(set(error_budget["hours_rounding_comparison"].values()))==1
    assert abs(error_budget["event_conditioning"]["linearized_time_difference_s"]-abs(spatial))<.001
    save_json(ROOT/"reports/q3_error_budget.json",error_budget)
    print(json.dumps(error_budget,ensure_ascii=False,indent=2))


if __name__=="__main__":
    main()
