#!/usr/bin/env python3
"""Four fixed/shrinking x appendix3/4 interventions with identical fresh states."""
import json
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import numpy as np
from scripts.question3 import read_case as read_q3
from scripts.question2 import read_case as read_q2
from scripts.question4 import fingerprint,read_case
from src.data_io import ROOT,digest,save_json
from validation.q4_convergence import get_case


def main():
    report={"source_fingerprint":fingerprint(),"script_sha256":digest(__file__),"complete":False,
        "interpretation":"controlled_model_scenarios_not_experimental_causal_identification",
        "A":"fixed radius + appendix3","B":"observed shrinkage + appendix3",
        "C":"fixed radius + appendix4","D":"observed shrinkage + appendix4",
        "geometry_support":"Moving scenarios B/D use observed R through 72 h with no extrapolation. Fixed scenarios define R=R0 for their entire counterfactual horizon.",
        "cases":[],"refinements":[]}
    results={}
    for letter,group,geometry,limit in (("A","q23","fixed",259200.),("B","q23","linear",259200.),
                                       ("C","q4","fixed",604800.),("D","q4","linear",259200.)):
        label="q4_final" if letter=="D" else f"factorial_{letter}"
        c,m=get_case(label,n=1280,group=group,geometry=geometry,max_step=30.,rtol=1e-10,limit=limit)
        assert m["statistics"]["event_found"]
        if geometry!="fixed":assert m["statistics"]["event_time_s"]<259200.
        results[letter]=(c,m)
        report["cases"].append({"letter":letter,"label":label,"group":group,"geometry":geometry,
            "event_time_s":m["statistics"]["event_time_s"],"event_h":m["statistics"]["event_time_s"]/3600,
            "event_radius_cm":100*m["statistics"]["radius_end_m"],
            "time_horizon_is_counterfactual_if_fixed_radius":geometry=="fixed",
            "at_72h_max_moisture":float(c["moisture_global_max"][np.flatnonzero(c["time_s"]==259200.)[0]]) if np.any(c["time_s"]==259200.) else None})
        if letter in ("B","C"):
            fine,fm=get_case(f"factorial_{letter}_n2560",n=2560,group=group,geometry=geometry,max_step=30.,rtol=1e-10,limit=limit)
            report["refinements"].append({"letter":letter,"fine_label":fm["label"],
                "event_delta_s":fm["statistics"]["event_time_s"]-m["statistics"]["event_time_s"]})
        save_json(ROOT/"reports/q4_factorial.json",report)
    q3,q3m=read_q3("q3_final");q2,q2m=read_q2("q2_final")
    whole_t=np.r_[q2["time_s"],q3["time_s"][1:]]
    whole={f:np.concatenate((q2[f],q3[f][1:])) for f in ("temperature","moisture")}
    A,Am=results["A"]
    ts,i,j=np.intersect1d(A["time_s"],whole_t,return_indices=True)
    report["A_fresh_fixed_recovery_vs_validated_Q23"]={"event_delta_s":Am["statistics"]["event_time_s"]-q3m["statistics"]["event_time_s"],
        "field_max_differences":{f:float(np.max(np.abs(A[f][i,:21]-whole[f][j]))) for f in whole}}
    t={k:m["statistics"]["event_time_s"] for k,(c,m) in results.items()}
    terms={"geometry_at_appendix3_s":t["B"]-t["A"],"properties_at_fixed_radius_s":t["C"]-t["A"],
           "interaction_s":t["D"]-t["B"]-t["C"]+t["A"],"total_D_minus_A_s":t["D"]-t["A"],
           "geometry_at_appendix4_s":t["D"]-t["C"],"properties_at_shrinking_radius_s":t["D"]-t["B"]}
    terms["decomposition_residual_s"]=terms["geometry_at_appendix3_s"]+terms["properties_at_fixed_radius_s"]+terms["interaction_s"]-terms["total_D_minus_A_s"]
    assert abs(terms["decomposition_residual_s"])<1e-8
    report["decomposition"]=terms;report["complete"]=True
    save_json(ROOT/"reports/q4_factorial.json",report);print(json.dumps(report),flush=True)


if __name__=="__main__":main()
