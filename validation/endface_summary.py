"""Assemble matched-grid end-face differences, never relabel them as errors."""
import csv
import json
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import numpy as np
from src.data_io import ROOT,digest,save_json,verify_sources
from validation.endface_run import fingerprint


def read(label):
    p=ROOT/"results/cache/endface"/(label+".npz");m=json.loads(p.with_suffix(".json").read_text())
    assert m["fingerprint"]==fingerprint() and m["cache_sha256"]==digest(p)
    return dict(np.load(p)),m


def contrast(group,nr,nz,rtol=1e-8,dt=120.):
    el=f"{group}_r{nr}_z{nz}_e1_tol{rtol:g}_dt{dt:g}"
    bl=f"{group}_r{nr}_z3_e0_tol{rtol:g}_dt{dt:g}"
    a,am=read(el);b,bm=read(bl)
    assert (am["nr"],bm["nr"],am["group"],bm["group"])==(nr,nr,group,group)
    assert (am["rtol"],bm["rtol"],am["max_step_s"],bm["max_step_s"])==(rtol,rtol,dt,dt)
    assert am["initialization"]==bm["initialization"]=="fresh"
    assert am["requested_limit"] is None and bm["requested_limit"] is None
    ts,ia,ib=np.intersect1d(a["time_s"],b["time_s"],return_indices=True)
    fields={}
    for f in ("temperature","moisture"):
        diff=a[f+"_mid"][ia]-b[f+"_mid"][ib]
        assert np.array_equal(np.isfinite(a[f+"_mid"][ia]),np.isfinite(b[f+"_mid"][ib]))
        where=np.unravel_index(np.nanargmax(abs(diff)),diff.shape)
        early=ts<=10800.
        enddiff=a[f+"_end"][ia,:20]-b[f+"_mid"][ib,:20]
        fields[f]={"midplane_max_abs":float(np.nanmax(abs(diff))),"time_s":float(ts[where[0]]),
                   "radius_column":int(where[1]),"midplane_first3h_max_abs":float(np.nanmax(abs(diff[early]))),
                   "end_vs_1D_max_abs_excluding_corner":float(np.nanmax(abs(enddiff)))}
    C=a["state"][:,1::2].reshape(len(a["time_s"]),nz,nr)
    radial_increase=float(np.max(np.diff(C,axis=2)));axial_increase=float(np.max(np.diff(C,axis=1)))
    late=a["time_s"]>=21600
    late_central=bool(np.all(a["argmax_node"][late]==0)) if np.any(late) else None
    # The coupled early field need not be monotone in either direction:
    # end heating changes the local diffusivity and surface replenishment.
    # Hence scan the full reconstructed 2D field; never infer a centre-only
    # event from geometric symmetry. The late central maximum is a result.
    if np.any(late):assert late_central
    time_records=[]
    for target in (1800.,10800.,21600.,172800.):
        if target not in ts:continue
        i=np.flatnonzero(ts==target)[0];aa=ia[i];bb=ib[i]
        time_records.append({"time_s":target,"radius_m":float(a["radius_m"][aa]),
            "mean_C_2D":float(a["mean_C"][aa]),"mean_C_1D":float(b["mean_C"][bb]),
            "mean_C_difference":float(a["mean_C"][aa]-b["mean_C"][bb]),
            "centre_C_difference":float(a["moisture_mid"][aa,0]-b["moisture_mid"][bb,0]),
            "end_centre_temperature_C":float(a["temperature_end"][aa,0]),
            "middle_centre_temperature_C":float(a["temperature_mid"][aa,0]),
            "end_centre_moisture":float(a["moisture_end"][aa,0])})
    delta=None if am["event_time_s"] is None else am["event_time_s"]-bm["event_time_s"]
    return {"group":group,"nr":nr,"nz":nz,"rtol":rtol,"max_step_s":dt,"exposed_label":el,"insulated_label":bl,
        "event_exposed_s":am["event_time_s"],"event_insulated_s":bm["event_time_s"],"paired_event_difference_s":delta,
        "fields":fields,"radial_C_monotonicity_violation":radial_increase,"axial_C_monotonicity_violation":axial_increase,
        "sampled_argmax_is_geometric_centre_for_t_ge_6h":late_central,"snapshots":time_records,
        "water_balance_relative":am["water_balance_relative"],"common_times":len(ts)}


def main():
    verify_sources();rows=[]
    for group in ("q1","q23","q4"):
        for nr,nz in ((32,32),(64,32),(64,64),(128,64)):
            rows.append(contrast(group,nr,nz))
    rows.append(contrast("q4",64,64,1e-10,30.))
    # Additional declared axial/time refinements, when present, remain visible.
    for group,nr,nz,rt,dt in (("q23",128,128,1e-8,120.),("q4",128,128,1e-8,120.),
                              ("q23",64,64,1e-9,60.),("q4",128,64,1e-10,30.)):
        a=ROOT/f"results/cache/endface/{group}_r{nr}_z{nz}_e1_tol{rt:g}_dt{dt:g}.json"
        if a.exists():rows.append(contrast(group,nr,nz,rt,dt))
    labels=set()
    for row in rows:labels.update((row["exposed_label"],row["insulated_label"]))
    artifact={}
    for label in labels:
        for ext in ("json","npz"):
            p=ROOT/f"results/cache/endface/{label}.{ext}";artifact[str(p.relative_to(ROOT))]=digest(p)
    report={"status":"completed_conditional_endface_assumption_audit_not_experimental_validation",
        "fingerprint":fingerprint(),"script_sha256":digest(__file__),"cases":rows,"cache_artifacts_sha256":artifact,
        "geometry":{"R0_m":.02,"half_length_m":.125,"axial_to_radial_time_scale_ratio_initial":(.125/.02)**2,
                    "end_fraction_initial":.02/(.25+.02),"end_fraction_Q4_event":.012/(.25+.012)},
        "assumptions":["same ambient on exposed side and ends","end h and hm equal side h and hm for exposed scenario",
                       "isotropic k and D","fixed length and symmetric ends","Q4 affine radial material motion, no axial shrinkage"],
        "limitations":["paired coarse-grid increments are not rigorous model-error bounds","corner field values are not exported as true surface data",
                       "unknown actual end exposure and anisotropy remain uncalibrated","official 1D workbooks unchanged"],"complete":True}
    save_json(ROOT/"reports/endface_audit.json",report)
    out=ROOT/"tables/endface_comparison.csv"
    with out.open("w",encoding="utf-8-sig",newline="") as f:
        w=csv.writer(f);w.writerow(["group","nr","nz","rtol","max_step_s","paired_event_difference_s","max_mid_T_difference_C","max_mid_C_difference","water_balance_relative"])
        for r in rows:w.writerow([r[k] for k in ("group","nr","nz","rtol","max_step_s","paired_event_difference_s")]+[r["fields"]["temperature"]["midplane_max_abs"],r["fields"]["moisture"]["midplane_max_abs"],r["water_balance_relative"]])
    print(json.dumps([{k:v for k,v in r.items() if k in ("group","nr","nz","rtol","paired_event_difference_s","fields","snapshots")} for r in rows],ensure_ascii=False,indent=2))


if __name__=="__main__":main()
