#!/usr/bin/env python3
"""Global maximum search, boundary split and checkpoint continuation tests."""
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))

import numpy as np
from scipy.interpolate import PchipInterpolator

from scripts.question2 import read_case as read_q2
from scripts.question3 import read_case, fingerprint
from src.boundary import Environment
from src.coupled.model import CoupledFVM
from src.coupled.integrate import integrate
from src.data_io import ROOT, digest, load_config, save_json
from src.events.integration import ScenarioEnvironment, global_maximum, continue_to_event
from src.fvm import RadialMesh


def main():
    cfg=load_config()
    env=ScenarioEnvironment(cfg)
    m=CoupledFVM(RadialMesh(40,.02,2.),cfg,env)
    C=.1+.2*np.exp(-((m.mesh.centres-.011)/.002)**2)
    state=np.empty(80); state[::2]=50.; state[1::2]=C
    x0,x1=m.mesh.centres[:2]
    axis=(x1*x1*C[0]-x0*x0*C[1])/(x1*x1-x0*x0)
    surf=m.surface(20000.,state[-2],state[-1])[1]
    nodes=np.r_[axis,C,surf]
    points=np.r_[0.,m.mesh.centres,.02]
    interpolator=PchipInterpolator(points,nodes)
    stationary=interpolator.derivative().roots(extrapolate=False)
    stationary=stationary[np.isfinite(stationary)&(stationary>=x0)&(stationary<=.02)]
    exact=float(np.max(interpolator(np.r_[points,stationary])))
    found=float(global_maximum(m,20000.,state))
    assert abs(found-exact)<1e-13 and found>.29 and axis<.15
    # Manufactured hump would be missed by a centre-only stopping condition.
    env.right_at_join=False
    left=[float(env.value(14400.,f)) for f in ("temperature","moisture")]
    env.right_at_join=True
    right=[float(env.value(14400.,f)) for f in ("temperature","moisture")]
    assert np.allclose(left,env.data[-1,1:]) and np.allclose(right,env.platform)
    rejected=[]
    for name,kwargs in (("reset_clock",{"start":0.}), ("wrong_shape",{"initial_state":state[:-2]}),
                        ("already_crossed",{"initial_state":np.tile([50.,.1],40)})):
        arguments={"model":m,"initial_state":state}; arguments.update(kwargs)
        try:
            continue_to_event(**arguments)
        except ValueError:
            rejected.append(name)
    assert len(rejected)==3
    # Independent start strategy: fresh integration versus a checked 3 h restart.
    tail,meta=read_case("q3_mesh320")
    parent,pm=read_q2(meta["q2_parent_label"])
    assert np.array_equal(tail["profile_state"][0],parent["final_state"])
    model=CoupledFVM(RadialMesh(320,.02,2.),cfg,Environment(cfg))
    fresh,stats=integrate(model,end=12600.,rtol=1e-10,atol_temperature=1e-11,
                         atol_moisture=1e-13,max_step=30.,output_step=60.)
    times,i,j=np.intersect1d(fresh["time_s"],tail["time_s"],return_indices=True)
    difference={f:float(np.max(np.abs(fresh[f][i]-tail[f][j]))) for f in ("temperature","moisture")}
    assert all(v<1e-6 for v in difference.values())
    baseline,bm=read_case("q3_final")
    st=bm["statistics"]
    assert abs(st["post_event_g"][1])<1e-11 and st["post_event_g"][0]>0 and max(st["post_event_g"][2:])<0
    assert np.max(np.diff(baseline["accepted_check_global_max"]))<1e-10
    report={"source_fingerprint":fingerprint(),"script_sha256":digest(__file__),
        "manufactured_interior_hump_global_max":found,"hump_axis_value":float(axis),
        "global_max_difference_from_PCHIP_stationary_point_search":abs(found-exact),
        "rejected_cases":rejected,"4h_left_boundary":left,"4h_right_boundary":right,
        "fresh_vs_checkpoint_3h_to_3p5h_max_differences":difference,
        "post_event_g":st["post_event_g"],"global_max_decreased_at_all_accepted_and_midpoint_checks":True,
        "status":"passed_Q3_global_event_and_timeline_checks_not_a_continuous_error_bound"}
    save_json(ROOT/"reports/q3_event_checks.json",report)
    print(report,flush=True)


if __name__=="__main__":
    main()
