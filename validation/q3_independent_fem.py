#!/usr/bin/env python3
"""Extend the independent P1 weak form, with its own BDF2 and event solve.

No primary RHS, event, surface reconstruction, or primary state is reused.
Its own verified Q2 nodal checkpoint is the only initial trajectory source.
"""
import argparse
import json
import sys
import time
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import numpy as np
from scipy.optimize import brentq

from src.data_io import ROOT, digest, load_config, save_json, verify_sources
from validation.q2_independent_fem import IndependentFEM, run as run_initial


def fingerprint():
    paths = [Path(__file__),ROOT/"validation/q2_independent_fem.py",ROOT/"config/model_config.json"]
    return {str(p.relative_to(ROOT)):digest(p) for p in paths}


def run(elements=640, max_step=30., initial_step=.9375):
    started = time.perf_counter()
    verify_sources()
    cfg = load_config()
    model = IndependentFEM(elements, cfg)
    parent_label = f"fem_bdf2_n{elements}_dt{initial_step:g}"
    parent_folder = ROOT/"results/cache/q2_independent"
    path, mp = parent_folder/(parent_label+".npz"), parent_folder/(parent_label+".json")
    try:
        im = json.loads(mp.read_text())
        assert im["script_sha256"] == digest(ROOT/"validation/q2_independent_fem.py")
        assert im["config_sha256"] == digest(ROOT/"config/model_config.json")
        assert im["cache_sha256"] == digest(path)
        with np.load(path) as raw:
            initial = {k:raw[k] for k in raw.files}
    except (FileNotFoundError, AssertionError):
        initial, im = run_initial(elements, initial_step, method="bdf2")
    assert np.array_equal(initial["nodes_m"],model.nodes) and initial["time_s"][-1] == 10800.
    T, C = initial["temperature_nodes"][-1].copy(), initial["moisture_nodes"][-1].copy()
    original_air = model.air.copy()
    platform = original_air[original_air[:,0]>=10800,1:].mean(axis=0)
    threshold = cfg["events"]["moisture_threshold"]
    t = 10800.
    older, previous_dt = None, None
    times, Ts, Cs = [t], [T.copy()], [C.copy()]
    step_count = iterations = 0
    minC = float(C.min())
    record_times = np.arange(14400.,604801.,1800.)
    cuts = np.unique(np.r_[original_air[(original_air[:,0]>10800),0],record_times])
    event_time = event_nodes = None
    event_bracket = None
    accepted_max_dt = 0.
    for stop in cuts:
        if t >= 14400.:
            # This array is local to the independent implementation. No call to
            # the primary boundary extension object is made.
            model.air = np.array([[14400.,*platform],[604800.,*platform]])
        while t < stop-1e-8:
            dt = min(max_step, stop-t)
            if previous_dt is None:
                dt = min(dt, .01)
            else:
                dt = min(dt, 1.8*previous_dt)
            newT,newC,count = model.step(t+dt,dt,T,C,older=older,previous_dt=previous_dt)
            accepted_max_dt = max(accepted_max_dt,dt)
            iterations += count; step_count += 1
            minC = min(minC,float(newC.min()))
            if np.max(newC) <= threshold:
                # Piecewise-linear FE maximum is exactly a nodal maximum.
                def residual(h):
                    if h == 0:
                        return float(np.max(C)-threshold)
                    return float(np.max(model.step(t+h,h,T,C,older=older,previous_dt=previous_dt)[1])-threshold)
                tau = brentq(residual,0.,dt,xtol=1e-6)
                event_time = t+tau
                newT,newC,_ = model.step(event_time,tau,T,C,older=older,previous_dt=previous_dt)
                event_nodes = (newT,newC)
                event_bracket = [t,t+dt]
                T,C = newT,newC; t = event_time
                break
            older,previous_dt = (T,C),dt
            T,C = newT,newC
            t += dt
        if event_time is not None:
            times.append(t); Ts.append(T.copy()); Cs.append(C.copy())
            break
        t = float(stop)
        if stop in record_times:
            times.append(t); Ts.append(T.copy()); Cs.append(C.copy())
        if stop == 14400.:
            # Reset multistep history at the boundary jump; BE startup is refined.
            older,previous_dt = None,None
        if stop % 43200 == 0:
            print(f"FEM n={elements}, dt={max_step}: {stop/3600:.0f} h, Cmax={C.max():.7f}",flush=True)
    if event_time is None:
        raise RuntimeError("Independent horizon exhausted before event")
    radii = np.arange(21)*.001
    cache = {"time_s":np.asarray(times),"radius_m":radii,"nodes_m":model.nodes,
        "temperature_nodes":np.asarray(Ts),"moisture_nodes":np.asarray(Cs),
        "temperature":np.array([np.interp(radii,model.nodes,row) for row in Ts]),
        "moisture":np.array([np.interp(radii,model.nodes,row) for row in Cs])}
    label = f"q3_fem_n{elements}_dt{max_step:g}_initial{initial_step:g}"
    folder = ROOT/"results/cache/q3_independent"
    folder.mkdir(parents=True,exist_ok=True)
    path = folder/(label+".npz")
    np.savez_compressed(path,**cache)
    meta = {"label":label,"source_fingerprint":fingerprint(),"cache_sha256":digest(path),
        "initial_parent_label":parent_label,"initial_parent_sha256":im["cache_sha256"],
        "elements":elements,"max_step_s":max_step,"initial_step_s":initial_step,
        "event_time_s":event_time,"event_bracket_s":event_bracket,
        "critical_maximum":float(C.max()),"maximum_radius_m":float(model.nodes[np.argmax(C)]),
        "minimum_accepted_moisture":minC,"steps":step_count,"picard_iterations":iterations,
        "accepted_max_step_s":accepted_max_dt,"elapsed_s":time.perf_counter()-started,
        "method":"independent_P1_consistent_mass_BDF2_Picard_and_substep_root",
        "status":"computed_pending_independent_convergence"}
    save_json(folder/(label+".json"),meta)
    print(json.dumps(meta),flush=True)
    return cache,meta


if __name__ == "__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--elements",type=int,default=640)
    parser.add_argument("--max-step",type=float,default=30.)
    parser.add_argument("--initial-step",type=float,default=.9375)
    run(**vars(parser.parse_args()))
