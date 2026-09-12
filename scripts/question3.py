#!/usr/bin/env python3
"""Q3 event from a provenance-checked Q2 full-precision checkpoint."""
import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import numpy as np

from scripts.question2 import read_case as read_q2, fingerprint as q2_fingerprint
from src.coupled.model import CoupledFVM
from src.data_io import ROOT, digest, load_config, save_json, verify_sources
from src.events.integration import ScenarioEnvironment, continue_to_event
from src.fvm import RadialMesh


def fingerprint():
    fp = q2_fingerprint()
    fp.update({str(p.relative_to(ROOT)): digest(p) for p in
               sorted((ROOT/"src/events").glob("*.py"))+[Path(__file__)]})
    return fp


def read_case(label):
    folder = ROOT/"results/cache/q3"
    path = folder/(label+".npz")
    meta = json.loads((folder/(label+".json")).read_text())
    if meta["source_fingerprint"] != fingerprint() or meta["cache_sha256"] != digest(path):
        raise ValueError("Stale Q3 case: "+label)
    _, parent = read_q2(meta["q2_parent_label"])
    if meta["q2_parent_sha256"] != parent["cache_sha256"]:
        raise ValueError("Q2 continuation checkpoint changed")
    with np.load(path) as raw:
        return {k: raw[k] for k in raw.files}, meta


def run_case(parent="q2_final", label="q3_trial", max_step=60., rtol=1e-10,
             method="BDF", output_step=60., scenario="main", limit=604800.):
    config = load_config()
    manifest = verify_sources()
    if config["problems"]["q3"] != {"properties":"q23", "radius":"fixed",
            "initialization":"q2_full_precision_same_config_or_fresh", "time_origin_s":0,
            "end":"global_moisture_event"}:
        raise ValueError("Q3 timeline contract changed")
    cache, meta = read_q2(parent)
    if (meta["property_set"], meta["initialization"], meta["end_s"], meta["time_origin_s"]) != ("q23", "fresh", 10800., 0.):
        raise ValueError("Wrong Q2 checkpoint timeline/properties")
    assert cache["time_s"][-1] == 10800 and np.array_equal(cache["profile_state"][-1], cache["final_state"])
    mesh = RadialMesh(meta["n"], config["geometry"]["reference_radius_m"], meta["grading"])
    assert np.array_equal(mesh.centres, cache["cell_centres_m"])
    model = CoupledFVM(mesh, config, ScenarioEnvironment(config, scenario), "q23", meta["perturbation"])
    result, stats = continue_to_event(model, cache["final_state"], max_step=max_step,
        rtol=rtol, method=method, output_step=output_step, limit=limit)
    folder = ROOT/"results/cache/q3"
    folder.mkdir(parents=True, exist_ok=True)
    path = folder/(label+".npz")
    np.savez_compressed(path, **result)
    metadata = {"label":label, "property_set":"q23", "physical_time_origin_s":0,
        "continuation_start_s":10800., "q2_parent_label":parent,
        "q2_parent_sha256":meta["cache_sha256"], "n":meta["n"], "grading":meta["grading"],
        "perturbation":meta["perturbation"], "rtol":rtol, "atol_temperature":1e-11,
        "atol_moisture":1e-13, "max_step_s":max_step, "method":method,
        "output_step_s":output_step, "boundary_scenario":scenario,
        "source_fingerprint":fingerprint(), "cache_sha256":digest(path),
        "source_inputs":manifest["files"], "statistics":stats,
        "status":"computed_pending_Q3_convergence_and_independent_verification"}
    save_json(folder/(label+".json"), metadata)
    print(json.dumps({"label":label, "statistics":stats}), flush=True)
    return result, metadata


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--parent", default="q2_final")
    parser.add_argument("--label", default="q3_trial")
    parser.add_argument("--max-step", type=float, default=60.)
    parser.add_argument("--rtol", type=float, default=1e-10)
    parser.add_argument("--method", choices=["BDF", "Radau"], default="BDF")
    parser.add_argument("--output-step", type=float, default=60.)
    parser.add_argument("--scenario", choices=["main", "last_observation", "mean_9000_14400"], default="main")
    parser.add_argument("--limit", type=float, default=604800.)
    run_case(**vars(parser.parse_args()))
