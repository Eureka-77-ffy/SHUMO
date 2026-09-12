#!/usr/bin/env python3
"""Q2: fresh appendix-three coupled trajectory; 3h milestone, not Q3 event."""
import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np

from src.boundary import Environment
from src.coupled.model import CoupledFVM
from src.coupled.integrate import integrate
from src.data_io import ROOT, digest, load_config, save_json, source_fingerprint, verify_sources
from src.fvm import RadialMesh


def fingerprint():
    result = source_fingerprint()
    paths = sorted((ROOT/"src/coupled").glob("*.py"))+[Path(__file__)]
    result.update({str(p.relative_to(ROOT)): digest(p) for p in paths})
    return result


def read_case(label):
    folder = ROOT/"results/cache/q23"
    path = folder/(label+".npz")
    metadata = json.loads((folder/(label+".json")).read_text())
    if metadata["source_fingerprint"] != fingerprint() or metadata["cache_sha256"] != digest(path):
        raise ValueError("Stale or altered Q2 cache: "+label)
    with np.load(path) as raw:
        cache = {key: raw[key] for key in raw.files}
    return cache, metadata


def run_case(n=320, grading=2., rtol=1e-9, atol_temperature=1e-10, atol_moisture=1e-12,
             max_step=10., method="BDF", label=None, perturbation=None):
    config = load_config()
    manifest = verify_sources()
    expected = {"properties": "q23", "radius": "fixed", "initialization": "fresh",
                "end": "q3_event", "paper_end_s": 10800}
    if config["problems"]["q2"] != expected:
        raise ValueError("Q2 timeline/property contract changed")
    env = Environment(config)
    mesh = RadialMesh(n, config["geometry"]["reference_radius_m"], grading)
    model = CoupledFVM(mesh, config, env, "q23", perturbation)
    cache, stats = integrate(model, end=10800., rtol=rtol, atol_temperature=atol_temperature,
                             atol_moisture=atol_moisture, max_step=max_step, method=method)
    label = label or f"q2_n{n}_g{grading:g}_{method}_r{rtol:g}_h{max_step:g}"
    folder = ROOT/"results/cache/q23"
    folder.mkdir(parents=True, exist_ok=True)
    path = folder/(label+".npz")
    np.savez_compressed(path, **cache)
    metadata = {"label": label, "problem": "q2", "property_set": "q23", "initialization": "fresh",
                "time_origin_s": 0., "start_s": 0., "end_s": 10800.,
                "scope": "first_three_hours_only; full_workbook_end_pending_Q3_event",
                "status": "computed_pending_comparative_validation", "n": n, "grading": grading,
                "rtol": rtol, "atol_temperature": atol_temperature, "atol_moisture": atol_moisture,
                "max_step_s": max_step, "method": method, "perturbation": dict(perturbation or {}),
                "boundary_half_cell": "integral_C_at_local_mean_T_and_k_at_local_mean_C",
                "source_fingerprint": fingerprint(), "source_inputs": manifest["files"],
                "cache_sha256": digest(path), "statistics": stats}
    save_json(folder/(label+".json"), metadata)
    print(json.dumps({"label": label, "statistics": stats,
                      "temperature_3h_centre_surface": cache["temperature"][-1, [0, -1]].tolist(),
                      "moisture_3h_centre_surface": cache["moisture"][-1, [0, -1]].tolist()}, indent=2), flush=True)
    return cache, metadata


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--n", type=int, default=320)
    parser.add_argument("--grading", type=float, default=2.)
    parser.add_argument("--rtol", type=float, default=1e-9)
    parser.add_argument("--max-step", type=float, default=10.)
    parser.add_argument("--method", choices=["BDF", "Radau"], default="BDF")
    parser.add_argument("--label")
    args = parser.parse_args()
    run_case(n=args.n, grading=args.grading, rtol=args.rtol, max_step=args.max_step,
             method=args.method, label=args.label)
